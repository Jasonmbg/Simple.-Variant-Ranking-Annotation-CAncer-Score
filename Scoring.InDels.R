scoring.func.indels.core <- function(rdata_dir, exp.genes=NULL,
                                     ref.genome=c("hg19","hg38"), 
                                     sample.name.output,
                                     w1=0.7, w2=0.8, w3=0.9, w4=1){
  
  
  # rdata_dir is the full path directory where the relative RData file has been stored
  # exp.genes is a character vector of expressed genes IDs as HUGO gene symbols
  # Load the necessary packages
  
  
  require(tidyverse)
  require(jsonlite)
  require(DT)
  require(htmlwidgets)
  
  ifr <- list.files(path=rdata_dir, ".RData",full.names = T)
  
  OC.input <-get(load(ifr)) 
  
  
  stopifnot(inherits(OC.input,"data.frame"))
  stopifnot(ncol(OC.input)>1)
  stopifnot(!is.null(sample.name.output))
  
  
  ref.genome=match.arg(ref.genome,several.ok = TRUE)
  
  if (ref.genome== "hg19") {
    
    nec_annot_cols <- c("chrom", "pos", "ref_base", "alt_base", "hugo", "transcript", "so", "cchange", "achange", "all_mappings", "samples", "numsample", "hg19__chrom", "hg19__pos","mutpred_indel__score","loftool__loftool_score",
                        "gnomad__af", "cgc__class", "cosmic__variant_count", "cgl__class", "cgl__class", "clinvar__sig")
    
    if(!all(nec_annot_cols%in%colnames(OC.input))) {
      stop("Not all necessary OpenCRAVAT annotator names for creating the relative scores !!")
    } 
    
    OC.input <- OC.input %>% dplyr::select(all_of(nec_annot_cols))
    
    
  }
  
  else if (ref.genome == "hg38") {
    
    nec_annot_cols <- c("chrom", "pos", "ref_base", "alt_base", "hugo", "transcript", "so", "cchange", "achange", "all_mappings", "samples", "numsample", "mutpred_indel__score", "loftool__loftool_score", "gnomad__af", "cgc__class", "cosmic__variant_count", "cgl__class", "cgl__class", "clinvar__sig") 
    
    if(!all(nec_annot_cols%in%colnames(OC.input))) {
      stop("Not all necessary OpenCRAVAT annotator names for creating the relative scores !!")
    }  
    
    OC.input <- OC.input %>% dplyr::select(all_of(nec_annot_cols))  
    
  } 
  
  
  if (is.null(exp.genes)){
    warning(paste0("No vector of gene symbols is provided, thus expression score would be 0"))
    
    d2 <- OC.input
    
  } else if (is.character(exp.genes)){
    
    d2 <- OC.input %>% mutate(Exp_score= 
                                if_else(hugo%in%exp.genes, 1, 0, missing = 0))
    
  } 
  
  
  # Initial creation of pathogenicity/effect prediction score
  
  d3 <- d2 %>% filter(so %in% c("inframe_deletion","frameshift_elongation", 
                                "frameshift_truncation","inframe_insertion")) %>%
    
    dplyr::filter((!gnomad__af > 0.05) %>% replace_na(TRUE)) %>% 
    
    mutate(Pathogenicity_Indel_score=case_when(
      
      so %in% c("frameshift_truncation","frameshift_elongation") & loftool__loftool_score >0.6 ~1,
      
      so %in% c("frameshift_truncation","frameshift_elongation") & loftool__loftool_score <0.6 ~0.8,
      
      so %in% c("inframe_deletion","inframe_insertion") & mutpred_indel__score >0.5 ~0.6,
      
      so %in% c("inframe_deletion","inframe_insertion") & mutpred_indel__score <0.5 ~0.5,
      
      TRUE~0)) 
  
  # Next we create the cancer evidence score
  
  d4 <- d3 %>%
    
    mutate(
      cancer_gene_ev_score=case_when(
        str_detect(cgc__class,"Oncogene|fusion|TSG") & cgl__class%in%c("Oncogene","TSG")~1,
        
        str_detect(cgc__class,"Oncogene|fusion|TSG") | cgl__class%in%c("Oncogene","TSG")~0.5,
        
        TRUE~0),
      
      cancer_freq_score= case_when(
        
        cosmic__variant_count > 10 ~1,
        
        cosmic__variant_count < 10 ~0.5,
        TRUE~0),
      
      clin_ev_score= case_when(
        
        str_detect(clinvar__sig, "Pathogenic|drug response|pathogenic$")~1,
        
        TRUE~0)) %>% mutate(Total_Cancer_score = rowMeans(select(., cancer_gene_ev_score,
                                                                 cancer_freq_score)))
  
  
  if (is.null(exp.genes)){
    
    
    d.final <- d4 %>% mutate(Final_Score= (Pathogenicity_Indel_score*w1 + Total_Cancer_score*w2 + clin_ev_score*w3)/2.4)
    
  } else if (!is.null(exp.genes)){
    
    d.final <- d4 %>% mutate(Final_Score= (Pathogenicity_Indel_score*w1 + Total_Cancer_score*w2 +clin_ev_score*w3 + Exp_score*w4)/3.4)
    
  } 
  
  
  write_csv(d.final,file=str_c("Finaloutput.Variant.Ranked.OpenCRAVAT",
                               sample.name.output,"csv",sep="."))
  
  dat.sel.top10 <- d.final %>% arrange(desc(Final_Score)) %>% top_n(10)
  
  mytable <- datatable(dat.sel.top10, filter = 'top', options = list(paging = FALSE))
  
  saveWidget(mytable, file=str_c("Top10RankedVariants",
                                 sample.name.output,"html",sep="."))
  
}