scoring.func.indels <- function(rdata_dir,exp.genes=NULL,sample.name.output,
                                w1=0.7,w2=0.8,w3=0.9,w4=1){
  
  
  # rdata_dir is the full path directory where the relative RData file has been stored
  # exp.genes is a character vector of expressed genes IDs as HUGO gene symbols
  # Load the necessary packages
  
  
  require(tidyverse)
  require(jsonlite)
  require(DT)
  require(htmlwidgets)
  
  ifr <- list.files(path=rdata_dir, ".RData", full.names = T)
  
  OC.input <-get(load(ifr)) 
  
  
  stopifnot(inherits(OC.input,"data.frame"))
  stopifnot(ncol(OC.input)>1)
  stopifnot(!is.null(sample.name.output))
  
  if(!any(colnames(OC.input)%in%c("chrom", "pos", "ref_base", "alt_base", "hugo,transcript", "so", "cchange", "achange",
                                  "all_mappings", "hg19__chrom", "hg19__pos","mutpred_indel__score",
                                  "loftool__loftool_score","Exp_score,gnomad__af",
                                  "cgc__class","cosmic__variant_count","cgl__class","cgl__class", 
                                  "civic__clinical_a_score","clinvar__sig"))) {
    
    stop("Not necessary OpenCRAVAT annotator names for creating the merged score !!")
  }  
  
  if (is.null(exp.genes)){
    warning(paste0("No vector of gene symbols is provided, thus expression score would be 0"))
    
    d2 <- OC.input %>% mutate(Exp_score=0)
    
  } else if (is.character(exp.genes)){
    
    d2 <- OC.input %>% mutate(Exp_score= if_else(hugo%in%c(exp.genes), 1, 0, missing = 0))
    
  } 
  
  
  # Initial creation of pathogenicity/effect prediction score
  
  d3 <- d2 %>% filter(so %in% c("inframe_deletion","frameshift_elongation", 
                                "frameshift_truncation","inframe_insertion")) %>%
    
    
    dplyr::select(chrom, pos, ref_base, alt_base, hugo,transcript, so, cchange, achange,
                  all_mappings, hg19__chrom, hg19__pos,mutpred_indel__score,
                  loftool__loftool_score,Exp_score,gnomad__af,
                  cgc__class,cosmic__variant_count,cgl__class,cgl__class, 
                  civic__clinical_a_score,clinvar__sig) %>%
    
    
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
        
        !is.na(civic__clinical_a_score) & str_detect(clinvar__sig, "Pathogenic|drug response|pathogenic$")~1,
        
        !is.na(civic__clinical_a_score) | str_detect(clinvar__sig, "Pathogenic|drug response|pathogenic$")~0.5,
        TRUE~0)) %>% mutate(Total_Cancer_score = rowMeans(select(., cancer_gene_ev_score, cancer_freq_score)))
  
  # Final merging of the created sub-scores
  
  weight.func <- function(a,b,c,d,w1=0.7,w2=0.8,w3=0.9,w4=1){
    weighted.score <- (a*w1 + b*w2 + c*w3 + d*w4)/3.4
    return(weighted.score)
  }
  
  # Final merging of the created sub-scores
  
  d.final <- d4 %>% mutate(Final_Score= weight.func(a=Pathogenicity_Indel_score, b=Total_Cancer_score, c=clin_ev_score, d=Exp_score))
  
  write_csv(d.final,file=str_c("Finaloutput.Variant.Ranked.OpenCRAVAT",sample.name.output,"csv",sep="."))
  
  dat.sel.top10 <- d.final %>% arrange(desc(Final_Score)) %>% top_n(10)
  
  mytable <- datatable(dat.sel.top10, filter = 'top', options = list(paging = FALSE))
  
  saveWidget(mytable, file=str_c("Top10.Ranked.OC", sample.name.output, "mytable.html",sep="."))
  
}