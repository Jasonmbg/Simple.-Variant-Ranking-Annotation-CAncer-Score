scoring.func.snvs.core <- function(rdata_dir, exp.genes=NULL, 
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
    
    nec_annot_cols <- c("chrom", "pos", "ref_base", "alt_base","cchange", "hg19__chrom", "hg19__pos", "all_mappings","achange","hugo","transcript","so","samples","numsample",
                        "spliceai__ds_dl","spliceai__ds_al","vest__score","gnomad__af","fathmm_xf_coding__fathmm_xf_coding_pred","siphy__logodds_rank","cscape_coding__score","phastcons__phastcons30_mamm_r",
                        "cgc__class","cosmic__variant_count","cgl__class","cancer_hotspots__samples", "civic__clinical_a_score","clinvar__sig")
    
    if(!all(nec_annot_cols%in%colnames(OC.input))) {
      stop("Not all necessary OpenCRAVAT annotator names for creating the relative scores !!")
    } 
    
    OC.input <- OC.input %>% dplyr::select(all_of(nec_annot_cols))
    
    
  }
  
  else if (ref.genome == "hg38") {
    
    nec_annot_cols <- c("chrom", "pos", "ref_base", "alt_base","cchange", 
                        "all_mappings","achange","hugo","transcript","so","samples","numsample","spliceai__ds_dl",
                        "spliceai__ds_al","vest__score","gnomad__af","fathmm_xf_coding__fathmm_xf_coding_pred","siphy__logodds_rank", "cscape_coding__score","phastcons__phastcons30_mamm_r","cgc__class",
                        "cosmic__variant_count","cgl__class","cancer_hotspots__samples", "civic__clinical_a_score","clinvar__sig")
    
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
  
  # modify the output of the cancer hotspots annotation db in order to sum the total variant count number
  
  ch_dt <- d2 %>%
    dplyr::select(cancer_hotspots__samples) %>%
    rownames_to_column() %>%
    group_by(rowname ) %>%
    mutate(cancer__hotspots__vc = ifelse(cancer_hotspots__samples == "[]",
                                         yes = "[{\"tissue\": \"none\", \"count\": 0}]",
                                         no = cancer_hotspots__samples ) %>%
             jsonlite::fromJSON(txt = .) %>%
             summarize(Count=sum(count)) %>%
             deframe())
  
  
  d2 <- d2 %>% mutate(Cancer_hotspots_vcount=ch_dt$cancer__hotspots__vc)
  
  # Initial creation of pathogenicity/effect prediction score
  
  d3 <- d2 %>% filter(so %in% c("start_lost","stop_lost", "stop_gained",
                                "missense_variant","splice_site_variant")) %>%
    
    dplyr::filter((!gnomad__af > 0.05) %>% replace_na(TRUE)) %>% 
    
    mutate(Pathogenicity_score= case_when(
      so=="splice_site_variant" &  (spliceai__ds_dl > 0.5 | spliceai__ds_al > 0.5) ~1,
      TRUE~0),
      
      Pathogenicity2_score=case_when(
        
        so %in% c("missense_variant","stop_gained","start_lost","stop_lost") & vest__score>0.6 & fathmm_xf_coding__fathmm_xf_coding_pred=="Damaging" & cscape_coding__score>0.6 ~1,
        
        so %in% c("missense_variant","stop_gained","start_lost","stop_lost") & vest__score >0.6 & fathmm_xf_coding__fathmm_xf_coding_pred=="Damaging" ~0.8,
        
        so %in% c("missense_variant","stop_gained","start_lost","stop_lost") & vest__score >0.6 & cscape_coding__score>0.6 ~0.8,
        
        so %in% c("missense_variant","stop_gained","start_lost","stop_lost") &fathmm_xf_coding__fathmm_xf_coding_pred=="Damaging" & cscape_coding__score>0.6 ~0.8,
        
        so %in% c("missense_variant","stop_gained","start_lost","stop_lost") & 
          vest__score>0.6 | fathmm_xf_coding__fathmm_xf_coding_pred=="Damaging" | cscape_coding__score>0.6 ~0.5,
        
        TRUE~0),
      
      evol_conserv_score= case_when(
        
        siphy__logodds_rank > 0.6 & phastcons__phastcons30_mamm_r > 0.6 ~1,
        
        siphy__logodds_rank > 0.6 | phastcons__phastcons30_mamm_r > 0.6 ~0.5,
        
        TRUE~0),
      
    ) %>% mutate(Merged_Pathogenicity_score = rowSums(select(., Pathogenicity_score,Pathogenicity2_score))) %>% 
    
    mutate(Total_predict_score=rowMeans(select(., Merged_Pathogenicity_score,evol_conserv_score)))  
  
  
  
  # Next we create the cancer evidence score
  
  d4 <- d3 %>%
    
    mutate(
      cancer_gene_ev_score=case_when(
        str_detect(cgc__class,"Oncogene|fusion|TSG") & cgl__class%in%c("Oncogene","TSG")~1,
        
        str_detect(cgc__class,"Oncogene|fusion|TSG") | cgl__class%in%c("Oncogene","TSG")~0.5,
        
        TRUE~0),
      
      cancer_freq_score= case_when(
        
        cosmic__variant_count > 10 & Cancer_hotspots_vcount > 10 ~1,
        
        cosmic__variant_count > 10 | Cancer_hotspots_vcount > 10 ~0.5,
        TRUE~0),
      
      clin_ev_score= case_when(
        
        civic__clinical_a_score > 1 & str_detect(clinvar__sig, "Pathogenic|drug response|pathogenic$")~1,
        
        civic__clinical_a_score > 1 | str_detect(clinvar__sig, "Pathogenic|drug response|pathogenic$")~0.5,
        TRUE~0)) %>% mutate(Total_Cancer_score = rowMeans(select(., cancer_gene_ev_score,
                                                                 cancer_freq_score)))
  
  # Final merging of the created sub-scores-follow this strategy in order not to have issues when the user does not provide any character vector of expressed genes
  
  if (is.null(exp.genes)){
    
    
    d.final <- d4 %>% mutate(Final_Score= (Total_predict_score*w1 + Total_Cancer_score*w2 +    clin_ev_score*w3)/2.4)
    
  } else if (!is.null(exp.genes)){
    
    d.final <- d4 %>% mutate(Final_Score= (Total_predict_score*w1 + Total_Cancer_score*w2 +
                                             clin_ev_score*w3 + Exp_score*w4)/3.4)
    
  } 
  
  
  write_csv(d.final,file=str_c("Finaloutput.Variant.Ranked.OpenCRAVAT",
                               sample.name.output,"csv",sep="."))
  
  dat.sel.top10 <- d.final %>% arrange(desc(Final_Score)) %>% top_n(10)
  
  mytable <- datatable(dat.sel.top10, filter = 'top', options = list(paging = FALSE))
  
  saveWidget(mytable, file=str_c("Top10RankedVariants",
                                 sample.name.output,"html",sep="."))
  
}