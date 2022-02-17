scoring.func.snvs <- function(rdata_dir,exp.genes=NULL,sample.name.output,
                              w1=0.7,w2=0.8,w3=0.9,w4=1){
  
  
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
  
  if(!any(colnames(OC.input)%in%c("chrom", "pos", "ref_base", "alt_base","cchange", 
                                  "all_mappings","achange", "hugo","transcript","so","samples","numsample","spliceai__ds_dl",
                                  "vest__score","gnomad__af","fathmm_xf_coding__fathmm_xf_coding_pred","siphy__logodds_rank",   
                                  "cscape_coding__score","phastcons__phastcons30_mamm_r","cgc__class",
                                  "cosmic__variant_count","cgl__class","Cancer_hotspots_vcount", "civic__clinical_a_score","clinvar__sig"))) {
    
    stop("Not necessary OpenCRAVAT annotator names for creating the merged score !!")
  }  
  
  if (is.null(exp.genes)){
    warning(paste0("No vector of gene symbols is provided, thus expression score would be 0"))
    
    d2 <- OC.input %>% mutate(Exp_score=0)
    
  } else if (is.character(exp.genes)){
    
    d2 <- OC.input %>% mutate(Exp_score= if_else(hugo%in%c(exp.genes), 1, 0, missing = 0))
    
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
  
  d3 <- d2 %>% filter(so %in% c("start_lost","stop_lost", "stop_gained", "missense_variant",
                                "splice_site_variant")) %>%
    
    
    dplyr::select(chrom, pos, ref_base, alt_base, hugo,transcript, so, samples,
                  numsample, cchange, achange,
                  all_mappings, spliceai__ds_dl, spliceai__ds_al, hg19__chrom, hg19__pos,
                  vest__score,Exp_score,gnomad__af,fathmm_xf_coding__fathmm_xf_coding_pred,
                  siphy__logodds_rank, phastcons__phastcons30_mamm_r,cscape_coding__score,
                  cgc__class,cosmic__variant_count,cgl__class,cgl__class, Cancer_hotspots_vcount, 
                  civic__clinical_a_score,clinvar__sig) %>%
    
    
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
  
  
  
  # Next we create the cancer evidence & the clinical scores
  
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
        
        !is.na(civic__clinical_a_score) & str_detect(clinvar__sig, "Pathogenic|drug response|pathogenic$")~1,
        
        !is.na(civic__clinical_a_score) | str_detect(clinvar__sig, "Pathogenic|drug response|pathogenic$")~0.5,
        TRUE~0)) %>% mutate(Total_Cancer_score = rowMeans(select(., cancer_gene_ev_score, cancer_freq_score)))
  
  
  # Final merging of the created sub-scores
  
  weight.func <- function(a,b,c,d,w1=0.7,w2=0.8,w3=0.9,w4=1){
    weighted.score <- (a*w1 + b*w2 + c*w3 + d*w4)/3.4
    return(weighted.score)
  }
  
  
  d.final <- d4 %>% mutate(Final_Score= weight.func(a=Total_predict_score, b=Total_Cancer_score, c=clin_ev_score, d=Exp_score))
  
  
  write_csv(d.final,file=str_c("Finaloutput.Variant.Ranked.OpenCRAVAT",sample.name.output,"csv", sep="."))
  
  dat.sel.top10 <- d.final %>% arrange(desc(Final_Score)) %>% top_n(10)
  
  mytable <- datatable(dat.sel.top10, filter = 'top', options = list(paging = FALSE))
  
  saveWidget(mytable, file=str_c("Top10.Ranked.OC", sample.name.output, "mytable.html",sep="."))
  
}
