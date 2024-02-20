## download reference genome and load required packages
genome_ref_dir = "/n/scratch3/users/m/mx40/working/harmonized_atlas_revision/GWAS_revision/reference_genome"
if(!file.exists(sprintf("%s/g1000_eur.bed",genome_ref_dir))){
  download.file("https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",destfile=sprintf("%s.zip",genome_ref_dir))
  unzip(sprintf("%s.zip",genome_ref_dir),exdir=genome_ref_dir)
}
genome_ref_path = sprintf("%s/g1000_eur",genome_ref_dir)

library(MAGMA.Celltyping)
library(EWCE)
library(One2One)
library(Seurat)
library(tidyverse)
library(rlist)
library(R.utils)
library(foreach)
library(doParallel)
library(tm)

numCores <- detectCores()
registerDoParallel(numCores) 


setwd("~")
## load quantile groups for celltypes
ctd_mouse_subtype <- readRDS("gene_quantile_40.Rds")

## headache gwas summary statistics file path
gwas_sumstats_path_headache = "3799.gwas.imputed_v3.both_sexes.tsv"
#gunzip(sprintf("%s.bgz",gwas_sumstats_path_headache),gwas_sumstats_path_headache)


## control gwas summary statistics file path
gwas_sumstats_path_control = "6159_100.gwas.imputed_v3.both_sexes.tsv"
#gunzip(sprintf("%s.bgz",gwas_sumstats_path_control),gwas_sumstats_path_control)

gwas_list <- c(
               gwas_sumstats_path_headache,gwas_sumstats_path_control
               
)

gwas_name_list <- c(
                    "headache",
                    "control")


gwas_output <- function(ctd_file,gwas_sumstats_path,genome_ref_path,gwas_name,counts_name,species,cluster,gwas_sumstats_path_control,temp_upstream_kb,temp_downstream_kb){
  
  
  ## check whether gwas is properly formatted
  gwas_sumstats_path_formatted = sprintf("%s.formatted",gwas_sumstats_path)
  if(!file.exists(gwas_sumstats_path_formatted)){
    col_headers = format_sumstats_for_magma(gwas_sumstats_path,1)
    file.copy(from=col_headers,to=gwas_sumstats_path_formatted,overwrite = TRUE)
  }
  
  ## Map SNPs to Genes
  genesOutPath = map.snps.to.genes(gwas_sumstats_path_formatted,genome_ref_path=genome_ref_path,genome_build = "GRCh37",upstream_kb = temp_upstream_kb,downstream_kb = temp_downstream_kb)
  
  gwas_sumstats_path_control_formatted = sprintf("%s.formatted",gwas_sumstats_path_control)
  if(!file.exists(gwas_sumstats_path_control_formatted)){
    col_headers_control = format_sumstats_for_magma(gwas_sumstats_path_control,1)
    file.copy(from=col_headers_control,to=gwas_sumstats_path_control_formatted,overwrite = TRUE)
  }
  genesOutPath_control = map.snps.to.genes(gwas_sumstats_path_control_formatted,genome_ref_path=genome_ref_path,genome_build = "GRCh37",upstream_kb = temp_upstream_kb,downstream_kb = temp_downstream_kb)
  
  
  
  l5_path <- paste0(gwas_name,"_",counts_name,"_",species,"_",cluster,"_upstream",temp_upstream_kb,"_downstream",temp_downstream_kb)
  dir.create(l5_path)
  setwd(l5_path)
  
  ## test GWAS enrichments that remain in a GWAS after controling for a second GWAS
  controlGenesOut = sprintf("%s.genes.out",get.magma.paths(gwas_sumstats_path_control_formatted,upstream_kb = temp_upstream_kb,downstream_kb = temp_downstream_kb)$filePathPrefix)
  
  ctAssocsLinear_control_for_second_gwas = calculate_celltype_associations(ctd_file,gwas_sumstats_path_formatted,
                                                                           genome_ref_path=genome_ref_path,specificity_species = "mouse",
                                                                           genesOutCOND=controlGenesOut,analysis_name = "Controlling",
                                                                           upstream_kb = temp_upstream_kb,downstream_kb = temp_downstream_kb)
  
  ## save GWAS enrichment results
  list.save(ctAssocsLinear_control_for_second_gwas,"ctAssocsLinear_control_for_second_gwas.rds")

  subtype <- ctAssocsLinear_control_for_second_gwas[[1]]$results
  
  write.csv(subtype,paste0(cluster,".csv"))
  
  setwd("../")
}




foreach (up_m = c(1,10),.combine ='c') %dopar% {
  foreach (down_n = c(1,10),.combine ='c') %dopar% {
    foreach(i = 1:length(gwas_list),.combine ='c') %dopar% {
      gwas_output(ctd_mouse_subtype,gwas_list[i],genome_ref_path,gwas_name_list[i],"mouse_TG_neuron_nonneuron","mouse","subtype",gwas_sumstats_path_control,up_m,down_n)
    }
  }
}
