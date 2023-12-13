# produce the input for analysis ( alph-diversity beta-diversity etc)

# 2020/10/17 

library(tidyverse)
library(fs)
library(data.table)
library(metagenomeSeq)
library(sva)
#BiocManager::install("sva")

#library(readr)
#library(vegan)



# list files 

workdir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/OTU_clustering_gg13_8/output_from_qiime2/"

workdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/OTU_clustering_gg13_8/output_from_qiime2/"

outdir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/OTU_clustering_gg13_8/output/Normalization_part/"

outdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/OTU_clustering_gg13_8/output/Normalization_part/"

files <- dir(workdir, pattern = "*taxonomy.tsv") # i did some manual fix on the file : first delete the first line; second separate the taxonomy column

table <- files


# meta 
metatable <- read_delim("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/20201007-metatable-simple.xlsx",delim = "\t") %>% 
  .[2:nrow(.),]
metatable <- read_delim("/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/20201007-metatable-simple.xlsx",delim = "\t") %>% 
  .[2:nrow(.),]




read_table_fun <- function(workdir, table){
  otu <- read.table(paste(workdir, table, sep = "/"), stringsAsFactors = F, sep = "\t", header = T)
  colnames(otu) <- str_remove_all(colnames(otu),"X")
  taxonomy <- str_split(otu$taxonomy,"; ",simplify = T)
  colnames(taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  otu <- otu %>%
    select(-taxonomy) %>%
    cbind(., taxonomy, stringsAsFactors = F)
  return(otu)
}


Normalization_metagenome <- function(OTUfile, col, outdir){
  normfile <- OTUfile[,2:col] %>%    
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    metagenomeSeq::newMRexperiment(.) %>%       ## creating MR object                      
    metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
    metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>%
    as_tibble(.) %>%
    add_column(OTU = OTUfile[,1], Kingdom = OTUfile$Kingdom, Phylum = OTUfile$Phylum, 
               Class = OTUfile$Class, Order = OTUfile$Order, 
               Family = OTUfile$Family, Genus = OTUfile$Genus,
               Species = OTUfile$Species) %>%
    as.data.table() 
  write_csv(normfile, paste(outdir, "normalized_otu_table.csv", sep = ""))
  return(normfile)
}



ten_percent_prevelance <- function(otu_norm, col){
  otu_norm_sel <- list()
  for (i in (col-7):col){
    print(i)
    otu_norm_sel[[colnames(otu_norm)[i]]] <- otu_norm %>%
      select(-OTU, -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
      cbind(., otu_norm[,i, with = FALSE])
    colnames(otu_norm_sel[[colnames(otu_norm)[i]]])[col-7] <- "tax"
    otu_norm_sel[[colnames(otu_norm)[i]]] <- otu_norm_sel[[colnames(otu_norm)[i]]] %>%
      pivot_longer(., -tax, names_to = "sample") %>%
      group_by(tax, sample) %>%
      summarise(., tax_sum = sum(value)) %>%  # check it alreay,correct : > dplyr::filter(otu_norm, Phylum ==" p__Bacteroidetes") %>% .[,1] %>% sum() [1] 1127.818 > dplyr::filter(otu_norm, Phylum ==" p__Bacteroidetes") %>% .[,2] %>% sum() [1] 554.3
      #filter(., !tax =="") %>%
      pivot_wider(., names_from = sample, values_from = tax_sum) %>% #values_fill = list(values = 0) 
      column_to_rownames(., var = "tax") %>%
        .[rowSums(. > 0) >= 0.1*158,]
  }
  return(otu_norm_sel)
 
}

Write_normalized_table <- function(outdir, otu_norm_sel, name_of_file){
  for (i in names(otu_norm_sel)){
    print(i)
    otu_norm_sel[[i]] <- otu_norm_sel[[i]] %>%
      rownames_to_column(., var = i) 
    write_csv(otu_norm_sel[[i]], paste(outdir, paste(i, name_of_file, sep = "_"), sep = ""))
  }
}



Normalization_pip <- function(workdir, table){
  otu <- read_table_fun(workdir, table)
  col <- ncol(otu)
  otu_norm <- Normalization_metagenome(otu, col, outdir)
  otu_norm_sel <- list()
  otu_norm_sel <- ten_percent_prevelance(otu_norm, col)
  Write_normalized_table(outdir, otu_norm_sel, "10_percent_normalized.csv")
  return(otu_norm_sel)
}





make_metatable_sel <- function(group_file, otu_norm_sel){
  metatable_sel <- dplyr::filter(group_file, SampleID %in% colnames(otu_norm_sel[[1]])) %>%
    add_column(., Group = ifelse(.$Class %in% c("1","2","3","4","5"), "H","C")) %>%
    mutate(SampleID = factor(.$SampleID, levels = colnames(otu_norm_sel[[1]]))) %>%
    add_column(Batch = str_sub(.$SampleID, 0,7)) %>%
    add_column(sub_class1 = ifelse(.$Class %in% c("1","2","3","4","5"), "H", .$Class)) %>%
    add_column(sub_class2 = ifelse(.$Class %in% c("9","10","12","13"), "C", .$Class)) %>%
    as.data.table()
  metatable_sel <- metatable_sel[order(SampleID)]
  return(metatable_sel)
}



combat_remove <- function(input_matrix, batch_file){
  cmb.model <- model.matrix(~1, data = batch_file)
  input_matrix_corr <- ComBat(dat = input_matrix,
                              batch = batch_file$batch,
                              mod = cmb.model,
                              par.prior = TRUE,
                              prior.plots = FALSE)
  return(input_matrix_corr)
}


log_transform_combat_exp <- function(input_table, batch.mat){
  input_matrix <- as.matrix(input_table)
  input_matrix[input_matrix==0] <- min(input_matrix[input_matrix>0])/2
  input_matrix_log <- as.matrix(log(input_matrix))
  input_matrix_corr <- combat_remove(input_matrix_log, batch.mat) 
  input_matrix_corr_exp <- as.data.frame(exp(input_matrix_corr))
  return(input_matrix_corr_exp)
}



Combat_normalization_pip <- function(otu_norm_sel, metatable_sel, outdir){ # combat from sva 
  # first combat remove 
  batch.mat <- cbind.data.frame("sample" = metatable_sel$SampleID, 
                                "batch" = metatable_sel$Batch)
  print(head(batch.mat))
  
  # second combat remove 
  batch.mat2 <- cbind.data.frame("sample" = metatable_sel$SampleID, 
                                "batch" = metatable_sel$sub_class1)
  # third combat remove 
  batch.mat3 <- cbind.data.frame("sample" = metatable_sel$SampleID, 
                                 "batch" = metatable_sel$sub_class2)
 
  # build list 
  otu_norm_sel_batch1 <- list()
  # loop - 1st 
  dir_create(paste(outdir, 'batch1', sep = "")) # from fs package
  
  for (i in names(otu_norm_sel)){
    print(i)
    if (nrow(otu_norm_sel[[i]])>1){
      otu_norm_sel_batch1[[i]] <- otu_norm_sel[[i]] %>%
        log_transform_combat_exp(., batch.mat) 
    }else{
      print(paste("only 1 row:",otu_norm_sel[[i]],sep = ":"))
    }
    
  }
  Write_normalized_table(outdir = paste(outdir, 'batch1/', sep = ""),otu_norm_sel = otu_norm_sel_batch1, name_of_file = "10_percent_normalized_batch1_combat.csv")
  
  # loop - 2nd and 3st : together 
  otu_norm_sel_batch2 <- list()
  for (i in names(otu_norm_sel)){
    print(i)
    dir_create(paste(outdir, 'batch2', sep = "")) # from fs package
    if (nrow(otu_norm_sel[[i]]) >1){
      otu_norm_sel_batch2[[i]] <- otu_norm_sel_batch1[[i]] %>%
        log_transform_combat_exp(., batch.mat2) %>%
        log_transform_combat_exp(., batch.mat3)
    }else{
      print(paste("only 1 row:",otu_norm_sel[[i]],sep = ":"))
    }
   
  }
  Write_normalized_table(outdir = paste(outdir, 'batch2/', sep = ""),otu_norm_sel = otu_norm_sel_batch2, name_of_file = "10_percent_normalized_batch2_combat.csv")
  
}


# main_command
#=======================================================================================================================
otu_norm_sel <- list()
otu_norm_sel <- Normalization_pip(workdir, table)
metatable_sel <- make_metatable_sel(metatable, otu_norm_sel)
write_csv(metatable_sel, paste(outdir, "metatable_sel.csv", sep = "/"))
Combat_normalization_pip(otu_norm_sel,metatable_sel, outdir)


# another file 
workdir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_from_qiime2/"

workdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_from_qiime2/"

outdir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output/normalization_part/"

outdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output/normalization_part/"

files <- dir(workdir, pattern = "*taxonomy.tsv") # i did some manual fix on the file : first delete the first line; second separate the taxonomy column

# meta 
metatable <- read_delim("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/20201007-metatable-simple-deblur.xlsx",delim = "\t") %>% 
  .[2:nrow(.),]
metatable <- read_delim("/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/20201007-metatable-simple-deblur.xlsx",delim = "\t") %>% 
  .[2:nrow(.),]


table2 <- files


otu_norm_sel2 <- list()
otu_norm_sel2 <- Normalization_pip(workdir, table2)
metatable_sel2 <- make_metatable_sel(metatable, otu_norm_sel2)
write_csv(metatable_sel2, paste(outdir, "metatable_sel.csv", sep = "/"))
Combat_normalization_pip(otu_norm_sel2,metatable_sel2, outdir)



