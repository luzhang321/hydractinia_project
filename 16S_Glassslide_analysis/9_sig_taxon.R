# sig.taxons of glassslide samples 

#######################################################################################################


#                                   load library                                                      # 


#######################################################################################################

library(metagenomeSeq)
library(tidyverse)
library(data.table)
library(ggpubr)
library(rstatix)
library(Rfit)
library(viridis)
library(hrbrthemes)
library(VennDiagram)
library(metagenomeSeq)

# this is run on local linux mint computer 
# 1. first get the whole count form the separte new matrix 
# 2. normalization & do sig OTU level & draw 3 heatmaps
# 3. normalization & do significant genus/families/orders 
# 3. draw heatmap separately for GS vs GS+HS  
# 4. draw heatmap separately for GS vs GS+HS*
# 5. draw heatmap separately for GS+HS vs GS+HS* 
#######################################################################################################


#                                   Prepare input                                                      # 


#######################################################################################################



setwd("/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/9_sig_taxon/")
files_original_counts <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/7_alphadiversity/glasslide_average_count_table.tsv" %>%
  read_tsv()

input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/7_alphadiversity/"
load(paste0(input_dir, "glassslide_alpha_diversity.RData"))

metadata <- alpha_table_bee

# make taxonomy table 

otu <- files_original_counts
colnames(otu)[1] <- "OTUID"

otu_tax <- select(otu, OTUID, Kingdom, Phylum, Class,
                  Order, Family, Genus, Species)


# form the format for normalization

otu <- otu %>%
  select(., -Kingdom, -Phylum, -Class,-Order, -Family, -Genus, -Species)


# form separate matrix 
# separate otu table 
#=======================================================================
input <- otu

select_target_class <- function(input = NA, key_word = NA, metatable_sel_biorep = NA){
  
  otu_separate <- input %>%
    pivot_longer(., -OTUID, names_to = "sample_id_merged", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("sample_id_merged","Group","Class")], by = "sample_id_merged") %>%
    filter(., Class %in% key_word) %>% 
    column_to_rownames(., var = "sample_id_merged") %>% 
    select(., -Class, -Group)
  
  return(otu_separate)
  
}
otu_GS <- select_target_class(input = input, key_word = c(16), metatable_sel_biorep = metadata)  
otu_GS_HS <- select_target_class(input = input, key_word = c(15), metatable_sel_biorep = metadata) 
otu_GS_HS_dead <- select_target_class(input = input, key_word = c(20), metatable_sel_biorep = metadata) 


#######################################################################################################


#                                  Normalization & sig.otu                                             # 


#######################################################################################################


# metagenomseq function 

metagenomeseq_difference_abundance <- function(Dataset_1 = NA,Dataset_2 = NA, prevCutoff = 0.1, Combined_metadata = NA) {
  
  
  dataset_name_1 <- deparse(substitute(Dataset_1)) # get the name 
  dataset_name_2 <- deparse(substitute(Dataset_2))
  # 
  ### combine them together 
  
  full_features = rbind(Dataset_1, Dataset_2)
  
  ## setting the metadata samples that exist in the current comparison
  Combined_metadata <- Combined_metadata %>%
    column_to_rownames(., var = "sample_id_merged")
  full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]
  
  
  # full_metadata$Groups = factor(c(rep(dataset_name_1,nrow(Dataset_1)),
  #                                 rep(dataset_name_2,nrow(Dataset_2))))
  
  ### Order them : same order 
  full_metadata = full_metadata[order(rownames(full_metadata)),] #%>% .[complete.cases(.), ]
  full_features = full_features[order(rownames(full_features)),] %>% .[which(rownames(.) %in% rownames(full_metadata)),]
  
  
  ### transoposed data frame of features
  tran = t(full_features) %>% as.data.frame(.)
  
  ### Taxonomy dataframe , item needed for MRExperiment object
  taxaData = data.frame(colnames(full_features))
  rownames(taxaData) = names(full_features)
  colnames(taxaData) = "OTU"
  taxaData <- AnnotatedDataFrame(taxaData)
  
  ### Metadata in annotated daa frame format with teh same order and rownames as the full features
  phenoData <- AnnotatedDataFrame(full_metadata)
  
  ##3 Creating MR experiment Object
  prof_mr <- newMRexperiment(counts = tran, phenoData = phenoData, featureData = taxaData)
  
  
  prevalence = ncol(prof_mr)*0.1
  ## Filtering for 1) Depth > 0.001 2) 10% prevalence
  #prof_mr <- filterData(prof_mr, depth = 0.001, present = prevalence)
  prof_mr <- filterData(prof_mr, present = prevalence)
  
  ## Cumulative Sum Scaling Normalisation
  p <- cumNormStatFast(prof_mr)
  tmpMR <- cumNorm(prof_mr, p)
  
  #NORMFACTOR <- normFactors(prof_mr) # lU ADDING 
  ## Creating the model 
  mod <-  model.matrix(~ Group, data = pData(tmpMR))
  
  
  ## FitZig =  zero inflated model ideal for metaphlan data
  results_zeller <-  fitZig(tmpMR, mod)
  
  ## Results
  result = MRtable(results_zeller, number = nrow(results_zeller))
  result = result[result$adjPvalues < 0.05,]
  
  return(result)
}  


diffotus_GS_GS_HS <- metagenomeseq_difference_abundance(Dataset_1 = otu_GS, Dataset_2 = otu_GS_HS, 
                                                          Combined_metadata = metadata) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = .$`GroupGS+HS`) # the default is GS+HS higher will be postivie 

diffotus_GS_GS_HS_dead <- metagenomeseq_difference_abundance(Dataset_1 = otu_GS, Dataset_2 = otu_GS_HS_dead, 
                                                                Combined_metadata = metadata) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = .$`GroupGS+HS*`) # the default is GS+HS+dead higher will be postivie 

diffotus_GS_HS_GS_HS_dead <- metagenomeseq_difference_abundance(Dataset_1 = otu_GS_HS, Dataset_2 = otu_GS_HS_dead, 
                                                             Combined_metadata = metadata) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = .$`GroupGS+HS*`) # the default is GS+HS+dead higher will be postivie 


#######################################################################################################


#                         Draw heatmap for sig.OTUs                                                   # 


#######################################################################################################


# draw 2 heatmaps for the OTU levels 
otudir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/9_sig_taxon/"

# normalization object based on the selected group 


Normalization_metagenome_selsamples <- function(OTUfile = NA, otudir = NA, removegroups = NA, combine_metadata = NA, name = NA){
  
  keep_samples <- combine_metadata[combine_metadata$Group != removegroups,]$sample_id_merged
  colnames(OTUfile)[1] <- "OTU"
  OTUfile <- OTUfile %>%
    select(., c("OTU", all_of(keep_samples), "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  col <- ncol(OTUfile)
  
  normfile <- OTUfile %>% 
    column_to_rownames(., var = "OTU") %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    metagenomeSeq::newMRexperiment(.) %>%       ## creating MR object                      
    metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
    metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>%
    as_tibble(.) %>%
    add_column(OTU = OTUfile$OTU, Kingdom = OTUfile$Kingdom, Phylum = OTUfile$Phylum,  # the order will not change tibble(data1 = rownames(normfile), data2 = OTUfile$OTU) %>% filter(., data1!=data2) : 0 
               Class = OTUfile$Class, Order = OTUfile$Order, 
               Family = OTUfile$Family, Genus = OTUfile$Genus,
               Species = OTUfile$Species) %>%
    as_tibble() 
  write_csv(normfile, paste(otudir, name, ".csv", sep = ""))
  return(normfile)
}


combine_metadata <- metadata

renormalized_otu_GS_GS_HS <- Normalization_metagenome_selsamples(OTUfile = files_original_counts, otudir = otudir, removegroups = "GS+HS*",
                                                              combine_metadata = metadata, name = "GS_GS_HS_renormalization_table")


renormalized_otu_GS_GS_HS_dead <- Normalization_metagenome_selsamples(OTUfile = files_original_counts, otudir = otudir, removegroups = "GS+HS",
                                                                   combine_metadata = metadata, name = "GS_GS_HS_dead_renormalization_table")

renormalized_otu_GS_HS_GS_HS_dead <- Normalization_metagenome_selsamples(OTUfile = files_original_counts, otudir = otudir, removegroups = "GS",
                                                                      combine_metadata = metadata, name = "GS_HS_GS_HS_dead_renormalization_table")






# heatmap drawing function 

draw_pheatmap <- function(renormalized_table = NA, combine_metadata = NA, differentOTU = NA, tag = NA){
  
  
  input_sel_outtable <- renormalized_table %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species, -OTU)
  
  
  group_info <- data.table(sample_id_merged = colnames(input_sel_outtable)) %>%
    left_join(., combine_metadata, by = "sample_id_merged") %>%
    select(., sample_id_merged, Group) %>%
    column_to_rownames(., var = 'sample_id_merged') %>% 
    arrange(., Group)
  
  
  input_pheatmap <- renormalized_table %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    column_to_rownames(., var = "OTU")
  
  if (tag != "nofilter"){
    input_pheatmap_sel <- input_pheatmap %>%
      .[rowSums(.) > 0,] %>%
      .[differentOTU$OTU,] %>%
      .[,rownames(group_info)]
  }else if (tag == "nofilter"){
    input_pheatmap_sel <- input_pheatmap %>%
      .[rowSums(.) > 0,] %>%
      .[,rownames(group_info)]
  }
  
  
  
  # draw heatmap 
  pheatmap_object <- pheatmap::pheatmap(log(input_pheatmap_sel+0.0001), 
                                        annotation_col = group_info,
                                        annotation_colors = ann_colors,
                                        #annotation_row = input_taxnomy,
                                        show_rownames = FALSE,
                                        cluster_cols = F,
                                        fontsize_col = 4)
  
  
  return(pheatmap_object)
  
  
}


combine_metadata <- metadata
color <- viridis::viridis(10)[c(4,6,10)]

ann_colors = list(
  Group = c("GS" = viridis::viridis(10)[c(4)], "GS+HS" = viridis::viridis(10)[c(6)], "GS+HS*" = viridis::viridis(10)[c(10)])
)



GS_GS_HS_OTU_pheatmap <- draw_pheatmap(renormalized_table = renormalized_otu_GS_GS_HS, combine_metadata = metadata, 
                                    differentOTU = diffotus_GS_GS_HS, tag = "filter")

GS_GS_HS_dead_OTU_pheatmap <- draw_pheatmap(renormalized_table = renormalized_otu_GS_GS_HS_dead, combine_metadata = metadata, 
                                       differentOTU = diffotus_GS_GS_HS_dead, tag = "filter")

GS_HS_GS_HS_dead_OTU_pheatmap <- draw_pheatmap(renormalized_table = renormalized_otu_GS_HS_GS_HS_dead, combine_metadata = metadata, 
                                       differentOTU = diffotus_GS_HS_GS_HS_dead, tag = "filter")


pdf(paste(otudir, "GS_GS_HS_OTU_pheatmap.pdf", sep = "/"),width = 10, height = 12)
GS_GS_HS_OTU_pheatmap
dev.off()

pdf(paste(otudir, "GS_GS_HS_dead_OTU_pheatmap.odf", sep = "/"),width = 10, height = 12)
GS_GS_HS_dead_OTU_pheatmap
dev.off()

pdf(paste(otudir, "GS_HS_GS_HS_dead_OTU_pheatmap.pdf", sep = "/"),width = 10, height = 12)
GS_HS_GS_HS_dead_OTU_pheatmap
dev.off()

#######################################################################################################


#  prepare the input & normalization for phylum, families, genus level separately                     # 


#######################################################################################################


# produce different levels table 
#====================
produce_different_taxa <- function(otu_norm = NA, col = NA){
  otu_norm_sel <- list()
  for (i in (col-7):col){
    print(i)
    otu_norm_sel[[colnames(otu_norm)[i]]] <- otu_norm %>%
      select(-OTU, -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
      cbind(., otu_norm[,i])
    colnames(otu_norm_sel[[colnames(otu_norm)[i]]])[col-7] <- "tax"
    otu_norm_sel[[colnames(otu_norm)[i]]] <- otu_norm_sel[[colnames(otu_norm)[i]]] %>%
      pivot_longer(., -tax, names_to = "sample") %>%
      group_by(tax, sample) %>%
      summarise(., tax_sum = sum(value)) %>%  # check it alreay,correct : > dplyr::filter(otu_norm, Phylum ==" p__Bacteroidetes") %>% .[,1] %>% sum() [1] 1127.818 > dplyr::filter(otu_norm, Phylum ==" p__Bacteroidetes") %>% .[,2] %>% sum() [1] 554.3
      #filter(., !tax =="") %>%
      pivot_wider(., names_from = sample, values_from = tax_sum) #values_fill = list(values = 0) 
    # for the ones with corresponding tax name is NA 
    otu_norm_sel[[colnames(otu_norm)[i]]]$tax[is.na(otu_norm_sel[[colnames(otu_norm)[i]]]$tax)] <- "NA"
    
    otu_norm_sel[[colnames(otu_norm)[i]]] <- otu_norm_sel[[colnames(otu_norm)[i]]] %>%
      column_to_rownames(., var = "tax") 
  }
  return(otu_norm_sel)
  
}



otu_combined_count <- list()
temp <- files_original_counts %>%
  add_column(OTU = .$OTUID) %>%
  select(.,-OTUID)

col <- ncol(temp)

otu_combined_count <- produce_different_taxa(otu_norm = temp, col = col)

# separate otu table depending on samples 
#====================


separate_subgroups <- function(table = NA, metatable_sel_biorep = NA, key_word = NA){
  
  colnames(metatable_sel_biorep)[4] <- "SampleID"
  
  table_sel <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% key_word) %>% 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  return(table_sel)
}




# separate main program 
otu_combined_count_sub <- list()
for (i in names(otu_combined_count)[names(otu_combined_count)!="Kingdom"]){
  print(i)
  otu_combined_count_sub[[i]][["GS"]] <- separate_subgroups(otu_combined_count[[i]],
                                                                    metatable_sel_biorep = metadata, key_word = c(16))
  otu_combined_count_sub[[i]][["GS_HS"]] <- separate_subgroups(otu_combined_count[[i]],
                                                            metatable_sel_biorep = metadata, key_word = c(15))
  otu_combined_count_sub[[i]][["GS_HS_dead"]] <- separate_subgroups(otu_combined_count[[i]],
                                                            metatable_sel_biorep = metadata, key_word = c(20))
  
}


#######################################################################################################


#  Differentially abundant Taxon                                                                      # 


#######################################################################################################



metagenomeseq_diff_taxa <- function(input_list = NA, metadata = NA){
  
  otu_GS_taxa <- input_list[["GS"]]
  otu_GS_HS_taxa <- input_list[["GS_HS"]]
  otu_GS_HS_dead_taxa <- input_list[["GS_HS_dead"]]
  
  diffotus_GS_GS_HS_taxa <- metagenomeseq_difference_abundance(Dataset_1 = otu_GS_taxa, Dataset_2 = otu_GS_HS_taxa, 
                                                                 Combined_metadata = metadata) %>%
    mutate(.,Taxa = rownames(.)) %>% mutate(., Group_Impact = .$`GroupGS+HS`)
  
  diffotus_GS_GS_HS_dead_taxa <- metagenomeseq_difference_abundance(Dataset_1 = otu_GS_taxa, Dataset_2 = otu_GS_HS_dead_taxa, 
                                                            Combined_metadata = metadata) %>%
    mutate(.,Taxa = rownames(.)) %>% mutate(., Group_Impact = .$`GroupGS+HS*`)
  
  diffotus_GS_HS_GS_HS_dead_taxa <- metagenomeseq_difference_abundance(Dataset_1 = otu_GS_HS_taxa, Dataset_2 = otu_GS_HS_dead_taxa, 
                                                                 Combined_metadata = metadata) %>%
    mutate(.,Taxa = rownames(.)) %>% mutate(., Group_Impact = .$`GroupGS+HS*`)
  
  result <- list("GS_GS_HS" = diffotus_GS_GS_HS_taxa,
                 "GS_GS_HS_dead" = diffotus_GS_GS_HS_dead_taxa,
                 "GS_HS_GS_HS_dead" = diffotus_GS_HS_GS_HS_dead_taxa)
  
  
  return(result)
}


otu_combined_count_sub_result <- list()
otudir
for (i in names(otu_combined_count_sub)){
  print(i)
  otu_combined_count_sub_result[[i]] <- metagenomeseq_diff_taxa(input_list = otu_combined_count_sub[[i]], 
                                                                metadata = metadata)
  write_csv(otu_combined_count_sub_result[[i]][["GS_GS_HS"]], paste(otudir, paste("Taxonomy/sig", i, "GS_GS_HS.csv", sep = "_"), sep = "/"))
  write_csv(otu_combined_count_sub_result[[i]][["GS_GS_HS_dead"]], paste(otudir, paste("Taxonomy/sig", i, "GS_GS_HS_dead.csv", sep = "_"), sep = "/"))
  write_csv(otu_combined_count_sub_result[[i]][["GS_HS_GS_HS_dead"]], paste(otudir, paste("Taxonomy/sig", i, "GS_HS_GS_HS_dead.csv", sep = "_"), sep = "/"))
  
}

#######################################################################################################


#  draw heatmap   - renormalization                                                                   # 


#######################################################################################################




Normalization_metagenome <- function(OTUfile = NA, otudir = NA, name = NA, removegroups = NA, metadata = NA){
  
  keep_samples <- metadata[metadata$Group != removegroups,]$sample_id_merged
  
  OTUfile <- OTUfile %>%
    select(., c(all_of(keep_samples)))
  
  
  
  normfile <- OTUfile %>%    
    metagenomeSeq::newMRexperiment(.) %>%       ## creating MR object                      
    metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
    metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(., var = "OTUID") %>% 
    as.data.table() 
  write_csv(normfile, paste(otudir, name,"_normalized_otu_table.csv", sep = ""))
  return(normfile)
}


Phylum_famiy_genus_norm <- list()
for (i in names(otu_combined_count)[names(otu_combined_count)!="Kingdom"]){
  
  Phylum_famiy_genus_norm[[i]][["GS_GS_HS"]] <- Normalization_metagenome(OTUfile = otu_combined_count[[i]],
                                                                              otudir = paste(otudir, "Taxonomy/heatmap/", sep = "/"),
                                                                              name = paste(i, "GSvsGS_HS", sep = "_"),
                                                                              removegroups = "GS+HS*",
                                                                              metadata = metadata)
  
  Phylum_famiy_genus_norm[[i]][["GS_GS_HS_dead"]] <- Normalization_metagenome(OTUfile = otu_combined_count[[i]],
                                                                         otudir = paste(otudir, "Taxonomy/heatmap/", sep = "/"),
                                                                         name = paste(i, "GSvsGS_HS_dead", sep = "_"),
                                                                         removegroups = "GS+HS",
                                                                         metadata = metadata)
  
  Phylum_famiy_genus_norm[[i]][["GS_HS_GS_HS_dead"]] <- Normalization_metagenome(OTUfile = otu_combined_count[[i]],
                                                                              otudir = paste(otudir, "Taxonomy/heatmap/", sep = "/"),
                                                                              name = paste(i, "GS_HSvsGS_HS_dead", sep = "_"),
                                                                              removegroups = "GS",
                                                                              metadata = metadata)
  
  
}

#######################################################################################################


#  draw heatmap   - renormalization                                                                   # 


#######################################################################################################

pheatmap_output_taxa <- list()
# c("#386cb0", "#f0027f")
for (j in c("Phylum", "Family", "Genus")){ # names(otu_combined_count)[names(otu_combined_count)!="Kingdom"]
  print(j)
  
  # GS vs GS+HS 
  renormalized_tablev2 <- Phylum_famiy_genus_norm[[j]][["GS_GS_HS"]]
  combine_metadatav2 <- metadata
  differentOTUv2 <- otu_combined_count_sub_result[[j]][["GS_GS_HS"]]
  
  tag <- "filter"
  ann_colors = list(
    Group = c("GS" = viridis::viridis(10)[c(4)], "GS+HS" = viridis::viridis(10)[c(6)])
  )
  
  
  pheatmap_output_taxa[[j]][["GS_GS_HS"]] <- draw_pheatmap_taxa_update(renormalized_tablev2 = renormalized_tablev2,
                                                                              combine_metadatav2 = combine_metadatav2, 
                                                                              differentOTUv2 = differentOTUv2, tag = "filter", range_value = 3, ann_colors = ann_colors)
  
  pdf(paste(otudir, "/Taxonomy/heatmap/figure/", "GS_GS_HS", "diff.", j, ".pdf", sep = ""), width = 9, height = 10*1.6)
  print(pheatmap_output_taxa[[j]][["GS_GS_HS"]])
  dev.off()
  
  # GS vs GS+HS* 
  renormalized_tablev2 <- Phylum_famiy_genus_norm[[j]][["GS_GS_HS_dead"]]
  combine_metadatav2 <- metadata
  differentOTUv2 <- otu_combined_count_sub_result[[j]][["GS_GS_HS_dead"]]
  
  tag <- "filter"
  ann_colors2 = list(
    Group = c("GS" = viridis::viridis(10)[c(4)], "GS+HS*" = viridis::viridis(10)[c(10)])
  )
  
  
  pheatmap_output_taxa[[j]][["GS_GS_HS_dead"]] <- draw_pheatmap_taxa_update(renormalized_tablev2 = renormalized_tablev2,
                                                                              combine_metadatav2 = combine_metadatav2, 
                                                                              differentOTUv2 = differentOTUv2, tag = "filter", range_value = 3, ann_colors = ann_colors2)
  
  pdf(paste(otudir, "/Taxonomy/heatmap/figure/", "GS_GS_HS_dead", "diff.", j, ".pdf", sep = ""), width = 9, height = 10*1.6)
  print(pheatmap_output_taxa[[j]][["GS_GS_HS_dead"]])
  dev.off()
  
  
  
  
  
  # GS+HS vs GS+HS* 
  renormalized_tablev2 <- Phylum_famiy_genus_norm[[j]][["GS_HS_GS_HS_dead"]]
  combine_metadatav2 <- metadata
  differentOTUv2 <- otu_combined_count_sub_result[[j]][["GS_HS_GS_HS_dead"]]
  
  tag <- "filter"
  ann_colors3 = list(
    Group = c("GS+HS" = viridis::viridis(10)[c(6)], "GS+HS*" = viridis::viridis(10)[c(10)])
  )
  
  
  pheatmap_output_taxa[[j]][["GS_HS_GS_HS_dead"]] <- draw_pheatmap_taxa_update(renormalized_tablev2 = renormalized_tablev2,
                                                                                   combine_metadatav2 = combine_metadatav2, 
                                                                                   differentOTUv2 = differentOTUv2, tag = "filter", range_value = 3, ann_colors = ann_colors3)
  
  pdf(paste(otudir, "/Taxonomy/heatmap/figure/", "GS_HS_GS_HS_dead", "diff.", j, ".pdf", sep = ""), width = 9, height = 10*1.6)
  print(pheatmap_output_taxa[[j]][["GS_HS_GS_HS_dead"]])
  dev.off()
  
  
  
}


draw_pheatmap_taxa_update <- function(renormalized_tablev2 = NA, combine_metadatav2 = NA, differentOTUv2 = NA, tag = NA, range_value = 3, ann_colors = NA){
  
  
  input_sel_outtable <- renormalized_tablev2 %>%
    column_to_rownames(., var = "OTUID")
  
  
  group_info <- data.table(sample_id_merged = colnames(input_sel_outtable)) %>%
    left_join(., combine_metadatav2, by = "sample_id_merged") %>%
    select(., sample_id_merged, Group) %>%
    column_to_rownames(., var = 'sample_id_merged') %>% 
    arrange(., Group)
  
  input_pheatmap <- input_sel_outtable
  
  if (tag != "nofilter"){
    input_pheatmap_sel <- input_pheatmap %>%
      .[rowSums(.) > 0,] %>%
      .[differentOTUv2$Taxa,] %>%
      .[,rownames(group_info)]
  }else if (tag == "nofilter"){
    input_pheatmap_sel <- input_pheatmap %>%
      .[rowSums(.) > 0,] %>%
      .[,rownames(group_info)]
  }
  
  input_pheatmap_sel <- input_pheatmap_sel %>% filter(rownames(.) != "NA") %>% filter(rownames(.) != "f__") %>% filter(rownames(.) != "g__") %>% filter(rownames(.) != "p__")
  # change the taxonomic level to the real name, removing p__,f__, g__ et al 
  rownames(input_pheatmap_sel) <- str_remove_all(rownames(input_pheatmap_sel), "p__") %>%
    str_remove_all(., "f__") %>%
    str_remove_all(., "g__")
  
  # italic the taxonomic name 
  newnames_HS <- lapply(
    rownames(input_pheatmap_sel),
    function(x) bquote(italic(.(x))))
  
 
  
  # scale it 
  input_pheatmap_sel_scale <- t(scale(t(input_pheatmap_sel)))
  # default scale is column scale, scale is generic function whose default method centers and/or scales the columns of a numeric matrix.
  # so here we used z-score to do normalization, then replace the extreme value with the same value and draw plot 
  
  # distribution of the value 
  p2 <- as.matrix(input_pheatmap_sel_scale) %>% as.numeric() %>% data.frame(values = .) %>% ggplot(., aes(values)) + geom_density(bw = "SJ")
  
  # set extreme value to a certain value 
  input_pheatmap_sel_scale[input_pheatmap_sel_scale > range_value] <- range_value
  
  # draw heatmap 
  pheatmap_object <- pheatmap::pheatmap(input_pheatmap_sel_scale, 
                                        annotation_col = group_info,
                                        annotation_colors = ann_colors,
                                        show_rownames = TRUE,
                                        cluster_cols = FALSE,
                                        # scale = "row",
                                        fontsize_col = 5,
                                        fontsize_row = 8,
                                        cellheight=12, cellwidth = 6,
                                        labels_row = as.expression(newnames_HS))
  
  
  
  
  return(pheatmap_object)
  
  
}

save.image(paste(otudir, "/glassslide_metagenomeseq.RData", sep = "/"))
