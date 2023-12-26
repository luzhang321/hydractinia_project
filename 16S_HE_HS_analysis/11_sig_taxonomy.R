# sig taxonomy 
# The aim of this script is to get significant species 

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
# 2. normalization & do sig OTU level & draw 2 heatmaps
# 3. normalization & do significant genus/families/orders 
# 3. draw heatmap separately for HS vs Control 
# 4. draw heatmap separately for HE vs HS 


# part1 :prepare the input & normalization 
#=================================================================================================================================

setwd("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/")
files_original_counts <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/average_count_table.csv" %>%
  read_csv()

metadata <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/metatable_biorep_unique.csv")

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
  
  otu_separate <- otu %>%
    pivot_longer(., -OTUID, names_to = "uniquebiorep", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("uniquebiorep","Class")], by = "uniquebiorep") %>%
    filter(., Class %in% key_word) %>% 
    column_to_rownames(., var = "uniquebiorep") %>% 
    select(., -Class)
  
  return(otu_separate)
  
}
otu_HScontrol <- select_target_class(input = otu_full_table, key_word = c(32), metatable_sel_biorep = metadata) #7 samples 
otu_HScontrolv2 <- select_target_class(input = otu, key_word = c(32), metatable_sel_biorep = metadata) #7 samples 
all.equal(otu_HScontrol, otu_HScontrolv2)
#[1] TRUE
otu_HS <- select_target_class(input = otu_full_table, key_word = c(23,33), metatable_sel_biorep = metadata) #15 samples 
otu_HSv2 <- select_target_class(input = otu_full_table, key_word = c(23,33), metatable_sel_biorep = metadata) #15 samples 
all.equal(otu_HS, otu_HSv2)
#[1] TRUE
otu_HE <- select_target_class(input = otu_full_table, key_word = c(1,2,3,4), metatable_sel_biorep = metadata) #26 samples 

# above function doesn't use input parameters 

# part2 :normalization & sig.taxnomy 
#=================================================================================================================================
# metagenomseq function 

metagenomeseq_difference_abundance <- function(Dataset_1,Dataset_2,prevCutoff = 0.1,Combined_metadata) {
  
  
  dataset_name_1 <- deparse(substitute(Dataset_1)) # get the name 
  dataset_name_2 <- deparse(substitute(Dataset_2))
  # 
  ### combine them together 
  
  full_features = rbind(Dataset_1, Dataset_2)
  
  ## setting the metadata samples that exist in the current comparison
  Combined_metadata <- Combined_metadata %>%
    column_to_rownames(., var = "uniquebiorep")
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


diffotus_HS_control <- metagenomeseq_difference_abundance(Dataset_1 = otu_HS, Dataset_2 = otu_HScontrol, 
                                                         Combined_metadata = metadata) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupHS_Shell)

diffotus_HE_HS <- metagenomeseq_difference_abundance(Dataset_1 = otu_HE, Dataset_2 = otu_HS, Combined_metadata = metadata) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupHS)

# draw 2 heatmaps for the OTU levels 
otudir <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/"

# normalization object based on the selected group 


Normalization_metagenome_selsamples <- function(OTUfile = NA, otudir = NA, removegroups = NA, combine_metadata = NA, name = NA){
  
  keep_samples <- metadata[metadata$Group != removegroups,]$uniquebiorep
  
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

renormalized_otu_HE_HS <- Normalization_metagenome_selsamples(OTUfile = files_original_counts, otudir = otudir, removegroups = "HS_Shell",
                                                              combine_metadata = metadata, name = "HS_and_HE_renormalization_table")


renormalized_otu_HS_Control <- Normalization_metagenome_selsamples(OTUfile = files_original_counts, otudir = otudir, removegroups = "H",
                                                              combine_metadata = metadata, name = "HS_and_control_renormalization_table")






#normalized_otu_table <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/normalized_otu_table.rds"

# annotation the group 
# #868686FF GRAY; #CD534CFF LIGHT RED 



  
draw_pheatmap <- function(renormalized_table = NA, combine_metadata = NA, differentOTU = NA, tag = NA){
  

  input_sel_outtable <- renormalized_table %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species, -OTU)
  
  
  group_info <- data.table(uniquebiorep = colnames(input_sel_outtable)) %>%
    left_join(., combine_metadata, by = "uniquebiorep") %>%
    select(., uniquebiorep, Group) %>%
    column_to_rownames(., var = 'uniquebiorep') %>% 
    arrange(., Group)
    
  
  #input_taxnomy <- renormalized_table %>%
  #  select(Kingdom, Phylum, Class, Order, Family, Genus, Species, OTU) %>%
  #  select(OTU, Phylum) %>%
  #  column_to_rownames(., var = "OTU")
  
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

ann_colors = list(
  Group = c("H" = "cornflowerblue", "HS" = "#CD534CFF")#, 
  #Class = c("1" = brewer.pal(n = 10, name = "Paired")[1], "2" = brewer.pal(n = 10, name = "Paired")[2], "3" = brewer.pal(n = 10, name = "Paired")[3], "4" = brewer.pal(n = 10, name = "Paired")[4], "23" = brewer.pal(n = 10, name = "Paired")[5], "32" = brewer.pal(n = 10, name = "Paired")[6], "33" = brewer.pal(n = 10, name = "Paired")[7])
)

HS_HE_allOTU_pheatmap <- draw_pheatmap(renormalized_table = renormalized_otu_HE_HS, combine_metadata = metadata, differentOTU = diffotus_HE_HS, tag = "nofilter")

HS_HE_OTU_pheatmap <- draw_pheatmap(renormalized_table = renormalized_otu_HE_HS, combine_metadata = metadata, 
                                       differentOTU = diffotus_HE_HS, tag = "filter")


pdf(paste(otudir, "HS_and_HE_differenetial_OTU.pdf", sep = "/"),width = 10, height = 12)
HS_HE_OTU_pheatmap
dev.off()
pdf(paste(otudir, "HS_and_HE_allOTU.pdf", sep = "/"),width = 10, height = 12)
HS_HE_allOTU_pheatmap
dev.off()


ann_colors = list(
  Group = c("HS" = "#CD534CFF", "HS_Shell" = "#868686FF")#, 
  #Class = c("1" = brewer.pal(n = 10, name = "Paired")[1], "2" = brewer.pal(n = 10, name = "Paired")[2], "3" = brewer.pal(n = 10, name = "Paired")[3], "4" = brewer.pal(n = 10, name = "Paired")[4], "23" = brewer.pal(n = 10, name = "Paired")[5], "32" = brewer.pal(n = 10, name = "Paired")[6], "33" = brewer.pal(n = 10, name = "Paired")[7])
)

HS_control_all_pheatmap <- draw_pheatmap(renormalized_table = renormalized_otu_HS_Control, combine_metadata = metadata,
                                     differentOTU = diffotus_HS_control, tag = "nofilter")
HS_control_pheatmap <- draw_pheatmap(renormalized_table = renormalized_otu_HS_Control, combine_metadata = metadata,
                                     differentOTU = diffotus_HS_control, tag = "filter")

pdf(paste(otudir, "HS_control_differenetial_OTU.pdf", sep = "/"),width = 10, height = 12)
HS_control_pheatmap
dev.off()

pdf(paste(otudir, "HS_control_all_OTU.pdf", sep = "/"),width = 10, height = 12)
HS_control_all_pheatmap
dev.off()


# part3 :prepare the input & normalization for phylum, families, genus level separately 
#=================================================================================================================================

# produce different levels table 
#====================
produce_different_taxa <- function(otu_norm, col){
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
  add_column(OTUID = .$OTU) %>%
  select(.,-OTU)
col <- ncol(temp)
colnames(temp)[col] <- "OTU"

otu_combined_count <- produce_different_taxa(temp, col)

# separate otu table depending on samples 
#====================

table <- otu_combined_count[[2]]
metatable_sel_biorep <- metadata
key_word <- c(32)


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
  otu_combined_count_sub[[i]][["HS_control"]] <- separate_subgroups(otu_combined_count[[i]],
                                                    metatable_sel_biorep = metadata, key_word = c(32))
  otu_combined_count_sub[[i]][["HS"]] <- separate_subgroups(otu_combined_count[[i]],
                                                          metatable_sel_biorep = metadata, key_word = c(23,33))
  otu_combined_count_sub[[i]][["HE"]] <- separate_subgroups(otu_combined_count[[i]],
                                                        metatable_sel_biorep = metadata, key_word = c(1,2,3,4))
  
}



# # differentially 
#====================
metatable_sel_biorep <- metadata

metagenomeseq_diff_taxa <- function(input_list = NA, metadata = NA){
 
  otu_HScontrol_taxa <- input_list[["HS_control"]]
  otu_HS_taxa <- input_list[["HS"]]
  otu_HE_taxa <- input_list[["HE"]]
  
  diffotus_HS_control_taxa <- metagenomeseq_difference_abundance(Dataset_1 = otu_HS_taxa, Dataset_2 = otu_HScontrol_taxa, 
                                                            Combined_metadata = metadata) %>%
    mutate(.,Taxa = rownames(.)) %>% mutate(.,Group_Impact = GroupHS_Shell)
  
  diffotus_HE_HS_taxa <- metagenomeseq_difference_abundance(Dataset_1 = otu_HE_taxa, Dataset_2 = otu_HS_taxa, 
                                                            Combined_metadata = metadata) %>%
    mutate(.,Taxa = rownames(.)) %>% mutate(.,Group_Impact = GroupHS)
  
  
  
  result <- list("HS_vs_control" = diffotus_HS_control_taxa,
                 "HE_vs_HS" = diffotus_HE_HS_taxa)
  
  
  return(result)
}


otu_combined_count_sub_result <- list()
otudir
for (i in names(otu_combined_count_sub)){
  print(i)
  otu_combined_count_sub_result[[i]] <- metagenomeseq_diff_taxa(input_list = otu_combined_count_sub[[i]], 
                                                                metadata = metadata)
  write_csv(otu_combined_count_sub_result[[i]][["HS_vs_control"]], paste(otudir, paste("Taxonomy/sig", i, "HS_vs_control.csv", sep = "_"), sep = "/"))
  write_csv(otu_combined_count_sub_result[[i]][["HE_vs_HS"]], paste(otudir, paste("Taxonomy/sig", i, "HS_vs_HE.csv", sep = "_"), sep = "/"))
}

# Part4 : draw heatmap 
#==========================================================================================

# renormalization
#====================


Normalization_metagenome <- function(OTUfile, otudir, name, removegroups = NA, metadata = NA){
  
  keep_samples <- metadata[metadata$Group != removegroups,]$uniquebiorep
  
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
  
  Phylum_famiy_genus_norm[[i]][["HS_vs_control"]] <- Normalization_metagenome(otu_combined_count[[i]],
                                                           otudir = paste(otudir, "Taxonomy/", sep = "/"),
                                                           name = paste(i, "HSvsctrl", sep = "_"),
                                                           removegroups = "H",
                                                           metadata = metadata)
  
  Phylum_famiy_genus_norm[[i]][["HE_vs_HS"]] <- Normalization_metagenome(otu_combined_count[[i]],
                                                                              otudir = paste(otudir, "Taxonomy/", sep = "/"),
                                                                              name = paste(i, "HEvsHS", sep = "_"),
                                                                              removegroups = "HS_Shell",
                                                                              metadata = metadata)
  
  
}


### draw heatmap 


draw_pheatmap_taxa <- function(renormalized_table = NA, combine_metadata = NA, differentOTU = NA, tag = NA, annotation_tag = NA){
  
  
  input_sel_outtable <- renormalized_table %>%
    column_to_rownames(., var = "OTUID")
  
  
  group_info <- data.table(uniquebiorep = colnames(input_sel_outtable)) %>%
    left_join(., combine_metadata, by = "uniquebiorep") %>%
    select(., uniquebiorep, Group) %>%
    column_to_rownames(., var = 'uniquebiorep') %>% 
    arrange(., Group)
  
  
  input_taxnomy <- input_sel_outtable %>%
     add_column(taxa = rownames(.)) %>%
     select(taxa)
  colnames(input_taxnomy) <- annotation_tag
 
  input_pheatmap <- input_sel_outtable
  
  if (tag != "nofilter"){
    input_pheatmap_sel <- input_pheatmap %>%
      .[rowSums(.) > 0,] %>%
      .[differentOTU$Taxa,] %>%
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
                                        show_rownames = TRUE,
                                        cluster_cols = F,
                                        fontsize_col = 4, fontsize_row = 4)
  
  
  if (annotation_tag == "OTU"){
    pheatmap_object <- pheatmap::pheatmap(log(input_pheatmap_sel+0.0001), 
                                          annotation_col = group_info,
                                          annotation_colors = ann_colors,
                                          #annotation_row = input_taxnomy,
                                          show_rownames = TRUE,
                                          cluster_cols = F,
                                          fontsize_col = 4, fontsize_row = 4)
  }
  
  
  return(pheatmap_object)
  
  
}


pheatmap_output_taxa <- list()

for (i in names(otu_combined_count)[names(otu_combined_count)!="Kingdom"]){
  print(i)
  
  renormalized_table <- Phylum_famiy_genus_norm[[i]][["HS_vs_control"]]
  combine_metadata <- metadata
  differentOTU <- otu_combined_count_sub_result[[i]][["HS_vs_control"]]
  tag <- "filter"
  ann_colors = list(
    Group = c("HS" = "#CD534CFF", "HS_Shell" = "#868686FF")#, 
    #Class = c("1" = brewer.pal(n = 10, name = "Paired")[1], "2" = brewer.pal(n = 10, name = "Paired")[2], "3" = brewer.pal(n = 10, name = "Paired")[3], "4" = brewer.pal(n = 10, name = "Paired")[4], "23" = brewer.pal(n = 10, name = "Paired")[5], "32" = brewer.pal(n = 10, name = "Paired")[6], "33" = brewer.pal(n = 10, name = "Paired")[7])
  )
  
  pheatmap_output_taxa[[i]][["HS_vs_control"]] <- draw_pheatmap_taxa(renormalized_table = renormalized_table,
                combine_metadata = combine_metadata, 
                differentOTU = differentOTU, tag = "filter", annotation_tag = i)
  
  pdf(paste(otudir, "/Taxonomy/", "HS_vs_control", "diff.", i, ".pdf", sep = ""), width = 8, height = 10)
  print(pheatmap_output_taxa[[i]][["HS_vs_control"]])
  dev.off()

  renormalized_table <- Phylum_famiy_genus_norm[[i]][["HE_vs_HS"]]
  combine_metadata <- metadata
  differentOTU <- otu_combined_count_sub_result[[i]][["HE_vs_HS"]]
  tag <- "filter"
  ann_colors = list(
    Group = c("H" = "cornflowerblue", "HS" = "#CD534CFF")#, 
    #Class = c("1" = brewer.pal(n = 10, name = "Paired")[1], "2" = brewer.pal(n = 10, name = "Paired")[2], "3" = brewer.pal(n = 10, name = "Paired")[3], "4" = brewer.pal(n = 10, name = "Paired")[4], "23" = brewer.pal(n = 10, name = "Paired")[5], "32" = brewer.pal(n = 10, name = "Paired")[6], "33" = brewer.pal(n = 10, name = "Paired")[7])
  )
  
  
  pheatmap_output_taxa[[i]][["HE_vs_HS"]] <- draw_pheatmap_taxa(renormalized_table = renormalized_table,
                                                                combine_metadata = combine_metadata, 
                                                                differentOTU = differentOTU, tag = "filter", annotation_tag = i)
  
  pdf(paste(otudir, "/Taxonomy/", "HE_vs_HS", "diff.", i, ".pdf", sep = ""), width = 8, height = 10)
  print(pheatmap_output_taxa[[i]][["HE_vs_HS"]])
  dev.off()
  
}


save.image(paste(otudir, "metagenomeseq.RData", sep = "/"))

#################### 2022.07 plot modification of heatmap #########################################
otudir <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/"
load(paste(otudir, "metagenomeseq.RData", sep = "/"))


pheatmap_output_taxa_update <- list()
# c("#386cb0", "#f0027f")
for (j in c("Phylum", "Family", "Genus")){ # names(otu_combined_count)[names(otu_combined_count)!="Kingdom"]
  print(j)
  
  renormalized_tablev2 <- Phylum_famiy_genus_norm[[j]][["HS_vs_control"]]
  combine_metadatav2 <- metadata
  differentOTUv2 <- otu_combined_count_sub_result[[j]][["HS_vs_control"]]
  tag <- "filter"
  ann_colors = list(
    Group = c("HS+SH" = "#f0027f", "SH" = "#386cb0")
  )
  
  pheatmap_output_taxa_update[[j]][["HS_vs_control"]] <- draw_pheatmap_taxa_update(renormalized_table = renormalized_tablev2,
                                                                     combine_metadata = combine_metadatav2, 
                                                                     differentOTU = differentOTUv2, tag = "filter", range_value = 2.5)
  
  pdf(paste(otudir, "/Taxonomy/", "HS_vs_control", "diff.", j, "-update.pdf", sep = ""), width = 9, height = 10*1.6)
  print(pheatmap_output_taxa_update[[j]][["HS_vs_control"]])
  dev.off()
  
  
}


draw_pheatmap_taxa_update <- function(renormalized_tablev2 = NA, combine_metadatav2 = NA, differentOTUv2 = NA, tag = NA, range_value = 3){
  
  
  input_sel_outtable <- renormalized_tablev2 %>%
    column_to_rownames(., var = "OTUID")
  
  
  group_info <- data.table(uniquebiorep = colnames(input_sel_outtable)) %>%
    left_join(., combine_metadatav2, by = "uniquebiorep") %>%
    select(., uniquebiorep, Group) %>%
    column_to_rownames(., var = 'uniquebiorep') %>% 
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
  
  # change group name 
  group_info_new <- group_info
  group_info_new$Group <- ifelse(group_info_new$Group == "HS", "HS+SH", "SH")

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
                                        annotation_col = group_info_new,
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

save.image(paste(otudir, "metagenomeseq_202207_figuremodification.RData", sep = "/"))


