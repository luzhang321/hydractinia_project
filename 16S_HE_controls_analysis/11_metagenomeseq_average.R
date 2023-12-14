# 10 metagenomseq 

# metagenomseq for differentially OTUs 


library(tidyverse)
library(fs)
library(vegan)
library(data.table)
library(ggpubr)
library(rstatix)
library(Rfit)
library(viridis)
library(hrbrthemes)
library(VennDiagram)
library(metagenomeSeq)
library(RColorBrewer)


# read my OTU file into the RData 
#===================================================================================================================================
# this is for whole count 


setwd("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/")
files_original_counts <- "/media/lu//Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/10_average_count_table.csv"

# read file 


otu <- read_csv(files_original_counts)
# otu - taxonomy 
#=======================================================================
colnames(otu)[1] <- "OTUID"

otu_tax <- select(otu, OTUID, Kingdom, Phylum, Class,
                  Order, Family, Genus, Species)

#otu_tax <- read_csv("otu_tax.csv")
colnames(otu_tax)[1] <- "OTUID"


# read metatable 
#======================================================================
metatable_sel_biorep <- read_csv("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/metatable_biorep_unique.csv")

head(metatable_sel_biorep)

# form the format of the function 

otu <- otu %>%
  select(., -Kingdom, -Phylum, -Class,-Order, -Family, -Genus, -Species)

# separate otu table 
#=======================================================================


otu_hydrac <- otu %>%
  pivot_longer(., -OTUID, names_to = "uniquebiorep", values_to = "count") %>%
  pivot_wider(., names_from = "OTUID", values_from = "count") %>%
  left_join(., metatable_sel_biorep[,c("uniquebiorep","Class")]) %>%
  filter(., Class %in% c(1,2,3,4,5)) %>% #90 samples 
  column_to_rownames(., var = "uniquebiorep") %>% 
  select(., -Class)


otu_control1 <- otu %>%
  pivot_longer(., -OTUID, names_to = "uniquebiorep", values_to = "count") %>%
  pivot_wider(., names_from = "OTUID", values_from = "count") %>%
  left_join(., metatable_sel_biorep[,c("uniquebiorep","Class")]) %>%
  filter(., Class %in% c(9,10)) %>% #14 samples 
  column_to_rownames(., var = "uniquebiorep") %>% 
  select(., -Class)

otu_control2 <- otu %>%
  pivot_longer(., -OTUID, names_to = "uniquebiorep", values_to = "count") %>%
  pivot_wider(., names_from = "OTUID", values_from = "count") %>%
  left_join(., metatable_sel_biorep[,c("uniquebiorep","Class")]) %>%
  filter(., Class %in% 12) %>% #11 samples 
  column_to_rownames(., var = "uniquebiorep") %>% 
  select(., -Class)


otu_control3 <- otu %>%
  pivot_longer(., -OTUID, names_to = "uniquebiorep", values_to = "count") %>%
  pivot_wider(., names_from = "OTUID", values_from = "count") %>%
  left_join(., metatable_sel_biorep[,c("uniquebiorep","Class")]) %>%
  filter(., Class %in% 13) %>% #5 samples 
  column_to_rownames(., var = "uniquebiorep") %>% 
  select(., -Class)


# metagenomseq function 

metagenomeseq_difference_abundance <- function(Dataset_1,Dataset_2,prevCutoff = 0.1,Combined_metadata) {
  
  #Dataset_1 = otu_hydrac
  #Dataset_2 = otu_control1
  #Combined_metadata =  metatable_sel_biorep
  dataset_name_1 = deparse(substitute(Dataset_1))
  dataset_name_2 = deparse(substitute(Dataset_2))
  # 
  ### Setting the sets
  #full_features = bind_rows(Dataset_1[,2:ncol(Dataset_1)], Dataset_2[,2:ncol(Dataset_2)])
  #full_features[is.na(full_features)] <- 0
  
  full_features = rbind(Dataset_1, Dataset_2)
  
  ## setting the metadata samples that exist in the current comparison
  Combined_metadata <- Combined_metadata %>%
    column_to_rownames(., var = "uniquebiorep")
  full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]
  
  
  # full_metadata$Groups = factor(c(rep(dataset_name_1,nrow(Dataset_1)),
  #                                 rep(dataset_name_2,nrow(Dataset_2))))
  
  ### Order them 
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
  mod <-  model.matrix(~ Group + Crab_presence + Location_Collection, data = pData(tmpMR))
  
  
  ## FitZig =  zero inflated model ideal for metaphlan data
  results_zeller <-  fitZig(tmpMR, mod)
  
  ## Results
  result = MRtable(results_zeller, number = nrow(results_zeller))
  result = result[result$adjPvalues < 0.05,]
  
  return(result)
}  



control1_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac, otu_control1, Combined_metadata = metatable_sel_biorep) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
control2_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac, otu_control2, Combined_metadata = metatable_sel_biorep) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
control3_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac, otu_control3, Combined_metadata = metatable_sel_biorep) %>%
  mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)




dim(control1_diff_otus)
dim(control2_diff_otus)
dim(control3_diff_otus)


## Venn plot 
#=======================================================================

# Load library

outdir <- "/media/lu//Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/"

# Generate 3 sets of 200 words
Shell <- control1_diff_otus$OTU
Sand <- control2_diff_otus$OTU
Stone <- control3_diff_otus$OTU

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
draw_venn_plot <- function(Shell, Sand, Stone, outdir, PngName){
  venn.diagram(
    x = list(Shell, Sand, Stone),
    category.names = c("Shell" , "Sand" , "Stone"),
    filename = paste(outdir,PngName,sep = "/"),
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
}

outdir <- "/media/lu//Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/"

outdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/"
draw_venn_plot(Shell, Sand, Stone, outdir, '#1-Venn-OTU.png')
# updates 2022,03 
draw_venn_plot_update(input = list(control1_diff_otus, control2_diff_otus, control3_diff_otus), 
               outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
               PngName = "#1-Venn-OTU.png", main = "OTU Level") # this function see below


# heatmap of the differetially OTU with matched genus information 
#=======================================================================

# core microbiome 

filedir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_normalized_otu_table.csv"

otu_norm <- read_csv(filedir) %>%
  add_column(OTUID = .$OTU)
otu_norm_sel <- otu_norm[colnames(otu)]


intersect_OTUs <- data.table(OTUID = intersect(Shell, Sand) %>% intersect(., Stone)) %>%
  left_join(., otu_tax) %>%
  left_join(., otu_norm_sel) %>%
  add_column(GenusOTU = paste(.$Genus,"(",.$OTUID,")",sep = "")) %>%
  select(., -Phylum, -OTUID, -Kingdom, -Class, -Order, -Family, -Genus, -Species) %>%
  column_to_rownames(., var = "GenusOTU")
# annotation the group 

ann_colors = list(
  Group = c(" H. echinata" = "cornflowerblue", "Control" = "darkgoldenrod1"),
  Class = c("1" = brewer.pal(n = 10, name = "Paired")[1], "2" = brewer.pal(n = 10, name = "Paired")[2], "3" = brewer.pal(n = 10, name = "Paired")[3], "4" = brewer.pal(n = 10, name = "Paired")[4], "5" = brewer.pal(n = 10, name = "Paired")[5], "9" = brewer.pal(n = 10, name = "Paired")[6], "10" = brewer.pal(n = 10, name = "Paired")[7], "12" = brewer.pal(n = 10, name = "Paired")[9], "13" = brewer.pal(n = 11, name = "Paired")[10])
)



group_info <- data.table(uniquebiorep = colnames(intersect_OTUs)) %>%
  left_join(., metatable_sel_biorep) %>%
  select(., uniquebiorep, Group) %>%
  column_to_rownames(., var = 'uniquebiorep')
#group_info$Class <- as.character(group_info$Class)


otu_info <- data.table(OTUID = intersect(Shell, Sand) %>% intersect(., Stone)) %>%
  left_join(., otu_tax) %>%
  left_join(., otu_norm_sel) %>%
  add_column(GenusOTU = paste(.$Genus,"(",.$OTUID,")",sep = "")) %>%
  select(., Phylum, GenusOTU) %>% # Phylum, OTUID, Kingdom, Class, Order, Family, Genus, Species, GenusOTU ,Family, Genus, 
  column_to_rownames(., var = "GenusOTU")




# draw heatmap 
pdf("#2differenetial_OTU.pdf",width = 10, height = 10)
pheatmap::pheatmap(log(intersect_OTUs+0.001), 
                   annotation_col = group_info,
                   annotation_colors = ann_colors,
                   annotation_row = otu_info,
                   show_rownames = FALSE,
                   fontsize_col = 4)

dev.off()


# Genus/Family/Phylum level 
#=======================================================================

# produce different levels 
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
temp <- read_csv(file = files_original_counts) %>%
  add_column(OTUID = .$OTU) %>%
  select(.,-OTU)
colnames(temp)[70] <- "OTU"

otu_combined_count <- produce_different_taxa(temp, 70)


# separate otu table 
#====================

separate_subgroups <- function(table, metatable_sel_biorep){

  colnames(metatable_sel_biorep)[4] <- "SampleID"
  
  otu_hydrac <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% c(1,2,3,4,5)) %>% 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  otu_control1 <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% c(9,10)) %>%
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  otu_control2 <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% 12) %>% 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  
  otu_control3 <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% 13) %>% 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  sep_table <- list()
  sep_table[["hydract"]] <- otu_hydrac
  sep_table[["control1"]] <- otu_control1
  sep_table[["control2"]] <- otu_control2
  sep_table[["control3"]] <- otu_control3
  
  return(sep_table)
}


# separate
otu_combined_count_sub <- list()
for (i in c("Phylum","Family","Genus")){
  otu_combined_count_sub[[i]] <- separate_subgroups(otu_combined_count[[i]],
                                                    metatable_sel_biorep = metatable_sel_biorep)
  
  
}


# differentially 
#====================

metagenomeseq_diff_taxa <- function(input_list, metatable_sel_biorep, outdir = outdir, PngName){
  otu_hydrac <- input_list[[1]]
  otu_control1 <- input_list[[2]]
  otu_control2 <- input_list[[3]]
  otu_control3 <- input_list[[4]]
  
  control1_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac, otu_control1, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
  control2_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac, otu_control2, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
  control3_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac, otu_control3, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
  
  
  
  Shell <- control1_diff_otus$OTU
  Sand <- control2_diff_otus$OTU
  Stone <- control3_diff_otus$OTU
  
  draw_venn_plot(Shell, Sand, Stone, outdir, PngName)
  
  OTUID = intersect(Shell, Sand) %>% intersect(., Stone)
  
  result <- list(control1_diff_otus,
                 control2_diff_otus,
                 control3_diff_otus,
                 OTUID)
  
  
  return(result)
}

otu_combined_count_sub_result <- list()
for (i in names(otu_combined_count_sub)){
  otu_combined_count_sub_result[[i]] <- metagenomeseq_diff_taxa(otu_combined_count_sub[[i]], metatable_sel_biorep, outdir = outdir, PngName = paste("#2-Venn-", i, ".png",sep = ""))
}


# write all this results into excel file 
# phylum wirte : shell sand stone order 
library(tidyverse)
write_csv(otu_combined_count_sub_result[[1]][[1]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_phylum_shell.csv")
write_csv(otu_combined_count_sub_result[[1]][[2]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_phylum_sand.csv")
write_csv(otu_combined_count_sub_result[[1]][[3]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_phylum_stone.csv")



# family write 
write_csv(otu_combined_count_sub_result[[2]][[1]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_family_shell.csv")
write_csv(otu_combined_count_sub_result[[2]][[2]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_family_sand.csv")
write_csv(otu_combined_count_sub_result[[2]][[3]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_family_stone.csv")


# genus write 

write_csv(otu_combined_count_sub_result[[3]][[1]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_genus_shell.csv")
write_csv(otu_combined_count_sub_result[[3]][[2]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_genus_sand.csv")
write_csv(otu_combined_count_sub_result[[3]][[3]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_genus_stone.csv")










# update on 2022.03, venn plot group name change 
library(VennDiagram)
# List of items

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
input <- otu_combined_count_sub_result[[1]]
outdir <- "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/"
PngName <- "#2_venn_testtest.png"

draw_venn_plot_update(input = otu_combined_count_sub_result[[1]], 
                      outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                      PngName = "#2-Venn-Phylum.png", main = "Phylum Level")

draw_venn_plot_update(input = otu_combined_count_sub_result[[2]], 
                      outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                      PngName = "#2-Venn-Family.png", main = "Family Level")

draw_venn_plot_update(input = otu_combined_count_sub_result[[3]], 
                      outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                      PngName = "#2-Venn-Genus.png", main = "Genus Level")




draw_venn_plot_update <- function(input = NA, outdir = NA, PngName = NA, main = NA){
  venn.diagram(
    x = list(input[[1]]$OTU, input[[2]]$OTU, input[[3]]$OTU),
    category.names = c("HE vs Shell" , "HE vs Sand" , "HE vs Stone"),
    filename = paste(outdir,PngName,sep = "/"),
    output=TRUE,
    main = main,
    main.cex = 0.5,
    main.pos = c(0.5,1),
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 1,
    #lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .4,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-10, 10, 5),
    cat.dist = c(0.055, 0.055, 0.055),
    cat.fontfamily = "sans",
    rotation = 1
  )
}


# 2022.07 updates Venn plot color 
"#dede18" - SH/SH+CR
"#a77be3" - SA
"#f5973d" - ST


draw_venn_plot_update_color_phylum(input = otu_combined_count_sub_result[[1]], 
                      outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                      PngName = "#2-Venn-Phylum-v2.svg", main = "Phylum Level")

draw_venn_plot_update_color_family(input = otu_combined_count_sub_result[[2]], 
                      outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                      PngName = "#2-Venn-Family-v2.svg", main = "Family Level")

draw_venn_plot_update_color_family(input = otu_combined_count_sub_result[[3]], 
                      outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                      PngName = "#2-Venn-Genus-v2.svg", main = "Genus Level")


load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/11_metagenomeseq.RData")
saveRDS(otu_combined_count_sub_result[[1]], "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/2_Phylum_venn_input.rds")
saveRDS(otu_combined_count_sub_result[[2]], "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/2_Family_venn_input.rds")
saveRDS(otu_combined_count_sub_result[[3]], "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/2_Genus_venn_input.rds")


draw_venn_plot_update_color_phylum <- function(input = NA, outdir = NA, PngName = NA, main = NA){
  p <- venn.diagram(
    x = list(input[[1]]$OTU, input[[2]]$OTU, input[[3]]$OTU),
    category.names = c("HE vs SH/SH+CR" , "HE vs SA" , "HE vs ST"),
    filename = NULL,
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5,1),
    #col=c('#dede18', '#a77be3', '#f5973d'),
    
    # Output features
    imagetype="svg" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 5,
    #lty = 'blank',
    fill = c('#dede18', '#a77be3', '#f5973d'),
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-10, 10, 5),
    cat.dist = c(0.055, 0.055, 0.055),
    cat.fontfamily = "sans",
    rotation = 1
  )
  
  ggsave(p, file=paste(outdir,PngName,sep = "/"), device = "svg")
  
}

# for family 
draw_venn_plot_update_color_family <- function(input = NA, outdir = NA, PngName = NA, main = NA){
  p <- venn.diagram(
    x = list(input[[1]]$OTU, input[[2]]$OTU, input[[3]]$OTU),
    category.names = c("HE vs SH/SH+CR" , "HE vs SA" , "HE vs ST"),
    filename = NULL,
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5,1),
    #col=c('#dede18', '#a77be3', '#f5973d'),
    
    # Output features
    imagetype="svg" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 5,
    #lty = 'blank',
    fill = c('#dede18', '#a77be3', '#f5973d'),
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-10, 10, 5),
    cat.dist = c(0.035, 0.055, 0.055),
    cat.fontfamily = "sans",
    rotation = 1
  )
  
  ggsave(p, file=paste(outdir,PngName,sep = "/"), device = "svg")
  
}

# i didn't save updated image for this 







# i also need all information of different phyla, not only the significant ones 
####

metagenomeseq_difference_abundance_nofilter <- function(Dataset_1,Dataset_2,prevCutoff = 0.1,Combined_metadata) {
  
  #Dataset_1 = otu_hydrac
  #Dataset_2 = otu_control1
  #Combined_metadata =  metatable_sel_biorep
  dataset_name_1 = deparse(substitute(Dataset_1))
  dataset_name_2 = deparse(substitute(Dataset_2))
  # 
  ### Setting the sets
  #full_features = bind_rows(Dataset_1[,2:ncol(Dataset_1)], Dataset_2[,2:ncol(Dataset_2)])
  #full_features[is.na(full_features)] <- 0
  
  full_features = rbind(Dataset_1, Dataset_2)
  
  ## setting the metadata samples that exist in the current comparison
  Combined_metadata <- Combined_metadata %>%
    column_to_rownames(., var = "uniquebiorep")
  full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]
  
  
  # full_metadata$Groups = factor(c(rep(dataset_name_1,nrow(Dataset_1)),
  #                                 rep(dataset_name_2,nrow(Dataset_2))))
  
  ### Order them 
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
  mod <-  model.matrix(~ Group + Crab_presence + Location_Collection, data = pData(tmpMR))
  
  
  ## FitZig =  zero inflated model ideal for metaphlan data
  results_zeller <-  fitZig(tmpMR, mod)
  
  ## Results
  result = MRtable(results_zeller, number = nrow(results_zeller))
  #result = result[result$adjPvalues < 0.05,]
  
  return(result)
}  

metagenomeseq_diff_taxa_no_filter <- function(input_list, metatable_sel_biorep, outdir = outdir, PngName){
  otu_hydrac <- input_list[[1]]
  otu_control1 <- input_list[[2]]
  otu_control2 <- input_list[[3]]
  otu_control3 <- input_list[[4]]
  
  control1_diff_otus <- metagenomeseq_difference_abundance_nofilter(otu_hydrac, otu_control1, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
  control2_diff_otus <- metagenomeseq_difference_abundance_nofilter(otu_hydrac, otu_control2, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
  control3_diff_otus <- metagenomeseq_difference_abundance_nofilter(otu_hydrac, otu_control3, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
  
  
  
  Shell <- control1_diff_otus$OTU
  Sand <- control2_diff_otus$OTU
  Stone <- control3_diff_otus$OTU
  
  #draw_venn_plot(Shell, Sand, Stone, outdir, PngName)
  
  OTUID = intersect(Shell, Sand) %>% intersect(., Stone)
  
  result <- list(control1_diff_otus,
                 control2_diff_otus,
                 control3_diff_otus,
                 OTUID)
  
  
  return(result)
}

otu_combined_count_sub_result_no_filter <- list()
for (i in names(otu_combined_count_sub)){
  otu_combined_count_sub_result_no_filter[[i]] <- metagenomeseq_diff_taxa_no_filter(otu_combined_count_sub[[i]], metatable_sel_biorep, outdir = outdir, PngName = paste("#2-Venn-", i, ".png",sep = ""))
}

# write the unfiltered information 
# phylum wirte : shell sand stone order 
library(tidyverse)
write_csv(otu_combined_count_sub_result_no_filter[[1]][[1]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_phylum_shell.csv")
write_csv(otu_combined_count_sub_result_no_filter[[1]][[2]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_phylum_sand.csv")
write_csv(otu_combined_count_sub_result_no_filter[[1]][[3]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_phylum_stone.csv")


# family write 
write_csv(otu_combined_count_sub_result_no_filter[[2]][[1]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_family_shell.csv")
write_csv(otu_combined_count_sub_result_no_filter[[2]][[2]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_family_sand.csv")
write_csv(otu_combined_count_sub_result_no_filter[[2]][[3]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_family_stone.csv")


# genus write 

write_csv(otu_combined_count_sub_result_no_filter[[3]][[1]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_genus_shell.csv")
write_csv(otu_combined_count_sub_result_no_filter[[3]][[2]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_genus_sand.csv")
write_csv(otu_combined_count_sub_result_no_filter[[3]][[3]] , "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_unfiltered_genus_stone.csv")
















# draw heatmap of core microbiome 
#====================

Normalization_metagenome <- function(OTUfile, otudir, name){
  normfile <- OTUfile %>%    
    metagenomeSeq::newMRexperiment(.) %>%       ## creating MR object                      
    metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
    metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>%
    as_tibble(.) %>%
    add_column(OTUID = rownames(OTUfile)) %>%
    as.data.table() 
  #write_csv(normfile, paste(otudir, name,"_normalized_otu_table.csv", sep = ""))
  return(normfile)
}


Phylum_famiy_genus_norm <- list()
for (i in c("Phylum","Family","Genus")){
  
  Phylum_famiy_genus_norm[[i]] <- Normalization_metagenome(otu_combined_count[[i]],
                                                           otudir = "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                                                           name = i)
}




#a <- read_csv(files[7]) #2 - family; 3- Genus; 7- Phylum

#draw_pheatmap_batch(Phylum_famiy_genus_norm[[1]], 
#                    otu_combined_count_sub_result[["Phylum"]][[4]], 
#                    pdfName  = paste(outdir, "#3-Differential-Phylum.pdf",sep = "/"))


#draw_pheatmap_batch(Phylum_famiy_genus_norm[[2]], 
#                    otu_combined_count_sub_result[["Family"]][[4]], 
#                    pdfName = paste(outdir, "#3-Differential-Family.pdf",sep = "/"))

#draw_pheatmap_batch(Phylum_famiy_genus_norm[[3]], 
#                    otu_combined_count_sub_result[["Genus"]][[4]], 
#                    pdfName = paste(outdir, "#3-Differential-Genus.pdf",sep = "/"))

draw_pheatmap_batch <- function(table, otulist, pdfName){
  
  table$OTUID[1] <- "NA"
  intersect_OTUs <- data.table(OTUID = otulist) %>%
    left_join(., table) %>%
    column_to_rownames(., var = "OTUID")
  #return(intersect_OTUs)
  
  # draw heatmap 
  pdf(pdfName,width = 10, height = 10)
  pheatmap::pheatmap(log(intersect_OTUs+0.001), 
                     annotation_col = group_info,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     fontsize_col = 4)
  
  dev.off()
}

# update 2022.03 

draw_pheatmap_batch_update <- function(table = NA, otulist = NA, pdfName = NA){
  
  table$OTUID[1] <- "NA"
  intersect_OTUs <- data.table(OTUID = otulist) %>%
    left_join(., table) %>%
    column_to_rownames(., var = "OTUID")
  #return(intersect_OTUs)
  group_info_order <- arrange(group_info, Group)
  intersect_OTUs_order <- intersect_OTUs[, rownames(group_info_order)]
  
  # draw heatmap 
  pdf(pdfName,width = 10, height = 10)
  pheatmap::pheatmap(log(intersect_OTUs_order+0.001), 
                     annotation_col = group_info_order,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     cluster_cols = FALSE,
                     fontsize_col = 4)
  
  dev.off()
}

draw_pheatmap_batch_update(Phylum_famiy_genus_norm[[1]], 
                    otu_combined_count_sub_result[["Phylum"]][[4]], 
                    pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/", 
                                     "#3-Differential-Phylum.pdf",sep = "/"))


draw_pheatmap_batch_update(Phylum_famiy_genus_norm[[2]], 
                    otu_combined_count_sub_result[["Family"]][[4]], 
                    pdfName = paste("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                                    "#3-Differential-Family.pdf",sep = "/"))

draw_pheatmap_batch_update(Phylum_famiy_genus_norm[[3]], 
                    otu_combined_count_sub_result[["Genus"]][[4]], 
                    pdfName = paste("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/", 
                                    "#3-Differential-Genus.pdf",sep = "/"))

# update heatmap figure 

# update 2022.05
pdfName <- "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Phylum-update.pdf"
pdfName <- "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Family-update.pdf"
otutable <- Phylum_famiy_genus_norm[[2]]
otulist <- otu_combined_count_sub_result[["Family"]][[4]] #intersect of the significantly different taxon
ann_colors = list(
  Group = c("H. echinata" = "#7fc97f", "Control" = "#f1e2cc")
)
#ann_colors2 = list(
#  Group = c("H" = "#7fc97f", "C" = "#fc8d62")
#)


library(RColorBrewer)
library("RColorBrewer")
library(circlize)

# phylum pdf 
draw_pheatmap_batch_update2(otutable = Phylum_famiy_genus_norm[[1]], 
                            otulist = otu_combined_count_sub_result[["Phylum"]][[4]],
                            pdfName = "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Phylum-update.pdf",
                            width_seletive = 9)

# family pdf 
draw_pheatmap_batch_update2(otutable = Phylum_famiy_genus_norm[[2]], 
                            otulist = otu_combined_count_sub_result[["Family"]][[4]],
                            pdfName = "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Family-update.pdf",
                            width_seletive = 9)

# genus pdf 
draw_pheatmap_batch_update2(otutable = Phylum_famiy_genus_norm[[3]], 
                            otulist = otu_combined_count_sub_result[["Genus"]][[4]],
                            pdfName = "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Genus-update.pdf",
                            width_seletive = 9)


width_seletive <- 9
draw_pheatmap_batch_update2 <- function(otutable = NA, otulist = NA, pdfName = NA, width_seletive = NA){
  
  otutable$OTUID[1] <- "NA"
  intersect_OTUs <- data.table(OTUID = otulist) %>%
    left_join(., otutable) %>%
    column_to_rownames(., var = "OTUID")
  #return(intersect_OTUs)
  group_info_order <- arrange(group_info, Group) %>%
    mutate(Group = ifelse(.$Group == "H", "H. echinata", "Control"))
  intersect_OTUs_order <- intersect_OTUs[, rownames(group_info_order)]
  rownames(intersect_OTUs_order) <- str_remove_all(rownames(intersect_OTUs_order), "p__") %>%
    str_remove_all(., "f__") %>%
    str_remove_all(., "g__")
  
  newnames <- lapply(
    rownames(intersect_OTUs_order),
    function(x) bquote(italic(.(x))))
  
  # draw heatmap 
  pdf(pdfName, width = width_seletive)
  pheatmap::pheatmap(log(intersect_OTUs_order+0.001), 
                     annotation_col = group_info_order,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     cluster_cols = FALSE,
                     #color = my_palette,
                     fontsize_col = 4,
                     fontsize_row = 8,
                     cellheight=12, cellwidth = 6,
                     labels_row = as.expression(newnames))
  
  pheatmap::pheatmap(intersect_OTUs_order, 
                     annotation_col = group_info_order,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     cluster_cols = FALSE,
                     scale = "row",
                     fontsize_col = 5,
                     fontsize_row = 8,
                     cellheight=12, cellwidth = 6,
                     labels_row = as.expression(newnames))
  
  
  dev.off()
}



# update.07.2022, change the color scale of the heatmap 
###############################################################################################################################


library(RColorBrewer)
library("RColorBrewer")
library(circlize)
ann_colors = list(
  Group = c("HE" = "#09ba09", "SH/SH+CR/SA/ST" = "#f1e2cc")
)
# phylum pdf 
phylum_pheatmap <- draw_pheatmap_batch_update3(otutable = Phylum_famiy_genus_norm[[1]], 
                            otulist = otu_combined_count_sub_result[["Phylum"]][[4]],
                            pdfName = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Phylum-updatev2.pdf",
                            width_seletive = 9)

##################### my test code below ##########################################################################
#dat <-as.matrix(Phylum_famiy_genus_norm[[1]]) %>% as.matrix() %>% as.numeric() %>% data.frame(values = .)
#mat <- Phylum_famiy_genus_norm[[1]] %>% column_to_rownames(., var = "OTUID") %>% as.matrix()
#ggplot(dat, aes(values)) + geom_density(bw = "SJ")
## you can see the data, the value is not very high there
#quantile_breaks <- function(xs, n = 10) {
#  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
#  breaks[!duplicated(breaks)]
#}

#color_test <- c("#5a89be", "#89b7d6", "#ebf7e2", "#fef8b5", "#fdcc82", "#f1734b")
library(viridis)
#mat_breaks <- quantile_breaks(mat, n = 11)
#pheatmap(
#  mat               = mat,
  #color             = inferno(length(mat_breaks) - 1),
#  color             = color_test,
#  breaks            = mat_breaks,
#  border_color      = NA,
#  show_colnames     = FALSE,
#  show_rownames     = TRUE,
#  annotation_col    = group_info_order,
#  annotation_colors = ann_colors,
#  drop_levels       = TRUE,
#  fontsize          = 14,
#  main              = "Quantile Color Scale"
#)
##################### my test code above ##########################################################################


# family pdf 
family_pheatmap <-draw_pheatmap_batch_update3(otutable = Phylum_famiy_genus_norm[[2]], 
                            otulist = otu_combined_count_sub_result[["Family"]][[4]],
                            pdfName = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Family--updatev2.pdf",
                            width_seletive = 9)

# genus pdf 
genus_pheatmap <-draw_pheatmap_batch_update3(otutable = Phylum_famiy_genus_norm[[3]], 
                            otulist = otu_combined_count_sub_result[["Genus"]][[4]],
                            pdfName = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#3-Differential-Genus--updatev2.pdf",
                            width_seletive = 9)


width_seletive <- 9
##################### my test code below ##########################################################################

#quantile_breaks <- function(xs, n = 10) {
#  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
#  breaks[!duplicated(breaks)]
#}

#color_test <- c("#5a89be", "#89b7d6", "#ebf7e2", "#fef8b5", "#fdcc82", "#f1734b")
#color_test2 <- c("#5a89be", "#ebf7e2", "#fef8b5", "#f1734b")
#mat_breaks <- quantile_breaks(mat, n = 11)
##################### my test code above ##########################################################################

draw_pheatmap_batch_update3 <- function(otutable = NA, otulist = NA, pdfName = NA, width_seletive = NA){
  
  otutable$OTUID[1] <- "NA"
  intersect_OTUs <- data.table(OTUID = otulist) %>%
    left_join(., otutable) %>%
    column_to_rownames(., var = "OTUID")
  #return(intersect_OTUs)
  group_info_order <- arrange(group_info, Group) %>%
    rownames_to_column(., var = "sample") %>% 
    mutate(Group = ifelse(.$Group == "H", "HE", "SH/SH+CR/SA/ST")) %>%
    column_to_rownames(., var = "sample")
  
  intersect_OTUs_order <- intersect_OTUs[, rownames(group_info_order)]
  rownames(intersect_OTUs_order) <- str_remove_all(rownames(intersect_OTUs_order), "p__") %>%
    str_remove_all(., "f__") %>%
    str_remove_all(., "g__")
  
  newnames <- lapply(
    rownames(intersect_OTUs_order),
    function(x) bquote(italic(.(x))))
  
  
  library(viridis)
  mat <- log(intersect_OTUs_order+0.001) %>% as.matrix()
  mat_breaks2 <- quantile_breaks(mat, n = 11)
  
  # below code doesn't work good 
  #pheatmap::pheatmap(
  #                   log(intersect_OTUs_order+0.001),
  #                   annotation_col = group_info_order,
  #                   annotation_colors = ann_colors,
  #                   show_rownames = TRUE,
  #                   cluster_cols = FALSE,
  #                   color = color_test2,
  #                   fontsize_col = 4,
  #                   fontsize_row = 8,
  #                   cellheight=12, cellwidth = 6,
  #                   breaks = mat_breaks2,
  #                   labels_row = as.expression(newnames))
  
  
  
  intersect_OTUs_order_scale <- t(scale(t(intersect_OTUs_order)))
  # default scale is column scale, scale is generic function whose default method centers and/or scales the columns of a numeric matrix.
  # so here we used z-score to do normalization, then replace the extreme value with the same value and draw plot 
  intersect_OTUs_order_scale[intersect_OTUs_order_scale > 3.5] <- 3.5
  
  
  p <- pheatmap::pheatmap(intersect_OTUs_order_scale, 
                          annotation_col = group_info_order,
                          annotation_colors = ann_colors,
                          show_rownames = TRUE,
                          cluster_cols = FALSE,
                          # scale = "row",
                          fontsize_col = 5,
                          fontsize_row = 8,
                          cellheight=12, cellwidth = 6,
                          labels_row = as.expression(newnames))
  # draw heatmap 
  pdf(pdfName, width = width_seletive)
  print(p)

  dev.off()
  return(p)
}

save.image("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/step11_metagenomeseq_average_plot_modification.RData")

###############################################################################################################################

# extract sig list 

# Phylum 
Phylum_list1 <- select(otu_combined_count_sub_result[["Phylum"]][[1]], OTU, adjPvalues, Group_Impact) %>%
  add_column(Control = "Shell")
Phylum_list2 <- select(otu_combined_count_sub_result[["Phylum"]][[2]], OTU, adjPvalues, Group_Impact) %>%
  add_column(Control = "Sand")
Phylum_list3 <- select(otu_combined_count_sub_result[["Phylum"]][[3]], OTU, adjPvalues, Group_Impact) %>%
  add_column(Control = "Stone")

phylum_list <- rbind(Phylum_list1, Phylum_list2, Phylum_list3)

# Family 
Family_list1 <-select(otu_combined_count_sub_result[["Family"]][[1]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Control = "Shell")
Family_list2 <- select(otu_combined_count_sub_result[["Family"]][[2]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Control = "Sand")
Family_list3 <-select(otu_combined_count_sub_result[["Family"]][[3]], OTU, adjPvalues, Group_Impact) %>%
  add_column(Control = "Stone")

Family_list <- rbind(Family_list1, Family_list2, Family_list3)

# Genus 
Genus_list1 <- select(otu_combined_count_sub_result[["Genus"]][[1]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Control = "Shell")
Genus_list2 <- select(otu_combined_count_sub_result[["Genus"]][[2]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Control = "Sand")
Genus_list3 <- select(otu_combined_count_sub_result[["Genus"]][[3]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Control = "Stone")

Genus_list <- rbind(Genus_list1, Genus_list2, Genus_list3)



save.image("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/11_metagenomeseq.RData")

load("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/11_metagenomeseq.RData")


# 19, Nov, 2022  figure modification on venn plot 
###########################################################################################################

load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/step11_metagenomeseq_average_plot_modification.RData")

library(tidyverse)
library(VennDiagram)
library(gridExtra)



"#dede18" - SH/SH+CR
"#a77be3" - SA
"#f5973d" - ST

# https://stackoverflow.com/questions/43324180/adding-legend-to-venn-diagram
x <- list(c(1,2,3,4,5),c(4,5,6,7,8,9,10))
diag <- venn.diagram(x,NULL,fill=c("#80b1d3","#b3de69"),
                     category.names=c("A","B"),height=500,width=500,res=150)


cols <- c("#80b1d3","#b3de69")
lg <- legendGrob(labels=c("A","B"), pch=rep(19,length(c("A","B"))),
                 gp=gpar(col=cols, fill="gray"),
                 byrow=TRUE)
g <- gTree(children = gList(diag))
gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4,1))


draw_venn_plot_update_color_phylum_with_legend <- function(input = NA, outdir = NA, PngName = NA, main = NA){
  p <- venn.diagram(
    x = list(input[[1]]$OTU, input[[2]]$OTU, input[[3]]$OTU),
    category.names = c("" , "" , ""),
    #category.names = c("HE vs SH/SH+CR" , "HE vs SA" , "HE vs ST"),
    filename = NULL,
    output = FALSE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5,1),
    #col=c('#dede18', '#a77be3', '#f5973d'),
    
    # Output features
    #imagetype="svg" ,
    #height = 480 , 
    #width = 480 , 
    #resolution = 600,
    #compression = "lzw",
    
    # Circles
    lwd = 5,
    #lty = 'blank',
    fill = c('#dede18', '#a77be3', '#f5973d'),
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "Helvetica",
    
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-10, 10, 5),
    cat.dist = c(0.055, 0.055, 0.055),
    cat.fontfamily = "Helvetica", # https://studiokayama.com/5-fonts-that-add-credibility-and-professionalism-to-scientific-research/
    rotation = 1
  )
  
  cols <- c('#dede18', '#a77be3', '#f5973d')
  lg <- legendGrob(labels = c("HE vs SH/SH+CR","HE vs SA", "HE vs ST"), 
                   pch = rep(15, length(c("HE vs SH/SH+CR","HE vs SA", "HE vs ST"))),
                   gp=gpar(col = cols, fill = "gray"),
                   byrow = TRUE)
  
  g <- gTree(children = gList(p))
  final_p <- gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4, 1))
  
  
  ggsave(final_p, file = paste(outdir, PngName, sep = "/"), device = "pdf", width = 9)
  return(final_p)

}# arial error 



phylum_venn <- draw_venn_plot_update_color_phylum_with_legend(input = otu_combined_count_sub_result[[1]], 
                                   outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                                   PngName = "#2-Venn-Phylum-v3.pdf", main = "Phylum Level")


# for family 
draw_venn_plot_update_color_family_with_legend <- function(input = NA, outdir = NA, PngName = NA, main = NA){
  p <- venn.diagram(
    x = list(input[[1]]$OTU, input[[2]]$OTU, input[[3]]$OTU),
    category.names = c("" , "" , ""),
    #category.names = c("HE vs SH/SH+CR" , "HE vs SA" , "HE vs ST"),
    filename = NULL,
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5,1.1),
    #col=c('#dede18', '#a77be3', '#f5973d'),
    
    
    # Circles
    lwd = 5,
    #lty = 'blank',
    fill = c('#dede18', '#a77be3', '#f5973d'),
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "Helvetica",
    
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-10, 10, 5),
    cat.dist = c(0.035, 0.055, 0.055),
    cat.fontfamily = "Helvetica",
    rotation = 1
  )
  
  
  cols <- c('#dede18', '#a77be3', '#f5973d')
  lg <- legendGrob(labels = c("HE vs SH/SH+CR","HE vs SA", "HE vs ST"), 
                   pch = rep(15, length(c("HE vs SH/SH+CR","HE vs SA", "HE vs ST"))),
                   gp=gpar(col = cols, fill = "gray"),
                   byrow = TRUE)
  
  g <- gTree(children = gList(p))
  final_p <- gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4, 1))
  
  
  ggsave(final_p, file = paste(outdir, PngName, sep = "/"), device = "pdf", width = 9)  
  
  return(final_p)

}


family_venn <- draw_venn_plot_update_color_family_with_legend(input = otu_combined_count_sub_result[[2]], 
                                   outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                                   PngName = "#2-Venn-Family-v3.pdf", main = "Family Level")

genus_venn <- draw_venn_plot_update_color_family_with_legend(input = otu_combined_count_sub_result[[3]], 
                                   outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/",
                                   PngName = "#2-Venn-Genus-v3.pdf", main = "Genus Level")


save.image("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/step11_metagenomeseq_average_plot_modificationv2.RData")

load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/step11_metagenomeseq_average_plot_modificationv2.RData")

combine_plot <- ggpubr::ggarrange(phylum_venn, family_venn, genus_venn, ncol = 3)

ggsave("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#2-combine_venn_plot.pdf", combine_plot, width = 28)
save.image("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/step11_metagenomeseq_average_plot_modificationv2.RData")



