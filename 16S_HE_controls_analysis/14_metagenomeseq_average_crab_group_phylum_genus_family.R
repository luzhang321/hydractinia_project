# separate otu table 
#====================

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
#===================================================================================================================================


# read otu table and metatable 
# =============================

setwd("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10")
load("#6-4-otu_combined_count.RData")

metatable_sel_biorep <- read_csv("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/metatable_biorep_unique.csv")


outdir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10"


# separate table into 4 groups 
# =============================

separate_subgroups <- function(table, metatable_sel_biorep){
  
  colnames(metatable_sel_biorep)[4] <- "SampleID"
  
  otu_hydrac_crab <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% c(1,3)) %>% #15 samples 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  otu_hydrac_no_crab <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% c(2,4)) %>% #41 samples 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  otu_control_crab <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% 9) %>% #5 samples 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  
  otu_control_no_crab <- table %>%
    rownames_to_column(., var = "OTUID") %>%
    pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
    pivot_wider(., names_from = "OTUID", values_from = "count") %>%
    left_join(., metatable_sel_biorep[,c("SampleID","Class")]) %>%
    filter(., Class %in% 10) %>% #15 samples 
    column_to_rownames(., var = "SampleID") %>% 
    select(., -Class)
  
  sep_table <- list()
  sep_table[["hydract-crab"]] <- otu_hydrac_crab #15 
  sep_table[["hydract-nocrab"]] <- otu_hydrac_no_crab #11
  sep_table[["control-crab"]] <- otu_control_crab #5
  sep_table[["control-nocrab"]] <- otu_control_no_crab #9
  
  return(sep_table)
}


# separate code 
otu_combined_count_sub <- list()
for (i in c("Phylum","Family","Genus")){
  otu_combined_count_sub[[i]] <- separate_subgroups(otu_combined_count[[i]],
                                                    metatable_sel_biorep = metatable_sel_biorep)
  
  
}

# differentially by using metagenomeseq 
#====================

# venn plot function
draw_venn_plot <- function(WithCrab, NoCrab, outdir, name){
  venn.diagram(
    x = list(WithCrab, NoCrab),
    category.names = c("WithCrab : \n Hydrac vs Control" , "NoCrab : \n Hydrac vs Control"),
    filename = paste(outdir, name, sep = "/"),
    output=TRUE, 
    # Output features
    #imagetype="pdf" ,
    resolution = 300,
    compression = "lzw",
    height = 1880 , 
    width = 1880 , 
    
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = brewer.pal(3, "Pastel2")[1:2],
    
    # Numbers
    cex = 1.6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.0,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-5, 5),
    cat.dist = c(0.035, 0.035),
    cat.fontfamily = "sans"
  )
}

metagenomeseq_difference_abundance <- function(Dataset_1,Dataset_2,prevCutoff = 0.1,Combined_metadata) {
  
  #Dataset_1 = otu_hydrac_crab 
  #Dataset_2 = otu_control_with_crab
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
  
  
  #full_metadata$Group = factor(c(rep(dataset_name_1,nrow(Dataset_1)),
  #                                 rep(dataset_name_2,nrow(Dataset_2))))
  full_metadata$Group <- factor(full_metadata$Group)
  
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
  prof_mr <- filterData(prof_mr, depth = 5, present = prevalence)
  
  ## Cumulative Sum Scaling Normalisation
  p <- cumNormStatFast(prof_mr)
  tmpMR <- cumNorm(prof_mr, p)
  
  #NORMFACTOR <- normFactors(prof_mr) # lU ADDING 
  ## Creating the model 
  mod <-  model.matrix(~ Group + Location_Collection, data = pData(tmpMR))
  
  
  ## FitZig =  zero inflated model ideal for metaphlan data
  results_zeller <-  fitZig(tmpMR, mod)
  
  ## Results
  result = MRtable(results_zeller, number = nrow(results_zeller))
  result = result[result$adjPvalues < 0.05,]
  
  return(result)
}  


metagenomeseq_diff_taxa <- function(input_list, metatable_sel_biorep, outdir = outdir, PngName){
  otu_hydrac_crab <- input_list[[1]]
  otu_hydrac_no_crab <- input_list[[2]]
  otu_control_crab <- input_list[[3]]
  otu_control_no_crab <- input_list[[4]]
  
  control1_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac_crab, otu_control_crab, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)
  control2_diff_otus <- metagenomeseq_difference_abundance(otu_hydrac_no_crab, otu_control_no_crab, Combined_metadata = metatable_sel_biorep) %>%
    mutate(.,OTU = rownames(.)) %>% mutate(.,Group_Impact = GroupH)

  
  
  
  crab <- control1_diff_otus$OTU
  nocrab <- control2_diff_otus$OTU
 
  
  draw_venn_plot(crab, nocrab, outdir, PngName)
  
  
  result <- list(control1_diff_otus,
                 control2_diff_otus,
                 crab,
                 nocrab)
  
  
  return(result)
}

otu_combined_count_sub_result <- list()
for (i in names(otu_combined_count_sub)){
  otu_combined_count_sub_result[[i]] <- metagenomeseq_diff_taxa(otu_combined_count_sub[[i]], metatable_sel_biorep, outdir = outdir, PngName = paste("#6-5-Venn-", i, ".pdf",sep = ""))
}

# write these sig taxa into csv file ( add : 2021.05.20 )

library(tidyverse)

# phylum 
write_csv(otu_combined_count_sub_result[[1]][[1]], "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/6_sig_crab_phylum.csv")
write_csv(otu_combined_count_sub_result[[1]][[2]], "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/6_sig_nocrab_phylum.csv")

# family 
write_csv(otu_combined_count_sub_result[[2]][[1]], "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/6_sig_crab_family.csv")
write_csv(otu_combined_count_sub_result[[2]][[2]], "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/6_sig_nocrab_family.csv")

# genus 
write_csv(otu_combined_count_sub_result[[3]][[1]], "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/6_sig_crab_genus.csv")
write_csv(otu_combined_count_sub_result[[3]][[2]], "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/6_sig_nocrab_genus.csv")




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
  write_csv(normfile, paste(otudir, "#6-6-",name,"_normalized_otu_table.csv", sep = ""))
  return(normfile)
}


Phylum_famiy_genus_norm <- list()
for (i in c("Phylum","Family","Genus")){
  
  Phylum_famiy_genus_norm[[i]] <- Normalization_metagenome(otu_combined_count[[i]],
                                                           otudir = "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/",
                                                           name = i)
}




#a <- read_csv(files[7]) #2 - family; 3- Genus; 7- Phylum


# prepare the metatable 
group_info_with_crab <- filter(metatable_sel_biorep, Class %in% c(1,3,9)) %>% 
  select(., uniquebiorep, Group, Class) %>%
  column_to_rownames(., var = 'uniquebiorep')

write_csv(group_info_with_crab, paste(outdir, "#6-8-metatable-crab.csv", sep = "/"))

group_info_without_crab <- filter(metatable_sel_biorep, Class %in% c(2,4,10)) %>% 
  select(., uniquebiorep, Group, Class) %>%
  column_to_rownames(., var = 'uniquebiorep')

write_csv(group_info_without_crab, paste(outdir, "#6-8-metatable-no-crab.csv", sep = "/"))


# draw heatmap 


draw_pheatmap_batch(Phylum_famiy_genus_norm[[1]], 
                    otu_combined_count_sub_result[["Phylum"]][[3]], 
                    pdfName  = paste(outdir, "#6-7-Differential-Phylum-crab.pdf",sep = "/"),
                    group_info = group_info_with_crab)

draw_pheatmap_batch(Phylum_famiy_genus_norm[[1]], 
                    otu_combined_count_sub_result[["Phylum"]][[4]], 
                    pdfName  = paste(outdir, "#6-7-Differential-Phylum-without-crab.pdf",sep = "/"),
                    group_info = group_info_without_crab)



draw_pheatmap_batch(Phylum_famiy_genus_norm[[2]], 
                    otu_combined_count_sub_result[["Family"]][[3]], 
                    pdfName = paste(outdir, "#6-7-Differential-Family-crab.pdf",sep = "/"),
                    group_info = group_info_with_crab)

draw_pheatmap_batch(Phylum_famiy_genus_norm[[2]], 
                    otu_combined_count_sub_result[["Family"]][[4]], 
                    pdfName = paste(outdir, "#6-7-Differential-Family-without-crab.pdf",sep = "/"),
                    group_info = group_info_without_crab)


draw_pheatmap_batch(Phylum_famiy_genus_norm[[3]], 
                    otu_combined_count_sub_result[["Genus"]][[3]], 
                    pdfName = paste(outdir, "#6-7-Differential-Genus-crab.pdf",sep = "/"),
                    group_info = group_info_with_crab)

draw_pheatmap_batch(Phylum_famiy_genus_norm[[3]], 
                    otu_combined_count_sub_result[["Genus"]][[4]], 
                    pdfName = paste(outdir, "#6-7-Differential-Genus-without.pdf",sep = "/"),
                    group_info = group_info_without_crab)


draw_pheatmap_batch <- function(table, otulist, pdfName, group_info){
  
    ann_colors = list(
      Group = c("H" = "cornflowerblue", "C" = "darkgoldenrod1"),
      Class = c("1" = brewer.pal(n = 10, name = "Paired")[1], "2" = brewer.pal(n = 10, name = "Paired")[2], "3" = brewer.pal(n = 10, name = "Paired")[3], "4" = brewer.pal(n = 10, name = "Paired")[4], "5" = brewer.pal(n = 10, name = "Paired")[5], "9" = brewer.pal(n = 10, name = "Paired")[6], "10" = brewer.pal(n = 10, name = "Paired")[7], "12" = brewer.pal(n = 10, name = "Paired")[9], "13" = brewer.pal(n = 11, name = "Paired")[10])
    )
  
  
  table$OTUID[1] <- "NA"
  intersect_OTUs <- data.table(OTUID = otulist) %>%
    left_join(., table) %>%
    column_to_rownames(., var = "OTUID") %>%
    .[rownames(group_info)]
  #return(intersect_OTUs)
  
  # group info order 
  
  group_info <- group_info %>%
    add_column(sample = rownames(group_info)) %>%
    left_join(data.table(sample = colnames(intersect_OTUs)), .) %>%
    column_to_rownames(., var = "sample")
  
  group_info <- group_info[order(group_info$Group),]
  
  group_info$Class <- as.character(group_info$Class)
  
  # order the intersect_OTUs 
  intersect_OTUs <- intersect_OTUs[rownames(group_info)]
  
  
  # draw heatmap 
  pdf(pdfName,width = 10, height = 10)
  pheatmap::pheatmap(log(intersect_OTUs+0.001), 
                     annotation_col = group_info,
                     annotation_colors = ann_colors,
                     show_rownames = TRUE,
                     fontsize_col = 4,
                     cluster_cols = FALSE)
  
  dev.off()
}



# extract sig list 

# Phylum 
Phylum_list1 <- select(otu_combined_count_sub_result[["Phylum"]][[1]], OTU, adjPvalues, Group_Impact) %>%
  add_column(Group = "crab")
Phylum_list2 <- select(otu_combined_count_sub_result[["Phylum"]][[2]], OTU, adjPvalues, Group_Impact) %>%
  add_column(Group = "no-crab")


phylum_list <- rbind(Phylum_list1, Phylum_list2)

# Family 
Family_list1 <-select(otu_combined_count_sub_result[["Family"]][[1]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Group = "crab")
Family_list2 <- select(otu_combined_count_sub_result[["Family"]][[2]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Group = "no-crab")


Family_list <- rbind(Family_list1, Family_list2)

# Genus 
Genus_list1 <- select(otu_combined_count_sub_result[["Genus"]][[1]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Group = "crab")
Genus_list2 <- select(otu_combined_count_sub_result[["Genus"]][[2]], OTU, adjPvalues, Group_Impact)%>%
  add_column(Group = "no-crab")


Genus_list <- rbind(Genus_list1, Genus_list2)







## update the figure 2022,07 ####################################################################################################

color_label <- c( "#f5f50c", "#09ba09")
H_label <- "HE+SH+CR"
control_label <- "SH+CR"

draw_pheatmap_batch_update <- function(otutable = NA, otulist = NA, pdfName = NA, width_seletive = NA, group_info = NA, crab_tag = NA){
  
  otutable$OTUID[1] <- "NA"
  intersect_OTUs <- data.table(OTUID = otulist) %>%
    left_join(., otutable) %>%
    column_to_rownames(., var = "OTUID") %>%
    .[rownames(group_info)]
  #return(intersect_OTUs)
  
  # group info order 
  
  group_info <- group_info %>%
    add_column(sample = rownames(group_info)) %>%
    left_join(data.table(sample = colnames(intersect_OTUs)), .) %>%
    column_to_rownames(., var = "sample")
  
  group_info_order <- group_info[order(group_info$Group),] %>% select(-Class)
  if (crab_tag == "yes"){
    group_info_order$Group <- ifelse(group_info_order$Group == "H", "HE+SH+CR", "SH+CR")
    ann_colors = list(
      Group = c("HE+SH+CR" = "#09ba09", "SH+CR" = "#f5f50c"))
  }
  
  if (crab_tag == "no"){
    group_info_order$Group <- ifelse(group_info_order$Group == "H", "HE+SH", "SH")
    ann_colors = list(
      Group = c("HE+SH" = "#09ba09", "SH" = "#f5f50c"))
  }

  # order the intersect_OTUs 
  intersect_OTUs <- intersect_OTUs[rownames(group_info_order)]
  rownames(intersect_OTUs) <- str_remove_all(rownames(intersect_OTUs), "p__") %>%
    str_remove_all(., "f__") %>%
    str_remove_all(., "g__")
  
  
  newnames <- lapply(
    rownames(intersect_OTUs),
    function(x) bquote(italic(.(x))))
  
  intersect_OTUs_scale <- t(scale(t(intersect_OTUs)))
  # default scale is column scale, scale is generic function whose default method centers and/or scales the columns of a numeric matrix.
  # so here we used z-score to do normalization, then replace the extreme value with the same value and draw plot 
  intersect_OTUs_scale[intersect_OTUs_scale > 3.5] <- 3.5
  
  p2 <- as.matrix(intersect_OTUs_scale) %>% as.numeric() %>% data.frame(values = .) %>% ggplot(., aes(values)) + geom_density(bw = "SJ")
  # for crab phylum, it is very even : -2 - 2 
  # draw heatmap 
  print(p2)

  p <- pheatmap::pheatmap(intersect_OTUs_scale, 
                          annotation_col = group_info_order,
                          annotation_colors = ann_colors,
                          show_rownames = TRUE,
                          cluster_cols = FALSE,
                          # scale = "row",
                          fontsize_col = 5,
                          fontsize_row = 8,
                          cellheight=12, cellwidth = 6,
                          labels_row = as.expression(newnames))
  
  pdf(pdfName, height = 10, width = width_seletive)
  print(p)
  dev.off()
  result <- list(p, p2)
  return(result)
}


# draw heatmap 

outdir_new <- "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/"
phylum_crab_heatmap <- draw_pheatmap_batch_update(otutable = Phylum_famiy_genus_norm[[1]], 
                           otulist  = otu_combined_count_sub_result[["Phylum"]][[3]], 
                           pdfName  = paste(outdir_new, "#6-7-Differential-Phylum-crab-update.pdf",sep = "/"),
                           width_seletive = 9,
                           group_info = group_info_with_crab,
                           crab_tag = "yes")
# -2 to 2 
phylum_nocrab_heatmap <- draw_pheatmap_batch_update(otutable = Phylum_famiy_genus_norm[[1]], 
                    otulist = otu_combined_count_sub_result[["Phylum"]][[4]], 
                    pdfName  = paste(outdir_new, "#6-7-Differential-Phylum-without-crab-update.pdf",sep = "/"),
                    width_seletive = 9,
                    group_info = group_info_without_crab,
                    crab_tag = "no")
#phylum_nocrab_heatmap[[2]]
# -2 to 2 


family_crab_heatmap <- draw_pheatmap_batch_update(otutable = Phylum_famiy_genus_norm[[2]], 
                    otulist = otu_combined_count_sub_result[["Family"]][[3]], 
                    pdfName = paste(outdir_new, "#6-7-Differential-Family-crab-update.pdf",sep = "/"),
                    width_seletive = 9,
                    group_info = group_info_with_crab,
                    crab_tag = "yes")
family_crab_heatmap[[2]]
#-2 to approaching 3 

family_nocrab_heatmap <- draw_pheatmap_batch_update(otutable = Phylum_famiy_genus_norm[[2]], 
                    otulist = otu_combined_count_sub_result[["Family"]][[4]], 
                    pdfName = paste(outdir_new, "#6-7-Differential-Family-without-crab-update.pdf",sep = "/"),
                    width_seletive = 9,
                    group_info = group_info_without_crab,
                    crab_tag = "no")
family_nocrab_heatmap[[2]]
# -2 to 2 

genus_crab_heatmap <- draw_pheatmap_batch_update(otutable = Phylum_famiy_genus_norm[[3]], 
                    otulist = otu_combined_count_sub_result[["Genus"]][[3]], 
                    pdfName = paste(outdir_new, "#6-7-Differential-Genus-crab-update.pdf",sep = "/"),
                    width_seletive = 9,
                    group_info = group_info_with_crab,
                    crab_tag = "yes")
genus_crab_heatmap[[2]]

# -2 to 2 
genus_nocrab_heatmap <- draw_pheatmap_batch_update(otutable = Phylum_famiy_genus_norm[[3]], 
                    otulist = otu_combined_count_sub_result[["Genus"]][[4]], 
                    pdfName = paste(outdir_new, "#6-7-Differential-Genus-without-crab-update.pdf",sep = "/"),
                    width_seletive = 9,
                    group_info = group_info_without_crab,
                    crab_tag = "no")

# -2 to ~ 3 

save.image(paste0(outdir_new, "14.metagenomseq_average_crab_phylum_07_2022_figuremodification.RData"))


color1 <- "#53ed53"
color2 <- "#09ba09"

############### Venn plot modification ############################################################ 2022.07
draw_venn_plot_update_color <- function(crab = NA, nocrab = NA, PngName = NA, outdir = NA, main = NA){
  
  p <- venn.diagram(
    x = list(crab, nocrab),
    category.names = c("HE+SH+CR vs SH+CR" , "HE+SH vs SH"),
    filename=NULL, 
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5,0.8),
    
    # Output features
    imagetype="svg" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 5,
    #lty = 'blank',
    fill = c('#069106', '#09ba09'),
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-5, 5),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    margin = 0.5
  )
  
  
  
  ggsave(p, file=paste(outdir,PngName,sep = "/"), device = "svg", width = 6)
  return(p)
  
}


venn_plot <- list()
for (i in names(otu_combined_count_sub)){
  venn_plot[[i]] <- draw_venn_plot_update_color(crab = otu_combined_count_sub_result[[i]][[3]], 
                                                nocrab = otu_combined_count_sub_result[[i]][[4]], 
                                                outdir = outdir_new, PngName = paste("#6-5-Venn-", i, "-update.svg",sep = ""), 
                                                main = i)
}
save.image(paste0(outdir_new, "14.metagenomseq_average_crab_phylum_07_2022_figuremodification.RData"))

# the figures are open in inkscape and slightly modify the potion of main and save it to pdf. 

save.image("14.metagenomseq_average_crab_phylum.RData")



# venn plot modification again Nov. 19th. 2022 
############################################################################################################
outdir_new <- "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/"
load(paste0(outdir_new, "14.metagenomseq_average_crab_phylum_07_2022_figuremodification.RData"))

library(VennDiagram)
library(tidyverse)

draw_venn_plot_update_color_with_legend <- function(crab = NA, nocrab = NA, PngName = NA, outdir = NA, main = NA, pos = 0.65, color = NA, do.lines = FALSE){
  
  p <- venn.diagram(
    x = list(crab, nocrab),
    category.names = c("" , ""),
    #category.names = c("HE+SH+CR vs SH+CR" , "HE+SH vs SH"),
    filename=NULL, 
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5, pos),
    
    # Output features
    #imagetype="svg" ,
    #height = 480 , 
    #width = 480 , 
    #resolution = 300,
    #compression = "lzw",
    
    # Circles
    lwd = 5,
    #lty = 'blank',

    #fill = c('#069106', '#99d8c9'), # '#09ba09' #e5f5f9
    fill = color,
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "Helvetica",
    
    # Set names
    
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-5, 5),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "Helvetica",
    margin = 0.5#,
    
    # outline https://rdrr.io/cran/VennDiagram/man/venn.diagram.html 
    
    #lty = 'blank'
  )
  
  
  #cols <- c('#069106', '#99d8c9')
  cols <- color 
  lg <- legendGrob(labels = c("HE+SH+CR vs SH+CR" , "HE+SH vs SH"), 
                   pch = rep(15, length(c("HE+SH+CR vs SH+CR" , "HE+SH vs SH"))),
                   gp=gpar(col = cols, fill = "gray"),
                   byrow = TRUE, do.lines = do.lines)
  
  g <- gTree(children = gList(p))
  final_p <- gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4, 1.2))
  
  
  ggsave(final_p, file = paste(outdir, PngName, sep = "/"), device = "pdf", width = 9)  
  
  return(final_p)
  
}


venn_plot <- list()
for (i in names(otu_combined_count_sub)){
  if (i != "Genus"){
    pos <- 0.67
  }else{
    pos <- 0.65
  }
  venn_plot[[i]] <- draw_venn_plot_update_color_with_legend(crab = otu_combined_count_sub_result[[i]][[3]], 
                                                nocrab = otu_combined_count_sub_result[[i]][[4]], 
                                                outdir = outdir_new, 
                                                PngName = paste("#6-5-Venn-", i, "-updatev2.pdf",sep = ""), 
                                                main = i,
                                                pos = pos,
                                                color = c('#069106', '#99d8c9'),
                                                do.lines = FALSE)
}


# remotes::install_github("trevorld/gridpattern")



# version 2 


draw_venn_plot_update_color_with_legend2 <- function(crab = NA, nocrab = NA, PngName = NA, outdir = NA, main = NA, pos = 0.65, color = NA, do.lines = FALSE){
  
  p <- venn.diagram(
    x = list(crab, nocrab),
    category.names = c("" , ""),
    #category.names = c("HE+SH+CR vs SH+CR" , "HE+SH vs SH"),
    filename=NULL, 
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5, pos),
    
    # Output features
    #imagetype="svg" ,
    #height = 480 , 
    #width = 480 , 
    #resolution = 300,
    #compression = "lzw",
    
    # Circles
    lwd = 5,
    #lty = 'blank',
    
    #fill = c('#069106', '#99d8c9'), # '#09ba09' #e5f5f9
    fill = color,
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "Helvetica",
    
    # Set names
    
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    #cat.pos = c(-5, 5),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "Helvetica",
    margin = 0.5#,
    
    # outline https://rdrr.io/cran/VennDiagram/man/venn.diagram.html 
    
    #lty = 'blank'
  )
  
  
  #cols <- c('#069106', '#99d8c9')
  cols <- color 
  lg <- legendGrob(labels = c("HE+SH+CR vs SH+CR" , "HE+SH vs SH"), 
                   pch = rep(15, length(c("HE+SH+CR vs SH+CR" , "HE+SH vs SH"))),
                   gp=gpar(col = cols, fill = cols, lwd = 0),
                   byrow = TRUE)
  
  g <- gTree(children = gList(p))
  final_p <- gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4, 1.2))
  
  
  ggsave(final_p, file = paste(outdir, PngName, sep = "/"), device = "pdf", width = 9)  
  
  return(final_p)
  
}



venn_plot_white <- list()
for (i in names(otu_combined_count_sub)){
  if (i != "Genus"){
    pos <- 0.67
  }else{
    pos <- 0.65
  }
  venn_plot_white[[i]] <- draw_venn_plot_update_color_with_legend2(crab = otu_combined_count_sub_result[[i]][[3]], 
                                                            nocrab = otu_combined_count_sub_result[[i]][[4]], 
                                                            outdir = outdir_new, 
                                                            PngName = paste("#6-5-Venn-", i, "-updatev3.pdf",sep = ""), 
                                                            main = i,
                                                            pos = pos,
                                                            color = c("#FFFFFF", "#FFFFFF"))
}
# manuall fill the pattern https://www.tutorviacomputer.com/inkscape/patterns/ 
# https://www.youtube.com/watch?v=YvmjNsS9FzA&t=75s 


#g$children$GRID.polygon.1310$x

#Diag <- data.frame(
#  x = c(1,1,1.45,1.45), # 1st 2 values dictate starting point of line. 2nd 2 dictate width. Each whole = one background grid
#  y = c(0,0,20,20))

#g + geom_path(data=Diag, aes(x=x, y=y),colour = "black")

outdir_new <- "/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/"
save.image(paste0(outdir_new, "14.metagenomseq_average_crab_phylum_11_2022_figuremodification.RData"))

load(paste0(outdir_new, "14.metagenomseq_average_crab_phylum_11_2022_figuremodification.RData"))

combine_plot <- ggpubr::ggarrange(venn_plot[[1]], venn_plot[[2]], venn_plot[[3]], ncol = 3)

ggsave(paste0(outdir_new, "/#6-5-Venn-updatev3-combine.pdf"), combine_plot, width = 28)

save.image(paste0(outdir_new, "14.metagenomseq_average_crab_phylum_11_2022_figuremodification.RData"))

# update 2023.09.07
############################################################################################################
# change the color of the heatmap. 

load("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/14.metagenomseq_average_crab_phylum_11_2022_figuremodification.RData")

# draw heamap with different colors in the heatmap 
# 
draw_pheatmap_batch_updatev2 <- function(otutable = NA, otulist = NA, pdfName = NA, width_seletive = NA, group_info = NA, crab_tag = NA){
  
  otutable$OTUID[1] <- "NA"
  intersect_OTUs <- data.table(OTUID = otulist) %>%
    left_join(., otutable) %>%
    column_to_rownames(., var = "OTUID") %>%
    .[rownames(group_info)]
  #return(intersect_OTUs)
  
  # group info order 
  
  group_info <- group_info %>%
    add_column(sample = rownames(group_info)) %>%
    left_join(data.table(sample = colnames(intersect_OTUs)), .) %>%
    column_to_rownames(., var = "sample")
  
  group_info_order <- group_info[order(group_info$Group),] %>% select(-Class)
  if (crab_tag == "yes"){
    group_info_order$Group <- ifelse(group_info_order$Group == "H", "HE+SH+CR", "SH+CR")
    ann_colors = list(
      Group = c("HE+SH+CR" = "#09ba09", "SH+CR" = "#f5f50c"))
  }
  # c('#069106', '#99d8c9')
  if (crab_tag == "no"){
    group_info_order$Group <- ifelse(group_info_order$Group == "H", "HE+SH", "SH")
    ann_colors = list(
      Group = c("HE+SH" = "#99d8c9", "SH" = "#fdfd81"))
  }
  
  # order the intersect_OTUs 
  intersect_OTUs <- intersect_OTUs[rownames(group_info_order)]
  rownames(intersect_OTUs) <- str_remove_all(rownames(intersect_OTUs), "p__") %>%
    str_remove_all(., "f__") %>%
    str_remove_all(., "g__")
  
  
  newnames <- lapply(
    rownames(intersect_OTUs),
    function(x) bquote(italic(.(x))))
  
  intersect_OTUs_scale <- t(scale(t(intersect_OTUs)))
  # default scale is column scale, scale is generic function whose default method centers and/or scales the columns of a numeric matrix.
  # so here we used z-score to do normalization, then replace the extreme value with the same value and draw plot 
  intersect_OTUs_scale[intersect_OTUs_scale > 3.5] <- 3.5
  
  p2 <- as.matrix(intersect_OTUs_scale) %>% as.numeric() %>% data.frame(values = .) %>% ggplot(., aes(values)) + geom_density(bw = "SJ")
  # for crab phylum, it is very even : -2 - 2 
  # draw heatmap 
  print(p2)
  
  p <- pheatmap::pheatmap(intersect_OTUs_scale, 
                          annotation_col = group_info_order,
                          annotation_colors = ann_colors,
                          show_rownames = TRUE,
                          cluster_cols = FALSE,
                          # scale = "row",
                          fontsize_col = 5,
                          fontsize_row = 8,
                          cellheight=12, cellwidth = 6,
                          labels_row = as.expression(newnames))
  
  pdf(pdfName, height = 10, width = width_seletive)
  print(p)
  dev.off()
  result <- list(p, p2)
  return(result)
}


# draw heatmap 

outdir_new <- "/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/"
phylum_crab_heatmapv2 <- draw_pheatmap_batch_updatev2(otutable = Phylum_famiy_genus_norm[[1]], 
                                                  otulist  = otu_combined_count_sub_result[["Phylum"]][[3]], 
                                                  pdfName  = paste(outdir_new, "#6-7-Differential-Phylum-crab-updatev2.pdf",sep = "/"),
                                                  width_seletive = 9,
                                                  group_info = group_info_with_crab,
                                                  crab_tag = "yes")
# -2 to 2 
phylum_nocrab_heatmapv2 <- draw_pheatmap_batch_updatev2(otutable = Phylum_famiy_genus_norm[[1]], 
                                                    otulist = otu_combined_count_sub_result[["Phylum"]][[4]], 
                                                    pdfName  = paste(outdir_new, "#6-7-Differential-Phylum-without-crab-updatev2.pdf",sep = "/"),
                                                    width_seletive = 9,
                                                    group_info = group_info_without_crab,
                                                    crab_tag = "no")
#phylum_nocrab_heatmapv2[[2]]
# -2 to 2 


family_crab_heatmapv2 <- draw_pheatmap_batch_updatev2(otutable = Phylum_famiy_genus_norm[[2]], 
                                                  otulist = otu_combined_count_sub_result[["Family"]][[3]], 
                                                  pdfName = paste(outdir_new, "#6-7-Differential-Family-crab-updatev2.pdf",sep = "/"),
                                                  width_seletive = 9,
                                                  group_info = group_info_with_crab,
                                                  crab_tag = "yes")
family_crab_heatmap[[2]]
#-2 to approaching 3 

family_nocrab_heatmapv2 <- draw_pheatmap_batch_updatev2(otutable = Phylum_famiy_genus_norm[[2]], 
                                                    otulist = otu_combined_count_sub_result[["Family"]][[4]], 
                                                    pdfName = paste(outdir_new, "#6-7-Differential-Family-without-crab-updatev2.pdf",sep = "/"),
                                                    width_seletive = 9,
                                                    group_info = group_info_without_crab,
                                                    crab_tag = "no")
family_nocrab_heatmapv2[[2]]
# -2 to 2 

genus_crab_heatmapv2 <- draw_pheatmap_batch_updatev2(otutable = Phylum_famiy_genus_norm[[3]], 
                                                 otulist = otu_combined_count_sub_result[["Genus"]][[3]], 
                                                 pdfName = paste(outdir_new, "#6-7-Differential-Genus-crab-updatev2.pdf",sep = "/"),
                                                 width_seletive = 9,
                                                 group_info = group_info_with_crab,
                                                 crab_tag = "yes")
genus_crab_heatmap[[2]]

# -2 to 2 
genus_nocrab_heatmapv2 <- draw_pheatmap_batch_updatev2(otutable = Phylum_famiy_genus_norm[[3]], 
                                                   otulist = otu_combined_count_sub_result[["Genus"]][[4]], 
                                                   pdfName = paste(outdir_new, "#6-7-Differential-Genus-without-crab-updatev2.pdf",sep = "/"),
                                                   width_seletive = 9,
                                                   group_info = group_info_without_crab,
                                                   crab_tag = "no")

# -2 to ~ 3 

save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/14.metagenomseq_average_crab_phylum_09_2023_figuremodification.RData")


