# functional comparison 
library(readr)
library("magrittr")
library(tidyr)
library(dplyr)
library(rstatix)
library(compositions)
library(textshape)
library(pheatmap)
library("tibble")
library(viridis)
library(stringr)
library(fossil)
library(data.table)
library(vegan)
library(ggpubr)
library(ggbeeswarm)
library(ggplot2)

pathway_input <- read_tsv("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv")


input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/9_diversity/"
load(paste0(input_dir, "9_alpha_diversity_analysis_HE_series_0909.RData"))

# metatable 

alpha_table_add_sel # 
#################################################################################################


#                       Average the counts                                                      #


#################################################################################################

# average OTU count table and do normalization again
HE_series_average_path <- pathway_input %>% 
  select(pathway, metadata_combine_all$SampleID) %>% 
  pivot_longer(., -pathway, names_to = "SampleID", values_to = "count") %>%
  left_join(., metadata_combine_all, by = "SampleID") %>%
  group_by(sample_id_merged, pathway) %>%
  summarize(mean_count = round(mean(count, na.rm = TRUE))) %>%
  pivot_wider(names_from = "sample_id_merged", values_from = "mean_count")  #65 samples 

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/"
write_tsv(HE_series_average_path, paste0(output_dir, "HE_series_pathway_average_count_table.tsv"))

#################################################################################################


#		relative abundance                                                                    			#


#################################################################################################


HE_series_average_path_sel <- HE_series_average_path %>% 
  select(pathway, alpha_table_add_sel$sample_id_merged) #63 samples 
HE_series_average_path_sel_long <- HE_series_average_path_sel %>% 
  pivot_longer(!pathway) 


HE_series_average_path_sel_sum <- HE_series_average_path_sel %>% 
  pivot_longer(!pathway) %>% 
  group_by(name) %>% 
  summarise(sum_value = sum(value))
  
HE_series_average_path_sel_relab <- HE_series_average_path_sel_long %>% 
  left_join(., HE_series_average_path_sel_sum, by = "name") %>% 
  mutate(relab = value/sum_value)


#################################################################################################


#	clr transformed data                                                                    			#


#################################################################################################

HE_series_average_path_sel_CLR <- HE_series_average_path_sel %>% 
  column_to_rownames(., var = "pathway") %>% 
  clr(.)




#################################################################################################


#		sig.pathways compared each comparisons - relative abundance                            			#


#################################################################################################

# wilcox test 
# compare fresh vs 7-19d 
HE_series_average_path_sel_relab_group <- HE_series_average_path_sel_relab %>% 
  left_join(., alpha_table_add_sel, by = c("name" = "sample_id_merged"))



# first picked the ones that are different among all the groups 

HE_series_average_path_sel_relab_group_statistic <- HE_series_average_path_sel_relab_group %>% 
  group_by(pathway) %>%
  kruskal_test(relab ~ Group)
  

# second only picked up the sig.ones and performed dunn test 
HE_series_average_path_sel_relab_group_statistic_sig <- HE_series_average_path_sel_relab_group_statistic %>% 
  mutate(q = p.adjust(p, method = "fdr")) %>% 
  filter(p < 0.05)

HE_series_average_path_sel_relab_group_statistic_q_sig <- HE_series_average_path_sel_relab_group_statistic_sig %>% 
  filter(q < 0.05)



HE_series_average_path_sel_relab_group_statistic_pairwise <- HE_series_average_path_sel_relab_group %>% 
  dplyr::filter(pathway %in% HE_series_average_path_sel_relab_group_statistic_q_sig$pathway) 


#################################################################################################


#		sig.pathways compared each comparisons - clr                                          			#


#################################################################################################



HE_series_average_path_sel_CLR_df <- HE_series_average_path_sel_CLR %>% 
  as.data.frame()

HE_series_average_path_sel_CLR_df$pathway <- rownames(HE_series_average_path_sel_CLR_df)

HE_series_average_path_sel_CLR_df_long <- HE_series_average_path_sel_CLR_df %>% 
  pivot_longer(!pathway) %>% 
  left_join(., alpha_table_add_sel, by = c("name" = "sample_id_merged"))

comparison_table <- combn(HE_series_average_path_sel_CLR_df_long$Group %>% unique(), 2)

i <- 1 
clr_comparison_result <- list()
clr_comparison_result_effectsize <- list()

for (i in 1:ncol(comparison_table)){
  print(i)
  compar1 <- comparison_table[1,i] %>% as.character()
  compar2 <- comparison_table[2,i] %>% as.character()
  compar_final <- paste(compar1, compar2, sep = "_")  
  
  n_sample <- length(HE_series_average_path_sel_CLR_df_long %>% 
                       filter(Group %in% c(compar1, compar2)) %>% 
                       .$name %>% unique())
  HE_series_average_path_sel_CLR_df_long_sel <- HE_series_average_path_sel_CLR_df_long %>% 
    filter(Group %in% c(compar1, compar2)) %>% 
    filter(value > 0) %>% 
    group_by(pathway) %>% 
    summarise(num_sample = n()) %>% 
    filter(num_sample > round(n_sample * 0.2))
    
  
  clr_comparison_result[[compar_final]] <- HE_series_average_path_sel_CLR_df_long %>% 
    filter(Group %in% c(compar1, compar2)) %>% 
    filter(pathway %in% HE_series_average_path_sel_CLR_df_long_sel$pathway) %>% 
    group_by(pathway) %>% 
    wilcox_test(value ~ Group, p.adjust.method = "none") %>% 
    adjust_pvalue(method = "fdr")

  
  
}

# calculate the wilcox effect size in rstatix 

i <- 1 
clr_comparison_result_effectsize <- list()

for (i in 1:ncol(comparison_table)){
  print(i)
  compar1 <- comparison_table[1,i] %>% as.character()
  compar2 <- comparison_table[2,i] %>% as.character()
  compar_final <- paste(compar1, compar2, sep = "_")  
  
  n_sample <- length(HE_series_average_path_sel_CLR_df_long %>% 
                       filter(Group %in% c(compar1, compar2)) %>% 
                       .$name %>% unique())
  HE_series_average_path_sel_CLR_df_long_sel <- HE_series_average_path_sel_CLR_df_long %>% 
    filter(Group %in% c(compar1, compar2)) %>% 
    filter(value > 0) %>% 
    group_by(pathway) %>% 
    summarise(num_sample = n()) %>% 
    filter(num_sample > round(n_sample * 0.2))
  
  clr_comparison_result_effectsize[[compar_final]] <- HE_series_average_path_sel_CLR_df_long %>% 
    filter(Group %in% c(compar1, compar2)) %>% 
    filter(pathway %in% HE_series_average_path_sel_CLR_df_long_sel$pathway) %>% 
    group_by(pathway) %>% 
    wilcox_effsize(value ~ Group, p.adjust.method = "none") %>% 
    adjust_pvalue(method = "fdr")
  
}




# heatmap for HE + HE3 
################################################################################################

#      phetamap for the OTUs                                                                   # 

#################################################################################################


# pathway and description 
pathway_description <- read_tsv("pathways_out/path_abun_unstrat_descrip.tsv")
pathway_description_sel <- pathway_description %>% select(pathway, description)


draw_pheatmap_update <- function(table = NA, otulist = NA, pdfName = NA, ann_color = NA){
  
  library(tidyr)
  intersect_OTUs <- otulist %>%
    filter(p < 0.05)
  tsvname <- str_replace(pdfName, "pdf", "tsv")
  
  write_tsv(intersect_OTUs, tsvname)
  
  group_info <- alpha_table_add_sel %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    select(sample_id_merged, Group) %>% 
    column_to_rownames(., var = "sample_id_merged")
  
  annotation_row <- intersect_OTUs %>% 
    select(pathway, p.adj) %>% 
    left_join(., pathway_description_sel) %>% 
    select(-pathway) %>% 
    rename(pathway = description) %>% 
    column_to_rownames(., var = "pathway") %>% 
    mutate(p.adj = case_when(p.adj >= 0.2 ~ "q>=0.2", p.adj < 0.05 ~ "q<0.05", p.adj < 0.1 ~ "q<0.1", p.adj < 0.2 ~ "q<0.2")) %>% 
    select(p.adj)
    
  table_wide <- table %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    filter(pathway %in% intersect_OTUs$pathway) %>% 
    select(pathway, name, value) %>% 
    pivot_wider(., names_from = name, values_from = value) %>% 
    left_join(., pathway_description_sel) %>% 
    select(-pathway) %>% 
    rename(pathway = description) %>% 
    column_to_rownames(., var = "pathway")
  
  
  
  # draw heatmap 
  
  pdf(pdfName,width = 20, height = 20)
  pheatmap::pheatmap(table_wide, 
                     annotation_col = group_info,
                     annotation_colors = ann_color,
                     annotation_row = annotation_row,
                     show_rownames = TRUE,
                     cluster_cols = FALSE,
                     fontsize_col = 5,
                     fontsize_row = 8,
                     cellheight=12, cellwidth = 6)
  
  dev.off()
}

color <- c("#079F07","#068406","#056A05","#034F03", "#023502") #011A01
color <- viridis::viridis(10)[c(4,6,8,9,10)]
#  "#31688EFF" "#1F9E89FF" "#6DCD59FF" "#B4DE2CFF" "#FDE725FF"
labels = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH")


# HE+SH vs HE3+SH
draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE+SH_HE3+SH"]], 
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE+SH_HE3+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE3+SH" = "#6DCD59FF")))

# HE+SH vs HE2+SH 
draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE+SH_HE2+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE+SH_HE2+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE2+SH" = "#1F9E89FF")))

# HE+SH vs HE4+SH 
draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE+SH_HE4+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE+SH_HE4+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE4+SH" = "#B4DE2CFF")))

# HE+SH vs HE5+SH 

draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE+SH_HE5+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE5+SH" = "#FDE725FF")))


# HE2+SH vs HE3+SH 
draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE3+SH_HE2+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE2+SH_HE3+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE3+SH" = "#6DCD59FF")))

# HE2+SH vs HE4+SH 

draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE2+SH_HE4+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE2+SH_HE4+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE4+SH" = "#B4DE2CFF")))


# HE2+SH vs HE5+SH 

draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE2+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE2+SH_HE5+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE5+SH" = "#FDE725FF")))


# HE3+SH vs HE4+SH 

draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE3+SH_HE4+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE3+SH_HE4+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE3+SH" = "#6DCD59FF", "HE4+SH" = "#B4DE2CFF")))


# HE3+SH vs HE5+SH 
draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE3+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE3+SH_HE5+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE3+SH" = "#6DCD59FF", "HE5+SH" = "#FDE725FF")))


# HE4+SH vs HE5+SH 
draw_pheatmap_update(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE4+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE4+SH_HE5+SH_differential_pathways_clr.pdf",sep = "/"),
                     ann_color = list(Group = c("HE4+SH" = "#B4DE2CFF", "HE5+SH" = "#FDE725FF")))



# there are too much functional output 

draw_pheatmap_update_q0.05 <- function(table = NA, otulist = NA, pdfName = NA, ann_color = NA){
  
  library(tidyr)
  intersect_OTUs <- otulist %>%
    filter(p.adj < 0.05)
  tsvname <- str_replace(pdfName, "pdf", "tsv")
  
  write_tsv(intersect_OTUs, tsvname)
  
  group_info <- alpha_table_add_sel %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    select(sample_id_merged, Group) %>% 
    column_to_rownames(., var = "sample_id_merged") %>% 
    arrange(Group)
  
  annotation_row <- intersect_OTUs %>% 
    select(pathway, p.adj) %>% 
    left_join(., pathway_description_sel) %>% 
    select(-pathway) %>% 
    rename(pathway = description) %>% 
    column_to_rownames(., var = "pathway") %>% 
    mutate(p.adj = case_when(p.adj >= 0.2 ~ "q>=0.2", p.adj < 0.05 ~ "q<0.05", p.adj < 0.1 ~ "q<0.1", p.adj < 0.2 ~ "q<0.2")) %>% 
    select(p.adj)
  
  table_wide <- table %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    filter(pathway %in% intersect_OTUs$pathway) %>% 
    select(pathway, name, value) %>% 
    pivot_wider(., names_from = name, values_from = value) %>% 
    left_join(., pathway_description_sel) %>% 
    select(-pathway) %>% 
    rename(pathway = description) 
  
  
  # draw heatmap 
  if (nrow(table_wide) > 0){
    
    table_wide <- table_wide %>% 
      select(pathway, rownames(group_info)) %>% 
      column_to_rownames(., var = "pathway")
    
    pdf(pdfName,width = 20, height = 20)
    pheatmap::pheatmap(table_wide, 
                       annotation_col = group_info,
                       annotation_colors = ann_color,
                       annotation_row = annotation_row,
                       show_rownames = TRUE,
                       cluster_cols = FALSE,
                       fontsize_col = 5,
                       fontsize_row = 8,
                       cellheight=12, cellwidth = 6)
    
    dev.off()
  }else{
    print("no result")
  }
 
}



# HE+SH vs HE3+SH: noresult 
draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE+SH_HE3+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE+SH_HE3+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE3+SH" = "#6DCD59FF")))

# HE+SH vs HE2+SH :yes 
draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE+SH_HE2+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE+SH_HE2+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE2+SH" = "#1F9E89FF")))

# HE+SH vs HE4+SH : no result 
draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE+SH_HE4+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE+SH_HE4+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE4+SH" = "#B4DE2CFF")))

# HE+SH vs HE5+SH : yes 

draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE+SH_HE5+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE5+SH" = "#FDE725FF")))


# HE2+SH vs HE3+SH : yes 
draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE3+SH_HE2+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE2+SH_HE3+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE3+SH" = "#6DCD59FF")))

# HE2+SH vs HE4+SH :yes 

draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE2+SH_HE4+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE2+SH_HE4+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE4+SH" = "#B4DE2CFF")))


# HE2+SH vs HE5+SH : yes 

draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE2+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE2+SH_HE5+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE5+SH" = "#FDE725FF")))


# HE3+SH vs HE4+SH : no result 

draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE3+SH_HE4+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE3+SH_HE4+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE3+SH" = "#6DCD59FF", "HE4+SH" = "#B4DE2CFF")))


# HE3+SH vs HE5+SH : yes 
draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE3+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE3+SH_HE5+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE3+SH" = "#6DCD59FF", "HE5+SH" = "#FDE725FF")))


# HE4+SH vs HE5+SH : yes 
draw_pheatmap_update_q0.05(table = HE_series_average_path_sel_CLR_df_long, 
                     otulist = clr_comparison_result[["HE4+SH_HE5+SH"]], 
                     pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                      "#HE4+SH_HE5+SH_differential_pathways_clr_q0.05.pdf",sep = "/"),
                     ann_color = list(Group = c("HE4+SH" = "#B4DE2CFF", "HE5+SH" = "#FDE725FF")))


# I need to order effect size 
library(colorspace)
library(scales)
test_color1 <- heat_hcl(n=9)
show_col(test_color1)
hcl_palettes(plot = TRUE)
#https://colorspace.r-forge.r-project.org/articles/colorspace.html
q4 <- qualitative_hcl(4, palette = "Dark 3")
plot(log(EuStockMarkets), plot.type = "single", col = q4, lwd = 2)

qvalue <- sequential_hcl(4, palette = "Light Grays") # 
show_col(sequential_hcl(4, palette = "Light Grays"))
# "#474747" "#8A8A8A" "#C2C2C2" "#E2E2E2"

hcl.pals("divergingx")
#[1] "ArmyRose" "Earth"    "Fall"     "Geyser"   "TealRose" "Temps"    "PuOr"     "RdBu"     "RdGy"     "PiYG"     "PRGn"     "BrBG"     "RdYlBu"   "RdYlGn"  
#[15] "Spectral" "Zissou 1" "Cividis"  "Roma" 


q100_roma <- divergingx_hcl(100, palette = "Roma") %>% rev
show_col(q100_roma)
q100 <- divergingx_hcl(100, palette = "ArmyRose")
show_col(q100)
q100 <- divergingx_hcl(100, palette = "Earth")
show_col(q100)

q100 <- divergingx_hcl(100, palette = "Fall")
show_col(q100)

q100 <- divergingx_hcl(100, palette = "Geyser")
show_col(q100)

q100_rdbu <- divergingx_hcl(100, palette = "RdBu") %>% rev()
show_col(q100_rdbu)


draw_pheatmap_update_q0.05_effect_size <- function(table = NA, otulist = NA, otu_effectsize_list = NA, pdfName = NA, ann_color = NA, target_color = NA){
  
  library(tidyr)
  intersect_OTUs <- otulist %>%
    left_join(., otu_effectsize_list %>% select(pathway, effsize, magnitude), by = "pathway") %>% 
    filter(p.adj < 0.05) %>% 
    filter(effsize > 0.25) %>% 
    arrange(effsize)
    
  tsvname <- str_replace(pdfName, "pdf", "tsv")
  
  write_tsv(intersect_OTUs, tsvname)
  
  group_info <- alpha_table_add_sel %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    select(sample_id_merged, Group) %>% 
    column_to_rownames(., var = "sample_id_merged") %>% 
    arrange(Group)
  
  annotation_row <- intersect_OTUs %>% 
    select(pathway, p.adj, effsize) %>%
    rename(qvalue = p.adj) %>%  
    left_join(., pathway_description_sel) %>% 
    select(-pathway) %>% 
    rename(pathway = description) %>% 
    column_to_rownames(., var = "pathway") %>% 
    mutate(qvalue = case_when(qvalue >= 0.2 ~ "q>=0.2", qvalue < 0.05 ~ "q<0.05", qvalue < 0.1 ~ "q<0.1", qvalue < 0.2 ~ "q<0.2")) %>% 
    select(qvalue, effsize)
  
  table_wide <- table %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    filter(pathway %in% intersect_OTUs$pathway) %>% 
    select(pathway, name, value) %>% 
    pivot_wider(., names_from = name, values_from = value) %>% 
    left_join(., pathway_description_sel) %>% 
    left_join(tibble(pathway = intersect_OTUs$pathway), .) %>% 
    select(-pathway) %>% 
    rename(pathway = description) 
  
  
  # draw heatmap 
  if (nrow(table_wide) > 0){
    
    table_wide <- table_wide %>% 
      select(pathway, rownames(group_info)) %>% 
      column_to_rownames(., var = "pathway") 

    #ann_color[["qvalue"]] <- c("q>=0.2" = "#E2E2E2", "q<0.2" = "#C2C2C2", "q<0.1" = "#8A8A8A", "q<0.05" = "#474747")
    ann_color[["qvalue"]] <- c("q<0.05" = "#8A8A8A")
    
    pdf(pdfName,width = 20, height = 20)
    pheatmap::pheatmap(table_wide, 
                       annotation_col = group_info,
                       annotation_colors = ann_color,
                       annotation_row = annotation_row,
                       show_rownames = TRUE,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       fontsize_col = 5,
                       fontsize_row = 8,
                       color = target_color,
                       cellheight=12, cellwidth = 6)
    
    dev.off()
  }else{
    print("no result")
  }
  
}


# the ones doesn't have result before, i dont run 

# HE+SH vs HE2+SH :yes 
draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE+SH_HE2+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE+SH_HE2+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE+SH_HE2+SH_differential_pathways_clr_q0.05_effsize.pdf",sep = "/"),
                           ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE2+SH" = "#1F9E89FF")),
                           target_color = q100_rdbu)

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                                       otulist = clr_comparison_result[["HE+SH_HE2+SH"]], 
                                       otu_effectsize_list = clr_comparison_result_effectsize[["HE+SH_HE2+SH"]],
                                       pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                                        "#HE+SH_HE2+SH_differential_pathways_clr_q0.05_effsize_roma.pdf",sep = "/"),
                                       ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE2+SH" = "#1F9E89FF")),
                                       target_color = q100_roma)



# HE+SH vs HE5+SH : yes 

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE+SH_HE5+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE+SH_HE5+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE+SH_HE5+SH_differential_pathways_clr_q0.05_effsize.pdf",sep = "/"),
                           ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE5+SH" = "#FDE725FF")),
                           target_color = q100_rdbu)

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE+SH_HE5+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE+SH_HE5+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE+SH_HE5+SH_differential_pathways_clr_q0.05_effsize_roma.pdf",sep = "/"),
                           ann_color = list(Group = c("HE+SH" = "#31688EFF", "HE5+SH" = "#FDE725FF")),
                           target_color = q100_roma)



# HE2+SH vs HE3+SH : yes 
draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE3+SH_HE2+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE3+SH_HE2+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE2+SH_HE3+SH_differential_pathways_clr_q0.05_effsize.pdf",sep = "/"),
                           ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE3+SH" = "#6DCD59FF")),
                           target_color = q100_rdbu)

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                                       otulist = clr_comparison_result[["HE3+SH_HE2+SH"]], 
                                       otu_effectsize_list = clr_comparison_result_effectsize[["HE3+SH_HE2+SH"]],
                                       pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                                        "#HE2+SH_HE3+SH_differential_pathways_clr_q0.05_effsize_roma.pdf",sep = "/"),
                                       ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE3+SH" = "#6DCD59FF")),
                                       target_color = q100_roma)


# HE2+SH vs HE4+SH :yes 

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE2+SH_HE4+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE2+SH_HE4+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE2+SH_HE4+SH_differential_pathways_clr_q0.05_effsize.pdf",sep = "/"),
                           ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE4+SH" = "#B4DE2CFF")),
                           target_color = q100_rdbu)


draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                                       otulist = clr_comparison_result[["HE2+SH_HE4+SH"]], 
                                       otu_effectsize_list = clr_comparison_result_effectsize[["HE2+SH_HE4+SH"]],
                                       pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                                        "#HE2+SH_HE4+SH_differential_pathways_clr_q0.05_effsize_roma.pdf",sep = "/"),
                                       ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE4+SH" = "#B4DE2CFF")),
                                       target_color = q100_roma)


# HE2+SH vs HE5+SH : yes 

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE2+SH_HE5+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE2+SH_HE5+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE2+SH_HE5+SH_differential_pathways_clr_q0.05_effsize.pdf",sep = "/"),
                           ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE5+SH" = "#FDE725FF")),
                           target_color = q100_rdbu)

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE2+SH_HE5+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE2+SH_HE5+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE2+SH_HE5+SH_differential_pathways_clr_q0.05_effsize_roma.pdf",sep = "/"),
                           ann_color = list(Group = c("HE2+SH" = "#1F9E89FF", "HE5+SH" = "#FDE725FF")),
                           target_color = q100_roma)




# HE3+SH vs HE5+SH : yes 
draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE3+SH_HE5+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE3+SH_HE5+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE3+SH_HE5+SH_differential_pathways_clr_q0.05_effsize.pdf",sep = "/"),
                           ann_color = list(Group = c("HE3+SH" = "#6DCD59FF", "HE5+SH" = "#FDE725FF")),
                           target_color = q100_rdbu)

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                                       otulist = clr_comparison_result[["HE3+SH_HE5+SH"]], 
                                       otu_effectsize_list = clr_comparison_result_effectsize[["HE3+SH_HE5+SH"]],
                                       pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                                        "#HE3+SH_HE5+SH_differential_pathways_clr_q0.05_effsize_roma.pdf",sep = "/"),
                                       ann_color = list(Group = c("HE3+SH" = "#6DCD59FF", "HE5+SH" = "#FDE725FF")),
                                       target_color = q100_roma)


# HE4+SH vs HE5+SH : yes 
draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                           otulist = clr_comparison_result[["HE4+SH_HE5+SH"]], 
                           otu_effectsize_list = clr_comparison_result_effectsize[["HE4+SH_HE5+SH"]],
                           pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                            "#HE4+SH_HE5+SH_differential_pathways_clr_q0.05_effsize.pdf",sep = "/"),
                           ann_color = list(Group = c("HE4+SH" = "#B4DE2CFF", "HE5+SH" = "#FDE725FF")),
                           target_color = q100_rdbu)

draw_pheatmap_update_q0.05_effect_size(table = HE_series_average_path_sel_CLR_df_long, 
                                       otulist = clr_comparison_result[["HE4+SH_HE5+SH"]], 
                                       otu_effectsize_list = clr_comparison_result_effectsize[["HE4+SH_HE5+SH"]],
                                       pdfName  = paste("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/", 
                                                        "#HE4+SH_HE5+SH_differential_pathways_clr_q0.05_effsize_roma.pdf",sep = "/"),
                                       ann_color = list(Group = c("HE4+SH" = "#B4DE2CFF", "HE5+SH" = "#FDE725FF")),
                                       target_color = q100_roma)










# heatmap with cluster version : haven't run 



draw_pheatmap_update_q0.05_effect_size_cluster <- function(table = NA, otulist = NA, otu_effectsize_list = NA, pdfName = NA, ann_color = NA){
  
  library(tidyr)
  intersect_OTUs <- otulist %>%
    left_join(., otu_effectsize_list %>% select(pathway, effsize, magnitude), by = "pathway") %>% 
    filter(p.adj < 0.05) %>% 
    filter(effsize > 0.25) %>% 
    arrange(effsize)
  
  tsvname <- str_replace(pdfName, "pdf", "tsv")
  
  write_tsv(intersect_OTUs, tsvname)
  
  group_info <- alpha_table_add_sel %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    select(sample_id_merged, Group) %>% 
    column_to_rownames(., var = "sample_id_merged") %>% 
    arrange(Group)
  
  annotation_row <- intersect_OTUs %>% 
    select(pathway, p.adj, effsize) %>%
    rename(qvalue = p.adj) %>%  
    left_join(., pathway_description_sel) %>% 
    select(-pathway) %>% 
    rename(pathway = description) %>% 
    column_to_rownames(., var = "pathway") %>% 
    mutate(qvalue = case_when(qvalue >= 0.2 ~ "q>=0.2", qvalue < 0.05 ~ "q<0.05", qvalue < 0.1 ~ "q<0.1", qvalue < 0.2 ~ "q<0.2")) %>% 
    select(qvalue, effsize)
  
  table_wide <- table %>% 
    filter(Group %in% c(unique(otulist$group1), unique(otulist$group2))) %>% 
    filter(pathway %in% intersect_OTUs$pathway) %>% 
    select(pathway, name, value) %>% 
    pivot_wider(., names_from = name, values_from = value) %>% 
    left_join(., pathway_description_sel) %>% 
    #left_join(tibble(pathway = intersect_OTUs$pathway), .) %>% 
    select(-pathway) %>% 
    rename(pathway = description) 
  
  
  # draw heatmap 
  if (nrow(table_wide) > 0){
    
    table_wide <- table_wide %>% 
      select(pathway, rownames(group_info)) %>% 
      column_to_rownames(., var = "pathway") 
    
    ann_color[["qvalue"]] <- c("q>=0.2" = "#E2E2E2", "q<0.2" = "#C2C2C2", "q<0.1" = "#8A8A8A", "q<0.05" = "#474747")
    
    pdf(pdfName,width = 20, height = 20)
    pheatmap::pheatmap(table_wide, 
                       annotation_col = group_info,
                       annotation_colors = ann_color,
                       annotation_row = annotation_row,
                       show_rownames = TRUE,
                       cluster_cols = FALSE,
                       cluster_rows = TRUE,
                       fontsize_col = 5,
                       fontsize_row = 8,
                       color = q100,
                       cellheight=12, cellwidth = 6)
    
    dev.off()
  }else{
    print("no result")
  }
  
}






#################################################################################################


#		beta diversity  of the functional analysis                                            			#


#################################################################################################


sta_beta_result <- function(otu_table = NA, metatable_sel = NA, name = NA){
  set.seed(5)
  distance.bray <- otu_table %>%
    as.matrix(.) %>%
    t(.) %>%
    vegdist(.,method = 'bray')
  
  
  metatable_sel_order <- data.table(sample_id_merged = names(otu_table)) %>%
    left_join(., metatable_sel) 
  
  
  # beta dispersion for group 
  beta_dispersion <- betadisper(distance.bray, group = metatable_sel_order$Group)
  beta_dispersion_pvalue <-  permutest(beta_dispersion)
  
  
  # anosim for group 
  anosim.result <- anosim(distance.bray, metatable_sel_order$Group, permutations = 999) # it will be affected by beta dispersion 
  anosim.r <- anosim.result$statistic
  anosim.p <- anosim.result$signif
  
  # adonis for group and class 
  set.seed(5)
  adonis_group <- adonis(distance.bray ~ Group, data = metatable_sel_order, distance = "bray") # the group should be in the end position
  
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[1]#last postion is group
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[1]#last postion is group
  
  
  
  stat_beta <- list()
  stat_beta[["disp"]] <- beta_dispersion_pvalue
  stat_beta[["anosim"]] <- anosim.result
  stat_beta[["adonisr"]] <- adonis_group
  stat_beta[["adonisp"]] <- adonis_group_p_value
  
  set.seed(1234)
  pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  colnames(points) <- c("x", "y", "z")
  eig <- pcoa$eig
  points <- cbind(points, metatable_sel_order) # correct this one 
  c <- PcoA_plot(points_df = points, eig = eig, statistics = stat_beta) 
  stat_beta[["plot"]] <- c
  
  stat_beta[["dis"]] <- distance.bray
  
  return(stat_beta)
}

#color <- c("#079F07","#068406","#056A05","#034F03", "#023502") #011A01
color <- viridis(10)[c(4,6,8,9,10)]

PcoA_plot <- function(points_df = NA, eig = NA, statistics = NA){
  rsquare_beta_diversity <- round(statistics[["adonisr"]]$aov.tab$R2[1], 3)
  mean_table <- points_df %>% group_by(Group) %>% summarise(x_mean = mean(x), y_mean = mean(y))
  p <- ggscatter(points_df, x = "x", y = "y",
                 color = "Group", palette = color,
                 shape = "Group",
                 ellipse = TRUE, 
                 ellipse.type = "convex"
  ) + 
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    annotate(geom="text", x = 0.2, y = 0.3, label = bquote(P == ~ .(statistics[["adonisp"]]))) + 
    annotate(geom="text", x = 0.2, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_beta_diversity))) +
    geom_point (data = mean_table, aes (x_mean, y_mean), size = 2, shape = "triangle", color = "gray") + 
    geom_segment (x = mean_table$x_mean[1], y = mean_table$y_mean[1], xend = mean_table$x_mean[2], yend = mean_table$y_mean[2], color = "gray") +
    geom_segment (x = mean_table$x_mean[2], y = mean_table$y_mean[2], xend = mean_table$x_mean[3], yend = mean_table$y_mean[3], color = "gray") +
    geom_segment (x = mean_table$x_mean[3], y = mean_table$y_mean[3], xend = mean_table$x_mean[4], yend = mean_table$y_mean[4], color = "gray") +
    geom_segment (x = mean_table$x_mean[4], y = mean_table$y_mean[4], xend = mean_table$x_mean[5], yend = mean_table$y_mean[5], color = "gray") 
  #stat_rows_center(fun.center = "mean", size = 3, shape = "rhombus")
  
  
  
  
  return(p)
}



beta_diversity_path_statistics <- sta_beta_result(otu_table = HE_series_average_path_sel %>% column_to_rownames(., var = "pathway"), metatable_sel = alpha_table_add_sel)
head(beta_diversity_path_statistics)


pdf(paste(output_dir, "HE_series_path_betadiversity.pdf", sep = "/"))
beta_diversity_path_statistics[["plot"]]
dev.off()

# pairwise adonis 


pairwise_adonis <- function(input = NA, meta = NA){
  meta$Group <- as.character(meta$Group)
  comparison <- combn(unique(meta$Group), 2)
  pairwise_adonis_tibble <- comparison %>% t() %>% as.data.frame() %>% .[1:10,] %>% as.tibble() %>%
    add_column(R2 = -1, p_value = -1)
  
  
  for (i in 1:ncol(comparison)){
    print(i)
    meta_sel <- meta %>% filter(Group %in% comparison[,i])
    meta_sel$Group <- meta_sel$Group %>% as.character()
    distance.bray.sel <- input %>%
      select(meta_sel$sample_id_merged) %>%
      add_column(., rowsum = rowSums(.)) %>%
      filter(rowsum > 0) %>% 
      select(-rowsum) %>%
      as.matrix(.) %>%
      t(.) %>%
      vegdist(.,method = 'bray')
    
    set.seed(1234)
    adonis_result_temp <- vegan::adonis(distance.bray.sel ~ Group, data = meta_sel, distance = "bray")
    pairwise_adonis_tibble$p_value[i] <- adonis_result_temp$aov.tab$`Pr(>F)`[1]
    pairwise_adonis_tibble$R2[i] <- adonis_result_temp$aov.tab$R2[1]
  }
  
  return(pairwise_adonis_tibble)
}


metatable_sel_path_order <- data.table(sample_id_merged = names(HE_series_average_path_sel %>% column_to_rownames(., var = "pathway"))) %>%
  left_join(., alpha_table_add_sel) 


input <- HE_series_average_path_sel %>% column_to_rownames(., var = "pathway")
input_t <- input %>%
  as.matrix(.) %>%
  t(.) 

pairwise_table <- pairwise_adonis(input = HE_series_average_path_sel %>% column_to_rownames(., var = "pathway"), 
                                  meta = metatable_sel_path_order)


#> pairwise_table
## A tibble: 10 × 4
#V1     V2         R2 p_value
#<chr>  <chr>   <dbl>   <dbl>
1 HE+SH  HE3+SH 0.0749   0.17 
2 HE+SH  HE2+SH 0.0828   0.132
3 HE+SH  HE4+SH 0.149    0.123
4 HE+SH  HE5+SH 0.0588   0.17 
5 HE3+SH HE2+SH 0.0778   0.119
6 HE3+SH HE4+SH 0.241    0.013
7 HE3+SH HE5+SH 0.121    0.022
8 HE2+SH HE4+SH 0.294    0.007
9 HE2+SH HE5+SH 0.127    0.019
10 HE4+SH HE5+SH 0.0355   0.357

#############################################################################################################
# change the abundance to relab and do the bray curtis calculation again 


#############################################################################################################


sta_beta_result_relab <- function(otu_table = NA, metatable_sel = NA, name = NA){
  set.seed(5)
  distance.bray <- otu_table %>%
    as.matrix(.) %>%
    t(.) %>%
    vegdist(.,method = 'bray')
  
  
  metatable_sel_order <- data.table(sample_id_merged = names(otu_table)) %>%
    left_join(., metatable_sel) 
  
  
  # beta dispersion for group 
  beta_dispersion <- betadisper(distance.bray, group = metatable_sel_order$Group)
  beta_dispersion_pvalue <-  permutest(beta_dispersion)
  
  
  # anosim for group 
  anosim.result <- anosim(distance.bray, metatable_sel_order$Group, permutations = 999) # it will be affected by beta dispersion 
  anosim.r <- anosim.result$statistic
  anosim.p <- anosim.result$signif
  
  # adonis for group and class 
  set.seed(5)
  adonis_group <- adonis(distance.bray ~ Group, data = metatable_sel_order, distance = "bray") # the group should be in the end position
  
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[1]#last postion is group
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[1]#last postion is group
  
  
  
  stat_beta <- list()
  stat_beta[["disp"]] <- beta_dispersion_pvalue
  stat_beta[["anosim"]] <- anosim.result
  stat_beta[["adonisr"]] <- adonis_group
  stat_beta[["adonisp"]] <- adonis_group_p_value
  
  set.seed(1234)
  pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  colnames(points) <- c("x", "y", "z")
  eig <- pcoa$eig
  points <- cbind(points, metatable_sel_order) # correct this one 
  c <- PcoA_plot_relab(points_df = points, eig = eig, statistics = stat_beta) 
  stat_beta[["plot"]] <- c
  
  stat_beta[["dis"]] <- distance.bray
  
  return(stat_beta)
}

#color <- c("#079F07","#068406","#056A05","#034F03", "#023502") #011A01
color <- viridis(10)[c(4,6,8,9,10)]

PcoA_plot_relab <- function(points_df = NA, eig = NA, statistics = NA){
  rsquare_beta_diversity <- round(statistics[["adonisr"]]$aov.tab$R2[1], 3)
  mean_table <- points_df %>% group_by(Group) %>% summarise(x_mean = mean(x), y_mean = mean(y))
  p <- ggscatter(points_df, x = "x", y = "y",
                 color = "Group", palette = color,
                 shape = "Group",
                 ellipse = TRUE, 
                 ellipse.type = "convex"
  ) + 
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    annotate(geom="text", x = 0.13, y = 0.1, label = bquote(P == ~ .(statistics[["adonisp"]]))) + 
    annotate(geom="text", x = 0.13, y = 0.12, label = bquote(R^2 ==  ~ .(rsquare_beta_diversity))) +
    geom_point (data = mean_table, aes (x_mean, y_mean), size = 2, shape = "triangle", color = "gray") + 
    geom_segment (x = mean_table$x_mean[1], y = mean_table$y_mean[1], xend = mean_table$x_mean[2], yend = mean_table$y_mean[2], color = "gray") +
    geom_segment (x = mean_table$x_mean[2], y = mean_table$y_mean[2], xend = mean_table$x_mean[3], yend = mean_table$y_mean[3], color = "gray") +
    geom_segment (x = mean_table$x_mean[3], y = mean_table$y_mean[3], xend = mean_table$x_mean[4], yend = mean_table$y_mean[4], color = "gray") +
    geom_segment (x = mean_table$x_mean[4], y = mean_table$y_mean[4], xend = mean_table$x_mean[5], yend = mean_table$y_mean[5], color = "gray") 
  #stat_rows_center(fun.center = "mean", size = 3, shape = "rhombus")
  
  
  
  
  return(p)
}

## relab transform 
otu_table_path_sum <- HE_series_average_path_sel %>% 
  pivot_longer(!pathway) %>% 
  group_by(name) %>% 
  summarise(sum_ab = sum(value)) %>% 
  ungroup()
otu_table_path <- HE_series_average_path_sel %>% 
  pivot_longer(!pathway) %>% 
  left_join(., otu_table_path_sum) %>% 
  mutate(relab = value/sum_ab) %>% 
  mutate(log10_ab = log10(relab + 1))

otu_table_path_relab <- otu_table_path %>% 
  select(-log10_ab, -value, -sum_ab) %>% 
  pivot_wider(names_from = name, values_from = relab)



beta_diversity_path_statistics_relab <- sta_beta_result_relab(otu_table = otu_table_path_relab %>% column_to_rownames(., var = "pathway"), metatable_sel = alpha_table_add_sel)
head(beta_diversity_path_statistics_relab)


pdf(paste(output_dir, "HE_series_path_betadiversity_relab.pdf", sep = "/"))
beta_diversity_path_statistics_relab[["plot"]]
dev.off()

pairwise_table_relab <- pairwise_adonis(input = otu_table_path_relab %>% column_to_rownames(., var = "pathway"), 
                                  meta = metatable_sel_path_order)
pairwise_table_relab
#<chr>  <chr>   <dbl>   <dbl>
#  1 HE+SH  HE3+SH 0.128    0.012
#2 HE+SH  HE2+SH 0.225    0.001
#3 HE+SH  HE4+SH 0.193    0.046
#4 HE+SH  HE5+SH 0.188    0.001
#5 HE3+SH HE2+SH 0.125    0.007
#6 HE3+SH HE4+SH 0.0598   0.36 
#7 HE3+SH HE5+SH 0.0573   0.113
#8 HE2+SH HE4+SH 0.103    0.085
#9 HE2+SH HE5+SH 0.152    0.001
#10 HE4+SH HE5+SH 0.0916   0.036







#######################################################################################################################

# change the distance to ecludiean distance with log10(X + 1): log transform of the count table                       # 

#######################################################################################################################


euclidean 
## euclidean distance 


sta_beta_result_euclidean <- function(otu_table = NA, metatable_sel = NA, name = NA){
  set.seed(5)
  distance.bray <- otu_table %>%
    as.matrix(.) %>%
    t(.) %>%
    vegdist(.,method = 'euclidean')
  
  
  metatable_sel_order <- data.table(sample_id_merged = names(otu_table)) %>%
    left_join(., metatable_sel) 
  
  
  # beta dispersion for group 
  beta_dispersion <- betadisper(distance.bray, group = metatable_sel_order$Group)
  beta_dispersion_pvalue <-  permutest(beta_dispersion)
  
  
  # anosim for group 
  anosim.result <- anosim(distance.bray, metatable_sel_order$Group, permutations = 999) # it will be affected by beta dispersion 
  anosim.r <- anosim.result$statistic
  anosim.p <- anosim.result$signif
  
  # adonis for group and class 
  set.seed(5)
  adonis_group <- adonis(distance.bray ~ Group, data = metatable_sel_order, distance = "bray") # the group should be in the end position
  
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[1]#last postion is group
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[1]#last postion is group
  
  
  
  stat_beta <- list()
  stat_beta[["disp"]] <- beta_dispersion_pvalue
  stat_beta[["anosim"]] <- anosim.result
  stat_beta[["adonisr"]] <- adonis_group
  stat_beta[["adonisp"]] <- adonis_group_p_value
  
  set.seed(1234)
  pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  colnames(points) <- c("x", "y", "z")
  eig <- pcoa$eig
  points <- cbind(points, metatable_sel_order) # correct this one 
  c <- PcoA_plot_euclidean(points_df = points, eig = eig, statistics = stat_beta) 
  stat_beta[["plot"]] <- c
  
  stat_beta[["dis"]] <- distance.bray
  
  return(stat_beta)
}



PcoA_plot_euclidean <- function(points_df = NA, eig = NA, statistics = NA){
  rsquare_beta_diversity <- round(statistics[["adonisr"]]$aov.tab$R2[1], 3)
  mean_table <- points_df %>% group_by(Group) %>% summarise(x_mean = mean(x), y_mean = mean(y))
  p <- ggscatter(points_df, x = "x", y = "y",
                 color = "Group", palette = color,
                 shape = "Group",
                 ellipse = TRUE, 
                 ellipse.type = "convex"
  ) + 
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    annotate(geom="text", x = 0.010, y = 0.013, label = bquote(P == ~ .(statistics[["adonisp"]]))) + 
    annotate(geom="text", x = 0.010, y = 0.015, label = bquote(R^2 ==  ~ .(rsquare_beta_diversity))) +
    geom_point (data = mean_table, aes (x_mean, y_mean), size = 2, shape = "triangle", color = "gray") + 
    geom_segment (x = mean_table$x_mean[1], y = mean_table$y_mean[1], xend = mean_table$x_mean[2], yend = mean_table$y_mean[2], color = "gray") +
    geom_segment (x = mean_table$x_mean[2], y = mean_table$y_mean[2], xend = mean_table$x_mean[3], yend = mean_table$y_mean[3], color = "gray") +
    geom_segment (x = mean_table$x_mean[3], y = mean_table$y_mean[3], xend = mean_table$x_mean[4], yend = mean_table$y_mean[4], color = "gray") +
    geom_segment (x = mean_table$x_mean[4], y = mean_table$y_mean[4], xend = mean_table$x_mean[5], yend = mean_table$y_mean[5], color = "gray") 
  #stat_rows_center(fun.center = "mean", size = 3, shape = "rhombus")
  
  
  
  
  return(p)
}

otu_table_path_log <- otu_table_path %>% 
  select(-relab, -value, -sum_ab) %>% 
  pivot_wider(names_from = name, values_from = log10_ab)



beta_diversity_path_statistics_log <- sta_beta_result_euclidean(otu_table = otu_table_path_log %>% column_to_rownames(., var = "pathway"), metatable_sel = alpha_table_add_sel)
head(beta_diversity_path_statistics_log)


pdf(paste(output_dir, "HE_series_path_betadiversity_relab_log.pdf", sep = "/"))
beta_diversity_path_statistics_log[["plot"]]
dev.off()



pairwise_adonis_euclidean <- function(input = NA, meta = NA){
  meta$Group <- as.character(meta$Group)
  comparison <- combn(unique(meta$Group), 2)
  pairwise_adonis_tibble <- comparison %>% t() %>% as.data.frame() %>% .[1:10,] %>% as.tibble() %>%
    add_column(R2 = -1, p_value = -1)
  
  
  for (i in 1:ncol(comparison)){
    print(i)
    meta_sel <- meta %>% filter(Group %in% comparison[,i])
    meta_sel$Group <- meta_sel$Group %>% as.character()
    distance.bray.sel <- input %>%
      select(meta_sel$sample_id_merged) %>%
      add_column(., rowsum = rowSums(.)) %>%
      filter(rowsum > 0) %>% 
      select(-rowsum) %>%
      as.matrix(.) %>%
      t(.) %>%
      vegdist(.,method = 'euclidean')
    
    set.seed(1234)
    adonis_result_temp <- vegan::adonis(distance.bray.sel ~ Group, data = meta_sel, distance = "bray")
    pairwise_adonis_tibble$p_value[i] <- adonis_result_temp$aov.tab$`Pr(>F)`[1]
    pairwise_adonis_tibble$R2[i] <- adonis_result_temp$aov.tab$R2[1]
  }
  
  return(pairwise_adonis_tibble)
}


  pairwise_table_log <- pairwise_adonis_euclidean(input = otu_table_path_log %>% column_to_rownames(., var = "pathway"), 
                                                                          meta = metatable_sel_path_order)

  
#  V1     V2         R2 p_value
#  <chr>  <chr>   <dbl>   <dbl>
#    1 HE+SH  HE3+SH 0.131    0.005
#  2 HE+SH  HE2+SH 0.227    0.001
#  3 HE+SH  HE4+SH 0.207    0.033
#  4 HE+SH  HE5+SH 0.181    0.001
#  5 HE3+SH HE2+SH 0.135    0.005
#  6 HE3+SH HE4+SH 0.0683   0.306
#  7 HE3+SH HE5+SH 0.0568   0.112
#  8 HE2+SH HE4+SH 0.0889   0.125
#  9 HE2+SH HE5+SH 0.146    0.001
#  10 HE4+SH HE5+SH 0.0927   0.034
  
  
#######################################################################################################################

# metagenomeSeq transformed ???                     # 

#######################################################################################################################
# so far i just keep using the relab 




#################################################################################################


#		alpha diversity of the functional analysis                                            			#


#################################################################################################



# alpha diversity 
alpha_diversity <- function(otu_table = NA, metatable_sel = NA, name = NA){
  
  otu <- otu_table %>%
    as.matrix(.) %>%
    t(.) %>%
    as.data.frame()
  
  otu_shannon <- vegan::diversity(otu,"shannon")
  
  otu_simpson <- vegan::diversity(otu, "simpson")
  
  otu_chao1 <- apply(otu,1,fossil::chao1)
  
  # make table 
  alpha_df <- data.table(shannon = otu_shannon, simpson = otu_simpson, chao1 = otu_chao1) %>%
    add_column(sample_id_merged = names(otu_shannon)) %>% 
    left_join(., metatable_sel, by = "sample_id_merged")
  
  return(alpha_df)
  
}

metadata_HE_series_unique_pathway <- alpha_table_add_sel %>% select(Group, sample_id_merged) %>% unique() 
alpha_table_path <- alpha_diversity(otu_table = HE_series_average_path_sel %>% column_to_rownames(., var = "pathway"), metatable_sel = metadata_HE_series_unique_pathway)

write_csv(alpha_table_path, paste(output_dir, "HE_series_path_average_alpha_table.csv",sep = "/"))

# alpha diversity figure 

alpha_diversity_statistic_shannon <- rstatix::kruskal_test(alpha_table_path, shannon ~ Group)
alpha_diversity_statistic_simpson <- rstatix::kruskal_test(alpha_table_path, simpson ~ Group)
alpha_diversity_statistic_chao1 <- rstatix::kruskal_test(alpha_table_path, chao1 ~ Group)

alpha_diversity_statistic_ks_all <- rbind(alpha_diversity_statistic_shannon,
                                          alpha_diversity_statistic_simpson,
                                          alpha_diversity_statistic_chao1)

alpha_diversity_statistic_ks_all
#> alpha_diversity_statistic_ks_all
## A tibble: 3 × 6
#.y.         n statistic    df       p method        
#<chr>   <int>     <dbl> <int>   <dbl> <chr>         
#  1 shannon    63      13.8     4 0.00803 Kruskal-Wallis
#2 simpson    63      11.9     4 0.0178  Kruskal-Wallis
#3 chao1      63      16.3     4 0.00269 Kruskal-Wallis



alpha_diversity_dunn_shannon <- rstatix::dunn_test(alpha_table_path, shannon ~ Group, p.adjust.method = "fdr")
alpha_diversity_dunn_simpson <- rstatix::dunn_test(alpha_table_path, simpson ~ Group, p.adjust.method = "fdr")
alpha_diversity_dunn_chao1 <- rstatix::dunn_test(alpha_table_path, chao1 ~ Group, p.adjust.method = "fdr")

alpha_diversity_dunn_all <- rbind(alpha_diversity_dunn_shannon,
                                  alpha_diversity_dunn_simpson,
                                  alpha_diversity_dunn_chao1)

alpha_diversity_dunn_all

write_tsv(alpha_diversity_dunn_all, paste0(output_dir, "/HE_series_path_dunn_test.tsv"))
write_tsv(alpha_diversity_statistic_ks_all, paste0(output_dir, "/HE_series_path_ks_test.tsv"))



#################################################################################################


#		alpha diversity - boxplot                                                             			#


#################################################################################################


color <- c("#079F07","#068406","#056A05","#034F03", "#023502") #011A01
color <- viridis::viridis(10)[c(4,6,8,9,10)]

alpha_table_path$Group <- factor(alpha_table_path$Group, levels = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH"))

alpha_table_bee <- alpha_table_path

#########################  chao1 
chao1_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(chao1 ~ Group, p.adjust.method = "fdr") 

chao1_pvalue_statistics <- chao1_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)

chao1_pvalue_statistics_sel <- chao1_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
chao1_pvalue_statistics_sel$p.adj <- format(chao1_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)


chao1_plot <- ggboxplot(alpha_table_bee, x = "Group", y = "chao1", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means() +
  geom_beeswarm(cex = 1.5) +
  #geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH"), 
                    values = color,
                    labels = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH")) +
  scale_x_discrete(labels=c("HE+SH" = "HE+SH\nFresh",
                            "HE2+SH" = "HE2+SH\nUp to 6days",  # (Up to 6days)
                            "HE3+SH" = "HE3+SH\n7-19d",  # (11-19days)
                            "HE4+SH" = "HE4+SH\n20-59d",  # (21-30days)
                            "HE5+SH" = "HE5+SH\n>59d") # # (>=61days)
  )+ 
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

chao1_plot




shannon_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(shannon ~ Group, p.adjust.method = "fdr") 

shannon_pvalue_statistics <- shannon_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)

shannon_pvalue_statistics_sel <- shannon_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
shannon_pvalue_statistics_sel$p.adj <- format(shannon_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)
shannon_pvalue_statistics_sel$y.position
#5.50664 5.52944

shannon_plot <- ggboxplot(alpha_table_bee, x = "Group", y = "shannon", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means() +
  geom_beeswarm(cex = 1.5) +
  #geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH"), 
                    values = color,
                    labels = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH")) +
  scale_x_discrete(labels=c("HE+SH" = "HE+SH\nFresh",
                            "HE2+SH" = "HE2+SH\nUp to 6days",  # (Up to 6days)
                            "HE3+SH" = "HE3+SH\n7-19d",  # (11-19days)
                            "HE4+SH" = "HE4+SH\n20-59d",  # (21-30days)
                            "HE5+SH" = "HE5+SH\n>59d") # # (>=61days)
  )+ 
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

shannon_plot

# Shannon–Wiener and Simpson index. However, while the Shannon–Wiener index is strongly influenced by species richness and by
# rare species, the Simpson index gives more weight to evenness and common species. 
# https://www.researchgate.net/post/What-is-the-deference-between-Shannon-Wiener-diversity-Index-and-Simpson-diversity-Index


simpson_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(simpson ~ Group, p.adjust.method = "fdr") 

simpson_pvalue_statistics <- simpson_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)

simpson_pvalue_statistics$y.position

simpson_pvalue_statistics_sel <- simpson_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
simpson_pvalue_statistics_sel$p.adj <- format(simpson_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)



simpson_plot <- ggboxplot(alpha_table_bee, x = "Group", y = "simpson", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means() +
  geom_beeswarm(cex = 1.5) +
  #geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH"), 
                    values = color,
                    labels = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH")) +
  scale_x_discrete(labels=c("HE+SH" = "HE+SH\nFresh",
                            "HE2+SH" = "HE2+SH\nUp to 6days",  # (Up to 6days)
                            "HE3+SH" = "HE3+SH\n7-19d",  # (11-19days)
                            "HE4+SH" = "HE4+SH\n20-59d",  # (21-30days)
                            "HE5+SH" = "HE5+SH\n>59d") # # (>=61days)
  )+ 
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

simpson_plot


pdf(paste0(output_dir, "/simpson_path.pdf"))
simpson_plot
dev.off()

pdf(paste0(output_dir, "/shannon_path.pdf"))
shannon_plot
dev.off()

pdf(paste0(output_dir, "/chao1_path.pdf"))
chao1_plot
dev.off()

pdf(paste0(output_dir, "/HE_series_alpha_diversity_combine_path.pdf"), width = 12)
ggpubr::ggarrange(shannon_plot, simpson_plot, chao1_plot, ncol = 3)
dev.off()


saveRDS(shannon_pvalue_statistics, paste0(output_dir, "/shannon_statistics_path.rds"))
saveRDS(simpson_pvalue_statistics, paste0(output_dir, "/simpson_statistics_path.rds"))
saveRDS(chao1_pvalue_statistics, paste0(output_dir, "/chao1_statistics_path.rds"))




save.image("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/functional_analysis/picrust2_out_pipeline/pathways_out/pathways.RData")



