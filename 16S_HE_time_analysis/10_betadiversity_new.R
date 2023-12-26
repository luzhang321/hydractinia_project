
#######################################################################################################


#                                   load library                                                      # 


#######################################################################################################

# script updated by remove agepoint 28, 40 point depending on the alpha diversity 
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(data.table)
library(vegan)

input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/9_diversity/"
load(paste0(input_dir, "9_alpha_diversity_analysis_HE_series_0909.RData"))

#######################################################################################################


#                prepare input                                                                        # 


#######################################################################################################

col - 7 
beta_diversity_input <- otu_table_average_norm[,1:66] %>%
  column_to_rownames(., var = "OTU") %>% 
  select(alpha_table_add_sel$sample_id_merged)

beta_diverstiy_metadata <- alpha_table_add_sel


########################################################################################################


#                    beta diversity statistics                                                         # 

########################################################################################################

metatable_sel <- beta_diverstiy_metadata 
otu_table <- beta_diversity_input
  
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
  p <- ggscatter(points_df, x = "x", y = "y",
                 color = "Group", palette = color,
                 shape = "Group",
                 ellipse = TRUE, 
                 ellipse.type = "convex"
                 ) + 
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    annotate(geom="text", x = 0.2, y = 0.3, label = bquote(P == ~ .(statistics[["adonisp"]]))) + 
    annotate(geom="text", x = 0.2, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_beta_diversity)))
  
  

  return(p)
}



beta_diversity_statistics <- sta_beta_result(otu_table = beta_diversity_input, metatable_sel = beta_diverstiy_metadata)
head(beta_diversity_statistics)


output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/10_beta_diversity/"
pdf(paste(output_dir, "HE_series_betadiversity.pdf", sep = "/"))
beta_diversity_statistics[["plot"]]
dev.off()



# pairwise adonis  

metatable_sel_order <- data.table(sample_id_merged = names(beta_diversity_input)) %>%
  left_join(., beta_diverstiy_metadata) 

#adonis(distance.bray ~ Group, data = metatable_sel_order, distance = "bray")

input <- beta_diversity_input
input_t <- input %>%
  as.matrix(.) %>%
  t(.) 
meta <- metatable_sel_order

pairwise_table <- pairwise_adonis(input = beta_diversity_input, meta = metatable_sel_order)
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

#remotes::install_github("Jtrachsel/funfuns")
pairwise_adonis_result <- funfuns::pairwise.adonis(input_t, meta$Group)
pairwise_adonis_result
pairs  F.Model         R2 p.value p.adjusted
1   HE+SH vs HE3+SH 4.277013 0.16920565   0.001      0.001
2   HE+SH vs HE2+SH 4.110002 0.15160463   0.001      0.001
3   HE+SH vs HE4+SH 2.546224 0.15388552   0.001      0.001
4   HE+SH vs HE5+SH 3.162991 0.09537712   0.001      0.001
5  HE3+SH vs HE2+SH 3.427799 0.12497534   0.001      0.001
6  HE3+SH vs HE4+SH 1.978778 0.11654421   0.031      0.031
7  HE3+SH vs HE5+SH 4.485817 0.12641154   0.001      0.001
8  HE2+SH vs HE4+SH 2.688626 0.13655732   0.001      0.001
9  HE2+SH vs HE5+SH 2.979453 0.08280985   0.001      0.001
10 HE4+SH vs HE5+SH 2.337696 0.08875857   0.001      0.001



write_tsv(pairwise_table, paste0(output_dir, "pairwise_table_myownfunction.tsv"))
write_tsv(pairwise_adonis_result, paste0(output_dir, "pairwise_table_package_funfuns.tsv"))
save.image(paste0(output_dir, "/10_betadiverstiy.RData"))



## figure modification - 2023.09.07

#############################################################
output_dir <- "/Volumes/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/10_beta_diversity/"
load(paste0(output_dir, "/10_betadiverstiy.RData"))
# add centroid to the beta diversity 
beta_diversity_statistics
set.seed(1234)
pcoa <- cmdscale(beta_diversity_statistics$dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) <- c("x", "y", "z")
eig <- pcoa$eig
points <- cbind(points, metatable_sel_order) # correct this one 

PcoA_plot_new <- function(points_df = NA, eig = NA, statistics = NA){
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

# remotes::install_github("corybrunson/ordr")


c <- PcoA_plot_new(points_df = points, eig = eig, statistics = beta_diversity_statistics) 

output_dir <- "/Volumes/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/10_beta_diversity/"
pdf(paste(output_dir, "HE_series_betadiversityv2.pdf", sep = "/"))
c
dev.off()

save.image(paste0(output_dir, "/10_betadiverstiy_figuremodification_0907.RData"))
