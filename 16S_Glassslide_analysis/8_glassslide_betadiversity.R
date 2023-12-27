# beta diversity of glassslide samples 

#######################################################################################################


#                                   load library                                                      # 


#######################################################################################################


library(tidyverse)
library(ggpubr)
library(ggplot2)
library(data.table)
library(vegan)

input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/7_alphadiversity/"
load(paste0(input_dir, "glassslide_alpha_diversity.RData"))


#######################################################################################################


#                prepare input                                                                        # 


#######################################################################################################

39 - 7 
beta_diversity_input <- otu_table_average_norm[,1:32] %>%
  column_to_rownames(., var = "OTU") %>% 
  select(alpha_table_bee$sample_id_merged)

beta_diverstiy_metadata <- alpha_table_bee


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
color <- viridis::viridis(10)[c(4,6,10)]

PcoA_plot <- function(points_df = NA, eig = NA, statistics = NA){
  rsquare_beta_diversity <- round(statistics[["adonisr"]]$aov.tab$R2[1], 3)
  p <- ggscatter(points_df, x = "x", y = "y",
                 color = "Group", palette = color,
                 shape = "Group",
                 ellipse = TRUE, 
                 ellipse.type = "norm"
  ) + 
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    annotate(geom="text", x = 0.2, y = 0.3, label = bquote(P == ~ .(statistics[["adonisp"]]))) + 
    annotate(geom="text", x = 0.2, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_beta_diversity)))
  
  
  
  return(p)
}



beta_diversity_statistics <- sta_beta_result(otu_table = beta_diversity_input, metatable_sel = beta_diverstiy_metadata)
head(beta_diversity_statistics)

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/8_betadiversity/"
pdf(paste(output_dir, "HE_series_betadiversity.pdf", sep = "/"))
beta_diversity_statistics [["plot"]]
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
#> pairwise_table
# A tibble: 3 x 4
#V1    V2        R2 p_value
#<fct> <fct>  <dbl>   <dbl>
#  1 GS+HS GS     0.222   0.001
#2 GS+HS GS+HS* 0.219   0.001
#3 GS    GS+HS* 0.211   0.001


pairwise_adonis <- function(input = NA, meta = NA){
  meta$Group <- as.character(meta$Group)
  comparison <- combn(unique(meta$Group), 2)
  pairwise_adonis_tibble <- comparison %>% t() %>% as.data.frame() %>% .[1:3,] %>% as.tibble() %>%
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

pairwise_adonis_result <- funfuns::pairwise.adonis(input_t, meta$Group)

pairwise_adonis_result
pairs  F.Model        R2 p.value p.adjusted
1     GS+HS vs GS 5.149234 0.2224365   0.001      0.001
2 GS+HS vs GS+HS* 5.606565 0.2189503   0.001      0.001
3    GS vs GS+HS* 4.826990 0.2114598   0.001      0.001


write_tsv(pairwise_table, paste0(output_dir, "glassslide_pairwise_table_myownfunction.tsv"))
write_tsv(pairwise_adonis_result, paste0(output_dir, "glassslide_pairwise_table_package_funfuns.tsv"))
save.image(paste0(output_dir, "/10_glassslide_betadiverstiy.RData"))
