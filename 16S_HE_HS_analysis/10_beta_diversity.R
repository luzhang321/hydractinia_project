# aim of this script is to do beta diversity [namely bray-curtis pcoa]
library(tidyverse)
library(vegan)

# steps :
#1. use the tidied table from alpha-diversity comparisons;
#2. do beta diversity comparisons 

otu_table_average_norm <- readRDS("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/normalized_otu_table.rds")
  
metadata_sel <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/metatable_biorep_unique.csv")
metatable_sel_biorep_unique <- metadata_sel

otu_table_HS_control <- otu_table_average_norm %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metadata_sel) %>%
  filter(., Class %in%  c(23,33,32)) %>%
  .[,1:7274] %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")
metadata_HS_control <- filter(metatable_sel_biorep_unique, Class %in% c(23,33,32))


otu_table_HS_HE <- otu_table_average_norm %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metadata_sel) %>%
  filter(., Class %in%  c(1,2,3,4,23,33)) %>%
  .[,1:7274] %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")
  
metadata_HS_HE <- filter(metatable_sel_biorep_unique, Class %in% c(1,2,3,4,23,33))

# pcoa 

# PCoA plot 

PcoA_plot <- function(points_df, eig){
  # p <- ggplot(points_df, aes(x=x, y=y, color=Group))
  # p <- p  +
  #   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
  #        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + #title = title_name
  #   #main_theme +
  #   #labs(shape = "Class") +
  #   scale_color_manual(values=c("#4B0082", "golden rod")) +
  #   geom_point() +
  #   scale_shape_manual(values=seq(0,15)) +
  #   stat_ellipse()
  # 
  p <- ggscatter(points_df, x = "x", y = "y",
                 color = "Group", palette = "jco",
                 shape = "Group",
                 ellipse = TRUE, ellipse.type = "convex") + 
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) 
  
  
  
  #annotate(geom = "text", label = paste("p-value",pvalue, sep = ":"),x=(max(ICU_pheno_sel$x)), y = min(ICU_pheno_sel$y),hjust=1,color = "black", size = 3) 
  #scale_color_brewer(palette = "Set1")
  return(p)
}




sta_beta_result <- function(otu_table, metatable_sel, name){
  set.seed(5)
  distance.bray <- otu_table %>%
    #column_to_rownames(., var = name) %>%
    .[,-1] %>%
    as.matrix(.) %>%
    t(.) %>%
    vegdist(.,method = 'bray')
  
  
  metatable_sel_order <- tibble(uniquebiorep = names(otu_table)[2:length(names(otu_table))]) %>%
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
  
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[1]#first postion is group
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[1]#first postion is group
  
  
  
  stat_beta <- list()
  stat_beta[["disp"]] <- beta_dispersion_pvalue
  stat_beta[["anosim"]] <- anosim.result
  stat_beta[["adonisr"]] <- adonis_group
  stat_beta[["adonisp"]] <- adonis_group_p_value
  
  
  pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  colnames(points) <- c("x", "y", "z")
  eig <- pcoa$eig
  points <- cbind(points, metatable_sel_order) # correct this one 
  c <- PcoA_plot(points, eig) + ggtitle(paste(name,"P:",adonis_group_p_value,"\n","R square:", 
                                              round(adonis_group_r_square,3), sep = " "))
  stat_beta[["plot"]] <- c
  
  stat_beta[["dis"]] <- distance.bray
  
  return(stat_beta)
}


beta_HS_control <- sta_beta_result(otu_table_HS_control, metadata_HS_control, "Hydractinia symbiolongicarpus vs Control")


beta_HS_HE <- sta_beta_result(otu_table_HS_HE, metadata_HS_HE, "hydractinia echinata vs Hydractinia symbiolongicarpus")


# write it into pdf 

pdf("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/HS_control_beta.pdf", width = 12)
beta_HS_control[[5]]
dev.off()

pdf("//media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/HE_HS_beta.pdf", width = 12)
beta_HS_HE[[5]]
dev.off()


saveRDS(beta_HS_control, "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/HS_control_beta.rds")
saveRDS(beta_HS_HE, "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/HE_HS_beta.rds")
save.image("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/Beta.RData")

### 2022.07 plot modification ####################################################################################################################
## 2022.05 modification of the plot + 2022.06 i change the color 

# shellwithcrab as control
rsquare_beta_HS_control <- round(beta_HS_control[["adonisr"]]$aov.tab$R2[1], 3)
library(ggpubr)
library(ggplot2)
beta_diversity_beta_HS_control <- beta_HS_control[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 4) +
  scale_color_manual(breaks = c("HS", "HS_Shell"), 
                     values = c( "#f0027f", "#386cb0"),
                     labels = c("HS+SH", "SH")) + 
  scale_fill_manual(breaks = c("HS", "HS_Shell"), 
                    values = c( "#f0027f", "#386cb0"),
                    labels = c("HS+SH", "SH")) + 
  scale_shape_manual(breaks = c("HS", "HS_Shell"), values = c(16,16), labels = c("HS+SH", "SH")) +
  annotate(geom="text", x = 0.4, y = 0.3, label = bquote(P == ~ .(beta_HS_control[["adonisp"]]))) + 
  annotate(geom="text", x = 0.4, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_beta_HS_control)))

pdf("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/HE_HS_beta_update.pdf")
beta_diversity_beta_HS_control
dev.off()
save.image("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/Beta_202207_figuremodification.RData")

#load("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/Beta.RData")


# figure modification 

load("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/Beta_202207_figuremodification.RData")

beta_diversity_beta_HS_control_big <- beta_diversity_beta_HS_control + theme(
  legend.position="top",
  plot.title = element_text(size=11),
  axis.text.x = element_text(face="bold"),
  axis.text.y = element_text(face="bold"),
  axis.title.y = element_text(face="bold"),
  axis.title.x = element_text(face="bold"),
  text = element_text(size=15)
) 

pdf("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/HE_HS_beta_updatev2.pdf")
beta_diversity_beta_HS_control_big
dev.off()
save.image("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/10_betadiversity/Beta_202211_figuremodification.RData")



