# the is to produce a combined boxplot for the significantly differntially abundant phylum/genus/family 


library(tidyverse)
library(ggpubr)

#################################################################################################


#				input file 							# 


#################################################################################################


# phylum level 
input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/9_sig_taxon/"
GS_GS_HS_sig_phylum <- read_csv(paste0(input_dir, "Taxonomy/sig_Phylum_GS_GS_HS.csv"))
GS_GS_HS_dead_sig_phylum <- read_csv(paste0(input_dir, "Taxonomy/sig_Phylum_GS_GS_HS_dead.csv"))
GS_HS_GS_HS_dead_sig_phylum <- read_csv(paste0(input_dir, "Taxonomy/sig_Phylum_GS_HS_GS_HS_dead.csv"))

Phylum_normalized_table <- read_tsv(paste0(input_dir, "Taxonomy_distribution/Phylum_combined_normalized.tsv"))

# family level 

input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/9_sig_taxon/"
GS_GS_HS_sig_family <- read_csv(paste0(input_dir, "Taxonomy/sig_Family_GS_GS_HS.csv"))
GS_GS_HS_dead_sig_family <- read_csv(paste0(input_dir, "Taxonomy/sig_Family_GS_GS_HS_dead.csv"))
GS_HS_GS_HS_dead_sig_family <- read_csv(paste0(input_dir, "Taxonomy/sig_Family_GS_HS_GS_HS_dead.csv"))

Family_normalized_table <- read_tsv(paste0(input_dir, "Taxonomy_distribution/Family_combined_normalized.tsv"))


# genera level 

input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/9_sig_taxon/"
GS_GS_HS_sig_genera <- read_csv(paste0(input_dir, "Taxonomy/sig_Genus_GS_GS_HS.csv"))
GS_GS_HS_dead_sig_genera <- read_csv(paste0(input_dir, "Taxonomy/sig_Genus_GS_GS_HS_dead.csv"))
GS_HS_GS_HS_dead_sig_genera <- read_csv(paste0(input_dir, "Taxonomy/sig_Genus_GS_HS_GS_HS_dead.csv"))

Genera_normalized_table <- read_tsv(paste0(input_dir, "Taxonomy_distribution/Genus_combined_normalized.tsv"))


# metadata 

input_dir2 <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/7_alphadiversity/"
load(paste0(input_dir2, "glassslide_alpha_diversity.RData"))

metadata <- alpha_table_bee
output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/9_sig_taxon/Taxonomy_boxplot/"
write_tsv(metadata, paste0(output_dir, "metadata.tsv"))




#################################################################################################


#				select the overlap of the 3 significant list and pick the overlap           						# 


#################################################################################################

# phylum 
union_phylum <- tibble(phylum = c(GS_GS_HS_sig_phylum$Taxa, GS_GS_HS_dead_sig_phylum$Taxa, GS_HS_GS_HS_dead_sig_phylum$Taxa)) %>%
  add_column(Group = c(rep("GS_GS_HS", length(GS_GS_HS_sig_phylum$Taxa)), rep("GS_GS_HS_dead", length(GS_GS_HS_dead_sig_phylum$Taxa)), rep("GS_HS_GS_HS_dead", length(GS_HS_GS_HS_dead_sig_phylum$Taxa)))) %>%
  add_column(num = 1) %>% 
  group_by(phylum) %>% 
  summarise(Group = paste(Group, collapse = ";"), num = sum(num)) %>%
  ungroup() %>% 
  arrange(desc(num))

nrow(union_phylum %>% filter(num >= 3)) #2
nrow(union_phylum %>% filter(num >= 2)) #15
nrow(union_phylum %>% filter(num >= 1)) #22
             
overlap_3_phylum <- union_phylum %>% filter(num >= 3)
overlap_2_phylum <- union_phylum %>% filter(num >= 2)
overlap_2_phylum_no_NA <- overlap_2_phylum %>% filter(!is.na(phylum)) %>% filter(phylum != "p__")
overlap_1_phylum <- union_phylum %>% filter(num >= 1)

# family 
union_family <- tibble(family = c(GS_GS_HS_sig_family$Taxa, GS_GS_HS_dead_sig_family$Taxa, GS_HS_GS_HS_dead_sig_family$Taxa)) %>%
  add_column(Group = c(rep("GS_GS_HS", length(GS_GS_HS_sig_family$Taxa)), rep("GS_GS_HS_dead", length(GS_GS_HS_dead_sig_family$Taxa)), rep("GS_HS_GS_HS_dead", length(GS_HS_GS_HS_dead_sig_family$Taxa)))) %>%
  add_column(num = 1) %>% 
  group_by(family) %>% 
  summarise(Group = paste(Group, collapse = ";"), num = sum(num)) %>%
  ungroup() %>% 
  arrange(desc(num))

nrow(union_family %>% filter(num >= 3)) #3
nrow(union_family %>% filter(num >= 2)) #50
nrow(union_family %>% filter(num >= 1)) #85

overlap_3_family <- union_family %>% filter(num >= 3)
overlap_2_family <- union_family %>% filter(num >= 2)
overlap_2_family_no_NA <- overlap_2_family %>% filter(!is.na(family)) %>% filter(family != "f__")
overlap_1_family <- union_family %>% filter(num >= 1)


# genera 
union_genera <- tibble(genera = c(GS_GS_HS_sig_genera$Taxa, GS_GS_HS_dead_sig_genera$Taxa, GS_HS_GS_HS_dead_sig_genera$Taxa)) %>%
  add_column(Group = c(rep("GS_GS_HS", length(GS_GS_HS_sig_genera$Taxa)), rep("GS_GS_HS_dead", length(GS_GS_HS_dead_sig_genera$Taxa)), rep("GS_HS_GS_HS_dead", length(GS_HS_GS_HS_dead_sig_genera$Taxa)))) %>%
  add_column(num = 1) %>% 
  group_by(genera) %>% 
  summarise(Group = paste(Group, collapse = ";"), num = sum(num)) %>%
  ungroup() %>% 
  arrange(desc(num))

nrow(union_genera %>% filter(num >= 3)) #11
nrow(union_genera %>% filter(num >= 2)) #85
nrow(union_genera %>% filter(num >= 1)) #133

overlap_3_genera <- union_genera %>% filter(num >= 3)
overlap_3_genera_no_NA <- overlap_3_genera %>% filter(!is.na(genera)) %>% filter(genera != "g__")

overlap_2_genera <- union_genera %>% filter(num >= 2)
overlap_1_genera <- union_genera %>% filter(num >= 1)



#################################################################################################


#			prepare	draw boxplot shown taxa : at 2 overlap  -phylum                        						# 


#################################################################################################


GS_GS_HS_sig_phylum_simple <- GS_GS_HS_sig_phylum %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_GS_HS")
GS_GS_HS_dead_sig_phylum_simple <- GS_GS_HS_dead_sig_phylum %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_GS_HS_dead")
GS_HS_GS_HS_dead_sig_phylum_simple <- GS_HS_GS_HS_dead_sig_phylum %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_HS_GS_HS_dead")

Phylum_all_simple_table <- rbind(GS_GS_HS_sig_phylum_simple,
                                 GS_GS_HS_dead_sig_phylum_simple,
                                 GS_HS_GS_HS_dead_sig_phylum_simple)


Phylum_all_simple_table_sel <- Phylum_all_simple_table %>% filter(Taxa %in% overlap_2_phylum_no_NA$phylum)

nrow(Phylum_all_simple_table_sel) #28
colnames(Phylum_all_simple_table_sel)[3] <- "taxon"
Phylum_all_simple_table_sel$group1 <- str_split(Phylum_all_simple_table_sel$Group, "_GS", simplify = T) %>% .[,1] %>% str_replace_all(., "_", "+")
Phylum_all_simple_table_sel$group2 <- str_remove(Phylum_all_simple_table_sel$Group, "GS_") %>% 
  str_remove(., "^HS_") %>% 
  str_replace_all(., "_", "+") %>%
  str_replace(., "\\+dead", "\\*")


#################################################################################################


#			draw boxplot         - phylum level                                                           						# 


#################################################################################################


Phylum_normalized_table_metadata <- Phylum_normalized_table %>%
  left_join(., metadata %>% select(-name), by = c("name" = "sample_id_merged")) 

Phylum_normalized_table_metadata_sta <- Phylum_normalized_table_metadata %>%
  filter(taxon %in% unique(overlap_2_phylum_no_NA$phylum)) %>% 
  group_by(taxon) %>%
  rstatix::t_test(abundance ~ Group) 

Phylum_normalized_table_metadata_sta <- Phylum_normalized_table_metadata_sta %>% rstatix::add_xy_position(x = "taxon", fun = "mean_sd", dodge = 0.8)


Phylum_normalized_table_metadata_sta_combine <- Phylum_normalized_table_metadata_sta %>%
  left_join(., Phylum_all_simple_table_sel, by = c("taxon", "group1", "group2")) %>% 
  filter(!is.na(pvalues))

Phylum_normalized_table_metadata_sta_combine$p <- Phylum_normalized_table_metadata_sta_combine$pvalues
Phylum_normalized_table_metadata_sta_combine$p.adj <- format(Phylum_normalized_table_metadata_sta_combine$adjPvalues, scientific = TRUE, digits = 3)
Phylum_normalized_table_metadata_sta_combine$y.position <- 3.5
Phylum_normalized_table_metadata$log_abundance <- log10(Phylum_normalized_table_metadata$abundance + 0.01)
Phylum_normalized_table_metadata_sta_combine$y.position[Phylum_normalized_table_metadata_sta_combine$group1 == "GS"] <- 2.3
Phylum_normalized_table_metadata_sta_combine$y.position[Phylum_normalized_table_metadata_sta_combine$group1 == "GS" & Phylum_normalized_table_metadata_sta_combine$group2 == "GS+HS*"] <- 2.5
Phylum_normalized_table_metadata_sta_combine$y.position[Phylum_normalized_table_metadata_sta_combine$group1 == "GS+HS"] <- 2.7

color <- viridis::viridis(10)[c(4,6,10)]

library(ggprism)
Phylum_normalized_table_metadata_sel <- Phylum_normalized_table_metadata %>%
  filter(taxon %in% overlap_2_phylum_no_NA$phylum)

phylum_boxplot <- ggboxplot(Phylum_normalized_table_metadata_sel, x = "taxon", y = "log_abundance", fill = "Group") +
  stat_pvalue_manual(Phylum_normalized_table_metadata_sta_combine, tip.length = 0.0, label = "p.adj") +
  scale_fill_manual(breaks = c("GS", "GS+HS", "GS+HS*"), 
                    values = color,
                    labels = c("GS", "GS+HS", "GS+HS*")) + 
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  ) +
  rotate_x_text(45) +
  ylab("Log10 Abundance") +
  xlab("Phylum")
  

ggsave(paste0(output_dir, "phylum_boxplot_glassslide.pdf"), phylum_boxplot, width = 15)




#################################################################################################


#			prepare	draw boxplot shown taxa : at 2 overlap  - family                      						# 


#################################################################################################




GS_GS_HS_sig_family_simple <- GS_GS_HS_sig_family %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_GS_HS")
GS_GS_HS_dead_sig_family_simple <- GS_GS_HS_dead_sig_family %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_GS_HS_dead")
GS_HS_GS_HS_dead_sig_family_simple <- GS_HS_GS_HS_dead_sig_family %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_HS_GS_HS_dead")

Family_all_simple_table <- rbind(GS_GS_HS_sig_family_simple,
                                 GS_GS_HS_dead_sig_family_simple,
                                 GS_HS_GS_HS_dead_sig_family_simple)


Family_all_simple_table_sel <- Family_all_simple_table %>% filter(Taxa %in% overlap_3_family$family)

nrow(Family_all_simple_table_sel) #103
colnames(Family_all_simple_table_sel)[3] <- "taxon"
Family_all_simple_table_sel$group1 <- str_split(Family_all_simple_table_sel$Group, "_GS", simplify = T) %>% .[,1] %>% str_replace_all(., "_", "+")
Family_all_simple_table_sel$group2 <- str_remove(Family_all_simple_table_sel$Group, "GS_") %>% 
  str_remove(., "^HS_") %>% 
  str_replace_all(., "_", "+") %>%
  str_replace(., "\\+dead", "\\*")


#################################################################################################


#			draw boxplot         - family level                                                           						# 


#################################################################################################


Family_normalized_table_metadata <- Family_normalized_table %>%
  left_join(., metadata %>% select(-name), by = c("name" = "sample_id_merged")) 

Family_normalized_table_metadata_sta <- Family_normalized_table_metadata %>%
  filter(taxon %in% unique(overlap_3_family$family)) %>% 
  group_by(taxon) %>%
  rstatix::t_test(abundance ~ Group) 

Family_normalized_table_metadata_sta <- Family_normalized_table_metadata_sta %>% rstatix::add_xy_position(x = "taxon", fun = "mean_sd", dodge = 0.8)


Family_normalized_table_metadata_sta_combine <- Family_normalized_table_metadata_sta %>%
  left_join(., Family_all_simple_table_sel, by = c("taxon", "group1", "group2")) %>% 
  filter(!is.na(pvalues))

Family_normalized_table_metadata_sta_combine$p <- Family_normalized_table_metadata_sta_combine$pvalues
Family_normalized_table_metadata_sta_combine$p.adj <- format(Family_normalized_table_metadata_sta_combine$adjPvalues, scientific = TRUE, digits = 3)
Family_normalized_table_metadata_sta_combine$y.position <- 3.5
Family_normalized_table_metadata$log_abundance <- log10(Family_normalized_table_metadata$abundance + 0.01)
Family_normalized_table_metadata_sta_combine$y.position[Family_normalized_table_metadata_sta_combine$group1 == "GS"] <- 2.3
Family_normalized_table_metadata_sta_combine$y.position[Family_normalized_table_metadata_sta_combine$group1 == "GS" & Family_normalized_table_metadata_sta_combine$group2 == "GS+HS*"] <- 2.5
Family_normalized_table_metadata_sta_combine$y.position[Family_normalized_table_metadata_sta_combine$group1 == "GS+HS"] <- 2.7


library(ggprism)
Family_normalized_table_metadata_sel <- Family_normalized_table_metadata %>%
  filter(taxon %in% overlap_3_family$family)

family_boxplot <- ggboxplot(Family_normalized_table_metadata_sel, x = "taxon", y = "log_abundance", fill = "Group") +
  stat_pvalue_manual(Family_normalized_table_metadata_sta_combine, tip.length = 0.0, label = "p.adj") +
  scale_fill_manual(breaks = c("GS", "GS+HS", "GS+HS*"), 
                    values = color,
                    labels = c("GS", "GS+HS", "GS+HS*")) + 
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  ) +
  rotate_x_text(45) +
  ylab("Log10 Abundance") +
  xlab("Family")


ggsave(paste0(output_dir, "family_boxplot_glassslide.pdf"), family_boxplot, width = 7)


#################################################################################################


#			prepare	draw boxplot shown taxa : at 2 overlap  - genera                     						# 


#################################################################################################




GS_GS_HS_sig_genera_simple <- GS_GS_HS_sig_genera %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_GS_HS")
GS_GS_HS_dead_sig_genera_simple <- GS_GS_HS_dead_sig_genera %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_GS_HS_dead")
GS_HS_GS_HS_dead_sig_genera_simple <- GS_HS_GS_HS_dead_sig_genera %>% select(pvalues, adjPvalues, Taxa, Group_Impact) %>% add_column(Group = "GS_HS_GS_HS_dead")

Genera_all_simple_table <- rbind(GS_GS_HS_sig_genera_simple,
                                 GS_GS_HS_dead_sig_genera_simple,
                                 GS_HS_GS_HS_dead_sig_genera_simple)


Genera_all_simple_table_sel <- Genera_all_simple_table %>% filter(Taxa %in% overlap_3_genera$genera)

nrow(Genera_all_simple_table_sel) #33
colnames(Genera_all_simple_table_sel)[3] <- "taxon"
Genera_all_simple_table_sel$group1 <- str_split(Genera_all_simple_table_sel$Group, "_GS", simplify = T) %>% .[,1] %>% str_replace_all(., "_", "+")
Genera_all_simple_table_sel$group2 <- str_remove(Genera_all_simple_table_sel$Group, "GS_") %>% 
  str_remove(., "^HS_") %>% 
  str_replace_all(., "_", "+") %>%
  str_replace(., "\\+dead", "\\*")



#################################################################################################


#			draw boxplot         - genera level                                                           						# 


#################################################################################################


Genera_normalized_table_metadata <- Genera_normalized_table %>%
  left_join(., metadata %>% select(-name), by = c("name" = "sample_id_merged")) 

Genera_normalized_table_metadata_sta <- Genera_normalized_table_metadata %>%
  filter(taxon %in% unique(overlap_3_genera$genera)) %>% 
  group_by(taxon) %>%
  rstatix::t_test(abundance ~ Group) 

Genera_normalized_table_metadata_sta <- Genera_normalized_table_metadata_sta %>% rstatix::add_xy_position(x = "taxon", fun = "mean_sd", dodge = 0.8)


Genera_normalized_table_metadata_sta_combine <- Genera_normalized_table_metadata_sta %>%
  left_join(., Genera_all_simple_table_sel, by = c("taxon", "group1", "group2")) %>% 
  filter(!is.na(pvalues))

Genera_normalized_table_metadata_sta_combine$p <- Genera_normalized_table_metadata_sta_combine$pvalues
Genera_normalized_table_metadata_sta_combine$p.adj <- format(Genera_normalized_table_metadata_sta_combine$adjPvalues, scientific = TRUE, digits = 3)
Genera_normalized_table_metadata_sta_combine$y.position <- 3.5
Genera_normalized_table_metadata$log_abundance <- log10(Genera_normalized_table_metadata$abundance + 0.01)
Genera_normalized_table_metadata_sta_combine$y.position[Genera_normalized_table_metadata_sta_combine$group1 == "GS"] <- 2.3
Genera_normalized_table_metadata_sta_combine$y.position[Genera_normalized_table_metadata_sta_combine$group1 == "GS" & Genera_normalized_table_metadata_sta_combine$group2 == "GS+HS*"] <- 2.5
Genera_normalized_table_metadata_sta_combine$y.position[Genera_normalized_table_metadata_sta_combine$group1 == "GS+HS"] <- 2.7


library(ggprism)
Genera_normalized_table_metadata_sel <- Genera_normalized_table_metadata %>%
  filter(taxon %in% overlap_3_genera$genera)

genera_boxplot <- ggboxplot(Genera_normalized_table_metadata_sel, x = "taxon", y = "log_abundance", fill = "Group") +
  stat_pvalue_manual(Genera_normalized_table_metadata_sta_combine, tip.length = 0.0, label = "p.adj") +
  scale_fill_manual(breaks = c("GS", "GS+HS", "GS+HS*"), 
                    values = color,
                    labels = c("GS", "GS+HS", "GS+HS*")) + 
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  ) +
  rotate_x_text(45) +
  ylab("Log10 Abundance") +
  xlab("Genus")


ggsave(paste0(output_dir, "genera_boxplot_glassslide.pdf"), genera_boxplot, width = 15)

save.image(paste0(output_dir, "11_glassslide_sig_taxon_boxplot.RData"))
