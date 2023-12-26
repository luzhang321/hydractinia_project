
#################################################################################################


#				loaded packages 						#


#################################################################################################

library(tidyverse)
library(dplyr)
library(tidyr)
library(metagenomeSeq)
library(stringr)
library(ggpubr)
library(ggbeeswarm)
library(data.table)

#################################################################################################


#				input file 							# 


#################################################################################################


#taxonomy_table_HE_series <- read_tsv("/sbidata/lzhang/201911_hydractinia/202207_16S_HE_time_series_analysis/8_taxonomy_export/exported/table-with-taxonomy-header-removed.tsv")
taxonomy_table_HE_series <- read_tsv("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/table-with-taxonomy.tsv", skip = 1)
#metadata_HE_series <- read_tsv("/sbidata/lzhang/201911_hydractinia/202207_16S_HE_time_series_analysis/7_metadata.tsv")
metadata_HE_series <- read_tsv("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/metadata_combine_merge_revised.tsv")
metadata_HE_series_simple <- metadata_HE_series %>%
  select(sample_id, subclass, sample_id_merged) %>%
  add_column(SampleID = str_replace_all(.$sample_id, "_", ".")) %>%
  select(-sample_id) %>% 
  dplyr::relocate(SampleID)

# the fresh samples metadata 
metadata_fresh <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/metatable_sel.csv") %>% 
  filter(Class %in% c(2,4))
metadata_fresh_with_biorepname <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/metatable_sel_biorep.csv") %>%
  filter(Class %in% c(2,4))
# otu_table_fresh <- read_csv("/sbidata/lzhang/201911_hydractinia/202207_16S_HE_time_series_analysis/10_average_count_table.csv")
metadata_fresh_with_biorepname_simple <- metadata_fresh_with_biorepname %>% 
  select(SampleID, uniquebiorep) %>%
  add_column(subclass = "Fresh") %>% 
  mutate(sample_id_merged = uniquebiorep) %>%
  select(-uniquebiorep)
  
metadata_combine_all <- rbind(metadata_HE_series_simple, metadata_fresh_with_biorepname_simple) %>% 
  filter(SampleID != "May2020.Class8.SClass2.B1.T2") #139 samples 
# i check the sample reads number, obviously the sample 
# May2020_Class8_SClass2_B1_T2 sample should be removed in the analysis 


#################################################################################################


#			separately the taxonomy 						#


#################################################################################################

otu_table <- taxonomy_table_HE_series
metatable <- metadata_combine_all
read_table_fun <- function(otu_table = NA, metatable = NA){
  # select the target samples and remove the asv with rowsum = 0 in this ASVs 
  colnames(otu_table)[1] <- "OTUID"
  otu_table_long <- otu_table %>%
    select("OTUID", metatable$SampleID, "taxonomy") %>%
    pivot_longer(cols = -c(OTUID, taxonomy))
  otu_table_long_keep <- otu_table_long %>%
    group_by(OTUID) %>%
    summarise(sum = sum(value)) %>% 
    filter(sum >0 )
  
  # separate taxonomy and only select the non-rowsum zero ASV
  otu_table_sel <- otu_table %>% 
    filter(OTUID %in% otu_table_long_keep$OTUID) %>%
    select("OTUID", metatable$SampleID, "taxonomy") 
  
  taxonomy <- str_split(otu_table_sel$taxonomy,"; ",simplify = T)
  colnames(taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  otu_table_sel <- otu_table_sel %>%
    select(-taxonomy) %>%
    cbind(., taxonomy, stringsAsFactors = F)
  
  
  return(otu_table_sel)
}

taxonomy_table_HE_series_clean <- read_table_fun(otu_table = taxonomy_table_HE_series, metatable = metadata_combine_all)


taxonomy_table_HE_series_tax <- select(taxonomy_table_HE_series_clean, OTUID, Kingdom, Phylum, Class,
                  Order, Family, Genus, Species)			        

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/9_diversity/"

write_tsv(taxonomy_table_HE_series_tax, paste0(output_dir, "HE_series_taxon.tsv"))#9702 asvs 

#################################################################################################


#			Check the total number of ASVs                                                						#


#################################################################################################

HE_sample_total_ASV <- taxonomy_table_HE_series_clean %>% 
  select(., -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
  filter(count > 0) %>%
  add_column(num = 1) %>%
  group_by(SampleID) %>%
  summarise(asv_number = sum(num), total_counts = sum(count))





#################################################################################################


#                       Average the counts                                                      #


#################################################################################################

# average OTU count table and do normalization again
HE_series_average <- taxonomy_table_HE_series_clean %>% 
  select(., -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
  left_join(., metadata_combine_all, by = "SampleID") %>%
  group_by(sample_id_merged, OTUID) %>%
  summarize(mean_count = round(mean(count, na.rm = TRUE))) %>%
  pivot_wider(names_from = "sample_id_merged", values_from = "mean_count") %>%
  left_join(., taxonomy_table_HE_series_tax, by = "OTUID")

write_tsv(HE_series_average, paste0(output_dir, "HE_series_average_count_table.tsv"))


#################################################################################################


#                       check the total ASV number                                              #


#################################################################################################

HE_series_average

HE_sample_average_total_ASV <- HE_series_average %>% 
  select(., -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
  filter(count > 0) %>%
  add_column(num = 1) %>%
  group_by(SampleID) %>%
  summarise(asv_number = sum(num), total_counts = sum(count))




#################################################################################################


#			Normalization 					                                                            			#


#################################################################################################

# normalization

Normalization_metagenome <- function(OTUfile = NA, col = NA, outdir = NA){
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
  write_csv(normfile, paste(outdir, "HE_series_average_normalized_otu_table.csv", sep = ""))
  return(normfile)
}


col <- ncol(HE_series_average)
otu_table_average_norm <- Normalization_metagenome(OTUfile= HE_series_average, col = col, outdir = output_dir) 
col - 7
otu_table_average_norm_notaxa <- otu_table_average_norm %>%
  .[,1:66] %>%
  column_to_rownames(., var = "OTU") %>%
  rownames_to_column(., var = "OTU")


HE_sample_normalized_total_ASV <- otu_table_average_norm_notaxa %>% 
  pivot_longer(., -OTU, names_to = "SampleID", values_to = "count") %>%
  filter(count > 0) %>%
  add_column(num = 1) %>%
  group_by(SampleID) %>%
  summarise(asv_number = sum(num), total_counts = sum(count))
View(HE_sample_normalized_total_ASV)



#################################################################################################


#			alpha diversity calculation 					                                                  	#


#################################################################################################

# alpha diversity 
alpha_diversity <- function(otu_table = NA, metatable_sel = NA, name = NA){
  
  otu <- otu_table %>%
    #column_to_rownames(., var = name) %>%
    .[,-1] %>%
    as.matrix(.) %>%
    t(.) %>%
    as.data.frame()
  
  otu_shannon <- vegan::diversity(otu,"shannon")
  
  otu_simpson <- vegan::diversity(otu, "simpson")
  
  otu_chao1 <- apply(otu,1,fossil::chao1)
  
  # make table 
  alpha_df <- data.table(shannon = otu_shannon, simpson = otu_simpson, chao1 = otu_chao1) %>%
    add_column(sample_id_merged = names(otu_shannon), name = name) %>% 
    left_join(., metatable_sel, by = "sample_id_merged", sort = FALSE)
  
  return(alpha_df)
  
}

metadata_HE_series_unique <- metadata_combine_all %>% select(subclass, sample_id_merged) %>% unique() #65samples 
alpha_table <- alpha_diversity(otu_table = otu_table_average_norm_notaxa, metatable_sel = metadata_HE_series_unique, name = "OTUtable")

write_csv(alpha_table, paste(output_dir, "HE_series_average_alpha_table.csv",sep = "/"))


#################################################################################################


#                      alpha diversity comparison                                               #


#################################################################################################
#load(paste0(output_dir, "/0_data.RData"))
#metadata_combine_merge# I checked all of them are without crab, and since they are all the same location aquarium. No need to adjust these 2 factors.

# I will use simple kruskal wallis test and posthoc dunn test for the comparison 



alpha_diversity_statistic_shannon <- rstatix::kruskal_test(alpha_table, shannon ~ subclass)
alpha_diversity_statistic_simpson <- rstatix::kruskal_test(alpha_table, simpson ~ subclass)
alpha_diversity_statistic_chao1 <- rstatix::kruskal_test(alpha_table, chao1 ~ subclass)

alpha_diversity_statistic_ks_all <- rbind(alpha_diversity_statistic_shannon,
                                          alpha_diversity_statistic_simpson,
                                          alpha_diversity_statistic_chao1)

alpha_diversity_statistic_ks_all

# A tibble: 3 x 6
#.y.         n statistic    df        p method        
#<chr>   <int>     <dbl> <int>    <dbl> <chr>         
#  1 shannon    54      22.5     4 0.000158 Kruskal-Wallis
#  2 simpson    54      23.1     4 0.000123 Kruskal-Wallis
#  3 chao1      54      22.1     4 0.000193 Kruskal-Wallis

alpha_diversity_dunn_shannon <- rstatix::dunn_test(alpha_table, shannon ~ subclass, p.adjust.method = "fdr")
alpha_diversity_dunn_simpson <- rstatix::dunn_test(alpha_table, simpson ~ subclass, p.adjust.method = "fdr")
alpha_diversity_dunn_chao1 <- rstatix::dunn_test(alpha_table, chao1 ~ subclass, p.adjust.method = "fdr")

alpha_diversity_dunn_all <- rbind(alpha_diversity_dunn_shannon,
                                  alpha_diversity_dunn_simpson,
                                  alpha_diversity_dunn_chao1)

alpha_diversity_dunn_all

write_tsv(alpha_diversity_dunn_all, paste0(output_dir, "/HE_series_dunn_test.tsv"))
write_tsv(alpha_diversity_statistic_ks_all, paste0(output_dir, "/HE_series_ks_test.tsv"))


#################################################################################################


#                      boxplot                                                                  #


#################################################################################################


# bee plot 

alpha_table_bee <- alpha_table
alpha_table_bee$Group <- as.character(alpha_table_bee$subclass)
alpha_table_bee$Group[alpha_table_bee$subclass == "Fresh"] <- "HE+SH"
alpha_table_bee$Group[alpha_table_bee$subclass %in% 0] <- "HE2+SH"
alpha_table_bee$Group[alpha_table_bee$subclass %in% 1] <- "HE3+SH"
alpha_table_bee$Group[alpha_table_bee$subclass %in% 2] <- "HE4+SH"
alpha_table_bee$Group[alpha_table_bee$subclass %in% 3] <- "HE5+SH"

# the group information 
# fresh : the class2+class4 control hydractinia with shell fresh samples as control
#8-0 Hydractinia- Shell HE NO HKIAquarium Up to 6 days - 
#8-1 Hydractinia- Shell HE NO HKIAquarium 7 to 19 days -
#8-2 Hydractinia- Shell HE NO HKIAquarium 20 to 59 days -
#8-3 Hydractinia- Shell HE NO HKIAquarium >= 60 days -
  
color <- c("#079F07","#068406","#056A05","#034F03", "#023502") #011A01
color <- viridis::viridis(10)[c(4,6,8,9,10)]

alpha_table_bee$Group <- factor(alpha_table_bee$Group, levels = c("HE+SH", "HE2+SH", "HE3+SH", "HE4+SH", "HE5+SH"))


#########################  chao1 
chao1_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(chao1 ~ Group, p.adjust.method = "fdr") 

chao1_pvalue_statistics <- chao1_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
#chao1_pvalue_statistics$y.position <- c(3100, 3300, 3500)
#chao1_pvalue_statistics$p.adj <- format(chao1_pvalue_statistics$p.adj, scientific = TRUE, digits = 3)
chao1_pvalue_statistics_sel <- chao1_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
chao1_pvalue_statistics_sel$p.adj <- format(chao1_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)




library(ggbeeswarm)   
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
#shannon_pvalue_statistics$p.adj <- format(shannon_pvalue_statistics$p.adj, scientific = TRUE, digits = 3)

shannon_pvalue_statistics_sel <- shannon_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
shannon_pvalue_statistics_sel$p.adj <- format(shannon_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)
shannon_pvalue_statistics_sel$y.position
#[1] 7.44248 7.88848 8.06688 8.15608
shannon_pvalue_statistics_sel$y.position <- c(7.44248, 7.88848, 8.06688, 8.21)

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



simpson_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(simpson ~ Group, p.adjust.method = "fdr") 

simpson_pvalue_statistics <- simpson_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)

simpson_pvalue_statistics$y.position
#[1] 1.00136 1.00176 1.00216 1.00256 1.00296 1.00336 1.00376 1.00416 1.00456 1.00496
simpson_pvalue_statistics$y.position <- c(1.001, 1.002, 1.003, 1.004, 1.005, 1.006, 1.007, 1.008, 1.009, 1.010)

#simpson_pvalue_statistics$p.adj <- format(simpson_pvalue_statistics$p.adj, scientific = TRUE, digits = 3)
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


pdf(paste0(output_dir, "/simpson.pdf"))
simpson_plot
dev.off()

pdf(paste0(output_dir, "/shannon.pdf"))
shannon_plot
dev.off()

pdf(paste0(output_dir, "/chao1.pdf"))
chao1_plot
dev.off()

pdf(paste0(output_dir, "/HE_series_alpha_diversity_combine.pdf"), width = 12)
ggpubr::ggarrange(shannon_plot, simpson_plot, chao1_plot, ncol = 3)
dev.off()


saveRDS(shannon_pvalue_statistics, paste0(output_dir, "/shannon_statistics.rds"))
saveRDS(simpson_pvalue_statistics, paste0(output_dir, "/simpson_statistics.rds"))
saveRDS(chao1_pvalue_statistics, paste0(output_dir, "/chao1_statistics.rds"))

############################################################################################################



#             go into details of age unit time and do the boxplot                                          # 

############################################################################################################
metadata_combine_merge_revised <- metadata_HE_series %>% select(sample_id_merged, age_unit) %>% unique()
alpha_table_add <- alpha_table_bee %>% 
  left_join(., metadata_combine_merge_revised, by = "sample_id_merged")

alpha_table_add$age_unit[is.na(alpha_table_add$age_unit)] <- "Fresh"
alpha_table_add$age_unit <- factor(alpha_table_add$age_unit, levels = c("Fresh", "6","10","19", "28", "30", "40", "60","70", "90"))

chao1_pvalue_statistics_age <- alpha_table_add %>%
  rstatix::dunn_test(chao1 ~ age_unit, p.adjust.method = "fdr") %>%
  filter(p.adj <= 0.05)

chao1_pvalue_statistics_age <- chao1_pvalue_statistics_age %>% 
  rstatix::add_xy_position(x = "age_unit", fun = "mean_sd", dodge = 0.8) 
chao1_pvalue_statistics_age$p.adj <- format(chao1_pvalue_statistics_age$p.adj, scientific = TRUE, digits = 3)

 
chao1_plot_age <- ggboxplot(alpha_table_add, x = "age_unit", y = "chao1", fill = "age_unit",outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics_age, tip.length = 0, label = "p.adj") +
  stat_compare_means() +
  geom_beeswarm(cex = 3) +
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

chao1_plot_age

simpson_pvalue_statistics_age <- alpha_table_add %>%
  rstatix::dunn_test(simpson ~ age_unit, p.adjust.method = "fdr") %>%
  filter(p.adj <= 0.05)

simpson_pvalue_statistics_age <- simpson_pvalue_statistics_age %>% 
  rstatix::add_xy_position(x = "age_unit", fun = "mean_sd", dodge = 0.8) 
simpson_pvalue_statistics_age$p.adj <- format(simpson_pvalue_statistics_age$p.adj, scientific = TRUE, digits = 3)


simpson_plot_age <- ggboxplot(alpha_table_add, x = "age_unit", y = "simpson", fill = "age_unit",outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics_age, tip.length = 0, label = "p.adj") +
  stat_compare_means() +
  geom_beeswarm(cex = 3) +
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

simpson_plot_age


shannon_pvalue_statistics_age <- alpha_table_add %>%
  rstatix::dunn_test(shannon ~ age_unit, p.adjust.method = "fdr") %>%
  filter(p.adj <= 0.05)

shannon_pvalue_statistics_age <- shannon_pvalue_statistics_age %>% 
  rstatix::add_xy_position(x = "age_unit", fun = "mean_sd", dodge = 0.8) 
shannon_pvalue_statistics_age$p.adj <- format(shannon_pvalue_statistics_age$p.adj, scientific = TRUE, digits = 3)

shannon_plot_age <- ggboxplot(alpha_table_add, x = "age_unit", y = "shannon", fill = "age_unit",outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics_age, tip.length = 0, label = "p.adj") +
  stat_compare_means() +
  geom_beeswarm(cex = 3) +
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

shannon_plot_age

pdf(paste0(output_dir, "/HE_series_alpha_diversity_combine_agetimeunit.pdf"), width = 30, height = 14)
ggarrange(shannon_plot_age, simpson_plot_age, chao1_plot_age, ncol = 3)
dev.off()

save.image(paste0(output_dir, "9_alpha_diversity_analysis_HE_series_new0808.RData"))


############################################################################################################



#             depending on the plot, we decide that we remove the 28 point sample, which are more likely to 
#             to be a outlier, also the outlier at 40d in  HE_series_alpha_diversity_combine_agetimeunit.pdf


############################################################################################################
# modification on september 9th 
output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/9_diversity/"


load(paste0(output_dir, "9_alpha_diversity_analysis_HE_series_new0808.RData"))

alpha_table_add_sel <- alpha_table_add %>%
  filter(age_unit != "28") %>%  #82 chao1 
  filter(sample_id_merged != "Aug2020_Class8_SClass3_B2") # only 126 chao1, it's also a outlier 

# 63 samples 


#########################  chao1 
chao1_pvalue_statistics <- alpha_table_add_sel %>%
  rstatix::dunn_test(chao1 ~ Group, p.adjust.method = "fdr") 

chao1_pvalue_statistics <- chao1_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
#chao1_pvalue_statistics$y.position <- c(3100, 3300, 3500)
#chao1_pvalue_statistics$p.adj <- format(chao1_pvalue_statistics$p.adj, scientific = TRUE, digits = 3)
chao1_pvalue_statistics_sel <- chao1_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
chao1_pvalue_statistics_sel$p.adj <- format(chao1_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)
chao1_pvalue_statistics_sel$y.position <- c(2500, 2600, 2700)



library(ggbeeswarm)   
chao1_plot <- ggboxplot(alpha_table_add_sel, x = "Group", y = "chao1", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means(label.y = 2800) +
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


shannon_pvalue_statistics <- alpha_table_add_sel %>%
  rstatix::dunn_test(shannon ~ Group, p.adjust.method = "fdr") 

shannon_pvalue_statistics <- shannon_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
#shannon_pvalue_statistics$p.adj <- format(shannon_pvalue_statistics$p.adj, scientific = TRUE, digits = 3)

shannon_pvalue_statistics_sel <- shannon_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
shannon_pvalue_statistics_sel$p.adj <- format(shannon_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)
shannon_pvalue_statistics_sel$y.position
#[1] 7.44248 7.88848 8.06688 8.15608
shannon_pvalue_statistics_sel$y.position <- c(7.44248, 7.58848, 7.66688)#, 8.21)

shannon_plot <- ggboxplot(alpha_table_add_sel, x = "Group", y = "shannon", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means(label.y = 7.76) +
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



simpson_pvalue_statistics <- alpha_table_add_sel %>%
  rstatix::dunn_test(simpson ~ Group, p.adjust.method = "fdr") 

simpson_pvalue_statistics <- simpson_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)

simpson_pvalue_statistics$y.position
#[1] 1.00136 1.00176 1.00216 1.00256 1.00296 1.00336 1.00376 1.00416 1.00456 1.00496
simpson_pvalue_statistics$y.position <- c(1.001, 1.002, 1.003, 1.004, 1.005, 1.006, 1.007, 1.008, 1.009, 1.010)

#simpson_pvalue_statistics$p.adj <- format(simpson_pvalue_statistics$p.adj, scientific = TRUE, digits = 3)
simpson_pvalue_statistics_sel <- simpson_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
simpson_pvalue_statistics_sel$y.position <- c(1.002,1.0025,1.003)
simpson_pvalue_statistics_sel$p.adj <- format(simpson_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)



simpson_plot <- ggboxplot(alpha_table_add_sel, x = "Group", y = "simpson", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means(label.y = 1.0035) +
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


pdf(paste0(output_dir, "/simpson.pdf"))
simpson_plot
dev.off()

pdf(paste0(output_dir, "/shannon.pdf"))
shannon_plot
dev.off()

pdf(paste0(output_dir, "/chao1.pdf"))
chao1_plot
dev.off()

pdf(paste0(output_dir, "/HE_series_alpha_diversity_combine.pdf"), width = 12)
ggpubr::ggarrange(shannon_plot, simpson_plot, chao1_plot, ncol = 3)
dev.off()


saveRDS(shannon_pvalue_statistics, paste0(output_dir, "/shannon_statistics.rds"))
saveRDS(simpson_pvalue_statistics, paste0(output_dir, "/simpson_statistics.rds"))
saveRDS(chao1_pvalue_statistics, paste0(output_dir, "/chao1_statistics.rds"))
save.image(paste0(output_dir, "9_alpha_diversity_analysis_HE_series_0909.RData"))
