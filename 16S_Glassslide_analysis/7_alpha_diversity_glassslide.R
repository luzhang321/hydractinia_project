
# for glassslide samples alpha diversity comparison 
library(dplyr)
library(tidyr)
library(metagenomeSeq)
library(stringr)
library(ggpubr)
library(ggbeeswarm)
library(data.table)
library(tidyverse)

#################################################################################################


#				input file 							# 


#################################################################################################


taxonomy_table_glasslide <- read_tsv("/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/7_alphadiversity/table-with-taxonomy.tsv", skip = 1)
metadata_glassslide <- read_tsv("/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/7_alphadiversity/0_metada_clean_classslide.tsv") # all no crab + location Aquarium
metadata_glassslide_simple <- metadata_glassslide %>%
  select(sample_id, Class, sample_id_merged) %>%
  add_column(SampleID = str_replace_all(.$sample_id, "_", ".")) %>%
  select(-sample_id) %>% 
  dplyr::relocate(SampleID)


metadata_combine_all <- metadata_glassslide_simple


#################################################################################################


#			separately the taxonomy 						#


#################################################################################################


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
  
  taxonomy_tibble <- str_split(otu_table_sel$taxonomy,"; ",simplify = T)
  colnames(taxonomy_tibble) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  otu_table_sel <- otu_table_sel %>%
    select(-taxonomy) %>%
    cbind(., taxonomy_tibble, stringsAsFactors = F)
  
  
  return(otu_table_sel)
}

otu_table <- taxonomy_table_glasslide
metatable <- metadata_glassslide_simple
taxonomy_table_glassslide_clean <- read_table_fun(otu_table = taxonomy_table_glasslide, metatable = metadata_glassslide_simple)


taxonomy_table_glassslide_tax <- select(taxonomy_table_glassslide_clean, OTUID, Kingdom, Phylum, Class,
                                       Order, Family, Genus, Species)			        

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022_glassslide_samples/7_alphadiversity/"

write_tsv(taxonomy_table_glassslide_tax, paste0(output_dir, "HS_glassslide_taxon.tsv"))#2186 asvs 

#################################################################################################


#			Check the total number of ASVs  in each sample                                 						#


#################################################################################################

glassslide_sample_total_ASV <- taxonomy_table_glassslide_clean %>% 
  select(., -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
  filter(count > 0) %>%
  add_column(num = 1) %>%
  group_by(SampleID) %>%
  summarise(asv_number = sum(num), total_counts = sum(count))

glassslide_sample_total_ASV %>% arrange(asv_number)
# some of the asv are very small : SampleID               asv_number total_counts
#<chr>                       <dbl>        <dbl>
#  1 Oct2020.Class20.B11.T1         21        28744
#2 May2020.Class15.BF.T1          30        28882
#3 May2020.Class15.BM1.T1         34        25495
#4 May2020.Class15.B2.T2          53        20792
#5 May2020.Class15.BM1.T2         65        23327 


#################################################################################################


#                       Average the counts                                                      #


#################################################################################################

# average OTU count table and do normalization again
glassslide_average <- taxonomy_table_glassslide_clean %>% 
  select(., -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
  left_join(., metadata_glassslide_simple, by = "SampleID") %>%
  group_by(sample_id_merged, OTUID) %>%
  summarize(mean_count = round(mean(count, na.rm = TRUE))) %>%
  pivot_wider(names_from = "sample_id_merged", values_from = "mean_count") %>%
  left_join(., taxonomy_table_glassslide_tax, by = "OTUID")

write_tsv(glassslide_average, paste0(output_dir, "glasslide_average_count_table.tsv")) #31 samples , 8 taxonomy 


#################################################################################################


#                       check the total ASV number                                              #


#################################################################################################

glassslide_average

glassslide_average_total_ASV <- glassslide_average %>% 
  select(., -Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  pivot_longer(., -OTUID, names_to = "SampleID", values_to = "count") %>%
  filter(count > 0) %>%
  add_column(num = 1) %>%
  group_by(SampleID) %>%
  summarise(asv_number = sum(num), total_counts = sum(count))

arrange(glassslide_average_total_ASV, asv_number)
#SampleID            asv_number total_counts
#<chr>                    <dbl>        <dbl>
#  1 Oct2020_Class20_B11         21        28744
#2 May2020_Class15_BF          30        28882
#3 May2020_Class15_B2          53        20792
#4 May2020_Class15_BM1         79        24415


#################################################################################################


#			Normalization 					                                                            			#


#################################################################################################

# normalization

Normalization_metagenome <- function(OTUfile = NA, col = NA, outdir = NA){
  normfile <- OTUfile[,2:col] %>%    
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    metagenomeSeq::newMRexperiment(.) %>%       ## creating MR object                      
    metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
    metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>% # the rownames will keep the same 
    as_tibble(.) %>%
    add_column(OTU = OTUfile[,1], Kingdom = OTUfile$Kingdom, Phylum = OTUfile$Phylum, 
               Class = OTUfile$Class, Order = OTUfile$Order, 
               Family = OTUfile$Family, Genus = OTUfile$Genus,
               Species = OTUfile$Species) %>%
    as.data.table() 
  write_csv(normfile, paste(outdir, "HS_glassslide_average_normalized_otu_table.csv", sep = ""))
  return(normfile)
}


col <- ncol(glassslide_average)
otu_table_average_norm <- Normalization_metagenome(OTUfile= glassslide_average, col = col, outdir = output_dir) 
col - 7
otu_table_average_norm_notaxa <- otu_table_average_norm %>%
  .[,1:32] %>%
  column_to_rownames(., var = "OTU") %>%
  rownames_to_column(., var = "OTU")

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

metadata_glassslide_unique <- metadata_combine_all %>% select(Class, sample_id_merged) %>% unique() #31samples 
alpha_table_glass_slide <- alpha_diversity(otu_table = otu_table_average_norm_notaxa, metatable_sel = metadata_glassslide_unique, name = "OTUtable")

write_csv(alpha_table_glass_slide, paste(output_dir, "glassslide_average_alpha_table.csv",sep = "/"))



#################################################################################################


#                      alpha diversity comparison                                               #


#################################################################################################


alpha_diversity_statistic_shannon <- rstatix::kruskal_test(alpha_table_glass_slide, shannon ~ Class)
alpha_diversity_statistic_simpson <- rstatix::kruskal_test(alpha_table_glass_slide, simpson ~ Class)
alpha_diversity_statistic_chao1 <- rstatix::kruskal_test(alpha_table_glass_slide, chao1 ~ Class)

alpha_diversity_statistic_ks_all <- rbind(alpha_diversity_statistic_shannon,
                                          alpha_diversity_statistic_simpson,
                                          alpha_diversity_statistic_chao1)

alpha_diversity_statistic_ks_all

# A tibble: 3 x 6
#.y.         n statistic    df      p method        
#<chr>   <int>     <dbl> <int>  <dbl> <chr>         
#  1 shannon    31      8.07     2 0.0177 Kruskal-Wallis
#2 simpson    31      7.81     2 0.0201 Kruskal-Wallis
#3 chao1      31      8.18     2 0.0168 Kruskal-Wallis

alpha_diversity_dunn_shannon <- rstatix::dunn_test(alpha_table_glass_slide, shannon ~ Class, p.adjust.method = "fdr")
alpha_diversity_dunn_simpson <- rstatix::dunn_test(alpha_table_glass_slide, simpson ~ Class, p.adjust.method = "fdr")
alpha_diversity_dunn_chao1 <- rstatix::dunn_test(alpha_table_glass_slide, chao1 ~ Class, p.adjust.method = "fdr")

alpha_diversity_dunn_all <- rbind(alpha_diversity_dunn_shannon,
                                  alpha_diversity_dunn_simpson,
                                  alpha_diversity_dunn_chao1)

alpha_diversity_dunn_all

write_tsv(alpha_diversity_dunn_all, paste0(output_dir, "/glass_slide_dunn_test.tsv"))
write_tsv(alpha_diversity_statistic_ks_all, paste0(output_dir, "/glass_slide_ks_test.tsv"))

#################################################################################################


#                      boxplot                                                                  #


#################################################################################################


# bee plot 

alpha_table_bee <- alpha_table_glass_slide
alpha_table_bee$Group <- as.character(alpha_table_glass_slide$Class)
alpha_table_bee$Group[alpha_table_bee$Class == "16"] <- "GS"
alpha_table_bee$Group[alpha_table_bee$Class == "15"] <- "GS+HS"
alpha_table_bee$Group[alpha_table_bee$Class == "20"] <- "GS+HS*"

color <- c("#079F07","#068406","#056A05") #011A01
color <- viridis::viridis(10)[c(4,6,10)]

alpha_table_bee$Group <- factor(alpha_table_bee$Group, levels = c("GS", "GS+HS", "GS+HS*"))

#########################  chao1 
chao1_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(chao1 ~ Group, p.adjust.method = "fdr") 

chao1_pvalue_statistics <- chao1_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)

chao1_pvalue_statistics_sel <- chao1_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
chao1_pvalue_statistics_sel$p.adj <- format(chao1_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)
chao1_pvalue_statistics_sel$y.position <- c(720, 780)




library(ggbeeswarm)   
chao1_plot <- ggboxplot(alpha_table_bee, x = "Group", y = "chao1", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means(label.y = 840) +
  geom_beeswarm(cex = 1.5) +
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
  )

chao1_plot

##### shannon 
shannon_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(shannon ~ Group, p.adjust.method = "fdr") 

shannon_pvalue_statistics <- shannon_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
#shannon_pvalue_statistics$p.adj <- format(shannon_pvalue_statistics$p.adj, scientific = TRUE, digits = 3)

shannon_pvalue_statistics_sel <- shannon_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
shannon_pvalue_statistics_sel$p.adj <- format(shannon_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)
shannon_pvalue_statistics_sel$y.position
#[1] 7.44248 7.88848 8.06688 8.15608
shannon_pvalue_statistics_sel$y.position <- c(6.2, 6.4)

shannon_plot <- ggboxplot(alpha_table_bee, x = "Group", y = "shannon", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means(label.y = 6.6) +
  geom_beeswarm(cex = 1.5) +
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
  )

shannon_plot

simpson_pvalue_statistics <- alpha_table_bee %>%
  rstatix::dunn_test(simpson ~ Group, p.adjust.method = "fdr") 

simpson_pvalue_statistics <- simpson_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)



simpson_pvalue_statistics_sel <- simpson_pvalue_statistics %>% dplyr::filter(., p.adj <= 0.05)
simpson_pvalue_statistics_sel$p.adj <- format(simpson_pvalue_statistics_sel$p.adj, scientific = TRUE, digits = 3)
simpson_pvalue_statistics_sel$y.position <- c(1.001,1.005)


simpson_plot <- ggboxplot(alpha_table_bee, x = "Group", y = "simpson", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics_sel, tip.length = 0, label = "p.adj") +
  stat_compare_means(label.y = 1.01) +
  geom_beeswarm(cex = 1.5) +
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
  )

simpson_plot


ggsave(paste0(output_dir, "glassslide_shannon.pdf"), shannon_plot)
ggsave(paste0(output_dir, "glassslide_simpson.pdf"), simpson_plot)
ggsave(paste0(output_dir, "glassslide_chao1.pdf"), chao1_plot)



combine_alpha_diversity <- ggpubr::ggarrange(shannon_plot, simpson_plot, chao1_plot, ncol = 3)
ggsave(paste0(output_dir, "/HE_series_alpha_diversity_combine.pdf"), combine_alpha_diversity, width = 12)


saveRDS(shannon_pvalue_statistics, paste0(output_dir, "/shannon_statistics.rds"))
saveRDS(simpson_pvalue_statistics, paste0(output_dir, "/simpson_statistics.rds"))
saveRDS(chao1_pvalue_statistics, paste0(output_dir, "/chao1_statistics.rds"))
save.image(paste0(output_dir, "/glassslide_alpha_diversity.RData"))

