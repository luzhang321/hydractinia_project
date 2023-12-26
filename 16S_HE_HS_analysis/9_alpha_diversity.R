# aim of this script is to do alpha diversity analysis


# steps :
# 1. I will add the class 23 class 33 and class 32 samples 
# 2. second, I will combine the otu table together and calculate average count 
# 3. combine with taxonomy & (don't remove any unknow, i can put them all into unknown.)
# 4. normalization
# 5. do alpha diversity comparison(shannon & simpson )


# first read the batch from the previous analysis class1,2,3,4 metadata 
# the file is copied from local : 
# '/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/metatable_sel_biorep.csv'

# /sbidata/lzhang/201911_hydractinia/202201_16S_HS_analysis/9_alpha_diversity/metatable_sel_biorep.csv

library(tidyverse)
library(ggpubr)
library(rstatix)

metadata_sel_average <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/metatable_sel_biorep.csv") %>%
  filter(., Class %in% c(1,2,3,4)) #73


# 1. add new class to the metadata 

HS_metadata <- read_tsv("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/HS_metadata.xlsx")

HS_metadata_class22_33_32 <- HS_metadata %>% filter(., Class %in% c("23","32","33")) %>%
  add_column(Group = ifelse(.$Organism == "Hydractinia", "HS", "HS_Shell"), .before = "sample") %>%
  add_column(Batch = .$SampleID %>% str_remove_all(., "\\.") %>% substr(., 1,6), .before = "sample") %>%
  add_column(sub_class1 = NA, sub_class2 = NA, .before = "sample") %>%
  add_column(sampling = .$SampleID %>% str_remove_all(., "\\.") %>% substr(., 1,7), .before = "sample") %>%
  select(-sample) %>%
  add_column(uniquebiorep = paste(.$sampling, .$Class, .$Bio_Rep, sep = "_"))#.$Group, 
colnames(HS_metadata_class22_33_32)[14] <- "Egg_release/Start"
colnames(HS_metadata_class22_33_32)[8] <- "Crab_presence"

metadata_sel_average$uniquebiorep <- paste(metadata_sel_average$uniquebiorep, metadata_sel_average$Group,sep = "_")

HS_final_metadata <- rbind(metadata_sel_average, HS_metadata_class22_33_32) #95 SAMPLES 


# 2. I will combine the otu table together and calculate average count 

otu_table <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/8_export/8_export/table-with-taxonomy-header-modify.tsv"
# read file 

read_table_fun <- function(otu_table){
  otu <- read.table(otu_table, stringsAsFactors = F, sep = "\t", header = T)
  colnames(otu) <- str_remove_all(colnames(otu),"X")
  taxonomy <- str_split(otu$taxonomy,"; ",simplify = T)
  colnames(taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  otu <- otu %>%
    select(-taxonomy) %>%
    cbind(., taxonomy, stringsAsFactors = F)
  return(otu)
}

otu <- read_table_fun(otu_table)
colnames(otu)[1] <- "OTUID"

# otu - taxonomy 

otu_tax <- select(otu, OTUID, Kingdom, Phylum, Class,
                  Order, Family, Genus, Species)

write_csv(otu_tax, "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/otu_tax.csv")
colnames(otu_tax)[1] <- "OTU"


# average OTU count table and do normalization 
otu_table_average <- otu %>% 
  select(., -Kingdom, -Phylum, -Class,
         -Order, -Family, -Genus, -Species) %>%
  mutate(OTU = OTUID) %>%
  select(.,-OTUID) %>%
  pivot_longer(., -OTU, names_to = "SampleID", values_to = "count") %>%
  left_join(., HS_final_metadata, by = "SampleID") %>%
  group_by(uniquebiorep, OTU) %>%
  summarize(mean_count = round(mean(count, na.rm = TRUE))) %>%
  pivot_wider(names_from = "uniquebiorep", values_from = "mean_count") %>%
  full_join(., otu_tax, by = "OTU")


write_csv(otu_table_average, "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/average_count_table.csv")

# normalization
library(metagenomeSeq)

Normalization_metagenome <- function(OTUfile, col, outdir){
  normfile <- OTUfile[,2:col] %>%    
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    metagenomeSeq::newMRexperiment(.) %>%       ## creating MR object                      
    metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
    metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>%
    as_tibble(.) %>%
    add_column(OTU = OTUfile$OTU, Kingdom = OTUfile$Kingdom, Phylum = OTUfile$Phylum, 
               Class = OTUfile$Class, Order = OTUfile$Order, 
               Family = OTUfile$Family, Genus = OTUfile$Genus,
               Species = OTUfile$Species) %>%
    as_tibble() 
  write_csv(normfile, paste(outdir, "average_normalized_otu_table.csv", sep = ""))
  return(normfile)
}

outdir <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/"




col <- ncol(otu_table_average)
otu_table_average_norm <- Normalization_metagenome(otu_table_average, col, outdir) %>%
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  select(OTU, everything())


write_delim(otu_table_average_norm, paste(outdir, "normalized_otu_table.txt", sep = "/"), delim = "\t")
saveRDS(otu_table_average_norm, paste(outdir, "normalized_otu_table.rds", sep = "/"))

# alpha diversity 
alpha_diversity <- function(otu_table, metatable_sel, name){
  
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
  alpha_df <- tibble(shannon = otu_shannon, simpson = otu_simpson, chao1 = otu_chao1) %>%
    add_column(uniquebiorep = names(otu_shannon), name = name) %>% 
    left_join(., metatable_sel, by = "uniquebiorep", sort = FALSE)
  
  return(alpha_df)
  
}

metatable_sel_biorep_unique <- HS_final_metadata %>%
  select(., "Location_Collection","Group","Crab_presence",
         "uniquebiorep","Class") %>%
  unique(.)

write_csv(metatable_sel_biorep_unique, paste(outdir, "metatable_biorep_unique.csv", sep = "/"))


alpha_table <- alpha_diversity(otu_table_average_norm, metatable_sel_biorep_unique, "OTUtable")

write_csv(alpha_table, paste(outdir, "10_average_alpha_table.csv",sep = "/"))

# separate alpha table 
# HS vs Control : Class23+class33 vs Class32
alpha_table_HS <- filter(alpha_table, Class %in% c(23,33,32))

# HS vs HE : 23+33 vs 1,2,3,4
alpha_table_HS_HE <- filter(alpha_table, Class %in% c(1,2,3,4,23,33))

lm(shannon ~ Group, data = alpha_table_HS)

# alpha diversity statistics : for the crab presence part, cause HS samples they are majority with crab, and shell are nocrab, adjust this parameter will eluminate the true difference.
# So I don't adjust this parameter 
# ==================================================================================================================

check_residual_normality <- function(alpha_table_control){
  
  # normality assumptions 
  model_shannon  <- lm(shannon ~ Group, data = alpha_table_control) # I can't add the unique biorepID
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_shannon))
  # Compute Shapiro-Wilk test of normality
  shapiro.test(residuals(model_shannon))
  
  # simpson 
  model_simpson  <- lm(simpson ~ Group, data = alpha_table_control)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_simpson))
  # Compute Shapiro-Wilk test of normality
  shapiro.test(residuals(model_simpson))
  
  # chao1 
  model_chao1  <- lm(chao1 ~ Group, data = alpha_table_control)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_chao1))
  # Compute Shapiro-Wilk test of normality
  shapiro.test(residuals(model_chao1))
  
  normality_result <- list()
  normality_result[[1]] <- model_shannon
  normality_result[[2]] <- ggqqplot(residuals(model_shannon))
  normality_result[[3]] <- shapiro.test(residuals(model_shannon))
  normality_result[[4]] <- model_simpson
  normality_result[[5]] <- ggqqplot(residuals(model_simpson))  
  normality_result[[6]] <- shapiro.test(residuals(model_simpson))
  normality_result[[7]] <- model_chao1
  normality_result[[8]] <- ggqqplot(residuals(model_chao1))
  normality_result[[9]] <- shapiro.test(residuals(model_chao1))  
  return(normality_result)
}


alpha_table_HS_normality <- check_residual_normality(alpha_table_HS) 
alpha_table_HS_HE_normality <- check_residual_normality(alpha_table_HS_HE) 


# simpson doesn't fit lm for alpha_table_HS_normality
# alpha_table_HS_HE_normality : all fit linear 

# the statistics comparisons 
# I will use rfit for simpson comparison in HS comparisons 



# change to rfit for shannon and simpson ; using lm for chao1 
#================================================================================================================================


alpha_diverstiy_comparison <- function(alpha_table_control, para = "yes", group = NA){
  # chao1 
  model_chao1  <- lm(chao1 ~ Group, data = alpha_table_control)
  summary_model_chao1 <- summary(model_chao1)
  model_chao1_sta <- broom::tidy(model_chao1)
  pvalue_chao1 <- model_chao1_sta %>% filter(., term == group) %>% .$p.value
  rsquare_chao1 <- summary_model_chao1$r.squared # http://www.r-tutor.com/elementary-statistics/simple-linear-regression/coefficient-determination
  
  
  # shannon and simpson 
  model_shannon  <- lm(shannon ~ Group, data = alpha_table_control)
  summary_model_shannon <- summary(model_shannon)
  
  model_shannon_sta <- broom::tidy(model_shannon)
  pvalue_shannon <- model_shannon_sta %>% filter(., term == group) %>% .$p.value
  rsquare_shannon <- summary_model_shannon$r.squared
  
  
  
  if (para == "yes"){
    model_simpson  <- lm(simpson ~ Group, data = alpha_table_control)
    summary_model_simpson <- summary(model_simpson)
    
    model_simpson_sta <- broom::tidy(model_simpson)
    pvalue_simpson <- model_simpson_sta %>% filter(., term == group) %>% .$p.value
    rsquare_simpson <- summary_model_simpson$r.squared
    
    
  }else{
    simpson_wil <- wilcox.test(simpson ~ Group, data = alpha_table_control)
    # https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/#effect-size-1
    pvalue_simpson <- simpson_wil$p.value
    simpson_summary <- alpha_table_control %>% wilcox_effsize(simpson ~ Group)
    rsquare_simpson <- simpson_summary$effsize[[1]]
    
    
  }
  
  
  # produce table 
  pvalue_df <- tibble(shannon = c(pvalue_shannon, rsquare_shannon),
                          simpson = c(pvalue_simpson, rsquare_simpson),
                          chao1 = c(pvalue_chao1, rsquare_chao1))
  pvalue_df$category <- c("pvalue", "effect_size")
  
  return(pvalue_df)
  
}

HS_Ctrl_comparison <- alpha_diverstiy_comparison(alpha_table_HS, para = "no", group = "GroupHS_Shell") # cause simpson doesn't fit normal distribution, then i need to use rfit instead 
HS_HE_comparsion <- alpha_diverstiy_comparison(alpha_table_HS_HE, para = "yes", group = "GroupHS")



## draw the plot 

# "jco" 
library("scales")
# show_col(pal_jco("default")(10))
# #0073C2FF BLUE; #EFC000FF YELLOW ; #868686FF GRAY; #CD534CFF LIGHT RED 

# plot for simpson 

alpha_boxplot <- function(alpha_table = NA, key_word = NA, color_break = NA, color_value = NA, color_label = NA, pvalue_table = NA, title = NA){
  alpha_table_sel <- alpha_table %>% select(., Group, key_word)
  pvalue_sel <- pvalue_table %>% select(key_word)
  colnames(alpha_table_sel)[2] <- "ID"
  options(digits = 3) 
  p_boxplot <- alpha_table_sel %>%
    ggplot( aes(x=Group, y=ID, fill=Group)) +
    geom_boxplot() +
    scale_fill_manual(breaks = color_break, values = color_value, labels = color_label) +
    scale_x_discrete(labels = color_label) +
    theme_bw() + 
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.6, alpha=0.5) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(title) + # , "Pvlaue :", round(pvalue_sel[1,], 3), "\n", "Effect size:", round(pvalue_sel[2,], 3), sep = "")
    xlab("") + ylab(key_word) +
    annotate("text", x=2, y = max(alpha_table_sel[,2])*1, label= paste("Pvlaue :", format(pvalue_sel[[1]][1], scientific = T), "\n", "Effect size:", round(pvalue_sel[[1]][2], 4), sep = "")) + 
    rremove("legend")  
    
  return(p_boxplot)
}


# draw plot for hs vs control 

color_break <- c("HS", "HS_Shell")
color_value <- c("#CD534CFF", "#868686FF")
color_label <- c("Hydractinia symbiolongicarpus", "Control")

title <- "Hydractinia Symbiolongicarpus VS Control"

p_HS_shannon <- alpha_boxplot(alpha_table = alpha_table_HS, key_word = "shannon", color_break = color_break, color_value = color_value, color_label = color_label,
              pvalue_table = HS_Ctrl_comparison, title = title)
p_HS_simpson <- alpha_boxplot(alpha_table = alpha_table_HS, key_word = "simpson", color_break = color_break, color_value = color_value, color_label = color_label,
              pvalue_table = HS_Ctrl_comparison, title = title)
p_HS_chao1 <- alpha_boxplot(alpha_table = alpha_table_HS, key_word = "chao1", color_break = color_break, color_value = color_value, color_label = color_label,
              pvalue_table = HS_Ctrl_comparison, title = title)

# draw plot for hs vs he 

color_break <- c("H", "HS")
color_value <- c("#EFC000FF", "#CD534CFF")
color_label <- c("hydractinia echinata", "Hydractinia symbiolongicarpus")

title <- "hydractinia Echinata VS Hydractinia Symbiolongicarpus"

p_HS_HE_shannon <- alpha_boxplot(alpha_table = alpha_table_HS_HE, key_word = "shannon", color_break = color_break, color_value = color_value, color_label = color_label,
                              pvalue_table = HS_HE_comparsion, title = title)
p_HS_HE_simpson <- alpha_boxplot(alpha_table = alpha_table_HS_HE, key_word = "simpson", color_break = color_break, color_value = color_value, color_label = color_label,
                              pvalue_table = HS_HE_comparsion, title = title)
p_HS_HE_chao1 <- alpha_boxplot(alpha_table = alpha_table_HS_HE, key_word = "chao1", color_break = color_break, color_value = color_value, color_label = color_label,
                            pvalue_table = HS_HE_comparsion, title = title)

## below not modified ##

pdf(paste(outdir, "alpha-diversity_HSvscontrol.pdf", sep = "/"), width = 12, height = 5)
gridExtra::grid.arrange(p_HS_shannon, p_HS_simpson, p_HS_chao1, ncol = 3)
dev.off()  

pdf(paste(outdir, "alpha-diversity_HSvsHE.pdf", sep = "/"), width = 12, height = 5)
gridExtra::grid.arrange(p_HS_HE_shannon, p_HS_HE_simpson, p_HS_HE_chao1, ncol=3)
dev.off()  


save.image(file = paste(outdir, "alpha-diversity.RData", sep = "/"))


## 2022, july modification of the plots 
outdir <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/9_alphadiversity/"
load(paste(outdir, "alpha-diversity.RData", sep = "/"))
# suggested color : see this webpage https://colorbrewer2.org/#type=qualitative&scheme=Accent&n=6 



library(ggbeeswarm)
#install.packages("ggprism")
library(ggprism)

alpha_table_control <- alpha_table_HS
pvalue_table <- HS_Ctrl_comparison
color_label <- c("#386cb0", "#f0027f")
H_label <- "HS+SH"
control_label <- "SH"

alphadiversity_HS_control_beeplot <- alpha_boxplot_advanced_update(alpha_table_control =  alpha_table_HS, pvalue_table = HS_Ctrl_comparison, color_label = color_label, control_label = control_label,
                              H_label = H_label)

alpha_boxplot_advanced_update <- function(alpha_table_control = NA, pvalue_table = NA, color_label = NA, control_label = NA, H_label = NA){
 
  alpha_table_control$Group <- factor(alpha_table_control$Group, levels = c("HS", "HS_Shell"))
  color <- color_label
  # manually add p value 
  shannon_pvalue_statistics <- alpha_table_control %>%
    rstatix::wilcox_test(shannon ~ Group) 
  shannon_pvalue_statistics$p <- format(pvalue_table$shannon[1], scientific = TRUE, digits = 3)
  shannon_pvalue_statistics <- shannon_pvalue_statistics %>% 
    rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
  two.means <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    shannon_pvalue_statistics$group1,    shannon_pvalue_statistics$group2,    shannon_pvalue_statistics$p, max(alpha_table_control$shannon)
  )
  #shannon_pvalue_statistics$y.position <- c(7.73,7.90,8.08)
  
  shannon_HS <- alpha_table_control %>%
    ggboxplot(x="Group", y="shannon", fill = "Group") +
    #geom_boxplot(outlier.shape = NA, color=Group) +
    geom_beeswarm(cex = 3) + 
    add_pvalue(two.means, tip.length = 0) + 
    #stat_pvalue_manual(shannon_pvalue_statistics, tip.length = 0, label = "p") + 
    scale_fill_manual(breaks = c("HS_Shell", "HS"), 
                      values = color_label,
                      labels = c("HS_Shell", "HS")) +
    scale_x_discrete(labels=c("HS" = H_label, "HS_Shell" = control_label)) +# change x-axis label
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    xlab("") 

  # manually add p value 
  simpson_pvalue_statistics <- alpha_table_control %>%
    rstatix::wilcox_test(simpson ~ Group) 
  simpson_pvalue_statistics$p <- format(pvalue_table$simpson[1], scientific = TRUE, digits = 3)
  simpson_pvalue_statistics <- simpson_pvalue_statistics %>% 
    rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
  two.means.simpson <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    simpson_pvalue_statistics$group1,    simpson_pvalue_statistics$group2,    simpson_pvalue_statistics$p, max(alpha_table_control$simpson)
  )
  
  simpson_HS <- alpha_table_control %>%
    ggboxplot(x="Group", y="simpson", fill = "Group") +
    geom_beeswarm(cex = 3) + 
    add_pvalue(two.means.simpson, tip.length = 0) + 
    scale_fill_manual(breaks = c("HS_Shell", "HS"), 
                      values = color_label,
                      labels = c("HS_Shell", "HS")) +
    scale_x_discrete(labels=c("HS" = H_label, "HS_Shell" = control_label)) +# change x-axis label
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    xlab("") 
  
  
  
  # manually add p value 
  chao1_pvalue_statistics <- alpha_table_control %>%
    rstatix::wilcox_test(chao1 ~ Group) 
  chao1_pvalue_statistics$p <- format(pvalue_table$chao1[1], scientific = TRUE, digits = 3)
  chao1_pvalue_statistics <- chao1_pvalue_statistics %>% 
    rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
  two.means.chao1 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    chao1_pvalue_statistics$group1,    chao1_pvalue_statistics$group2,    chao1_pvalue_statistics$p, max(alpha_table_control$chao1)
  )
  
  chao1_HS <- alpha_table_control %>%
    ggboxplot(x = "Group", y = "chao1", fill = "Group") +
    geom_beeswarm(cex = 3) + 
    add_pvalue(two.means.chao1, tip.length = 0) + 
    scale_fill_manual(breaks = c("HS_Shell", "HS"), 
                      values = color_label,
                      labels = c("HS_Shell", "HS")) +
    scale_x_discrete(labels=c("HS" = H_label, "HS_Shell" = control_label)) +# change x-axis label
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    xlab("") 
  
  p <- ggarrange(shannon_HS, simpson_HS, chao1_HS, ncol = 3, legend = NULL)
  return(p)
  
}

# save it to plot 

pdf(paste(outdir, "alpha-diversity_HSvscontrol_update.pdf", sep = "/"), width = 10, height = 4)
alphadiversity_HS_control_beeplot
dev.off()  

save.image(file = paste(outdir, "alpha-diversity-plot-modification-202207.RData", sep = "/"))








############# not run below ############################
# later for sig.genus, i need filtering 
#10% prevelance 
otu_table_average_norm_10pre <- otu_table_average_norm %>%
  column_to_rownames(., var = "OTU") %>%
  .[rowSums(. > 0) >= 0.1*62,] %>%
  rownames_to_column(., var = "OTU")
#write_csv(otu_table_average_norm_10pre, paste(outdir, "average_normalized_otu_table_10per_prevelance.csv", sep = "")) # not used this table for diversity


