# for average table : calculate alpha-diversity


library(tidyverse)
library(fs)
library(vegan)
library(data.table)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(rstatix)
library(Rfit)

# metatable 


metatable_sel_biorep <- read_csv("/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/metatable_sel_biorep.csv")


# this is for whole count 

# read file 
metatable_sel_biorep <- read_csv("/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/metatable_sel_biorep.csv")


# this is for whole count 
workdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_from_qiime2/"

files <- dir(workdir, pattern = "*taxonomy.tsv") # i did some manual fix on the file : first delete the first line; second separate the taxonomy column

read_table_fun <- function(workdir, table){
  otu <- read.table(paste(workdir, table, sep = "/"), stringsAsFactors = F, sep = "\t", header = T)
  colnames(otu) <- str_remove_all(colnames(otu),"X")
  taxonomy <- str_split(otu$taxonomy,"; ",simplify = T)
  colnames(taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  otu <- otu %>%
    select(-taxonomy) %>%
    cbind(., taxonomy, stringsAsFactors = F)
  return(otu)
}

otu <- read_table_fun(workdir, files)

# otu - taxonomy 

otu_tax <- select(otu, OTUID, Kingdom, Phylum, Class,
                  Order, Family, Genus, Species)

write_csv(otu_tax, "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/otu_tax.csv")
colnames(otu_tax)[1] <- "OTU"


# average OTU count table and do normalization again
otu_table_average <- otu %>% 
  select(., -Kingdom, -Phylum, -Class,
         -Order, -Family, -Genus, -Species) %>%
  mutate(OTU = OTUID) %>%
  select(.,-OTUID) %>%
  #.[colnames(otu_table)] %>%
  pivot_longer(., -OTU, names_to = "SampleID", values_to = "count") %>%
  left_join(., metatable_sel_biorep, by = "SampleID") %>%
  group_by(uniquebiorep, OTU) %>%
  summarize(mean_count = round(mean(count, na.rm = TRUE))) %>%
  pivot_wider(names_from = "uniquebiorep", values_from = "mean_count") %>%
  left_join(., otu_tax, by = "OTU")
write_csv(otu_table_average, "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/10_average_count_table.csv")

# normalization

Normalization_metagenome <- function(OTUfile, col, outdir){
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
  write_csv(normfile, paste(outdir, "average_normalized_otu_table.csv", sep = ""))
  return(normfile)
}

outdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/"

col <- ncol(otu_table_average)
otu_table_average_norm <- Normalization_metagenome(otu_table_average, col, outdir) %>%
  .[,1:63] %>%
  column_to_rownames(., var = "OTU") %>%
  rownames_to_column(., var = "OTU")

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
  alpha_df <- data.table(shannon = otu_shannon, simpson = otu_simpson, chao1 = otu_chao1) %>%
    add_column(uniquebiorep = names(otu_shannon), name = name) %>% 
    left_join(., metatable_sel, by = "uniquebiorep", sort = FALSE)
  
  return(alpha_df)
  
}

metatable_sel_biorep_unique <- metatable_sel_biorep %>%
  select(., "Location_Collection","Group","Crab_presence",
         "uniquebiorep","Class") %>%
  unique(.)

write_csv(metatable_sel_biorep_unique, "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/metatable_biorep_unique.csv")
#10% prevelance 
otu_table_average_norm_10pre <- otu_table_average_norm %>%
  column_to_rownames(., var = "OTU") %>%
  .[rowSums(. > 0) >= 0.1*62,] %>%
  rownames_to_column(., var = "OTU")
write_csv(otu_table_average_norm_10pre, paste(outdir, "average_normalized_otu_table_10per_prevelance.csv", sep = "")) # not used this table 


# alpha diversity statistics 
# ==================================================================================================================

alpha_table <- alpha_diversity(otu_table_average_norm, metatable_sel_biorep_unique, "OTUtable")
alpha_table_sel <- alpha_diversity(otu_table_average_norm_10pre, metatable_sel_biorep_unique, "test")

write_csv(alpha_table, paste(outdir, "10_average_alpha_table.csv",sep = "/"))


# separate alpha table 
# control1 : shell without hydractinia 
alpha_table_control1 <- filter(alpha_table, Class %in% c(1,2,3,4,5,9,10))
#alpha_table_control1_sel <- filter(alpha_table_sel, Class %in% c(1,2,3,4,5,9,10))

# control2 : sand-> class12 
alpha_table_control2 <- filter(alpha_table, Class %in% c(1,2,3,4,5,12))
#alpha_table_control2_sel <- filter(alpha_table_sel, Class %in% c(1,2,3,4,5,12))

# control3 : stone -> class13 
alpha_table_control3 <- filter(alpha_table, Class %in% c(1,2,3,4,5,13))
#alpha_table_control3_sel <- filter(alpha_table_sel, Class %in% c(1,2,3,4,5,13))

# lm 
# ==================================================================================================================

# check normality 

# hydractinia <- filter(alpha_table_control1, !Class %in% c(9,10)) 
# ggqqplot(hydractinia$shannon)

check_residual_normality <- function(alpha_table_control){
  
  # normality assumptions 
  model_shannon  <- lm(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_control) # I can't add the unique biorepID
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_shannon))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model_shannon))
  
  # simpson 
  model_simpson  <- lm(simpson ~ Group + Crab_presence + Location_Collection, data = alpha_table_control)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_simpson))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model_simpson))
  
  # chao1 
  model_chao1  <- lm(chao1 ~ Group + Crab_presence + Location_Collection, data = alpha_table_control)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_chao1))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model_chao1))
  
  normality_result <- list()
  normality_result[[1]] <- model_shannon
  normality_result[[2]] <- ggqqplot(residuals(model_shannon))
  normality_result[[3]] <- shapiro_test(residuals(model_shannon))
  normality_result[[4]] <- model_simpson
  normality_result[[5]] <- ggqqplot(residuals(model_simpson))  
  normality_result[[6]] <- shapiro_test(residuals(model_simpson))
  normality_result[[7]] <- model_chao1
  normality_result[[8]] <- ggqqplot(residuals(model_chao1))
  normality_result[[9]] <- shapiro_test(residuals(model_chao1))  
  return(normality_result)
}


alpha_table_control1_normality <- check_residual_normality(alpha_table_control1) 
alpha_table_control2_normality <- check_residual_normality(alpha_table_control2) 
alpha_table_control3_normality <- check_residual_normality(alpha_table_control3) 

# in lm, the order doesnt affect the result 
# > model_simpson  <- lm(simpson ~ Group + Crab_presence + Location_Collection, data = alpha_table_control)
# > ggqqplot(residuals(model_simpson))
# > shapiro_test(residuals(model_simpson))
# # A tibble: 1 x 3
# variable                 statistic p.value
# <chr>                        <dbl>   <dbl>
#   1 residuals(model_simpson)     0.912 0.00198
# > model_simpson  <- lm(simpson ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
# > ggqqplot(residuals(model_simpson))
# > shapiro_test(residuals(model_simpson))





# conclusion : chao1/shannon fits normality in all three comparions : 
# conclusion : simpson the qqplot is just okay, but the shapiro_test, < 0.05
# 0.00198; 0.00813; 0.0243 ; 


# for 10% prevelance 
#alpha_table_control1_normality_sel <- check_residual_normality(alpha_table_control1_sel) 
#alpha_table_control2_normality_sel <- check_residual_normality(alpha_table_control2_sel) 
#alpha_table_control3_normality_sel <- check_residual_normality(alpha_table_control3_sel) 

# all shannon/chao1, non-sig 
# all simpson, sig ; plot not good 



# change to rfit for shannon and simpson ; using lm for chao1 
#================================================================================================================================


alpha_diverstiy_comparison <- function(alpha_table_control, para = "yes"){
  # chao1 
  model_chao1  <- lm(chao1 ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
  #model_chao1  <- lm(chao1 ~ Group + Crab_presence + Location_Collection, data = alpha_table_control)
  # after confirmation with Bastian, it would be better if I used directly summary(), in this function, 
  # order doesn't matter, anova if i put group in first position, it will produce typeI error, which is sequencial
  #pvalue_chao1 <- anova(model_chao1)$`Pr(>F)`[1]
  model_chao1_sta <- broom::tidy(model_chao1)
  pvalue_chao1 <- model_chao1_sta %>% filter(., term == "GroupH") %>% .$p.value
  #print(anova(model_chao1))
  
  # shannon and simpson 
  model_shannon  <- lm(shannon ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
  #pvalue_shannon <- anova(model_shannon)$`Pr(>F)`[1]
  #print(anova(model_shannon))
  model_shannon_sta <- broom::tidy(model_shannon)
  pvalue_shannon <- model_shannon_sta %>% filter(., term == "GroupH") %>% .$p.value
  
  if (para == "yes"){
    model_simpson  <- lm(simpson ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
    #pvalue_simpson <- anova(model_simpson)$`Pr(>F)`[1]
    #print(anova(model_simpson))
    model_simpson_sta <- broom::tidy(model_simpson)
    pvalue_simpson <- model_simpson_sta %>% filter(., term == "GroupH") %>% .$p.value
    
  }else{
    
    rfit.simpson <- rfit(simpson ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
    pvalue_rfit.simpson <- summary(rfit.simpson)
    pvalue_simpson <- pvalue_rfit.simpson$coefficients["GroupH",4]
    #print(pvalue_simpson)
  }
  
  
  # produce table 
  pvalue_df <- data.table(shannon = pvalue_shannon,
                          simpson = pvalue_simpson,
                          chao1 = pvalue_chao1)
  
  return(pvalue_df)
  
}

control1 <- alpha_diverstiy_comparison(alpha_table_control1, para = "no") # cause simpson doesn't fit normal distribution, then i need to use rfit instead 
control2 <- alpha_diverstiy_comparison(alpha_table_control2, para = "no")
control3 <- alpha_diverstiy_comparison(alpha_table_control3, para = "no")
# > alpha_diverstiy_comparison(alpha_table_control3, para = "no")
# shannon    simpson      chao1
# 1: 0.04024729 0.07554251 0.05701211


# control1_sel <- alpha_diverstiy_comparison(alpha_table_control1_sel, para = "non") # simpson not accurate 
# control2_sel <- alpha_diverstiy_comparison(alpha_table_control2_sel, para = "non") # simpson not accurate 
# control3_sel <- alpha_diverstiy_comparison(alpha_table_control3_sel, para = "non") # simpson not accurate 
alpha_diverstiy_comparison(alpha_table_control1, para = "yes")
alpha_diverstiy_comparison(alpha_table_control3, para = "yes")
alpha_diverstiy_comparison(alpha_table_control2)



# details statistics for each one 

alpha_diverstiy_comparison_return_statistic <- function(alpha_table_control, para = "yes"){
  detailed_sta <- list()
  # chao1 
  model_chao1  <- lm(chao1 ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
  
  detailed_sta[[3]] <- broom::tidy(model_chao1)
  
  # shannon and simpson 
  model_shannon  <- lm(shannon ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
  detailed_sta[[1]] <- broom::tidy(model_shannon)
  
  if (para == "yes"){
    model_simpson  <- lm(simpson ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
    detailed_sta[[2]] <- broom::tidy(model_simpson)
  }else{
    
    rfit.simpson <- rfit(chao1 ~ Crab_presence + Location_Collection + Group, data = alpha_table_control)
    pvalue_rfit.simpson <- summary(rfit.simpson)
    detailed_sta[[2]] <- pvalue_rfit.simpson
  }
  
  
  
  return(detailed_sta)
  
}


alpha_sta1 <- alpha_diverstiy_comparison_return_statistic(alpha_table_control1, para = "no")
alpha_sta2 <- alpha_diverstiy_comparison_return_statistic(alpha_table_control2, para = "no")
alpha_sta3 <- alpha_diverstiy_comparison_return_statistic(alpha_table_control3, para = "no")

# t value 
alpha_sta1[[1]]$statistic[5] # shannon - group 
alpha_sta2[[1]]$statistic[5] 
alpha_sta3[[1]]$statistic[5]
alpha_sta1[[2]]$coefficients["GroupH","t.value"] # simpson - group 
alpha_sta2[[2]]$coefficients["GroupH","t.value"]
alpha_sta3[[2]]$coefficients["GroupH","t.value"]
alpha_sta1[[3]]$statistic[5] # chao1 - group 
alpha_sta2[[3]]$statistic[5] 
alpha_sta3[[3]]$statistic[5] 

# partial r square 
R_square_output <- function(alpha_table_r_square_input){
  
  # shannon & chao1 in this 3 alpha table 
  mA_shannon = lm(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_r_square_input)
  mB_shannon = lm(shannon ~ Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  Eta_shannon = modelCompare(mB_shannon, mA_shannon) %>% .$PRE
  
  mA_chao1 = lm(chao1 ~ Group + Crab_presence + Location_Collection, data = alpha_table_r_square_input)
  mB_chao1 = lm(chao1 ~ Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  Eta_chao1 = modelCompare(mB_chao1, mA_chao1) %>% .$PRE
  
  # for simpson, rfit 
  mA_rfit = rfit(simpson ~ Group + Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  mB_rfit = rfit(simpson ~ Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  drop_test_rfit <- drop.test(mA_rfit, mB_rfit)
  Rsquare <- effectsize::F_to_eta2(drop_test_rfit$F, df = drop_test_rfit$df1, df_error = drop_test_rfit$df2) %>% .$Eta2_partial
  
  alpha_r_square_table <- tibble(
    shannon_r_square = Eta_shannon,
    simpson_r_square = Rsquare, 
    chao1_r_square = Eta_chao1
  )
  
  return(alpha_r_square_table)
}


R_square_output_control1 <- R_square_output(alpha_table_control1)
R_square_output_control2 <- R_square_output(alpha_table_control2)
R_square_output_control3 <- R_square_output(alpha_table_control3)
# > R_square_output_control1
# # A tibble: 1 x 3
# shannon_r_square simpson_r_square chao1_r_square
# <dbl>            <dbl>          <dbl>
#   1            0.134            0.157          0.122
# > R_square_output_control2
# # A tibble: 1 x 3
# shannon_r_square simpson_r_square chao1_r_square
# <dbl>            <dbl>          <dbl>
#   1            0.112            0.128         0.0647
# > R_square_output_control3
# # A tibble: 1 x 3
# shannon_r_square simpson_r_square chao1_r_square
# <dbl>            <dbl>          <dbl>
#   1            0.125            0.117          0.109


# semi-partial r square 



Part_R_square_output <- function(alpha_table_r_square_input){
  
  # shannon & chao1 in this 3 alpha table 
  mA_shannon = lm(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_r_square_input)
  mB_shannon = lm(shannon ~ Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  Eta_shannon = modelCompare(mB_shannon, mA_shannon) %>% .$DeltaR2
  
  mA_chao1 = lm(chao1 ~ Group + Crab_presence + Location_Collection, data = alpha_table_r_square_input)
  mB_chao1 = lm(chao1 ~ Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  Eta_chao1 = modelCompare(mB_chao1, mA_chao1) %>% .$DeltaR2
  
  # for simpson, rfit 
  mA_rfit = rfit(simpson ~ Group + Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  mB_rfit = rfit(simpson ~ Crab_presence + Location_Collection, data = alpha_table_r_square_input) 
  Rsquare <- summary(mA_rfit)$R2[1,1] - summary(mB_rfit)$R2[1,1]
  
  alpha_r_square_table <- tibble(
    shannon_r_square = Eta_shannon,
    simpson_r_square = Rsquare, 
    chao1_r_square = Eta_chao1
  )
  
  return(alpha_r_square_table)
}

Part_R_square_output_control1 <- Part_R_square_output(alpha_table_control1)
Part_R_square_output_control2 <- Part_R_square_output(alpha_table_control2)
Part_R_square_output_control3 <- Part_R_square_output(alpha_table_control3)



# "jco" 
library("scales")
show_col(pal_jco("default")(10))
# #0073C2FF BLUE; #EFC000FF YELLOW ; #868686FF GRAY; #CD534CFF LIGHT RED 


alpha_boxplot <- function(alpha_table_control, pvalue_table, name, alpha_sta){
  label1 <- paste("P = ",format(pvalue_table$shannon, scientific = TRUE, digits = 2), sep = "")
  t_value_shannon <- paste("t value = ", format(alpha_sta[[1]]$statistic[5], digits = 3), sep = "")
  shannon <- alpha_table_control %>%
    ggplot( aes(x=Group, y=shannon, fill=Group)) +
    geom_boxplot() +
    scale_fill_manual(breaks = c("C", "H"), 
                      values=c( "#0073C2FF", "#EFC000FF"),
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "Hydractinia", "C" = "Control")) +# change x-axis label
    theme_bw() + 
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.6, alpha=0.5) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste(name, " ", label1, "\n", t_value_shannon, sep = "")) +
    xlab("") #+
    #annotate("text", x = 1.5, y = 4, label = , scientific = TRUE), sep = ""))
  
  label1 <- paste("P = ",format(pvalue_table$simpson, scientific = TRUE,digits = 2), sep = "")
  t_value_simpson <- paste("t value = ", format(alpha_sta[[2]]$coefficients["GroupH","t.value"],digits = 3), sep = "")
  simpson <- alpha_table_control %>%
    ggplot( aes(x=Group, y=simpson, fill=Group)) +
    geom_boxplot() +
    scale_fill_manual(breaks = c("C", "H"), 
                      values=c( "#0073C2FF", "#EFC000FF"),
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "Hydractinia", "C" = "Control")) +
    theme_bw() + 
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.6, alpha=0.5) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste(name, " ",  label1, "\n", t_value_simpson, sep = "")) + 
    xlab("") 
  
  label1 <- paste("P = ",format(pvalue_table$chao1, scientific = TRUE,digits = 2), sep = "")
  t_value_chao1 <- paste("t value = ", format(alpha_sta[[3]]$statistic[5], digits = 3), sep = "")
  chao1 <- alpha_table_control %>%
    ggplot( aes(x=Group, y=chao1, fill=Group)) +
    geom_boxplot() +
    scale_fill_manual(breaks = c("C", "H"), 
                      values=c( "#0073C2FF", "#EFC000FF"),
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "Hydractinia", "C" = "Control")) + 
    theme_bw() + 
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.8, alpha=0.5) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste(name, " ",  label1, "\n", t_value_chao1, sep = "")) +
    xlab("") #+
    #annotate("text", x = 1.5, y = 5, label = paste("P = ",format(pvalue_table$chao1, scientific = TRUE,digits = 2), sep = ""))
  
  p <- ggarrange(shannon, simpson, chao1, ncol = 3, legend = NULL)
  return(p)
  
}

a <- alpha_boxplot(alpha_table_control1, control1, "Shell as Control", alpha_sta = alpha_sta1)
b <- alpha_boxplot(alpha_table_control2, control2, "Sand as Control", alpha_sta = alpha_sta2)
c <- alpha_boxplot(alpha_table_control3, control3, "Stone as Control", alpha_sta = alpha_sta3)

alpha_sig_table <- rbind(control1, control2, control3)


pdf(paste("/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/", "average_table/#10-alpha-diversity_average_boxplot.pdf", sep = ""), width = 10, height = 9)
pdf(paste(outdir, "average_table/#10-alpha-diversity_average_boxplot.pdf", sep = ""), width = 10, height = 9)
gridExtra::grid.arrange(a, b, c)
dev.off()  


#alpha_boxplot(alpha_table_control1_sel, control1_sel, "Shell as Control")
#alpha_boxplot(alpha_table_control2_sel, control2_sel, "Sand as Control")
#alpha_boxplot(alpha_table_control3_sel, control3_sel, "Stone as Control")

# change the color of the boxplot 


alpha_boxplot_advanced <- function(alpha_table_control = NA, pvalue_table = NA, name = NA, alpha_sta = NA, color_label = NA, control_label = NA){
  label1 <- paste("p = ",format(pvalue_table$shannon, scientific = TRUE, digits = 2), sep = "")
  #t_value_shannon <- paste("t value = ", format(alpha_sta[[1]]$statistic[5], digits = 3), sep = "")
  alpha_table_control$Group <- factor(alpha_table_control$Group, levels = c("H", "C"))
  color <- color_label
  shannon <- alpha_table_control %>%
    ggplot(aes(x=Group, y=shannon, fill=Group)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(breaks = c("C", "H"), 
                      values = color_label,
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "H. echinata", "C" = control_label)) +# change x-axis label
    theme_bw() + 
    geom_jitter(color="black", size=1) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    #ggtitle(paste(name, " ", label1, "\n", t_value_shannon, sep = "")) +
    xlab("") +
    annotate("text", x = 1.5, y = max(alpha_table_control$shannon), label = label1, fontface = 'bold.italic')
  
  label1 <- paste("p = ",format(pvalue_table$simpson, scientific = TRUE,digits = 2), sep = "")
  
  simpson <- alpha_table_control %>%
    ggplot( aes(x=Group, y=simpson, fill=Group)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(breaks = c("C", "H"), 
                      values=color_label,
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "H. echinata", "C" = control_label)) +
    theme_bw() + 
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=1) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    #ggtitle(paste(name, " ",  label1, "\n", t_value_simpson, sep = "")) + 
    xlab("") + 
    annotate("text", x = 1.5, y = max(alpha_table_control$simpson), label = label1, fontface = 'bold.italic')
  
  label1 <- paste("p = ",format(pvalue_table$chao1, scientific = TRUE,digits = 2), sep = "")
  chao1 <- alpha_table_control %>%
    ggplot( aes(x=Group, y=chao1, fill=Group)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(breaks = c("C", "H"), 
                      values= color_label,
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "H. echinata", "C" = control_label)) + 
    theme_bw() + 
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=1) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    #ggtitle(paste(name, " ",  label1, "\n", t_value_chao1, sep = "")) +
    xlab("") + 
    annotate("text", x = 1.5, y = max(alpha_table_control$chao1), label = label1, fontface = 'bold.italic')
  
  
  p <- ggarrange(shannon, simpson, chao1, ncol = 3, legend = NULL)
  return(p)
  
}

alpha_boxplot_advanced

color_label <- c( "#ffff99", "#7fc97f")
control_label <- "Shell"
shell_alpha_diversity <- alpha_boxplot_advanced(alpha_table_control1, control1, "Shell as Control", alpha_sta = alpha_sta1, color_label = color_label, control_label = "Shell")
color_label <- c( "#beaed4", "#7fc97f")
sand_alpha_diversity <- alpha_boxplot_advanced(alpha_table_control2, control2, "Sand as Control", alpha_sta = alpha_sta2, color_label = color_label, control_label = "Sand")
color_label <- c( "#fdc086", "#7fc97f")
stone_alpha_diversity <- alpha_boxplot_advanced(alpha_table_control3, control3, "Stone as Control", alpha_sta = alpha_sta3, color_label = color_label, control_label = "Stone")


pdf(paste("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/", "average_table/#10-alpha-diversity_average_boxplot_modified_color.pdf", sep = ""), width = 10, height = 9)
gridExtra::grid.arrange(shell_alpha_diversity, sand_alpha_diversity, stone_alpha_diversity)
dev.off() 


#7fc97f
#beaed4
#fdc086
#ffff99
#386cb0
#f0027f

## combined alpha diversity plots 
alpha_table_control1_revised <- alpha_table_control1 
alpha_table_control1_revised[alpha_table_control1_revised$Group == "C",]$Group <- "Shell"

alpha_table_control2_revised <- alpha_table_control2 
alpha_table_control2_revised[alpha_table_control2_revised$Group == "C",]$Group <- "Sand"

alpha_table_control3_revised <- alpha_table_control3
alpha_table_control3_revised[alpha_table_control3_revised$Group == "C",]$Group <- "Stone"

alpha_table_combine <- rbind(alpha_table_control1_revised, 
                             alpha_table_control2_revised, 
                             alpha_table_control3_revised) %>% unique()
alpha_table_combine$Group <- factor(alpha_table_combine$Group, levels = c("H", "Shell", "Sand", "Stone"))

shannon_pvalue_statistics <- alpha_table_combine %>%
  rstatix::wilcox_test(shannon ~ Group) %>% 
  filter(group1=="H")

shannon_pvalue_statistics <- shannon_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
shannon_pvalue_statistics$y.position <- c(7.73,7.90,8.08)

color <- c("#ffff99", "#beaed4", "#fdc086", "#7fc97f")
color <- c("#f5f50c", "#a77be3", "#f5973d", "#09ba09")
shannon_pvalue_statistics$p[shannon_pvalue_statistics$group2 == "Shell"] <- control1$shannon
shannon_pvalue_statistics$p[shannon_pvalue_statistics$group2 == "Sand"] <- control2$shannon
shannon_pvalue_statistics$p[shannon_pvalue_statistics$group2 == "Stone"] <- control3$shannon
shannon_pvalue_statistics$p <- format(shannon_pvalue_statistics$p, scientific = TRUE, digits = 2)


shannon_combine_plot <- ggboxplot(alpha_table_combine, x = "Group", y = "shannon", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics, tip.length = 0, label = "p") +
  geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
    theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

shannon_combine_plot

# simpson 
simpson_pvalue_statistics <- alpha_table_combine %>%
  rstatix::wilcox_test(simpson ~ Group) %>% 
  filter(group1=="H")

simpson_pvalue_statistics <- simpson_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
simpson_pvalue_statistics$y.position <- c(simpson_pvalue_statistics$y.position[1] - sd(simpson_pvalue_statistics$y.position), simpson_pvalue_statistics$y.position[2], simpson_pvalue_statistics$y.position[3] + sd(simpson_pvalue_statistics$y.position))

color <- c("#ffff99", "#beaed4", "#fdc086", "#7fc97f")
color <- c("#f5f50c", "#a77be3", "#f5973d", "#09ba09")
simpson_pvalue_statistics$p[simpson_pvalue_statistics$group2 == "Shell"] <- control1$simpson
simpson_pvalue_statistics$p[simpson_pvalue_statistics$group2 == "Sand"] <- control2$simpson
simpson_pvalue_statistics$p[simpson_pvalue_statistics$group2 == "Stone"] <- control3$simpson
simpson_pvalue_statistics$p <- format(simpson_pvalue_statistics$p, scientific = TRUE, digits = 2)


simpson_combine_plot <- ggboxplot(alpha_table_combine, x = "Group", y = "simpson", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics, tip.length = 0, label = "p") +
  geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

simpson_combine_plot

# chao1 
chao1_pvalue_statistics <- alpha_table_combine %>%
  rstatix::wilcox_test(chao1 ~ Group) %>% 
  filter(group1=="H")

chao1_pvalue_statistics$Groups
chao1_pvalue_statistics <- chao1_pvalue_statistics %>% 
  rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
chao1_pvalue_statistics$y.position <- c(3100, 3300, 3500)


color <- c("#ffff99", "#beaed4", "#fdc086", "#7fc97f")
color <- c("#f5f50c", "#a77be3", "#f5973d", "#09ba09")

chao1_pvalue_statistics$p[chao1_pvalue_statistics$group2 == "Shell"] <- control1$chao1
chao1_pvalue_statistics$p[chao1_pvalue_statistics$group2 == "Sand"] <- control2$chao1
chao1_pvalue_statistics$p[chao1_pvalue_statistics$group2 == "Stone"] <- control3$chao1
chao1_pvalue_statistics$p <- format(chao1_pvalue_statistics$p, scientific = TRUE, digits = 2)


chao1_combine_plot <- ggboxplot(alpha_table_combine, x = "Group", y = "chao1", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics, tip.length = 0, label = "p") +
  geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

chao1_combine_plot



pdf(paste("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/", "average_table/#10-alpha-diversity_average_boxplot_modified_color_combine.pdf", sep = ""), width = 10, height = 4)
ggarrange(shannon_combine_plot, simpson_combine_plot, chao1_combine_plot, nrow= 1)
dev.off() 
pdf(paste("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/", "average_table/#10-alpha-diversity_average_boxplot_modified_color_combinev2.pdf", sep = ""), width = 10, height = 4)
ggarrange(shannon_combine_plot, simpson_combine_plot, chao1_combine_plot, nrow= 1)
dev.off() 
save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-2022_alpha_diversity_plot_modificationv2.RData")
# 06.26 : i modified the figure group name again and change it already the code above, only excecute the plot part , and save it to pdfv2
# bee plot : 06.27 update 
load("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-2022_alpha_diversity_plot_modificationv2.RData")

alpha_table_combine_bee <- alpha_table_combine
alpha_table_combine_bee$group <- as.character(alpha_table_combine_bee$Group)
alpha_table_combine_bee$group[alpha_table_combine_bee$Class %in% c(1,3)] <- "HE+SH+CR"
alpha_table_combine_bee$group[alpha_table_combine_bee$Class %in% 5] <- "HE1+SH+CR"
alpha_table_combine_bee$group[alpha_table_combine_bee$Class %in% c(2,4)] <- "HE+SH"

alpha_table_combine_bee$group[alpha_table_combine_bee$Class %in% 9] <- "SH+CR"
alpha_table_combine_bee$group[alpha_table_combine_bee$Class %in% 10] <- "SH"
alpha_table_combine_bee$group[alpha_table_combine_bee$Class %in% 12] <- "SA"
alpha_table_combine_bee$group[alpha_table_combine_bee$Class %in% 13] <- "ST"
alpha_table_combine_bee$class <- factor(alpha_table_combine_bee$group, levels = c("HE+SH+CR", "HE+SH", "HE1+SH+CR", "CH+CR", "SH+CR", "SH", "SA", "ST"))

chao1_combine_plot
# https://www.r-statistics.com/2011/03/beeswarm-boxplot-and-plotting-it-with-r/
#install.packages("beeswarm")
library(beeswarm)
color <- c("#f5f50c", "#a77be3", "#f5973d", "#09ba09")
#alpha_table_combine_bee <- as.data.table(alpha_table_combine)
#alpha_table_combine_bee$Group <- as.character(alpha_table_combine_bee$Group)
#beeswarm(chao1 ~ Group, data = alpha_table_combine_bee,
#         pch = 19,
#         col = c("#f5f50c", "#a77be3", "#f5973d", "#09ba09"),
#         method = "swarm")
#boxplot(chao1 ~ Group, data = alpha_table_combine_bee, add = T, names = c("",""), col="#0000ff22")


#beeswarm(alpha_table_combine_bee$Group ~ alpha_table_combine_bee$chao1)
library(ggbeeswarm)   
chao1_combine_plotv2 <- ggboxplot(alpha_table_combine_bee, x = "Group", y = "chao1", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics, tip.length = 0, label = "p") +
  geom_beeswarm(cex = 3) +
  #geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

chao1_combine_plotv2

simpson_combine_plotv2 <- ggboxplot(alpha_table_combine_bee, x = "Group", y = "simpson", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics, tip.length = 0, label = "p") +
  #geom_point(aes(color = factor(Class))) + 
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  ) + 
  geom_beeswarm(cex = 3) 
simpson_combine_plotv2

shannon_combine_plotv2 <- ggboxplot(alpha_table_combine, x = "Group", y = "shannon", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics, tip.length = 0, label = "p") +
  geom_beeswarm(cex = 3) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )


pdf(paste("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/", "average_table/#10-alpha-diversity_average_boxplot_modified_color_combinev3.pdf", sep = ""), width = 13, height = 5)
ggarrange(shannon_combine_plotv2, simpson_combine_plotv2, chao1_combine_plotv2, nrow= 1)
dev.off() 

library(ggbeeswarm)   
chao1_combine_plotv3 <- ggboxplot(alpha_table_combine_bee, x = "Group", y = "chao1", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics, tip.length = 0, label = "p") +
  geom_beeswarm(cex = 3, aes(color = class), show.legend = TRUE) + 
  #geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("HE", "SH/SH+CR", "SA", "ST")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    #legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  )

chao1_combine_plotv3

simpson_combine_plotv3 <- ggboxplot(alpha_table_combine_bee, x = "Group", y = "simpson", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics, tip.length = 0, label = "p") +
  #geom_point(aes(color = factor(Class))) + 
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("HE", "SH/SH+CR", "SA", "ST")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    #legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  ) + 
  geom_beeswarm(cex = 3, aes(color = class), show.legend = TRUE) 
simpson_combine_plotv3

shannon_combine_plotv3 <- ggboxplot(alpha_table_combine_bee, x = "Group", y = "shannon", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics, tip.length = 0, label = "p") +
  geom_beeswarm(cex = 3, aes(color = class), show.legend = TRUE) + 
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("HE", "SH/SH+CR", "SA", "ST")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    #legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  ) 

shannon_combine_plotv3
pdf(paste("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/", "average_table/#10-alpha-diversity_average_boxplot_modified_color_combinev3_legend.pdf", sep = ""), width = 16, height = 4)
ggarrange(shannon_combine_plotv3, simpson_combine_plotv3, chao1_combine_plotv3, nrow= 1)
dev.off() 

save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-2022_alpha_diversity_plot_modificationv3.RData")


# 11.18 : here I should increase the font size of the figure 
load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-2022_alpha_diversity_plot_modificationv3.RData")

ggarrange(shannon_combine_plotv2, simpson_combine_plotv2, chao1_combine_plotv2, nrow = 1)

library(ggbeeswarm)   
chao1_combine_plotv2_modified <- ggboxplot(alpha_table_combine_bee, x = "Group", y = "chao1", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(chao1_pvalue_statistics, tip.length = 0, label = "p", size = 5.5) +
  geom_beeswarm(cex = 3) +
  #geom_jitter(color="black", size=1) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    text = element_text(size=20)
  )

chao1_combine_plotv2_modified


simpson_combine_plotv2_modified <- simpson_combine_plotv2 <- ggboxplot(alpha_table_combine_bee, x = "Group", y = "simpson", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(simpson_pvalue_statistics, tip.length = 0, label = "p", size = 5.5) +
  #geom_point(aes(color = factor(Class))) + 
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    text = element_text(size=20)
  ) + 
  geom_beeswarm(cex = 3) 
simpson_combine_plotv2_modified


shannon_combine_plotv2_modified <- ggboxplot(alpha_table_combine, x = "Group", y = "shannon", fill = "Group", outlier.shape = NA) +
  xlab("") + 
  stat_pvalue_manual(shannon_pvalue_statistics, tip.length = 0, label = "p", size = 5.5) +
  geom_beeswarm(cex = 3) +
  scale_fill_manual(breaks = c("Shell", "Sand", "Stone", "H"), 
                    values = color,
                    labels = c("Shell", "Sand", "Stone", "H")) +
  scale_x_discrete(labels=c("H" = "HE", "Shell" = "SH/SH+CR", "Sand" = "SA", "Stone" = "ST"))+
  theme_bw() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    text = element_text(size=20)
  ) 


shannon_combine_plotv2_modified

ggarrange(shannon_combine_plotv2_modified, simpson_combine_plotv2_modified, chao1_combine_plotv2_modified, nrow = 1)
pdf(paste("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/", "average_table/#10-alpha-diversity_average_boxplot_modified_color_combinev4.pdf", sep = ""), width = 13, height = 6)
ggarrange(shannon_combine_plotv2_modified, simpson_combine_plotv2_modified, chao1_combine_plotv2_modified, nrow= 1)
dev.off() 

save.image("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-2022_alpha_diversity_plot_modificationv4.RData")








# beta diversity 
# separate samples 
#================================================================================================================================
otu_table_control1 <- otu_table_average_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(1,2,3,4,5,9,10)) %>%
  .[,1:12738] %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")


metatable_control1 <- filter(metatable_sel_biorep_unique, Class %in% c(1,2,3,4,5,9,10))

otu_table_control2 <- otu_table_average_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(1,2,3,4,5,12)) %>%
  .[,1:12738] %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")

metatable_control2 <- filter(metatable_sel_biorep_unique, Class %in% c(1,2,3,4,5,12))


otu_table_control3 <- otu_table_average_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(1,2,3,4,5,13)) %>%
  .[,1:12738] %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")

metatable_control3 <- filter(metatable_sel_biorep_unique, Class %in% c(1,2,3,4,5,13))

# for 10% prevelance - separate samples 
# ================================================================================================================================
#otu_table_control1_sel <- otu_table_average_norm_10pre[,1:63] %>%
#  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
#  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
#  left_join(., metatable_sel_biorep_unique) %>%
#  filter(., Class %in% c(1,2,3,4,5,9,10)) %>%
#  .[,1:4531] %>%
#  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
#  pivot_wider(., names_from = "uniquebiorep", values_from = "value")



#otu_table_control2_sel <- otu_table_average_norm_10pre[,1:63] %>%
#  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
#  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
#  left_join(., metatable_sel_biorep_unique) %>%
#  filter(., Class %in% c(1,2,3,4,5,12)) %>%
#  .[,1:4531] %>%
#  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
#  pivot_wider(., names_from = "uniquebiorep", values_from = "value")


#otu_table_control3_sel <- otu_table_average_norm_10pre[,1:63] %>%
#  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
#  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
#  left_join(., metatable_sel_biorep_unique) %>%
#  filter(., Class %in% c(1,2,3,4,5,13)) %>%
#  .[,1:4531] %>%
#  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
#  pivot_wider(., names_from = "uniquebiorep", values_from = "value")


# ================================================================================================================================



sta_beta_result <- function(otu_table, metatable_sel, name){
  set.seed(5)
  distance.bray <- otu_table %>%
    #column_to_rownames(., var = name) %>%
    .[,-1] %>%
    as.matrix(.) %>%
    t(.) %>%
    vegdist(.,method = 'bray')
  
  
  metatable_sel_order <- data.table(uniquebiorep = names(otu_table)[2:length(names(otu_table))]) %>%
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
  adonis_group <- adonis(distance.bray ~ Crab_presence + Location_Collection + Group, data = metatable_sel_order, distance = "bray") # the group should be in the end position
  
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[3]#last postion is group
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[3]#last postion is group
  
  
  
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


beta_control1 <- sta_beta_result(otu_table_control1, metatable_control1, "Shell as Control")


beta_control2 <- sta_beta_result(otu_table_control2, metatable_control2, "Sand as Control")


beta_control3 <- sta_beta_result(otu_table_control3, metatable_control3, "Stone as Control")

pdf("deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity.pdf", width = 12)

ggarrange(beta_control1[[5]], beta_control2[[5]], beta_control3[[5]], ncol = 3, widths = 5, heights = 3)
#ggarrange(beta_control1_correct[[5]], beta_control2_correct[[5]], beta_control3_correct[[5]], ncol = 3, widths = 5, heights = 3)

dev.off()

## 2022.05 modification of the plot 

# shell as control
rsquare_shell_beta <- round(beta_control1[["adonisr"]]$aov.tab$R2[3], 3)
 
beta_diversity_shell <- beta_control1[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 4) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c("#ffff99", "#7fc97f"),
                     labels = c("Shell", expression(italic("H. echinata")))) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c("#ffff99", "#7fc97f"),
                    labels = c("Shell", expression(italic("H. echinata")))) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("Shell", expression(italic("H. echinata")))) +
  #annotate(geom="text", x = 0.3, y = 0.3, label = paste0("P = ", beta_control1[["adonisp"]])) +
  annotate(geom="text", x = 0.3, y = 0.3, label = bquote(P == ~ .(beta_control1[["adonisp"]]))) + 
  annotate(geom="text", x = 0.3, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_shell_beta)))
  

# sand as control
rsquare_sand_beta <- round(beta_control2[["adonisr"]]$aov.tab$R2[3], 3)

beta_diversity_sand <- beta_control2[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 4) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c("#beaed4", "#7fc97f"),
                     labels = c("Sand", expression(italic("H. echinata")))) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c("#beaed4", "#7fc97f"),
                    labels = c("Sand", expression(italic("H. echinata")))) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("Sand", expression(italic("H. echinata")))) +
  annotate(geom="text", x = 0.3, y = 0.3, label = bquote(P == ~ .(beta_control2[["adonisp"]]))) + 
  annotate(geom="text", x = 0.3, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_sand_beta)))


# stone as control 
rsquare_stone_beta <- round(beta_control3[["adonisr"]]$aov.tab$R2[3], 3)

beta_diversity_stone <- beta_control3[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 4) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c("#fdc086", "#7fc97f"),
                     labels = c("Stone", expression(italic("H. echinata")))) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c("#fdc086", "#7fc97f"),
                    labels = c("Stone", expression(italic("H. echinata")))) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("Stone", expression(italic("H. echinata")))) +
  annotate(geom="text", x = -0.4, y = 0.3, label = bquote(P == ~ .(beta_control3[["adonisp"]]))) + 
  annotate(geom="text", x = -0.4, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_stone_beta)))

pdf("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-updated.pdf", height = 4, width = 12)

ggarrange(beta_diversity_shell, beta_diversity_sand, beta_diversity_stone, ncol = 3, widths = 5, heights = 3)
#ggarrange(beta_control1_correct[[5]], beta_control2_correct[[5]], beta_control3_correct[[5]], ncol = 3, widths = 5, heights = 3)

dev.off()
save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-plots-modification.RData")

# modification again on 06.26th 
load("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-plots-modification.RData")
# shell as control
rsquare_shell_beta <- round(beta_control1[["adonisr"]]$aov.tab$R2[3], 3)
color <- c("#ffff99", "#beaed4", "#fdc086", "#7fc97f")
color <- c("#f5f50c", "#a77be3", "#f5973d", "#09ba09")

beta_diversity_shellv2 <- beta_control1[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 2.5) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c("#dede18", "#09ba09"),
                     #labels = c("Shell", expression(italic("H. echinata")))) +
                     labels = c("SH/SH+CR", "HE")) +
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c("#dede18", "#09ba09"), # #f5f50c
                    labels = c("SH/SH+CR", "HE")) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("SH/SH+CR", "HE")) +
  #annotate(geom="text", x = 0.3, y = 0.3, label = paste0("P = ", beta_control1[["adonisp"]])) +
  annotate(geom="text", x = 0.3, y = 0.3, label = bquote(P == ~ .(beta_control1[["adonisp"]]))) + 
  annotate(geom="text", x = 0.3, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_shell_beta)))


# sand as control
rsquare_sand_beta <- round(beta_control2[["adonisr"]]$aov.tab$R2[3], 3)

beta_diversity_sandv2 <- beta_control2[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 2.5) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c("#a77be3", "#09ba09"),
                     labels = c("SA", "HE")) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c("#a77be3", "#09ba09"),
                    labels = c("SA", "HE")) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("SA", "HE")) +
  annotate(geom="text", x = 0.3, y = 0.3, label = bquote(P == ~ .(beta_control2[["adonisp"]]))) + 
  annotate(geom="text", x = 0.3, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_sand_beta)))


# stone as control 
rsquare_stone_beta <- round(beta_control3[["adonisr"]]$aov.tab$R2[3], 3)

beta_diversity_stonev2 <- beta_control3[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 2.5) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c("#f5973d", "#09ba09"),
                     labels = c("ST","HE")) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c("#f5973d", "#09ba09"),
                    labels = c("ST", "HE")) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("ST", "HE")) +
  annotate(geom="text", x = -0.4, y = 0.3, label = bquote(P == ~ .(beta_control3[["adonisp"]]))) + 
  annotate(geom="text", x = -0.4, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_stone_beta)))


pdf("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-updatedv3.pdf", height = 5, width = 12)

ggarrange(beta_diversity_shellv2, beta_diversity_sandv2, beta_diversity_stonev2, ncol = 3, widths = 6, heights = 4)
#ggarrange(beta_control1_correct[[5]], beta_control2_correct[[5]], beta_control3_correct[[5]], ncol = 3, widths = 5, heights = 3)

dev.off()

save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-plots-modification-0626.RData")


## figure modification on 18, Nov 
load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-plots-modification-0626.RData")
ggarrange(beta_diversity_shellv2, beta_diversity_sandv2, beta_diversity_stonev2, ncol = 3, widths = 6, heights = 4)


beta_diversity_shellv2_modified <- beta_diversity_shellv2 + theme(
  legend.position="top",
  plot.title = element_text(size=11),
  axis.text.x = element_text(face="bold"),
  axis.text.y = element_text(face="bold"),
  axis.title.y = element_text(face="bold"),
  axis.title.x = element_text(face="bold"),
  text = element_text(size=15)
) 

beta_diversity_sandv2_modified <- beta_diversity_sandv2 + 
  theme(
    legend.position="top",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    axis.title.x = element_text(face="bold"),
    text = element_text(size=15)
  ) 


beta_diversity_stonev2_modified <- beta_diversity_stonev2 + 
  theme(
    legend.position="top",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    axis.title.x = element_text(face="bold"),
    text = element_text(size=15)
  ) 

pdf("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-updatedv4.pdf", height = 5, width = 14)

ggarrange(beta_diversity_shellv2_modified, 
          NULL,
          beta_diversity_sandv2_modified, 
          NULL,
          beta_diversity_stonev2_modified, 
          ncol = 5, widths = c(1,0.05,1,0.05, 1), heights = 4)
#ggarrange(beta_control1_correct[[5]], beta_control2_correct[[5]], beta_control3_correct[[5]], ncol = 3, widths = 5, heights = 3)
# https://github.com/kassambara/ggpubr/issues/151
dev.off()

save.image("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/#10-beta-diversity-plots-modification-1118.RData")


# below not used 
# beta_control1_sel <- sta_beta_result(otu_table_control1_sel, metatable_control1, "control1") # after I did it, the significance is not sig anymore 
# 
# beta_control2_sel <- sta_beta_result(otu_table_control2_sel, metatable_control2, "control2") # after I did it, the significance is not sig anymore 
# 
# beta_control3_sel <- sta_beta_result(otu_table_control3_sel, metatable_control3, "control3") # after I did it, the significance is not sig anymore 


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









#below not used. 
# the weird thing is that the original count table, doesn't affect the result, beta dispersion is okay.
# however result got from metagenomseq, there is a big dispersion. 
# =================================================================================================================================
dev.off()

produce_distance_matrix <- function(distance, metatable_sel){
  distance <- as.matrix(distance) 
  distance_df <- as.data.table(distance) %>%
    add_column(., UniqueID1 = rownames(distance)) %>%
    pivot_longer(., -UniqueID1, names_to = "uniquebiorep", values_to = "distance") %>%
    unique(.) %>%
    filter(UniqueID1!=uniquebiorep) %>%
    left_join(., metatable_sel) %>%
    select(., -Location_Collection, -Crab_presence, -Class) 
  
  colnames(distance_df)[2] <- "UniqueID2"
  colnames(distance_df)[4] <- "Group2"
  
  
  colnames(distance_df)[1] <- "uniquebiorep"
  
  
  distance_df <- distance_df %>%
    left_join(., metatable_sel)
  
  
  colnames(distance_df)[1] <- "UniqueID1"
  colnames(distance_df)[6] <- "Group1"
  
  within_dis <- distance_df %>%
    filter(., Group1 == "H",Group2 == "H")
  
  within_dis <- distance_df %>%
    filter(., Group1 == "C",Group2 == "C") %>%
    rbind(., within_dis)
  
  
  between_dis <- distance_df %>%
    filter(., Group1 == "H",Group2 == "C")
  
  between_dis <- distance_df %>%
    filter(., Group1 == "C",Group2 == "H") %>%
    rbind(., between_dis)
  
  dis_plot <- rbind(within_dis, between_dis) %>%
    add_column(Distype = c(rep("Within_distance", nrow(within_dis)),
                           rep("Between_distance", nrow(between_dis))))
  
  pvalue <-wilcox.test(within_dis$distance, between_dis$distance)$p.value
  
  p <- ggplot(dis_plot, aes(x=Distype, y=distance, color=Distype)) +
    geom_boxplot() +
    scale_color_brewer(palette="Dark2") +   
    ggtitle(paste("Wilcox P value :",format(pvalue,scientific = TRUE,digits = 2) , sep = " " ))
  
  # my_comparisons <- list( c("Within_distance_H", "Within_distance_C"), 
  #                         c("Within_distance_H", "Between_distance"), 
  #                         c("Within_distance_C", "Between_distance") )
  
  # p <- p + stat_compare_means(comparisons = my_comparisons)
  
  return(p)
}
# this is not correct, since distance are depedent, not correct using wilcox
#pdf("deblur_gg13_8/output_multiIV/average_table/#10-distance-between-within.pdf", width = 10)
#ggarrange(produce_distance_matrix(beta_control1[["dis"]], metatable_sel_biorep_unique),
#produce_distance_matrix(beta_control2[["dis"]], metatable_sel_biorep_unique),
#produce_distance_matrix(beta_control3[["dis"]], metatable_sel_biorep_unique),ncol = 3, common#.legend = T)

#dev.off()

# Rdata is useful 
save.image("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_alpha.RData")
load("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_alpha.RData")
save.image("/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_alpha.RData")
















# i did the sum up to see  if there is any difference (not used below )
# =================================================================================================================================

# average OTU count table and do normalization again
otu_table_sum <- otu %>% 
  select(., -Kingdom, -Phylum, -Class,
         -Order, -Family, -Genus, -Species) %>%
  mutate(OTU = OTUID) %>%
  select(.,-OTUID) %>%
  .[colnames(otu_table)] %>%
  pivot_longer(., -OTU, names_to = "SampleID", values_to = "count") %>%
  left_join(., metatable_sel_biorep, by = "SampleID") %>%
  group_by(uniquebiorep, OTU) %>%
  summarize(mean_count = round(sum(count, na.rm = TRUE))) %>%
  pivot_wider(names_from = "uniquebiorep", values_from = "mean_count") %>%
  left_join(., otu_tax, by = "OTU")


# normalization



outdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/"

col <- ncol(otu_table_sum)
otu_table_sum_norm <- Normalization_metagenome(otu_table_sum, col, outdir) %>%
  .[,1:63] %>%
  column_to_rownames(., var = "OTU") %>%
  rownames_to_column(., var = "OTU")


otu_table_control1_sum <- otu_table_sum_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(1,2,3,4,5,9,10)) %>%
  select(., -Location_Collection, -Group, -Crab_presence, -Class) %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")



otu_table_control2_sum <- otu_table_sum_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(1,2,3,4,5,12)) %>%
  select(., -Location_Collection, -Group, -Crab_presence, -Class) %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")


otu_table_control3_sum <- otu_table_sum_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(1,2,3,4,5,13)) %>%
  select(., -Location_Collection, -Group, -Crab_presence, -Class) %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")



beta_control1_sum <- sta_beta_result(otu_table_control1_sum, metatable_control1, "control1") # Group not sig 

beta_control2_sum <- sta_beta_result(otu_table_control2_sum, metatable_control2, "control2") # Group not sig 

beta_control3_sum <- sta_beta_result(otu_table_control3_sum, metatable_control3, "control3") # Group not sig 





