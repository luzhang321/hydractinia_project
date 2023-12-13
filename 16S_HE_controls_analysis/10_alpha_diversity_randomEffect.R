# treat techinical replicates as random effect 
# 2020/12/17  


library(tidyverse)
library(fs)
library(vegan)
library(data.table)
library(ggpubr)
library(rstatix)
library(Rfit)
library(viridis)
library(hrbrthemes)
library(nlme)


# list files 
setwd("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/random_effect/")

workdir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/random_effect/"


#inputdir <- "/Volumes/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output/normalization_part/"
inputdir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output/normalization_part/"


files <- dir_ls(inputdir, glob = "*percent*normalized.csv",recurse = TRUE, type = "file")



# metatable standarization 
# ================================================================================================================================
metatable_sel <- read_csv("/media/lu//Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/metatable_sel.csv") 
colnames(metatable_sel)[8] <- "Crab_presence"


metatable_sel$Crab_presence

table(metatable_sel$Crab_presence)
metatable_sel$Crab_presence <- ifelse(metatable_sel$Crab_presence %in% c("Yes", "yes", "YES"), 
                                      "CrabYes", metatable_sel$Crab_presence) 
metatable_sel$Crab_presence <- ifelse(metatable_sel$Crab_presence %in% c("NO", "Stone", "-","No"), 
                                      "CrabNo", metatable_sel$Crab_presence) 


table(metatable_sel$Location_Collection) # this is already standard 

# metatable add unique idd 
sampling_time_id <- readxl::read_excel("/media/lu//Lucy/documents/201911_hydractina/metadata/Untitled 2.xlsx") %>%
  .[,3:4]

metatable_sel_biorep <- left_join(metatable_sel, sampling_time_id, by = "SampleID") 
metatable_sel_biorep$uniquebiorep <- paste(metatable_sel_biorep$sampling, 
                                           metatable_sel_biorep$Class, 
                                           metatable_sel_biorep$Bio_Rep, sep = "_")


write_csv(metatable_sel_biorep, "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/random_effect/metatable_sel_biorep.csv")

# read file 

read_otu_table <- function(input){
  otu_table <- read_csv(input)
}

otu_table <- read_otu_table(files[5]) # 5 is the OTU table 

# alpha-diverstiy 
# ==================================================================================================================

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
    add_column(SampleID = names(otu_shannon), name = name) %>% 
    left_join(., metatable_sel, by = "SampleID", sort = FALSE)
  
  return(alpha_df)
  
}


# alpha diversity statistics 
# ==================================================================================================================

alpha_table <- alpha_diversity(otu_table, metatable_sel_biorep, "OTUtable")


# separate alpha table 
# control1 : shell without hydractinia 
alpha_table_control1 <- filter(alpha_table, Class %in% c(1,2,3,4,5,9,10))

# control2 : sand-> class12 
alpha_table_control2 <- filter(alpha_table, Class %in% c(1,2,3,4,5,12))

# control3 : stone -> class13 
alpha_table_control3 <- filter(alpha_table, Class %in% c(1,2,3,4,5,13))



# using linear mixed model for duplicate sampleID) 
# ==================================================================================================================

lmm_compare <- function(alpha_table){
 
  mix_model_simpson <- nlme::lme(
    fixed = simpson ~ Group + Crab_presence + Location_Collection,
    data = alpha_table,
    random = ~ 1|uniquebiorep
  )
  anova(mix_model_simpson)
  ggqqplot(residuals(mix_model_simpson))
  shapiro_test(residuals(mix_model_simpson))
  
  
  
  mix_model_shannon <- nlme::lme(
    fixed = shannon ~ Group + Crab_presence + Location_Collection,
    data = alpha_table,
    random = ~ 1|uniquebiorep
  )
  anova(mix_model_shannon)
  ggqqplot(residuals(mix_model_shannon))
  shapiro_test(residuals(mix_model_shannon))
  
  
  mix_model_chao1 <- nlme::lme(
    fixed = chao1 ~ Group + Crab_presence + Location_Collection,
    data = alpha_table,
    random = ~ 1|uniquebiorep
  )
  anova(mix_model_chao1)
  ggqqplot(residuals(mix_model_chao1))
  shapiro_test(residuals(mix_model_chao1))
  
  result <- list()
  result["model_simpson"] <- mix_model_simpson 
  result["model_shannon"] <- mix_model_shannon
  result["model_chao1"] <- mix_model_chao1 
  result["anova_simpson"] <- anova(mix_model_simpson) 
  result["anova_shannon"] <- anova(mix_model_shannon) 
  result["anova_chao1"] <- anova(mix_model_chao1) 
  result["qqplot_simpson"] <- ggqqplot(residuals(mix_model_simpson))
  result["qqplot_shannon"] <- ggqqplot(residuals(mix_model_shannon)) 
  result["qqplot_chao1"] <- ggqqplot(residuals(mix_model_chao1)) 
  result["shapiro_simpson"] <- shapiro_test(residuals(mix_model_simpson)) 
  result["shapiro_shannon"] <- shapiro_test(residuals(mix_model_shannon))
  result["shapiro_chao1"] <- shapiro_test(residuals(mix_model_chao1)) 
  
  return(result)
}


alpha_diversity_sta1 <- lmm_compare(alpha_table = alpha_table_control1)
  
alpha_diversity_sta2 <- 
  
alpha_diversity_sta3 <- 


# lme4 
# re-fit model
require(lmerTest)
library(lmtest)
library(lme4)
#install.packages("lmerTest")
summary(lmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1))
#m.semTest <- lmer(Semantic.error ~ TestTime * Diagnosis + (TestTime | SubjectID),
#                  data = NamingRecovery, REML = FALSE)

m.sem <- lmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1)
coefs <- data.frame(coef(summary(m.sem)))

m.semTest <- lmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1,
                  REML = FALSE)

# get Satterthwaite-approximated degrees of freedom
coefs$df.Satt <- coef(summary(m.semTest))[, 3]
# get approximate p-values
coefs$p.Satt <- coef(summary(m.semTest))[, 5]
coefs
lme4::nlmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1)

library(nlme)
mix_model <- nlme::lme(
  fixed = chao1 ~ Group + Crab_presence + Location_Collection,
  data = alpha_table_control1,
  random = ~ 1|uniquebiorep
)
anova(mix_model)


ggqqplot(residuals(mix_model))
ggqqplot(residuals(m.sem))
shapiro_test(residuals(mix_model))

















######### below not used ##########################################################################################

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
# conclusion : only chao1 fits normality in all three comparions : 
# conclusion : control2 -> 0.0339 but plot is very good;
# conclusion : control3 -> 0.0497 but plot is very good;


# transform of log 

log_transform <- function(alpha_table_control){
  alpha_table_control$shannon <- log10(alpha_table_control$shannon)
  alpha_table_control$simpson <- log10(alpha_table_control$simpson)
  return(alpha_table_control)
}

alpha_table_control1_log <- log_transform(alpha_table_control1)
alpha_table_control2_log <- log_transform(alpha_table_control2)
alpha_table_control3_log <- log_transform(alpha_table_control3)

alpha_table_control1_log_normality <- check_residual_normality(alpha_table_control1_log) 

alpha_table_control2_log_normality <- check_residual_normality(alpha_table_control2_log) 

alpha_table_control3_log_normality <- check_residual_normality(alpha_table_control3_log) 

# conclusion : still not fit normality ; 


# change to rfit for shannon and simpson ; using lm for chao1 
#================================================================================================================================


alpha_diverstiy_comparison <- function(alpha_table_control){
  # chao1 
  model_chao1  <- lm(chao1 ~ Group + Crab_presence + Location_Collection, data = alpha_table_control)
  pvalue_chao1 <- anova(model_chao1)$`Pr(>F)`[1]
  
  # shannon and simpson 
  rfit.shannon <- rfit(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_control)
  pvalue_rfit.shannon <- summary(rfit.shannon)
  pvalue_rfit.shannon$coefficients[2,4]
  
  rfit.simpson <- rfit(simpson ~ Group + Crab_presence + Location_Collection, data = alpha_table_control)
  pvalue_rfit.simpson <- summary(rfit.simpson)
  pvalue_rfit.simpson$coefficients[2,4]
  
  # produce table 
  pvalue_df <- data.table(shannon = pvalue_rfit.shannon$coefficients[2,4],
                          simpson = pvalue_rfit.simpson$coefficients[2,4],
                          chao1 = pvalue_chao1)
  
  return(pvalue_df)
  
}

control1 <- alpha_diverstiy_comparison(alpha_table_control1)
control2 <- alpha_diverstiy_comparison(alpha_table_control2)
control3 <- alpha_diverstiy_comparison(alpha_table_control3)


alpha_boxplot <- function(alpha_table_control, pvalue_table, name){
  shannon <- alpha_table_control %>%
    ggplot( aes(x=Group, y=shannon, fill=Group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.6, alpha=0.5) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(name) +
    xlab("") +
    annotate("text", x = 1.5, y = 4, label = paste("P = ",format(pvalue_table$shannon, scientific = TRUE), sep = ""))
  
  simpson <- alpha_table_control %>%
    ggplot( aes(x=Group, y=simpson, fill=Group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.6, alpha=0.5) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    xlab("") +
    annotate("text", x = 1.5, y = 0.96, label = paste("P = ",format(pvalue_table$simpson, scientific = TRUE), sep = ""))
  
  chao1 <- alpha_table_control %>%
    ggplot( aes(x=Group, y=chao1, fill=Group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.8, alpha=0.5) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    #ggtitle("chao1") +
    xlab("") +
    annotate("text", x = 1.5, y = 5, label = paste("P = ",format(pvalue_table$chao1, scientific = TRUE,digits = 2), sep = ""))
  
  p <- ggarrange(shannon, simpson, chao1, ncol = 3)
  return(p)
  
}


a <- alpha_boxplot(alpha_table_control1, control1, "Shell as Control")
b <- alpha_boxplot(alpha_table_control2, control2)
c <- alpha_boxplot(alpha_table_control3, control3)

ggarrange(a,b,c,ncol = 1)

# what about try with the non-filtered data?
#================================================================================================================================

filedir <- "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output/Normalization_part/"

otu_norm <- read_csv(paste(filedir,"normalized_otu_table.csv",sep = "")) 
otu_norm_sel <- otu_norm[colnames(otu_table)]

otu_norm_alpha_table <- alpha_diversity(otu_norm_sel, metatable_sel_biorep, "OTUtable")

# separate alpha table 
# control1 : shell without hydractinia 
otu_norm_alpha_table_control1 <- filter(otu_norm_alpha_table, Class %in% c(1,2,3,4,5,9,10))

# control2 : sand-> class12 
otu_norm_alpha_table_control2 <- filter(otu_norm_alpha_table, Class %in% c(1,2,3,4,5,12))

# control3 : stone -> class13 
otu_norm_alpha_table_control3 <- filter(otu_norm_alpha_table, Class %in% c(1,2,3,4,5,13))

#check normality 
otu_norm_alpha_table_control1_normality <- check_residual_normality(otu_norm_alpha_table_control1) 
otu_norm_alpha_table_control2_normality <- check_residual_normality(otu_norm_alpha_table_control2) 
otu_norm_alpha_table_control3_normality <- check_residual_normality(otu_norm_alpha_table_control3) 

# it is the same, except chao1, none of the residuel fits normal distribution 


# beta diversity with adonis adjustment 
#================================================================================================================================

# separate samples 
#================================================================================================================================
otu_table_control1 <- otu_table %>%
  pivot_longer(., -OTU, names_to = "SampleID",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep) %>%
  filter(., Class %in% c(1,2,3,4,5,9,10)) %>%
  .[,1:2799] %>%
  pivot_longer(.,-SampleID, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "SampleID", values_from = "value")


metatable_control1 <- filter(metatable_sel_biorep, Class %in% c(1,2,3,4,5,9,10))

otu_table_control2 <- otu_table %>%
  pivot_longer(., -OTU, names_to = "SampleID",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep) %>%
  filter(., Class %in% c(1,2,3,4,5,12)) %>%
  .[,1:2799] %>%
  pivot_longer(.,-SampleID, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "SampleID", values_from = "value")

metatable_control2 <- filter(metatable_sel_biorep, Class %in% c(1,2,3,4,5,12))


otu_table_control3 <- otu_table %>%
  pivot_longer(., -OTU, names_to = "SampleID",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep) %>%
  filter(., Class %in% c(1,2,3,4,5,13)) %>%
  .[,1:2799] %>%
  pivot_longer(.,-SampleID, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "SampleID", values_from = "value")

metatable_control3 <- filter(metatable_sel_biorep, Class %in% c(1,2,3,4,5,13))

sta_beta_result <- function(otu_table, metatable_sel, name){
  set.seed(5)
  distance.bray <- otu_table %>%
    #column_to_rownames(., var = name) %>%
    .[,-1] %>%
    as.matrix(.) %>%
    t(.) %>%
    vegdist(.,method = 'bray')
  
  
  # beta dispersion for group 
  beta_dispersion <- betadisper(distance.bray, group = metatable_sel$Group)
  beta_dispersion_pvalue <-  permutest(beta_dispersion)
  
  
  # anosim for group 
  anosim.result <- anosim(distance.bray, metatable_sel$Group, permutations = 999) # it will be affected by beta dispersion 
  anosim.r <- anosim.result$statistic
  anosim.p <- anosim.result$signif
  
  # adonis for group and class 
  adonis_group <- adonis(distance.bray ~ Group + Crab_presence + Location_Collection + uniquebiorep, data = metatable_sel, distance = "bray")
  
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[1]
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[1]
  
  stat_beta <- list()
  stat_beta[["disp"]] <- beta_dispersion
  stat_beta[["anosim"]] <- anosim.result
  stat_beta[["adonisr"]] <- adonis_group
  stat_beta[["adonisp"]] <- adonis_group_p_value
  return(stat_beta)
}


beta_control1 <- sta_beta_result(otu_table_control1, metatable_control1, "control1")

beta_control2 <- sta_beta_result(otu_table_control2, metatable_control2, "control2")

beta_control3 <- sta_beta_result(otu_table_control3, metatable_control3, "control3")


# PCoA plot for beta diverstiy 








# ANOVA Function 

anova_alpha_function <- function(alpha_table){  # ref : https://www.datanovia.com/en/lessons/anova-in-r/
  
  # control1 : shell without hydractinia 
  alpha_table_control1 <- filter(alpha_table, Class %in% c(1,2,3,4,5,9,10))
  
  # control2 : sand-> class12 
  alpha_table_control2 <- filter(alpha_table, Class %in% c(1,2,3,4,5,12))
  
  # control3 : stone -> class13 
  alpha_table_control3 <- filter(alpha_table, Class %in% c(1,2,3,4,5,13))
  
  # simple visualization of the data considering of three factors 
  bxp1 <- boxplot_three_factors(alpha_table_control = alpha_table_control1) 
  bxp2 <- boxplot_three_factors(alpha_table_control = alpha_table_control2) 
  bxp3 <- boxplot_three_factors(alpha_table_control = alpha_table_control3) 
  
  
  
}

# ANOVA summary statistics visualization

boxplot_three_factors <- function(alpha_table_control){
  bxp_shannon <- ggboxplot(
    alpha_table_control, x = "Group", y = "shannon", 
    color = "Location_Collection", palette = "jco", facet.by = "Crab_presence"
  ) 
  bxp_simpson <- ggboxplot(
    alpha_table_control, x = "Group", y = "simpson", 
    color = "Location_Collection", palette = "jco", facet.by = "Crab_presence"
  )
  bxp_chao1 <- ggboxplot(
    alpha_table_control, x = "Group", y = "chao1", 
    color = "Location_Collection", palette = "jco", facet.by = "Crab_presence"
  ) 
  bxp_alpha <- ggarrange(bxp_shannon, bxp_simpson, bxp_chao1,ncol=3,common.legend = T)
  return(bxp_alpha)
}



# identify outliers # Do I need to remove outlier? https://www.rdocumentation.org/packages/rstatix/versions/0.6.0/topics/identify_outliers

# https://www.rdocumentation.org/packages/rstatix/versions/0.6.0/topics/identify_outliers

##############################################################################################################

a <- alpha_table_control1 %>% 
  group_by(Group, Crab_presence, Location_Collection) %>%
  identify_outliers(shannon)



# check assumptions 
####################################################################################################
# normality assumptions 
# model  <- lm(shannon ~ Group*Crab_presence*Location_Collection, data = alpha_table_control1)
# # Create a QQ plot of residuals
# ggqqplot(residuals(model))
# # Compute Shapiro-Wilk test of normality
# shapiro_test(residuals(model))
# 


check_anova_asuumptions <- function(alpha_table_control){
  # check normality 
  # check shapiro 
  ###### shannon ########################################
  table_normality_shannon <- alpha_table_control1 %>%
    group_by(Group, Crab_presence, Location_Collection) %>%
    shapiro_test(shannon)
  # check qq plot 
  qqplot_shannon <- ggqqplot(alpha_table_control1, "shannon", ggtheme = theme_bw()) +
    facet_grid(Crab_presence + Location_Collection ~ Group, labeller = "label_both")
  
  ###### simpson ########################################
  
  table_normality_simpson <- alpha_table_control1 %>%
    group_by(Group, Crab_presence, Location_Collection) %>%
    shapiro_test(simpson)
  # check qq plot 
  qqplot_simpson <- ggqqplot(alpha_table_control1, "simpson", ggtheme = theme_bw()) +
    facet_grid(Crab_presence + Location_Collection ~ Group, labeller = "label_both")
  
  
  ###### chao1 ########################################
  table_normality_chao1 <- alpha_table_control1 %>%
    group_by(Group, Crab_presence, Location_Collection) %>%
    shapiro_test(chao1)
  # check qq plot 
  qqplot_chao1 <- ggqqplot(alpha_table_control1, "chao1", ggtheme = theme_bw()) +
    facet_grid(Crab_presence + Location_Collection ~ Group, labeller = "label_both")
  
  
  # check homogeneity of variance assumption 
  # if not significant, we can assume the homogeneity of variances in the different groups 
  Levenes_test_shannon <- alpha_table_control1 %>% levene_test(shannon ~ Crab_presence*Location_Collection*Group)
  Levenes_test_simpson <- alpha_table_control1 %>% levene_test(simpson ~ Crab_presence*Location_Collection*Group)
  Levenes_test_chao1 <- alpha_table_control1 %>% levene_test(chao1 ~ Crab_presence*Location_Collection*Group)
  
  # > table(paste(alpha_table_control1$Group, alpha_table_control1$Crab_presence, alpha_table_control1$Location_Collection))
  
  #C CrabNo Ship     C CrabYes Ship      H CrabNo Ship      H CrabNo Watt H CrabYes Aquarium     H CrabYes Ship     H CrabYes Watt 
  #27                 14                 25                  6                 17                 37                  5 
  
}


# ANOVA statistical method 
#########################################################################################

three_way_anova <- function(alpha_table_control){
  res.aov1 <- alpha_table_control1 %>% anova_test(shannon ~ Crab_presence + Location_Collection + Group) %>%
    as.data.frame() %>%
    add_column(index = "shannon")
  
  res.aov2 <- alpha_table_control1 %>% anova_test(simpson ~ Crab_presence + Location_Collection + Group) %>%
    as.data.frame() %>%
    add_column(index = "simpson")
  
  res.aov3 <- alpha_table_control1 %>% anova_test(chao1 ~ Crab_presence + Location_Collection + Group) %>%
    as.data.frame() %>%
    add_column(index = "chao1")
  
  res.aov_table <- rbind(res.aov1, res.aov2, res.aov3)
  return(res.aov_table)
}

#??? 
#summary(aov(simpson ~ Crab_presence + Location_Collection + Group, data = alpha_table_control1)) why this 2 results are so different ?? 


# lm 
#############################################################################################

# check normality 

# hydractinia <- filter(alpha_table_control1, !Class %in% c(9,10)) 
# ggqqplot(hydractinia$shannon)

check_residual_normality <- function(alpha_table_control){
  # normality assumptions 
  model  <- lm(shannon ~ Group*Crab_presence*Location_Collection, data = alpha_table_control1)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model))
  
  # simpson 
  model_simpson  <- lm(simpson ~ Group*Crab_presence*Location_Collection, data = alpha_table_control1)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_simpson))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model_simpson))
  
  # chao1 
  model_chao1  <- lm(chao1 ~ Group + Crab_presence + Location_Collection, data = alpha_table_control1)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_chao1))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model_chao1))
  
}
# hydractinia <- filter(alpha_table_control1, !Class %in% c(9,10)) 
# ggqqplot(hydractinia$shannon)


# linear regression test 
lm_test <- function(alpha_table_control1){
  summary(glm(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_control1))
  summary(glm(simpson ~ Group + Crab_presence + Location_Collection, data = alpha_table_control1))
  summary(glm(chao1 ~ Group + Crab_presence + Location_Collection, data = alpha_table_control1))
}

glm.D93 <- glm(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_control1)
anova(glm.D93)
lm.D93 <- lm(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_control1)
summary(lm.D93)
summary(glm.D93)
summary(anova(glm.D93))

# only chao1 fits the glm 
# for shannon and simpson, which one i should use? 
# how to slove the problem? 
install.packages("Rfit")
library(Rfit)
rfit.shannon <- rfit(shannon ~ Group + Crab_presence + Location_Collection, data = alpha_table_control1)
summary(rfit.shannon)

anova(glm.D93, test = "F")



# adonis 

sta_beta_alpha_result <- function(otu_table, metatable_sel, name){
  set.seed(5)
  distance.bray <- otu_table %>%
    #column_to_rownames(., var = name) %>%
    .[,-1] %>%
    as.matrix(.) %>%
    t(.) %>%
    vegdist(.,method = 'bray')
  
  # beta dispersion for group 
  beta_dispersion <- betadisper(distance.bray, group = metatable_sel$Group)
  beta_dispersion_pvalue <-  permutest(beta_dispersion)
  
  
  # beta dispersion for class 
  beta_dispersion_class <- betadisper(distance.bray, group = metatable_sel$Class)
  beta_dispersion_class_pvalue <- permutest(beta_dispersion_class)
  # boxplot(beta_dispersion_class)
  
  
  # anosim for group 
  anosim.result <- anosim(distance.bray, metatable_sel$Group, permutations = 999) # it will be affected by beta dispersion 
  anosim.r <- anosim.result$statistic
  anosim.p <- anosim.result$signif
  
  # adonis for group and class 
  adonis_group <- adonis(distance.bray ~ Group + Crab_presence + Location_Collection, data = metatable_sel, distance = "bray")
  adonis_class <- adonis(distance.bray ~ Class, data = metatable_sel, distance = "bray")
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[1]
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[1]
  
  adonis_class_r_square <- adonis_class$aov.tab$R2[1]
  adonis_class_p_value <- adonis_class$aov.tab$`Pr(>F)`[1]
}




# main program 
# ===========================================================================================================================

otu_table_list <- list()

files_new <- files 
files_new[1] <- files[]


for (i in files){
  print(i)
  print(basename(i))
  print(dirname(i))
  name <- str_remove(basename(i),".csv")
  otu_table_list[[basename(i)]] <- read_otu_table(i)
  
}




# lme4 
# re-fit model
require(lmerTest)
library(lmtest)
library(lme4)
#install.packages("lmerTest")
summary(lmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1))
#m.semTest <- lmer(Semantic.error ~ TestTime * Diagnosis + (TestTime | SubjectID),
#                  data = NamingRecovery, REML = FALSE)

m.sem <- lmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1)
coefs <- data.frame(coef(summary(m.sem)))

m.semTest <- lmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1,
                  REML = FALSE)

# get Satterthwaite-approximated degrees of freedom
coefs$df.Satt <- coef(summary(m.semTest))[, 3]
# get approximate p-values
coefs$p.Satt <- coef(summary(m.semTest))[, 5]
coefs
lme4::nlmer(shannon ~ Group + Crab_presence + Location_Collection + (1|uniquebiorep), data = alpha_table_control1)

library(nlme)
mix_model <- nlme::lme(
  fixed = chao1 ~ Group + Crab_presence + Location_Collection,
  data = alpha_table_control1,
  random = ~ 1|uniquebiorep
)
anova(mix_model)


ggqqplot(residuals(mix_model))
ggqqplot(residuals(m.sem))
shapiro_test(residuals(mix_model))



