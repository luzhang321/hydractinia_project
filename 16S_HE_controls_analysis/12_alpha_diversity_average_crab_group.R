# crab presence effect on alpha diversity 
# author 2021.02.16 

library(tidyverse)

# A. read metatable 
#========================================

setwd("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/")

metatable_sel_biorep_unique <- read_csv("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/metatable_biorep_unique.csv")
# original counts table 
# otu_table_average <- read_csv("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/10_average_count_table.csv")

# B. normalization table 
#========================================

norm_table <- read_csv("/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_normalized_otu_table.csv")

###### tips : just follow the script in 10_alpha_diversity_average.R 


# C. alpha diversity 
# =========================================

alpha_table <- read_csv("../10_average_alpha_table.csv")

head(alpha_table)

# 1. separate alpha table 
# =================
# 1+3 vs 9 (shell with crab with hydractinia VS shell with crab without hydractinia)
alpha_table_with_crab <- filter(alpha_table, Class %in% c(1,3,9))

# control2 : sand-> class12 
alpha_table_without_crab <- filter(alpha_table, Class %in% c(2,4,10))


# 2. check normality 
# =================

# made a new function called check normality 

check_residual_normality <- function(alpha_table_control){
  
  # normality assumptions 
  model_shannon  <- lm(shannon ~ Group + Location_Collection, data = alpha_table_control) # I can't add the unique biorepID
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_shannon))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model_shannon))
  
  # simpson 
  model_simpson  <- lm(simpson ~ Group + Location_Collection, data = alpha_table_control)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model_simpson))
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model_simpson))
  
  # chao1 
  model_chao1  <- lm(chao1 ~ Group + Location_Collection, data = alpha_table_control)
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



alpha_table_with_crab_normality <- check_residual_normality(alpha_table_with_crab) 

rbind(alpha_table_with_crab_normality[[3]], alpha_table_with_crab_normality[[6]], alpha_table_with_crab_normality[[9]])
# variable                 statistic p.value
# <chr>                        <dbl>   <dbl>
# 1 residuals(model_shannon)     0.943 0.269  
# 2 residuals(model_simpson)     0.865 0.00977 : this plot is okay, even though not perfectly match
# 3 residuals(model_chao1)       0.993 1.00   

alpha_table_without_crab_normlity <- check_residual_normality(alpha_table_without_crab)

rbind(alpha_table_without_crab_normlity[[3]], alpha_table_without_crab_normlity[[6]], alpha_table_without_crab_normlity[[9]])

# variable                 statistic p.value
# <chr>                        <dbl>   <dbl>
#   1 residuals(model_shannon)     0.959  0.521 
# 2 residuals(model_simpson)     0.960  0.534 
# 3 residuals(model_chao1)       0.910  0.0628

pdf("#1-3-normality-with-crab.pdf")
# order : shannon; simpson; chao1 

alpha_table_with_crab_normality[[2]]
alpha_table_with_crab_normality[[5]]# non-normality 
alpha_table_with_crab_normality[[8]]

dev.off()


pdf("#1-4-normality-without-crab.pdf")
# order : shannon; simpson; chao1 

alpha_table_without_crab_normlity[[2]]
alpha_table_without_crab_normlity[[5]]
alpha_table_without_crab_normlity[[8]]

dev.off()




# 3. compare alpha diversity 
# ===========================

alpha_diverstiy_comparison <- function(alpha_table_control, para = "yes"){
  # chao1 
  #model_chao1  <- lm(chao1 ~ Group + Location_Collection, data = alpha_table_control)
  model_chao1  <- lm(chao1 ~ Location_Collection + Group, data = alpha_table_control)
  #pvalue_chao1 <- anova(model_chao1)$`Pr(>F)`[1] : this sequencial, not correct 
  model_chao1_sta <- broom::tidy(model_chao1)
  pvalue_chao1 <- model_chao1_sta %>% filter(., term == "GroupH") %>% .$p.value
  
  # shannon and simpson 
  #model_shannon  <- lm(shannon ~ Group + Location_Collection, data = alpha_table_control)
  model_shannon  <- lm(shannon ~ Location_Collection + Group, data = alpha_table_control)
  #pvalue_shannon <- anova(model_shannon)$`Pr(>F)`[1]
  model_shannon_sta <- broom::tidy(model_shannon)
  pvalue_shannon <- model_shannon_sta %>% filter(., term == "GroupH") %>% .$p.value
  
  
  if (para == "yes"){
    #model_simpson  <- lm(simpson ~ Group + Location_Collection, data = alpha_table_control)
    model_simpson  <- lm(simpson ~  Location_Collection + Group, data = alpha_table_control)
    #pvalue_simpson <- anova(model_simpson)$`Pr(>F)`[1]
    model_simpson_sta <- broom::tidy(model_simpson)
    pvalue_simpson <- model_simpson_sta %>% filter(., term == "GroupH") %>% .$p.value
    
  }else{
    
    rfit.simpson <- Rfit::rfit(simpson ~  Location_Collection + Group, data = alpha_table_control)
    pvalue_rfit.simpson <- summary(rfit.simpson)
    pvalue_simpson <- pvalue_rfit.simpson$coefficients["GroupH",4]
    
  }
  
  
  # produce table 
  pvalue_df <- data.table(shannon = pvalue_shannon,
                          simpson = pvalue_simpson,
                          chao1 = pvalue_chao1)
  
  return(pvalue_df)
  
}


alpha_table_with_crab_compare <- alpha_diverstiy_comparison(alpha_table_with_crab, para = "no")

# shannon   simpson     chao1
# 1: 0.2391163 0.1987562 0.1843356
#:  C  H 
#   5 15

alpha_table_without_crab_compare <- alpha_diverstiy_comparison(alpha_table_without_crab, para = "yes")


# shannon    simpson      chao1
# 1: 0.0273521 0.05264617 0.03483778
# C  H 
# 9 11 

# detailed statistics of alpha diversity 

alpha_diverstiy_comparison_return_statistic <- function(alpha_table_control, para = "yes"){
  
  detailed_sta <- list()
  # chao1 
  #model_chao1  <- lm(chao1 ~ Group + Location_Collection, data = alpha_table_control)
  model_chao1  <- lm(chao1 ~ Location_Collection + Group, data = alpha_table_control)
  #pvalue_chao1 <- anova(model_chao1)$`Pr(>F)`[1] : this sequencial, not correct 
  detailed_sta[[3]] <- broom::tidy(model_chao1)
  
  # shannon and simpson 
  #model_shannon  <- lm(shannon ~ Group + Location_Collection, data = alpha_table_control)
  model_shannon  <- lm(shannon ~ Location_Collection + Group, data = alpha_table_control)
  #pvalue_shannon <- anova(model_shannon)$`Pr(>F)`[1]
  detailed_sta[[1]] <- broom::tidy(model_shannon)
  
  
  if (para == "yes"){
    #model_simpson  <- lm(simpson ~ Group + Location_Collection, data = alpha_table_control)
    model_simpson  <- lm(simpson ~  Location_Collection + Group, data = alpha_table_control)
    #pvalue_simpson <- anova(model_simpson)$`Pr(>F)`[1]
    detailed_sta[[2]] <- broom::tidy(model_simpson)
    
  }else{
    
    rfit.simpson <- Rfit::rfit(simpson ~  Location_Collection + Group, data = alpha_table_control)
    pvalue_rfit.simpson <- summary(rfit.simpson)
    detailed_sta[[2]] <- pvalue_rfit.simpson
    
  }
  
  
  # produce table 
  return(detailed_sta)
  
}

alpha_sta1 <- alpha_diverstiy_comparison_return_statistic(alpha_table_with_crab, para = "no")
alpha_sta2 <- alpha_diverstiy_comparison_return_statistic(alpha_table_without_crab, para = "no")

## t value 

# t value 
alpha_sta1[[1]]$statistic[3] # shannon - group 
alpha_sta2[[1]]$statistic[3] 

alpha_sta1[[2]]$coefficients["GroupH","t.value"] # simpson - group 
alpha_sta2[[2]]$coefficients["GroupH","t.value"]

alpha_sta1[[3]]$statistic[3] # chao1 - group 
alpha_sta2[[3]]$statistic[3] 

# rsquare 

R_square_output <- function(alpha_table_r_square_input){
  
  # shannon & chao1 in this 3 alpha table 
  mA_shannon = lm(shannon ~ Group + Location_Collection, data = alpha_table_r_square_input)
  mB_shannon = lm(shannon ~ Location_Collection, data = alpha_table_r_square_input) 
  Eta_shannon = modelCompare(mB_shannon, mA_shannon) %>% .$PRE
  
  mA_chao1 = lm(chao1 ~ Group + Location_Collection, data = alpha_table_r_square_input)
  mB_chao1 = lm(chao1 ~ Location_Collection, data = alpha_table_r_square_input) 
  Eta_chao1 = modelCompare(mB_chao1, mA_chao1) %>% .$PRE
  
  # for simpson, rfit 
  mA_rfit = rfit(simpson ~ Group + Location_Collection, data = alpha_table_r_square_input) 
  mB_rfit = rfit(simpson ~ Location_Collection, data = alpha_table_r_square_input) 
  drop_test_rfit <- drop.test(mA_rfit, mB_rfit)
  Rsquare <- effectsize::F_to_eta2(drop_test_rfit$F, df = drop_test_rfit$df1, df_error = drop_test_rfit$df2) %>% .$Eta2_partial
  
  alpha_r_square_table <- tibble(
    shannon_r_square = Eta_shannon,
    simpson_r_square = Rsquare, 
    chao1_r_square = Eta_chao1
  )
  
  return(alpha_r_square_table)
}


R_square_output_with_crab <- R_square_output(alpha_table_with_crab)
R_square_output_without_crab <- R_square_output(alpha_table_without_crab)

# part r square 

Part_R_square_output <- function(alpha_table_r_square_input){
  
  # shannon & chao1 in this 3 alpha table 
  mA_shannon = lm(shannon ~ Group + Location_Collection, data = alpha_table_r_square_input)
  mB_shannon = lm(shannon ~ Location_Collection, data = alpha_table_r_square_input) 
  Eta_shannon = modelCompare(mB_shannon, mA_shannon) %>% .$DeltaR2
  
  mA_chao1 = lm(chao1 ~ Group + Location_Collection, data = alpha_table_r_square_input)
  mB_chao1 = lm(chao1 ~ Location_Collection, data = alpha_table_r_square_input) 
  Eta_chao1 = modelCompare(mB_chao1, mA_chao1) %>% .$DeltaR2 
  
  # for simpson, rfit 
  mA_rfit = rfit(simpson ~ Group + Location_Collection, data = alpha_table_r_square_input) 
  mB_rfit = rfit(simpson ~ Location_Collection, data = alpha_table_r_square_input) 
  Rsquare <- summary(mA_rfit)$R2[1,1] - summary(mB_rfit)$R2[1,1]
  
  alpha_r_square_table <- tibble(
    shannon_r_square = Eta_shannon,
    simpson_r_square = Rsquare, 
    chao1_r_square = Eta_chao1
  )
  
  return(alpha_r_square_table)
}
Part_R_square_output_with_crab <- Part_R_square_output(alpha_table_with_crab)
Part_R_square_output_without_crab <- Part_R_square_output(alpha_table_without_crab)


# 4. draw the boxplot 
# =====================

library(viridis)

alpha_boxplot <- function(alpha_table_control, pvalue_table, name, alpha_sta){
  label1 <- paste("P = ",format(pvalue_table$shannon, scientific = TRUE,digits = 2), sep = "")
  t_value_shannon <- paste("t value = ", format(alpha_sta[[1]]$statistic[3], digits = 3), sep = "")
  
  shannon <- alpha_table_control %>%
    ggplot( aes(x=Group, y=shannon, fill=Group)) +
    geom_boxplot() +
    scale_fill_manual(breaks = c("C", "H"), 
                      values=c( "#0073C2FF", "#EFC000FF"),
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "Hydractinia", "C" = "Control")) + # change x-axis label
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
  t_value_chao1 <- paste("t value = ", format(alpha_sta[[3]]$statistic[3], digits = 3), sep = "")
  
  chao1 <- alpha_table_control %>%
    ggplot( aes(x=Group, y=chao1, fill=Group)) +
    geom_boxplot() +
    scale_fill_manual(breaks = c("C", "H"), 
                      values=c( "#0073C2FF", "#EFC000FF"),
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = "Hydractinia", "C" = "Control")) + 
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.8, alpha=0.5) +
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste(name, " ",  label1, "\n", t_value_chao1, sep = "")) +
    xlab("") #+
  #annotate("text", x = 1.5, y = 5, label = paste("P = ",format(pvalue_table$chao1, scientific = TRUE,digits = 2), sep = ""))
  
  p <- ggarrange(shannon, simpson, chao1, ncol = 3)
  return(p)
  
}

alpha_boxplot_with_crab <- alpha_boxplot(alpha_table_with_crab, alpha_table_with_crab_compare, "Shell with Crab", alpha_sta = alpha_sta1)

alpha_boxplot_without_crab <- alpha_boxplot(alpha_table_without_crab, alpha_table_without_crab_compare, "Shell without Crab", alpha_sta = alpha_sta2)

pdf("#1-1.alpha_boxplot_with_crab.pdf", width = 9, height = 6)
alpha_boxplot_with_crab
dev.off()

pdf("#1-2.alpha_boxplot_without_crab.pdf", width = 9, height = 6)
alpha_boxplot_without_crab
dev.off()

####### 2022.05 update plots 


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
    scale_x_discrete(labels=c("H" = expression(italic("H. echinata")), "C" = control_label)) +# change x-axis label
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
    scale_x_discrete(labels=c("H" = expression(italic("H. echinata")), "C" = control_label)) +
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
    scale_x_discrete(labels=c("H" = expression(italic("H. echinata")), "C" = control_label)) + 
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



color_label <- c( "#ffff99", "#7fc97f")
control_label <- "S_C(Shell_with_crab)"
alpha_boxplot_with_crab_update <- alpha_boxplot_advanced(alpha_table_with_crab, alpha_table_with_crab_compare, "Shell with Crab", alpha_sta = alpha_sta1, color_label = color_label, control_label = control_label)

color_label <- c( "#ffff99", "#7fc97f")
control_label <- "S(Shell_without_crab)"
alpha_boxplot_without_crab_update <- alpha_boxplot_advanced(alpha_table_without_crab, alpha_table_without_crab_compare, "Shell without Crab", alpha_sta = alpha_sta2, color_label = color_label, control_label = control_label)

pdf("#1-1.alpha_boxplot_with_crab_update.pdf", width = 10, height = 4)
alpha_boxplot_with_crab_update
dev.off()

pdf("#1-2.alpha_boxplot_without_crab_update.pdf", width = 10, height = 4)
alpha_boxplot_without_crab_update
dev.off()

# i didn't save rdata, just use previous one and exceute the plot funciton 

################  update2 plot again, change the name & color : 2022.06.29.######################
load("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/12_alpha_beta_diversity_crab.RData")

#install.packages("ggprism")
library(ggprism)
library(ggbeeswarm)

alpha_boxplot_advanced_update <- function(alpha_table_control = NA, pvalue_table = NA, name = NA, alpha_sta = NA, color_label = NA, control_label = NA, H_label = NA){
  
  alpha_table_control$Group <- factor(alpha_table_control$Group, levels = c("H", "C"))
  color <- color_label
  # manually add p value from lm result 
  shannon_pvalue_statistics <- alpha_table_control %>%
    rstatix::wilcox_test(shannon ~ Group) 
  shannon_pvalue_statistics$p <- format(pvalue_table$shannon, scientific = TRUE, digits = 3)
  shannon_pvalue_statistics <- shannon_pvalue_statistics %>% 
    rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
  two.means <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    shannon_pvalue_statistics$group1,    shannon_pvalue_statistics$group2,    shannon_pvalue_statistics$p, max(alpha_table_control$shannon)
  )
  #shannon_pvalue_statistics$y.position <- c(7.73,7.90,8.08)
  
  shannon <- alpha_table_control %>%
    ggboxplot(x="Group", y="shannon", fill = "Group") +
    #geom_boxplot(outlier.shape = NA, color=Group) +
    geom_beeswarm(cex = 3) + 
    add_pvalue(two.means, tip.length = 0) + 
    #stat_pvalue_manual(shannon_pvalue_statistics, tip.length = 0, label = "p") + 
    scale_fill_manual(breaks = c("C", "H"), 
                      values = color_label,
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = H_label, "C" = control_label)) +# change x-axis label
    theme_bw() + 
    #geom_jitter(color="black", size=1) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    #ggtitle(paste(name, " ", label1, "\n", t_value_shannon, sep = "")) +
    xlab("") #+
    #annotate("text", x = 1.5, y = max(alpha_table_control$shannon), label = label1, fontface = 'bold.italic')
  
  label1 <- paste("p = ",format(pvalue_table$simpson, scientific = TRUE,digits = 2), sep = "")
  # manually add p value 
  simpson_pvalue_statistics <- alpha_table_control %>%
    rstatix::wilcox_test(simpson ~ Group) 
  simpson_pvalue_statistics$p <- format(pvalue_table$simpson, scientific = TRUE, digits = 2)
  simpson_pvalue_statistics <- simpson_pvalue_statistics %>% 
    rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
  two.means.simpson <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    simpson_pvalue_statistics$group1,    simpson_pvalue_statistics$group2,    simpson_pvalue_statistics$p, max(alpha_table_control$simpson)
  )
  
  simpson <- alpha_table_control %>%
    ggboxplot(x="Group", y="simpson", fill = "Group") +
    geom_beeswarm(cex = 3) + 
    add_pvalue(two.means.simpson, tip.length = 0) + 
    scale_fill_manual(breaks = c("C", "H"), 
                      values = color_label,
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = H_label, "C" = control_label)) +
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
  chao1_pvalue_statistics$p <- format(pvalue_table$chao1, scientific = TRUE, digits = 2)
  chao1_pvalue_statistics <- chao1_pvalue_statistics %>% 
    rstatix::add_xy_position(x = "Group", fun = "mean_sd", dodge = 0.8)
  two.means.chao1 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    chao1_pvalue_statistics$group1,    chao1_pvalue_statistics$group2,    chao1_pvalue_statistics$p, max(alpha_table_control$chao1)
  )
  
  label1 <- paste("p = ",format(pvalue_table$chao1, scientific = TRUE,digits = 2), sep = "")
  chao1 <- alpha_table_control %>%
    ggboxplot(x = "Group", y = "chao1", fill = "Group") +
    geom_beeswarm(cex = 3) + 
    add_pvalue(two.means.chao1, tip.length = 0) + 
    scale_fill_manual(breaks = c("C", "H"), 
                      values= color_label,
                      labels = c("Control", "Hydractinia")) +
    scale_x_discrete(labels=c("H" = H_label, "C" = control_label)) + 
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(face="bold"),
      axis.text.y = element_text(face="bold"),
      axis.title.y = element_text(face="bold")
    ) +
    xlab("") 
  
  p <- ggarrange(shannon, simpson, chao1, ncol = 3, legend = NULL)
  return(p)
  
}


#"#f5973d", "#09ba09"
color_label <- c( "#f5f50c", "#09ba09")
H_label <- "HE+SH+CR"
control_label <- "SH+CR"
alpha_boxplot_with_crab_updatev2 <- alpha_boxplot_advanced_update(alpha_table_with_crab, alpha_table_with_crab_compare, "Shell with Crab", alpha_sta = alpha_sta1, color_label = color_label, control_label = control_label,  H_label = H_label)

color_label <- c( "#f5f50c", "#09ba09")
H_label <- "HE+SH"
control_label <- "SH"
alpha_boxplot_without_crab_updatev2 <- alpha_boxplot_advanced_update(alpha_table_without_crab, alpha_table_without_crab_compare, "Shell without Crab", alpha_sta = alpha_sta2, color_label = color_label, control_label = control_label,  H_label = H_label)

pdf("#1-1.alpha_boxplot_with_crab_updatev2.pdf", width = 10, height = 4)
alpha_boxplot_with_crab_updatev2
dev.off()

pdf("#1-2.alpha_boxplot_without_crab_updatev2.pdf", width = 10, height = 4)
alpha_boxplot_without_crab_updatev2
dev.off()


save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/12_alpha_beta_diversity_crab_plot_modification.RData")




# D. beta diversity 
# ===========================================

# 1. separate samples 
#==================

otu_table_average_norm <- norm_table %>%
  .[,1:63] %>%
  column_to_rownames(., var = "OTU") %>%
  rownames_to_column(., var = "OTU")



otu_table_withcrab <- otu_table_average_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(1,3,9)) %>%
  .[,1:12738] %>% # after this column are all empty
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")


metatable_otu_table_withcrab <- filter(metatable_sel_biorep_unique, Class %in% c(1,3,9))

otu_table_withoutcrab <- otu_table_average_norm[,1:63] %>%
  pivot_longer(., -OTU, names_to = "uniquebiorep",values_to = "Value") %>%
  pivot_wider(., names_from = "OTU", values_from = "Value") %>%
  left_join(., metatable_sel_biorep_unique) %>%
  filter(., Class %in% c(2,4,10)) %>%
  .[,1:12738] %>%
  pivot_longer(.,-uniquebiorep, names_to = "OTU", values_to = "value") %>%
  pivot_wider(., names_from = "uniquebiorep", values_from = "value")

metatable_otu_table_withoutcrab <- filter(metatable_sel_biorep_unique, Class %in% c(2,4,10))


# 2. statistical of beta-diversity
#==================

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
  set.seed(5)
  anosim.result <- anosim(distance.bray, metatable_sel_order$Group, permutations = 999) # it will be affected by beta dispersion 
  anosim.r <- anosim.result$statistic
  anosim.p <- anosim.result$signif
  
  # adonis for group and class 
  set.seed(5)
  adonis_group <- adonis(distance.bray ~ Location_Collection + Group, data = metatable_sel_order, distance = "bray")
  
  
  adonis_group_r_square <- adonis_group$aov.tab$R2[2] # i modified this one, since the position changed
  adonis_group_p_value <- adonis_group$aov.tab$`Pr(>F)`[2] # i modified this one, since the position changed 
  
  
  
  stat_beta <- list()
  stat_beta[["disp"]] <- beta_dispersion_pvalue
  stat_beta[["anosim"]] <- anosim.result
  stat_beta[["adonisr"]] <- adonis_group
  stat_beta[["adonisp"]] <- adonis_group_p_value
  
  set.seed(5)
  pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  colnames(points) <- c("x", "y", "z")
  eig <- pcoa$eig
  points <- cbind(points, metatable_sel_order)
  c <- PcoA_plot(points, eig) + ggtitle(paste(name,"P:",adonis_group_p_value,"\n","R square:", 
                                              round(adonis_group_r_square,3), sep = " "))
  stat_beta[["plot"]] <- c
  
  stat_beta[["dis"]] <- distance.bray
  
  return(stat_beta)
}

PcoA_plot <- function(points_df, eig){
  
  p <- ggscatter(points_df, x = "x", y = "y",
                 color = "Group", palette = "jco",
                 shape = "Group",
                 ellipse = TRUE, ellipse.type = "convex") + 
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) 
  
  return(p)
}

beta_with_crab <- sta_beta_result(otu_table_withcrab, metatable_otu_table_withcrab, "Shell with Crab")

beta_without_crab <- sta_beta_result(otu_table_withoutcrab, metatable_otu_table_withoutcrab, "Shell without Crab")


## 2022.05 modification of the plot + 2022.06 i change the color 

# shellwithcrab as control
rsquare_shell_beta_withcrab <- round(beta_with_crab[["adonisr"]]$aov.tab$R2[2], 3)
library(ggpubr)
library(ggplot2)
beta_diversity_withcrab <- beta_with_crab[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 4) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c( "#dede18", "#09ba09"),
                     labels = c("SH+CR", "HE+SH+CR")) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c( "#dede18", "#09ba09"),
                    labels = c("SH+CR", "HE+SH+CR")) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("SH+CR", "HE+SH+CR")) +
  annotate(geom="text", x = 0.2, y = 0.3, label = bquote(P == ~ .(beta_with_crab[["adonisp"]]))) + 
  annotate(geom="text", x = 0.2, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_shell_beta_withcrab)))

# shellwithoutcrab as control
rsquare_shell_beta_withoutcrab <- round(beta_without_crab[["adonisr"]]$aov.tab$R2[2], 3)
library(ggpubr)
library(ggplot2)
beta_diversity_withoutcrab <- beta_without_crab[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 4) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c( "#dede18", "#09ba09"),
                     labels = c("SH", "HE+SH")) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c( "#dede18", "#09ba09"),
                    labels = c("SH", "HE+SH")) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("SH", "HE+SH")) +
  annotate(geom="text", x = 0.2, y = 0.3, label = bquote(P == ~ .(beta_without_crab[["adonisp"]]))) + 
  annotate(geom="text", x = 0.2, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_shell_beta_withoutcrab)))


pdf("#2-2.beta-without-crab-updatev2.pdf")
beta_diversity_withoutcrab
dev.off()

pdf("#2-1.beta-with-crab-updatev2.pdf")
beta_diversity_withcrab
dev.off()

save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/12_alpha_beta_diversity_crab_plot_modification.RData")


## 2022.11.19 modification of the plot: make the ylim and xlim bigger and also fill the color 
load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/12_alpha_beta_diversity_crab_plot_modification.RData")


beta_diversity_withoutcrabv2 <- beta_diversity_withoutcrab + theme(
  legend.position="top",
  plot.title = element_text(size=11),
  axis.text.x = element_text(face="bold"),
  axis.text.y = element_text(face="bold"),
  axis.title.y = element_text(face="bold"),
  axis.title.x = element_text(face="bold"),
  text = element_text(size=15)
) 

beta_diversity_withcrabv2 <- beta_diversity_withcrab + theme(
  legend.position="top",
  plot.title = element_text(size=11),
  axis.text.x = element_text(face="bold"),
  axis.text.y = element_text(face="bold"),
  axis.title.y = element_text(face="bold"),
  axis.title.x = element_text(face="bold"),
  text = element_text(size=15)
) 

pdf("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/#2-2.beta-without-crab-updatev3.pdf")
beta_diversity_withoutcrabv2
dev.off()

pdf("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/#2-1.beta-with-crab-updatev3.pdf")
beta_diversity_withcrabv2
dev.off()


pdf("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/#1-1.alpha_boxplot_with_crab_updatev3.pdf", width = 10, height = 4)
alpha_boxplot_with_crab_updatev3 <- alpha_boxplot_with_crab_updatev2 +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    text = element_text(size=20)
  ) 
alpha_boxplot_with_crab_updatev3
dev.off()

pdf("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/#1-2.alpha_boxplot_without_crab_updatev3.pdf", width = 10, height = 4)
alpha_boxplot_without_crab_updatev3 <- alpha_boxplot_without_crab_updatev2 + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    text = element_text(size=20)
  ) 
alpha_boxplot_without_crab_updatev3
dev.off()

# tip: for alpha diversity, i already changed the font size several months before v2 and v3 are same 

save.image("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/12_alpha_beta_diversity_crab_plot_modificationv2.RData")

# modification plot again: 2023.09.07
# change the color of the beta and alpha diversity 

load("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/12_alpha_beta_diversity_crab_plot_modificationv2.RData")

#color_label <- c( "#f5f50c", "#09ba09")
color_label <- c("#fdfd81", "#99d8c9")
H_label <- "HE+SH"
control_label <- "SH"
alpha_boxplot_without_crab_updatev4 <- alpha_boxplot_advanced_update(alpha_table_without_crab, alpha_table_without_crab_compare, "Shell without Crab", alpha_sta = alpha_sta2, color_label = color_label, control_label = control_label,  H_label = H_label)

pdf("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/#1-2.alpha_boxplot_without_crab_updatev4.pdf", width = 10, height = 4)

alpha_boxplot_without_crab_updatev4 <- alpha_boxplot_without_crab_updatev4 + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    text = element_text(size=20)
  ) 
alpha_boxplot_without_crab_updatev4
dev.off()




pdf("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/#2-2.beta-without-crab-updatev4.pdf")

beta_diversity_withoutcrabv3 <- beta_without_crab[["plot"]] + ggtitle("")+
  geom_point(aes(colour = factor(Group)), size = 4) +
  scale_color_manual(breaks = c("C", "H"), 
                     values = c( "#fdfd81", "#99d8c9"), # c("#fdfd81", "#99d8c9")
                     labels = c("SH", "HE+SH")) + 
  scale_fill_manual(breaks = c("C", "H"), 
                    values = c( "#fdfd81", "#99d8c9"),
                    labels = c("SH", "HE+SH")) + 
  scale_shape_manual(breaks = c("C", "H"), values = c(16,16), labels = c("SH", "HE+SH")) +
  annotate(geom="text", x = 0.2, y = 0.3, label = bquote(P == ~ .(beta_without_crab[["adonisp"]]))) + 
  annotate(geom="text", x = 0.2, y = 0.35, label = bquote(R^2 ==  ~ .(rsquare_shell_beta_withoutcrab))) +
  theme(
    legend.position="top",
    plot.title = element_text(size=11),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    axis.title.x = element_text(face="bold"),
    text = element_text(size=15)
  ) 
beta_diversity_withoutcrabv3
dev.off()

save.image("/Volumes/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/1_3vs9_2_4vs10/12_alpha_beta_diversity_crab_plot_modificationv3.RData")



# 3. distance ( within and between distance )

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


distance_compare <- ggarrange(produce_distance_matrix(beta_with_crab[["dis"]], metatable_otu_table_withcrab),
          produce_distance_matrix(beta_without_crab[["dis"]], metatable_otu_table_withoutcrab),
          ncol = 2, common.legend = T)


pdf("#2.beta.pdf", width = 10, height = 7)
ggarrange(beta_with_crab[[5]], beta_without_crab[[5]])
dev.off()

pdf("#2-2.beta-without-crab.pdf")
beta_without_crab[[5]]
dev.off()

pdf("#2-1.beta-with-crab.pdf")
beta_with_crab[[5]]
dev.off()

pdf("#3.distance-compare.pdf", width = 8, height = 6)
# left is with crab; right is without crab
distance_compare
dev.off()

#table 

data.table(adonis_rsquare = beta_without_crab[["adonisr"]]$aov.tab$R2[2], adonis_p = beta_without_crab[["adonisr"]]$aov.tab$`Pr(>F)`[2], disper = beta_without_crab[["disp"]]$tab$`Pr(>F)`[1], anosim_r = beta_without_crab[["anosim"]]$statistic[1], anosim_p = beta_without_crab[["anosim"]]$signif)
#    adonis_rsquare adonis_p disper  anosim_r anosim_p
#1:      0.2067709    0.001  0.003 0.6157176    0.002


data.table(adonis_rsquare = beta_with_crab[["adonisr"]]$aov.tab$R2[2], adonis_p = beta_with_crab[["adonisr"]]$aov.tab$`Pr(>F)`[2], disper = beta_with_crab[["disp"]]$tab$`Pr(>F)`[1], anosim_r = beta_with_crab[["anosim"]]$statistic[1], anosim_p = beta_with_crab[["anosim"]]$signif)

# adonis_rsquare adonis_p disper  anosim_r anosim_p
  # 1:      0.1261014    0.001   0.05 0.5642899    0.001

ggarrange(beta_with_crab[[5]], beta_without_crab[[5]])

save.image("12_alpha_beta_diversity_crab.RData")

load("12_alpha_beta_diversity_crab.RData")

