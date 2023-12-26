
# the aim of this script is to compare sig.taxa from HE vs Control and sig.taxa from HS vs Control and try to see whether they are some difference
# between the two
# for HE sig.taxa 
#=================================================================================================================================

load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/11_metagenomeseq.RData")

# phylum wirte : shell sand stone order 
library(tidyverse)
#otu_combined_count_sub_result[[1]][[1]] :  "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_phylum_shell.csv"
#otu_combined_count_sub_result[[1]][[2]] :  "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_phylum_sand.csv"
#otu_combined_count_sub_result[[1]][[3]] :  "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_phylum_stone.csv"



# family write 
#otu_combined_count_sub_result[[2]][[1]] : "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_family_shell.csv"
#otu_combined_count_sub_result[[2]][[2]] : "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_family_sand.csv"
#otu_combined_count_sub_result[[2]][[3]] : "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_family_stone.csv"


# genus write 

#otu_combined_count_sub_result[[3]][[1]] : "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_genus_shell.csv"
#otu_combined_count_sub_result[[3]][[2]] : "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_genus_sand.csv"
#otu_combined_count_sub_result[[3]][[3]] : "/media/lu/Lucy/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/10_sig_genus_stone.csv"


# the overlap of the HS sig.taxa 
HE_sig_phylum <- otu_combined_count_sub_result[[1]][[4]]
HE_sig_family <- otu_combined_count_sub_result[[2]][[4]]
HE_sig_genus <- otu_combined_count_sub_result[[3]][[4]]


# For HS vs Control taxa 
# ==================================================================================================================================

HS_sig_phylum <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/Taxonomy/sig_Phylum_HS_vs_control.csv")
HS_sig_Family <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/Taxonomy/sig_Family_HS_vs_control.csv")
HS_sig_Genus <- read_csv("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/Taxonomy/sig_Genus_HS_vs_control.csv")



# Compare them and draw venn plto 
# ==================================================================================================================================
library(VennDiagram)



library(RColorBrewer)
myCol <- c("cornflowerblue", "#CD534CFF")

draw_venn_plot_update_HE_HS <- function(input_HE = NA, input_HS = NA, outdir = NA, PngName = NA, main = NA){
  venn.diagram(
    x = list(input_HE, input_HS),
    category.names = c("HE vs Control" , "HS vs Control"),
    filename = paste(outdir,PngName,sep = "/"),
    output=TRUE,
    main = main,
    main.cex = 0.5,
    main.pos = c(0.5,1),
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 1,
    #lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .4,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    cat.pos = c(-10, 10),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"#
  )
}

draw_venn_plot_update_HE_HS(input_HE = HE_sig_phylum, input_HS = HS_sig_phylum$Taxa, 
                            outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                            PngName = "Phylum_sig_vennplot.png", main = "Phylum Level")

HS_sig_Family$Taxa[which(is.na(HS_sig_Family$Taxa))] <- "NA"
draw_venn_plot_update_HE_HS(input_HE = HE_sig_family, input_HS = as.character(HS_sig_Family$Taxa), 
                            outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                            PngName = "Family_sig_vennplot.png", main = "Family Level")


draw_venn_plot_update_HE_HS(input_HE = HE_sig_genus, input_HS = HS_sig_Genus$Taxa, 
                            outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                            PngName = "Genus_sig_vennplot.png", main = "Genus Level")


# extract the overlapped taxa and observed whether they are belongs to the same direction 
# =====================================================================================================================================

# 1, phylum 
overlap_HE_sig_phylum <- intersect(HE_sig_phylum, HS_sig_phylum$Taxa)
# shell control - HE 
overlap_HE_sig_phylum_shell <- otu_combined_count_sub_result[[1]][[1]] %>% 
  filter(OTU %in% overlap_HE_sig_phylum) %>% add_column(control = "HE_vs_shell")

overlap_HE_sig_phylum_sand <- otu_combined_count_sub_result[[1]][[2]] %>% 
  filter(OTU %in% overlap_HE_sig_phylum) %>% add_column(control = "HE_vs_sand")

overlap_HE_sig_phylum_stone <- otu_combined_count_sub_result[[1]][[3]] %>% 
  filter(OTU %in% overlap_HE_sig_phylum) %>% add_column(control = "HE_vs_stone")

overlap_HS_sig_phylum <- HS_sig_phylum %>% 
  filter(Taxa %in% overlap_HE_sig_phylum) %>% add_column(control = "HS_vs_control")

overlap_phylum_HE_table <- rbind(overlap_HE_sig_phylum_shell, overlap_HE_sig_phylum_sand, overlap_HE_sig_phylum_stone)

dir <- "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/"
write_csv(overlap_phylum_HE_table, 
          paste(dir, "overlap_HE_sig_phylum.csv", sep = "/"))
write_csv(overlap_HS_sig_phylum, 
          paste(dir, "overlap_HS_sig_phylum.csv", sep = "/"))

# regarding to the direction of HS and HS_control 
load("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/metagenomeseq.RData")
otu_combined_count_sub[["Phylum"]][["HS_control"]][,"p__Firmicutes"] %>% sum() # 5788, 5 samples 
otu_combined_count_sub[["Phylum"]][["HS"]][,"p__Firmicutes"] %>% sum() # 444, 10 samples 
# here by alphabetic, the group 0 is HS group 1 is control, So even though the scaling factor is negative, this means the control is higher 
overlap_HS_sig_phylum %>% filter(., Taxa == "p__Firmicutes")


# 2. family
load("/media/lu/Seagate Por/documents/201911_hydractina/202010_qiime2_16S_analysis/1_hydarctinia_microbiome/deblur_gg13_8/output_multiIV/average_table/11_metagenomeseq.RData")
overlap_HE_sig_family <- intersect(HE_sig_family, HS_sig_Family$Taxa)
# shell control - HE 
overlap_HE_sig_family_shell <- otu_combined_count_sub_result[[2]][[1]] %>% 
  filter(OTU %in% overlap_HE_sig_family) %>% add_column(control = "HE_vs_shell")

overlap_HE_sig_family_sand <- otu_combined_count_sub_result[[2]][[2]] %>% 
  filter(OTU %in% overlap_HE_sig_family) %>% add_column(control = "HE_vs_sand")

overlap_HE_sig_family_stone <- otu_combined_count_sub_result[[2]][[3]] %>% 
  filter(OTU %in% overlap_HE_sig_family) %>% add_column(control = "HE_vs_stone")
overlap_family_HE_table <- rbind(overlap_HE_sig_family_shell, overlap_HE_sig_family_sand, overlap_HE_sig_family_stone)

overlap_HS_sig_family <- HS_sig_Family %>% 
  filter(Taxa %in% overlap_HE_sig_family) %>% add_column(control = "HS_vs_control")

write_csv(overlap_family_HE_table, 
          paste(dir, "overlap_HE_sig_family.csv", sep = "/"))
write_csv(overlap_HS_sig_family, 
          paste(dir, "overlap_HS_sig_family.csv", sep = "/"))


# 3. genus 

overlap_HE_sig_genus <- intersect(HE_sig_genus, HS_sig_Genus$Taxa)
# shell control - HE 
overlap_HE_sig_genus_shell <- otu_combined_count_sub_result[[3]][[1]] %>% 
  filter(OTU %in% overlap_HE_sig_genus) %>% add_column(control = "HE_vs_shell")

overlap_HE_sig_genus_sand <- otu_combined_count_sub_result[[3]][[2]] %>% 
  filter(OTU %in% overlap_HE_sig_genus) %>% add_column(control = "HE_vs_sand")

overlap_HE_sig_genus_stone <- otu_combined_count_sub_result[[3]][[3]] %>% 
  filter(OTU %in% overlap_HE_sig_genus) %>% add_column(control = "HE_vs_stone")
overlap_genus_HE_table <- rbind(overlap_HE_sig_genus_shell, overlap_HE_sig_genus_sand, overlap_HE_sig_genus_stone)

overlap_HS_sig_genus <- HS_sig_Genus %>% 
  filter(Taxa %in% overlap_HE_sig_genus) %>% add_column(control = "HS_vs_control")

write_csv(overlap_genus_HE_table, 
          paste(dir, "overlap_HE_sig_genus.csv", sep = "/"))
write_csv(overlap_HS_sig_genus, 
          paste(dir, "overlap_HS_sig_genus.csv", sep = "/"))







save.image("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/HE_3controls_HS_control_sig_taxa_compare.RData")


############ 2022.07 figure modification only for venn plot ##############################################################################################################

load("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/HE_3controls_HS_control_sig_taxa_compare.RData")


library(VennDiagram)



library(RColorBrewer)
myCol <- c("#068f06", "#c90e70")

draw_venn_plot_update_HE_HS_v2 <- function(input_HE = NA, input_HS = NA, outdir = NA, PngName = NA, main = NA){
  p <- venn.diagram(
    x = list(input_HE, input_HS),
    category.names = c("HE vs \n SH/SH+CR/ST/SA" , "HS+SH vs SH"),
    filename=NULL, 
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5,1),
    
    # Output features
    imagetype="svg" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    # Circles
    lwd = 5,
    fill = myCol,
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    cat.pos = c(-10, 10),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
    #cat.pos = c(-5, 5),
    #cat.fontfamily = "sans",
    #margin = 0.5
  )
  
  ggsave(p, file=paste(outdir,PngName,sep = "/"), device = "svg", width = 6)
  return(p)
}

phylum_overlap <- draw_venn_plot_update_HE_HS_v2(input_HE = HE_sig_phylum, input_HS = HS_sig_phylum$Taxa, 
                            outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                            PngName = "Phylum_sig_vennplot_update.svg", main = "Phylum Level")

HS_sig_Family$Taxa[which(is.na(HS_sig_Family$Taxa))] <- "NA"
family_overlap <- draw_venn_plot_update_HE_HS_v2(input_HE = HE_sig_family, input_HS = as.character(HS_sig_Family$Taxa), 
                            outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                            PngName = "Family_sig_vennplot_update.svg", main = "Family Level")


genera_overlap <- draw_venn_plot_update_HE_HS_v2(input_HE = HE_sig_genus, input_HS = HS_sig_Genus$Taxa, 
                            outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                            PngName = "Genus_sig_vennplot_update.svg", main = "Genus Level")


save.image("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/HE_3controls_HS_control_sig_taxa_compare_202207_plotmodification.RData")


###### venn plot figure update, by using legend : Nov. 2022 


load("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/HE_3controls_HS_control_sig_taxa_compare_202207_plotmodification.RData")

library(VennDiagram)
myCol <- c("#068f06", "#c90e70")

draw_venn_plot_update_HE_HS_v3 <- function(input_HE = NA, input_HS = NA, outdir = NA, PngName = NA, main = NA, pos = NA){
  p <- venn.diagram(
    x = list(input_HE, input_HS),
    category.names = c("", ""),
    #category.names = c("HE vs \n SH/SH+CR/ST/SA" , "HS+SH vs SH"),
    filename=NULL, 
    output=TRUE,
    main = main,
    main.cex = 2,
    main.pos = c(0.5, pos),
    
    # Circles
    # Circles
    lwd = 5,
    fill = myCol,
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = c("text"),
    cat.pos = c(-10, 10),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
    #cat.pos = c(-5, 5),
    #cat.fontfamily = "sans",
    #margin = 0.5
  )
  
  
  cols <- c("#068f06", "#c90e70")
  lg <- legendGrob(labels =c("HE vs \n SH/SH+CR/ST/SA" , "HS+SH vs SH"), 
                   pch = rep(15, 2),
                   gp=gpar(col = cols, fill = "gray"),
                   byrow = TRUE)
  
  g <- gTree(children = gList(p))
  final_p <- gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4, 1))
  
  
  ggsave(final_p, file = paste(outdir, PngName, sep = "/"), device = "pdf", width = 9)  
  
  
  return(final_p)
  
  
}


phylum_overlap <- draw_venn_plot_update_HE_HS_v3(input_HE = HE_sig_phylum, input_HS = HS_sig_phylum$Taxa, 
                                                 outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                                                 PngName = "Phylum_sig_vennplot_updatev2.pdf", main = "Phylum Level", pos = 0.95)

HS_sig_Family$Taxa[which(is.na(HS_sig_Family$Taxa))] <- "NA"
family_overlap <- draw_venn_plot_update_HE_HS_v3(input_HE = HE_sig_family, input_HS = as.character(HS_sig_Family$Taxa), 
                                                 outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                                                 PngName = "Family_sig_vennplot_update2.pdf", main = "Family Level", pos = 0.95)


genera_overlap <- draw_venn_plot_update_HE_HS_v3(input_HE = HE_sig_genus, input_HS = HS_sig_Genus$Taxa, 
                                                 outdir = "/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/",
                                                 PngName = "Genus_sig_vennplot_updatev2.pdf", main = "Genus Level", pos = 0.95)


combine_venn <- ggpubr::ggarrange(phylum_overlap, family_overlap, genera_overlap, ncol = 3)
ggsave("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/combine_venn.pdf", combine_venn, width = 28)

save.image("/media/lu/Seagate Por/documents/201911_hydractina/202201_HS/11_significant_spes/compare_HSvscontrol_and_HEvs3controls_sig_taxa/HE_3controls_HS_control_sig_taxa_compare_202211_plotmodification.RData")
