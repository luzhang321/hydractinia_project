# the aim of this script is to do taxonomic distribution of the different time series 
# Lu Zhang 
# 08.08

#######################################################################################################


#                                   load library                                                      # 


#######################################################################################################

library(tidyverse)
library(ggpubr)
library(ggplot2)



input_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/9_diversity/"
load(paste0(input_dir, "9_alpha_diversity_analysis_HE_series_0909.RData"))

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile/"



#######################################################################################################


#                                   prepare input                                                    # 


#######################################################################################################

metadata_taxon <- alpha_table_add_sel # the updated samples 63 samples 
input_normalized_taxon <- otu_table_average_norm #73columns, 8taxon, 2 samples needed to be removed 


########################################################################################################


#                    combined the taxon table for each taxon                                           # 

########################################################################################################

input_normalized_taxon

taxon_table <- list()
taxon_table_long_combined <- list()
for (i in c("Phylum", "Class", "Order", "Family", "Genus")){
  print(i)
  taxon_table[[i]] <- input_normalized_taxon %>% select(metadata_taxon$sample_id_merged, all_of(i)) #63 samples 
  taxon_table_long <- taxon_table[[i]] %>% pivot_longer(., !all_of(i)) 
  colnames(taxon_table_long)[1] <- "taxon"
  
  taxon_table_long$taxon[taxon_table_long$taxon == ""] <- "unassigned"
  taxon_table_long$taxon[taxon_table_long$taxon == "p__"] <- "unassigned"
  taxon_table_long$taxon[taxon_table_long$taxon == "c__"] <- "unassigned"
  taxon_table_long$taxon[taxon_table_long$taxon == "f__"] <- "unassigned"
  taxon_table_long$taxon[taxon_table_long$taxon == "o__"] <- "unassigned"
  taxon_table_long$taxon[taxon_table_long$taxon == "g__"] <- "unassigned"
  
  keep_taxon <- taxon_table_long %>% group_by(taxon) %>% summarise(sum = sum(value)) %>% filter(sum > 0) # only keep the non-zero taxon 
  
  taxon_table_long <- taxon_table_long %>% filter(taxon %in% keep_taxon$taxon)
  
 
  taxon_table_long_combined[[i]] <- taxon_table_long %>%
    group_by(taxon, name) %>%
    summarise(abundance = sum(value))
  
  write_tsv(taxon_table_long_combined[[i]], paste0(output_dir, i, "_combined_normalized.tsv"))

}


########################################################################################################


#                    select the top10 taxon and make other as name unassigned                          # 

########################################################################################################

taxon_table_long_combined_mean <- list()
for (i in names(taxon_table_long_combined)){
  
  taxon_table <- taxon_table_long_combined[[i]]
  taxon_table_long_combined_mean[[i]] <- select_top10(taxon_table = taxon_table, input_metadata = metadata_taxon)
  
}

#input_metadata <- metadata_taxon
select_top10 <- function(taxon_table = NA, input_metadata = NA){
  
  taxon_table_combine_metadata <- taxon_table %>%
    left_join(., input_metadata, by = c("name" = "sample_id_merged")) %>% 
    group_by(taxon, Group) %>%
    summarise(mean_abundance = mean(abundance), num = n()) %>%
    ungroup()
  
  top_10 <- taxon_table %>%
    left_join(., input_metadata, by = c("name" = "sample_id_merged")) %>% 
    group_by(taxon) %>%
    summarise(mean_abundance = mean(abundance)) %>%
    arrange(desc(mean_abundance)) %>%
    .[1:10,] 
    
  
  taxon_table_combine_metadata_rename <- taxon_table_combine_metadata
  taxon_table_combine_metadata_rename$taxon[!taxon_table_combine_metadata_rename$taxon %in% top_10$taxon] <- "Others"
  
  taxon_table_combine_metadata_rename <- taxon_table_combine_metadata_rename %>% 
    group_by(taxon, Group) %>%
    summarise(mean_abundance = sum(mean_abundance)) %>%
    ungroup()
  
  result <- list()
  result[["combine_table"]] <- taxon_table_combine_metadata
  result[["taxon"]] <- top_10
  result[["combine_table_top10"]] <- taxon_table_combine_metadata_rename
  return(result)
  
}

########################################################################################################


#                    Visualization of the taxon distribution                                           # 

########################################################################################################



#taxon <- taxon_table_long_combined_mean[["Phylum"]][["combine_table_top10"]]
manualcolors <- c( ' forestgreen ' , ' red2 ' , ' orange ' , ' cornflowerblue ' ,
                 ' magenta ' , ' darkolivegreen4 ' ,
                 ' indianred1 ' , ' tan4 ' , ' darkblue ' ,
                 ' mediumorchid1 ' , ' firebrick4 ' , ' yellowgreen ' , ' lightsalmon ' , ' tan3 ' ,
                 "tan1", ' darkgray ' , ' wheat4 ' , ' #DDAD4B ' , ' chartreuse ' , ' seagreen1 ' ,
                 ' moccasin ' , ' mediumvioletred ' , ' seagreen ' , ' cadetblue1 ' ,
                 "darkolivegreen1" ,"tan2" , "tomato3" , "#7CE3D8","gainsboro","black")


draw_barplot <- function(taxon = NA, tag = NA){
  
  
  barchart1 <- ggplot(taxon,aes(x = Group, y = mean_abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual(name = tag,values = manualcolors) +
    labs(x = "", y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12)
    )
  
  return(barchart1)
  
}


barchart <- list()
for (i in names(taxon_table_long_combined_mean)){
  
  barchart[[i]] <- draw_barplot(taxon = taxon_table_long_combined_mean[[i]][["combine_table_top10"]], tag = i)
  ggsave(filename = paste0(output_dir, i, "_barchart.pdf", sep = ""), barchart[[i]], width = 8)
  
}
save.image(paste0(output_dir, "11_taxonomic_profile.RData"))

########################################################################################################


#        Visualization of the taxon distribution  - modified fixed color for fixed taxon               # 

########################################################################################################


# modification 09.23(fixed color for fixed taxon) 

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile/"
load(paste0(output_dir, "11_taxonomic_profile.RData"))


barchart <- list()
for (i in names(taxon_table_long_combined_mean)){
  
  barchart[[i]] <- draw_barplot(taxon = taxon_table_long_combined_mean[[i]][["combine_table_top10"]], tag = i)
  ggsave(filename = paste0(output_dir, i, "_barchart.pdf", sep = ""), barchart[[i]], width = 8)
  
}


# phylum manually color values=c("4" = "red", "r" = "blue", "f" = "darkgreen")


i <- "Phylum"
taxon_table_long_combined_mean[[i]][["combine_table_top10"]]$taxon %>% unique()

barchart[[i]] <- draw_barplot_phylum(taxon = taxon_table_long_combined_mean[[i]][["combine_table_top10"]], tag = i)
#[1] "Others"             "p__Acidobacteria"   "p__Actinobacteria"  "p__Bacteroidetes"   "p__Chloroflexi"     "p__Cyanobacteria"   "p__Firmicutes"      "p__Planctomycetes" 
#[9] "p__Proteobacteria"  "p__Tenericutes"     "p__Verrucomicrobia"

draw_barplot_phylum <- function(taxon = NA, tag = NA){

  barchart1 <- ggplot(taxon,aes(x = Group, y = mean_abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual(name = tag, values = c("Others" = "grey", 
                                             "p__Acidobacteria" = "forestgreen", 
                                             "p__Actinobacteria" = "red2",
                                             "p__Bacteroidetes" = "orange",
                                             "p__Chloroflexi" = "cornflowerblue",
                                             "p__Cyanobacteria" = "magenta",
                                             "p__Firmicutes" = "darkolivegreen4",
                                             "p__Planctomycetes" = "indianred1",
                                             "p__Proteobacteria" = "tan4",
                                             "p__Tenericutes" = "darkblue",
                                             "p__Verrucomicrobia" = "mediumorchid1")) +
    labs(x = "", y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12)
    )
  
  return(barchart1)
  
}



## Family 
i <- "Family"
taxon_table_long_combined_mean[[i]][["combine_table_top10"]]$taxon %>% unique()

barchart[[i]] <- draw_barplot_family(taxon = taxon_table_long_combined_mean[[i]][["combine_table_top10"]], tag = i)
#[1] "f__Alteromonadaceae"    "f__Flammeovirgaceae"    "f__Flavobacteriaceae"   "f__Hyphomicrobiaceae"   "f__Phyllobacteriaceae"  "f__Piscirickettsiaceae" "f__Pseudanabaenaceae"
#[8] "f__Rhodobacteraceae"    "f__Saprospiraceae"      "Others"                 "unassigned"          

draw_barplot_family <- function(taxon = NA, tag = NA){
  
  barchart1 <- ggplot(taxon,aes(x = Group, y = mean_abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual(name = tag, values = c("Others" = "grey", 
                                             "f__Alteromonadaceae" = "forestgreen", 
                                             "f__Flammeovirgaceae" = "red2",
                                             "f__Flavobacteriaceae" = "orange",
                                             "f__Hyphomicrobiaceae" = "cornflowerblue",
                                             "f__Phyllobacteriaceae" = "magenta",
                                             "f__Piscirickettsiaceae" = "darkolivegreen4",
                                             "f__Pseudanabaenaceae" = "indianred1",
                                             "f__Rhodobacteraceae" = "tan4",
                                             "f__Saprospiraceae" = "darkblue",
                                             "unassigned" = "mediumorchid1")) +
    labs(x = "", y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12)
    )
  
  return(barchart1)
  
}


i <- "Genus"
taxon_table_long_combined_mean[[i]][["combine_table_top10"]]$taxon %>% unique()

barchart[[i]] <- draw_barplot_genus(taxon = taxon_table_long_combined_mean[[i]][["combine_table_top10"]], tag = i)
#[1] "g__Aquimarina"     "g__Arcobacter"     "g__Loktanella"     "g__Maribacter"     "g__Mycoplasma"     "g__Octadecabacter" "g__Phaeobacter"    "g__Tenacibaculum" 
#[9] "g__Vibrio"         "Others"            "unassigned"       

draw_barplot_genus <- function(taxon = NA, tag = NA){
  
  barchart1 <- ggplot(taxon,aes(x = Group, y = mean_abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual(name = tag, values = c("Others" = "grey", 
                                             "g__Aquimarina" = "forestgreen", 
                                             "g__Arcobacter" = "red2",
                                             "g__Loktanella" = "orange",
                                             "g__Maribacter" = "cornflowerblue",
                                             "g__Mycoplasma" = "magenta",
                                             "g__Octadecabacter" = "darkolivegreen4",
                                             "g__Phaeobacter" = "indianred1",
                                             "g__Tenacibaculum" = "tan4",
                                             "g__Vibrio" = "darkblue",
                                             "unassigned" = "mediumorchid1")) +
    labs(x = "", y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12)
    )
  
  return(barchart1)
  
}

for (i in names(taxon_table_long_combined_mean)){
  
 
  ggsave(filename = paste0(output_dir, i, "_barchart.pdf", sep = ""), barchart[[i]], width = 8)
  
}
save.image(paste0(output_dir, "11_taxonomic_profile_modification0923.RData"))




########################################################################################################


#        check the biggest values of unclassified asv                                                  # 

########################################################################################################

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile/"
load(paste0(output_dir, "11_taxonomic_profile_modification0923.RData"))



head(input_normalized_taxon)

input_normalized_taxon_asv <- input_normalized_taxon %>%
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Species) %>%
  pivot_longer(!c(OTU, Genus))

input_normalized_taxon_asv$Genus[input_normalized_taxon_asv$Genus == ""] <- "unassigned"
input_normalized_taxon_asv$Genus[input_normalized_taxon_asv$Genus == "g__"] <- "unassigned"

  
  
input_normalized_taxon_asv_unclassfied <- input_normalized_taxon_asv %>%
  filter(Genus == "unassigned") %>%
  group_by(OTU) %>%
  summarise(value = sum(value))

# percentage 
input_normalized_taxon_asv_unclassfied_percent <- input_normalized_taxon_asv_unclassfied %>%
  add_column(percent = .$value/sum(.$value)) %>%
  arrange(desc(percent))

top20 <- head(input_normalized_taxon_asv_unclassfied_percent,n = 20)

filter(input_normalized_taxon, OTU == "34835afbd7b2e3ce535d53436d8b9a86") %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species)
#        Kingdom            Phylum                  Class                Order                Family Genus Species
1: k__Bacteria p__Proteobacteria c__Gammaproteobacteria o__Enterobacteriales f__Enterobacteriaceae              
(base) lu@sbi-A1754:/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile/unclassfiied/72368f24-0da1-47e8-8947-e7cf98943164/data$ grep -A 1 '34835afbd7b2e3ce535d53436d8b9a86' dna-sequences.fasta 
>34835afbd7b2e3ce535d53436d8b9a86
TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCA

filter(input_normalized_taxon, OTU == "a199468b017e0ae5750066cbb047560d") %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species)
> filter(input_normalized_taxon, OTU == "a199468b017e0ae5750066cbb047560d") %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species)
Kingdom            Phylum                  Class          Order               Family Genus Species
1: k__Bacteria p__Proteobacteria c__Alphaproteobacteria o__Rhizobiales f__Hyphomicrobiaceae   g__     s__
grep -A 1 'a199468b017e0ae5750066cbb047560d' dna-sequences.fasta 
>a199468b017e0ae5750066cbb047560d
TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTGTGATGAAGGCCTTAGGGTTGTAAAACACTTTCATCGGTGAAGATAATGACGGTAGCCGAAGAAGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACTGATCAGTTGGGGGTGAAATCCCGAGGCTTAACCTCGGAACTGCCTTCAATACTGTCAGTCTTGAGATCGGAAGAGGTGAGTGGAATTCCTAGTGTAGAGGTGAAATTCGTAGATATTAGGAAGAACACCAGTGGCGAAGGCGGCTCACTGGTCCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCGAACA
filter(input_normalized_taxon, OTU == "62a553bf8e561b989f56cb38b943b9ba") %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species)
Kingdom            Phylum                  Class          Order Family Genus Species
1: k__Bacteria p__Proteobacteria c__Alphaproteobacteria o__Rhizobiales    f__   g__     s__

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile/"
save.image(paste0(output_dir, "11_taxonomic_profile_modification0928_checktaxon.RData"))


########################################################################################################


#        figure modification 10.16            # 

########################################################################################################

# modification 10.16 (change the yaxis to 0-100 not demical number)
output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile/"
load(paste0(output_dir, "11_taxonomic_profile_modification0923.RData"))

i <- "Phylum"
barchart[[i]] <- barchart[[i]] + scale_y_continuous(labels = scales::percent)
i <- "Family"
barchart[[i]] <- barchart[[i]] + scale_y_continuous(labels = scales::percent)
i <- "Genus"
barchart[[i]] <- barchart[[i]] + scale_y_continuous(labels = scales::percent)


for (i in c("Phylum", "Family", "Genus")){
  
  
  ggsave(filename = paste0(output_dir, i, "_barchart.pdf", sep = ""), barchart[[i]], width = 8)
  
}
save.image(paste0(output_dir, "11_taxonomic_profile_modification1016.RData"))


