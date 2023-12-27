# do the taxonomic distribution based on silva database classifier prediction 
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

setwd("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile_silva_classifier/")

#################################################################################################


#				input file 							# 


#################################################################################################


load(paste0("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/9_diversity/",
            "9_alpha_diversity_analysis_HE_series_0909.RData"))

write_tsv(alpha_table_add_sel, "alpha_diversity_metadata_average.tsv")
write_tsv(metadata_combine_all, "metadata_combine_all.tsv")
rm(list = ls())

metadata_combine_all <- read_tsv("metadata_combine_all.tsv")
taxonomy_table_HE_series <- read_tsv("/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile_silva_classifier/table-with-taxonomy-silva-classifer.tsv", skip = 1)


#################################################################################################


#			separately the taxonomy 						#


#################################################################################################

otu_table <- taxonomy_table_HE_series #173 samples 
metatable <- metadata_combine_all # 139 samples 
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

output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile_silva_classifier//"

write_tsv(taxonomy_table_HE_series_tax, paste0(output_dir, "HE_series_taxon.tsv"))#9702 asvs 


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

write_tsv(HE_series_average, paste0(output_dir, "HE_series_average_count_table_silva.tsv"))



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

# May2020_Class8_SClass2_B2 84  453.6170
# Aug2020_Class8_SClass3_B2 128 661.0689
# removed : see 9_alpha_diversity_analysis_new 



#######################################################################################################


#                                   prepare input                                                    # 


#######################################################################################################

metadata_taxon <- read_tsv("alpha_diversity_metadata_average.tsv") # the updated samples 63 samples 
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
  taxon_table_long$taxon[taxon_table_long$taxon == "Ambiguous_taxa"] <- "unassigned"
  
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



# phylum manually color values=c("4" = "red", "r" = "blue", "f" = "darkgreen")


i <- "Phylum"
taxon_table_long_combined_mean[[i]][["combine_table_top10"]]$taxon %>% unique()

barchart[[i]] <- draw_barplot_phylum(taxon = taxon_table_long_combined_mean[[i]][["combine_table_top10"]], tag = i)
#greengenes result: 
#[1] "Others"             "p__Acidobacteria"   "p__Actinobacteria"  "p__Bacteroidetes"   "p__Chloroflexi"     "p__Cyanobacteria"   "p__Firmicutes"      "p__Planctomycetes" 
#[9] "p__Proteobacteria"  "p__Tenericutes"     "p__Verrucomicrobia"

#[1] "D_1__Acidobacteria"      "D_1__Actinobacteria"     "D_1__Bacteroidetes"      "D_1__Chloroflexi"        "D_1__Cyanobacteria"      "D_1__Epsilonbacteraeota" "D_1__Planctomycetes"    
#[8] "D_1__Proteobacteria"     "D_1__Tenericutes"        "D_1__Verrucomicrobia"    "Others"         

draw_barplot_phylum <- function(taxon = NA, tag = NA){
  
  barchart1 <- ggplot(taxon,aes(x = Group, y = mean_abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual(name = tag, values = c("Others" = "grey", 
                                             "D_1__Acidobacteria" = "forestgreen", 
                                             "D_1__Actinobacteria" = "red2",
                                             "D_1__Bacteroidetes" = "orange",
                                             "D_1__Chloroflexi" = "cornflowerblue",
                                             "D_1__Cyanobacteria" = "magenta",
                                             "D_1__Firmicutes" = "darkolivegreen4",
                                             "D_1__Planctomycetes" = "indianred1",
                                             "D_1__Proteobacteria" = "tan4",
                                             "D_1__Tenericutes" = "darkblue",
                                             "D_1__Verrucomicrobia" = "mediumorchid1",
                                             "D_1__Epsilonbacteraeota" = "lightsalmon"
                                             )) +
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
#greengenes classifer result 
#[1] "f__Alteromonadaceae"    "f__Flammeovirgaceae"    "f__Flavobacteriaceae"   "f__Hyphomicrobiaceae"   "f__Phyllobacteriaceae"  "f__Piscirickettsiaceae" "f__Pseudanabaenaceae"
#[8] "f__Rhodobacteraceae"    "f__Saprospiraceae"      "Others"                 "unassigned"          

# current result 
#[1] "D_4__Cyclobacteriaceae"    "D_4__Flavobacteriaceae"    "D_4__Phycisphaeraceae"     "D_4__Rhizobiaceae"         "D_4__Rhodobacteraceae"     "D_4__Saprospiraceae"      
#[7] "D_4__Sphingomonadaceae"    "D_4__uncultured"           "D_4__uncultured bacterium" "Others"                    "unassigned"    

draw_barplot_family <- function(taxon = NA, tag = NA){
  
  barchart1 <- ggplot(taxon,aes(x = Group, y = mean_abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual(name = tag, values = c("Others" = "grey", 
                                             "D_4__Alteromonadaceae" = "forestgreen", 
                                             "D_4__Flammeovirgaceae" = "red2",
                                             "D_4__Flavobacteriaceae" = "orange",
                                             "D_4__Hyphomicrobiaceae" = "cornflowerblue",
                                             "D_4__Phyllobacteriaceae" = "magenta",
                                             "D_4__Piscirickettsiaceae" = "darkolivegreen4",
                                             "D_4__Pseudanabaenaceae" = "indianred1",
                                             "D_4__Rhodobacteraceae" = "tan4",
                                             "D_4__Saprospiraceae" = "darkblue",
                                             "unassigned" = "mediumorchid1",
                                             "D_4__Hyphomonadaceae" = "firebrick4",
                                             "D_4__Oceanospirillaceae" = "yellowgreen",
                                             "D_4__Vibrionaceae" = "lightsalmon",
                                             "D_4__Cyclobacteriaceae" = "tan3",
                                             "D_4__Phycisphaeraceae" = "tan1",
                                             "D_4__Rhizobiaceae" = "wheat4",
                                             "D_4__Sphingomonadaceae" = "#DDAD4B",
                                             "D_4__uncultured" = "chartreuse",
                                             "D_4__uncultured bacterium" = "seagreen1")) +
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
#greengenes result 
#[1] "g__Aquimarina"      "g__Arcobacter"      "g__Loktanella"      "g__Maribacter"      "g__Mycoplasma"      "g__Phaeobacter"     "g__Tenacibaculum"   "g__Vibrio"          "g__Winogradskyella"
#[10] "Others"             "unassigned" 

#[1] "D_5__Aquimarina"           "D_5__Arcobacter"           "D_5__Maribacter"           "D_5__Sphingorhabdus"       "D_5__Sulfitobacter"        "D_5__uncultured"          
#[7] "D_5__uncultured bacterium" "D_5__uncultured organism"  "D_5__Vibrio"               "Others"                    "unassigned"   

draw_barplot_genus <- function(taxon = NA, tag = NA){
  
  barchart1 <- ggplot(taxon,aes(x = Group, y = mean_abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual(name = tag, values = c("Others" = "grey", 
                                             "D_5__Aquimarina" = "forestgreen", 
                                             "D_5__Arcobacter" = "red2",
                                             "D_5__Loktanella" = "orange",
                                             "D_5__Maribacter" = "cornflowerblue",
                                             "D_5__Mycoplasma" = "magenta",
                                             "D_5__Octadecabacter" = "darkolivegreen4",
                                             "D_5__Phaeobacter" = "indianred1",
                                             "D_5__Tenacibaculum" = "tan4",
                                             "D_5__Vibrio" = "darkblue",
                                             "unassigned" = "mediumorchid1",
                                             "D_5__Bdellovibrio" = "firebrick4",
                                             "D_5__Devosia" = "yellowgreen",
                                             "D_5__Glaciecola" = "lightsalmon",
                                             "D_5__Pseudoalteromonas" = "tan3",
                                             "D_5__Roseovarius" = "tan1",
                                             "D_5__Ruegeria" = "wheat4",
                                             "D_5__Sphingorhabdus" = "#DDAD4B",
                                             "D_5__Sulfitobacter" = "moccasin",
                                             "D_5__uncultured" = "chartreuse",
                                             "D_5__uncultured bacterium" = "seagreen1",
                                             "D_5__uncultured organism" = "mediumvioletred")) +
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


save.image(paste0(output_dir, "11_taxonomic_profile.RData"))


########################################################################################################


#        figure modification 10.16            # 

########################################################################################################

# modification 10.16 (change the yaxis to 0-100 not demical number)
output_dir <- "/media/lu/Seagate Por/documents/201911_hydractina/2022HE_aquarium_time_series/new_analysis_with_class2_class4/11_taxonomy_profile_silva_classifier/"
load(paste0(output_dir, "11_taxonomic_profile.RData"))

i <- "Phylum"
new_label <- taxon_table_long_combined_mean[[i]][["combine_table_top10"]]$taxon %>% unique() 
for (pos in which(str_detect(new_label, "uncultured"))){
  new_label[pos] <- paste(new_label[pos], "*", sep = "")
  
}

barchart[[i]] <- barchart[[i]] + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(name = i, values = c("Others" = "grey", 
                                           "D_1__Acidobacteria" = "forestgreen", 
                                           "D_1__Actinobacteria" = "red2",
                                           "D_1__Bacteroidetes" = "orange",
                                           "D_1__Chloroflexi" = "cornflowerblue",
                                           "D_1__Cyanobacteria" = "magenta",
                                           "D_1__Firmicutes" = "darkolivegreen4",
                                           "D_1__Planctomycetes" = "indianred1",
                                           "D_1__Proteobacteria" = "tan4",
                                           "D_1__Tenericutes" = "darkblue",
                                           "D_1__Verrucomicrobia" = "mediumorchid1",
                                           "D_1__Epsilonbacteraeota" = "lightsalmon"
  ), labels = new_label)


i <- "Family"

new_label_family <- taxon_table_long_combined_mean[[i]][["combine_table_top10"]]$taxon %>% unique() 
for (pos in which(str_detect(new_label_family, "uncultured"))){
  new_label_family[pos] <- paste(new_label_family[pos], "*", sep = "")
  
}

barchart[[i]] <- barchart[[i]] + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(name = i, values = c("Others" = "grey", 
                                           "D_4__Alteromonadaceae" = "forestgreen", 
                                           "D_4__Flammeovirgaceae" = "red2",
                                           "D_4__Flavobacteriaceae" = "orange",
                                           "D_4__Hyphomicrobiaceae" = "cornflowerblue",
                                           "D_4__Phyllobacteriaceae" = "magenta",
                                           "D_4__Piscirickettsiaceae" = "darkolivegreen4",
                                           "D_4__Pseudanabaenaceae" = "indianred1",
                                           "D_4__Rhodobacteraceae" = "tan4",
                                           "D_4__Saprospiraceae" = "darkblue",
                                           "unassigned" = "mediumorchid1",
                                           "D_4__Hyphomonadaceae" = "firebrick4",
                                           "D_4__Oceanospirillaceae" = "yellowgreen",
                                           "D_4__Vibrionaceae" = "lightsalmon",
                                           "D_4__Cyclobacteriaceae" = "tan3",
                                           "D_4__Phycisphaeraceae" = "tan1",
                                           "D_4__Rhizobiaceae" = "wheat4",
                                           "D_4__Sphingomonadaceae" = "#DDAD4B",
                                           "D_4__uncultured" = "chartreuse",
                                           "D_4__uncultured bacterium" = "seagreen1"), labels = new_label_family)






i <- "Genus"

new_label_genus <- taxon_table_long_combined_mean[[i]][["combine_table_top10"]]$taxon %>% unique() 
for (pos in which(str_detect(new_label_genus, "uncultured"))){
  new_label_genus[pos] <- paste(new_label_genus[pos], "*", sep = "")
  
}

barchart[[i]] <- barchart[[i]] + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(name = i, values = c("Others" = "grey", 
                                           "D_5__Aquimarina" = "forestgreen", 
                                           "D_5__Arcobacter" = "red2",
                                           "D_5__Loktanella" = "orange",
                                           "D_5__Maribacter" = "cornflowerblue",
                                           "D_5__Mycoplasma" = "magenta",
                                           "D_5__Octadecabacter" = "darkolivegreen4",
                                           "D_5__Phaeobacter" = "indianred1",
                                           "D_5__Tenacibaculum" = "tan4",
                                           "D_5__Vibrio" = "darkblue",
                                           "unassigned" = "mediumorchid1",
                                           "D_5__Bdellovibrio" = "firebrick4",
                                           "D_5__Devosia" = "yellowgreen",
                                           "D_5__Glaciecola" = "lightsalmon",
                                           "D_5__Pseudoalteromonas" = "tan3",
                                           "D_5__Roseovarius" = "tan1",
                                           "D_5__Ruegeria" = "wheat4",
                                           "D_5__Sphingorhabdus" = "#DDAD4B",
                                           "D_5__Sulfitobacter" = "moccasin",
                                           "D_5__uncultured" = "chartreuse",
                                           "D_5__uncultured bacterium" = "seagreen1",
                                           "D_5__uncultured organism" = "mediumvioletred"), labels = new_label_genus)





for (i in c("Phylum", "Family", "Genus")){
  
  
  ggsave(filename = paste0(output_dir, i, "_barchart.pdf", sep = ""), barchart[[i]], width = 8)
  
}
save.image(paste0(output_dir, "11_taxonomic_profile_modification1016.RData"))



