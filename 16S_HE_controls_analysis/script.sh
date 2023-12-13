# fix the header 

# copy the data here in this directory and rename the header 
# make a new manifest 

head -n1 20201007-manifest-12345.txt > 20201011-manifest.txt

cat 20201007-manifest.txt|grep -v "sample-id" | while read line 
do
	#echo $line 
	file1=`echo $line | awk '{print $2}'`
	filename=`echo $line | awk '{print $1}'` 
	file2=`echo $line | awk '{print $3}'` 
	zcat $file1 | bioawk -c fastx '{split($name,a,"#");print "@"a[1]"/1 #"a[2]"\n"$seq"\n""+""\n"$qual}' - > data/$filename"_1.fastq"
	zcat $file2 | bioawk -c fastx '{split($name,a,"#");print "@"a[1]"/2 #"a[2]"\n"$seq"\n""+""\n"$qual}' - > data/$filename"_2.fastq"
	echo -e $filename"\t""/sbidata/lzhang/201911_hydractinia/202010_16S_analysis/qiime2_deblur/1_Hydractinia_microbiome/Hydractinia/data/"$filename"_1.fastq""\t""/sbidata/lzhang/201911_hydractinia/202010_16S_analysis/qiime2_deblur/1_Hydractinia_microbiome/Hydractinia/data/"$filename"_2.fastq" >> 20201011-manifest.txt
		
done  


# import data 

qiime tools import \
--input-path 20201011-manifest.txt \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-format PairedEndFastqManifestPhred64V2 \
--output-path hydractinia-paried-end-demux.qza

# visualization samples 
#### 
qiime demux summarize \
  --i-data ./hydractinia-paried-end-demux.qza \
  --o-visualization ./hydractinia-paried-end-demux.qzv 


# cut adapt : my ref : https://forum.qiime2.org/t/q2-cutadapt-of-primer-sequences-flag-verification/5805 
# https://docs.qiime2.org/2019.10/plugins/available/cutadapt/trim-paired/ I should make sure the primer is 3 prime end 
# see test.script.sh : I have already checked it. 

#Here is the primer information:
#V3-V4 16S    
#341F:ACTCCTACGGGAGGCAGCAG    
#806R:GGACTACHVGGGTWTCTAAT


qiime cutadapt trim-paired --p-cores 32 --p-front-f ACTCCTACGGGAGGCAGCAG --p-front-r GGACTACHVGGGTWTCTAAT \
 --i-demultiplexed-sequences ./hydractinia-paried-end-demux.qza \
 --o-trimmed-sequences ./output/cutadapt/hydractinia-paried-end-demux-primerTrim.qza --verbose --p-discard-untrimmed
# reason i used the --p-discard-untrimmed: https://forum.qiime2.org/t/joined-reads-before-deblur/15777/2 he recommends 

# view quality 

qiime demux summarize \
  --i-data ./output/cutadapt/hydractinia-paried-end-demux-primerTrim.qza \
  --o-visualization ./output/cutadapt/hydractinia-paried-end-demux-primerTrim.qzv 


# joined the reads 

mkdir -p ./output/joined/

qiime vsearch join-pairs \
  --i-demultiplexed-seqs ./output/cutadapt/hydractinia-paried-end-demux-primerTrim.qza \
  --o-joined-sequences ./output/joined/hydractinia-demux-joined.qza

qiime demux summarize \
  --i-data ./output/joined/hydractinia-demux-joined.qza \
  --o-visualization ./output/joined/hydractinia-demux-joined.qzv


# quality filter and deblur denoise 

qiime quality-filter q-score-joined \
 --i-demux ./output/joined/hydractinia-demux-joined.qza \
 --o-filtered-sequences ./output/quality-filter/hydractinia-demux-joined-filtered.qza \
 --o-filter-stats ./output/quality-filter/hydractinia-demux-joined-filtered-stats.qza

# deblur denoise  I have looked into my hydractinia-demux-joined.qzv,most of my reads are of good quality, but 500+ with low quality. So I trimmed with 404(2nd,25% length) 

qiime deblur denoise-16S \
  --i-demultiplexed-seqs ./output/quality-filter/hydractinia-demux-joined‐filtered.qza \
  --p-trim-length 404 \
  --p-sample-stats \
  --o-representative-sequences ./output/deblur/hydractinia-demux-joined‐filtered-deblur-rep-seqs.qza \
  --o-table ./output/deblur/hydractinia-demux-joined‐filtered-deblur-table.qza \
  --o-stats ./output/deblur/hydractinia-demux-joined‐filtered-deblur-stats.qza \
  --p-jobs-to-start 36


qiime feature-table tabulate-seqs \
 --i-data ./output/deblur/hydractinia-demux-joined‐filtered-deblur-rep-seqs.qza \
 --o-visualization ./output/deblur/hydractinia-demux-joined‐filtered-deblur-rep-seqs.qzv 

qiime feature-table summarize \
 --i-table ./output/deblur/hydractinia-demux-joined‐filtered-deblur-table.qza \
 --m-sample-metadata-file 20201007-metatable-simple-deblur.xlsx \
 --o-visualization ./output/deblur/hydractinia-demux-joined‐filtered-deblur-table.qzv

# analysis part 1 - phylogenetic tree 

mkdir -p ./output/phylogenetic/

wget \
  -O "sepp-refs-gg-13-8.qza" \
  "https://data.qiime2.org/2019.10/common/sepp-refs-gg-13-8.qza"

mv sepp-refs-gg-13-8.qza ./output/phylogenetic/

qiime fragment-insertion sepp \
  --i-representative-sequences ./output/deblur/hydractinia-demux-joined‐filtered-deblur-rep-seqs.qza \
  --i-reference-database ./output/phylogenetic/sepp-refs-gg-13-8.qza \
  --o-tree ./output/phylogenetic/tree.qza \
  --o-placements ./output/phylogenetic/tree_placements.qza \
  --p-threads 36

# analysis part2 - taxonomy 


mkdir -p ./output/taxonomy-train/
# taxonomy classfier training 

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /sbidata/lzhang/201911_hydractinia/202010_16S_analysis/qiime2_vsearch/1_Hydractinia_microbiome/Hydractinia/gg_13_8_otus/rep_set/99_otus.fasta \
  --output-path ./output/taxonomy-train/99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /sbidata/lzhang/201911_hydractinia/202010_16S_analysis/qiime2_vsearch/1_Hydractinia_microbiome/Hydractinia/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
  --output-path ./output/taxonomy-train/99-ref-taxonomy.qza

qiime feature-classifier extract-reads \
  --i-sequences ./output/taxonomy-train/99_otus.qza \
  --p-f-primer ACTCCTACGGGAGGCAGCAG \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --o-reads ./output/taxonomy-train/ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ./output/taxonomy-train/ref-seqs.qza \
  --i-reference-taxonomy ./output/taxonomy-train/99-ref-taxonomy.qza \
  --o-classifier ./output/taxonomy-train/gg13_8_99_classifier.qza

# taxonomy-unit prediction 
mkdir -p ./output/taxonomy/

qiime feature-classifier classify-sklearn \
  --i-classifier ./output/taxonomy-train/gg13_8_99_classifier.qza \
  --i-reads ./output/deblur/hydractinia-demux-joined‐filtered-deblur-rep-seqs.qza \
  --o-classification ./output/taxonomy/taxonomy.qza

qiime metadata tabulate \
  --m-input-file ./output/taxonomy/taxonomy.qza \
  --o-visualization ./output/taxonomy/taxonomy.qzv

qiime taxa barplot \
  --i-table ./output/deblur/hydractinia-demux-joined‐filtered-deblur-table.qza \
  --i-taxonomy ./output/taxonomy/taxonomy.qza \
  --m-metadata-file /sbidata/lzhang/201911_hydractinia/202010_16S_analysis/qiime2_vsearch/1_Hydractinia_microbiome/Hydractinia/20201007-metatable-simple-deblur.xlsx \
  --o-visualization ./output/taxonomy/taxa-bar-plots-gg13_5_99-classfier.qzv

# output - export 

mkdir -p ./output/taxonomy_export 
# Creating a BIOM table with taxonomy annotations
qiime tools export --input-path ./output/deblur/hydractinia-demux-joined‐filtered-deblur-table.qza --output-path ./output/taxonomy_export/exported
qiime tools export --input-path ./output/taxonomy/taxonomy.qza --output-path ./output/taxonomy_export/exported


# below has to be done manually 
#---------------------------------------------------
cd ./output/taxonomy_export
cp exported/taxonomy.tsv biom-taxonomy.tsv

#Change the first line of biom-taxonomy.tsv (i.e. the header) to this:
# remember they are tab-separated and starts with # 
#OTUID  taxonomy        confidence

biom add-metadata -i exported/feature-table.biom -o exported/table-with-taxonomy.biom --observation-metadata-fp biom-taxonomy.tsv --sc-separated taxonomy
biom convert -i exported/table-with-taxonomy.biom -o exported/table-with-taxonomy.tsv --to-tsv --table-type "OTU table" --header-key taxonomy
biom summarize-table -i exported/table-with-taxonomy.biom -o exported/otu_table_tax.sum 



