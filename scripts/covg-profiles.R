library(tidyverse)
library(viridis)

# Coverage and relative abundance calculations and plots

# Coverage stats for 10 metagenomes sequenced from specific time points for SNV profiling of these genomes that meet cutoffs

files_path <- "results/covg_breadth"
samples <- read.table(file.path(files_path, "samples.txt"), header=FALSE)
files <- file.path(files_path, samples$V1, "genomeCoverage.csv")

flanking_stats <- data_frame(filename = files) %>% 
  mutate(file_contents = map(filename, ~ read.csv(file.path(.)))) %>% 
  unnest() 

covg_stats <- flanking_stats %>% 
  select(filename, genome, coverage) %>% 
  mutate(sample = gsub("results/covg_breadth/", "", filename)) %>% 
  select(sample, genome, coverage) %>% 
  mutate(sample = gsub("-historicalRefs-quick-profile/genomeCoverage.csv", "", sample)) %>% 
  pivot_wider(names_from = sample, values_from = coverage)

covg_filtered <- flanking_stats %>% 
  mutate(sample = gsub("results/covg_breadth/", "", filename)) %>% 
  select(sample, genome, coverage, breadth) %>% 
  mutate(sample = gsub("-historicalVirusRefs-quick-profile/genomeCoverage.csv", "", sample)) %>% 
  filter(coverage > 10)

covg_filtered %>% group_by(sample) %>% count()

covg_filtered %>% 
  ggplot(aes(x=sample, y=coverage)) + geom_boxplot() + geom_jitter()

# write out list of good genomes that passed coverage filtes to pass for queue with inStrain to only process those genome-sample pairs 
write_tsv(covg_filtered, "results/covg_breadth/good-genomes-covg-list.txt", col_names=FALSE)

# Relative abundance of MAGs in the 10 samples 

relative_abundance <- read_tsv("results/covg_breadth/coverM-final-rel-abund.txt")
flanking_relative_abundance <- relative_abundance[,c(1, 7:16)]
colnames(flanking_relative_abundance) <- c("Genome", "2010-02-08", "2010-08-31", "2010-12-31", "2011-04-26", "2011-08-01", "2011-09-30", "2012-02-03", "2012-04-16", "2012-09-23", "2013-04-05")
flanking_relative_abundance$`2010-12-31` <- as.numeric(flanking_relative_abundance$`2010-12-31`)
flanking_relative_abundance[is.na(flanking_relative_abundance)] <- 0

flanking_table <- flanking_relative_abundance %>% 
  filter(Genome != c("unmapped")) %>% 
  filter(Genome != c("TKFM")) %>% 
  pivot_longer(!Genome, names_to = "sample", values_to = "relative_abundance")

metadata <- read.csv("results/binning/bin-codes.csv")
colnames(metadata) <- c("Genome", "Code")
rel_abund_metadata <- left_join(flanking_table, metadata)

# combine filtered genomes with sufficient coverage with genome metadata
covg_filtered$genome <- gsub(".fa", "", covg_filtered$genome)
colnames(covg_filtered) <- c("Sample", "Genome", "Coverage", "Breadth")
covg_filtered_info <- left_join(covg_filtered, metadata)

metagenome_info <- read.csv("metadata/metagenome_information/R1_Flanking_Metageomes_Information.csv") %>% 
  select(sample, operation_day)

rel_abundance_info <- left_join(rel_abund_metadata, metagenome_info) %>% 
  select(operation_day, Genome, relative_abundance, Code)

rel_abundance_info %>% ggplot(aes(x=as_factor(operation_day), y=relative_abundance, fill=Code)) + geom_bar(stat="identity", color="black", size=0.3, width=0.8) + scale_fill_brewer(palette="Set3", labels=c("Actinobacteria", "Alphaproteobacteria", "Bacteroidetes", "CAPIA", "CAPIIA", "Gammaproteobacteria", "Other Lineages", "EPV1 Phage")) + theme_classic() + scale_y_continuous(expand=c(0,0), limits=c(0,100), breaks=seq(0,100,25)) + ylab("% Relative Abundance") + xlab("Operation Day") + labs(fill="Genome Lineage") + theme(axis.text.y=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.title.y=element_text(size=12, face="bold"), axis.title.x=element_text(size=12, face="bold"), legend.title=element_text(size=10, face="bold"))

rel_abund_metadata %>% ggplot(aes(x=sample, fill=Code)) + geom_area(stat="count")
