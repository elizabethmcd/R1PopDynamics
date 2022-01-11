library(tidyverse)
library(viridis)
library(RColorBrewer)
library(gridGraphics)
library(cowplot)

# Coverage, relative abundance calculations, and within sample diversity and plots

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

full_metadata <- read.csv("results/R1-Flanking-final-bins-info.csv") %>% 
  select(bin, Code, classification)
colnames(full_metadata)[1] <- c("Genome")

relative_abundance_table_info <- left_join(full_metadata, flanking_relative_abundance) %>% 
  select(-Genome)
write.csv(relative_abundance_table_info, "results/binning/bins-relative-abundance-table.csv", quote=FALSE, row.names = FALSE)

manual_brewer_palette <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

flanking_abundance <- rel_abundance_info %>% ggplot(aes(x=as_factor(operation_day), y=relative_abundance, fill=Code)) + geom_bar(stat="identity", color="black", size=0.3, width=0.8) + scale_fill_manual(values=manual_brewer_palette, labels=c("Actinobacteria", "Alphaproteobacteria", "Bacteroidetes", "CAPIA", "CAPIIA", "Gammaproteobacteria", "Other Lineages", "EPV1 Phage")) + theme_classic() + scale_y_continuous(expand=c(0,0), limits=c(0,100), breaks=seq(0,100,25)) + ylab("% Relative Abundance") + xlab("Operation Day") + labs(fill="Genome Lineage") + theme(axis.text.y=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y=element_text(size=12, face="bold"), axis.title.x=element_blank(), legend.title=element_text(face="bold"))
flanking_abundance

rel_abund_metadata %>% ggplot(aes(x=sample, fill=Code)) + geom_area(stat="count")

final_rel_abund_codes <- read.csv("results/binning/bins-relative-abundance-table-modified.csv") %>% 
  select(-X, -GTDB.Classification) %>% 
  pivot_longer(!Genome, names_to="sample", values_to="relative_abundance")

final_rel_abund_codes %>% ggplot(aes(x=fct_rev(Genome), y=relative_abundance)) + geom_boxplot() + coord_flip() + scale_y_continuous(expand=c(0,0)) + theme_classic()

# Within-sample nucleotide diversity and r2 for genomes in samples that are above 10X coverage 

files_path <- "results/SNV_diversity/diversity"
files <- dir(files_path, pattern="*_genome_info.tsv")

diversity_files <- data_frame(filename = files) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(files_path, .)))) %>% 
  unnest() %>% 
  select(filename, genome, coverage, nucl_diversity, r2_mean, d_prime_mean)

diversity_table <- separate(diversity_files, filename, into=c("reference", "sample"), sep="-R1-") %>% 
  select(genome, sample, coverage, nucl_diversity, r2_mean, d_prime_mean) 

diversity_table$sample <- gsub("-inStrain.IS_genome_info.tsv", "", diversity_table$sample)
diversity_table$genome <- gsub(".fa", "", diversity_table$genome)

diversity_table_filtered <- diversity_table %>% 
  filter(coverage > 10)
colnames(diversity_table_filtered)[1] <- "Genome"

diversity_table_info <- left_join(diversity_table_filtered, metadata)
diversity_table_final <- left_join(diversity_table_info, metagenome_info)

flanking_nucleotide_diversity <- diversity_table_final %>% ggplot(aes(x=as_factor(operation_day), y=nucl_diversity, color=Code)) + geom_point(size=3, alpha=0.8) + scale_color_manual(values=manual_brewer_palette, labels=c("Actinobacteria", "Alphaproteobacteria", "Bacteroidetes", "CAPIA", "CAPIIA", "Gammaproteobacteria", "Other Lineages", "EPV1 Phage")) + scale_y_continuous(expand=c(0,0), limits=c(0,.015), breaks=seq(0,.015,.0025)) + ylab("Nucleotide Diveristy π") + xlab("Operation Day") + labs(color="Genome Lineage") + theme_classic()  + theme(axis.text.y=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y=element_text(size=12, face="bold"), axis.title.x=element_text(size=12, face="bold"), legend.title=element_text(face="bold"))
flanking_nucleotide_diversity

diversity_abundance <- left_join(diversity_table_final, rel_abundance_info, by=c("Genome", "operation_day")) %>% 
  select(Genome, sample, operation_day, relative_abundance, nucl_diversity, coverage, r2_mean, Code.x)

diversity_abund_plot <- diversity_abundance %>% ggplot(aes(x=nucl_diversity, y=relative_abundance, color=Code.x)) + geom_point() + scale_color_manual(values=manual_brewer_palette, labels=c("Actinobacteria", "Alphaproteobacteria", "Bacteroidetes", "CAPIA", "CAPIIA", "Gammaproteobacteria", "Other Lineages", "EPV1 Phage")) + scale_x_continuous(expand=c(0,0), limits=c(0,.015), breaks=seq(0,.015,.0025)) + labs(color="Genome Lineage") + ylab("% Relative Abundance") + xlab("Nucleotide Diversity π") + theme_classic() + theme(axis.title.y=element_text(face="bold"), axis.title.x=element_text(face="bold"), legend.title=element_text(face="bold"), legend.position=c(.9, .55))

diversity_faceted <- diversity_table_final %>% ggplot(aes(x=as_factor(operation_day), y=nucl_diversity, color=Genome)) + geom_point() + facet_wrap(~ Code, nrow=2) + theme_bw() + theme(legend.position = "none") + xlab("Operation Day") + ylab("Nucleotide Diversity π") + theme(axis.title.y=element_text(face="bold"), axis.title.x=element_text(face="bold"))

flanking_abund_div_grid <- plot_grid(flanking_abundance, flanking_nucleotide_diversity, labels=c("A", "B"), ncol=1)
flanking_abund_div_grid

flanking_abund_div_grid2 <- plot_grid(flanking_abundance, flanking_nucleotide_diversity, labels=c("C", "D"), ncol=1)

supp_div_abund <- plot_grid(diversity_abund_plot, diversity_faceted, ncol=1)
supp_div_abund

recombination_plot <- diversity_table_final %>% ggplot(aes(x=d_prime_mean, y=r2_mean)) + geom_point(aes(color=Code)) + geom_smooth(method="lm", se=TRUE, color="black") + ylab(expression(r ^2)) + xlab("D'") + labs(color="Genome Lineage") + scale_color_manual(values=manual_brewer_palette) + theme_bw() + theme(legend.position=c(.15, .7))
recombination_plot

ggsave("figs/flanking-abund-div-grid.png", flanking_abund_div_grid, width=10, height=6, units=c("in"))

ggsave("figs/flanking-abund-div-facets-supp.png", supp_div_abund, width=12, height=6, units=c("in"))

ggsave("figs/recombination-comparisons.png", recombination_plot, width=10, height=6, units=c("in"))

flanking_abund_div_grid


# full metadata of bins
# completeness, contamination, rRNAs, tRNAs 

complete_metadata <- read.csv("results/binning/R1-Flanking-final-bins-rRNA-info.csv")

genome_qual_plot <- complete_metadata %>% ggplot(aes(x=completeness, y=contamination, color=Lineage, size=rRNA)) + geom_point() + theme_classic()

complete_metadata %>% 
  filter(completeness > 90 & contamination < 10) %>% 
  count() # 136

complete_metadata %>% 
  filter(rRNA > 2) %>% 
  filter(completeness > 90 & contamination < 10) %>% 
  filter(X5S > 0 & X16S > 0 & X23S > 0) %>% 
  count() # 47

complete_metadata %>% 
  filter(X16S == 1) %>% 
  filter(completeness > 90 & contamination < 10) %>% 
  filter(Notes != "partial 16S") %>% 
  filter(Notes != "partial 16S, 23S") %>% 
  count() # 55

ggsave("figs/bins-quality-plot.png", genome_qual_plot, width=8, height=4, units=c("in"))

# Plot grid of chemical dynamics/clade dynamics/and the metagenome dynamics all in one figure

dynamics_grid <- plot_grid(p_data, clade_data, labels=c("A", "B"), label_x = 0, label_y=1, ncol=1, align="v")
dynamics_grid
flanking_abund_div_grid2 <- plot_grid(flanking_abundance, flanking_nucleotide_diversity, labels=c("C", "D"), label_x=0, label_y=1, ncol=1, align="v")
flanking_abund_div_grid2

flanking_complete_grid <- plot_grid(dynamics_grid, flanking_abund_div_grid2, ncol=2, rel_widths=c(2.5,1))
flanking_complete_grid
ggsave("figs/Flanking-chemical-diversity-complete-grid.png", flanking_complete_grid, width=24, height=12, units=c("in"))


new_flanking_grid <- plot_grid(flanking_abundance, flanking_nucleotide_diversity, ncol=1, align="v", labels=c("B", "C"))
new_flanking_grid
flanking_grid_v2 <- plot_grid(p_data, new_flanking_grid, labels=c("A"), ncol=1, rel_heights = c(1.3,2))

ggsave("~/Desktop/flanking-figure.png", flanking_grid_v2, width=17, height=11, units=c("in"))
