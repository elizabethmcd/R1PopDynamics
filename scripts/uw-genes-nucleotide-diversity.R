library(tidyverse)
library(cowplot)
library(gridGraphics)
library(RColorBrewer)

#################################
# UW1 and UW3 gene-wise nucleotide diversity
#################################

uw_genes_path <- "results/SNV_diversity/UW_SNVs"
uw_genes_files <- dir(uw_genes_path, pattern="*_gene_info.tsv")

uw_genes <- data_frame(filename = uw_genes_files) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(uw_genes_path, .)))) %>% 
  unnest() %>% 
  select(filename, gene, coverage, SNV_count,  nucl_diversity)

uw_genes_table <- separate(uw_genes, filename, into=c("reference", "sample"), sep="-R1-") %>% 
  select(reference, sample, gene, coverage, SNV_count, nucl_diversity) %>% 
  filter(coverage > 5) # only keep genes with a coverage of 10X, gets rid of ~5000 genes for both refs across all samples

uw_genes_table$sample <- gsub("-inStrain.IS_gene_info.tsv", "", uw_genes_table$sample)

uw_genes_table %>% 
  group_by(reference, sample) %>% 
  count() # group by ref then sample to count how many genes for each ref in each sample meet the threshold > 5X coverage

aggregate(SNV_count ~ sample, uw_genes_table, sum)
uw_genes_table %>% 
  group_by(reference, sample) %>% 
  summarise(sum_SNV = sum(SNV_count))
           
uw_nucl_div_table %>% ggplot(aes(x=sample, y=nucl_diversity)) + geom_jitter(alpha=0.2) + facet_wrap(~ reference) # nucl div of all genes not separated by category

uw_genes_table$gene <- gsub("gnl\\|X\\|", "", uw_genes_table$gene)
colnames(uw_genes_table) <- c("reference", "sample", "locus_tag", "coverage", "SNV_count", "nucl_diversity")

#################################
# COGs CSV from AcDiv project for defining groups of genes that are core, core among Accumulibacter, and accessory reciprocally in UW1 and UW3
#################################

cog_table <- read.csv("results/pangenomics/UW1_UW3_COGS_TAGS_TABLE.csv")

uw_cog_diversity_table <- left_join(uw_genes_table, cog_table, by="locus_tag") %>% 
  select(reference.x, sample, locus_tag, coverage, SNV_count, nucl_diversity, label) %>% 
  drop_na()

colnames(uw_cog_diversity_table)[1] <- c("reference")

uw_div_genes_table <- left_join(uw_cog_diversity_table, metagenome_info) %>% 
  select(reference, operation_day, locus_tag, coverage, SNV_count, nucl_diversity, label)

uw_div_genes_table$label <- factor(uw_div_genes_table$label, levels=c("core", "acc_core", "accessory"))

uw.labs <- c("UW1 IIA", "UW3 IA")
names(uw.labs) <- c("UW1", "UW3")
genes.labs <- c("Core", "Accumulibacter Core", "Accessory")
names(genes.labs) <- c("core", "acc_core", "accessory")

uw_div_genes_plot <- uw_div_genes_table %>% ggplot(aes(x=factor(operation_day), y=nucl_diversity)) + geom_jitter(alpha=0.2) + facet_grid(rev(vars(label)), switch='y', vars(reference), labeller=labeller(label=genes.labs, reference=uw.labs)) + scale_y_continuous(breaks=seq(0, 0.05, 0.005), position="right") + ylab("Nucleotide Diversity π \n") + xlab("\n Operation Day") + theme_bw() + theme(strip.placement = "outside")
uw_div_genes_plot

ggsave("figs/UW1-UW3-genes-diversity-time-series.png", uw_div_genes_plot, width=9, height=7, units=c("in"))

#################################
# CRISPR Diveristy in UW3 
#################################

uw3_crispr <- uw_genes_table %>% 
  filter(reference == 'UW3') %>% 
  filter(locus_tag == 'PEBOPBNB_2_544' | locus_tag == 'PEBOPBNB_2_545' | locus_tag == 'PEBOPBNB_2_546' | locus_tag == 'PEBOPBNB_2_547' | locus_tag == 'PEBOPBNB_2_548' | locus_tag == 'PEBOPBNB_2_549' | locus_tag == 'PEBOPBNB_2_551')

uw3_crispr %>% ggplot(aes(x=locus_tag, y=coverage)) + geom_point() + facet_wrap(~ sample, scales = "free", ncol=2)

uw3_crispr %>% ggplot(aes(x=locus_tag, y=nucl_diversity)) + geom_point() + facet_wrap(~ sample, ncol=2)
