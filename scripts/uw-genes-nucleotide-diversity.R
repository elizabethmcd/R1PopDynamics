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
