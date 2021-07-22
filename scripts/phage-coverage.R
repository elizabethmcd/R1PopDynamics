library(tidyverse)

phage_genes <- read_tsv("results/SNV_diversity/EPV1/EPV1v-R1-2011-09-30-inStrain.IS_gene_info.tsv") %>% 
  select(gene, coverage, SNV_count, nucl_diversity)
phage_genes$index <- seq.int(nrow(phage_genes))

phage_genes_low <- read_tsv("results/SNV_diversity/EPV1/EPV1v-R1-2011-04-26-inStrain.IS_gene_info.tsv") %>% 
  select(gene, coverage, SNV_count, nucl_diversity)
phage_genes_low$index <- seq.int(nrow(phage_genes_low))

write.csv(phage_genes, "results/SNV_diversity/EPV1/phage_genes_index_1.csv", row.names = FALSE, quote=FALSE)
write.csv(phage_genes_low, "results/SNV_diversity/EPV1/phage_genes_index_2.csv", row.names = FALSE, quote=FALSE)

all_phage <- read.csv("results/SNV_diversity/EPV1/epv1-dates-stats.csv")

all_phage %>% 
  ggplot(aes(x=index, y=coverage)) + geom_area() + geom_line() + facet_wrap(~ sample, ncol=1, scales = "free") + scale_x_continuous(limits=c(1,55), breaks=seq(1,55,1), expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_classic() + xlab("Gene Index") + ylab("Coverage in Sample") + theme(axis.text.x=element_text(size=6))


all_phage %>% 
  filter(coverage > avg_covg) %>% 
  ggplot(aes(x=index, y=nucl_diversity, color=sample)) + geom_point() + scale_x_continuous(limits=c(1,55), breaks=seq(1,55,1))

p1 <- phage_genes_low %>% 
  ggplot(aes(x=index, y=coverage)) + geom_area() + geom_line() + scale_x_continuous(limits=c(1,51), breaks=seq(1,51,1), expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_classic() 

phage_genes %>% 
  ggplot(aes(x=gene, y=coverage)) + geom_area() + geom_line() + scale_x_discrete(limits=c(1,54), breaks=seq(1,54,1), expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + theme_classic()
