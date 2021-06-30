library(tidyverse)

# merge checkm, ANI, and classification for clusters that are very similar

# Checkm
checkm_stats <- read.delim("results/binning/dereplicated_checkm_stats.tsv", header=FALSE, sep="\t")
colnames(checkm_stats) <- c("bin1", "lineage", "completeness", "contamination", "size", "contigs", "gc")
# ANI
similar_clusters <- read.delim("results/binning/similar-clusters-ANI.tsv", header=FALSE, sep="\t") %>% 
  select(V1, V2, V3)
colnames(similar_clusters) <- c("bin1", "bin2", "ANI")
similar_clusters$bin1 <- gsub("/home/GLBRCORG/emcdaniel/EBPR/Flanking/metagenomes/binningResults/checkm_stats/dRep/dereplicated_genomes/", "", similar_clusters$bin1)
similar_clusters$bin2 <- gsub("/home/GLBRCORG/emcdaniel/EBPR/Flanking/metagenomes/binningResults/checkm_stats/dRep/dereplicated_genomes/", "", similar_clusters$bin2)
similar_clusters$bin1 <- gsub(".fa", "", similar_clusters$bin1)
similar_clusters$bin2 <- gsub(".fa", "", similar_clusters$bin2)
# Classifications 
gtdb <- read.delim("results/binning/ebpr-flanking-mags-gtdb.tsv", sep="\t")
colnames(gtdb) <- c("bin1", "classification")

all_bin_stats <- left_join(gtdb, checkm_stats) %>% select(-lineage)

similar_cluster_info <- left_join(similar_clusters, all_bin_stats)
colnames(similar_cluster_info) <- c("bin1", "bin2", "ANI", "classification", "comp1", "cont1", "size1", "contigs1", "gc1")
checkm_stats_2 <- checkm_stats
colnames(checkm_stats_2)[1] <- c("bin2")
similar_clusters_info <- left_join(similar_cluster_info, checkm_stats_2) %>% 
  select(-lineage)
colnames(similar_clusters_info) <- c("bin1", "bin2", "ANI", "classification", "comp1", "cont1", "size1", "contigs1", "gc1", "comp2", "cont2", "size2", "contigs2", "gc2")

write.csv(similar_clusters_info, "results/binning/similar-clusters-info.csv", row.names = FALSE, quote = FALSE)
write.csv(all_bin_stats, "results/binning/all-bins-info.csv", row.names = FALSE, quote = FALSE)
