library(dada2)
library(phyloseq)
library(ggplot2)
library(ampvis2)
library(grid)
library(gridExtra)
library(vegan)
library(tidyverse)
library(patchwork)
library(cowplot)

# before starting the dada2 workflow, samples must be demultiplexed, primers/adapters are removed, and the F and R files contain reads in matching order
# this preprocessing script works through a time-series of samples from engineered bioreactors, amplified the 16S region using the V3-V4 primers/region, and was sequenced with the Illumina 2x300 chemistry (incorrectly, was supposed to be 2x250, but will be more to throw out)
# forward primer: CCTACGGGNGGCWGCAG
# reverse primer: GACTACHVGGGTATCTAATCC
# also check that all fastq files and databases aren't in the cloud, will need to pull them down

# primer lengths
primerF <- "CCTACGGGNGGCWGCAG"
primerR <- "GACTACHVGGGTATCTAATCC"
lengthF <- nchar(primerF)
lengthR <- nchar(primerR)

# path and files setup
path <- "raw_data/FlankingTags2010/raw_data"
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# inspect quality profiles of a couple forward/reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# filter and quality trim
# the forward reads look really good up until the end
# the reverse reads really crap out around 210, which is expected for reverse Illumina reads and with the 2x250/300 chemistry
# sequencing done on V3-V4 region with 2x300 chemistry instead of 2x250
# cut forward to 290, cut reverse to 215, need to check overlap length

filtFs <- file.path("raw_data/FlankingTags2010/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("raw_data/FlankingTags2010/filtered", paste0(sample.names,"_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,215), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE, trimLeft = c(lengthF, lengthR))

# check quality again
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

# error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# sample inference for obtaining sequence variants
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# can use the dada function to pool sequences to inform from multiple samples, but can address that at a later time

# merge paired reads to get the full denoised dataasets
# most reads should merge, if not getting a lot of merge paired-reads check trimming parameters
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))
# can get specific length of expected amplicon length
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
# have to check the size with V3/V4 for overlap amount, since V4 pretty much completely overlaps

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# looking at number of reads that went through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
# if a lot of reads are removed as chimeric, have to check the removal of primers, ambiguous nucleotides in unremoved primers can interfere with ID'ing chimeras
# when account for quality and trimming off left sides for primers, gets the most sequences merged with a good matching overlap and least amount of chimeric sequences

# taxa assigned with silva, general and exact species
taxa <- assignTaxonomy(seqtab.nochim, "databases/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
taxa <- addSpecies(taxa, "databases/silva_species_assignment_v138.1.fa.gz")
# taxa assigned with GTDB taxonomy, general and exact species
# gtdb_taxa <- assignTaxonomy(seqtab.nochim, "databases/GTDB_bac-arc_ssu_r86.fa.gz", multithread=FALSE)
# gtdb_taxa <- addSpecies(gtdb_taxa, "databases/GTDB_dada2_assignment_species.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
taxa.print

####################################
# create a phyloseq object
####################################
samples.out <- as.data.frame(rownames(seqtab.nochim))
colnames(samples.out) <- c("timepoint")
rownames(samples.out) <- rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samples.out), tax_table(taxa))

# merge with sample information
# P date, operation day, SRT
sample_info <- read.csv("metadata/reactor_metadata/r1_subset_chemical_data_2010_2013.csv") %>% 
  select(Date, P_release, P_uptake, percent_P_removal, operation_day, SRT)
colnames(sample_info)[1] <- c("timepoint")

samples.out$timepoint <- gsub("R1-", "", samples.out$timepoint)
sample_metadata <- left_join(samples.out, sample_info)
sample_metadata$timepoint <- Map(paste, 'R1-', sample_metadata$timepoint)
row.names(sample_metadata) <- sample_metadata$timepoint
row.names(sample_metadata) <- gsub(" ", "", row.names(sample_metadata))
sample_metadata$timepoint <- gsub(" ", "", sample_metadata$timepoint)

# phyloseq object with metadata
ps2 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(sample_metadata), tax_table(taxa))

# Alpha diversity / Shannon diversity richness within samples

shannon_timeseries <- plot_richness(ps2, x="operation_day", measure="Shannon") + scale_x_continuous(expand=c(0,0), limits=c(0,1200), breaks=seq(0,1175,30)) + xlab(label="Operation Day") + ylab('Shannon\n Alpha Diversity\n') + theme_bw() + theme(axis.title.x=element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=7), strip.background=element_blank(), strip.text.x=element_blank(), plot.title=element_text(size=12, face="bold"), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))
shannon_timeseries

ggsave("figs/R1-shannon-timeseries.png", shannon_timeseries, width=14, height=4, units=c("in"))

# subset taxa to only include Accumulibacter ASVs
accumulibacter_ps <- subset_taxa(ps2, Genus=="Candidatus_Accumulibacter")
plot_richness(accumulibacter_ps, x="operation_day", measures="Shannon") + scale_x_continuous(expand=c(0,0), limits=c(0,1200), breaks=seq(0,1175,30)) + xlab(label="Operation Day") + ylab('Shannon\n Alpha Diversity\n') + theme_bw() + theme(axis.title.x=element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=7), strip.background=element_blank(), strip.text.x=element_blank(), plot.title=element_text(size=12, face="bold"), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))

# ampvis
#source the phyloseq_to_ampvis2() function from the gist
#devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")
ampvis2_obj <- phyloseq_to_ampvis2(ps2)

r1_genus_plot <- amp_heatmap(ampvis2_obj, tax_aggregate = "Genus", group_by="operation_day", tax_add="Phylum", tax_show=8, plot_values = FALSE, plot_legendbreaks=c(0.1,1,10,60)) + scale_y_discrete(labels=c("Actinobacteria; Tetrasphaera", "Actinobacteria; Leucobacter", "Proteobacteria; Brevundimonas", "Proteobacteria; Diaphorobacter", "Proteobacteria; Pseudoxanthomonas", "Proteobacteria; Gemmobacter", "Bacteroidetes; Chryseobacterium", "Candidatus Accumulibacter")) + theme(axis.text.y = element_text(face="bold.italic", size=6), legend.position="right") + theme(axis.text.x=element_text(angle=80, size=6, vjust=1))
r1_genus_plot

acc_asv_heatmap <- amp_heatmap(ampvis2_obj, tax_aggregate="OTU", tax_show=2, plot_values=FALSE, group_by="operation_day", plot_legendbreaks=c(0.1,1,10,60)) + scale_y_discrete(labels=c("Ca. Accumulibacter ASV2", "Ca. Accumulibacter ASV1"), expand=c(0,0)) + theme(axis.text.x=element_text(angle=80, size=6, vjust=1), axis.text.y=element_text(size=8))

amp_boxplot(ampvis2_obj, tax_show=8)

ggsave(filename="figs/r1-ampvis2-genus-plot-heatmap.png", r1_genus_plot, width=16, height=3, units=c("in"))

ggsave(filename="figs/acc-asvs-heatmap.png", acc_asv_heatmap, width=9, height=3, units=c("in"))

p1 <- plot_grid(r1_genus_plot, shannon_timeseries, ncol=1, labels=c("A", "B"), label_size=10, vjust=1)

ggsave(filename="figs/16S-heatmap-shannon-grid.png", p1, height=4, width=10, units=c("in"))
