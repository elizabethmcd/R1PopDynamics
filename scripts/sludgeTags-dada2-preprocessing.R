library(dada2)
library(phyloseq)
library(ggplot2)

# before starting the dada2 workflow, samples must be demultiplexed, primers/adapters are removed, and the F and R files contain reads in matching order
# this preprocessing script works through a time-series of samples from engineered bioreactors, amplified the 16S region using the V3-V4 primers/region, and was sequenced with the Illumina 2x300 chemistry (incorrectly, was supposed to be 2x250, but will be more to throw out)
# forward primer: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG
# reverse primer: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC

# primer lengths
primerF <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
primerR <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
lengthF <- nchar(primerF)
lengthR <- nchar(primerR)

# path and files setup
path <- "SludgeTags/raw"
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# inspect quality profiles of a couple forward/reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# filter and quality trim
# the forward reads look really good up until the end
# the reverse reads really crap out around 210, which is expected for reverse Illumina reads and with the 2x250/300 chemistry
# sequencing done on V3-V4 region with 2x300 chemistry instead of 2x250
# cut forward to 275, cut reverse to 200, need to check overlap length

filtFs <- file.path("SludgeTags/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("SludgeTags/filtered", paste0(sample.names,"_R_filt.fastq.gz"))
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
# if a lot of reads are removed as chimeric, have to check the removal of primers, ambiguous nucleotides in unremoved primers can interfere with ID'ing chimeras
# when account for quality and trimming off left sides for primers, gets the most sequences merged with a good matching overlap and least amount of chimeric sequences

# taxa assigned with silva, general and exact species
taxa <- assignTaxonomy(seqtab.nochim, "databases/silva_nr_v132_train_set.fa.gz", multithread=FALSE)
taxa <- addSpecies(taxa, "databases/silva_species_assignment_v132.fa.gz")
# taxa assigned with GTDB taxonomy, general and exact species
# gtdb_taxa <- assignTaxonomy(seqtab.nochim, "databases/GTDB_bac-arc_ssu_r86.fa.gz", multithread=FALSE)
# gtdb_taxa <- addSpecies(gtdb_taxa, "databases/GTDB_dada2_assignment_species.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL

# create a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))

# store ASV name
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# some terrible bar charts
# top 20
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
family <- plot_bar(ps.top20, x="Sample", fill="Family")
genus <- plot_bar(ps.top20, x="Sample", fill="Genus")
# top 50
top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
family50 <- plot_bar(ps.top50, x="Sample", fill="Family")
genus50 <- plot_bar(ps.top50, x="Sample", fill="Genus")
# top 100 for fun
top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
family100 <- plot_bar(ps.top100, x="Sample", fill="Family")
genus100 <- plot_bar(ps.top100, x="Sample", fill="Genus")
phyla100 <- plot_bar(ps.top100, x="Sample", fill="Phylum")

# save as a dataframe for now
OTU1 <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
df <- as.data.frame(OTU1)
write.csv(df, "data/otu_tables/2019-08-23-DADA2-asv-table.csv", quote=FALSE)
family50
genus50
