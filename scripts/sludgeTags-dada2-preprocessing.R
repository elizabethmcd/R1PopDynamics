library(dada2)
library(phyloseq)
library(ggplot2)
library(ampvis2)
library(grid)
library(gridExtra)

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
samples.out <- as.data.frame(rownames(seqtab.nochim))
colnames(samples.out) <- c("timepoint")
rownames(samples.out) <- rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samples.out), tax_table(taxa))
write.csv(samples.out, "~/Desktop/sludge_samples.csv", quote=FALSE)
metadata = read.csv("~/Desktop/metadata.csv", row.names = 'X')
ps2 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa))

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
family
genus <- plot_bar(ps.top20, x="Sample", fill="Genus")

# save as a dataframe for now
OTU1 <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
sludgeOTUs <- as.data.frame(OTU1)
write.csv(df, "data/otu_tables/2019-08-23-DADA2-asv-table.csv", quote=FALSE)

# try all
all <- names(sort(taxa_sums(ps), decreasing=TRUE))
ps.all <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
# all levels plots
allfamily <- plot_bar(ps.all, x="Sample", fill="Family") + theme(legend.position="none")


# ampvis
#source the phyloseq_to_ampvis2() function from the gist
#devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")
ampvis2_obj <- phyloseq_to_ampvis2(ps2)
amp_heatmap(ampvis2_obj, tax_aggregate = "Phylum", tax_show=4, plot_colorscale = "sqrt", plot_values = FALSE) + theme(axis.text.y = element_text(size=8), legend.position="right")

genus <- amp_heatmap(ampvis2_obj, tax_aggregate = "Genus", tax_show=7, plot_colorscale = "sqrt", plot_values = FALSE) + theme(axis.text.y = element_text(size=8), legend.position="right")
cleaned.genus <- genus + theme(legend.position="none", axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

accumulibacter <- amp_heatmap(ampvis2_obj, tax_aggregate="OTU", tax_show=2, plot_colorscale="sqrt", plot_values = FALSE) + theme(axis.text.y=element_blank(), legend.position="right")
cleaned.accum <- accumulibacter + theme(legend.position="none", axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# qPCR data
ppk = read.table("data/raw_data/Accum-qPCR-Timeseries.csv", sep=',', header=TRUE)
ppk.df = as.data.frame(ppk)
ppk.df$Date <- as.Date(df$Date, "%m/%d/%y")
 # subset of dates
t1 <- as.Date("2009-01-02")
t2 <- as.Date("2011-12-02")
sub <- ppk.df %>% filter(Date>=t1 & Date<=t2)
sub.m <- melt(sub,id.vars="Date",measure.vars=c("Clade.IA", "Clade.IIA"))
sub.m$Date <- as.character(sub.m$Date)
# full timeseries
full.m <- melt(ppk.df, vars="Date", measure.vars=c("Clade.IA", "Clade.IIA"))
full.m$Date <- as.character(full.m$Date)
# subset plot
p1 <- ggplot(sub.m, aes(x=Date, y=variable, fill=value)) + geom_tile(color="white") + scale_x_discrete(labels=sub.m$Date, position="bottom", expand=c(0,0)) + scale_fill_gradientn(colors = rev(viridis_pal()(9)), limits=c(0,16000000)) + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(axis.text.x= element_text(angle=75, hjust=1))
# clean plot
p2 <- p1 + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank()) + theme(legend.position="none") + labs(x=NULL, y=NULL) + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
# full plot
f1 <- ggplot(full.m, aes(x=Date, y=variable, fill=value)) + geom_tile(color="white") + scale_x_discrete(labels=full.m$Date, position="bottom", expand=c(0,0)) + scale_fill_gradientn(colors = rev(viridis_pal()(9)), limits=c(0,16000000)) + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(axis.text.x= element_text(angle=75, hjust=1))
f2 <- f1 + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank()) + theme(legend.position="none") + labs(x=NULL, y=NULL) + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))

# stitch 16S and Accumulibacer ASV together
genus.tmp = ggplot_build(cleaned.genus)
ppk.tmp = ggplot_build(cleaned.accum)
g1 = ggplot_gtable(ppk.tmp) ; g2 = ggplot_gtable(genus.tmp)
n1 = length(ppk.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
n2 = length(genus.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
g = rbind(g1, g2, size="first")
ggsave(file="~/Desktop/test.png", plot=g, width=100, height=30, units=c("cm"))
    # can do this when I finish qPCR of the extra 16S samples to put directly on top with ppk1
    # ppk1 data instead of the ASV of Accumulibacter IA and IIA abundant ones
p1

# save qPCR plots
ggsave(file="~/Desktop/ppk1.png", plot=p1, width=30, height=10, units=c('cm'))
ggsave(file="~/Desktop/full-ppk1.png", plot=f2, width=60, height=5, units=c('cm'))
ggsave(file="~/Desktop/full-ppk1-legend.png", plot=f1, width=60, height=5, units=c('cm'))
