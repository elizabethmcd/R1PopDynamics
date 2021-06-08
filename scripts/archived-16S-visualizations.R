####################
# Old visualization script for 16S data with phyloseq and ampvis objects
####################

# store ASV name
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
theme_set(theme_bw())
plot_richness(ps, x="timepoint", measures=c("Shannon"))

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
top10 <- names(sort(taxa_sums(ps), decreasing=TRUE)[1:10])
ps.top10 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top10 <- prune_taxa(top10, ps.top10)
plot_bar(ps.top10, x="timepoint", fill="Genus")
####################################
# basic statistics
####################################
# average Bray-Curtis dissimilarity for each sample
pairBray <- distance(ps, "bray")
braycalc <- as.matrix(pairBray)
avg <- colMeans(braycalc)
result <- rbind(braycalc, avg)
# average bray curtis dissim over time
avg_bc <- rownames_to_column(as.data.frame(avg))
colnames(avg_bc) <- c("Sample", "BC")
avg_bc %>% ggplot(aes(x=Sample, y=BC, group=1)) + geom_point() + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle=85, hjust=1))
braydf <- rownames_to_column(as.data.frame(braycalc))
colnames(braydf) <- c("")
# select first subset of samples before reseeding event
#
ps_pruned <- prune_samples(samples, ps)
prunedBray <- distance(ps_pruned, "bray")
prunedcalc <- as.matrix(prunedBray)
prunedavg <- colMeans(prunedcalc)
prunedresult <- rbind(prunedcalc, prunedavg)
avg_pruned_bc <- rownames_to_column(as.data.frame(prunedavg))
colnames(avg_pruned_bc) <- c("sample", "BC")
avg_pruned_bc %>% ggplot(aes(x=sample, y=BC, group=1)) + geom_point() + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle=85, hjust=1))

####################################
# visualization with bar charts and heatmaps
####################################

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
write.csv(df, "results/otu_tables/2021-06-08-DADA2-flanking-OTU-table.csv", quote=FALSE)

# try all
all <- names(sort(taxa_sums(ps), decreasing=TRUE))
ps.all <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
# all levels plots
allfamily <- plot_bar(ps.all, x="Sample", fill="Family")
allfamily

# ampvis
#source the phyloseq_to_ampvis2() function from the gist
#devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")
ampvis2_obj <- phyloseq_to_ampvis2(ps2)
amp_heatmap(ampvis2_obj, tax_aggregate = "Phylum", tax_show=4, plot_colorscale = "sqrt", plot_values = FALSE) + theme(axis.text.y = element_text(size=8), legend.position="right")

genus <- amp_heatmap(ampvis2_obj, tax_aggregate = "Genus", tax_show=7, plot_colorscale = "sqrt", plot_values = FALSE) + theme(axis.text.y = element_text(size=8), legend.position="right")
genus
cleaned.genus <- genus + theme(legend.position="none", axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
cleaned.genus

amp_heatmap(ampvis2_obj, tax_aggregate="Genus", tax_show=10, plot_colorscale="sqrt", plot_values=FALSE) + theme(axis.text.y=element_blank())


accumulibacter <- amp_heatmap(ampvis2_obj, tax_aggregate="OTU", tax_show=2, plot_colorscale="sqrt", plot_values = FALSE) + theme(axis.text.y=element_blank(), legend.position="right")
cleaned.accum <- accumulibacter + theme(legend.position="none", axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
accumulibacter
cleaned.accum

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
