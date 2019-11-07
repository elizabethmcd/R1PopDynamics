library(phyloseq)
library(vegan)
library(tidyverse)

# otu table
sludge <- read.csv("data/otu_tables/2019-08-23-DADA2-asv-table.csv")
# after running cohesion, comparing turnover rates with connectivity and sample cohesion over time

# average Bray-Curtis dissimilarity for each sample
pairBray <- distance(ps, "bray")
braycalc <- as.matrix(pairBray)
avg <- colMeans(braycalc)
result <- rbind(braycalc, avg)
# looking at avg BC and cohesion values
pos <- as.matrix(cohesion.pos)
calc <- as.matrix(avg)
braypos <- cbind(pos, calc)
colnames(braypos) <- c("cohesionpos", "BC")
poscorr <- as.data.frame(braypos)
colnames(poscorr) <- c("cohesionpos", "BC")
poscorr %>% ggplot(aes(x=cohesionpos, y=BC)) + geom_point()
neg <- as.matrix(cohesion.neg)
brayneg <- cbind(neg, calc)
negcorr <- as.data.frame(brayneg)
colnames(negcorr) <- c("cohesionneg", "BC")
negcorr %>% ggplot(aes(x=cohesionneg, y=BC)) + geom_point()
# cohesion over time
cohesionPos <- as.data.frame(pos)
cohesionPos <- rownames_to_column(cohesionPos)
colnames(cohesionPos) <- c("Sample", "cohesion")
cohesionPos %>% ggplot(aes(x=Sample, y=cohesion, group=1)) + geom_point() + geom_line()
cohesionNeg <- as.data.frame(neg)
cohesionNeg <- rownames_to_column(cohesionNeg)
colnames(cohesionNeg) <- c("Sample", "cohesion")
cohesionNeg %>% ggplot(aes(x=Sample, y=cohesion, group=1)) + geom_point() + geom_line()
# average bray curtis dissim over time
avg_bc <- rownames_to_column(as.data.frame(avg))
colnames(avg_bc) <- c("Sample", "BC")
avg_bc %>% ggplot(aes(x=Sample, y=BC, group=1)) + geom_point() + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle=85, hjust=1))
