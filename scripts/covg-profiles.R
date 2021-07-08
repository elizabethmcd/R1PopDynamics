library(tidyverse)

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
  mutate(sample = gsub("-historicalRefs-quick-profile/genomeCoverage.csv", "", sample)) %>% 
  filter(coverage > 10)

covg_filtered %>% group_by(sample) %>% count()

covg_filtered %>% 
  ggplot(aes(x=sample, y=coverage)) + geom_boxplot() + geom_jitter()
