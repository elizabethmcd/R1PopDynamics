library(tidyverse)
library(readxl)
library(reshape2)
library(viridis)
library(patchwork)
library(cowplot)

# R1R2 metadata

####################################################
# R1
####################################################

r1_meta <- read_excel("metadata/r1-2005-2013-metadata-summaries.xlsx")
r1_clean <- r1_meta %>% 
  select(Date, `P release`, `P uptake`, TSS, VSS) %>%
  filter(!is.na(`P uptake`)) %>% 
  mutate(VSS_TSS = (VSS / TSS) * 100)

colnames(r1_clean) <- c("Date", "P_release", "P_uptake", "TSS", "VSS", "VSS_TSS")

r1_clean$P_uptake <- ifelse(r1_clean$P_uptake < 0, 0, r1_clean$P_uptake)
r1_clean$P_release <- ifelse(r1_clean$P_release < 0, 0, r1_clean$P_release)

r1_table <- r1_clean %>% 
  mutate(percent_P_removal = ((P_release - P_uptake) / P_release) * 100) %>% 
  filter(!is.na(P_release)) %>% 
  filter(P_release != 0)

r1_table$percent_P_removal <- ifelse(r1_table$percent_P_removal < 0, 0, r1_table$percent_P_removal)
r1_table$Date <- as.Date(r1_table$Date, "YYYY-MM-DD")
r1_table$VSS_TSS <- ifelse(r1_table$VSS_TSS < 0, 0, r1_table$VSS_TSS)
r1_table$VSS_TSS <- ifelse(r1_table$VSS_TSS > 100, 100, r1_table$VSS_TSS)

r1_table %>% ggplot(aes(x=Date, y=percent_P_removal)) + geom_line()

r1_full <- r1_table %>% 
  filter_all(all_vars(!is.na(.)))

write.csv(r1_full, "cleaned_data/r1_cleaned_chemical_data.csv", quote=FALSE, row.names = FALSE)

r1_full %>% ggplot(aes(x=Date)) + geom_line(aes(y=percent_P_removal), colour="#4D3F83") + geom_line(aes(y=VSS_TSS), color="#D8820F") + theme_classic()

t1 <- as.Date("2010-01-28")
t2 <- as.Date("2013-05-30")
sub_r1 <- r1_table %>% filter(Date>=t1 & Date<=t2)
sub_r1_clean <- r1_full %>% filter(Date>=t1 & Date<=t2)

write.csv(sub_r1_clean, "cleaned_data/r1_subset_chemical_data_2010_2013.csv", quote=FALSE, row.names = FALSE)

sub_r1_clean %>% ggplot(aes(x=Date)) + geom_line(aes(y=percent_P_removal), colour="#4D3F83") + geom_line(aes(y=VSS_TSS), color="#D8820F") + theme_classic()

start_date <- as.Date("2010-01-28", format=("%Y-%m-%d"))
sub_r1_clean$operation_day <- difftime(sub_r1_clean$Date, start_date, units="days")
sub_r1_clean$SRT <- sub_r1_clean$operation_day / 4

r1_2011_sub_plot <- sub_r1_clean %>% ggplot(aes(x=operation_day)) + geom_line(aes(y=percent_P_removal), colour="#4D3F83", size=1) + geom_line(aes(y=VSS_TSS), color="#D8820F", size=1) + scale_x_continuous(breaks=seq(0,1220, by=50), expand=c(0,0), limits=c(0,1220)) + theme_classic()

r1_sub_SRT_plot <- sub_r1_clean %>% ggplot(aes(x=SRT)) + geom_line(aes(y=percent_P_removal), colour="#4D3F83", size=1.5) + geom_line(aes(y=VSS_TSS), color="#D8820F", size=1.5) + scale_x_continuous(breaks=seq(0,305, by=20), expand=c(0,0), limits=c(0,305)) + theme_classic()
r1_sub_SRT_plot

# p release vs p remaining values
sub_r1_clean %>% ggplot(aes(x=operation_day)) + geom_line(aes(y=P_release), colour="#4D3F83", size=1.5) + geom_line(aes(y=P_uptake), colour="#0FA1D8", size=1.5) + scale_x_continuous(breaks=seq(0,1220, by=50), expand=c(0,0), limits=c(0,1220)) + theme_classic()



####################################################
# R2
####################################################

r2_meta <- read_excel("metadata/r2-metadata-summary.xlsx")

r2_clean <- r2_meta %>% 
  select(Date, `P release`, `P uptake`, TSS, VSS) %>%
  filter(!is.na(`P uptake`)) %>% 
  mutate(VSS_TSS = (VSS / TSS) * 100)

colnames(r2_clean) <- c("Date", "P_release", "P_uptake", "TSS", "VSS", "VSS_TSS")
r2_clean$P_uptake <- ifelse(r2_clean$P_uptake < 0, 0, r2_clean$P_uptake)
r2_clean$P_release <- ifelse(r2_clean$P_release < 0, 0, r2_clean$P_release)
r2_table <- r2_clean %>% 
  mutate(percent_P_removal = ((P_release - P_uptake) / P_release) * 100) %>% 
  filter(!is.na(P_release)) %>% 
  filter(P_release != 0)

r2_table$percent_P_removal <- ifelse(r2_table$percent_P_removal < 0, 0, r2_table$percent_P_removal)
r2_table$Date <- as.Date(r2_table$Date, "YYYY-MM-DD")
r2_table$VSS_TSS <- ifelse(r2_table$VSS_TSS < 0, 0, r2_table$VSS_TSS)
r2_table$VSS_TSS <- ifelse(r2_table$VSS_TSS > 100, 100, r2_table$VSS_TSS)

r2_table %>% ggplot(aes(x=Date, y=percent_P_removal)) + geom_line()

r2_full <- r2_table %>% 
  filter_all(all_vars(!is.na(.)))

r2_full %>% ggplot(aes(x=Date)) + geom_line(aes(y=percent_P_removal), colour="#4D3F83") + geom_line(aes(y=VSS_TSS), color="#D8820F") + theme_classic()

####################################################
# qPCR Data
####################################################

qPCR <- read.csv("raw_data/ppk1_qPCR_data/Accum-qPCR-Timeseries.csv")
qPCR$Date <- as.Date(qPCR$Date, "%m/%d/%y")
t1 <- as.Date("2010-01-28")
t2 <- as.Date("2013-05-30")
sub_qPCR <- qPCR %>% filter(Date>=t1 & Date<=t2)
sub.m <- melt(sub_qPCR,id.vars="Date",measure.vars=c("Clade.IA", "Clade.IIA"))
sub.m$Date <- as.character(sub.m$Date)

qPCR_plot <- sub.m %>% ggplot(aes(x=Date, y=variable, fill=value)) + geom_tile(color="white") + scale_x_discrete(labels=sub.m$Date) + scale_fill_gradientn(colors = rev(viridis_pal()(9)), limits=c(0,16000000))
qPCR_plot
qPCR_clean <- qPCR_plot + theme(axis.line=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    legend.position="none",
                    panel.background=element_blank(),
                    panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank())
fig <- r1_2011_sub_plot / qPCR_clean
fig

ggsave("figs/R1_chemical_clade_dynamics.png", fig, width=20, height=8, units=c("cm"))


####################################################
# R1 sample inventories
####################################################

inventory <- read.csv("metadata/EBPR-SBR-R1-Sample-Inventory.csv", stringsAsFactors = FALSE)
inventory$Date <- as.Date(inventory$Date, format="%m/%d/%y")

colnames(inventory) <- c("Date", "Replicates", "Box", "Reactor", "Rack", "Existing", "Number_Replicates")

inventory_cleaned <- inventory %>% select(Date, Existing, Number_Replicates)

r1_inventory_data <- left_join(inventory_cleaned, r1_full) %>% filter_all(all_vars(!is.na(.)))

sub_r1_inventory_data <- r1_inventory_data %>%  filter(Date>=t1 & Date<=t2)

####################################################
# R1 qPCR and P data
####################################################

sub_r1_clean$P_release <- ifelse(sub_r1_clean$P_release > 100, 100, sub_r1_clean$P_release)

r1_clean_gathered <- sub_r1_clean %>% 
  select(operation_day, P_release, P_uptake) %>% 
  gather(key="variable", value="value", -operation_day)

qPCR_metadata <- left_join(sub_qPCR, sub_r1_clean) %>% filter_all(all_vars(!is.na(.)))

p_data <- r1_clean_gathered %>% ggplot(aes(x=operation_day, y=value)) + geom_line(aes(color=variable), size=1.4) + scale_color_manual(values=c("#4D3F83", "#0FA1D8"), labels=c("Total P Release at \n end of Anaerobic Phase \n", "Total P Remaining at \n end of Aerobic Phase")) + scale_x_continuous(breaks=seq(0,1220, by=30), expand=c(0,0), limits=c(0,1220)) + xlab("Operation Day") + ylab("Phosphorus (mg/L)") + theme_classic() + theme(axis.text.y=element_text(size=10), axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold"), legend.title=element_blank(), legend.position=c(0.15,0.9), legend.text=element_text(size=12),axis.text.x=element_text(size=11), axis.text.y.left=element_text(size=12))
p_data

clade_operation <- qPCR_metadata %>% select(Clade.IA, Clade.IIA, operation_day)

operation.m <- melt(clade_operation,id.vars="operation_day",measure.vars=c("Clade.IA", "Clade.IIA"))
operation.m$operation_day <- as.numeric(operation.m$operation_day)

clade_data <- operation.m %>% ggplot(aes(x=operation_day, y=value, group=variable)) + geom_point(size=1.4) + geom_line(aes(linetype=variable), size=1.4) + scale_linetype_manual(values=c("solid", "twodash"), labels=c("Clade IA", "Clade IIA")) + scale_x_continuous(breaks=seq(0,1220, by=30), expand=c(0,0), limits=c(0,1220)) + scale_y_continuous(limits=c(0,1.5e7), breaks=seq(0,1.5e7, 2.5e6), labels=c("0","2.5e+06", "5.0e+06", "7.5e+06", "1.0e+07", "1.25e+07", "1.5e+07")) + xlab("\nOperation Day") + ylab("Copies/ng DNA") + theme_classic() + theme(axis.title.x=element_text(size=14, face="bold"), axis.title.y=element_text(size=14, face="bold"), legend.title=element_blank(), legend.position=c(0.15,0.9), legend.text=element_text(size=12),axis.text.x=element_text(size=11), axis.text.y=element_text(size=12))

clade_data

r1_dynamics <- p_data / clade_data
r2 <- r1_dynamics + plot_annotation(tag_levels="A", tag_sep=".") & theme(plot.tag=element_text(size=16, face="bold", hjust=0.1, vjust=0.2))

ggsave("figs/R1_2010_2013_p_clade_data_modified.png", r2, width=17, height=10, units=c("in"))
