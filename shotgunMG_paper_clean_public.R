library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(vegan)
library(pheatmap)
library(tools)
library(lubridate)
library(stringr)
library(ggpubr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(ape)
library(ggVennDiagram)
library(ggpubr)
library(scales)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dendextend)

# Notes: In the following R code, variables labeled *china* refers to the human gut SM dataset
# which comes from Chinese participants. Variables labeled *rock* refers to the Antarctic SM dataset
# which comes from rocky environment.
# Soil refers to agricultural soil SM dataset. 
# Mock refers to both mock even and staggered communities SM dataset. 
# Mock even refers to the evenly staggered mock community SM dataset.
# Mock staggered refers to the staggered spiked-in mock community SM dataset.

# Define the directory where the zenodo archive has been downloaded and extracted.
# Change the following directory path to where the extracted zenodo archive lives:
# the working directory of the script is relative to the china_microbiome (gut microbiome) dataset, so please set your
# working directory accordingly.
setwd("/fs/vnas_Hnrc/eme/jut001/work/projects/shotgunMG/china_microbiome/")
# Set directory where output files will be written.
# The script expect the two R code files (shotgunMG_"paper_clean.R and shotgunMG_utils.R to be located in this (output) directory)
root = "~/projects/shotgunMG_paper/"
outdir = paste0(root)

options(stringsAsFactors = FALSE)

vColors = c(
  "#0000CD", "#00FF00", "#FF0000", "#808080", "#000000", "#B22222", "#DAA520", 
  "#DDA0DD", "#FF00FF", "#00FFFF", "#4682B4", "#E6E6FA", "#FF8C00", "#80008B", 
  "#8FBC8F", "#00BFFF", "#FFFF00", "#808000", "#FFCCCC", "#FFE5CC", "#FFFFCC", "#E5FFCC", 
  "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF", "#CCCCFF", "#E5CCFF", "#FFCCFF", "#FFCCE5", 
  "#FFFFFF", "#990000", "#666600", "#006666", "#330066", "#A0A0A0", "#99004C"
)
vColorsA = c(
  "#0000CD",  "#FF0000", "#808080", "#000000", "#B22222", "#DAA520", 
  "#DDA0DD", "#FF00FF", "#00FFFF", "#4682B4", "#E6E6FA", "#FF8C00", "#80008B", 
  "#8FBC8F", "#00BFFF", "#FFFF00", "#808000", "#FFCCCC", "#FFE5CC", "#FFFFCC", "#E5FFCC", 
  "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF", "#CCCCFF", "#E5CCFF", "#FFCCFF", "#FFCCE5", 
  "#FFFFFF", "#990000", "#666600", "#006666", "#330066", "#A0A0A0", "#99004C"
)
vColors2 = c(
  #red
  "#cc6155",
  "#e74c3c",
  "#a93225",
  #purple
  "#af7ac4",
  "#8d44ad",
  "#633974",
  #blue
  "#5599c7",
  "#2874a6",
  "#1f618d",
  #blue green
  "#48c9b0",
  "#169f85",
  "#148f77",
  #green
  "#26ae60",
  "#2fcb71",
  "#249b56",
  #yellow orange
  "#f5b041",
  "#e57e23",
  "#ca6f1d",
  "#9c640d",
  "#7d6608",
  # Gray blue-gray
  "#edbb99",
  "#5d6d7d",
  "#515a5a",
  "#abb2b9"
)

# This file contains functions that create data frames the multiple output
# files inherent to this study's design.
source(paste0(root, "./shotgunMG_utils.R"))

######################################
# Assembly stats                     #
# key metrics in assembly stats      #
# FIGURE S1                          #
######################################
mapping_file  = "./mapping_file.tsv"
mapping = data.frame(fread(mapping_file), check.names=FALSE)
row.names(mapping) = mapping$`#SampleID`
head(mapping)

df_stats = data.frame(fread("./assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = hms(df_stats$`compute time (hrs)`)
df_stats = df_stats[,c(1,2,4,8,9,10,11,12,13,14,15,17,18)]

colnames(df_stats)[4] = "N50 (bases)"
colnames(df_stats)[5] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("Number of bases (Gb)", "mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs", "Number of genes", "Number of assembled bases"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Human gut"
tmp = df_stats_melted_1[1:6,]
tmp$variable = "Number of sequencing clusters"
tmp$value = c(0.1, 0.5, 1, 4, 8, 12)
df_stats_melted_1 = rbind(df_stats_melted_1, tmp)

df_stats_melted_china = df_stats_melted_1

# Rock/Antarctic microbiome
df_stats = data.frame(fread("../rock_microbiome/assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = hms(df_stats$`compute time (hrs)`)
df_stats = df_stats[,c(1,2,4,8,9,10,11,12,13,14,15,17,18,19,20)]

colnames(df_stats)[4] = "N50 (bases)"
colnames(df_stats)[5] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs", "Number of genes", "Number of assembled bases"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Antarctic"
tmp = df_stats_melted_1[1:7,]
tmp$variable = "Number of sequencing clusters"
tmp$value = c(0.1, 0.5, 1, 4, 8, 12, 100)
df_stats_melted_rock = df_stats_melted_1
df_stats_melted_final = rbind(df_stats_melted_china, df_stats_melted_rock)
df_stats_melted_final = df_stats_melted_final[!df_stats_melted_final$variable %in% c("Number of sequencing clusters", "Number of bases (Gb)"), ]
unique(df_stats_melted_final$variable)

p_stats_1 <- ggplot(data=df_stats_melted_final, aes(x=.data[["Number of bases (Gb)"]], y=value, fill=Dataset, group=Dataset, color=Dataset)) + 
  facet_wrap(variable ~ ., ncol=3,  scales="free", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) + #, scales="free_x", space="free_x") + 
  geom_point(pch=21) + geom_line() + 
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title.y=element_text(family="Helvetica", size=0),
    axis.title.x=element_text(family="Helvetica", size=11),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=9, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=9, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA) +
  scale_y_continuous(labels = label_number(prefix = "", suffix = "", big.mark = ",")) +
  scale_x_continuous(labels = label_number(prefix = "", suffix = "", big.mark = ","))
p_stats_1

# Alternative version of the figure.
df_stats = data.frame(fread("./assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = as.numeric(hms(df_stats$`compute time (hrs)`))/60/60
df_stats = df_stats[,c(1,2,4,7,8,9,10,11,12,13,14,15,17,18)]
colnames(df_stats)[5] = "N50 (bases)"
colnames(df_stats)[6] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs", 
                                                                    "Number of genes", "Number of assembled bases", "compute time (hrs)"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Human gut"
df_stats_melted_china = df_stats_melted_1

# Rock/Antarctic microbiome
df_stats = data.frame(fread("../rock_microbiome/assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = as.numeric(hms(df_stats$`compute time (hrs)`))/60/60
df_stats = df_stats[,c(1,2,4,7,8,9,10,11,12,13,14,15,17,18,19,20)]
colnames(df_stats)[5] = "N50 (bases)"
colnames(df_stats)[6] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs", 
                                                                    "Number of genes", "Number of assembled bases", "compute time (hrs)"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Antarctic"
df_stats_melted_rock = df_stats_melted_1

# soil microbiome
df_stats = data.frame(fread("../soil/assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = as.numeric(hm(df_stats$`compute time (hrs)`))/60/60
df_stats = df_stats[,c(1,2,4,7,8,9,10,11,12,13,14,15,17,18)]
colnames(df_stats)[5] = "N50 (bases)"
colnames(df_stats)[6] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs", 
                                                                    "Number of genes", "Number of assembled bases", "compute time (hrs)"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Agricultural soil"
df_stats_melted_soil = df_stats_melted_1

# mock all
df_stats = data.frame(fread("../mock/assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = as.numeric(hm(df_stats$`compute time (hrs)`))/60/60
df_stats = df_stats[,c(1,2,4,7,8,9,10,11,12,13,14,15,17,18)]
colnames(df_stats)[5] = "N50 (bases)"
colnames(df_stats)[6] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs",
                                                                    "Number of genes", "Number of assembled bases", "compute time (hrs)"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Combined mock communities"
df_stats_melted_mock = df_stats_melted_1

# mock even
df_stats = data.frame(fread("../mock_even/assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = as.numeric(hm(df_stats$`compute time (hrs)`))/60/60
df_stats = df_stats[,c(1,2,4,7,8,9,10,11,12,13,14,15,17,18)]
colnames(df_stats)[5] = "N50 (bases)"
colnames(df_stats)[6] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs", 
                                                                    "Number of genes", "Number of assembled bases", "compute time (hrs)"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Even mock community"
df_stats_melted_mock_even = df_stats_melted_1

# mock staggered
df_stats = data.frame(fread("../mock_staggered/assembly_stats.csv", sep="\t", quote=FALSE), check.names=F)
df_stats$`mem (GB)` = as.numeric(gsub(",", "", df_stats$`mem (GB)`))
df_stats$`Number of bases (Gb)` = as.numeric(gsub(",", "", df_stats$`Number of bases (Gb)`))
df_stats$N50 = as.numeric(gsub(",", "", df_stats$N50))
df_stats$`max contig length` = as.numeric(gsub(",", "", df_stats$`max contig length`))
df_stats$`Number of assembled bases` = as.numeric(gsub(",", "", df_stats$`number of assembled bases`))
df_stats$`number of assembled bases` = NULL
df_stats$`Number of contigs` = as.numeric(gsub(",", "", df_stats$`Number of contigs`))
df_stats$`Number of genes` = as.numeric(gsub(",", "", df_stats$`Number of genes`))
df_stats$gt10 = as.numeric(gsub(",", "", df_stats$gt10))
df_stats$gt20 = as.numeric(gsub(",", "", df_stats$gt20))
df_stats$gt40 = as.numeric(gsub(",", "", df_stats$gt40))
df_stats$gt80 = as.numeric(gsub(",", "", df_stats$gt80))
df_stats$gt160 = as.numeric(gsub(",", "", df_stats$gt160))
df_stats$gt320 = as.numeric(gsub(",", "", df_stats$gt320))
df_stats$`compute time (hrs)` = as.numeric(hm(df_stats$`compute time (hrs)`))/60/60
df_stats = df_stats[,c(1,2,4,7,8,9,10,11,12,13,14,15,17,18)]
colnames(df_stats)[5] = "N50 (bases)"
colnames(df_stats)[6] = "max contigs length (bases)"
df_stats_melted = melt(df_stats, id.vars=c("Number of bases (Gb)"), measure.vars=colnames(df_stats[,2:(ncol(df_stats))]))
df_stats_melted_1 = df_stats_melted[df_stats_melted$variable %in% c("mem (GB)", "N50 (bases)", "max contig length (bases)", "Number of contigs", 
                                                                    "Number of genes", "Number of assembled bases", "compute time (hrs)"),]
df_stats_melted_1$value = as.numeric(df_stats_melted_1$value)
df_stats_melted_1$Dataset = "Staggered mock community"
df_stats_melted_mock_staggered = df_stats_melted_1

df_stats_melted_final = rbind(df_stats_melted_china, df_stats_melted_rock, 
                              df_stats_melted_soil, df_stats_melted_mock,
                              df_stats_melted_mock_even, df_stats_melted_mock_staggered)
#df_stats_melted_final = df_stats_melted_final[df_stats_melted_final$variable != "Number of bases (Gb)", ]
unique(df_stats_melted_final$variable)
#df_stats_melted_final$`Number of clusters (M)` = factor(df_stats_melted_final$`Number of clusters (M)`, levels=c("0.1", "0.5", "1", "4", "8", "12", "all data"))


df_stats_melted_final_high = df_stats_melted_final[df_stats_melted_final$Dataset %in% c("Human gut", "Antarctic", "Agricultural soil"),]
df_stats_melted_final_low = df_stats_melted_final[df_stats_melted_final$Dataset %in% c("Combined mock communities", "Even mock community", "Staggered mock community"),]

p_stats_3a <- ggplot(data=df_stats_melted_final_high, aes(x=.data[["Number of bases (Gb)"]], y=value, group=Dataset)) + 
  facet_grid(variable ~ Dataset,  scales="free", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) + #, scales="free_x", space="free_x") + 
  geom_point(pch=21) + geom_line() + 
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title.y=element_text(family="Helvetica", size=0),
    axis.title.x=element_text(family="Helvetica", size=11),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, hjust=0, vjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA) +
  scale_y_continuous(labels = label_number(prefix = "", suffix = "", big.mark = ","))

p_stats_3b <- ggplot(data=df_stats_melted_final_low, aes(x=.data[["Number of bases (Gb)"]], y=value, group=Dataset)) + 
  facet_grid(variable ~ Dataset,  scales="free", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) + #, scales="free_x", space="free_x") + 
  geom_point(pch=21) + geom_line() + 
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title.y=element_text(family="Helvetica", size=0),
    axis.title.x=element_text(family="Helvetica", size=11),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, hjust=0, vjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA) +
  scale_y_continuous(labels = label_number(prefix = "", suffix = "", big.mark = ","))

# Assemble the figure(s):
figure2a <- ggarrange(
  p_stats_3a, p_stats_3b,
  labels = NULL,
  ncol = 2, nrow = 1, widths = c(0.5, 0.5),
  common.legend = TRUE, legend = "none",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f2a = annotate_figure(figure2a, 
                     fig.lab = "\n    ", 
                     fig.lab.pos = "top.left",
                     top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)


final_list = annotate_figure(
  ggarrange(f2a
            #labels=c("A", "B"),#, "C", "D", "E"),
            #nrow = 2,
            #heights=c(0.45, 0.55)
  ),
  top = text_grob("Figure S1\n", color = "black", face = "plain", size = 11),
)
final_list
pdf( file=paste0(outdir, "./figures/FIGURE_S1.pdf"), height=8, width=10)
final_list
dev.off()


################################################
# qc mapping stats                             #
# For each qc mapping stats file, put          # 
# properly aligned reads in a df and boxplot.  #
# Figure S2                                    #
################################################    
df_qc = get_qc_as_df(my_qc_files_china)
curr_order_facet = as.character(unique(df_qc$Type))

df_qc$Type = factor(df_qc$Type, levels=curr_order_facet)
df_qc$dummy = "dummy"
df_qc$num_clusters = df_qc$Type
df_qc$num_clusters = gsub("_all", "", df_qc$Type)
df_qc$num_clusters = factor(df_qc$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))
df_qc$Type = gsub("_all", " - arm", df_qc$Type)
df_qc$properlyPaired = gsub(",", "", df_qc$properlyPaired)
df_qc$properlyPaired = as.numeric(df_qc$properlyPaired)
df_qc$`properlyPaired%` = as.numeric(df_qc$`properlyPaired%`)

p_qc1 <- ggplot(data=df_qc, aes(x=Type, y=properlyPaired)) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") +
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("Number of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

p_qc2 <- ggplot(data=df_qc, aes(x=Type, y=.data[["properlyPaired%"]])) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("% of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

figure3a <- ggarrange(
                    p_qc1, p_qc2,
                    labels = NULL,
                    ncol = 2, nrow = 1,
                    common.legend = TRUE, legend = "bottom",
                    align = "hv", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f3a = annotate_figure(figure3a, 
                     fig.lab = "     Human gut", 
                     fig.lab.pos = "top.left",
                     top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# Rock microbiome
df_qc = NULL
df_qc = get_qc_as_df(my_qc_files_rock)
# order
curr_order_facet = as.character(unique(df_qc$Type))
df_qc$Type = factor(df_qc$Type, levels=curr_order_facet)
df_qc$dummy = "dummy"
df_qc$num_clusters = df_qc$Type
df_qc$num_clusters = gsub("_all", "", df_qc$Type)
df_qc$num_clusters = gsub("complete", "", df_qc$num_clusters)
df_qc$num_clusters = factor(df_qc$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", ""))
df_qc$Type = gsub("_all", " - arm", df_qc$Type)
df_qc$properlyPaired = gsub(",", "", df_qc$properlyPaired)
df_qc$properlyPaired = as.numeric(df_qc$properlyPaired)
df_qc$`properlyPaired%` = as.numeric(df_qc$`properlyPaired%`)

p_qc3 <- ggplot(data=df_qc, aes(x=Type, y=properlyPaired)) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("Number of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

p_qc4 <- ggplot(data=df_qc, aes(x=Type, y=.data[["properlyPaired%"]])) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("% of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

figure3b <- ggarrange(
  p_qc3, p_qc4,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f3b = annotate_figure(figure3b, 
                     fig.lab = "     Antarctic", 
                     fig.lab.pos = "top.left",
                     top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# Agricultural soil
df_qc = NULL
df_qc = get_qc_as_df(my_qc_files_soil)
# order
curr_order_facet = as.character(unique(df_qc$Type))
df_qc$Type = factor(df_qc$Type, levels=curr_order_facet)
df_qc$dummy = "dummy"
df_qc$num_clusters = df_qc$Type
df_qc$num_clusters = gsub("_all", "", df_qc$Type)
df_qc$num_clusters = gsub("complete", "", df_qc$num_clusters)
df_qc$num_clusters = factor(df_qc$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", ""))
df_qc$Type = gsub("_all", " - arm", df_qc$Type)
df_qc$properlyPaired = gsub(",", "", df_qc$properlyPaired)
df_qc$properlyPaired = as.numeric(df_qc$properlyPaired)
df_qc$`properlyPaired%` = as.numeric(df_qc$`properlyPaired%`)

p_qc5 <- ggplot(data=df_qc, aes(x=Type, y=properlyPaired)) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("Number of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

p_qc6 <- ggplot(data=df_qc, aes(x=Type, y=.data[["properlyPaired%"]])) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("% of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

figure3c <- ggarrange(
  p_qc5, p_qc6,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f3c = annotate_figure(figure3c, 
                      fig.lab = "     Agricultural soil", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# mock microbiome
df_qc = NULL
df_qc = get_qc_as_df(my_qc_files_mock)
# order
curr_order_facet = as.character(unique(df_qc$Type))
df_qc$Type = factor(df_qc$Type, levels=curr_order_facet)
df_qc$dummy = "dummy"
df_qc$num_clusters = df_qc$Type
df_qc$num_clusters = gsub("_all", "", df_qc$Type)
df_qc$num_clusters = gsub("complete", "", df_qc$num_clusters)
df_qc$num_clusters = factor(df_qc$num_clusters, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", ""))
df_qc$Type = gsub("_all", " - arm", df_qc$Type)
df_qc$properlyPaired = gsub(",", "", df_qc$properlyPaired)
df_qc$properlyPaired = as.numeric(df_qc$properlyPaired)
df_qc$`properlyPaired%` = as.numeric(df_qc$`properlyPaired%`)

p_qc7 <- ggplot(data=df_qc, aes(x=Type, y=properlyPaired)) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("Number of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

p_qc8 <- ggplot(data=df_qc, aes(x=Type, y=.data[["properlyPaired%"]])) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("% of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

figure3d <- ggarrange(
  p_qc7, p_qc8,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f3d = annotate_figure(figure3d, 
                      fig.lab = "     Combined mock communities", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# mock microbiome even
df_qc = NULL
df_qc = get_qc_as_df(my_qc_files_mock_even)
# order
curr_order_facet = as.character(unique(df_qc$Type))
df_qc$Type = factor(df_qc$Type, levels=curr_order_facet)
df_qc$dummy = "dummy"
df_qc$num_clusters = df_qc$Type
df_qc$num_clusters = gsub("_all", "", df_qc$Type)
df_qc$num_clusters = gsub("complete", "", df_qc$num_clusters)
df_qc$num_clusters = factor(df_qc$num_clusters, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", ""))
df_qc$Type = gsub("_all", " - arm", df_qc$Type)
df_qc$properlyPaired = gsub(",", "", df_qc$properlyPaired)
df_qc$properlyPaired = as.numeric(df_qc$properlyPaired)
df_qc$`properlyPaired%` = as.numeric(df_qc$`properlyPaired%`)

p_qc9 <- ggplot(data=df_qc, aes(x=Type, y=properlyPaired)) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("Number of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

p_qc10 <- ggplot(data=df_qc, aes(x=Type, y=.data[["properlyPaired%"]])) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("% of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

figure3e <- ggarrange(
  p_qc9, p_qc10,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f3e = annotate_figure(figure3e, 
                      fig.lab = "     Even mock community", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# mock microbiome staggered
df_qc = NULL
df_qc = get_qc_as_df(my_qc_files_mock_staggered)
# order
curr_order_facet = as.character(unique(df_qc$Type))
df_qc$Type = factor(df_qc$Type, levels=curr_order_facet)
df_qc$dummy = "dummy"
df_qc$num_clusters = df_qc$Type
df_qc$num_clusters = gsub("_all", "", df_qc$Type)
df_qc$num_clusters = gsub("complete", "", df_qc$num_clusters)
df_qc$num_clusters = factor(df_qc$num_clusters, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", ""))
df_qc$Type = gsub("_all", " - arm", df_qc$Type)
df_qc$properlyPaired = gsub(",", "", df_qc$properlyPaired)
df_qc$properlyPaired = as.numeric(df_qc$properlyPaired)
df_qc$`properlyPaired%` = as.numeric(df_qc$`properlyPaired%`)

p_qc11 <- ggplot(data=df_qc, aes(x=Type, y=properlyPaired)) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("Number of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

p_qc12 <- ggplot(data=df_qc, aes(x=Type, y=.data[["properlyPaired%"]])) + 
  facet_grid(. ~ num_clusters, scales="free_x", space="free_x") + 
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, notch=FALSE, color="black", fill="gray") + 
  geom_jitter(size=0.1, position=position_jitter(0.2), alpha=0.4) +
  ylab(str_wrap("% of properly aligned paired reads", width = 20)) + xlab(paste0("")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))

figure3f <- ggarrange(
  p_qc11, p_qc12,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f3f = annotate_figure(figure3f, 
                      fig.lab = "     Staggered mock community", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

final_list = annotate_figure(
  ggarrange(f3a,f3b,f3c,f3d,f3e,f3f,
            labels=c("A", "B", "C", "D", "E", "F"),
            nrow = 6,
            heights=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  ),
  top = text_grob("Figure S2\n", color = "black", face = "plain", size = 11),
)
final_list

# Labels of the scatter plot
pdf( file=paste0(outdir, "./figures/FIGURE_S2.pdf"), height=14, width=6.8)
final_list
dev.off()

##############################################################
# Beta diversity                                             #
#                                                            #
# Generate Mantel correlation between all distance matrices  #
# The goal is to see if there are difference overall between #
# all the different pipeline methods.                        #
# Then plot the correlations on a heatmap for visual support.#
# Also plot PcoA plot with alpha diversity quitilies         #
# Figure 2                                                   #
##############################################################
## China microbiome
df_list = get_dists_as_list(my_dist_files, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files), ncol = length(my_dist_files)))
colnames(df2) = names(my_dist_files)
row.names(df2) = names(my_dist_files)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
   
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

# Write the table because it takes several hours to generate all the spearman correlation values.
#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_betadiv_df.tsv"), sep="\t"), check.names=FALSE)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

p_beta1 = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) + 
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta1)

## Rock microbiome
df_list = get_dists_as_list(my_dist_files_rock, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_rock), ncol = length(my_dist_files_rock)))
colnames(df2) = names(my_dist_files_rock)
row.names(df2) = names(my_dist_files_rock)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta2 = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta2)

## Agricultural soil
df_list = get_dists_as_list(my_dist_files_soil, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_rock), ncol = length(my_dist_files_soil)))
colnames(df2) = names(my_dist_files_soil)
row.names(df2) = names(my_dist_files_soil)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta3 = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta3)

## mock all
df_list = get_dists_as_list(my_dist_files_mock, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_mock), ncol = length(my_dist_files_mock)))
colnames(df2) = names(my_dist_files_mock)
row.names(df2) = names(my_dist_files_mock)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M", y="2M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M", y="3M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M_all", y="2M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M_all", y="3M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta4 = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta4)

## mock even
df_list = get_dists_as_list(my_dist_files_mock_even, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_mock_even), ncol = length(my_dist_files_mock_even)))
colnames(df2) = names(my_dist_files_mock_even)
row.names(df2) = names(my_dist_files_mock_even)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M", y="2M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M", y="3M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M_all", y="2M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M_all", y="3M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta5 = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta5)

## mock staggered
df_list = get_dists_as_list(my_dist_files_mock_staggered, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_mock_staggered), ncol = length(my_dist_files_mock_staggered)))
colnames(df2) = names(my_dist_files_mock_staggered)
row.names(df2) = names(my_dist_files_mock_staggered)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 3, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M", y="2M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M", y="3M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M_all", y="2M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M_all", y="3M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta6 = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta6)

#splitted mock community too small to generate mantel analyses

##################################################
# Beta div tile maps with genes                  #
#                                                #
##################################################
# Do similar analyses, but based on gene abundance distance matrices.
## China microbiome
df_list = get_dists_as_list(my_dist_files_genes, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files), ncol = length(my_dist_files)))
colnames(df2) = names(my_dist_files)
row.names(df2) = names(my_dist_files)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

# Write the table because it takes several hours to generate all the spearman correlation values.
#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_genes_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_betadiv_genes_df.tsv"), sep="\t"), check.names=FALSE)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

p_beta1g = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) + 
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta1g)

## Rock microbiome
df_list = get_dists_as_list(my_dist_files_rock_genes, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_rock_genes), ncol = length(my_dist_files_rock_genes)))
colnames(df2) = names(my_dist_files_rock_genes)
row.names(df2) = names(my_dist_files_rock_genes)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta2g = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta2g)

## Agricultural soil
df_list = get_dists_as_list(my_dist_files_soil_genes, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_soil_genes), ncol = length(my_dist_files_soil_genes)))
colnames(df2) = names(my_dist_files_soil_genes)
row.names(df2) = names(my_dist_files_soil_genes)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_soil_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta3g = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta3g)

## mock all
df_list = get_dists_as_list(my_dist_files_mock_genes, debug=FALSE)
df = NULL
df2 = data.frame(matrix(1, nrow = length(my_dist_files_mock_genes), ncol = length(my_dist_files_mock_genes)))
colnames(df2) = names(my_dist_files_mock_genes)
row.names(df2) = names(my_dist_files_mock_genes)
df3 = NULL
x = 1
for(i in 1:length(df_list)){
  if(i == length(df_list)){ break; }
  for(j in (i+1):length(df_list)){
    df_x = df_list[[i]]
    df_y = df_list[[j]]
    name_x = names(df_list[i])
    name_y = names(df_list[j])
    res = mantel(df_x, df_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    df2[[name_x, name_y]] = res$statistic
    df2[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M", y="2M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M", y="3M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M_all", y="2M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M_all", y="3M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm"))

#write.table(df3, paste0(root, "/melted_style_spearman_betadiv_df_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 4M, 0.1M - arm, 0.5M - arm and 1M - arm
p_beta4g = ggplot(df3, aes(x=y, y=x,label=round(r,2))) +
  geom_tile(aes(fill=r)) +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_beta4g)



# Figure 2 will be the tiled heatmaps and Fig S3 will be the ordinations.

# Tiled heatmaps
#A)
p_beta1_ann = annotate_figure(p_beta1, 
                      fig.lab = "                           Contigs", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
p_beta1g_ann = annotate_figure(p_beta1g, 
                           fig.lab = "                       Genes", 
                           fig.lab.pos = "top.left",
                           top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

figure4a <- annotate_figure(
  ggarrange(
    p_beta1_ann, p_beta1g_ann,
    labels = c("", ""),
    ncol = 2, nrow = 1,
    widths = c(0.5,0.5),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n       Human gut", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

#B)
p_beta2_ann = annotate_figure(p_beta2, 
                              fig.lab = "                           Contigs", 
                              fig.lab.pos = "top.left",
                              top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
p_beta2g_ann = annotate_figure(p_beta2g, 
                               fig.lab = "                       Genes", 
                               fig.lab.pos = "top.left",
                               top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

figure4b <- annotate_figure(
  ggarrange(
    p_beta2_ann, p_beta2g_ann,
    labels = c("", ""),
    ncol = 2, nrow = 1,
    widths = c(0.5,0.5),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n       Antarctic", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

#C)
p_beta3_ann = annotate_figure(p_beta3, 
                              fig.lab = "                           Contigs", 
                              fig.lab.pos = "top.left",
                              top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
p_beta3g_ann = annotate_figure(p_beta3g, 
                               fig.lab = "                       Genes", 
                               fig.lab.pos = "top.left",
                               top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
figure4c <- annotate_figure(
  ggarrange(
    p_beta3_ann, p_beta3g_ann,
    labels = c("", ""),
    ncol = 2, nrow = 1,
    widths = c(0.5,0.5),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n       Agricultural soil", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

#D)
p_beta4_ann = annotate_figure(p_beta4, 
                              fig.lab = "                           Contigs", 
                              fig.lab.pos = "top.left",
                              top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
p_beta4g_ann = annotate_figure(p_beta4g, 
                               fig.lab = "                       Genes", 
                               fig.lab.pos = "top.left",
                               top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
figure4d <- annotate_figure(
  ggarrange(
    p_beta4_ann, p_beta4g_ann,
    labels = c(""),
    ncol = 2, nrow = 1,
    widths = c(0.5,0.5),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n       Combined mock communities", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

## Final
final_list = annotate_figure(
  ggarrange(figure4a, figure4b, figure4c, figure4d,
            labels=c("A", "B", "C", "D"),#, "C", "D", "E"),
            nrow = 4,ncol=1,
            heights=c(0.5, 0.5, 0.5, 0.5)
  ),
  top = text_grob("Figure 2\n", color = "black", face = "plain", size = 11),
)
#final_list
pdf( file=paste0(outdir, "./figures/F_FIGURE_2.pdf"), height=18, width=12)
final_list
dev.off()

# pcoas
##################################
# PCoAs                          #
#                                #
##################################

################################
# Human gut CONTIGS
mapping = data.frame(fread("./mapping_file3.tsv", sep="\t"), check.names=FALSE)

# Beta div contigs
# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diversity
df_coords = get_coords_as_df(my_coords_files)
colnames(df_coords)[5] = "workflow"
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files, mapping=mapping)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all")

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
df_coords$D2 = ifelse(df_coords$workflow %in% c("4M", "0.1M_all", "0.5M_all", "1M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bd = list()
# for simplicity sake, only show quartile 1 and 5
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  curr_ad = df_alpha[df_alpha$Type == curr_workflow,]
  
  # Divide by quintiles:
  curr_ad$quantiles = NULL
  curr_ad$quantiles = cut(curr_ad$value, 
                          quantile(
                            curr_ad$value, 
                            prob = seq(0, 1, length = 6), 
                            type = 5, names=TRUE, 
                            labels=FALSE, 
                            ordered_result=TRUE
                          )
  )
  
  df_link = data.frame(quantiles=levels(curr_ad$quantiles), rank=seq(1,(length(unique(curr_ad$quantiles))-1), by=1))
  df_tmp = curr_ad[,c("variable", "quantiles")]
  df_tmp = merge(df_tmp, df_link, by="quantiles")
  curr_mapping = merge(mapping, df_tmp, by.x="#SampleID", by.y="variable")
  colnames(curr_mapping)[ncol(curr_mapping)] = "Diversity_quartile"
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  curr_bd$Diversity_quartile = as.character(curr_bd$Diversity_quartile)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Diversity_quartile)) +
    geom_point(size=0.95, shape=21, stroke=0) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) +
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa) +
    guides(fill = guide_legend(title="Alpha diversity\n(richness) quintiles",override.aes = list(size=3)))
  ps_bd[[i]] = p
  
}

# Human gut GENES
mapping = data.frame(fread("./mapping_file3.tsv", sep="\t"), check.names=FALSE)

# Beta div contigs
# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diversity
df_coords = get_coords_as_df(my_coords_files_genes)
colnames(df_coords)[5] = "workflow"
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files, mapping=mapping)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all")

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
df_coords$D1 = ifelse(df_coords$workflow %in% c("0.1M", "0.5M", "0.1M_all"), df_coords$D1 * -1, df_coords$D1)
df_coords$D2 = ifelse(df_coords$workflow %in% c("0.1M", "1M", "0.5M_all", "1M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bd_genes = list()
# for simplicity sake, only show quartile 1 and 5
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  curr_ad = df_alpha[df_alpha$Type == curr_workflow,]
  
  # Divide by quintiles:
  curr_ad$quantiles = NULL
  curr_ad$quantiles = cut(curr_ad$value, 
                          quantile(
                            curr_ad$value, 
                            prob = seq(0, 1, length = 6), 
                            type = 5, names=TRUE, 
                            labels=FALSE, 
                            ordered_result=TRUE
                          )
  )
  
  df_link = data.frame(quantiles=levels(curr_ad$quantiles), rank=seq(1,(length(unique(curr_ad$quantiles))-1), by=1))
  df_tmp = curr_ad[,c("variable", "quantiles")]
  df_tmp = merge(df_tmp, df_link, by="quantiles")
  curr_mapping = merge(mapping, df_tmp, by.x="#SampleID", by.y="variable")
  colnames(curr_mapping)[ncol(curr_mapping)] = "Diversity_quartile"
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  curr_bd$Diversity_quartile = as.character(curr_bd$Diversity_quartile)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Diversity_quartile)) +
    geom_point(size=0.95, shape=21, stroke=0) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) +
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa) +
    guides(fill = guide_legend(title="Alpha diversity\n(richness) quintiles",override.aes = list(size=3)))
  ps_bd_genes[[i]] = p
  
}

#########################################
# Do PCoA plots for Antarctic microbiome CONTIGS
mapping = data.frame(fread("../rock_microbiome/0.1M_clusters/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
df_coords = get_coords_as_df(my_coords_files_rock)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
vColorsPcoa2 = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00",
                 "#B22222", "#DAA520", "#DDA0DD", "#FF00FF",
                 "#00FFFF", "#4682B4", "#80008B")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all")
curr_mapping = mapping

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
colnames(df_coords)[5] = "workflow"
df_coords$D1 = ifelse(df_coords$workflow %in% c("complete", "8M_all", "12M_all"), df_coords$D1 * -1, df_coords$D1)
df_coords$D2 = ifelse(df_coords$workflow %in% c("0.5M", "8M", "complete", "4M_all", "8M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bdr = list()

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 0.5M, 8M, complete, 4M - arm, 8M - arm
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_bd$Treatment = gsub(".nord", "", curr_bd$Treatment)
  curr_bd$Treatment = gsub(".sud", "", curr_bd$Treatment)
  
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Treatment)) +
    geom_point(size=2.0, shape=21, stroke=0.5) + #geom_point(colour="grey90", size = 1.5) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) + 
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa2) +
    guides(fill = guide_legend(title="Sampling location",override.aes = list(size=3)))
  ps_bdr[[i]] = p
  
}

# Do PCoA plots for Antarctic microbiome GENES
mapping = data.frame(fread("../rock_microbiome/0.1M_clusters/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
df_coords = get_coords_as_df(my_coords_files_rock_genes)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
vColorsPcoa2 = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00",
                 "#B22222", "#DAA520", "#DDA0DD", "#FF00FF",
                 "#00FFFF", "#4682B4", "#80008B")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all")
curr_mapping = mapping

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
colnames(df_coords)[5] = "workflow"
#df_coords$D1 = ifelse(df_coords$workflow %in% c("0.1M", "8M", "complete", "4M_all", "8M_all"), df_coords$D1 * -1, df_coords$D1)
df_coords$D2 = ifelse(df_coords$workflow %in% c("0.5M", "1M", "4M", "12M", "complete", "0.1M_all", "0.5M_all", "4M_all", "8M_all", "12M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bdr_genes = list()

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 0.5M, 8M, complete, 4M - arm, 8M - arm
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_bd$Treatment = gsub(".nord", "", curr_bd$Treatment)
  curr_bd$Treatment = gsub(".sud", "", curr_bd$Treatment)
  
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Treatment)) +
    geom_point(size=2.0, shape=21, stroke=0.5) + #geom_point(colour="grey90", size = 1.5) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) + 
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa2) +
    guides(fill = guide_legend(title="Sampling location",override.aes = list(size=3)))
  ps_bdr_genes[[i]] = p
  
}

#########################################
# Do PCoA plots for Agr soil CONTIGS

mapping = data.frame(fread("../soil/0.1M_clusters/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
df_coords = get_coords_as_df(my_coords_files_soil)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
vColorsPcoa2 = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00",
                 "#B22222", "#DAA520", "#DDA0DD", "#FF00FF",
                 "#00FFFF", "#4682B4", "#80008B")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all")
curr_mapping = mapping

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
colnames(df_coords)[5] = "workflow"
df_coords$D1 = ifelse(df_coords$workflow %in% c("4M_all"), df_coords$D1 * -1, df_coords$D1)
df_coords$D2 = ifelse(df_coords$workflow %in% c("0.1M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bds = list()

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 0.5M, 8M, complete, 4M - arm, 8M - arm
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_bd$Nitrogen = gsub(".nord", "", curr_bd$Nitrogen)
  curr_bd$Nitrogen = gsub(".sud", "", curr_bd$Nitrogen)
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Nitrogen)) +
    geom_point(size=2.0, shape=21, stroke=0.5) + #geom_point(colour="grey90", size = 1.5) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) + 
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa2) +
    guides(fill = guide_legend(title="Nitrogen",override.aes = list(size=3)))
  ps_bds[[i]] = p
  
}

# Do PCoA plots for Agr soil GENES

mapping = data.frame(fread("../soil/0.1M_clusters/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
df_coords = get_coords_as_df(my_coords_files_soil_genes)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
vColorsPcoa2 = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00",
                 "#B22222", "#DAA520", "#DDA0DD", "#FF00FF",
                 "#00FFFF", "#4682B4", "#80008B")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all")
curr_mapping = mapping

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
colnames(df_coords)[5] = "workflow"
df_coords$D2 = ifelse(df_coords$workflow %in% c("0.1M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bds_genes = list()

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 0.5M, 8M, complete, 4M - arm, 8M - arm
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_bd$Nitrogen = gsub(".nord", "", curr_bd$Nitrogen)
  curr_bd$Nitrogen = gsub(".sud", "", curr_bd$Nitrogen)
  
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Nitrogen)) +
    geom_point(size=2.0, shape=21, stroke=0.5) + #geom_point(colour="grey90", size = 1.5) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) + 
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa2) +
    guides(fill = guide_legend(title="Nitrogen",override.aes = list(size=3)))
  ps_bds_genes[[i]] = p
  
}

#########################################
# Do PCoA plots for combined mock community CONTIGS

mapping = data.frame(fread("../mock/0.1M_clusters/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
df_coords = get_coords_as_df(my_coords_files_mock)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
vColorsPcoa2 = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00",
                 "#B22222", "#DAA520", "#DDA0DD", "#FF00FF",
                 "#00FFFF", "#4682B4", "#80008B")
my_workflows = c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M_all", "0.5M_all", "1M_all", "2M_all", "3M_all", "4M_all")
curr_mapping = mapping

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
colnames(df_coords)[5] = "workflow"
#df_coords$D2 = ifelse(df_coords$workflow %in% c("0.5M", "8M", "complete", "4M_all", "8M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bdm = list()

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 0.5M, 8M, complete, 4M - arm, 8M - arm
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  #curr_bd$Nitrogen = gsub(".nord", "", curr_bd$Nitrogen)
  #curr_bd$Nitrogen = gsub(".sud", "", curr_bd$Nitrogen)
  
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Treatment)) +
    geom_point(size=2.0, shape=21, stroke=0.5) + #geom_point(colour="grey90", size = 1.5) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) + 
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa2) +
    guides(fill = guide_legend(title="Community type",override.aes = list(size=3)))
  ps_bdm[[i]] = p
  
}

# Do PCoA plots for combined mock community GENES
mapping = data.frame(fread("../mock/0.1M_clusters/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
df_coords = get_coords_as_df(my_coords_files_mock_genes)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
vColorsPcoa2 = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00",
                 "#B22222", "#DAA520", "#DDA0DD", "#FF00FF",
                 "#00FFFF", "#4682B4", "#80008B")
my_workflows = c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete", "0.1M_all", "0.5M_all", "1M_all", "2M_all", "3M_all", "4M_all")
curr_mapping = mapping

# Invert axes of 4M, 0.1M -arm, 0.5M -arm and 1M -arm
colnames(df_coords)[5] = "workflow"
#df_coords$D2 = ifelse(df_coords$workflow %in% c("0.5M", "8M", "complete", "4M_all", "8M_all"), df_coords$D2 * -1, df_coords$D2)

ps_bdm_genes = list()

# Since signs of axes do not mean anything in PcoAs, invert y-axis of 0.5M, 8M, complete, 4M - arm, 8M - arm
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_bd = df_coords[df_coords[,5] == curr_workflow,]
  
  #pcoa
  curr_bd = merge(curr_bd, curr_mapping, by.x="variable", by.y="#SampleID")
  percent1 = unique(curr_bd$PCo1)
  percent2 = unique(curr_bd$PCo2)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  #curr_bd$Nitrogen = gsub(".nord", "", curr_bd$Nitrogen)
  #curr_bd$Nitrogen = gsub(".sud", "", curr_bd$Nitrogen)
  
  
  p <- ggplot(data=curr_bd, aes(x=D1, y=D2, fill=Treatment)) +
    geom_point(size=2.0, shape=21, stroke=0.5) + #geom_point(colour="grey90", size = 1.5) +
    ggtitle(curr_workflow2) + 
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_text(size=8, colour="black"),
      axis.text.y=element_text(size=8, colour="black"),
      axis.title=element_text(size=8),
      axis.ticks.length=unit(0.2,"cm"),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.key=element_blank(),
      legend.margin = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5, size=8, vjust=0),
      strip.background =  element_blank(),
      strip.text.x = element_text(size=8)
    ) + 
    geom_hline(aes(yintercept=0), size=0.2) +
    xlab(paste0("PCo1 (", percent1,"%)")) +
    ylab(paste0("PCo2 (", percent2,"%)")) + 
    scale_fill_manual(values=vColorsPcoa2) +
    guides(fill = guide_legend(title="Community type",override.aes = list(size=3)))
  ps_bdm_genes[[i]] = p
  
}

#A)
figure_pcoa1 = do.call("grid_arrange_shared_legend", c(ps_bd, ncol=6, nrow=2, position="right"))
figure_pcoa1 <- annotate_figure(
  ggarrange(
    figure_pcoa1,
    labels = c("A\n"),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Human gut; contigs", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
figure_pcoa1g = do.call("grid_arrange_shared_legend", c(ps_bd_genes, ncol=6, nrow=2, position="right"))
figure_pcoa1g <- annotate_figure(
  ggarrange(
    figure_pcoa1g,
    labels = c(""),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Human gut; genes", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# B)
figure_pcoa2 = do.call("grid_arrange_shared_legend", c(ps_bdr, ncol=7, nrow=2, position="right"))
figure_pcoa2 <- annotate_figure(
  ggarrange(
    figure_pcoa2,
    labels = c("B\n"),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Arctic soil; contigs", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
figure_pcoa2g = do.call("grid_arrange_shared_legend", c(ps_bdr_genes, ncol=7, nrow=2, position="right"))
figure_pcoa2g <- annotate_figure(
  ggarrange(
    figure_pcoa2g,
    labels = c(""),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Arctic soil; genes", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# C)
figure_pcoa3 = do.call("grid_arrange_shared_legend", c(ps_bds, ncol=7, nrow=2, position="right"))
figure_pcoa3 <- annotate_figure(
  ggarrange(
    figure_pcoa3,
    labels = c("C\n"),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Agricultural soil; contigs", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
figure_pcoa3g = do.call("grid_arrange_shared_legend", c(ps_bds_genes, ncol=7, nrow=2, position="right"))
figure_pcoa3g <- annotate_figure(
  ggarrange(
    figure_pcoa3g,
    labels = c(""),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Agricultural soil; genes", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# D)
figure_pcoa4 = do.call("grid_arrange_shared_legend", c(ps_bdm, ncol=7, nrow=2, position="right"))
figure_pcoa4 <- annotate_figure(
  ggarrange(
    figure_pcoa4,
    labels = c("D\n"),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Combined mock communities; contigs", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
figure_pcoa4g = do.call("grid_arrange_shared_legend", c(ps_bdm_genes, ncol=7, nrow=2, position="right"))
figure_pcoa4g <- annotate_figure(
  ggarrange(
    figure_pcoa4g,
    labels = c(""),
    ncol = 1, nrow = 1,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "\n     Agricultural soil; genes", 
  fig.lab.pos = "top.left",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

figure_pcoas <- annotate_figure(
  ggarrange(
    figure_pcoa1, figure_pcoa1g,
    figure_pcoa2, figure_pcoa2g,
    figure_pcoa3, figure_pcoa3g,
    figure_pcoa4, figure_pcoa4g,
    #labels = c(""),
    ncol = 1, nrow = 8,
    widths = c(1),
    common.legend = FALSE, legend = NULL,
    align = "hv", 
    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
  ),
  fig.lab = "Figure S3", 
  fig.lab.pos = "top",
  fig.lab.face = "plain",
  top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)
pdf( file=paste0(outdir, "./figures/F_FIGURE_S3.pdf"), height=26, width=11)
figure_pcoas
dev.off()


###########################################################################################
# Alpha diversity                                                                         #
# For each alpha div table, find the minimal value (baseline) common to all samples       #
# and generate a boxplot for each. This way we can observe if different pipeline methods  #
# affect resulting alpha diversity                                                        #
# Results are being split into 5 quantiles (quintiles)                                    #
###########################################################################################
### China microbiome
mapping = data.frame(fread(mapping_file, header=TRUE, sep="\t"), check.names=FALSE)
df_alpha = get_alphadiv_as_df(my_alpha_div_contigs_files, mapping=mapping)

curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diveristy
df_tmp = df_alpha[df_alpha$Type == "12M",]
df_tmp$quantiles = NULL
df_tmp$quantiles = cut(df_tmp$value, quantile(df_tmp$value, prob = seq(0, 1, length = 6), type = 5, names=TRUE, labels=FALSE, ordered_result=TRUE), include.lowest=TRUE)
head(df_tmp)
unique(df_tmp$quantiles)
df_link = data.frame(quantiles=levels(df_tmp$quantiles), rank=seq(1,(length(unique(df_tmp$quantiles))), by=1))
df_tmp2 = df_tmp[,c("variable", "quantiles")]
df_tmp2 = merge(df_tmp2, df_link, by="quantiles", all=TRUE)
# then relink with df_alpha
df_alpha = merge(df_alpha, df_tmp2[c("variable", "rank")], by="variable")
df_alpha$rank_type2 = paste0("quantile ", df_alpha$rank)
df_alpha$rank_type2 = factor(df_alpha$rank_type2, c("quantile 1",
                                                    "quantile 2",
                                                    "quantile 3",
                                                    "quantile 4",
                                                    "quantile 5"
                                                    ))
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

# split panels by quantiles. So 5 panels in total. In each panels will be all workflows.
# This way differences between workflows will be better put in evidence.
p_alpha1 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value, fill=workflow_type)) + 
  facet_grid(. ~ rank_type2, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.05, pch=19,position=position_jitterdodge(0.3), alpha=0.5) +
  xlab("") + ylab(paste0("Observed contigs")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_blank(),
    legend.spacing = unit(1, "cm"),
    legend.position="top",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) + scale_fill_manual(values=(vColorsA)) + 
  scale_color_manual(values=(vColorsA))
p_alpha1

df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files, mapping=mapping)
curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diveristy
# Keep the same quantiles define for contigs for consistency between plots.

# then relink with df_alpha
df_alpha = merge(df_alpha, df_tmp2[c("variable", "rank")], by="variable")
df_alpha$rank_type2 = paste0("quantile ", df_alpha$rank)
df_alpha$rank_type2 = factor(df_alpha$rank_type2, c("quantile 1",
                                                    "quantile 2",
                                                    "quantile 3",
                                                    "quantile 4",
                                                    "quantile 5"
))
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

p_alpha2 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value)) + 
  facet_grid(. ~ rank_type2, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.05, pch=19,position=position_jitterdodge(0.3), alpha=0.5) +
  xlab("") + ylab(paste0("Observed genes")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA)
p_alpha2

#################################
# Rock microbiome alpha div 
mapping_rock = data.frame(fread("../rock_microbiome/0.1M_clusters/mapping_file.tsv"), check.names=FALSE)
row.names(mapping_rock) = mapping_rock$`#SampleID`
df_alpha = get_alphadiv_as_df(my_alpha_div_contigs_files_rock, mapping=mapping_rock)

# order
curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diveristy
df_tmp = df_alpha[df_alpha$Type == "complete",]
df_tmp$quantiles = NULL
df_tmp$quantiles = cut(df_tmp$value, quantile(df_tmp$value, prob = seq(0, 1, length = 6), type = 5, names=TRUE, labels=FALSE, ordered_result=TRUE))
head(df_tmp)
unique(df_tmp$quantiles)
df_link = data.frame(quantiles=levels(df_tmp$quantiles), rank=seq(1,(length(unique(df_tmp$quantiles))-1), by=1))
df_tmp2 = df_tmp[,c("variable", "quantiles")]
df_tmp2 = merge(df_tmp2, df_link, by="quantiles")
# then relink with df_alpha
df_alpha = merge(df_alpha, df_tmp2[c("variable", "rank")], by="variable")
df_alpha$rank_type2 = paste0("quantile ", df_alpha$rank)
df_alpha$rank_type2 = factor(df_alpha$rank_type2, c("quantile 1",
                                                    "quantile 2",
                                                    "quantile 3",
                                                    "quantile 4",
                                                    "quantile 5"
))
df_alpha$workflow_type = ifelse(grepl("complete", df_alpha$Type), "complete", "")
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

p_alpha3 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value)) + 
  facet_grid(. ~ rank_type2, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.5, pch=19,position=position_jitterdodge(0.3), alpha=1) +
  xlab("") + ylab(paste0("Observed contigs")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
   ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA)
p_alpha3

# genes
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files_rock, mapping=mapping_rock)
curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diveristy
# use the values previously computed

# then relink with df_alpha
df_alpha = merge(df_alpha, df_tmp2[c("variable", "rank")], by="variable")
df_alpha$rank_type2 = paste0("quantile ", df_alpha$rank)
df_alpha$rank_type2 = factor(df_alpha$rank_type2, c("quantile 1",
                                                    "quantile 2",
                                                    "quantile 3",
                                                    "quantile 4",
                                                    "quantile 5"
))
df_alpha$workflow_type = ifelse(grepl("complete", df_alpha$Type), "complete", "")
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

p_alpha4 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value)) + 
  facet_grid(. ~ rank_type2, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.5, pch=19,position=position_jitterdodge(0.3), alpha=1) +
  xlab("") + ylab(paste0("Observed genes")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
    ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA)

p_alpha4

#################################
# Agricultural soil alpha div 
mapping_agr = data.frame(fread("../soil/0.1M_clusters/mapping_file.tsv"), check.names=FALSE)
row.names(mapping_agr) = mapping_agr$`#SampleID`
df_alpha = get_alphadiv_as_df(my_alpha_div_contigs_files_soil, mapping=mapping_agr)

# order
curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diveristy
df_tmp = df_alpha[df_alpha$Type == "complete",]
df_tmp$quantiles = NULL
df_tmp$quantiles = cut(df_tmp$value, quantile(df_tmp$value, prob = seq(0, 1, length = 6), type = 5, names=TRUE, labels=FALSE, ordered_result=TRUE))
head(df_tmp)
unique(df_tmp$quantiles)
df_link = data.frame(quantiles=levels(df_tmp$quantiles), rank=seq(1,(length(unique(df_tmp$quantiles))-1), by=1))
df_tmp2 = df_tmp[,c("variable", "quantiles")]
df_tmp2 = merge(df_tmp2, df_link, by="quantiles")
# then relink with df_alpha
df_alpha = merge(df_alpha, df_tmp2[c("variable", "rank")], by="variable")
df_alpha$rank_type2 = paste0("quantile ", df_alpha$rank)
df_alpha$rank_type2 = factor(df_alpha$rank_type2, c("quantile 1",
                                                    "quantile 2",
                                                    "quantile 3",
                                                    "quantile 4",
                                                    "quantile 5"
))
df_alpha$workflow_type = ifelse(grepl("complete", df_alpha$Type), "complete", "")
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

p_alpha5 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value)) + 
  facet_grid(. ~ rank_type2, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.5, pch=19,position=position_jitterdodge(0.3), alpha=1) +
  xlab("") + ylab(paste0("Observed contigs")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA)
p_alpha5

# genes
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files_soil, mapping=mapping_agr)
curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diveristy
# use the values previously computed

# then relink with df_alpha
df_alpha = merge(df_alpha, df_tmp2[c("variable", "rank")], by="variable")
df_alpha$rank_type2 = paste0("quantile ", df_alpha$rank)
df_alpha$rank_type2 = factor(df_alpha$rank_type2, c("quantile 1",
                                                    "quantile 2",
                                                    "quantile 3",
                                                    "quantile 4",
                                                    "quantile 5"
))
df_alpha$workflow_type = ifelse(grepl("complete", df_alpha$Type), "complete", "")
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

p_alpha6 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value)) + 
  facet_grid(. ~ rank_type2, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.5, pch=19,position=position_jitterdodge(0.3), alpha=1) +
  xlab("") + ylab(paste0("Observed genes")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA)

p_alpha6

#################################
# Mock community alpha div 
mapping_mock = data.frame(fread("../mock/0.1M_clusters/mapping_file.tsv"), check.names=FALSE)
row.names(mapping_mock) = mapping_mock$`#SampleID`
df_alpha = get_alphadiv_as_df(my_alpha_div_contigs_files_mock, mapping=mapping_mock)

# order
curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

df_alpha$workflow_type = ifelse(grepl("complete", df_alpha$Type), "complete", "")
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

p_alpha7 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value)) + 
  facet_grid(. ~ Treatment, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.5, pch=19,position=position_jitterdodge(0.3), alpha=1) +
  xlab("") + ylab(paste0("Observed contigs")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA)
p_alpha7

# genes
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files_mock, mapping=mapping_mock)
curr_order_facet = as.character(unique(df_alpha$Type))
df_alpha$Type = factor(df_alpha$Type, levels=curr_order_facet)
df_alpha$dummy = "dummy"
df_alpha$num_clusters = df_alpha$Type
df_alpha$num_clusters = gsub("_all", "", df_alpha$Type)
df_alpha$num_clusters = factor(df_alpha$num_clusters, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete"))
df_alpha$Type = gsub("_all", " - arm", df_alpha$Type)

df_alpha$workflow_type = ifelse(grepl("complete", df_alpha$Type), "complete", "")
df_alpha$workflow_type = ifelse(grepl("arm", df_alpha$Type), "arm", "standard")

p_alpha8 <- ggplot(data=df_alpha, aes(x=num_clusters, y=value)) + 
  facet_grid(. ~ Treatment, scales="free_x", space="free_x") + 
  geom_boxplot(aes(fill=workflow_type), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow_type),size=0.5, pch=19,position=position_jitterdodge(0.3), alpha=1) +
  xlab("") + ylab(paste0("Observed genes")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=(11), face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="no",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=11, face="plain"),
    strip.text.y = element_blank(),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColorsA)) + scale_color_manual(values=vColorsA)

p_alpha8

# assemble figure
figure5a <- ggarrange(
  p_alpha1 + rremove("xlab") + rremove("x.text"), p_alpha2,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = TRUE, legend = NULL,
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f5a = annotate_figure(figure5a, 
                     fig.lab = "\n       Human gut", 
                     fig.lab.pos = "top.left",
                     top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

figure5b <- ggarrange(
  p_alpha3 + rremove("xlab") + rremove("x.text"), p_alpha4,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f5b = annotate_figure(figure5b, 
                     fig.lab = "\n     Antarctic", 
                     fig.lab.pos = "top.left",
                     top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

figure5c <- ggarrange(
  p_alpha5 + rremove("xlab") + rremove("x.text"), p_alpha6,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f5c = annotate_figure(figure5c, 
                      fig.lab = "\n     Agricultural soil", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

figure5d <- ggarrange(
  p_alpha7 + rremove("xlab") + rremove("x.text"), p_alpha8,
  labels = NULL,
  ncol = 1, nrow = 2,
  #widths = c(0.33,0.33,0.33),
  common.legend = FALSE, legend = NULL,
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f5d = annotate_figure(
  ggarrange(
    figure5d, nrow=1,ncol=3,widths= c(0.33,0.33,0.33)
  ),
            fig.lab = "\n     Mock communities",  c(0.33,0.33,0.33),
            fig.lab.pos = "top.left",
            top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
  
)

final_list = annotate_figure(
  ggarrange(f5a,f5b,f5c,f5d,
            labels=c("A", "B", "C", "D"),
            nrow = 4,
            heights=c(0.5, 0.5,0.5,0.5)
  ),
  top = text_grob("Figure 3\n", color = "black", face = "plain", size = 11),
)

pdf( file=paste0(outdir, "./figures/F_FIGURE_3.pdf"), height=16, width=9)
final_list
dev.off()

###########################################################################
# TAXONOMY China human gut                                                #
# Plot taxonomy at say class level or order level                         #
# for all samples. compact/squish/squeeze them and hide sample id         #
# so we can get nice summarized profiles to compare.                      #
# Figure S4                                                               #
###########################################################################
df_tax_counts = get_tax_as_df_unmelted(my_tax_files_china_counts, mapping_file="./mapping_file.tsv")
df_tax =        get_tax_as_df_unmelted(my_tax_files_china_counts, convert_to_rel=TRUE, mapping_file="./mapping_file.tsv")

mapping = data.frame(fread("./mapping_file.tsv"), check.names=FALSE)

taxa = unique(df_tax$Taxon)
taxa = taxa[taxa != "NULL"]

types = unique(df_tax$Type)

x = 1
final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_tax[df_tax$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$Taxon
    df_x$Taxon = NULL
    df_y = df_tax[df_tax$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$Taxon
    df_y$Taxon = NULL
    
    curr_taxa_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_taxa_x = unique(row.names(df_x))
    curr_taxa_y = unique(row.names(df_y))
    
    missing_in_x = curr_taxa_all[!curr_taxa_all %in% curr_taxa_x]
    missing_in_y = curr_taxa_all[!curr_taxa_all %in% curr_taxa_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)

# It is time intensive to generate the spearman correlations. Write the file as backup.
#write.table(df3, paste0(root, "/melted_style_spearman_taxonomy_china_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_taxonomy_china_df.tsv"), sep="\t"), check.names=FALSE)

df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

p_tax_china = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) + 
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    #axis.title=element_text(family="Helvetica", size=14),
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_tax_china)

##########################
# Venn diagram
##########################
taxa = unique(df_tax$Taxon)

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_tax$Type)

for(i in 1:length(types)){
  df_x = df_tax[df_tax$Type == types[i],]
  curr_taxa_x = unique(df_x$Taxon)
  
  missing_in_x = taxa[!taxa %in% curr_taxa_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_taxa = unique(df_x$Taxon)
  unique_taxa = unique_taxa[!unique_taxa %in% c("NULL", "0")]
  unique_taxa = data.frame(Taxon=unique_taxa, Type=types[i])
  freq_taxa = df_x[, c("Taxon", "sum")]
  freq_taxa = freq_taxa[!freq_taxa$Taxon %in% c("NULL", "0"),]
  freq_taxa$Type = types[i]
  
  if(x == 1){
    df_freq = freq_taxa
    df_unique = unique_taxa
  }else{
    df_freq = rbind(df_freq, freq_taxa)
    df_unique = rbind(df_unique, unique_taxa)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>% 
  mutate(Number_of_unique_taxa_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))

p_tax_china_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_taxa_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Workflow configuration") + ylab(str_wrap("Number of unique taxa recovered (genus level)", width = 20, indent = 0, exdent = 0)) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_tax_china_2

my_order = df_freq[df_freq$Type == "12M",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$Taxon
df_freq$Taxon = factor(df_freq$Taxon, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "4M", "4M_all", "8M", "8M_all", "12M", "12M_all"))

tax_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "4M", "8M", "12M"),]
tax_rest = tax_rest[tax_rest$sum > 0,]
tax_rest_12M = tax_rest[tax_rest$Type %in% c("12M"),]
tax_rest2 = tax_rest[tax_rest$Type != "12M",]
tax_diff = tax_rest_12M[!tax_rest_12M$Taxon %in% tax_rest2$Taxon,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$Taxon)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$Taxon)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$Taxon)))
M1 = M1[!is.na(M1)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$Taxon)))
M4 = M4[!is.na(M4)]
M8 = unique(as.character((df_freq[df_freq$Type == "8M",]$Taxon)))
M8 = M8[!is.na(M8)]
M12 = unique(as.character((df_freq[df_freq$Type == "12M",]$Taxon)))
M12 = M12[!is.na(M12)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "4M"=M4, "8M"=M8, "12M"=M12)
p_china_tax_venn = ggVennDiagram(x, label_size=3, label="count") +
  scale_fill_gradientn(name = "Taxonomic lineage (genus) count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  theme(
    text=element_text(size=12,family="Helvetica")
  )

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 368,]$item
groupB = venn_data2[venn_data2$count == 473,]$item
groupC = venn_data2[venn_data2$count == 342,]$item
groupD = venn_data2[venn_data2$count == 360,]$item
groupE = venn_data2[venn_data2$count == 526,]$item
groupF = venn_data2[venn_data2$count == 631,]$item

abundance = data.frame(fread("./export/12M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"), check.names=FALSE)

groups = c(groupA,groupB,groupC,groupD,groupE,groupF)

df = NULL
group_names = c("A - 368","B - 473","C - 342","D - 360","E - 526","F - 631")
for(i in 1:length(groups)){
  group = groups[[i]]
  abundance2 = abundance[abundance$Taxon %in% group,]
  tmp = data.frame(colSums(abundance2[,2:913]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the 12M reads configuration", width=30), size=8)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(8), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=8)
tab2_tax_china = tab2


###########################################
#  Rock/Antarctic microbiome TAXONOMY     #
###########################################
df_tax_counts = get_tax_as_df_unmelted(my_tax_files_rock_counts, mapping_file="../rock_microbiome/export/mapping_file.tsv")
df_tax =        get_tax_as_df_unmelted(my_tax_files_rock_counts, convert_to_rel=TRUE, mapping_file="../rock_microbiome/export/mapping_file.tsv")

mapping = data.frame(fread("../rock_microbiome/export/mapping_file.tsv"), check.names=FALSE)

taxa = unique(df_tax$Taxon)
taxa = taxa[taxa != "NULL"]

types = unique(df_tax$Type)

x = 1
final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_tax[df_tax$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$Taxon
    df_x$Taxon = NULL
    df_y = df_tax[df_tax$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$Taxon
    df_y$Taxon = NULL
    
    curr_taxa_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_taxa_x = unique(row.names(df_x))
    curr_taxa_y = unique(row.names(df_y))
    
    missing_in_x = curr_taxa_all[!curr_taxa_all %in% curr_taxa_x]
    missing_in_y = curr_taxa_all[!curr_taxa_all %in% curr_taxa_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
   
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)

# Write file as checkpoint.
#write.table(df3, paste0(root, "/melted_style_spearman_taxonomy_rock_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_taxonomy_rock_df.tsv"), sep="\t"), check.names=FALSE)
#
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm","complete"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm","complete"))
df3[58,]$x = "0.1M - arm"
df3[59,]$x = "0.5M - arm"
df3[60,]$x = "1M - arm"
df3[61,]$x = "4M - arm"
df3[62,]$x = "8M - arm"
df3[63,]$x = "12M - arm"

df3[58,]$y = "complete"
df3[59,]$y = "complete"
df3[60,]$y = "complete"
df3[61,]$y = "complete"
df3[62,]$y = "complete"
df3[63,]$y = "complete"

p_tax_rock = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) + 
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_tax_rock)

##########################
# Venn diagram rock
##########################
taxa = unique(df_tax$Taxon)

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_tax$Type)

for(i in 1:length(types)){
  df_x = df_tax[df_tax$Type == types[i],]
  curr_taxa_x = unique(df_x$Taxon)
 
  missing_in_x = taxa[!taxa %in% curr_taxa_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_taxa = unique(df_x$Taxon)
  unique_taxa = unique_taxa[!unique_taxa %in% c("NULL", "0")]
  unique_taxa = data.frame(Taxon=unique_taxa, Type=types[i])
  freq_taxa = df_x[, c("Taxon", "sum")]
  freq_taxa = freq_taxa[!freq_taxa$Taxon %in% c("NULL", "0"),]
  freq_taxa$Type = types[i]
  
  if(x == 1){
    df_freq = freq_taxa
    df_unique = unique_taxa
  }else{
    df_freq = rbind(df_freq, freq_taxa)
    df_unique = rbind(df_unique, unique_taxa)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>%
  mutate(Number_of_unique_taxa_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))

p_tax_rock_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_taxa_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Workflow configuration") + ylab(str_wrap("Number of unique taxa recovered (genus level)", width = 20, indent = 0, exdent = 0)) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_tax_rock_2

my_order = df_freq[df_freq$Type == "complete",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$Taxon
df_freq$Taxon = factor(df_freq$Taxon, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "4M", "4M_all", "8M", "8M_all", "12M", "12M_all", "complete"))

tax_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"),]
tax_rest = tax_rest[tax_rest$sum > 0,]
tax_rest_12M = tax_rest[tax_rest$Type %in% c("complete"),]
tax_rest2 = tax_rest[tax_rest$Type != "complete",]
tax_diff = tax_rest_12M[!tax_rest_12M$Taxon %in% tax_rest2$Taxon,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$Taxon)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$Taxon)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$Taxon)))
M1 = M1[!is.na(M1)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$Taxon)))
M4 = M4[!is.na(M4)]
M8 = unique(as.character((df_freq[df_freq$Type == "8M",]$Taxon)))
M8 = M8[!is.na(M8)]
M12 = unique(as.character((df_freq[df_freq$Type == "12M",]$Taxon)))
M12 = M12[!is.na(M12)]
Mcomplete = unique(as.character((df_freq[df_freq$Type == "complete",]$Taxon)))
Mcomplete = Mcomplete[!is.na(Mcomplete)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "4M"=M4, "8M"=M8, "12M"=M12, "complete"=Mcomplete)
p_rock_tax_venn = ggVennDiagram(x, label_size=3, label="count") +
  scale_fill_gradientn(name = "Taxonomic lineage (genus) count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  theme(
    text=element_text(size=12,family="Helvetica")
)

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 549,]$item
groupB = venn_data2[venn_data2$count == 234,]$item
groupC = venn_data2[venn_data2$count == 176,]$item
groupD = venn_data2[venn_data2$count == 571,]$item
groupE = venn_data2[venn_data2$count == 459,]$item
groupF = venn_data2[venn_data2$count == 933,]$item
groupG = venn_data2[venn_data2$count == 386,]$item

abundance = data.frame(fread("../rock_microbiome/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"), check.names=FALSE)

groups = c(groupA,groupB,groupC,groupD,groupE,groupF, groupG)

df = NULL
group_names = c("A - 549","B - 234","C - 176","D - 571","E - 459","F - 933", "G - 386")
for(i in 1:length(groups)){
  group = groups[[i]]
  abundance2 = abundance[abundance$Taxon %in% group,]
  tmp = data.frame(colSums(abundance2[,2:19]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the complete dataset configuration", width=30), size=8)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(8), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=8)
tab2_tax_rock = tab2


###########################################
#   Agricultural soil.                    #
###########################################
df_tax_counts = get_tax_as_df_unmelted(my_tax_files_soil_counts, mapping_file="../soil/export/mapping_file.tsv")
df_tax =        get_tax_as_df_unmelted(my_tax_files_soil_counts, convert_to_rel=TRUE, mapping_file="../soil/export/mapping_file.tsv")

mapping = data.frame(fread("../soil/export/mapping_file.tsv"), check.names=FALSE)

taxa = unique(df_tax$Taxon)
taxa = taxa[taxa != "NULL"]

types = unique(df_tax$Type)

x = 1
final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_tax[df_tax$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$Taxon
    df_x$Taxon = NULL
    df_y = df_tax[df_tax$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$Taxon
    df_y$Taxon = NULL
    
    curr_taxa_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_taxa_x = unique(row.names(df_x))
    curr_taxa_y = unique(row.names(df_y))
    
    missing_in_x = curr_taxa_all[!curr_taxa_all %in% curr_taxa_x]
    missing_in_y = curr_taxa_all[!curr_taxa_all %in% curr_taxa_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)

# Write file as checkpoint.
#write.table(df3, paste0(root, "/melted_style_spearman_taxonomy_rock_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_taxonomy_rock_df.tsv"), sep="\t"), check.names=FALSE)
#
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm","complete"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm","complete"))
df3[58,]$x = "0.1M - arm"
df3[59,]$x = "0.5M - arm"
df3[60,]$x = "1M - arm"
df3[61,]$x = "4M - arm"
df3[62,]$x = "8M - arm"
df3[63,]$x = "12M - arm"

df3[58,]$y = "complete"
df3[59,]$y = "complete"
df3[60,]$y = "complete"
df3[61,]$y = "complete"
df3[62,]$y = "complete"
df3[63,]$y = "complete"

p_tax_soil = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) + 
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_tax_soil)

##########################
# Venn diagram soil
##########################
taxa = unique(df_tax$Taxon)

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_tax$Type)

for(i in 1:length(types)){
  df_x = df_tax[df_tax$Type == types[i],]
  curr_taxa_x = unique(df_x$Taxon)
  
  missing_in_x = taxa[!taxa %in% curr_taxa_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_taxa = unique(df_x$Taxon)
  unique_taxa = unique_taxa[!unique_taxa %in% c("NULL", "0")]
  unique_taxa = data.frame(Taxon=unique_taxa, Type=types[i])
  freq_taxa = df_x[, c("Taxon", "sum")]
  freq_taxa = freq_taxa[!freq_taxa$Taxon %in% c("NULL", "0"),]
  freq_taxa$Type = types[i]
  
  if(x == 1){
    df_freq = freq_taxa
    df_unique = unique_taxa
  }else{
    df_freq = rbind(df_freq, freq_taxa)
    df_unique = rbind(df_unique, unique_taxa)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>%
  mutate(Number_of_unique_taxa_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))

p_tax_soil_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_taxa_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Workflow configuration") + ylab(str_wrap("Number of unique taxa recovered (genus level)", width = 20, indent = 0, exdent = 0)) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_tax_soil_2

my_order = df_freq[df_freq$Type == "complete",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$Taxon
df_freq$Taxon = factor(df_freq$Taxon, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "4M", "4M_all", "8M", "8M_all", "12M", "12M_all", "complete"))

tax_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"),]
tax_rest = tax_rest[tax_rest$sum > 0,]
tax_rest_12M = tax_rest[tax_rest$Type %in% c("complete"),]
tax_rest2 = tax_rest[tax_rest$Type != "complete",]
tax_diff = tax_rest_12M[!tax_rest_12M$Taxon %in% tax_rest2$Taxon,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$Taxon)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$Taxon)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$Taxon)))
M1 = M1[!is.na(M1)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$Taxon)))
M4 = M4[!is.na(M4)]
M8 = unique(as.character((df_freq[df_freq$Type == "8M",]$Taxon)))
M8 = M8[!is.na(M8)]
M12 = unique(as.character((df_freq[df_freq$Type == "12M",]$Taxon)))
M12 = M12[!is.na(M12)]
Mcomplete = unique(as.character((df_freq[df_freq$Type == "complete",]$Taxon)))
Mcomplete = Mcomplete[!is.na(Mcomplete)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "4M"=M4, "8M"=M8, "12M"=M12, "complete"=Mcomplete)
p_soil_tax_venn = ggVennDiagram(x, label_size=3, label="count") +
  scale_fill_gradientn(name = "Taxonomic lineage (genus) count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  theme(
    text=element_text(size=12,family="Helvetica")
  )

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 322,]$item
groupB = venn_data2[venn_data2$count == 86,]$item
groupC = venn_data2[venn_data2$count == 707,]$item
groupD = venn_data2[venn_data2$count == 420,]$item
groupE = venn_data2[venn_data2$count == 486,]$item
groupF = venn_data2[venn_data2$count == 1246,]$item
groupG = venn_data2[venn_data2$count == 400,]$item
groupH = venn_data2[venn_data2$count == 25,]$item

abundance = data.frame(fread("../soil/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"), check.names=FALSE)

groups = c(groupA,groupB,groupC,groupD,groupE,groupF,groupG,groupH)

df = NULL
group_names = c("A - 322","B - 86","C - 707","D - 420", "E - 486","F - 1246", "G - 400", "H - 25")
for(i in 1:length(groups)){
  group = groups[[i]]
  abundance2 = abundance[abundance$Taxon %in% group,]
  tmp = data.frame(colSums(abundance2[,2:73]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the complete dataset configuration", width=30), size=8)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F,G,H letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(10), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=10)
tab2_tax_soil = tab2


################################
# Mock communities             #
################################
df_tax_counts = get_tax_as_df_unmelted(my_tax_files_mock_counts, mapping_file="../mock/export/mapping_file.tsv")
df_tax =        get_tax_as_df_unmelted(my_tax_files_mock_counts, convert_to_rel=TRUE, mapping_file="../mock/export/mapping_file.tsv")

mapping = data.frame(fread("../mock/export/mapping_file.tsv"), check.names=FALSE)

taxa = unique(df_tax$Taxon)
taxa = taxa[taxa != "NULL"]

types = unique(df_tax$Type)

x = 1
final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_tax[df_tax$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$Taxon
    df_x$Taxon = NULL
    df_y = df_tax[df_tax$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$Taxon
    df_y$Taxon = NULL
    
    curr_taxa_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_taxa_x = unique(row.names(df_x))
    curr_taxa_y = unique(row.names(df_y))
    
    missing_in_x = curr_taxa_all[!curr_taxa_all %in% curr_taxa_x]
    missing_in_y = curr_taxa_all[!curr_taxa_all %in% curr_taxa_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M", y="2M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M", y="3M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M_all", y="2M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M_all", y="3M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)
df3$r = round(df3$r,3)
df3$r[is.na(df3$r)]<-0

# Write file as checkpoint.
#write.table(df3, paste0(root, "/melted_style_spearman_taxonomy_rock_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_taxonomy_rock_df.tsv"), sep="\t"), check.names=FALSE)
#
df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm","complete"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm","complete"))
df3[58,]$x = "0.1M - arm"
df3[59,]$x = "0.5M - arm"
df3[60,]$x = "1M - arm"
df3[61,]$x = "2M - arm"
df3[62,]$x = "3M - arm"
df3[63,]$x = "4M - arm"

df3[58,]$y = "complete"
df3[59,]$y = "complete"
df3[60,]$y = "complete"
df3[61,]$y = "complete"
df3[62,]$y = "complete"
df3[63,]$y = "complete"

p_tax_mock = ggplot(df3, aes(x=y, y=x,label=r)) +
  geom_tile(aes(fill=r)) + 
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_tax_mock)

##########################
# Venn diagram rock
##########################
taxa = unique(df_tax$Taxon)

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_tax$Type)

for(i in 1:length(types)){
  df_x = df_tax[df_tax$Type == types[i],]
  curr_taxa_x = unique(df_x$Taxon)
  
  missing_in_x = taxa[!taxa %in% curr_taxa_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_taxa = unique(df_x$Taxon)
  unique_taxa = unique_taxa[!unique_taxa %in% c("NULL", "0")]
  unique_taxa = data.frame(Taxon=unique_taxa, Type=types[i])
  freq_taxa = df_x[, c("Taxon", "sum")]
  freq_taxa = freq_taxa[!freq_taxa$Taxon %in% c("NULL", "0"),]
  freq_taxa$Type = types[i]
  
  if(x == 1){
    df_freq = freq_taxa
    df_unique = unique_taxa
  }else{
    df_freq = rbind(df_freq, freq_taxa)
    df_unique = rbind(df_unique, unique_taxa)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>%
  mutate(Number_of_unique_taxa_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete"))

p_tax_mock_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_taxa_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Workflow configuration") + ylab(str_wrap("Number of unique taxa recovered (genus level)", width = 20, indent = 0, exdent = 0)) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_tax_mock_2

my_order = df_freq[df_freq$Type == "complete",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$Taxon
df_freq$Taxon = factor(df_freq$Taxon, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "2M", "2M_all", "3M", "3M_all", "4M", "4M_all", "complete"))

tax_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete"),]
tax_rest = tax_rest[tax_rest$sum > 0,]
tax_rest_12M = tax_rest[tax_rest$Type %in% c("complete"),]
tax_rest2 = tax_rest[tax_rest$Type != "complete",]
tax_diff = tax_rest_12M[!tax_rest_12M$Taxon %in% tax_rest2$Taxon,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$Taxon)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$Taxon)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$Taxon)))
M1 = M1[!is.na(M1)]
M2 = unique(as.character((df_freq[df_freq$Type == "2M",]$Taxon)))
M2 = M2[!is.na(M2)]
M3 = unique(as.character((df_freq[df_freq$Type == "3M",]$Taxon)))
M3 = M3[!is.na(M3)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$Taxon)))
M4 = M4[!is.na(M4)]
Mcomplete = unique(as.character((df_freq[df_freq$Type == "complete",]$Taxon)))
Mcomplete = Mcomplete[!is.na(Mcomplete)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "2M"=M2, "3M"=M3, "4M"=M4, "complete"=Mcomplete)
p_mock_tax_venn = ggVennDiagram(x, label_size=3, label="count") +
  scale_fill_gradientn(name = "Taxonomic lineage (genus) count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  theme(
    text=element_text(size=12,family="Helvetica")
  )

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 31,]$item
groupB = venn_data2[venn_data2$name == "2M..3M..4M..complete",]$item
groupC = venn_data2[venn_data2$name == "0.5M..1M..2M..3M..complete",]$item
groupD = venn_data2[venn_data2$name == "0.1M..0.5M..1M..2M..3M..4M..complete",]$item

abundance = data.frame(fread("../mock/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"), check.names=FALSE)

groups = c(groupA,groupB,groupC,groupD)

df = NULL
group_names = c("A - 31","B - 1","C - 1","D - 1")
for(i in 1:length(groups)){
  group = groups[[i]]
  abundance2 = abundance[abundance$Taxon %in% group,]
  tmp = data.frame(colSums(abundance2[,2:7]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the complete dataset configuration", width=30), size=8)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(6), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=6)
tab2_tax_mock = tab2


# Assemble figure together

# plots china
figure6a1 <- ggarrange(
  p_tax_china_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1,0.9),
  widths = c(1, 1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6a2 <- ggarrange(
  p_tax_china, figure6a1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6a3 = ggarrange(
  p_china_tax_venn, tab2_tax_china,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6a4 = ggarrange(
  figure6a2, figure6a3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f6a = annotate_figure(figure6a4, 
                      fig.lab = "\n       Human gut", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# plots Rock
figure6b1 <- ggarrange(
  p_tax_rock_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1,0.9),
  widths = c(1,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6b2 <- ggarrange(
  p_tax_rock, figure6b1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6b3 = ggarrange(
  p_rock_tax_venn, tab2_tax_rock,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6b4 = ggarrange(
  figure6b2, figure6b3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f6b = annotate_figure(figure6b4, 
                      fig.lab = "\n       Antarctic", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# plots soil
figure6c1 <- ggarrange(
  p_tax_soil_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1,0.9),
  widths = c(1,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6c2 <- ggarrange(
  p_tax_soil, figure6c1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6c3 = ggarrange(
  p_soil_tax_venn, tab2_tax_soil,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6c4 = ggarrange(
  figure6c2, figure6c3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f6c = annotate_figure(figure6c4, 
                      fig.lab = "\n       Agricultural soil", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

#plots mock
figure6d1 <- ggarrange(
  p_tax_mock_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1,0.9),
  widths = c(1,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6d2 <- ggarrange(
  p_tax_mock, figure6d1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6d3 = ggarrange(
  p_mock_tax_venn, tab2_tax_mock,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure6d4 = ggarrange(
  figure6d2, figure6d3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f6d = annotate_figure(figure6d4, 
                      fig.lab = "\n       Mock communities", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)


# Assemble final figure.
final_list = annotate_figure(
  ggarrange(f6a,f6b,f6c,f6d,
            labels=c("A", "B", "C", "D"),
            nrow = 4,
            heights=c(1.2, 1.2,1.2,1.2)
  ),
  top = text_grob("Figure S4\n", color = "black", face = "plain", size = 11),
)

pdf( file=paste0(outdir, "./figures/FIGURE_S4.pdf"), height=28, width=12)
final_list
dev.off()


######################################
## MAGs quality                      #
## Summarize MAGs                    #
## Figure 6                          #
######################################
df_mags_china = get_mags_as_df(my_mags_files_china, contam=100, complete=0, print_high_qual=TRUE)
df_mags_china$status = "very_low"
df_mags_china[df_mags_china$Completeness >=  50 & df_mags_china$Contamination <= 50,]$status = "low"
df_mags_china[df_mags_china$Completeness >= 75 & df_mags_china$Contamination <= 25,]$status = "average"
df_mags_china[df_mags_china$Completeness >= 90 & df_mags_china$Contamination <= 10,]$status = "good"
df_mags_china[df_mags_china$Completeness >= 95 & df_mags_china$Contamination <= 5,]$status = "very_good"
df_mags_china[df_mags_china$Completeness >= 98 & df_mags_china$Contamination <= 2,]$status = "excellent"

df_mags_china_melted = melt(df_mags_china)
df_mags_china_melted = df_mags_china_melted[df_mags_china_melted$variable %in% c("Completeness", "Contamination"),]
df_mags_china_melted$number_of_clusters = gsub("_all", "", df_mags_china_melted$Type)
#df$var <- ifelse(df$var == " ?", " Private", as.character(df$var))   
df_mags_china_melted$Type2 = ifelse(grepl("^.*_all", df_mags_china_melted$Type), "arm", "standard")
df_mags_china_melted$number_of_clusters = factor(df_mags_china_melted$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M" , "12M"))

p_mags1 <- ggplot(data=df_mags_china_melted[df_mags_china_melted$variable == "Completeness",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + #, space="free") +
  xlab("Completeness (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) + scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" ))
p_mags1

p_mags2 <- ggplot(data=df_mags_china_melted[df_mags_china_melted$variable == "Contamination",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + #, space="free") +
  xlab("Contamination (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )+ scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" ))  + xlim(c(0, NA))
p_mags2

##################
# MAGs Rock      #
##################
df_mags_rock = get_mags_as_df(my_mags_files_rock, contam=100, complete=0, print_high_qual=TRUE)

df_mags_rock_melted = melt(df_mags_rock)
df_mags_rock_melted = df_mags_rock_melted[df_mags_rock_melted$variable %in% c("Completeness", "Contamination"),]
df_mags_rock_melted$number_of_clusters = gsub("_all", "", df_mags_rock_melted$Type)
df_mags_rock_melted$Type2 = ifelse(grepl("^.*_all", df_mags_rock_melted$Type), "arm", "standard")
df_mags_rock_melted$number_of_clusters = factor(df_mags_rock_melted$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M" , "12M", "complete"))

p_mags3 <- ggplot(data=df_mags_rock_melted[df_mags_rock_melted$variable == "Completeness",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + 
  xlab("Completeness (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) + scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" ))
p_mags3

p_mags4 <- ggplot(data=df_mags_rock_melted[df_mags_rock_melted$variable == "Contamination",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + #, space="free") +
  xlab("Contamination (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) +  scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" ))  + xlim(c(0, NA))
p_mags4

##################
# MAGs agri soil #
##################
df_mags_soil = get_mags_as_df(my_mags_files_soil, contam=100, complete=0, print_high_qual=TRUE)

df_mags_soil_melted = melt(df_mags_soil)
df_mags_soil_melted = df_mags_soil_melted[df_mags_soil_melted$variable %in% c("Completeness", "Contamination"),]
df_mags_soil_melted$number_of_clusters = gsub("_all", "", df_mags_soil_melted$Type)
df_mags_soil_melted$Type2 = ifelse(grepl("^.*_all", df_mags_soil_melted$Type), "arm", "standard")
df_mags_soil_melted$number_of_clusters = factor(df_mags_soil_melted$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M" , "12M", "complete"))

p_mags5 <- ggplot(data=df_mags_soil_melted[df_mags_soil_melted$variable == "Completeness",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + 
  xlab("Completeness (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) + scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" ))
p_mags5

p_mags6 <- ggplot(data=df_mags_soil_melted[df_mags_soil_melted$variable == "Contamination",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + #, space="free") +
  xlab("Contamination (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) +  scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" ))  + xlim(c(0, NA))
p_mags6

#######################
# MAGs mock community #
#######################
df_mags_mock = get_mags_as_df(my_mags_files_mock, contam=100, complete=0, print_high_qual=TRUE)

df_mags_mock_melted = melt(df_mags_mock)
df_mags_mock_melted = df_mags_mock_melted[df_mags_mock_melted$variable %in% c("Completeness", "Contamination"),]
df_mags_mock_melted$number_of_clusters = gsub("_all", "", df_mags_mock_melted$Type)
df_mags_mock_melted$Type2 = ifelse(grepl("^.*_all", df_mags_mock_melted$Type), "arm", "standard")
df_mags_mock_melted$number_of_clusters = factor(df_mags_mock_melted$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "2M", "3M" , "4M", "complete"))

p_mags7 <- ggplot(data=df_mags_mock_melted[df_mags_mock_melted$variable == "Completeness",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + 
  xlab("Completeness (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) + scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" ))
p_mags7

p_mags8 <- ggplot(data=df_mags_mock_melted[df_mags_mock_melted$variable == "Contamination",], aes(x=value, color=Type2, fill=Type2)) + 
  geom_histogram(binwidth=2) + 
  facet_wrap(. ~ number_of_clusters, scales="free_x", nrow=1) + #, space="free") +
  xlab("Contamination (%)") + 
  ylab("Count") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=6.5, colour="black"), 
    axis.title=element_text(family="Helvetica", size=10),
    plot.title = element_text(lineheight=1.2, face="bold", size=16),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size=10, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.05, "cm") ,
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  ) +  scale_fill_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) +
  scale_color_manual(values=c("arm" = "#0000CD", "standard" = "#FF0000" )) + xlim(c(0, NA))
p_mags8

# assemble figure
figure8a <- ggarrange(
  p_mags1 + rremove("legend"), p_mags2 + rremove("legend"),
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend=NULL,
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f8a = annotate_figure(figure8a,
                     fig.lab = "   Human gut",
                     fig.lab.pos = "top.left",
                     top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.2)
)

#f8a
figure8b <- ggarrange(
  p_mags3, p_mags4,
  ncol = 1, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f8b = annotate_figure(figure8b, 
                     fig.lab = "    Antarctic",
                     fig.lab.pos = "top.left",
                     top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.5)
)

#f8b
figure8c <- ggarrange(
  p_mags5, p_mags6,
  ncol = 1, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f8c = annotate_figure(figure8c, 
                      fig.lab = "    Agricultural soil",
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.5)
)

#f8c
figure8d <- ggarrange(
  p_mags7, p_mags8,
  ncol = 1, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f8d = annotate_figure(figure8d, 
                      fig.lab = "    Mock communities",
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.5)
)

#f8d
final_list = annotate_figure(
  ggarrange(f8a,f8b,f8c,f8d,
            labels=c("A", "B", "C", "D"),
            nrow = 4,
            heights=c(0.46, 0.52, 0.52, 0.52)
  ),
  top = text_grob("Figure 6\n", color = "black", face = "plain", size = 11),
)
#final_list

# Labels of the scatter plot
pdf( file=paste0(outdir, "./figures/F_FIGURE_6.pdf"), height=14, width=6.3)
final_list
dev.off()

#15 pixel (X) = 15 × 0.0104166667 in = 0.15625 in
800 *  0.0104166667
559 *  0.0104166667


#################################
# KOs                           #
# Figure S5                     #
#################################

############## CHINA KO ######################
# Normalize with recA (from COG annotations) #
##############################################
#df_KO_all = get_KO_abundance_as_df(my_KO_files_china, my_gene_abundance_files_china, my_COG_files_china, mapping_file="./mapping_file.tsv")
# write file because long to load...
#fwrite(df_KO_all, paste0(outdir, "./df_KO_china.tsv"), sep="\t", row.names=TRUE)
df_KO_all = data.frame(fread(paste0(outdir, "./df_KO_china.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO_all) = df_KO_all$V1; df_KO_all$V1=NULL;
head(df_KO_all)

mapping = data.frame(fread("./mapping_file.tsv"), check.names=FALSE)

KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

types = unique(df_KO_all$Type)

x = 1
final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_KO_all[df_KO_all$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$KO
    df_x$KO = NULL
    df_y = df_KO_all[df_KO_all$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$KO
    df_y$KO = NULL
    
    curr_KOs_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_KOs_all = curr_KOs_all[curr_KOs_all != "NULL"]
    curr_KOs_x = unique(row.names(df_x))
    curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
    curr_KOs_y = unique(row.names(df_y))
    curr_KOs_y = curr_KOs_y[curr_KOs_y != "NULL"]
    
    missing_in_x = curr_KOs_all[!curr_KOs_all %in% curr_KOs_x]
    missing_in_y = curr_KOs_all[!curr_KOs_all %in% curr_KOs_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)

#write.table(df3, paste0(root, "/melted_style_spearman_KO_china_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
df3 = data.frame(fread(paste0(root, "/melted_style_spearman_KO_china_df.tsv"), sep="\t"), check.names=FALSE)

df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm"))

p_KO_china = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) + 
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_KO_china)

##########################
# Investigate rare KOs
##########################
KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_KO_all$Type)

for(i in 1:length(types)){
  df_x = df_KO_all[df_KO_all$Type == types[i],]
  curr_KOs_x = unique(df_x$KO)
  curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
  
  missing_in_x = KOs[!KOs %in% curr_KOs_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_KOs = unique(df_x$KO)
  unique_KOs = unique_KOs[!unique_KOs %in% c("NULL", "0")]
  unique_KOs = data.frame(KO=unique_KOs, Type=types[i])
  freq_KOs = df_x[, c("KO", "sum")]
  freq_KOs = freq_KOs[!freq_KOs$KO %in% c("NULL", "0"),]
  freq_KOs$Type = types[i]
  
  if(x == 1){
    df_freq = freq_KOs
    df_unique = unique_KOs
  }else{
    df_freq = rbind(df_freq, freq_KOs)
    df_unique = rbind(df_unique, unique_KOs)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>% 
  mutate(Number_of_unique_KOs_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))

p_KO_china_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_KOs_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Analysis configuration") + ylab(paste0("Number of unique KOs recovered")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_KO_china_2

my_order = df_freq[df_freq$Type == "12M",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$KO
df_freq$KO = factor(df_freq$KO, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "4M", "4M_all", "8M", "8M_all", "12M", "12M_all"))

KO_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "4M", "8M", "12M"),]
KO_rest = KO_rest[KO_rest$sum > 0,]
KO_rest_12M = KO_rest[KO_rest$Type %in% c("12M"),]
KO_rest2 = KO_rest[KO_rest$Type != "12M",]
KO_diff = KO_rest_12M[!KO_rest_12M$KO %in% KO_rest2$KO,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$KO)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$KO)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$KO)))
M1 = M1[!is.na(M1)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$KO)))
M4 = M4[!is.na(M4)]
M8 = unique(as.character((df_freq[df_freq$Type == "8M",]$KO)))
M8 = M8[!is.na(M8)]
M12 = unique(as.character((df_freq[df_freq$Type == "12M",]$KO)))
M12 = M12[!is.na(M12)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "4M"=M4, "8M"=M8, "12M"=M12)
p_china_KO_venn = ggVennDiagram(x, label_size=3, label="count") +
    scale_fill_gradientn(name = "KO count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
    theme(
      text=element_text(size=12,family="Helvetica")
    )

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 202,]$item
groupB = venn_data2[venn_data2$count == 1527,]$item
groupC = venn_data2[venn_data2$count == 1583,]$item
groupD = venn_data2[venn_data2$count == 675,]$item
groupE = venn_data2[venn_data2$count == 292,]$item
groupF = venn_data2[venn_data2$count == 6082,]$item

abundance = data.frame(fread(my_gene_abundance_files_china[["12M"]]))
KO_df = data.frame(fread(my_KO_files_china[["12M"]], header=FALSE), check.names=FALSE)
colnames(KO_df)[1] = "gene_id"
colnames(KO_df)[3] = "KO"
groups = c(groupA,groupB,groupC,groupD,groupE,groupF)

df = NULL
group_names = c("A - 202","B - 1527","C - 1583","D - 675","E - 292","F - 6082")
for(i in 1:length(groups)){
  group = groups[[i]]
  KO_df2 = KO_df[KO_df$KO %in% group,]
  abundance2 = abundance[abundance$feature_id %in% KO_df2$gene_id,]
  tmp = data.frame(colSums(abundance2[,2:913]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the 12M reads configuration", width=30), size=8)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(8), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=8)

tab2_china = tab2

################################################
# ROCK KO                                      #
# Normalize with recA (from COG annotations)   #
################################################

#df_KO_all = get_KO_abundance_as_df(my_KO_files_rock, my_gene_abundance_files_rock, my_COG_files_rock, mapping_file="../rock_microbiome/export/mapping_file.tsv")
#fwrite(df_KO_all, paste0(outdir, "./df_KO_rock.tsv"), sep="\t", row.names=TRUE)
df_KO_all = data.frame(fread(paste0(outdir, "./df_KO_rock.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO_all) = df_KO_all$V1; df_KO_all$V1=NULL;

head(df_KO_all)

mapping = data.frame(fread("../rock_microbiome/export/mapping_file.tsv"), check.names=FALSE)

KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

types = unique(df_KO_all$Type)

x = 1

final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_KO_all[df_KO_all$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$KO
    df_x$KO = NULL
    df_y = df_KO_all[df_KO_all$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$KO
    df_y$KO = NULL
    
    curr_KOs_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_KOs_all = curr_KOs_all[curr_KOs_all != "NULL"]
    curr_KOs_x = unique(row.names(df_x))
    curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
    curr_KOs_y = unique(row.names(df_y))
    curr_KOs_y = curr_KOs_y[curr_KOs_y != "NULL"]
    
    missing_in_x = curr_KOs_all[!curr_KOs_all %in% curr_KOs_x]
    missing_in_y = curr_KOs_all[!curr_KOs_all %in% curr_KOs_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
   
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)

#write.table(df3, paste0(root, "/melted_style_spearman_KO_rock_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_KO_rock_df.tsv"), sep="\t"), check.names=FALSE)
df3[58,]$x = "0.1M - arm"
df3[59,]$x = "0.5M - arm"
df3[60,]$x = "1M - arm"
df3[61,]$x = "4M - arm"
df3[62,]$x = "8M - arm"
df3[63,]$x = "12M - arm"

df3[58,]$y = "complete"
df3[59,]$y = "complete"
df3[60,]$y = "complete"
df3[61,]$y = "complete"
df3[62,]$y = "complete"
df3[63,]$y = "complete"

df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm", "complete"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm", "complete"))

p_KO_rock = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) + #, width=0.9, height=0.9, color="gray")) +
  #facet_grid(LocationDepth ~ ., scales="free", space="free") +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  #ggtitle("Mantel test (Spearman coefficient) between Bray-Curtis dissimilarity of all datasets") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    #axis.title=element_text(family="Helvetica", size=14),
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_KO_rock)

##########################
# Investigate rare KOs   #
##########################
KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_KO_all$Type)

for(i in 1:length(types)){
  df_x = df_KO_all[df_KO_all$Type == types[i],]
  curr_KOs_x = unique(df_x$KO)
  curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
 
  missing_in_x = KOs[!KOs %in% curr_KOs_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_KOs = unique(df_x$KO)
  unique_KOs = unique_KOs[!unique_KOs %in% c("NULL", "0")]
  unique_KOs = data.frame(KO=unique_KOs, Type=types[i])
  freq_KOs = df_x[, c("KO", "sum")]
  freq_KOs = freq_KOs[!freq_KOs$KO %in% c("NULL", "0"),]
  freq_KOs$Type = types[i]
  
  if(x == 1){
    df_freq = freq_KOs
    df_unique = unique_KOs
  }else{
    df_freq = rbind(df_freq, freq_KOs)
    df_unique = rbind(df_unique, unique_KOs)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>% 
  mutate(Number_of_unique_KOs_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))

p_KO_rock_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_KOs_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Analysis configuration") + ylab(paste0("Number of unique KOs recovered")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_KO_rock_2

my_order = df_freq[df_freq$Type == "complete",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$KO
df_freq$KO = factor(df_freq$KO, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "4M", "4M_all", "8M", "8M_all", "12M", "12M_all", "complete"))

KO_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"),]
KO_rest = KO_rest[KO_rest$sum > 0,]
KO_rest_complete = KO_rest[KO_rest$Type %in% c("complete"),]
KO_rest2 = KO_rest[KO_rest$Type != "complete",]
KO_diff = KO_rest_complete[!KO_rest_complete$KO %in% KO_rest2$KO,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$KO)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$KO)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$KO)))
M1 = M1[!is.na(M1)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$KO)))
M4 = M4[!is.na(M4)]
M8 = unique(as.character((df_freq[df_freq$Type == "8M",]$KO)))
M8 = M8[!is.na(M8)]
M12 = unique(as.character((df_freq[df_freq$Type == "12M",]$KO)))
M12 = M12[!is.na(M12)]
Mcomplete = unique(as.character((df_freq[df_freq$Type == "complete",]$KO)))
Mcomplete = Mcomplete[!is.na(Mcomplete)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "4M"=M4, "8M"=M8, "12M"=M12, "complete"=Mcomplete)
p_rock_KO_venn = ggVennDiagram(x, label_size=3, label="count") +
  scale_fill_gradientn(name = "KO count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  theme(
    text=element_text(size=12,family="Helvetica")
  )
p_rock_KO_venn

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 767,]$item
groupB = venn_data2[venn_data2$count == 328,]$item
groupC = venn_data2[venn_data2$count == 202,]$item
groupD = venn_data2[venn_data2$count == 1648,]$item
groupE = venn_data2[venn_data2$count == 1114,]$item
groupF = venn_data2[venn_data2$count == 5001,]$item
groupG = venn_data2[venn_data2$count == 858,]$item

abundance = data.frame(fread(my_gene_abundance_files_rock[["complete"]]))
KO_df = data.frame(fread(my_KO_files_rock[["complete"]], header=FALSE), check.names=FALSE)
colnames(KO_df)[1] = "gene_id"
colnames(KO_df)[3] = "KO"
groups = c(groupA,groupB,groupC,groupD,groupE,groupF,groupG)

df = NULL
group_names = c("A - 767","B - 328","C - 202","D - 1648","E - 1114","F - 5001","G - 858")
for(i in 1:length(groups)){
  group = groups[[i]]
  KO_df2 = KO_df[KO_df$KO %in% group,]
  abundance2 = abundance[abundance$gene_id %in% KO_df2$gene_id,]
  tmp = data.frame(colSums(abundance2[,2:19]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the 'complete' reads configuration", width=27), size=7)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F,G letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(9), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=9)
tab2_rock = tab2

################################################
# AGR soil KO                                  #
# Normalize with recA (from COG annotations)   #
################################################

#df_KO_all = get_KO_abundance_as_df(my_KO_files_soil, my_gene_abundance_files_soil, 
#                                   my_COG_files_soil, mapping_file="../soil/export/mapping_file.tsv")
#fwrite(df_KO_all, paste0(outdir, "./df_KO_soil.tsv"), sep="\t", row.names=TRUE)
df_KO_all = data.frame(fread(paste0(outdir, "./df_KO_rock.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO_all) = df_KO_all$V1; df_KO_all$V1=NULL;

head(df_KO_all)

mapping = data.frame(fread("../rock_microbiome/export/mapping_file.tsv"), check.names=FALSE)

KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

types = unique(df_KO_all$Type)

x = 1

final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_KO_all[df_KO_all$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$KO
    df_x$KO = NULL
    df_y = df_KO_all[df_KO_all$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$KO
    df_y$KO = NULL
    
    curr_KOs_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_KOs_all = curr_KOs_all[curr_KOs_all != "NULL"]
    curr_KOs_x = unique(row.names(df_x))
    curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
    curr_KOs_y = unique(row.names(df_y))
    curr_KOs_y = curr_KOs_y[curr_KOs_y != "NULL"]
    
    missing_in_x = curr_KOs_all[!curr_KOs_all %in% curr_KOs_x]
    missing_in_y = curr_KOs_all[!curr_KOs_all %in% curr_KOs_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M", y="8M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M", y="12M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="8M_all", y="8M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="12M_all", y="12M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)

#write.table(df3, paste0(root, "/melted_style_spearman_KO_rock_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_KO_rock_df.tsv"), sep="\t"), check.names=FALSE)
df3[58,]$x = "0.1M - arm"
df3[59,]$x = "0.5M - arm"
df3[60,]$x = "1M - arm"
df3[61,]$x = "4M - arm"
df3[62,]$x = "8M - arm"
df3[63,]$x = "12M - arm"

df3[58,]$y = "complete"
df3[59,]$y = "complete"
df3[60,]$y = "complete"
df3[61,]$y = "complete"
df3[62,]$y = "complete"
df3[63,]$y = "complete"

df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm", "complete"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M - arm", "0.5M - arm", "1M - arm", "4M - arm", "8M - arm", "12M - arm", "complete"))

p_KO_soil = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) + #, width=0.9, height=0.9, color="gray")) +
  #facet_grid(LocationDepth ~ ., scales="free", space="free") +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  #ggtitle("Mantel test (Spearman coefficient) between Bray-Curtis dissimilarity of all datasets") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    #axis.title=element_text(family="Helvetica", size=14),
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_KO_soil)

##########################
# Investigate rare KOs   #
##########################
KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_KO_all$Type)

for(i in 1:length(types)){
  df_x = df_KO_all[df_KO_all$Type == types[i],]
  curr_KOs_x = unique(df_x$KO)
  curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
  
  missing_in_x = KOs[!KOs %in% curr_KOs_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_KOs = unique(df_x$KO)
  unique_KOs = unique_KOs[!unique_KOs %in% c("NULL", "0")]
  unique_KOs = data.frame(KO=unique_KOs, Type=types[i])
  freq_KOs = df_x[, c("KO", "sum")]
  freq_KOs = freq_KOs[!freq_KOs$KO %in% c("NULL", "0"),]
  freq_KOs$Type = types[i]
  
  if(x == 1){
    df_freq = freq_KOs
    df_unique = unique_KOs
  }else{
    df_freq = rbind(df_freq, freq_KOs)
    df_unique = rbind(df_unique, unique_KOs)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>% 
  mutate(Number_of_unique_KOs_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))

p_KO_soil_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_KOs_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Analysis configuration") + ylab(paste0("Number of unique KOs recovered")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_KO_soil_2

my_order = df_freq[df_freq$Type == "complete",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$KO
df_freq$KO = factor(df_freq$KO, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "4M", "4M_all", "8M", "8M_all", "12M", "12M_all", "complete"))

KO_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"),]
KO_rest = KO_rest[KO_rest$sum > 0,]
KO_rest_complete = KO_rest[KO_rest$Type %in% c("complete"),]
KO_rest2 = KO_rest[KO_rest$Type != "complete",]
KO_diff = KO_rest_complete[!KO_rest_complete$KO %in% KO_rest2$KO,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$KO)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$KO)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$KO)))
M1 = M1[!is.na(M1)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$KO)))
M4 = M4[!is.na(M4)]
M8 = unique(as.character((df_freq[df_freq$Type == "8M",]$KO)))
M8 = M8[!is.na(M8)]
M12 = unique(as.character((df_freq[df_freq$Type == "12M",]$KO)))
M12 = M12[!is.na(M12)]
Mcomplete = unique(as.character((df_freq[df_freq$Type == "complete",]$KO)))
Mcomplete = Mcomplete[!is.na(Mcomplete)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "4M"=M4, "8M"=M8, "12M"=M12, "complete"=Mcomplete)
p_soil_KO_venn = ggVennDiagram(x, label_size=3, label="count") +
  scale_fill_gradientn(name = "KO count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  theme(
    text=element_text(size=12,family="Helvetica")
  )
p_soil_KO_venn

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 858,]$item
groupB = venn_data2[venn_data2$count == 328,]$item
groupC = venn_data2[venn_data2$count == 202,]$item
groupD = venn_data2[venn_data2$count == 1648,]$item
groupE = venn_data2[venn_data2$count == 1114,]$item
groupF = venn_data2[venn_data2$count == 5001,]$item
groupG = venn_data2[venn_data2$count == 25,]$item
groupH = venn_data2[venn_data2$count == 767,]$item

abundance = data.frame(fread(my_gene_abundance_files_soil[["complete"]]))
colnames(abundance)[1] = "gene_id"
KO_df = data.frame(fread(my_KO_files_soil[["complete"]], header=FALSE), check.names=FALSE)
colnames(KO_df)[1] = "gene_id"
colnames(KO_df)[3] = "KO"
groups = NULL
groups = c(groupA,groupB,groupC,groupD,groupE,groupF,groupG,groupH)

df = NULL
group_names = c("A - 858","B - 328","C - 202","D - 1648","E - 1114","F - 5001","G - 25", "H - 767")
for(i in 1:(length(groups)-0)){
  print(i)
  group = groups[[i]]
  KO_df2 = KO_df[KO_df$KO %in% group,]
  abundance2 = abundance[abundance$gene_id %in% KO_df2$gene_id,]
  tmp = data.frame(colSums(abundance2[,2:73]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the 'complete' reads configuration", width=27), size=7)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F,G letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(9), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=9)
tab2_soil = tab2

################################################
# Mock communities KO                          #
# Normalize with recA (from COG annotations)   #
################################################

df_KO_all = get_KO_abundance_as_df(my_KO_files_mock, my_gene_abundance_files_mock, 
                                   my_COG_files_mock, mapping_file="../mock/export/mapping_file.tsv")
#fwrite(df_KO_all, paste0(outdir, "./df_KO_mock.tsv"), sep="\t", row.names=TRUE)
#df_KO_all = data.frame(fread(paste0(outdir, "./df_KO_mock.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO_all) = df_KO_all$V1; df_KO_all$V1=NULL;

head(df_KO_all)

mapping = data.frame(fread("../mock/export/mapping_file.tsv"), check.names=FALSE)

KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

types = unique(df_KO_all$Type)

x = 1

final_df = data.frame(matrix(1, nrow = length(types), ncol = length(types)))
colnames(final_df) = types
row.names(final_df) = types
df=NULL
df3=NULL
tmp_df=NULL
tmp_df3=NULL
for(i in 1:length(types)){
  if(i == length(types)){ break; }
  for(j in (i+1):length(types)){
    df_x = df_KO_all[df_KO_all$Type == types[i],]
    df_x$Type = NULL
    row.names(df_x) = df_x$KO
    df_x$KO = NULL
    df_y = df_KO_all[df_KO_all$Type == types[j],]
    df_y$Type = NULL
    row.names(df_y) = df_y$KO
    df_y$KO = NULL
    
    curr_KOs_all = c(unique(row.names(df_x)), unique(row.names(df_y)))
    curr_KOs_all = curr_KOs_all[curr_KOs_all != "NULL"]
    curr_KOs_x = unique(row.names(df_x))
    curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
    curr_KOs_y = unique(row.names(df_y))
    curr_KOs_y = curr_KOs_y[curr_KOs_y != "NULL"]
    
    missing_in_x = curr_KOs_all[!curr_KOs_all %in% curr_KOs_x]
    missing_in_y = curr_KOs_all[!curr_KOs_all %in% curr_KOs_y]
    
    add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
    colnames(add_df_x) = colnames(df_x)
    add_df_x[is.na(add_df_x)] <- 0
    df_x = rbind(df_x, add_df_x)
    
    add_df_y = data.frame(matrix(ncol=ncol(df_y), nrow=length(missing_in_y)), row.names=missing_in_y)
    colnames(add_df_y) = colnames(df_y)
    add_df_y[is.na(add_df_y)] <- 0
    df_y = rbind(df_y, add_df_y)
    
    df_x = df_x[order(row.names(df_x)), ]
    df_y = df_y[order(row.names(df_y)), ]
    
    if(!identical(row.names(df_x), row.names(df_y))){ stop("Something wrong with ordering of data.frames.") }
    
    name_x = types[i]
    name_y = types[j]
    dist_x = vegdist(t(df_x), method="bray")
    dist_y = vegdist(t(df_y), method="bray")
    res = mantel(dist_x, dist_y, method = "spearman", permutations = 999, na.rm = TRUE)
    
    tmp_df = data.frame(comparison=paste0(name_x, " vs ", name_y), r=res$statistic, signif=res$signif)
    tmp_df3 = data.frame(x=name_x, y=name_y, r=res$statistic, signif=res$signif)
    print(paste0(name_x, " - ", name_y, " - ", "no:", x))
    
    if(x == 1){
      df = tmp_df
      df3 = tmp_df3
    }else{
      df = rbind(df, tmp_df)
      df3 = rbind(df3, tmp_df3)
    }
    
    #Then populate df to create a dist matrix like df.
    final_df[[name_x, name_y]] = res$statistic
    final_df[[name_y, name_x]] = res$statistic
    x = x + 1
  }
}
df3 = rbind(df3, data.frame(x="0.1M", y="0.1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M", y="0.5M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M", y="1M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M", y="2M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M", y="3M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M", y="4M", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.1M_all", y="0.1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="0.5M_all", y="0.5M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="1M_all", y="1M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="2M_all", y="2M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="3M_all", y="3M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="4M_all", y="4M_all", r=1, signif=0))
df3 = rbind(df3, data.frame(x="complete", y="complete", r=1, signif=0))

# Then heatmap
df3$x = gsub("_all", " - arm", df3$x)
df3$y = gsub("_all", " - arm", df3$y)

#write.table(df3, paste0(root, "/melted_style_spearman_KO_rock_df.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#df3 = data.frame(fread(paste0(root, "/melted_style_spearman_KO_rock_df.tsv"), sep="\t"), check.names=FALSE)
df3[58,]$x = "0.1M - arm"
df3[59,]$x = "0.5M - arm"
df3[60,]$x = "1M - arm"
df3[61,]$x = "2M - arm"
df3[62,]$x = "3M - arm"
df3[63,]$x = "4M - arm"

df3[58,]$y = "complete"
df3[59,]$y = "complete"
df3[60,]$y = "complete"
df3[61,]$y = "complete"
df3[62,]$y = "complete"
df3[63,]$y = "complete"

df3$x = factor(df3$x, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm", "complete"))
df3$y = factor(df3$y, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "0.1M - arm", "0.5M - arm", "1M - arm", "2M - arm", "3M - arm", "4M - arm", "complete"))

p_KO_mock = ggplot(df3, aes(x=y, y=x,label=round(r,3))) +
  geom_tile(aes(fill=r)) + #, width=0.9, height=0.9, color="gray")) +
  #facet_grid(LocationDepth ~ ., scales="free", space="free") +
  scale_fill_gradientn(colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  ylab("") + xlab("") + labs(fill="r statistic") +
  geom_text(size=3,color="white", fontface=2) +
  #ggtitle("Mantel test (Spearman coefficient) between Bray-Curtis dissimilarity of all datasets") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=10, colour="black", angle=as.numeric(90), hjust=1, vjust=0.5),
    axis.text.y=element_text(size=10, colour="black"), 
    #axis.title=element_text(family="Helvetica", size=14),
    plot.title = element_text(lineheight=1.2, face="plain", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size=10, face="plain"),
    legend.title = element_text(size=12, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=9, face="bold"),
    strip.text.y = element_text(angle=0, hjust=0, size=7, face="bold"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
print(p_KO_mock)

##########################
# Investigate rare KOs   #
##########################
KOs = unique(df_KO_all$KO)
KOs = KOs[KOs != "NULL"]

x = 1
df_freq = NULL
df_unique = NULL
types = unique(df_KO_all$Type)

for(i in 1:length(types)){
  df_x = df_KO_all[df_KO_all$Type == types[i],]
  curr_KOs_x = unique(df_x$KO)
  curr_KOs_x = curr_KOs_x[curr_KOs_x != "NULL"]
  
  missing_in_x = KOs[!KOs %in% curr_KOs_x]
  
  add_df_x = data.frame(matrix(ncol=ncol(df_x), nrow=length(missing_in_x)), row.names=missing_in_x)
  colnames(add_df_x) = colnames(df_x)
  add_df_x[is.na(add_df_x)] <- 0
  df_x = rbind(df_x, add_df_x)
  df_x$Type = NULL
  df_x$sum = rowSums(df_x[,1:(ncol(df_x)-1)])
  unique_KOs = unique(df_x$KO)
  unique_KOs = unique_KOs[!unique_KOs %in% c("NULL", "0")]
  unique_KOs = data.frame(KO=unique_KOs, Type=types[i])
  freq_KOs = df_x[, c("KO", "sum")]
  freq_KOs = freq_KOs[!freq_KOs$KO %in% c("NULL", "0"),]
  freq_KOs$Type = types[i]
  
  if(x == 1){
    df_freq = freq_KOs
    df_unique = unique_KOs
  }else{
    df_freq = rbind(df_freq, freq_KOs)
    df_unique = rbind(df_unique, unique_KOs)
  }
  x=x+1
}

df_unique2 = df_unique %>% 
  count(Type) %>%
  dplyr::group_by(Type) %>% 
  mutate(Number_of_unique_KOs_recovered = sum(n)) %>% 
  as.data.frame()

df_unique2 = df_unique2[!grepl("_all", df_unique2$Type),]
df_unique2$Type = factor(df_unique2$Type, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete"))

p_KO_mock_2 <- ggplot(data=df_unique2, aes(x=Type, y=Number_of_unique_KOs_recovered)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Analysis configuration") + ylab(paste0("Number of unique KOs recovered")) +
  theme(
    text=element_text(family="Helvetica"),
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=0.5),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title=element_text(family="Helvetica", size=(11)),
    plot.title = element_text(lineheight=1.2, size=11, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.key = element_rect(fill = "NA"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=11, face="plain"),
    legend.spacing = unit(1, "cm"),
    legend.position="right",
    strip.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=0, face="plain"),
    strip.text.y = element_text(angle=90, hjust=1, size=11, face="plain"),
    strip.background =  element_blank()
  ) +  scale_fill_manual(values=(vColors))
p_KO_mock_2

my_order = df_freq[df_freq$Type == "complete",]
my_order = my_order[order(my_order$sum),]
my_order = my_order$KO
df_freq$KO = factor(df_freq$KO, levels=my_order)
df_freq$Type = factor(df_freq$Type, levels=c("0.1M", "0.1M_all","0.5M", "0.5M_all", "1M", "1M_all", "2M", "2M_all", "3M", "3M_all", "4M", "4M_all", "complete"))

KO_rest = df_freq[df_freq$Type %in% c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete"),]
KO_rest = KO_rest[KO_rest$sum > 0,]
KO_rest_complete = KO_rest[KO_rest$Type %in% c("complete"),]
KO_rest2 = KO_rest[KO_rest$Type != "complete",]
KO_diff = KO_rest_complete[!KO_rest_complete$KO %in% KO_rest2$KO,]

M01 = unique(as.character((df_freq[df_freq$Type == "0.1M",]$KO)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((df_freq[df_freq$Type == "0.5M",]$KO)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((df_freq[df_freq$Type == "1M",]$KO)))
M1 = M1[!is.na(M1)]
M2 = unique(as.character((df_freq[df_freq$Type == "2M",]$KO)))
M2 = M4[!is.na(M2)]
M3 = unique(as.character((df_freq[df_freq$Type == "3M",]$KO)))
M3 = M3[!is.na(M3)]
M4 = unique(as.character((df_freq[df_freq$Type == "4M",]$KO)))
M4 = M4[!is.na(M4)]
Mcomplete = unique(as.character((df_freq[df_freq$Type == "complete",]$KO)))
Mcomplete = Mcomplete[!is.na(Mcomplete)]

x = list("0.1M"=M01, "0.5M"=M05, "1M"=M1, "2M"=M2, "3M"=M3, "4M"=M4, "complete"=Mcomplete)
p_mock_KO_venn = ggVennDiagram(x, label_size=3, label="count") +
  scale_fill_gradientn(name = "KO count", colors=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100)) +
  theme(
    text=element_text(size=12,family="Helvetica")
  )
p_mock_KO_venn

# extract the KOs
venn_data <- process_data(Venn(x))
venn_data2 = data.frame(venn_region(venn_data))
groupA = venn_data2[venn_data2$count == 2494,]$item
groupB = venn_data2[venn_data2$count == 2723,]$item
groupC = venn_data2[venn_data2$count == 475,]$item
groupD = venn_data2[venn_data2$count == 99,]$item
groupE = venn_data2[venn_data2$name == "0.5M..2M..3M..4M..complete",]$item
groupF = venn_data2[venn_data2$name == "1M..3M..4M..complete",]$item
groupG = venn_data2[venn_data2$name == "2M..4M..complete",]$item
groupH = venn_data2[venn_data2$name == "1M..complete",]$item

abundance = data.frame(fread(my_gene_abundance_files_mock[["complete"]]))
colnames(abundance)[1] = "gene_id"
KO_df = data.frame(fread(my_KO_files_mock[["complete"]], header=FALSE), check.names=FALSE)
colnames(KO_df)[1] = "gene_id"
colnames(KO_df)[3] = "KO"
groups=NULL
groups = c(groupA,groupB,groupC,groupD,groupE,groupF,groupG,groupH)

df = NULL
group_names = c("A - 2494","B - 2723","C - 475","D - 99","E - 1","F - 1","G - 1","H - 1")
for(i in 1:(length(groups)-1)){
  print(i)
  group = groups[[i]]
  KO_df2 = KO_df[KO_df$KO %in% group,]
  abundance2 = abundance[abundance$gene_id %in% KO_df2$gene_id,]
  tmp = data.frame(colSums(abundance2[,2:7]))
  tmp = tmp[order(row.names(tmp)),,drop=FALSE]
  colnames(tmp) = group_names[i]
  if(i == 1){
    df = tmp
  }else{
    df = cbind(df, tmp)
  }
}

tmp = data.frame(colSums(df))
tmp2 = data.frame(prop.table(data.frame(colSums(df))))
tmp2 = format(tmp2*100, digits=1, scientific=FALSE)
tmp = cbind(tmp, tmp2)
colnames(tmp) = c("Number of reads", "%")
tmp$`Number of reads` = format(tmp$`Number of reads`, digits=0, big.mark=",", scientific=FALSE)
tmp

# then put this in a table as C in the figure.
# Summary table plot, medium orange theme
tab <- ggtexttable(tmp, theme = ttheme(base_size=7, base_style="blank", padding = unit(c(1, 1), "mm"),))
tab <- tab_add_title(tab, str_wrap("Distribution of reads for the 'complete' reads configuration", width=27), size=7)
tab = tab_add_footnote(tab, str_wrap("A,B,C,D,E,F,G letters correspond to the Venn diagram zones.", width=30), size=6)
tab2 = tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(2, 3), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(9), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 3, linetype = 2, to.row=9)
tab2_mock = tab2


#
# Then let's generate the figure.
#
# plots china
figure7a1 <- ggarrange(
  p_KO_china_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1,0.9),
  widths = c(1, 1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7a2 <- ggarrange(
  p_KO_china, figure7a1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7a3 = ggarrange(
  p_china_KO_venn, tab2_china,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7a4 = ggarrange(
  figure7a2, figure7a3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f7a = annotate_figure(figure7a4, 
                      fig.lab = "\n       Human gut", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# plots Rock
figure7b1 <- ggarrange(
  p_KO_rock_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  
  heights = c(1,0.9),
  widths = c(1,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7b2 <- ggarrange(
  p_KO_rock, figure7b1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7b3 = ggarrange(
  p_rock_KO_venn, tab2_rock,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7b4 = ggarrange(
  figure7b2, figure7b3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f7b = annotate_figure(figure7b4, 
                      fig.lab = "\n       Antarctic", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

# plots agr soil
figure7c1 <- ggarrange(
  p_KO_soil_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  
  heights = c(1,0.9),
  widths = c(1,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7c2 <- ggarrange(
  p_KO_soil, figure7c1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7c3 = ggarrange(
  p_soil_KO_venn, tab2_soil,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7c4 = ggarrange(
  figure7c2, figure7c3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f7c = annotate_figure(figure7c4, 
                      fig.lab = "\n       Agricultural soil", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

#plots mock
figure7d1 <- ggarrange(
  p_KO_mock_2,
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = FALSE, legend = NULL,
  
  heights = c(1,0.9),
  widths = c(1,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7d2 <- ggarrange(
  p_KO_mock, figure7d1,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(1.5, 1.23),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7d3 = ggarrange(
  p_mock_KO_venn, tab2_mock,
  labels = NULL,
  ncol = 2, nrow = 1,
  common.legend = FALSE, legend = NULL,
  widths = c(2.5,1),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

figure7d4 = ggarrange(
  figure7d2, figure7d3,
  labels = NULL,
  ncol = 1, nrow = 2,
  common.legend = FALSE, legend = NULL,
  heights = c(1, 1.6),
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f7d = annotate_figure(figure7d4, 
                      fig.lab = "\n       Mock communities", 
                      fig.lab.pos = "top.left",
                      top = textGrob("\n", gp = gpar(cex = 0.9), hjust=1.8)
)

final_list = annotate_figure(
  ggarrange(f7a,f7b,f7c,f7d,
            labels=c("A", "B", "C", "D"),
            nrow = 4,
            heights=c(1.2, 1.2, 1.2, 1.2)
  ),
  top = text_grob("Figure S5\n", color = "black", face = "plain", size = 11),
)

pdf( file=paste0(outdir, "./figures/FIGURE_S5.pdf"), height=26, width=12)
final_list
dev.off()



################################################
# Do taxonomic summaries for human gut         #
# + antarctic microbiome                       #
# Do selected taxa only                        #
#  and put complete tax summaries in figure S5 #
# Alpha div is already reported in figure 5a   #
# Here a) pcoa BC figures, b) selected taxa    #
# profiles, c) KO heatmap, d) KO+taxa barplot  #  
# Figure 4.                                    #
################################################

##################################
## Human gut                    ##
##################################
mapping = data.frame(fread("./mapping_file3.tsv", sep="\t"), check.names=FALSE)

# Beta div
# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diversity
df_coords = get_coords_as_df(my_coords_files)
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files, mapping=mapping)
# include taxa profile stacked bar plots based on alpha diversity quintiles.
df_taxa = get_tax_as_df_unmelted(my_tax_files_china_species, convert_to_rel=FALSE, mapping_file="./mapping_file.tsv")

selected_taxa = c(
  "g__Coprococcus$",
  "g__Bifidobacterium$",
  "s__[Eubacterium] rectale$",
  "g__Megamonas$",
  "g__Clostridium$",
  "s__Prevotella copri$",
  "g__Bacteroides$",
  "g__Prevotella$",
  "o__Clostridiales$"
)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all")

ps_ts = list()
curr_taxa2 = NULL
# for simplicity sake, only show quartile 1 and 5
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_ad = df_alpha[df_alpha$Type == curr_workflow,]
  
  # Divide by quintiles:
  curr_ad$quantiles = NULL
  curr_ad$quantiles = cut(curr_ad$value, 
                          quantile(
                            curr_ad$value, 
                            prob = seq(0, 1, length = 6), 
                            type = 5, names=TRUE, 
                            labels=FALSE, 
                            ordered_result=TRUE
                          )
  )
  
  df_link = data.frame(quantiles=levels(curr_ad$quantiles), rank=seq(1,(length(unique(curr_ad$quantiles))-1), by=1))
  df_tmp = curr_ad[,c("variable", "quantiles")]
  df_tmp = merge(df_tmp, df_link, by="quantiles")
  curr_mapping = merge(mapping, df_tmp, by.x="#SampleID", by.y="variable")
  colnames(curr_mapping)[ncol(curr_mapping)] = "Diversity_quartile"
  
  # Taxa
  curr_taxa = df_taxa[df_taxa$Type == curr_workflow,]
  curr_taxa$Type = NULL
  # clean taxa lineages:
  curr_taxa$Taxon = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL", "k__Bacteria", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("p__NULL;c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";c__NULL", ";c__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("k__Bacteria;", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";p__NULL;", ";p__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";s__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";f__NULL;g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";$", "", curr_taxa$Taxon)  
  
  # Then only keep selection
  curr_taxa = curr_taxa[grepl(paste(selected_taxa, collapse="|"), curr_taxa$Taxon),]
  # reorder by total abundance
  curr_taxa = cbind(curr_taxa, (rowSums(curr_taxa[1:(ncol(curr_taxa)-1)]))/ncol(curr_taxa) )
  curr_taxa = curr_taxa[order(-curr_taxa[, ncol(curr_taxa)]),,drop=FALSE]
  curr_taxa[,ncol(curr_taxa)] = NULL
  taxa_order = curr_taxa$Taxon
  
  curr_taxa = melt(curr_taxa)
  curr_taxa = merge(curr_taxa, curr_mapping, by.x="variable", by.y="#SampleID")
  
  # sort order of samples based on taxa order.
  
  tmp2 = dcast(curr_taxa, variable ~ Taxon, value.var="value")
  tmp2$variable = as.character(tmp2$variable)
  tmp3 = tmp2[with(tmp2, order(-`p__Firmicutes;c__Clostridia;o__Clostridiales`,
                               #-`p__Firmicutes`,
                               -`p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Selenomonadaceae;g__Megamonas`
                               #-`p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia`
  )),]
  sorted_variables = as.character(tmp3$variable)
  curr_taxa$variable = factor(curr_taxa$variable, levels=sorted_variables)
  
  # manually specify the colors of the bars
  color_list = c(
    "p__Firmicutes;c__Clostridia;o__Clostridiales" = "#0000CD",
    "p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium" = "#00FF00", 
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides" = "#FF0000",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella" = "#000000",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella;s__Prevotella copri" = "gray",
    "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium" = "#B22222",                       
    "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus" = "#FFA500",                                                        
    "p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Selenomonadaceae;g__Megamonas" = "#FF00FF"
  )
  
  # override taxa order
  taxa_order = c(
    "p__Firmicutes;c__Clostridia;o__Clostridiales",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella;s__Prevotella copri",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides",
    "p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Selenomonadaceae;g__Megamonas",
    "p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium", 
    "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium",                       
    "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus"
  )
  abundant_in_low_diversity = c(
    "p__Firmicutes;c__Clostridia;o__Clostridiales",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella;s__Prevotella copri",
    "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides",
    "p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Selenomonadaceae;g__Megamonas",
    "p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium"
  )
  abundant_in_high_diversity = c(
    "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium",                       
    "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus"
  )
  curr_taxa$high_low = ifelse(curr_taxa$Taxon %in% abundant_in_high_diversity, "Abundant in\nhigh diversity\n (quintile #1)", "Abundant in\nlow diversity\n (quintile #5)")
  curr_taxa$Taxon = factor(curr_taxa$Taxon, levels=(taxa_order))
  curr_taxa$value = curr_taxa$value * 100
  curr_taxa = curr_taxa[curr_taxa$Diversity_quartile %in% c("1","5"),]
  curr_taxa$Diversity_quartile = gsub("1", "Quintile #1\n(Participants harboring\n low diversity microbiota)", curr_taxa$Diversity_quartile)
  curr_taxa$Diversity_quartile = gsub("5", "Quintile #5\n(Participants harboring\n high diversity microbiota)", curr_taxa$Diversity_quartile)
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_taxa$workflow = curr_workflow2
  curr_taxa$workflow2 = ifelse(grepl("arm", curr_taxa$workflow), "arm", "standard")
  curr_taxa$number_of_clusters = gsub("^(.*M).*", "\\1", curr_taxa$workflow)
  # populate filtered taxa df for further usage.
  if(i == 1){
    curr_taxa2 = curr_taxa
  }else{
    curr_taxa2 = rbind(curr_taxa2, curr_taxa)
  }
  
  p <- ggplot(data=curr_taxa, aes(x=Taxon, y=value, fill=Taxon)) + 
    facet_grid(high_low ~ Diversity_quartile, scales="free_y", space="free_x") + 
    geom_boxplot(aes(fill=Taxon), outlier.shape = NA, position=position_dodge()) + 
    geom_jitter(aes(color=Taxon),size=0.005, pch=19,position=position_jitterdodge(1.45), alpha=1) +
    xlab("") + 
    ylab("abundance (%)") +
    scale_fill_manual(values=color_list) +
    scale_color_manual(values=color_list) +
    guides(fill=guide_legend(ncol=2)) + ggtitle(curr_workflow2) +
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=10, colour="black"), 
      axis.title=element_text(family="Helvetica", size=11),
      plot.title = element_text(face="plain", size=14, hjust=0.5),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.07, 'cm'),
      legend.position="bottom",
      strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
      strip.background =  element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.3, "lines")
    ) 
  
  if(i == 3 | i == 6 | i == 9 | i == 12){ 
    p = p + theme(strip.text.y = element_text(angle=90, vjust=0, hjust=0.5, size=11, face="plain"))
  }else{
    p = p + theme(strip.text.y = element_text(size=0))
  }
  ps_ts[[i]] = p
  
}

color_list = c("arm" = "#0000CD", "standard" = "#FF0000")
taxa_order = c(
  "p__Firmicutes;c__Clostridia;o__Clostridiales",
  "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella",
  "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella;s__Prevotella copri",
  "p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides",
  "p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Selenomonadaceae;g__Megamonas",
  "p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium", 
  "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium",                       
  "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus"
)
curr_taxa2$Taxon = gsub("p__Firmicutes;c__Clostridia;o__Clostridiales", "p__Firmicutes;c__Clostridia;\no__Clostridiales", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella", "p__Bacteroidetes;c__Bacteroidia;\no__Bacteroidales;f__Prevotellaceae;\ng__Prevotella", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella;s__Prevotella copri", "p__Bacteroidetes;c__Bacteroidia;\no__Bacteroidales;f__Prevotellaceae;\ng__Prevotella;s__Prevotella copri", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides", "p__Bacteroidetes;c__Bacteroidia;\no__Bacteroidales;f__Bacteroidaceae;\ng__Bacteroides", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Firmicutes;c__Negativicutes;o__Selenomonadales;f__Selenomonadaceae;g__Megamonas", "p__Firmicutes;c__Negativicutes;\no__Selenomonadales;f__Selenomonadaceae;\ng__Megamonas", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium", "p__Actinobacteria;c__Actinobacteria;\no__Bifidobacteriales;f__Bifidobacteriaceae;\ng__Bifidobacterium", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Firmicutes;c__Clostridia;\no__Clostridiales;f__Clostridiaceae;g__Clostridium", "p__Firmicutes;c__Clostridia;\no__Clostridiales;f__Clostridiaceae;\ng__Clostridium", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Firmicutes;c__Clostridia;\no__Clostridiales;f__Lachnospiraceae;g__Coprococcus", "p__Firmicutes;c__Clostridia;\no__Clostridiales;f__Lachnospiraceae;\ng__Coprococcus", curr_taxa2$Taxon)

curr_taxa2$number_of_clusters = factor(curr_taxa2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))
# Try to do all at once:
p_selected_tax_china <- ggplot(data=curr_taxa2, aes(x=number_of_clusters, y=value, fill=workflow2)) + 
  facet_grid(Taxon ~ Diversity_quartile, scales="free_y", space="free_x") + 
  geom_boxplot(aes(fill=workflow2), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow2),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  xlab("Number of clusters") + 
  ylab("abundance (%)") +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  guides(fill=guide_legend(ncol=2)) + ggtitle("") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title=element_text(family="Helvetica", size=11),
    plot.title = element_text(face="plain", size=14, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.07, 'cm'),
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
    strip.background =  element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(angle=0, vjust=0.5, hjust=0, size=10, face="plain")
  ) 
p_selected_tax_china

#################################
## Antarctic                    #
#################################
mapping = data.frame(fread("../rock_microbiome/0.1M_clusters/mapping_file.tsv", sep="\t"), check.names=FALSE)
mapping$Location = gsub(".nord", "", mapping$Treatment)
mapping$Location = gsub(".sud", "", mapping$Location)

# Beta div
# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diversity
df_coords = get_coords_as_df(my_coords_files_rock)
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files_rock, mapping=mapping)
# include taxa profile stacked bar plots based on alpha diversity quintiles.
df_taxa = get_tax_as_df_unmelted(my_tax_files_rock_species, convert_to_rel=FALSE, mapping_file="../rock_microbiome/0.1M_clusters/mapping_file.tsv")

# Keep most abundant 20 taxa.
df_taxa2 = melt(df_taxa)
df_taxa2 = df_taxa2 %>% 
  dplyr::group_by(Taxon) %>% 
  dplyr::summarise_at(vars(value), list(SUM = ~ sum(.), 
                                        SD = ~sd(.))) %>% 
  as.data.frame()
df_taxa2 = df_taxa2[order(-df_taxa2$SUM),]
selected_taxa = rev(unique(df_taxa2$Taxon[1:20]))

selected_taxa = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL", "k__Bacteria", selected_taxa)
selected_taxa = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL", "k__Bacteria", selected_taxa)
selected_taxa = gsub("p__NULL;c__NULL;o__NULL", "", selected_taxa)
selected_taxa = gsub("c__NULL;o__NULL", "", selected_taxa)
selected_taxa = gsub(";o__NULL", "", selected_taxa)
selected_taxa = gsub(";c__NULL", ";c__undefined", selected_taxa)
selected_taxa = gsub("k__Bacteria;", "", selected_taxa)
selected_taxa = gsub(";p__NULL;", ";p__undefined", selected_taxa)
selected_taxa = gsub(";s__NULL", "", selected_taxa)
selected_taxa = gsub(";f__NULL;g__NULL", "", selected_taxa)
selected_taxa = gsub(";g__NULL", "", selected_taxa)
selected_taxa = gsub(";$", "", selected_taxa)  

my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all", "complete")

ps_ts = list()
curr_taxa2 = NULL
# for simplicity sake, only show quartile 1 and 5
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
 
  # Taxa
  curr_taxa = df_taxa[df_taxa$Type == curr_workflow,]
  curr_taxa$Type = NULL
  # clean taxa lineages:
  curr_taxa$Taxon = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL;f__NULL;g__NULL;s__NULL", "k__Bacteria", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL", "k__Bacteria", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("p__NULL;c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";c__NULL", ";c__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("k__Bacteria;", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";p__NULL;", ";p__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";s__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";f__NULL;g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";$", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("p__Chloroflexi;;s__Chloroflexi bacterium", "p__Chloroflexi", curr_taxa$Taxon)
  
  # Then only keep selection
  curr_taxa = curr_taxa[grepl(paste(selected_taxa, collapse="|"), curr_taxa$Taxon),]
  # reorder by total abundance
  curr_taxa = cbind(curr_taxa, (rowSums(curr_taxa[1:(ncol(curr_taxa)-1)]))/ncol(curr_taxa) )
  curr_taxa = curr_taxa[order(-curr_taxa[, ncol(curr_taxa)]),,drop=FALSE]
  curr_taxa[,ncol(curr_taxa)] = NULL
  taxa_order = curr_taxa$Taxon
  
  curr_taxa = melt(curr_taxa)
  curr_taxa = merge(curr_taxa, mapping, by.x="variable", by.y="#SampleID")
  curr_taxa = curr_taxa[curr_taxa$Location %in% c("University.Valley", "Siegfried.Peak") ,]
  
  # sort order of samples based on taxa order.
  taxa_order = c(
    "k__Fungi;p__Ascomycota",                                                                                                   
    "k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Umbilicariales;f__Umbilicariaceae;g__Lasallia;s__Lasallia pustulata",         
    "p__Chloroflexi",                                                                                                           
    "k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Lecanorales",                                                                 
    "k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Lecanorales;f__Parmeliaceae;g__Bryoria;s__Bryoria tenuis",                    
    "p__Chloroflexi;c__Ktedonobacteria;o__Ktedonobacterales;f__Ktedonobacteraceae;g__Ktedonobacter;s__Ktedonobacter racemifer",
    "p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Intrasporangiaceae;g__Janibacter;s__Janibacter limosus"
  )
  
  curr_taxa$Taxon = factor(curr_taxa$Taxon, levels=(taxa_order))
  curr_taxa$value = curr_taxa$value * 100
  
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_taxa$workflow = curr_workflow2
  curr_taxa$workflow2 = ifelse(grepl("arm", curr_taxa$workflow), "arm", "standard")
  curr_taxa$number_of_clusters = gsub("^(.*M).*", "\\1", curr_taxa$workflow)
  # populate filtered taxa df for further usage.
  if(i == 1){
    curr_taxa2 = curr_taxa
  }else{
    curr_taxa2 = rbind(curr_taxa2, curr_taxa)
  }
  
  p <- ggplot(data=curr_taxa, aes(x=Taxon, y=value, fill=Taxon)) + 
    facet_grid(. ~ Location, scales="free_y", space="free_x") + 
    geom_boxplot(aes(fill=Taxon), outlier.shape = NA, position=position_dodge()) + 
    geom_jitter(aes(color=Taxon),size=0.005, pch=19,position=position_jitterdodge(1.45), alpha=1) +
    xlab("") + 
    ylab("abundance (%)") +
    scale_fill_manual(values=color_list) +
    scale_color_manual(values=color_list) +
    guides(fill=guide_legend(ncol=2)) + ggtitle(curr_workflow2) +
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=10, colour="black"), 
      axis.title=element_text(family="Helvetica", size=11),
      plot.title = element_text(face="plain", size=14, hjust=0.5),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.07, 'cm'),
      legend.position="bottom",
      strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
      strip.background =  element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.3, "lines")
    ) 
  
  if(i == 3 | i == 6 | i == 9 | i == 12){ 
    p = p + theme(strip.text.y = element_text(angle=90, vjust=0, hjust=0.5, size=11, face="plain"))
  }else{
    p = p + theme(strip.text.y = element_text(size=0))
  }
  ps_ts[[i]] = p
  
}

color_list = c("arm" = "#0000CD", "standard" = "#FF0000")
taxa_order = c(
  "k__Fungi;p__Ascomycota",                                                                                                   
  "k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Umbilicariales;f__Umbilicariaceae;g__Lasallia;s__Lasallia pustulata",         
  "p__Chloroflexi",                                                                                                           
  "k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Lecanorales",                                                                 
  "k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Lecanorales;f__Parmeliaceae;g__Bryoria;s__Bryoria tenuis",                    
  "p__Chloroflexi;c__Ktedonobacteria;o__Ktedonobacterales;f__Ktedonobacteraceae;g__Ktedonobacter;s__Ktedonobacter racemifer", 
  "p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Intrasporangiaceae;g__Janibacter;s__Janibacter limosus"
)
curr_taxa2$Taxon = factor(curr_taxa2$Taxon, levels=taxa_order)
curr_taxa2$number_of_clusters = factor(curr_taxa2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))
curr_taxa2 = na.omit(curr_taxa2)
curr_taxa2$Taxon = gsub("k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Lecanorales;f__Parmeliaceae;g__Bryoria;s__Bryoria tenuis", 
                        "k__Fungi;p__Ascomycota;\nc__Lecanoromycetes;o__Lecanorales;\nf__Parmeliaceae;g__Bryoria;\ns__Bryoria tenuis", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Umbilicariales;f__Umbilicariaceae;g__Lasallia;s__Lasallia pustulata", 
                        "k__Fungi;p__Ascomycota;\nc__Lecanoromycetes;o__Umbilicariales;\nf__Umbilicariaceae;g__Lasallia;\ns__Lasallia pustulata", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("k__Fungi;p__Ascomycota;c__Lecanoromycetes;o__Lecanorales", 
                        "k__Fungi;p__Ascomycota;\nc__Lecanoromycetes;o__Lecanorales", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Chloroflexi;c__Ktedonobacteria;o__Ktedonobacterales;f__Ktedonobacteraceae;g__Ktedonobacter;s__Ktedonobacter racemifer", 
                        "p__Chloroflexi;c__Ktedonobacteria;\no__Ktedonobacterales;f__Ktedonobacteraceae;\ng__Ktedonobacter;s__Ktedonobacter racemifer", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Intrasporangiaceae;g__Janibacter;s__Janibacter limosus", 
                        "p__Actinobacteria;c__Actinobacteria\n;o__Micrococcales;f__Intrasporangiaceae;\ng__Janibacter;s__Janibacter limosus", curr_taxa2$Taxon)


# Try to do all at once:
curr_taxa2$Location = gsub("\\.", " ", curr_taxa2$Location)
p_selected_tax_rock <- ggplot(data=curr_taxa2, aes(x=number_of_clusters, y=value, fill=workflow2)) + 
  facet_grid(Taxon ~ Location, scales="free_y", space="free_x") +
  geom_boxplot(aes(fill=workflow2), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow2),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=1) +
  xlab("Number of clusters") + 
  ylab("abundance (%)") +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  guides(fill=guide_legend(ncol=2)) + ggtitle("") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title=element_text(family="Helvetica", size=9),
    plot.title = element_text(face="plain", size=14, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.07, 'cm'),
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
    strip.background =  element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(angle=0, vjust=0.5, hjust=0, size=10, face="plain")
  ) 
p_selected_tax_rock

#################################
## Agricultural soil            #
#################################
mapping = data.frame(fread("../soil/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diversity
df_coords = get_coords_as_df(my_coords_files_soil)
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files_soil, mapping=mapping)
# include taxa profile stacked bar plots based on alpha diversity quintiles.
df_taxa = get_tax_as_df_unmelted(my_tax_files_soil_species, convert_to_rel=FALSE, mapping_file="../soil/mapping_file.tsv")

tmp = data.frame(fread("../soil/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv", header=T), check.names=FALSE)
tmp = cbind(tmp, (rowSums(tmp[2:(ncol(tmp)-1)]))/(ncol(tmp)-1) )
tmp = tmp[order(-tmp[, ncol(tmp)]),,drop=FALSE]
tmp[,ncol(tmp)] = NULL

selected_taxa = c(
  "p__Actinobacteria;c__Actinobacteria",
  "p__Acidobacteria",
  "p__Proteobacteria;c__Alphaproteobacteria",
  "p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae",
  "p__Proteobacteria",
  "p__Acidobacteria;c__Acidobacteriia"
)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
my_workflows = c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "0.1M_all", "0.5M_all", "1M_all", "4M_all", "8M_all", "12M_all", "complete")

ps_ts = list()
curr_taxa2 = NULL
# for simplicity sake, only show quartile 1 and 5
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]

  
  # Taxa
  curr_taxa = df_taxa[df_taxa$Type == curr_workflow,]
  curr_taxa$Type = NULL
  # clean taxa lineages:
  curr_taxa$Taxon = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL", "k__Bacteria", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("p__NULL;c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";c__NULL", ";c__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("k__Bacteria;", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";p__NULL;", ";p__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";s__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";f__NULL;g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";$", "", curr_taxa$Taxon)  
  curr_taxa$Taxon = gsub("^k__NULL$", "k__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";;s__Acidobacteria bacterium", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("^f__NULL$", "f__undefined", curr_taxa$Taxon)
  
  # Then only keep selection
  curr_taxa = curr_taxa[grepl(paste(selected_taxa, collapse="|"), curr_taxa$Taxon),]
  # reorder by total abundance
  curr_taxa = cbind(curr_taxa, (rowSums(curr_taxa[1:(ncol(curr_taxa)-1)]))/ncol(curr_taxa) )
  curr_taxa = curr_taxa[order(-curr_taxa[, ncol(curr_taxa)]),,drop=FALSE]
  curr_taxa[,ncol(curr_taxa)] = NULL
  taxa_order = curr_taxa$Taxon
  
  curr_taxa = melt(curr_taxa)
  curr_mapping = mapping
  curr_taxa = merge(curr_taxa, curr_mapping, by.x="variable", by.y="#SampleID")
  
  # sort order of samples based on taxa order.
  
  tmp2 = dcast(curr_taxa, variable ~ Taxon, value.var="value")
  tmp2$variable = as.character(tmp2$variable)
  sorted_variables = as.character(tmp3$variable)
  curr_taxa$variable = factor(curr_taxa$variable, levels=sorted_variables)
  
  abundant_in_low_diversity = c(
    "p__Acidobacteria",
    "p__Proteobacteria;c__Alphaproteobacteria",
    "p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae",
    "p__Proteobacteria",
    "p__Acidobacteria;c__Acidobacteriia"
  )
  abundant_in_high_diversity = c("p__Actinobacteria;c__Actinobacteria")
  curr_taxa = curr_taxa[curr_taxa$Taxon %in% selected_taxa,]
  curr_taxa$high_low = ifelse(curr_taxa$Taxon %in% abundant_in_high_diversity, "Abundant in\nN fertilized fields\n (Nitrogen)", 
                              "Abundant in\nunfertilized fields\n (no Nitrogen)")
  curr_taxa$Taxon = factor(curr_taxa$Taxon, levels=unique(taxa_order))
  curr_taxa$value = curr_taxa$value * 100
  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_taxa$workflow = curr_workflow2
  curr_taxa$workflow2 = ifelse(grepl("arm", curr_taxa$workflow), "arm", "standard")
  curr_taxa$number_of_clusters = gsub("^(.*M).*", "\\1", curr_taxa$workflow)
  # populate filtered taxa df for further usage.
  if(i == 1){
    curr_taxa2 = curr_taxa
  }else{
    curr_taxa2 = rbind(curr_taxa2, curr_taxa)
  }
  
  p <- ggplot(data=curr_taxa, aes(x=Taxon, y=value, fill=Taxon)) + 
    facet_grid(high_low ~ Nitrogen, scales="free_y", space="free_x") + 
    geom_boxplot(aes(fill=Taxon), outlier.shape = NA, position=position_dodge()) + 
    geom_jitter(aes(color=Taxon),size=0.005, pch=19,position=position_jitterdodge(1.45), alpha=1) +
    xlab("") + 
    ylab("abundance (%)") +
    scale_fill_manual(values=color_list) +
    scale_color_manual(values=color_list) +
    guides(fill=guide_legend(ncol=2)) + ggtitle(curr_workflow2) +
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=8, colour="black"), 
      axis.title=element_text(family="Helvetica", size=11),
      plot.title = element_text(face="plain", size=14, hjust=0.5),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.07, 'cm'),
      legend.position="bottom",
      strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
      strip.background =  element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.3, "lines")
    ) 
  
  if(i == 3 | i == 6 | i == 9 | i == 12){ 
    p = p + theme(strip.text.y = element_text(angle=90, vjust=0, hjust=0.5, size=11, face="plain"))
  }else{
    p = p + theme(strip.text.y = element_text(size=0))
  }
  ps_ts[[i]] = p
  
}

color_list = c("arm" = "#0000CD", "standard" = "#FF0000")
taxa_order = selected_taxa

curr_taxa2 = curr_taxa2[grep("c__Acidobacteriia*", curr_taxa2$Taxon),]
curr_taxa2$number_of_clusters = factor(curr_taxa2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))

curr_taxa2$Taxon = gsub("p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae", 
                        "p__Acidobacteria;c__Acidobacteriia\n;o__Acidobacteriales;f__Acidobacteriaceae", curr_taxa2$Taxon)
# Try to do all at once:
p_selected_tax_soil <- ggplot(data=curr_taxa2, aes(x=number_of_clusters, y=value, fill=workflow2)) + 
  facet_grid(Taxon ~ Nitrogen, scales="free_y", space="free_x") + 
  geom_boxplot(aes(fill=workflow2), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow2),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  xlab("Number of clusters") + 
  ylab("abundance (%)") +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  guides(fill=guide_legend(ncol=2)) + ggtitle("") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title=element_text(family="Helvetica", size=11),
    plot.title = element_text(face="plain", size=14, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.07, 'cm'),
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
    strip.background =  element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(angle=0, vjust=0.5, hjust=0, size=10, face="plain")
  ) 
p_selected_tax_soil

#################################
# Mock communities              # 
#################################
mapping = data.frame(fread("../mock/mapping_file.tsv", sep="\t"), check.names=FALSE)

# Beta div
# To highlight the fact that level of diversity is maintained through the different 
# types of assembly config, split/bin/categorize the samples according to their 
# overall level of diversity
df_coords = get_coords_as_df(my_coords_files_mock)
df_alpha = get_alphadiv_as_df(my_alpha_div_genes_files_mock, mapping=mapping)
# include taxa profile stacked bar plots based on alpha diversity quintiles.
df_taxa = get_tax_as_df_unmelted(my_tax_files_mock_species, convert_to_rel=FALSE, mapping_file="../mock/mapping_file.tsv")

tmp = data.frame(fread("../mock/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv", header=T), check.names=FALSE)
tmp = cbind(tmp, (rowSums(tmp[2:(ncol(tmp)-1)]))/(ncol(tmp)-1) )
tmp = tmp[order(-tmp[, ncol(tmp)]),,drop=FALSE]
tmp[,ncol(tmp)] = NULL

selected_taxa = c(
  "p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus mutans",                   
  "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium beijerinckii",                  
  "p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Rhodobacter;s__Rhodobacter sphaeroides",  
  "p__Deinococcus-Thermus;c__Deinococci;o__Deinococcales;f__Deinococcaceae;g__Deinococcus;s__Deinococcus radiodurans",          
  "p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus gasseri",                  
  "p__Firmicutes;c__Bacilli;o__Bacillales;f__Listeriaceae;g__Listeria;s__Listeria monocytogenes",                               
  "p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Schaalia;s__Schaalia odontolytica",
  "k__undefined"
)

vColorsPcoa = c("#0000CD", "#FF0000","#808080","#000000","#FF8C00")
my_workflows = c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "0.1M_all", "0.5M_all", "1M_all", "2M_all", "3M_all", "4M_all", "complete")

ps_ts = list()
curr_taxa2 = NULL
# for simplicity sake, only show quartile 1 and 5
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_mapping = mapping
  
  # Taxa
  curr_taxa = df_taxa[df_taxa$Type == curr_workflow,]
  curr_taxa$Type = NULL
  # clean taxa lineages:
  curr_taxa$Taxon = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL", "k__Bacteria", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("p__NULL;c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("c__NULL;o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";o__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";c__NULL", ";c__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub("k__Bacteria;", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";p__NULL;", ";p__undefined", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";s__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";f__NULL;g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";g__NULL", "", curr_taxa$Taxon)
  curr_taxa$Taxon = gsub(";$", "", curr_taxa$Taxon)  
  curr_taxa$Taxon = gsub("^k__NULL$", "k__undefined", curr_taxa$Taxon)

  # Then only keep selection
  curr_taxa = curr_taxa[grepl(paste(selected_taxa, collapse="|"), curr_taxa$Taxon),]
  # reorder by total abundance
  curr_taxa = cbind(curr_taxa, (rowSums(curr_taxa[1:(ncol(curr_taxa)-1)]))/ncol(curr_taxa) )
  curr_taxa = curr_taxa[order(-curr_taxa[, ncol(curr_taxa)]),,drop=FALSE]
  curr_taxa[,ncol(curr_taxa)] = NULL
  taxa_order = curr_taxa$Taxon
  
  curr_taxa = melt(curr_taxa)
  curr_taxa = merge(curr_taxa, curr_mapping, by.x="variable", by.y="#SampleID")
  
  # sort order of samples based on taxa order.
  tmp2 = dcast(curr_taxa, variable ~ Taxon, value.var="value")
  tmp2$variable = as.character(tmp2$variable)
  sorted_variables = as.character(tmp3$variable)
  curr_taxa$variable = factor(curr_taxa$variable, levels=sorted_variables)

  
  
  curr_taxa$Taxon = factor(curr_taxa$Taxon, levels=unique(taxa_order))
  curr_taxa$value = curr_taxa$value * 100

  curr_workflow2 = gsub("_all", " - arm", curr_workflow)
  curr_taxa$workflow = curr_workflow2
  curr_taxa$workflow2 = ifelse(grepl("arm", curr_taxa$workflow), "arm", "standard")
  curr_taxa$number_of_clusters = gsub("^(.*M).*", "\\1", curr_taxa$workflow)
  # populate filtered taxa df for further usage.
  if(i == 1){
    curr_taxa2 = curr_taxa
  }else{
    curr_taxa2 = rbind(curr_taxa2, curr_taxa)
  }
  
  p <- ggplot(data=curr_taxa, aes(x=Taxon, y=value, fill=Taxon)) + 
    facet_grid(high_low ~ Nitrogen, scales="free_y", space="free_x") + 
    geom_boxplot(aes(fill=Taxon), outlier.shape = NA, position=position_dodge()) + 
    geom_jitter(aes(color=Taxon),size=0.005, pch=19,position=position_jitterdodge(1.45), alpha=1) +
    xlab("") + 
    ylab("abundance (%)") +
    scale_fill_manual(values=color_list) +
    scale_color_manual(values=color_list) +
    guides(fill=guide_legend(ncol=2)) + ggtitle(curr_workflow2) +
    theme(
      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=10, colour="black"), 
      axis.title=element_text(family="Helvetica", size=11),
      plot.title = element_text(face="plain", size=14, hjust=0.5),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      legend.key.size = unit(0.45, "cm"),
      legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.07, 'cm'),
      legend.position="bottom",
      strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
      strip.background =  element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.3, "lines")
    ) 
  
  if(i == 3 | i == 6 | i == 9 | i == 12){ 
    p = p + theme(strip.text.y = element_text(angle=90, vjust=0, hjust=0.5, size=11, face="plain"))
  }else{
    p = p + theme(strip.text.y = element_text(size=0))
  }
  ps_ts[[i]] = p
  
}

color_list = c("arm" = "#0000CD", "standard" = "#FF0000")
taxa_order = selected_taxa
selected_taxa2 = c(
  "p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus mutans",                   
  "p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium beijerinckii",                  
  "p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Rhodobacter;s__Rhodobacter sphaeroides",  
  "p__Deinococcus-Thermus;c__Deinococci;o__Deinococcales;f__Deinococcaceae;g__Deinococcus;s__Deinococcus radiodurans",          
  "p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus gasseri",                  
  "p__Firmicutes;c__Bacilli;o__Bacillales;f__Listeriaceae;g__Listeria;s__Listeria monocytogenes",                               
  "p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Schaalia;s__Schaalia odontolytica"
)

curr_taxa2 = curr_taxa2[curr_taxa2$Taxon %in% selected_taxa2,]
curr_taxa2$number_of_clusters = factor(curr_taxa2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "2M", "3M", "4M", "complete"))

curr_taxa2$Taxon = gsub("p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus mutans", 
                        "p__Firmicutes;c__Bacilli;\no__Lactobacillales;f__Streptococcaceae;\ng__Streptococcus;s__Streptococcus mutans", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium beijerinckii","p__Firmicutes;c__Clostridia;\no__Clostridiales;f__Clostridiaceae;\ng__Clostridium;s__Clostridium beijerinckii", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Rhodobacter;s__Rhodobacter sphaeroides", 
                        "p__Proteobacteria;c__Alphaproteobacteria;\no__Rhodobacterales;f__Rhodobacteraceae;\ng__Rhodobacter;s__Rhodobacter sphaeroides", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Deinococcus-Thermus;c__Deinococci;o__Deinococcales;f__Deinococcaceae;g__Deinococcus;s__Deinococcus radiodurans", 
                        "p__Deinococcus-Thermus;c__Deinococci;\no__Deinococcales;f__Deinococcaceae;\ng__Deinococcus;s__Deinococcus radiodurans", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus gasseri", 
                        "p__Firmicutes;c__Bacilli;\no__Lactobacillales;f__Lactobacillaceae;\ng__Lactobacillus;s__Lactobacillus gasseri", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Firmicutes;c__Bacilli;o__Bacillales;f__Listeriaceae;g__Listeria;s__Listeria monocytogenes", 
                        "p__Firmicutes;c__Bacilli\n;o__Bacillales;f__Listeriaceae;\ng__Listeria;s__Listeria monocytogenes", curr_taxa2$Taxon)
curr_taxa2$Taxon = gsub("p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Schaalia;s__Schaalia odontolytica", 
                        "p__Actinobacteria;c__Actinobacteria;\no__Actinomycetales;f__Actinomycetaceae;\ng__Schaalia;s__Schaalia odontolytica", curr_taxa2$Taxon)
# Try to do all at once:

# Try to do all at once:
p_selected_tax_mock <- ggplot(data=curr_taxa2, aes(x=number_of_clusters, y=value, fill=workflow2)) + 
  facet_grid(Taxon ~ Treatment, scales="free_y", space="free_x") + 
  geom_boxplot(aes(fill=workflow2, color=workflow2), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow2),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  xlab("Number of clusters") + 
  ylab("abundance (%)") +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  guides(fill=guide_legend(ncol=2)) + ggtitle("") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title=element_text(family="Helvetica", size=11),
    plot.title = element_text(face="plain", size=14, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.07, 'cm'),
    legend.position="bottom",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=10, face="plain"),
    strip.background =  element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(angle=0, vjust=0.5, hjust=0, size=10, face="plain")
  ) 
p_selected_tax_mock

## Finally,
# Build figure
# Finally, generate figures.
fig_soil = ggarrange(
  p_selected_tax_soil + rremove("legend"), 
  nrow = 2, ncol= 1, heights=c(1,1,1)
)
f9a = annotate_figure(
  ggarrange(
    p_selected_tax_china + rremove("legend"), 
    p_selected_tax_rock + rremove("legend"),
    fig_soil, 
    p_selected_tax_mock,
    labels=c("A", "B", "C", "D"),
    nrow = 2, ncol= 2
  ),
  top = text_grob("Figure 4\n", color = "black", face = "plain", size = 11),
)

pdf( file=paste0(outdir, "./figures/F_FIGURE_4.pdf"), height=14, width=14)
f9a
dev.off()


#################################
# Finally find KO that are sda  #
#################################
#df_KO = get_KO_abundance_as_df(my_KO_files_china, my_gene_abundance_files_china, my_COG_files_china, mapping_file)
df_KO = data.frame(fread(paste0(outdir, "./df_KO_china.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO) = df_KO$V1; df_KO$V1=NULL;
head(df_KO)
mapping = data.frame(fread("./mapping_file.tsv"), check.names=FALSE)

my_workflows = unique(df_KO$Type)
ps_ko = list()
anova_df3 = NULL
# for simplicity sake, only show quartile 1 and 5 
for(i in 1:length(my_workflows)){
  print(paste0("i: ", i))
  curr_workflow = my_workflows[i]
  curr_KO = df_KO[df_KO$Type == curr_workflow,]
  KOs = unique(curr_KO$KO)
  KOs = KOs[KOs != "NULL"]
  curr_ad = df_alpha[df_alpha$Type == curr_workflow,]
  
  # Divide by quintiles:
  curr_ad$quantiles = NULL
  curr_ad$quantiles = cut(curr_ad$value, 
                          quantile(
                            curr_ad$value, 
                            prob = seq(0, 1, length = 6), 
                            type = 5, names=TRUE, 
                            labels=FALSE, 
                            ordered_result=TRUE
                          ), include.lowest=TRUE
  )
  
  df_link = data.frame(quantiles=levels(curr_ad$quantiles), rank=seq(1,(length(unique(curr_ad$quantiles))), by=1))
  df_tmp = curr_ad[,c("variable", "quantiles")]
  df_tmp = merge(df_tmp, df_link, by="quantiles")
  curr_mapping = merge(mapping, df_tmp, by.x="#SampleID", by.y="variable")
  colnames(curr_mapping)[ncol(curr_mapping)] = "Diversity_quartile"
  
  curr_KO = melt(curr_KO)
  curr_KO = merge(curr_KO, curr_mapping, by.x="variable", by.y="#SampleID")
  
  anova_df = NULL
  anova_df2 = NULL
  for(j in 1:length(KOs)){
    print(paste0("   j: ", j))
    curr_KO2 = curr_KO[curr_KO$KO == KOs[j],]
    if(nrow(curr_KO2) == 0){next;}
    res = aov(value ~ Diversity_quartile, data=curr_KO2)
    pvalue = summary(res)[[1]][["Pr(>F)"]][1]
    
    mean_1 = mean(curr_KO2[curr_KO2$Diversity_quartile == "1",]$value)
    mean_5 = mean(curr_KO2[curr_KO2$Diversity_quartile == "5",]$value)
    
    if(j == 1){
      anova_df = data.frame(KO=KOs[j], pvalue=pvalue, mean_1=mean_1, mean_5=mean_5)
    }else{
      anova_df = rbind(anova_df, data.frame(KO=KOs[j], pvalue=pvalue, mean_1=mean_1, mean_5=mean_5))
    }
    
  }
  anova_df$BH = p.adjust(anova_df$pvalue, method = "bonferroni")
  anova_df2 = anova_df[anova_df$BH < 100 & abs(anova_df$mean_1/anova_df$mean_5) >= 0,] # use lenient value, but use stringent later (say fc 11)
  anova_df2 = na.omit(anova_df2)
  anova_df2$Type = curr_workflow
  
  if(i == 1){
    anova_df3 = anova_df2
  }else{
    anova_df3 = rbind(anova_df3, anova_df2)
  }
}

# Then filter results from the raw table (also write table to avoid regenerating anovas)
#write.table(anova_df3, paste0(outdir, "/anova_df3_china.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
anova_df3 = data.frame(fread(paste0(outdir, "/anova_df3_china.tsv"), header=TRUE, sep="\t"), check.names=FALSE)

anova_df4 = anova_df3[anova_df3$BH < 0.05 & anova_df3$mean_1/anova_df3$mean_5 >= 5 | anova_df3$mean_5/anova_df3$mean_1 >= 5,]

# Add KOs info
keggs = data.frame(fread("~/databases/kegg_pathways_modules.tsv", header=TRUE), check.names=FALSE)
anova_df4 = merge(anova_df4, keggs[,c("KO", "KO_desc", "plevel2", "ko_desc")], by.x="KO", by.y="KO", all.x=TRUE, all.y=FALSE)

# Now that we have the interesting KOs in hand, lets generate a heatmap
df_KO2 = df_KO
df_KO2$KO_Type = paste0(df_KO2$KO, "_", df_KO2$Type)
anova_df4$KO_Type = paste0(anova_df4$KO, "_", anova_df4$Type)
df_KO2 = merge(anova_df4, df_KO2, by.x="KO_Type", by.y="KO_Type", all.x=FALSE, all.y=FALSE)
colnames(df_KO2)[2] = "KO"

# normalize by recA
## Then compute recA abundance
## get recA reads.  Here I follow the procedure described in acinas et al., biology comm.
## So basically, first sum all recA reads for each sample (i.e. column).
## Then divide the sum of each cell against their recA summed number.
## But then, row-wise, also divide each cell against the max of their row.
## The reason to normalize data against each KO max is that some KO are way more abundant than others
## , so the signal gets lost to the profit of the high abundant ones.
## Already normalized
col_order = unique(mapping$`#SampleID`)
abundance_KO = df_KO2[,c(col_order, "KO_Type")]
row.names(abundance_KO) = abundance_KO$KO_Type; abundance_KO$KO_Type = NULL;

############################
#row_max = as.numeric(apply(abundance_KO, 1, max))
#abundance_KO = abundance_KO / row_max
sample_order = row.names(abundance_KO)
row.names(abundance_KO) = gsub("_all", " - arm", row.names(abundance_KO))
uncombined = data.frame(do.call(rbind, strsplit(row.names(abundance_KO), "_", fixed=TRUE) ) )#[,c(2)])
colnames(uncombined) = c("KO", "workflow")
abundance_KO = cbind(abundance_KO, uncombined)
abundance_KO = abundance_KO[abundance_KO$KO %in% c("K22607","K01593","K11444","K18148","K22225",  "K16149", "K01178"),]
# melt + ggplot geom_tile
abundance_KO2 = melt(abundance_KO)
mapping4 = data.frame(fread("./mapping_file4.tsv", header=TRUE), check.names=FALSE)
sample_list = colnames(df_KO)[1:912]

# mapping file 4 has the diversity quintiles for 12M clusters workflow
abundance_KO2 = merge(abundance_KO2, mapping4[c("#SampleID", "Diversity_quartile")], by.x="variable", by.y="#SampleID")
tmp = abundance_KO2[abundance_KO2$workflow == "12M",]

abundance_KO2$workflow = factor(abundance_KO2$workflow, levels=c(
  "0.1M", "0.1M - arm", "0.5M", "0.5M - arm", "1M", "1M - arm", "4M",  "4M - arm",
  "8M", "8M - arm", "12M", "12M - arm"
))

KO_order = unique(abundance_KO2$KO)
abundance_KO2$KO = factor(abundance_KO2$KO, levels=KO_order)
keggs = data.frame(fread("~/databases/kegg_pathways_modules.tsv", header=TRUE), check.names=FALSE)
keggs$plevel2 = ifelse(keggs$KO %in% c("K22607","K01593","K11444","K18148","K22225",  "K16149", "K01178"), keggs$KO_desc, keggs$plevel2) 
abundance_KO2 = merge(abundance_KO2, keggs[,c("KO", "KO_desc", "plevel2", "ko_desc")], by.x="KO", by.y="KO", all.x=TRUE, all.y=FALSE)
abundance_KO2$KO_plevel2 = paste0(abundance_KO2$KO, " - ", abundance_KO2$plevel2)
unique(abundance_KO2$KO_plevel2[grepl(paste(KO_order, collapse="|"), unique(abundance_KO2$KO_plevel2))])
abundance_KO2[grepl("K01593", abundance_KO2$KO_plevel2),]$KO_plevel2 = "K01593 - Biosynthesis of other secondary metabolites"

KO_plevel2_order = unique(unlist(lapply(KO_order, grep,abundance_KO2$KO_plevel2, value = TRUE)))
abundance_KO2$KO_plevel2 = factor(abundance_KO2$KO_plevel2, levels=KO_plevel2_order)

# Genreate barplot
M01 = unique(as.character((abundance_KO2[abundance_KO2$workflow == "0.1M",]$KO)))
M01 = M01[!is.na(M01)]
M05 = unique(as.character((abundance_KO2[abundance_KO2$workflow == "0.5M",]$KO)))
M05 = M05[!is.na(M05)]
M1 = unique(as.character((abundance_KO2[abundance_KO2$workflow == "1M",]$KO)))
M1 = M1[!is.na(M1)]
M4 = unique(as.character((abundance_KO2[abundance_KO2$workflow == "4M",]$KO)))
M4 = M4[!is.na(M4)]
M8 = unique(as.character((abundance_KO2[abundance_KO2$workflow == "8M",]$KO)))
M8 = M8[!is.na(M8)]
M12 = unique(as.character((abundance_KO2[abundance_KO2$workflow == "12M",]$KO)))
M12 = M12[!is.na(M12)]

M01_arm = unique(as.character((abundance_KO2[abundance_KO2$workflow == "0.1M - arm",]$KO)))
M01_arm = M01[!is.na(M01_arm)]
M05_arm = unique(as.character((abundance_KO2[abundance_KO2$workflow == "0.5M - arm",]$KO)))
M05_arm = M05_arm[!is.na(M05_arm)]
M1_arm = unique(as.character((abundance_KO2[abundance_KO2$workflow == "1M - arm",]$KO)))
M1_arm = M1_arm[!is.na(M1_arm)]
M4_arm = unique(as.character((abundance_KO2[abundance_KO2$workflow == "4M - arm",]$KO)))
M4_arm = M4_arm[!is.na(M4_arm)]
M8_arm = unique(as.character((abundance_KO2[abundance_KO2$workflow == "8M - arm",]$KO)))
M8_arm = M8_arm[!is.na(M8)]
M12_arm = unique(as.character((abundance_KO2[abundance_KO2$workflow == "12M - arm",]$KO)))
M12_arm = M12_arm[!is.na(M12)]
all_vecs=NULL
all_vecs = list("0.1M"=M01,
                "0.1M_arm"=M01_arm, 
                "0.5M"=M05,
                "0.5M_arm"= M05_arm,
                "1M"=M1,
                "1M_arm"=M1_arm,
                "4M"=M4,
                "4M_arm"=M4_arm,
                "8M"=M8,
                "8M_arm"=M8_arm, 
                "12M"=M12,
                "12M_arm"=M12_arm)

data_bp = NULL
for(i in 1:length(all_vecs)){
  curr_vec = all_vecs[[i]]
  curr_name = names(all_vecs)[i]
  
  if(i == 1){
    data_bp = data.frame(workflow = curr_name, number_of_unique_KO = length(curr_vec))
  }else{
    data_bp = rbind(data_bp, data.frame(workflow = curr_name, number_of_unique_KO = length(curr_vec)))
  }
  
}  
data_bp
data_bp$workflow = gsub("_arm", " - arm", data_bp$workflow)
data_bp$number_of_clusters = data_bp$workflow
data_bp$number_of_clusters = gsub(" - arm", "", data_bp$number_of_clusters)
data_bp$workflow2 = ifelse(grepl(" - arm", data_bp$workflow), "arm", "normal")
data_bp$number_of_clusters = factor(data_bp$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))


p_bp = ggplot(data=data_bp, aes(x=number_of_clusters, y=number_of_unique_KO, fill=workflow2)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  xlab("Number of clusters") + ylab("Unique KOs recovered") +
  guides(fill=guide_legend("")) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.5),
    axis.text.x=element_text(size=9, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=9, colour="black"), 
    axis.title.y=element_text(size=9),
    axis.title.x=element_text(size=9),
    plot.title = element_text(lineheight=1.2, face="bold", size=11),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_text(size=9, face="plain"),
    legend.position="right",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=6, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=6, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
p_bp

# Then write plot for 12 M
abundance_KO2$KO_plevel2 = gsub("K22225 - Metabolism", "K22225 - Metabolism of cofactors and vitamins", abundance_KO2$KO_plevel2)
selected_KOs = c(
  "K22607","K01593","K11444","K18148","K22225",  "K16149", "K01178"
)
abundance_KO3 = abundance_KO2[abundance_KO2$KO %in% selected_KOs,]
abundance_KO3 = abundance_KO3[abundance_KO3$Diversity_quartile %in% c(1,5),]

abundance_KO3$workflow2 = gsub("_all", " - arm", abundance_KO3$workflow)
abundance_KO3$workflow3 = ifelse(grepl("arm", abundance_KO3$workflow), "arm", "standard")
abundance_KO3$number_of_clusters = gsub("^(.*M).*", "\\1", abundance_KO3$workflow)
abundance_KO3$Diversity_quartile = as.character(abundance_KO3$Diversity_quartile)
abundance_KO3$Diversity_quartile = gsub("1", "Quintile #1\n(Participants harboring\n low diversity microbiota)", abundance_KO3$Diversity_quartile)
abundance_KO3$Diversity_quartile = gsub("5", "Quintile #5\n(Participants harboring\n high diversity microbiota)", abundance_KO3$Diversity_quartile)
abundance_KO3$number_of_clusters = factor(abundance_KO3$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M"))
abundance_KO3$KO_plevel2 = factor(abundance_KO3$KO_plevel2, levels=unique(abundance_KO3$KO_plevel2)[c(2,3,5,6,7,8,1,4)])
abundance_KO3$KO_plevel2 = gsub("K11444 - wspR; two-component system, chemotaxis family, response regulator WspR",
                                "K11444 - wspR; two-component system, chemotaxis family,\n response regulator WspR", abundance_KO3$KO_plevel2)

color_list = c("arm" = "#0000CD", "standard" = "#FF0000")
p_KO_china_3 = ggplot(abundance_KO3, aes(y=value, x=number_of_clusters, fill=workflow3)) +
  geom_boxplot(aes(fill=workflow3), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow3),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  facet_grid(KO_plevel2 ~ Diversity_quartile, scales="free", space="free_x") +
  xlab("Number of clusters") + ylab("reads/recA reads") + 
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.5),
    axis.text.x=element_text(size=7, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=7, colour="black"), 
    axis.title.y=element_text(size=8),
    axis.title.x=element_text(size=8),
    axis.ticks.x = element_blank(),
    plot.title = element_text(lineheight=1.2, face="bold", size=11),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_blank(),
    legend.position="top",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
p_KO_china_3

###############################
# Do KO for Antarctic         #
###############################
#df_KO = get_KO_abundance_as_df(my_KO_files_china, my_gene_abundance_files_china, my_COG_files_china, mapping_file)
df_KO = data.frame(fread(paste0(outdir, "./df_KO_rock.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO) = df_KO$V1; df_KO$V1=NULL;
head(df_KO)
mapping = data.frame(fread("../rock_microbiome/0.1M_clusters/mapping_file.tsv"), check.names=FALSE)
mapping$Location = gsub(".nord", "", mapping$Treatment)
mapping$Location = gsub(".sud", "", mapping$Location)
mapping = mapping[mapping$Location %in% c("Siegfried.Peak", "University.Valley"),]
samples = unique(mapping$`#SampleID`)
curr_mapping = mapping

my_workflows = unique(df_KO$Type)
df_KO = df_KO[,c(samples, "KO", "Type")]
ps_ko = list()
anova_df3 = NULL
# for simplicity sake, only show quartile 1 and 5 
for(i in 1:length(my_workflows)){
  print(paste0("i: ", i))
  curr_workflow = my_workflows[i]
  curr_KO = df_KO[df_KO$Type == curr_workflow,]
  KOs = unique(curr_KO$KO)
  KOs = KOs[KOs != "NULL"]
  
  curr_KO = melt(curr_KO)
  curr_KO = merge(curr_KO, curr_mapping, by.x="variable", by.y="#SampleID")
  
  anova_df = NULL
  anova_df2 = NULL
  for(j in 1:length(KOs)){
    print(paste0("   j: ", j))
    curr_KO2 = curr_KO[curr_KO$KO == KOs[j],]
    if(nrow(curr_KO2) == 0){next;}
    res = aov(value ~ Location, data=curr_KO2)
    pvalue = summary(res)[[1]][["Pr(>F)"]][1]
    
    anova_1 = mean(curr_KO2[curr_KO2$Location == "Siegfried.Peak",]$value)
    anova_2 = mean(curr_KO2[curr_KO2$Location == "University.Valley",]$value)
    
    if(j == 1){
      anova_df = data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2)
    }else{
      anova_df = rbind(anova_df, data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2))
    }
    
  }
  anova_df$BH = p.adjust(anova_df$pvalue, method = "bonferroni")
  anova_df2 = anova_df[anova_df$pvalue < 0.05 & abs(anova_df$anova_1/anova_df$anova_2) >= 0,] # use lenient value, but use stringent later (say fc 11)
  anova_df2 = na.omit(anova_df2)
  anova_df2$Type = curr_workflow
  
  if(i == 1){
    anova_df3 = anova_df2
  }else{
    anova_df3 = rbind(anova_df3, anova_df2)
  }
 
}

# Then filter results from the raw table (also write table to avoid regenerating anovas)
#write.table(anova_df3, paste0(outdir, "/anova_df3_rock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
anova_df3 = data.frame(fread(paste0(outdir, "/anova_df3_rock.tsv"), header=TRUE, sep="\t"), check.names=FALSE)

anova_df4 = anova_df3[anova_df3$pvalue < 0.05 & abs(anova_df3$anova_1/anova_df3$anova_2) >= 2,]
mapping2 = mapping
row.names(mapping2) = mapping2$`#SampleID`;  mapping2$`#SampleID` = NULL;
## exploratory hm first
df_KO3 = df_KO[df_KO$KO %in% unique(anova_df4$KO),]
row.names(df_KO3) = paste0(df_KO3$KO, " - ", df_KO3$Type); df_KO3$KO = NULL; df_KO3$Type = NULL;
pheatmap::pheatmap(
  df_KO3, 
  file=paste0(outdir, "figures/tax_KO_test_rock_heatmap.pdf"),
  fontsize_row=9, 
  fontsize_col=3, 
  cellwidth=5,  
  cellheight=10,
  annotation=mapping2, 
  color=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100),
  clustering_method="average"
)
selected_KOs = c("K01915", "K01955", "K03556", "K07487", "K00341")
# Now that we have the interesting KOs in hand, lets generate a heatmap
df_KO2 = df_KO
df_KO2$KO_Type = paste0(df_KO2$KO, "_", df_KO2$Type)
df_KO2 = df_KO2[df_KO2$KO %in% selected_KOs,]
## Already normalized

# melt + ggplot geom_tile
abundance_KO2 = melt(df_KO2)
abundance_KO2 = merge(abundance_KO2, mapping[c("#SampleID", "Location")], by.x="variable", by.y="#SampleID")
abundance_KO2$workflow = gsub("_all", " - arm", abundance_KO2$Type)
abundance_KO2$workflow = factor(abundance_KO2$workflow, levels=c(
  "0.1M", "0.1M - arm", "0.5M", "0.5M - arm", "1M", "1M - arm", "4M",  "4M - arm",
  "8M", "8M - arm", "12M", "12M - arm", "complete"
))
KO_order = unique(abundance_KO2$KO)
abundance_KO2$KO = factor(abundance_KO2$KO, levels=KO_order)
keggs = data.frame(fread("~/databases/kegg_pathways_modules.tsv", header=TRUE), check.names=FALSE)
keggs$plevel2 = ifelse(keggs$KO %in% selected_KOs, keggs$KO_desc, keggs$plevel2) 
abundance_KO2 = merge(abundance_KO2, keggs[,c("KO", "KO_desc", "plevel2", "ko_desc")], by.x="KO", by.y="KO", all.x=TRUE, all.y=FALSE)
abundance_KO2$KO_plevel2 = paste0(abundance_KO2$KO, " - ", abundance_KO2$plevel2)
unique(abundance_KO2$KO_plevel2[grepl(paste(KO_order, collapse="|"), unique(abundance_KO2$KO_plevel2))])
KO_plevel2_order = unique(unlist(lapply(KO_order, grep,abundance_KO2$KO_plevel2, value = TRUE)))
abundance_KO2$KO_plevel2 = factor(abundance_KO2$KO_plevel2, levels=KO_plevel2_order)
abundance_KO2$workflow2 = gsub("_all", " - arm", abundance_KO2$workflow)
abundance_KO2$workflow3 = ifelse(grepl("arm", abundance_KO2$workflow2), "arm", "standard")
abundance_KO2$number_of_clusters = gsub("^(.*M).*", "\\1", abundance_KO2$workflow)
abundance_KO2$number_of_clusters = factor(abundance_KO2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))
abundance_KO2$Location  = gsub("\\.", " ", abundance_KO2$Location)
abundance_KO2$KO_plevel2 = gsub("malT; LuxR family transcriptional regulator, maltose regulon positive regulatory protein",
                                "malT; LuxR family transcriptional regulator, \nmaltose regulon positive regulatory protein", abundance_KO2$KO_plevel2)

p_KO_rock_3 = ggplot(abundance_KO2, aes(y=value, x=number_of_clusters, fill=workflow3)) +
  geom_boxplot(aes(fill=workflow3), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow3),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  facet_grid(KO_plevel2 ~ Location, scales="free", space="free_x") +
  xlab("Number of clusters") + ylab("reads/recA reads") + 
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.5),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title.y=element_text(size=9),
    axis.title.x=element_text(size=9),
    axis.ticks.x = element_blank(),
    plot.title = element_text(lineheight=1.2, face="bold", size=11),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_blank(),
    legend.position="top",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
p_KO_rock_3

########################
# Agricultural soil KO #
########################
mapping_file = "../soil/0.1M_clusters/mapping_file.tsv"
#df_KO = get_KO_abundance_as_df(my_KO_files_soil, my_gene_abundance_files_soil, my_COG_files_soil, mapping_file)
#fwrite(df_KO, paste0(outdir, "./df_KO_soil.tsv"), sep="\t", row.names=TRUE)
df_KO = data.frame(fread(paste0(outdir, "./df_KO_rock.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO) = df_KO$V1; df_KO$V1=NULL;
head(df_KO)
mapping = data.frame(fread("../soil/0.1M_clusters/mapping_file.tsv"), check.names=FALSE)
#mapping$Location = gsub(".nord", "", mapping$Treatment)
#mapping$Location = gsub(".sud", "", mapping$Location)
#mapping = mapping[mapping$Location %in% c("Siegfried.Peak", "University.Valley"),]
samples = unique(mapping$`#SampleID`)
colnames(df_KO) = gsub("\\.", "-", colnames(df_KO))
curr_mapping = mapping

my_workflows = unique(df_KO$Type)
df_KO = df_KO[,c(samples, "KO", "Type")]
ps_ko = list()
anova_df3 = NULL
# for simplicity sake, only show quartile 1 and 5 
for(i in 1:length(my_workflows)){
  print(paste0("i: ", i))
  curr_workflow = my_workflows[i]
  curr_KO = df_KO[df_KO$Type == curr_workflow,]
  KOs = unique(curr_KO$KO)
  KOs = KOs[KOs != "NULL"]
  
  curr_KO = melt(curr_KO)
  curr_KO = merge(curr_KO, curr_mapping, by.x="variable", by.y="#SampleID")
  
  anova_df = NULL
  anova_df2 = NULL
  for(j in 1:length(KOs)){
    print(paste0("   j: ", j))
    curr_KO2 = curr_KO[curr_KO$KO == KOs[j],]
    if(nrow(curr_KO2) == 0){next;}
    res = aov(value ~ Nitrogen, data=curr_KO2)
    pvalue = summary(res)[[1]][["Pr(>F)"]][1]
    
    anova_1 = mean(curr_KO2[curr_KO2$Nitrogen == "N",]$value)
    anova_2 = mean(curr_KO2[curr_KO2$Nitrogen == "noN",]$value)
    
    if(j == 1){
      anova_df = data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2)
    }else{
      anova_df = rbind(anova_df, data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2))
    }
    
  }
  anova_df$BH = p.adjust(anova_df$pvalue, method = "bonferroni")
  anova_df2 = anova_df[anova_df$pvalue < 0.05 & abs(anova_df$anova_1/anova_df$anova_2) >= 0,] # use lenient value, but use stringent later (say fc 11)
  anova_df2 = na.omit(anova_df2)
  anova_df2$Type = curr_workflow
  
  if(i == 1){
    anova_df3 = anova_df2
  }else{
    anova_df3 = rbind(anova_df3, anova_df2)
  }
}

# Then filter results from the raw table (also write table to avoid regenerating anovas)
#write.table(anova_df3, paste0(outdir, "/anova_df3_soil.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
anova_df3 = data.frame(fread(paste0(outdir, "/anova_df3_soil.tsv"), header=TRUE, sep="\t"), check.names=FALSE)

anova_df4 = anova_df3[anova_df3$BH < 0.05 & abs(anova_df3$anova_1/anova_df3$anova_2) >= 20,]
mapping2 = mapping
row.names(mapping2) = mapping2$`#SampleID`;  mapping2$`#SampleID` = NULL;
## exploratory hm first
df_KO3 = df_KO[df_KO$KO %in% unique(anova_df4$KO),]
row.names(df_KO3) = paste0(df_KO3$KO, " - ", df_KO3$Type); df_KO3$KO = NULL; df_KO3$Type = NULL;
pheatmap::pheatmap(
  df_KO3, 
  file=paste0(outdir, "figures/tax_KO_test_rock_heatmap.pdf"),
  fontsize_row=9, 
  fontsize_col=3, 
  cellwidth=5,  
  cellheight=10,
  annotation=mapping2, 
  color=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100),
  clustering_method="average"
)
selected_KOs = c("K09001", "K22020", "K09001", "K00363", "K18829")
# Now that we have the interesting KOs in hand, lets generate a heatmap
df_KO2 = df_KO
df_KO2$KO_Type = paste0(df_KO2$KO, "_", df_KO2$Type)
df_KO2 = df_KO2[df_KO2$KO %in% selected_KOs,]
## Already normalized

# melt + ggplot geom_tile
abundance_KO2 = melt(df_KO2)
abundance_KO2 = merge(abundance_KO2, mapping[c("#SampleID", "Nitrogen")], by.x="variable", by.y="#SampleID")
abundance_KO2$workflow = gsub("_all", " - arm", abundance_KO2$Type)
abundance_KO2$workflow = factor(abundance_KO2$workflow, levels=c(
  "0.1M", "0.1M - arm", "0.5M", "0.5M - arm", "1M", "1M - arm", "4M",  "4M - arm",
  "8M", "8M - arm", "12M", "12M - arm", "complete"
))
KO_order = unique(abundance_KO2$KO)
abundance_KO2$KO = factor(abundance_KO2$KO, levels=KO_order)
keggs = data.frame(fread("~/databases/kegg_pathways_modules.tsv", header=TRUE), check.names=FALSE)
keggs$plevel2 = ifelse(keggs$KO %in% selected_KOs, keggs$KO_desc, keggs$plevel2) 
abundance_KO2 = merge(abundance_KO2, keggs[,c("KO", "KO_desc", "plevel2", "ko_desc")], by.x="KO", by.y="KO", all.x=TRUE, all.y=FALSE)
abundance_KO2$KO_plevel2 = paste0(abundance_KO2$KO, " - ", abundance_KO2$plevel2)
unique(abundance_KO2$KO_plevel2[grepl(paste(KO_order, collapse="|"), unique(abundance_KO2$KO_plevel2))])
KO_plevel2_order = unique(unlist(lapply(KO_order, grep,abundance_KO2$KO_plevel2, value = TRUE)))
abundance_KO2$KO_plevel2 = factor(abundance_KO2$KO_plevel2, levels=KO_plevel2_order)
abundance_KO2$workflow2 = gsub("_all", " - arm", abundance_KO2$workflow)
abundance_KO2$workflow3 = ifelse(grepl("arm", abundance_KO2$workflow2), "arm", "standard")
abundance_KO2$number_of_clusters = gsub("^(.*M).*", "\\1", abundance_KO2$workflow)
abundance_KO2$number_of_clusters = factor(abundance_KO2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))


p_KO_soil_3 = ggplot(abundance_KO2, aes(y=value, x=number_of_clusters, fill=workflow3)) +
  geom_boxplot(aes(fill=workflow3), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow3),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  facet_grid(KO_plevel2 ~ Nitrogen, scales="free", space="free_x") +
  xlab("Number of clusters") + ylab("reads/recA reads") + 
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.5),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title.y=element_text(size=9),
    axis.title.x=element_text(size=9),
    axis.ticks.x = element_blank(),
    plot.title = element_text(lineheight=1.2, face="bold", size=11),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_blank(),
    legend.position="top",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
p_KO_soil_3

########################
# Agricultural soil KO #
########################
mapping_file = "../soil/0.1M_clusters/mapping_file.tsv"
#df_KO = get_KO_abundance_as_df(my_KO_files_soil, my_gene_abundance_files_soil, my_COG_files_soil, mapping_file)
#fwrite(df_KO, paste0(outdir, "./df_KO_soil.tsv"), sep="\t", row.names=TRUE)
df_KO = data.frame(fread(paste0(outdir, "./df_KO_soil.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO) = df_KO$V1; df_KO$V1=NULL;
head(df_KO)
mapping = data.frame(fread("../soil/0.1M_clusters/mapping_file.tsv"), check.names=FALSE)
samples = unique(mapping$`#SampleID`)
colnames(df_KO) = gsub("\\.", "-", colnames(df_KO))
curr_mapping = mapping

my_workflows = unique(df_KO$Type)
df_KO = df_KO[,c(samples, "KO", "Type")]
ps_ko = list()
anova_df3 = NULL
# for simplicity sake, only show quartile 1 and 5 
for(i in 1:length(my_workflows)){
  print(paste0("i: ", i))
  curr_workflow = my_workflows[i]
  curr_KO = df_KO[df_KO$Type == curr_workflow,]
  KOs = unique(curr_KO$KO)
  KOs = KOs[KOs != "NULL"]
  
  curr_KO = melt(curr_KO)
  curr_KO = merge(curr_KO, curr_mapping, by.x="variable", by.y="#SampleID")
  
  anova_df = NULL
  anova_df2 = NULL
  for(j in 1:length(KOs)){
    print(paste0("   j: ", j))
    curr_KO2 = curr_KO[curr_KO$KO == KOs[j],]
    if(nrow(curr_KO2) == 0){next;}
    res = aov(value ~ Nitrogen, data=curr_KO2)
    pvalue = summary(res)[[1]][["Pr(>F)"]][1]
    
    anova_1 = mean(curr_KO2[curr_KO2$Nitrogen == "N",]$value)
    anova_2 = mean(curr_KO2[curr_KO2$Nitrogen == "noN",]$value)
    
    if(j == 1){
      anova_df = data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2)
    }else{
      anova_df = rbind(anova_df, data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2))
    }
    
  }
  anova_df$BH = p.adjust(anova_df$pvalue, method = "bonferroni")
  anova_df2 = anova_df[anova_df$pvalue < 0.05 & abs(anova_df$anova_1/anova_df$anova_2) >= 0,] # use lenient value, but use stringent later (say fc 11)
  anova_df2 = na.omit(anova_df2)
  anova_df2$Type = curr_workflow
  
  if(i == 1){
    anova_df3 = anova_df2
  }else{
    anova_df3 = rbind(anova_df3, anova_df2)
  }
  
}

# Then filter results from the raw table (also write table to avoid regenerating anovas)
#write.table(anova_df3, paste0(outdir, "/anova_df3_soil.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
anova_df3 = data.frame(fread(paste0(outdir, "/anova_df3_soil.tsv"), header=TRUE, sep="\t"), check.names=FALSE)

anova_df4 = anova_df3[anova_df3$BH < 0.05 & abs(anova_df3$anova_1/anova_df3$anova_2) >= 20,]
mapping2 = mapping
row.names(mapping2) = mapping2$`#SampleID`;  mapping2$`#SampleID` = NULL;
## exploratory hm first
df_KO3 = df_KO[df_KO$KO %in% unique(anova_df4$KO),]
row.names(df_KO3) = paste0(df_KO3$KO, " - ", df_KO3$Type); df_KO3$KO = NULL; df_KO3$Type = NULL;
pheatmap::pheatmap(
  df_KO3, 
  file=paste0(outdir, "figures/tax_KO_test_rock_heatmap.pdf"),
  fontsize_row=9, 
  fontsize_col=3, 
  cellwidth=5,  
  cellheight=10,
  annotation=mapping2, 
  color=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100),
  clustering_method="average"
)
selected_KOs = c("K09001", "K22020", "K09001", "K00363", "K18829")
# Now that we have the interesting KOs in hand, lets generate a heatmap
df_KO2 = df_KO
df_KO2$KO_Type = paste0(df_KO2$KO, "_", df_KO2$Type)
df_KO2 = df_KO2[df_KO2$KO %in% selected_KOs,]
## Already normalized

# melt + ggplot geom_tile
abundance_KO2 = melt(df_KO2)
abundance_KO2 = merge(abundance_KO2, mapping[c("#SampleID", "Nitrogen")], by.x="variable", by.y="#SampleID")
abundance_KO2$workflow = gsub("_all", " - arm", abundance_KO2$Type)
abundance_KO2$workflow = factor(abundance_KO2$workflow, levels=c(
  "0.1M", "0.1M - arm", "0.5M", "0.5M - arm", "1M", "1M - arm", "4M",  "4M - arm",
  "8M", "8M - arm", "12M", "12M - arm", "complete"
))
KO_order = unique(abundance_KO2$KO)
abundance_KO2$KO = factor(abundance_KO2$KO, levels=KO_order)
keggs = data.frame(fread("~/databases/kegg_pathways_modules.tsv", header=TRUE), check.names=FALSE)
keggs$plevel2 = ifelse(keggs$KO %in% selected_KOs, keggs$KO_desc, keggs$plevel2) 
abundance_KO2 = merge(abundance_KO2, keggs[,c("KO", "KO_desc", "plevel2", "ko_desc")], by.x="KO", by.y="KO", all.x=TRUE, all.y=FALSE)
abundance_KO2$KO_plevel2 = paste0(abundance_KO2$KO, " - ", abundance_KO2$plevel2)
unique(abundance_KO2$KO_plevel2[grepl(paste(KO_order, collapse="|"), unique(abundance_KO2$KO_plevel2))])
KO_plevel2_order = unique(unlist(lapply(KO_order, grep,abundance_KO2$KO_plevel2, value = TRUE)))
abundance_KO2$KO_plevel2 = factor(abundance_KO2$KO_plevel2, levels=KO_plevel2_order)
abundance_KO2$workflow2 = gsub("_all", " - arm", abundance_KO2$workflow)
abundance_KO2$workflow3 = ifelse(grepl("arm", abundance_KO2$workflow2), "arm", "standard")
abundance_KO2$number_of_clusters = gsub("^(.*M).*", "\\1", abundance_KO2$workflow)
abundance_KO2$number_of_clusters = factor(abundance_KO2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))

abundance_KO2 = abundance_KO2[abundance_KO2$KO %in% c("K00363", "K09001", "K18829"),]
p_KO_soil_3 = ggplot(abundance_KO2, aes(y=value, x=number_of_clusters, fill=workflow3)) +
  geom_boxplot(aes(fill=workflow3), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow3),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  facet_grid(KO_plevel2 ~ Nitrogen, scales="free", space="free_x") +
  xlab("Number of clusters") + ylab("reads/recA reads") + 
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.5),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title.y=element_text(size=9),
    axis.title.x=element_text(size=9),
    axis.ticks.x = element_blank(),
    plot.title = element_text(lineheight=1.2, face="bold", size=11),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_blank(),
    legend.position="top",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
p_KO_soil_3

########################
# Mock communities  KO #
########################
mapping_file = "../mock/0.1M_clusters/mapping_file.tsv"
df_KO = get_KO_abundance_as_df(my_KO_files_mock, my_gene_abundance_files_mock, my_COG_files_mock, mapping_file)
#df_KO = data.frame(fread(paste0(outdir, "./df_KO_mock.tsv"), sep="\t", header=TRUE), check.names=TRUE); row.names(df_KO) = df_KO$V1; df_KO$V1=NULL;
head(df_KO)
mapping = data.frame(fread("../mock/0.1M_clusters/mapping_file.tsv"), check.names=FALSE)
samples = unique(mapping$`#SampleID`)
curr_mapping = mapping

my_workflows = unique(df_KO$Type)
df_KO = df_KO[,c(samples, "KO", "Type")]
ps_ko = list()
anova_df3 = NULL
# for simplicity sake, only show quartile 1 and 5 
for(i in 1:length(my_workflows)){
  print(paste0("i: ", i))
  curr_workflow = my_workflows[i]
  curr_KO = df_KO[df_KO$Type == curr_workflow,]
  KOs = unique(curr_KO$KO)
  KOs = KOs[KOs != "NULL"]
  
  curr_KO = melt(curr_KO)
  curr_KO = merge(curr_KO, curr_mapping, by.x="variable", by.y="#SampleID")
  
  anova_df = NULL
  anova_df2 = NULL
  for(j in 1:length(KOs)){
    print(paste0("   j: ", j))
    curr_KO2 = curr_KO[curr_KO$KO == KOs[j],]
    if(nrow(curr_KO2) == 0){next;}
    res = aov(value ~ Treatment, data=curr_KO2)
    pvalue = summary(res)[[1]][["Pr(>F)"]][1]
    
    anova_1 = mean(curr_KO2[curr_KO2$Treatment == "Even",]$value)
    anova_2 = mean(curr_KO2[curr_KO2$Treatment == "Staggered",]$value)
    
    if(j == 1){
      anova_df = data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2)
    }else{
      anova_df = rbind(anova_df, data.frame(KO=KOs[j], pvalue=pvalue, anova_1=anova_1, anova_2=anova_2))
    }
    
  }
  anova_df$BH = p.adjust(anova_df$pvalue, method = "bonferroni")
  anova_df2 = anova_df[anova_df$pvalue < 0.05 & abs(anova_df$anova_1/anova_df$anova_2) >= 0,] # use lenient value, but use stringent later (say fc 11)
  anova_df2 = na.omit(anova_df2)
  anova_df2$Type = curr_workflow
  
  if(i == 1){
    anova_df3 = anova_df2
  }else{
    anova_df3 = rbind(anova_df3, anova_df2)
  }
  
}

# Then filter results from the raw table (also write table to avoid regenerating anovas)
#write.table(anova_df3, paste0(outdir, "/anova_df3_mock.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
anova_df3 = data.frame(fread(paste0(outdir, "/anova_df3_mock.tsv"), header=TRUE, sep="\t"), check.names=FALSE)

anova_df4 = anova_df3[anova_df3$BH < 0.001 & abs(anova_df3$anova_1/anova_df3$anova_2) >= 200,]
mapping2 = mapping
row.names(mapping2) = mapping2$`#SampleID`;  mapping2$`#SampleID` = NULL;
## exploratory hm first
df_KO3 = df_KO[df_KO$KO %in% unique(anova_df4$KO),]
row.names(df_KO3) = paste0(df_KO3$KO, " - ", df_KO3$Type); df_KO3$KO = NULL; df_KO3$Type = NULL;
pheatmap::pheatmap(
  df_KO3, 
  file=paste0(outdir, "figures/tax_KO_test_mock_heatmap.pdf"),
  fontsize_row=9, 
  fontsize_col=3, 
  cellwidth=5,  
  cellheight=10,
  annotation=mapping2, 
  color=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100),
  clustering_method="average"
)
selected_KOs = c("K21572", "K11636", "K02118", "K18197", "K21298")
# Now that we have the interesting KOs in hand, lets generate a heatmap
df_KO2 = df_KO
df_KO2$KO_Type = paste0(df_KO2$KO, "_", df_KO2$Type)
df_KO2 = df_KO2[df_KO2$KO %in% selected_KOs,]
## Already normalized

# melt + ggplot geom_tile
abundance_KO2 = melt(df_KO2)
abundance_KO2 = merge(abundance_KO2, mapping[c("#SampleID", "Treatment")], by.x="variable", by.y="#SampleID")
abundance_KO2$workflow = gsub("_all", " - arm", abundance_KO2$Type)
abundance_KO2$workflow = factor(abundance_KO2$workflow, levels=c(
  "0.1M", "0.1M - arm", "0.5M", "0.5M - arm", "1M", "1M - arm", "2M",  "2M - arm",
  "3M", "3M - arm", "4M", "4M - arm", "complete"
))
KO_order = unique(abundance_KO2$KO)
abundance_KO2$KO = factor(abundance_KO2$KO, levels=KO_order)
keggs = data.frame(fread("~/databases/kegg_pathways_modules.tsv", header=TRUE), check.names=FALSE)
keggs$plevel2 = ifelse(keggs$KO %in% selected_KOs, keggs$KO_desc, keggs$plevel2) 
abundance_KO2 = merge(abundance_KO2, keggs[,c("KO", "KO_desc", "plevel2", "ko_desc")], by.x="KO", by.y="KO", all.x=TRUE, all.y=FALSE)
abundance_KO2$KO_plevel2 = paste0(abundance_KO2$KO, " - ", abundance_KO2$plevel2)
unique(abundance_KO2$KO_plevel2[grepl(paste(KO_order, collapse="|"), unique(abundance_KO2$KO_plevel2))])
KO_plevel2_order = unique(unlist(lapply(KO_order, grep,abundance_KO2$KO_plevel2, value = TRUE)))
abundance_KO2$KO_plevel2 = factor(abundance_KO2$KO_plevel2, levels=KO_plevel2_order)
abundance_KO2$workflow2 = gsub("_all", " - arm", abundance_KO2$workflow)
abundance_KO2$workflow3 = ifelse(grepl("arm", abundance_KO2$workflow2), "arm", "standard")
abundance_KO2$number_of_clusters = gsub("^(.*M).*", "\\1", abundance_KO2$workflow)
abundance_KO2$number_of_clusters = factor(abundance_KO2$number_of_clusters, levels=c("0.1M", "0.5M", "1M", "4M", "8M", "12M", "complete"))
#abundance_KO2$Location  = gsub("\\.", " ", abundance_KO2$Location)
abundance_KO2$KO_plevel2 = gsub("K02118 - ATPVB, ntpB, atpB; V/A-type H\\+/Na\\+-transporting ATPase subunit B",
                                "K02118 - ATPVB ntpB atpB; V/A-type \n H Na -transporting ATPase subunit B", abundance_KO2$KO_plevel2)
abundance_KO2$KO_plevel2 = gsub("K11636 - yxdM; putative ABC transport system permease protein",
                                "K11636 - yxdM; putative ABC transport\n system permease protein", abundance_KO2$KO_plevel2)
abundance_KO2$KO_plevel2 = gsub("K21572 - susD; starch-binding outer membrane protein, SusD/RagB family",
                                "K21572 - susD; starch-binding outer membrane protein,\n SusD/RagB family", abundance_KO2$KO_plevel2)

abundance_KO2 = abundance_KO2[abundance_KO2$KO %in% c("K02118", "K11636", "K21572"),]
p_KO_mock_3 = ggplot(abundance_KO2, aes(y=value, x=number_of_clusters, fill=workflow3)) +
  geom_boxplot(aes(fill=workflow3), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow3),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  facet_grid(KO_plevel2 ~ Treatment, scales="free", space="free_x") +
  xlab("Number of clusters") + ylab("reads/recA reads") + 
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.5),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(90), hjust=1),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title.y=element_text(size=9),
    axis.title.x=element_text(size=9),
    axis.ticks.x = element_blank(),
    plot.title = element_text(lineheight=1.2, face="bold", size=11),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=9, face="plain"),
    legend.title = element_blank(),
    legend.position="top",
    strip.text.x = element_text(angle=0, vjust=0,hjust=0.5, size=9, face="plain"),
    strip.text.y = element_text(angle=0, hjust=0, size=9, face="plain"),
    strip.background =  element_blank(),
    panel.spacing = unit(0.3, "lines")
  )
p_KO_mock_3

# Finally, generate figures.
fig_rock = ggarrange(
  p_KO_rock_3,
  nrow = 2, ncol= 1, heights=c(5, 1)
)


fig_soil = ggarrange(
  p_KO_soil_3 + rremove("legend"), 
  nrow = 2, ncol= 1, heights=c(1.5,1)
)

fig_mock = ggarrange(
  p_KO_mock_3,
  nrow = 2, ncol= 1, heights=c(1.85, 1)
)

f10a = annotate_figure(
  ggarrange(
    p_KO_china_3 + rremove("legend"), 
    fig_rock,# + rremove("legend"),
    fig_soil, 
    fig_mock,
    labels=c("A", "B", "C", "D"),
    nrow = 2, ncol= 2
  ),
  top = text_grob("Figure 5\n", color = "black", face = "plain", size = 11),
)

pdf( file=paste0(outdir, "./figures/F_FIGURE_5.pdf"), height=10, width=14)
f10a
dev.off()

##############################################
# Reproducibility investigations             #
##############################################
my_0.1M_tax = c(
  "1" =   "../china_microbiome/reproducibility/0.1M_clusters/1/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2" =   "../china_microbiome/reproducibility/0.1M_clusters/2/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3" =   "../china_microbiome/reproducibility/0.1M_clusters/3/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)
my_0.5M_tax = c(
  "1" =   "../china_microbiome/reproducibility/0.5M_clusters/1/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2" =   "../china_microbiome/reproducibility/0.5M_clusters/2/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3" =   "../china_microbiome/reproducibility/0.5M_clusters/3/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)
my_1M_tax = c(
  "1" =   "../china_microbiome/reproducibility/1M_clusters/1/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2" =   "../china_microbiome/reproducibility/1M_clusters/2/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3" =   "../china_microbiome/reproducibility/1M_clusters/3/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)

my_0.1M_alpha = c(
  "1" =   "../china_microbiome/reproducibility/0.1M_clusters/1/alpha_div/gene_abundance/alpha_richness.tsv",
  "2" =   "../china_microbiome/reproducibility/0.1M_clusters/2/alpha_div/gene_abundance/alpha_richness.tsv",
  "3" =   "../china_microbiome/reproducibility/0.1M_clusters/3/alpha_div/gene_abundance/alpha_richness.tsv"
)
my_0.5M_alpha = c(
  "1" =   "../china_microbiome/reproducibility/0.5M_clusters/1/alpha_div/gene_abundance/alpha_richness.tsv",
  "2" =   "../china_microbiome/reproducibility/0.5M_clusters/2/alpha_div/gene_abundance/alpha_richness.tsv",
  "3" =   "../china_microbiome/reproducibility/0.5M_clusters/3/alpha_div/gene_abundance/alpha_richness.tsv"
)
my_1M_alpha = c(
  "1" =   "../china_microbiome/reproducibility/1M_clusters/1/alpha_div/gene_abundance/alpha_richness.tsv",
  "2" =   "../china_microbiome/reproducibility/1M_clusters/2/alpha_div/gene_abundance/alpha_richness.tsv",
  "3" =   "../china_microbiome/reproducibility/1M_clusters/3/alpha_div/gene_abundance/alpha_richness.tsv"
)

mapping = data.frame(fread("./mapping_file3.tsv", header=T), check.names=F)
row.names(mapping) = mapping$`#SampleID`

df_alpha_0.1M = get_alphadiv_as_df(my_0.1M_alpha, mapping)
df_alpha_0.1M$subsampling = "0.1M"
df_alpha_0.5M = get_alphadiv_as_df(my_0.5M_alpha, mapping)
df_alpha_0.5M$subsampling = "0.5M"
df_alpha_1M = get_alphadiv_as_df(my_1M_alpha, mapping)
df_alpha_1M$subsampling = "1M"

df_alpha = rbind(df_alpha_0.1M, df_alpha_0.5M, df_alpha_1M)
colnames(df_alpha)[5] = "number_of_clusters"
df_alpha$workflow2 = "standard"
df_alpha$Type_subsampling = paste0(df_alpha$Type, "_", df_alpha$subsampling)

tax_df_0.1M = get_tax_as_df(my_0.1M_tax)
tax_df_0.1M$subsampling = "0.1M"
tax_df_0.5M = get_tax_as_df(my_0.5M_tax)
tax_df_0.5M$subsampling = "0.5M"
tax_df_1M = get_tax_as_df(my_1M_tax)
tax_df_1M$subsampling = "1M"

tax_df = rbind(tax_df_0.1M, tax_df_0.5M, tax_df_1M)
colnames(tax_df)[5] = "number_of_clusters"
tax_df$workflow2 = "standard"
tax_df2 = tax_df[grep(".*Prevotella.*", tax_df$Taxon),]

tax_df2$Taxon = gsub("k__Bacteria;p__NULL;c__NULL;o__NULL", "k__Bacteria", tax_df2$Taxon)
tax_df2$Taxon = gsub("p__NULL;c__NULL;o__NULL", "", tax_df2$Taxon)
tax_df2$Taxon = gsub("c__NULL;o__NULL", "", tax_df2$Taxon)
tax_df2$Taxon = gsub(";o__NULL", "", tax_df2$Taxon)
tax_df2$Taxon = gsub(";c__NULL", ";c__undefined", tax_df2$Taxon)
tax_df2$Taxon = gsub("k__Bacteria;", "", tax_df2$Taxon)
tax_df2$Taxon = gsub(";p__NULL;", ";p__undefined", tax_df2$Taxon)
tax_df2$Taxon = gsub(";s__NULL", "", tax_df2$Taxon)
tax_df2$Taxon = gsub(";f__NULL;g__NULL", "", tax_df2$Taxon)
tax_df2$Taxon = gsub(";g__NULL", "", tax_df2$Taxon)
tax_df2$Taxon = gsub(";$", "", tax_df2$Taxon)  

tax_df3 = tax_df2[grep(".*Prevotella copri$", tax_df2$Taxon),]
tax_df4 = tax_df2[grep(".*g__Prevotella$", tax_df2$Taxon),]
tax_df5 = rbind(tax_df3, tax_df4)
tax_df5$Type_subsampling = paste0(tax_df5$Type, "_", tax_df5$number_of_clusters)
df_taxa = tax_df5

curr_taxa2 = NULL
my_workflows = unique(df_alpha$Type_subsampling)
for(i in 1:length(my_workflows)){
  curr_workflow = my_workflows[i]
  curr_ad = df_alpha[df_alpha$Type_subsampling == curr_workflow,]
  
  # Divide by quintiles:
  curr_ad$quantiles = NULL
  curr_ad$quantiles = cut(curr_ad$value, 
                          quantile(
                            curr_ad$value, 
                            prob = seq(0, 1, length = 6), 
                            type = 5, names=TRUE, 
                            labels=FALSE, 
                            ordered_result=TRUE
                          )
  )
  
  df_link = data.frame(quantiles=levels(curr_ad$quantiles), rank=seq(1,(length(unique(curr_ad$quantiles))-1), by=1))
  df_tmp = curr_ad[,c("variable", "quantiles")]
  df_tmp = merge(df_tmp, df_link, by="quantiles")
  
  # Taxa
  curr_taxa = df_taxa[df_taxa$Type_subsampling == curr_workflow,]
  curr_taxa = merge(curr_taxa, df_tmp[, c("variable", "rank")], by="variable")
 
  
  # populate filtered taxa df for further usage.
  if(i == 1){
    curr_taxa2 = curr_taxa
  }else{
    curr_taxa2 = rbind(curr_taxa2, curr_taxa)
  }
}
  
colnames(curr_taxa2)[8] = "Diversity_quartile"
curr_taxa2$Diversity_quartile_rep = paste0(curr_taxa2$Diversity_quartile, " - rep #", curr_taxa2$Type)
color_list = c("arm" = "#0000CD", "standard" = "#FF0000")
curr_taxa2$value = curr_taxa2$value * 100
curr_taxa2 = curr_taxa2[curr_taxa2$Diversity_quartile %in% c(1,5),]
curr_taxa2$number_of_clusters_rep = paste0(curr_taxa2$number_of_clusters, "_", curr_taxa2$Type)
curr_taxa2$rank = as.character(curr_taxa2$Diversity_quartile)
curr_taxa2$rank_number_of_clusters = paste0("Diversity quantile #", curr_taxa2$rank, " - ", curr_taxa2$number_of_clusters, " clusters")

p_selected_tax <- ggplot(data=curr_taxa2, aes(x=Type, y=value, fill=workflow2)) + 
  facet_grid(Taxon ~ rank_number_of_clusters, scales="free_y", space="free_x") + 
  geom_boxplot(aes(fill=workflow2), outlier.shape = NA, position=position_dodge()) + 
  geom_jitter(aes(color=workflow2),size=0.005, pch=19,position=position_jitterdodge(0.1), alpha=0.5) +
  xlab("In silico replicates") + 
  ylab("abundance (%)") +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  guides(fill=guide_legend(ncol=2)) + ggtitle("") +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=0.65),
    axis.text.x=element_text(size=8, colour="black", angle=as.numeric(0), hjust=0.5),
    axis.text.y=element_text(size=8, colour="black"), 
    axis.title=element_text(family="Helvetica", size=11),
    plot.title = element_text(face="plain", size=14, hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size=8.5, face="plain", margin = margin(r = 0.5, unit = 'cm')),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.07, 'cm'),
    legend.position="bottom",
    strip.text.x = element_text(angle=90, vjust=0.5,hjust=0, size=10, face="plain"),
    strip.background =  element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.text.y = element_text(angle=0, vjust=0.5, hjust=0, size=10, face="plain")
  ) 
p_selected_tax

pdf( file=paste0(outdir, "./figures/FIGURE_S7.pdf"), height=6, width=10)
p_selected_tax
dev.off()
