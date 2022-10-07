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

options(stringsAsFactors = FALSE)

# Define color palettes.
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
#######################################################
# Define end results files of the shotgunMG pipeline  #
# used to generate figures.                           #
#######################################################

## QC files.
my_qc_files_china = c(
  "0.1M" = "./export/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M" = "./export/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M" = "./export/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M" = "./export/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "8M" = "./export/8M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "12M" = "./export/12M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/contig_abundance/qc_mapping_stats.tsv"
)

my_qc_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/contig_abundance/qc_mapping_stats.tsv"
)

my_qc_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M" = "../mock/export/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "2M" = "../mock/export/2M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "3M" = "../mock/export/3M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M" = "../mock/export/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "complete" = "../mock/export/all_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/contig_abundance/qc_mapping_stats.tsv"
)

my_qc_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M" = "../mock_even/export/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "2M" = "../mock_even/export/2M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "3M" = "../mock_even/export/3M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M" = "../mock_even/export/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "complete" = "../mock_even/export/all_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/contig_abundance/qc_mapping_stats.tsv"
)

my_qc_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "complete" = "../mock_staggered/export/all_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/contig_abundance/qc_mapping_stats.tsv"
)

my_qc_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M" = "../soil/export/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M" = "../soil/export/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "8M" = "../soil/export/8M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "12M" = "../soil/export/12M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "complete" = "../soil/export/all_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/contig_abundance/qc_mapping_stats.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/contig_abundance/qc_mapping_stats.tsv"
)

## Coords files contigs
my_coords_files = c(
  "0.1M" = "./export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M" = "./export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M" = "./export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M" = "./export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "8M" = "./export/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "12M" = "./export/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv"
)

my_coords_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv"
)

my_coords_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M" = "../mock/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "2M" = "../mock/export/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "3M" = "../mock/export/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M" = "../mock/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "complete" = "../mock/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv"
)

my_coords_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M" = "../mock_even/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "2M" = "../mock_even/export/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "3M" = "../mock_even/export/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M" = "../mock_even/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "complete" = "../mock_even/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv"
)

my_coords_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "complete" = "../mock_staggered/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv"
)

my_coords_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M" = "../soil/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M" = "../soil/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "8M" = "../soil/export/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "12M" = "../soil/export/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "complete" = "../soil/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance_coords.tsv"
)

## Coords files genes
my_coords_files_genes = c(
  "0.1M" = "./export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M" = "./export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M" = "./export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M" = "./export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "8M" = "./export/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "12M" = "./export/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv"
)

my_coords_files_rock_genes = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv"
)

my_coords_files_mock_genes = c(
  "0.1M" = "../mock/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M" = "../mock/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "2M" = "../mock/export/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "3M" = "../mock/export/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M" = "../mock/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "complete" = "../mock/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv"
)

my_coords_files_mock_even_genes = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M" = "../mock_even/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "2M" = "../mock_even/export/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "3M" = "../mock_even/export/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M" = "../mock_even/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "complete" = "../mock_even/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv"
)

my_coords_files_mock_staggered_genes = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "complete" = "../mock_staggered/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv"
)

my_coords_files_soil_genes = c(
  "0.1M" = "../soil/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M" = "../soil/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M" = "../soil/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "8M" = "../soil/export/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "12M" = "../soil/export/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "complete" = "../soil/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance_coords.tsv"
)

# Dist files
my_dist_files = c(
  "0.1M" = "./export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M" = "./export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M" = "./export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M" = "./export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "8M" = "./export/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "12M" = "./export/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv"
)

my_dist_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv"
)

my_dist_files_mock = c(
  "0.1M" =     "../mock/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M" =     "../mock/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M" =       "../mock/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "2M" =       "../mock/export/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "3M" =       "../mock/export/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M" =       "../mock/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "complete" = "../mock/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M_all" =   "../mock/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "2M_all" =   "../mock/export/map_all_reads/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "3M_all" =   "../mock/export/map_all_reads/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M_all" =   "../mock/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv"
)

my_dist_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M" = "../mock_even/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "2M" = "../mock_even/export/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "3M" = "../mock_even/export/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M" = "../mock_even/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "complete" = "../mock_even/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv"
)

my_dist_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "complete" = "../mock_staggered/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv"
)

my_dist_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M" = "../soil/export/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M" = "../soil/export/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "8M" = "../soil/export/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "12M" = "../soil/export/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "complete" = "../soil/export/all_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/beta_div/bray_curtis_contig_abundance/bray_curtis_contig_abundance.tsv"
)

# dist files genes
my_dist_files_genes = c(
  "0.1M" = "./export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M" = "./export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M" = "./export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M" = "./export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "8M" = "./export/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "12M" = "./export/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv"
)

my_dist_files_rock_genes = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv"
)

my_dist_files_mock_genes = c(
  "0.1M" =     "../mock/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M" =     "../mock/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M" =       "../mock/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "2M" =       "../mock/export/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "3M" =       "../mock/export/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M" =       "../mock/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "complete" = "../mock/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M_all"   = "../mock/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "2M_all"   = "../mock/export/map_all_reads/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "3M_all"   = "../mock/export/map_all_reads/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M_all"   = "../mock/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv"
)

my_dist_files_mock_even_genes = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M" = "../mock_even/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "2M" = "../mock_even/export/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "3M" = "../mock_even/export/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M" = "../mock_even/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "complete" = "../mock_even/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv"
)

my_dist_files_mock_staggered_genes = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "complete" = "../mock_staggered/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv"
)

my_dist_files_soil_genes = c(
  "0.1M" = "../soil/export/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M" = "../soil/export/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M" = "../soil/export/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "8M" = "../soil/export/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "12M" = "../soil/export/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "complete" = "../soil/export/all_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/beta_div/bray_curtis_gene_abundance/bray_curtis_gene_abundance.tsv"
)

# Alpha richness contigs
my_alpha_div_contigs_files = c(
  "0.1M" = "./export/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M" = "./export/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M" = "./export/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M" = "./export/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "8M" = "./export/8M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "12M" = "./export/12M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/alpha_div/contig_abundance/alpha_richness.tsv"
)

my_alpha_div_contigs_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/alpha_div/contig_abundance/alpha_richness.tsv"
)

my_alpha_div_contigs_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M" = "../mock/export/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "2M" = "../mock/export/2M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "3M" = "../mock/export/3M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M" = "../mock/export/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "complete" = "../mock/export/all_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv"
)

my_alpha_div_contigs_files_mock_even = c(
  "0.1M" = "../moc_even/export/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M" = "../mock_even/export/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "2M" = "../mock_even/export/2M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "3M" = "../mock_even/export/3M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M" = "../mock_even/export/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "complete" = "../mock_even/export/all_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M_all" = "../mock/export_even/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "2M_all" = "../mock/export_even/map_all_reads/2M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "3M_all" = "../mock/export_even/map_all_reads/3M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M_all" = "../mock/export_even/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv"
)

my_alpha_div_contigs_files_mock_staggered = c(
  "0.1M" = "../moc_staggered/export/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "complete" = "../mock_staggered/export/all_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M_all" = "../mock/export_staggered/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "2M_all" = "../mock/export_staggered/map_all_reads/2M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "3M_all" = "../mock/export_staggered/map_all_reads/3M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M_all" = "../mock/export_staggered/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv"
)

my_alpha_div_contigs_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M" = "../soil/export/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M" = "../soil/export/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "8M" = "../soil/export/8M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "12M" = "../soil/export/12M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "complete" = "../soil/export/all_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/alpha_div/contig_abundance/alpha_richness.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/alpha_div/contig_abundance/alpha_richness.tsv"
)

# alpha richness genes
my_alpha_div_genes_files = c(
  "0.1M" = "./export/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M" = "./export/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M" = "./export/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M" = "./export/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "8M" = "./export/8M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "12M" = "./export/12M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/alpha_div/gene_abundance/alpha_richness.tsv"
)

my_alpha_div_genes_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/alpha_div/gene_abundance/alpha_richness.tsv"
)

my_alpha_div_genes_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M" = "../mock/export/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "2M" = "../mock/export/2M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "3M" = "../mock/export/3M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M" = "../mock/export/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "complete" = "../mock/export/all_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv"
)

my_alpha_div_genes_files_mock_even = c(
  "0.1M" = "../moc_even/export/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M" = "../mock_even/export/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "2M" = "../mock_even/export/2M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "3M" = "../mock_even/export/3M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M" = "../mock_even/export/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "complete" = "../mock_even/export/all_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M_all" = "../mock/export_even/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "2M_all" = "../mock/export_even/map_all_reads/2M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "3M_all" = "../mock/export_even/map_all_reads/3M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M_all" = "../mock/export_even/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv"
)

my_alpha_div_genes_files_mock_staggered = c(
  "0.1M" = "../moc_staggered/export/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "complete" = "../mock_staggered/export/all_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M_all" = "../mock/export_staggered/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "2M_all" = "../mock/export_staggered/map_all_reads/2M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "3M_all" = "../mock/export_staggered/map_all_reads/3M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M_all" = "../mock/export_staggered/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv"
)

my_alpha_div_genes_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M" = "../soil/export/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M" = "../soil/export/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "8M" = "../soil/export/8M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "12M" = "../soil/export/12M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "complete" = "../soil/export/all_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/alpha_div/gene_abundance/alpha_richness.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/alpha_div/gene_abundance/alpha_richness.tsv"
)

# # Alpha chao1 contigsd
# my_alpha_chao1_div_contigs_files = c(
#   "0.1M" = "./export/0.1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "0.5M" = "./export/0.5M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "1M" = "./export/1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "4M" = "./export/4M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "8M" = "./export/8M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "12M" = "./export/12M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "0.1M_all" = "./export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "0.5M_all" = "./export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "1M_all" = "./export/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "4M_all" = "./export/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "8M_all" = "./export/map_all_reads/8M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "12M_all" = "./export/map_all_reads/12M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv"
# )
# 
# my_alpha_chao1_div_genes_files = c(
#   "0.1M" = "./export/0.1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "0.5M" = "./export/0.5M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "1M" = "./export/1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "4M" = "./export/4M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "8M" = "./export/8M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "12M" = "./export/12M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "0.1M_all" = "./export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "0.5M_all" = "./export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "1M_all" = "./export/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "4M_all" = "./export/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "8M_all" = "./export/map_all_reads/8M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "12M_all" = "./export/map_all_reads/12M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv"
# )
# 
# my_alpha_chao1_div_contigs_files_rock = c(
#   "0.1M" = "../rock_microbiome/export/0.1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "0.5M" = "../rock_microbiome/export/0.5M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "1M" = "../rock_microbiome/export/1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "4M" = "../rock_microbiome/export/4M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "8M" = "../rock_microbiome/export/8M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "12M" = "../rock_microbiome/export/12M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "complete" = "../rock_microbiome/export/all_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv",
#   "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/alpha_div/contig_abundance/alpha_chao1.tsv"
# )
# 
# my_alpha_chao1_div_genes_files_rock = c(
#   "0.1M" = "../rock_microbiome/export/0.1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "0.5M" = "../rock_microbiome/export/0.5M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "1M" = "../rock_microbiome/export/1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "4M" = "../rock_microbiome/export/4M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "8M" = "../rock_microbiome/export/8M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "12M" = "../rock_microbiome/export/12M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "complete" = "../rock_microbiome/export/all_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv",
#   "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/alpha_div/gene_abundance/alpha_chao1.tsv"
# )

# tax files L6 relative
my_tax_files_china = c(
  "0.1M" = "./export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M" = "./export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M" = "./export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M" = "./export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "8M" = "./export/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "12M" = "./export/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv"
)

my_tax_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv"
)

my_tax_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M" = "../mock/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "2M" = "../mock/export/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "3M" = "../mock/export/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M" = "../mock/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "complete" = "../mock/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv"
)

my_tax_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M" = "../mock_even/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "2M" = "../mock_even/export/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "3M" = "../mock_even/export/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M" = "../mock_even/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "complete" = "../mock_even/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv"
)

my_tax_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "complete" = "../mock_staggered/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv"
)

my_tax_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M" = "../soil/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M" = "../soil/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "8M" = "../soil/export/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "12M" = "../soil/export/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "complete" = "../soil/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L6.tsv"
)

# Tax files relative L7
my_tax_files_china_species = c(
  "0.1M" = "./export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M" = "./export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M" = "./export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M" = "./export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "8M" = "./export/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "12M" = "./export/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)

my_tax_files_rock_species = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)

my_tax_files_mock_species = c(
  "0.1M" = "../mock/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M" = "../mock/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2M" = "../mock/export/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3M" = "../mock/export/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M" = "../mock/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "complete" = "../mock/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)

my_tax_files_mock_even_species = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M" = "../mock_even/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2M" = "../mock_even/export/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3M" = "../mock_even/export/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M" = "../mock_even/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "complete" = "../mock_even/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)

my_tax_files_mock_staggered_species = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "complete" = "../mock_staggered/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)

my_tax_files_soil_species = c(
  "0.1M" = "../soil/export/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M" = "../soil/export/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M" = "../soil/export/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "8M" = "../soil/export/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "12M" = "../soil/export/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "complete" = "../soil/export/all_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/relative/feature_table_L7.tsv"
)


## Tax files absolute L6
my_tax_files_china_counts = c(
  "0.1M" = "./export/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M" = "./export/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M" = "./export/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M" = "./export/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "8M" = "./export/8M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "12M" = "./export/12M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"
)

my_tax_files_rock_counts = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"
)

my_tax_files_mock_counts = c(
  "0.1M" = "../mock/export/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M" = "../mock/export/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "2M" = "../mock/export/2M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "3M" = "../mock/export/3M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M" = "../mock/export/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "complete" = "../mock/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"
)

my_tax_files_mock_even_counts = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M" = "../mock_even/export/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "2M" = "../mock_even/export/2M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "3M" = "../mock_even/export/3M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M" = "../mock_even/export/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "complete" = "../mock_even/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"
)

my_tax_files_mock_staggered_counts = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "complete" = "../mock_staggered/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"
)

my_tax_files_soil_counts = c(
  "0.1M" = "../soil/export/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M" = "../soil/export/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M" = "../soil/export/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "8M" = "../soil/export/8M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "12M" = "../soil/export/12M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "complete" = "../soil/export/all_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/all/absolute/feature_table_L6.tsv"
)

## MAGs files
my_mags_files_china = c(
  "0.1M" = "./export/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M" = "./export/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M" = "./export/1M_clusters/binning/metabat2/out_checkm.txt",
  "4M" = "./export/4M_clusters/binning/metabat2/out_checkm.txt",
  "8M" = "./export/8M_clusters/binning/metabat2/out_checkm.txt",
  "12M" = "./export/12M_clusters/binning/metabat2/out_checkm.txt",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M_all" = "./export/map_all_reads/1M_clusters/binning/metabat2/out_checkm.txt",
  "4M_all" = "./export/map_all_reads/4M_clusters/binning/metabat2/out_checkm.txt",
  "8M_all" = "./export/map_all_reads/8M_clusters/binning/metabat2/out_checkm.txt",
  "12M_all" = "./export/map_all_reads/12M_clusters/binning/metabat2/out_checkm.txt"
)

my_mags_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M" = "../rock_microbiome/export/1M_clusters/binning/metabat2/out_checkm.txt",
  "4M" = "../rock_microbiome/export/4M_clusters/binning/metabat2/out_checkm.txt",
  "8M" = "../rock_microbiome/export/8M_clusters/binning/metabat2/out_checkm.txt",
  "12M" = "../rock_microbiome/export/12M_clusters/binning/metabat2/out_checkm.txt",
  "complete" = "../rock_microbiome/export/all_clusters/binning/metabat2/out_checkm.txt",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/binning/metabat2/out_checkm.txt",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/binning/metabat2/out_checkm.txt",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/binning/metabat2/out_checkm.txt",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/binning/metabat2/out_checkm.txt"
)

my_mags_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M" = "../mock/export/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M" = "../mock/export/1M_clusters/binning/metabat2/out_checkm.txt",
  "2M" = "../mock/export/2M_clusters/binning/metabat2/out_checkm.txt",
  "3M" = "../mock/export/3M_clusters/binning/metabat2/out_checkm.txt",
  "4M" = "../mock/export/4M_clusters/binning/metabat2/out_checkm.txt",
  "complete" = "../mock/export/all_clusters/binning/metabat2/out_checkm.txt",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/binning/metabat2/out_checkm.txt",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/binning/metabat2/out_checkm.txt",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/binning/metabat2/out_checkm.txt",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/binning/metabat2/out_checkm.txt"
)

my_mags_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M" = "../mock_even/export/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M" = "../mock_even/export/1M_clusters/binning/metabat2/out_checkm.txt",
  "2M" = "../mock_even/export/2M_clusters/binning/metabat2/out_checkm.txt",
  "3M" = "../mock_even/export/3M_clusters/binning/metabat2/out_checkm.txt",
  "4M" = "../mock_even/export/4M_clusters/binning/metabat2/out_checkm.txt",
  "complete" = "../mock_even/export/all_clusters/binning/metabat2/out_checkm.txt",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/binning/metabat2/out_checkm.txt",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/binning/metabat2/out_checkm.txt",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/binning/metabat2/out_checkm.txt",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/binning/metabat2/out_checkm.txt"
)

my_mags_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M" = "../mock_staggered/export/1M_clusters/binning/metabat2/out_checkm.txt",
  "2M" = "../mock_staggered/export/2M_clusters/binning/metabat2/out_checkm.txt",
  "3M" = "../mock_staggered/export/3M_clusters/binning/metabat2/out_checkm.txt",
  "4M" = "../mock_staggered/export/4M_clusters/binning/metabat2/out_checkm.txt",
  "complete" = "../mock_staggered/export/all_clusters/binning/metabat2/out_checkm.txt",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/binning/metabat2/out_checkm.txt",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/binning/metabat2/out_checkm.txt",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/binning/metabat2/out_checkm.txt",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/binning/metabat2/out_checkm.txt"
)

my_mags_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M" = "../soil/export/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M" = "../soil/export/1M_clusters/binning/metabat2/out_checkm.txt",
  "4M" = "../soil/export/4M_clusters/binning/metabat2/out_checkm.txt",
  "8M" = "../soil/export/8M_clusters/binning/metabat2/out_checkm.txt",
  "12M" = "../soil/export/12M_clusters/binning/metabat2/out_checkm.txt",
  "complete" = "../soil/export/all_clusters/binning/metabat2/out_checkm.txt",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/binning/metabat2/out_checkm.txt",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/binning/metabat2/out_checkm.txt",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/binning/metabat2/out_checkm.txt",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/binning/metabat2/out_checkm.txt",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/binning/metabat2/out_checkm.txt",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/binning/metabat2/out_checkm.txt"
)

## KO files
my_KO_files_china = c(
  "0.1M" = "./export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M" = "./export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M" = "./export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M" = "./export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "8M" = "./export/8M_clusters/annotations/blastp_kegg_parsed.tsv",
  "12M" = "./export/12M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.1M_all" = "./export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M_all" = "./export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M_all" = "./export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M_all" = "./export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "8M_all" = "./export/8M_clusters/annotations/blastp_kegg_parsed.tsv",
  "12M_all" = "./export/12M_clusters/annotations/blastp_kegg_parsed.tsv"
)

my_KO_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/annotations/blastp_kegg_parsed.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/annotations/blastp_kegg_parsed.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.1M_all" = "../rock_microbiome/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M_all" = "../rock_microbiome/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M_all" = "../rock_microbiome/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M_all" = "../rock_microbiome/export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "8M_all" = "../rock_microbiome/export/8M_clusters/annotations/blastp_kegg_parsed.tsv",
  "12M_all" = "../rock_microbiome/export/12M_clusters/annotations/blastp_kegg_parsed.tsv"
)

my_KO_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M" = "../mock/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "2M" = "../mock/export/2M_clusters/annotations/blastp_kegg_parsed.tsv",
  "3M" = "../mock/export/3M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M" = "../mock/export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "complete" = "../mock/export/all_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.1M_all" = "../mock/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M_all" = "../mock/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M_all" = "../mock/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "2M_all" = "../mock/export/2M_clusters/annotations/blastp_kegg_parsed.tsv",
  "3M_all" = "../mock/export/3M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M_all" = "../mock/export/4M_clusters/annotations/blastp_kegg_parsed.tsv"
)

my_KO_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M" = "../mock_even/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "2M" = "../mock_even/export/2M_clusters/annotations/blastp_kegg_parsed.tsv",
  "3M" = "../mock_even/export/3M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M" = "../mock_even/export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "complete" = "../mock_even/export/all_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.1M_all" = "../mock_even/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M_all" = "../mock_even/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M_all" = "../mock_even/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "2M_all" = "../mock_even/export/2M_clusters/annotations/blastp_kegg_parsed.tsv",
  "3M_all" = "../mock_even/export/3M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M_all" = "../mock_even/export/4M_clusters/annotations/blastp_kegg_parsed.tsv"
)

my_KO_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/annotations/blastp_kegg_parsed.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "complete" = "../mock_staggered/export/all_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.1M_all" = "../mock_staggered/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M_all" = "../mock_staggered/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M_all" = "../mock_staggered/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "2M_all" = "../mock_staggered/export/2M_clusters/annotations/blastp_kegg_parsed.tsv",
  "3M_all" = "../mock_staggered/export/3M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M_all" = "../mock_staggered/export/4M_clusters/annotations/blastp_kegg_parsed.tsv"
)

my_KO_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M" = "../soil/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M" = "../soil/export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "8M" = "../soil/export/8M_clusters/annotations/blastp_kegg_parsed.tsv",
  "12M" = "../soil/export/12M_clusters/annotations/blastp_kegg_parsed.tsv",
  "complete" = "../soil/export/all_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.1M_all" = "../soil/export/0.1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "0.5M_all" = "../soil/export/0.5M_clusters/annotations/blastp_kegg_parsed.tsv",
  "1M_all" = "../soil/export/1M_clusters/annotations/blastp_kegg_parsed.tsv",
  "4M_all" = "../soil/export/4M_clusters/annotations/blastp_kegg_parsed.tsv",
  "8M_all" = "../soil/export/8M_clusters/annotations/blastp_kegg_parsed.tsv",
  "12M_all" = "../soil/export/12M_clusters/annotations/blastp_kegg_parsed.tsv"
)

# COG files
my_COG_files_china = c(
  "0.1M" = "./export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M" = "./export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M" = "./export/1M_clusters/annotations/rpsblast_cog.tsv",
  "4M" = "./export/4M_clusters/annotations/rpsblast_cog.tsv",
  "8M" = "./export/8M_clusters/annotations/rpsblast_cog.tsv",
  "12M" = "./export/12M_clusters/annotations/rpsblast_cog.tsv",
  "0.1M_all" = "./export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M_all" = "./export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M_all" = "./export/1M_clusters/annotations/rpsblast_cog.tsv",
  "4M_all" = "./export/4M_clusters/annotations/rpsblast_cog.tsv",
  "8M_all" = "./export/8M_clusters/annotations/rpsblast_cog.tsv",
  "12M_all" = "./export/12M_clusters/annotations/rpsblast_cog.tsv"
)

my_COG_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/annotations/rpsblast_cog.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/annotations/rpsblast_cog.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/annotations/rpsblast_cog.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/annotations/rpsblast_cog.tsv",
  "0.1M_all" = "../rock_microbiome/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M_all" = "../rock_microbiome/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M_all" = "../rock_microbiome/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "4M_all" = "../rock_microbiome/export/4M_clusters/annotations/rpsblast_cog.tsv",
  "8M_all" = "../rock_microbiome/export/8M_clusters/annotations/rpsblast_cog.tsv",
  "12M_all" = "../rock_microbiome/export/12M_clusters/annotations/rpsblast_cog.tsv"
)

my_COG_files_mock = c(
  "0.1M" =     "../mock/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M" =     "../mock/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M" =       "../mock/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "2M" =       "../mock/export/2M_clusters/annotations/rpsblast_cog.tsv",
  "3M" =       "../mock/export/3M_clusters/annotations/rpsblast_cog.tsv",
  "4M" =       "../mock/export/4M_clusters/annotations/rpsblast_cog.tsv",
  "complete" = "../mock/export/all_clusters/annotations/rpsblast_cog.tsv",
  "0.1M_all" = "../mock/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M_all" = "../mock/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M_all" =   "../mock/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "2M_all" =   "../mock/export/2M_clusters/annotations/rpsblast_cog.tsv",
  "3M_all" =   "../mock/export/3M_clusters/annotations/rpsblast_cog.tsv",
  "4M_all" =   "../mock/export/4M_clusters/annotations/rpsblast_cog.tsv"
)

my_COG_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M" = "../mock_even/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "2M" = "../mock_even/export/2M_clusters/annotations/rpsblast_cog.tsv",
  "3M" = "../mock_even/export/3M_clusters/annotations/rpsblast_cog.tsv",
  "4M" = "../mock_even/export/4M_clusters/annotations/rpsblast_cog.tsv",
  "complete" = "../mock_even/export/all_clusters/annotations/rpsblast_cog.tsv",
  "0.1M_all" = "../mock_even/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M_all" = "../mock_even/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M_all" = "../mock_even/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "2M_all" = "../mock_even/export/2M_clusters/annotations/rpsblast_cog.tsv",
  "3M_all" = "../mock_even/export/3M_clusters/annotations/rpsblast_cog.tsv",
  "4M_all" = "../mock_even/export/4M_clusters/annotations/rpsblast_cog.tsv"
)

my_COG_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/annotations/rpsblast_cog.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/annotations/rpsblast_cog.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/annotations/rpsblast_cog.tsv",
  "complete" = "../mock_staggered/export/all_clusters/annotations/rpsblast_cog.tsv",
  "0.1M_all" = "../mock_staggered/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M_all" = "../mock_staggered/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M_all" = "../mock_staggered/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "2M_all" = "../mock_staggered/export/2M_clusters/annotations/rpsblast_cog.tsv",
  "3M_all" = "../mock_staggered/export/3M_clusters/annotations/rpsblast_cog.tsv",
  "4M_all" = "../mock_staggered/export/4M_clusters/annotations/rpsblast_cog.tsv"
)

my_COG_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M" = "../soil/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "4M" = "../soil/export/4M_clusters/annotations/rpsblast_cog.tsv",
  "8M" = "../soil/export/8M_clusters/annotations/rpsblast_cog.tsv",
  "12M" = "../soil/export/12M_clusters/annotations/rpsblast_cog.tsv",
  "complete" = "../soil/export/all_clusters/annotations/rpsblast_cog.tsv",
  "0.1M_all" = "../soil/export/0.1M_clusters/annotations/rpsblast_cog.tsv",
  "0.5M_all" = "../soil/export/0.5M_clusters/annotations/rpsblast_cog.tsv",
  "1M_all" = "../soil/export/1M_clusters/annotations/rpsblast_cog.tsv",
  "4M_all" = "../soil/export/4M_clusters/annotations/rpsblast_cog.tsv",
  "8M_all" = "../soil/export/8M_clusters/annotations/rpsblast_cog.tsv",
  "12M_all" = "../soil/export/12M_clusters/annotations/rpsblast_cog.tsv"
)

# Gene abundance
my_gene_abundance_files_china = c(
  "0.1M" = "./export/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M" = "./export/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M" = "./export/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M" = "./export/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "8M" = "./export/8M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "12M" = "./export/12M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/gene_abundance/merged_gene_abundance.tsv"
)

my_gene_abundance_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/gene_abundance/merged_gene_abundance.tsv"
)

my_gene_abundance_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M" = "../mock/export/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "2M" = "../mock/export/2M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "3M" = "../mock/export/3M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M" = "../mock/export/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "complete" = "../mock/export/all_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/gene_abundance/merged_gene_abundance.tsv"
)

my_gene_abundance_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M" = "../mock_even/export/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "2M" = "../mock_even/export/2M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "3M" = "../mock_even/export/3M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M" = "../mock_even/export/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "complete" = "../mock_even/export/all_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/gene_abundance/merged_gene_abundance.tsv"
)

my_gene_abundance_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "complete" = "../mock_staggered/export/all_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/gene_abundance/merged_gene_abundance.tsv"
)

my_gene_abundance_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M" = "../soil/export/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M" = "../soil/export/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "8M" = "../soil/export/8M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "12M" = "../soil/export/12M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "complete" = "../soil/export/all_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/gene_abundance/merged_gene_abundance.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/gene_abundance/merged_gene_abundance.tsv"
)

# Feature tables
my_ft_files_china = c(
  "0.1M" = "./export/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M" = "./export/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M" = "./export/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M" = "./export/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "8M" = "./export/8M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "12M" = "./export/12M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/feature_table.tsv"
)

my_ft_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/feature_table.tsv"
)

my_ft_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M" = "../mock/export/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "2M" = "../mock/export/2M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "3M" = "../mock/export/3M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M" = "../mock/export/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "complete" = "../mock/export/all_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv"
)

my_ft_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M" = "../mock_even/export/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "2M" = "../mock_even/export/2M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "3M" = "../mock_even/export/3M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M" = "../mock_even/export/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "complete" = "../mock_even/export/all_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv"
)

my_ft_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "complete" = "../mock_staggered/export/all_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv"
)

my_ft_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M" = "../soil/export/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M" = "../soil/export/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "8M" = "../soil/export/8M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "12M" = "../soil/export/12M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "complete" = "../soil/export/all_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/annotations/taxonomy_consensus/feature_table.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/annotations/taxonomy_consensus/feature_table.tsv"
)

## Link bins
my_link_files_china = c(
  "0.1M" = "./export/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M" = "./export/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M" = "./export/1M_clusters/annotations/link_contig_to_gene.tsv",
  "4M" = "./export/4M_clusters/annotations/link_contig_to_gene.tsv",
  "8M" = "./export/8M_clusters/annotations/link_contig_to_gene.tsv",
  "12M" = "./export/12M_clusters/annotations/link_contig_to_gene.tsv",
  "0.1M_all" = "./export/map_all_reads/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M_all" = "./export/map_all_reads/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M_all" = "./export/map_all_reads/1M_clusters/annotations/link_contig_to_gene.tsv",
  "4M_all" = "./export/map_all_reads/4M_clusters/annotations/link_contig_to_gene.tsv",
  "8M_all" = "./export/map_all_reads/8M_clusters/annotations/link_contig_to_gene.tsv",
  "12M_all" = "./export/map_all_reads/12M_clusters/annotations/link_contig_to_gene.tsv"
)

my_link_files_rock = c(
  "0.1M" = "../rock_microbiome/export/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M" = "../rock_microbiome/export/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M" = "../rock_microbiome/export/1M_clusters/annotations/link_contig_to_gene.tsv",
  "4M" = "../rock_microbiome/export/4M_clusters/annotations/link_contig_to_gene.tsv",
  "8M" = "../rock_microbiome/export/8M_clusters/annotations/link_contig_to_gene.tsv",
  "12M" = "../rock_microbiome/export/12M_clusters/annotations/link_contig_to_gene.tsv",
  "complete" = "../rock_microbiome/export/all_clusters/annotations/link_contig_to_gene.tsv",
  "0.1M_all" = "../rock_microbiome/export/map_all_reads/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M_all" = "../rock_microbiome/export/map_all_reads/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M_all" = "../rock_microbiome/export/map_all_reads/1M_clusters/annotations/link_contig_to_gene.tsv",
  "4M_all" = "../rock_microbiome/export/map_all_reads/4M_clusters/annotations/link_contig_to_gene.tsv",
  "8M_all" = "../rock_microbiome/export/map_all_reads/8M_clusters/annotations/link_contig_to_gene.tsv",
  "12M_all" = "../rock_microbiome/export/map_all_reads/12M_clusters/annotations/link_contig_to_gene.tsv"
)

my_link_files_mock = c(
  "0.1M" = "../mock/export/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M" = "../mock/export/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M" = "../mock/export/1M_clusters/annotations/link_contig_to_gene.tsv",
  "2M" = "../mock/export/2M_clusters/annotations/link_contig_to_gene.tsv",
  "3M" = "../mock/export/3M_clusters/annotations/link_contig_to_gene.tsv",
  "4M" = "../mock/export/4M_clusters/annotations/link_contig_to_gene.tsv",
  "complete" = "../mock/export/all_clusters/annotations/link_contig_to_gene.tsv",
  "0.1M_all" = "../mock/export/map_all_reads/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M_all" = "../mock/export/map_all_reads/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M_all" = "../mock/export/map_all_reads/1M_clusters/annotations/link_contig_to_gene.tsv",
  "2M_all" = "../mock/export/map_all_reads/2M_clusters/annotations/link_contig_to_gene.tsv",
  "3M_all" = "../mock/export/map_all_reads/3M_clusters/annotations/link_contig_to_gene.tsv",
  "4M_all" = "../mock/export/map_all_reads/4M_clusters/annotations/link_contig_to_gene.tsv"
)

my_link_files_mock_even = c(
  "0.1M" = "../mock_even/export/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M" = "../mock_even/export/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M" = "../mock_even/export/1M_clusters/annotations/link_contig_to_gene.tsv",
  "2M" = "../mock_even/export/2M_clusters/annotations/link_contig_to_gene.tsv",
  "3M" = "../mock_even/export/3M_clusters/annotations/link_contig_to_gene.tsv",
  "4M" = "../mock_even/export/4M_clusters/annotations/link_contig_to_gene.tsv",
  "complete" = "../mock_even/export/all_clusters/annotations/link_contig_to_gene.tsv",
  "0.1M_all" = "../mock_even/export/map_all_reads/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M_all" = "../mock_even/export/map_all_reads/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M_all" = "../mock_even/export/map_all_reads/1M_clusters/annotations/link_contig_to_gene.tsv",
  "2M_all" = "../mock_even/export/map_all_reads/2M_clusters/annotations/link_contig_to_gene.tsv",
  "3M_all" = "../mock_even/export/map_all_reads/3M_clusters/annotations/link_contig_to_gene.tsv",
  "4M_all" = "../mock_even/export/map_all_reads/4M_clusters/annotations/link_contig_to_gene.tsv"
)

my_link_files_mock_staggered = c(
  "0.1M" = "../mock_staggered/export/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M" = "../mock_staggered/export/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M" = "../mock_staggered/export/1M_clusters/annotations/link_contig_to_gene.tsv",
  "2M" = "../mock_staggered/export/2M_clusters/annotations/link_contig_to_gene.tsv",
  "3M" = "../mock_staggered/export/3M_clusters/annotations/link_contig_to_gene.tsv",
  "4M" = "../mock_staggered/export/4M_clusters/annotations/link_contig_to_gene.tsv",
  "complete" = "../mock_staggered/export/all_clusters/annotations/link_contig_to_gene.tsv",
  "0.1M_all" = "../mock_staggered/export/map_all_reads/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M_all" = "../mock_staggered/export/map_all_reads/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M_all" = "../mock_staggered/export/map_all_reads/1M_clusters/annotations/link_contig_to_gene.tsv",
  "2M_all" = "../mock_staggered/export/map_all_reads/2M_clusters/annotations/link_contig_to_gene.tsv",
  "3M_all" = "../mock_staggered/export/map_all_reads/3M_clusters/annotations/link_contig_to_gene.tsv",
  "4M_all" = "../mock_staggered/export/map_all_reads/4M_clusters/annotations/link_contig_to_gene.tsv"
)


my_link_files_soil = c(
  "0.1M" = "../soil/export/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M" = "../soil/export/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M" = "../soil/export/1M_clusters/annotations/link_contig_to_gene.tsv",
  "4M" = "../soil/export/4M_clusters/annotations/link_contig_to_gene.tsv",
  "8M" = "../soil/export/8M_clusters/annotations/link_contig_to_gene.tsv",
  "12M" = "../soil/export/12M_clusters/annotations/link_contig_to_gene.tsv",
  "complete" = "../soil/export/all_clusters/annotations/link_contig_to_gene.tsv",
  "0.1M_all" = "../soil/export/map_all_reads/0.1M_clusters/annotations/link_contig_to_gene.tsv",
  "0.5M_all" = "../soil/export/map_all_reads/0.5M_clusters/annotations/link_contig_to_gene.tsv",
  "1M_all" = "../soil/export/map_all_reads/1M_clusters/annotations/link_contig_to_gene.tsv",
  "4M_all" = "../soil/export/map_all_reads/4M_clusters/annotations/link_contig_to_gene.tsv",
  "8M_all" = "../soil/export/map_all_reads/8M_clusters/annotations/link_contig_to_gene.tsv",
  "12M_all" = "../soil/export/map_all_reads/12M_clusters/annotations/link_contig_to_gene.tsv"
)

## rpoB richness files.
my_rpob_richness_files_china = c(
  "0.1M"  = "./export/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M"  = "./export/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M"  = "./export/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M"  = "./export/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "8M"  = "./export/8M_cluster/alpha_div/rpob_abundance/alpha_richness.tsv",
  "12M"  = "./export/12M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.1M_all"  = "./export/map_all_reads/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M_all"  = "./export/map_all_reads/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M_all"  = "./export/map_all_reads/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M_all"  = "./export/map_all_reads/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "8M_all"  = "./export/map_all_reads/8M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "12M_all"  = "./export/map_all_reads/12M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv"
)

my_rpob_richness_files_rock_microbiome = c(
  "0.1M"  = "../rock_microbiome/export/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M"  = "../rock_microbiome/export/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M"  = "../rock_microbiome/export/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M"  = "../rock_microbiome/export/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "8M"  = "../rock_microbiome/export/8M_cluster/alpha_div/rpob_abundance/alpha_richness.tsv",
  "12M"  = "../rock_microbiome/export/12M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "complete"  = "../rock_microbiome/export/all_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.1M_all"  = "../rock_microbiome/export/map_all_reads/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M_all"  = "../rock_microbiome/export/map_all_reads/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M_all"  = "../rock_microbiome/export/map_all_reads/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M_all"  = "../rock_microbiome/export/map_all_reads/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "8M_all"  = "../rock_microbiome/export/map_all_reads/8M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "12M_all"  = "../rock_microbiome/export/map_all_reads/12M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv"
)

my_rpob_richness_files_mock = c(
  "0.1M"  = "../mock/export/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M"  = "../mock/export/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M"  = "../mock/export/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "2M"  = "../mock/export/2M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "3M"  = "../mock/export/3M_cluster/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M"  = "../mock/export/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "complete" = "../mock/export/all_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.1M_all"  = "../mock/export/map_all_reads/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M_all"  = "../mock/export/map_all_reads/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M_all"  = "../mock/export/map_all_reads/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "2M_all"  = "../mock/export/map_all_reads/2M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "3M_all"  = "../mock/export/map_all_reads/3M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M_all"  = "../mock/export/map_all_reads/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv"
)

my_rpob_richness_files_mock_even = c(
  "0.1M"  = "../mock_even/export/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M"  = "../mock_even/export/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M"  = "../mock_even/export/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "2M"  = "../mock_even/export/2M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "3M"  = "../mock_even/export/3M_cluster/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M"  = "../mock_even/export/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "complete" = "../mock_even/export/all_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.1M_all"  = "../mock_even/export/map_all_reads/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M_all"  = "../mock_even/export/map_all_reads/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M_all"  = "../mock_even/export/map_all_reads/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "2M_all"  = "../mock_even/export/map_all_reads/2M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "3M_all"  = "../mock_even/export/map_all_reads/3M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M_all"  = "../mock_even/export/map_all_reads/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv"
)

my_rpob_richness_files_mock_staggered = c(
  "0.1M"  = "../mock_staggered/export/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M"  = "../mock_staggered/export/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M"  = "../mock_staggered/export/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "2M"  = "../mock_staggered/export/2M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "3M"  = "../mock_staggered/export/3M_cluster/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M"  = "../mock_staggered/export/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "complete" = "../mock_staggered/export/all_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.1M_all"  = "../mock_staggered/export/map_all_reads/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M_all"  = "../mock_staggered/export/map_all_reads/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M_all"  = "../mock_staggered/export/map_all_reads/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "2M_all"  = "../mock_staggered/export/map_all_reads/2M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "3M_all"  = "../mock_staggered/export/map_all_reads/3M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M_all"  = "../mock_staggered/export/map_all_reads/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv"
)

my_rpob_richness_files_soil = c(
  "0.1M"  = "../soil/export/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M"  = "../soil/export/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M"  = "../soil/export/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M"  = "../soil/export/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "8M"  = "../soil/export/8M_cluster/alpha_div/rpob_abundance/alpha_richness.tsv",
  "12M"  = "../soil/export/12M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "complete" = "../soil/export/all_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.1M_all"  = "../soil/export/map_all_reads/0.1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "0.5M_all"  = "../soil/export/map_all_reads/0.5M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "1M_all"  = "../soil/export/map_all_reads/1M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "4M_all"  = "../soil/export/map_all_reads/4M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "8M_all"  = "../soil/export/map_all_reads/8M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv",
  "12M_all"  = "../soil/export/map_all_reads/12M_clusters/alpha_div/rpob_abundance/alpha_richness.tsv"
)

get_rpobtax_abundance_as_df = function(gene_abundance_files, cog_files, ft_table_files, link_files, mapping_file=NULL){
  final_df = NULL
  mapping = NULL
  mapping = data.frame(fread(mapping_file, sep="\t"), check.names=FALSE)
  print(head(mapping))
  col_order = unique(mapping$`#SampleID`)
  
  for(y in 1:length(cog_files)){
    gene_abundance_file = gene_abundance_files[y]
    cog_file = cog_files[y]
    link_file = link_files[y]
    ft_file = ft_table_files[y]
    prefix = names(cog_file)
    print(cog_file)
    tCOG = data.frame(fread(cog_file, sep="\t", header=F, fill=TRUE), check.names=FALSE)
   
    tAbundance = data.frame(fread(gene_abundance_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    colnames(tAbundance)[1] = "gene_id"
    rpob = tCOG[grepl(".*COG0085.*", tCOG$V13),]$V1
    
    tLink = data.frame(fread(link_file, sep="\t", header=F, fill=FALSE), check.names=FALSE)
    tLink = tLink[tLink$V1 %in% rpob,]
    contigs = tLink$V2
    
    ft = data.frame(fread(ft_file, header=TRUE, sep="\t"), check.names=FALSE) 
    ft = ft[ft$`#FEATURE_ID` %in% contigs,]
    ft$taxonomy = gsub(";s__.*", "", ft$taxonomy)
    
    tTax = ft %>% 
      dplyr::group_by(taxonomy) %>% # means that the following summarize function will be done by the kegg_entry variable.
      dplyr::summarise_at(.vars = as.character(unique(mapping$`#SampleID`)),
                          .funs = c(sum="sum")) %>% 
      as.data.frame()
    colnames(tTax) = gsub("_sum", "", colnames(tTax))
    
    tTax$Type = prefix
    
    if(y == 1){
      final_df = tTax
    }else{
      final_df = rbind(final_df, tTax)
    }
  }
 
  return(final_df)
}

get_coords_as_list = function(coords_files){
  
  df_list = list()
  for(i in 1:length(coords_files)){
    coords_file = coords_files[i][[1]]
    curr_name = names(coords_files[i])
    
    curr_file = file(coords_file)
    open(curr_file)
    percentVar = read.table(coords_file, skip=4, nrow=1)
    percent1 = percentVar[[1]] * 100 
    percent2 = percentVar[[2]] * 100 
    percent1 = formatC(round(percent1, 2), big.mark=",",drop0trailing=TRUE, format="f")
    percent2 = formatC(round(percent2, 2), big.mark=",",drop0trailing=TRUE, format="f")
    close(curr_file)
    
    # Then :
    curr_file = file(coords_file)
    open(curr_file)
    # Read starting from beginning of coordinates, fill last rows with NAs
    tData = read.table(curr_file, skip=9, fill=TRUE)
    # Then remove NAs
    tData = na.omit(tData)
    # Keep 3 first rows only.
    tData2 = tData[,1:4]
    for(j in 2:ncol(tData2)){
      tData2[,j] = as.numeric(tData2[,j]) 
    }   
    tData3 = tData2
    
    colnames(tData3) = c("variable", "D1", "D2", "D3")
    close(curr_file)
    tData3$workflow = curr_name
    tData3$PCo1 = percent1
    tData3$PCo2 = percent2
   
    df_list[[i]] = tData3
  }
  
  return(df_list)
}

get_coords_as_df = function(coords_files){
  df = NULL

  for(i in 1:length(coords_files)){
    coords_file = coords_files[i][[1]]
    curr_name = names(coords_files[i])
    print(curr_name)
    print(coords_file)
    
    curr_file = file(coords_file)
    open(curr_file)
    percentVar = read.table(coords_file, skip=4, nrow=1)
    percent1 = percentVar[[1]] * 100 
    percent2 = percentVar[[2]] * 100 
    percent1 = formatC(round(percent1, 2), big.mark=",",drop0trailing=TRUE, format="f")
    percent2 = formatC(round(percent2, 2), big.mark=",",drop0trailing=TRUE, format="f")
    close(curr_file)
    
    # Then :
    curr_file = file(coords_file)
    open(curr_file)
    # Read starting from beginning of coordinates, fill last rows with NAs
    tData = read.table(curr_file, skip=9, fill=TRUE)
    # Then remove NAs
    tData = na.omit(tData)
    # Keep 3 first rows only.
    tData2 = tData[,1:4]
    for(j in 2:ncol(tData2)){
      tData2[,j] = as.numeric(tData2[,j]) 
    }   
    tData3 = tData2
   
    colnames(tData3) = c("variable", "D1", "D2", "D3")
    close(curr_file)
    tData3$workflow = curr_name
    tData3$PCo1 = percent1
    tData3$PCo2 = percent2
    if(i == 1){
      df = tData3
    }else{
      df = rbind(df, tData3)
    }
  }
  
  return(df)
}

get_dists_as_list = function(dists_files, debug=FALSE){
  
  df_list = list()
  for(i in 1:length(dists_files)){
    dist_file = dists_files[i][[1]]
    curr_name = names(dists_files[i])
    print(paste0("processing ", curr_name))
    
    tData = data.frame(fread(dist_file, header=T), check.names=F)
    row.names(tData) = tData$V1
    tData$V1 = NULL
    # make sure matrice is in the same order and contains equals amount of row and cols.
    if(isFALSE(debug)){
      samples_col = colnames(tData)
      samples_row = row.names(tData)
      if(all(samples_col != samples_row)){
        stop("row names is not identical to col names...")
      }
      tData = tData[,samples_col]
      tData = tData[samples_col,]
    }else{
      samples_col = colnames(tData)[1:20]
      samples_row = row.names(tData)[1:20]
      if(all(samples_col != samples_row)){
        stop("row names is not identical to col names...")
      }
      tData = tData[,colnames(tData) %in% samples_col]
      tData = tData[row.names(tData) %in% samples_col,]
      tData = tData[,samples_col]
      tData = tData[samples_col,]
    }
    
    for(j in 1:ncol(tData)){
      tData[,j] = as.numeric(tData[,j]) 
    }   
    
    df_list[[curr_name]] = tData
  }
  
  return(df_list)
}

get_alphadiv_as_list = function(alpha_div_files){
  depth = NULL
  df_list = list()
  
  for(y in 1:length(alpha_div_files)){
    alpha_div_file = alpha_div_files[y]
    prefix = names(alpha_div_file)
    tAlpha = data.frame(fread(alpha_div_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    tAlpha = tAlpha[!tAlpha$depth %in% c("taxonomy"),]
    
    j = 2
    for(i in 2:ncol(tAlpha)){
      if(anyNA(tAlpha[,i])){
        print(colnames(tAlpha)[i])
        print("broke at i:");print(i);
        break
      }
      j = j + 1
    }
    j = j - 1
    
    read_depth = as.character(colnames(tAlpha)[j])
    if(!is.null(depth)){
      read_depth = depth
    }
    
    tAlpha = tAlpha[,colnames(tAlpha) %in% c("depth", read_depth)]
    tAlpha = data.frame(t(tAlpha))
    colnames(tAlpha) = tAlpha[1,]
    tAlpha = tAlpha[-1,]
    tAlpha[,1] = NULL
    
    tAlpha = tAlpha[,colnames(tAlpha) %in% row.names(mapping)]
    
    tAlpha = tAlpha[(nrow(tAlpha)-9):nrow(tAlpha),]
    for(k in 1:ncol(tAlpha)){
      tAlpha[,k] = as.numeric(tAlpha[,k])
      
    }
    
    tAlpha = Filter(function(x)!all(is.na(x)), tAlpha)
    
    
    tAlpha2 = data.frame(colMeans(tAlpha), check.names=FALSE)
    
    index_type = basename(file_path_sans_ext(as.vector(alpha_div_file)))
    
    df = data.frame(tAlpha2, check.names=FALSE)
    
    colnames(df)[1] = "variable"
    colnames(df)[2] = "value"
    df$workflow = prefix
    df$read_depth = read_depth
    
    df_list[[y]] = df
  }
  
  return(df_list)
}

get_alphadiv_as_df = function(alpha_div_files, mapping=mapping){
  depth = NULL 
  final_df = NULL
  
  for(y in 1:length(alpha_div_files)){
    alpha_div_file = alpha_div_files[y]
    prefix = names(alpha_div_file)
    tAlpha = data.frame(fread(alpha_div_file, fill=TRUE, sep="\t", header=T), check.names=FALSE)
    tAlpha = tAlpha[!tAlpha$depth %in% c("taxonomy"),]
    
    j = 2
    for(i in 2:ncol(tAlpha)){
      if(anyNA(tAlpha[,i])){
        print("max depth:"); print(colnames(tAlpha)[i])
        print("broke at i:"); print(i);
        break
      }
      j = j + 1
    }
    j = j - 1 # -1 because the row before did not have NAs.
    
    read_depth = as.character(colnames(tAlpha)[j])
    if(!is.null(depth)){
      read_depth = depth
    }
    
    tAlpha = tAlpha[,colnames(tAlpha) %in% c("depth", read_depth)]
    tAlpha = data.frame(t(tAlpha), check.names=FALSE)
    colnames(tAlpha) = tAlpha[1,]
    tAlpha = tAlpha[-1,]
    #tAlpha[,1] = NULL
    
    tAlpha = tAlpha[,colnames(tAlpha) %in% mapping[,1]]
    
    tAlpha = tAlpha[(nrow(tAlpha)-9):nrow(tAlpha),]
    for(k in 1:ncol(tAlpha)){
      tAlpha[,k] = as.numeric(tAlpha[,k])
    }
    
    tAlpha = Filter(function(x)!all(is.na(x)), tAlpha)
    
    tAlpha2 = data.frame(colMeans(tAlpha), check.names=FALSE)
    
    index_type = basename(file_path_sans_ext(as.vector(alpha_div_file)))
    
    df = data.frame(tAlpha2, check.names=FALSE)
    
    df = merge(df, mapping, by.x="row.names", by.y="#SampleID")
    colnames(df)[1] = "variable"
    colnames(df)[2] = "value"
    df$Type = prefix
    df$read_depth = read_depth
    
    if(y == 1){
      final_df = df
    }else{
      final_df = rbind(final_df, df)
    }
  }
  return(final_df)
}

get_qc_as_df = function(qc_files){
  final_df = NULL
  
  for(y in 1:length(qc_files)){
    qc_file = qc_files[y]
    prefix = names(qc_file)
    tQc = data.frame(fread(qc_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    tQc$`properlyPaired%` = gsub("%", "", tQc$`properlyPaired%`)
    tQc = tQc[,c("sampleName", "properlyPaired", "properlyPaired%")]
    tQc$Type = prefix
    
    if(y == 1){
      final_df = tQc
    }else{
      final_df = rbind(final_df, tQc)
    }
  }
  return(final_df)
}

get_mags_as_df = function(mags_files, contam=90, complete=85, print_high_qual=FALSE){
  final_df = NULL
  
  for(y in 1:length(mags_files)){
    mags_file = mags_files[y]
    prefix = names(mags_file)
    tMAGs = data.frame(fread(mags_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    if(nrow(tMAGs) == 0){ next;}
    tMAGs = tMAGs[tMAGs$Completeness >= complete & tMAGs$Contamination <= contam,,drop=FALSE]
    
    tMAGs$Type = prefix
    
    if(isTRUE(print_high_qual)){
      print(paste0(prefix, " - ", nrow(tMAGs)))
    }
    
    if(y == 1){
      final_df = tMAGs
    }else{
      final_df = rbind(final_df, tMAGs)
    }
  }
  return(final_df)
}

get_tax_as_df = function(tax_files, convert_to_rel=FALSE){
  final_df = NULL
  
  for(y in 1:length(tax_files)){
    tax_file = tax_files[y]
    #prefix = prefixes[y]
    prefix = names(tax_file)
    print(tax_file)
    tTax = data.frame(fread(tax_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    
    if(isTRUE(convert_to_rel)){
      tTaxPerc = prop.table(data.matrix(tTax[,2:ncol(tTax)]), margin=2)*100
      tTaxPerc = data.frame(Taxon=tTax$Taxon, tTaxPerc, check.names=FALSE)
      tTax = tTaxPerc
    }
    tTax = reshape2::melt(tTax)
    
    tTax$Type = prefix
    
    if(y == 1){
      final_df = tTax
    }else{
      final_df = rbind(final_df, tTax)
    }
  }
  return(final_df)
}

get_tax_as_df_unmelted = function(tax_files, convert_to_rel=FALSE, mapping_file=NULL){
  final_df = NULL
  mapping = NULL
  mapping = data.frame(fread(mapping_file, sep="\t"), check.names=FALSE)
  print(head(mapping))
  col_order = unique(mapping$`#SampleID`)
  
  for(y in 1:length(tax_files)){
    tax_file = tax_files[y]
    #prefix = prefixes[y]
    prefix = names(tax_file)
    print(tax_file)
    tTax = data.frame(fread(tax_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    
    if(isTRUE(convert_to_rel)){
      tTaxPerc = prop.table(data.matrix(tTax[,2:ncol(tTax)]), margin=2)*100
      tTaxPerc = data.frame(Taxon=tTax$Taxon, tTaxPerc, check.names=FALSE)
      tTax = tTaxPerc
    }
    tTax = tTax[,c(col_order, "Taxon")]
    tTax$Type = prefix
    
    if(y == 1){
      final_df = tTax
    }else{
      final_df = rbind(final_df, tTax)
    }
  }
  return(final_df)
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

get_KO_abundance_as_df = function(KO_files, gene_abundance_files, cog_files, mapping_file=NULL){
  final_df = NULL
  mapping = NULL
  mapping = data.frame(fread(mapping_file, sep="\t"), check.names=FALSE)
  print(head(mapping))
  col_order = unique(mapping$`#SampleID`)
  
  for(y in 1:length(KO_files)){
    KO_file = KO_files[y]
    gene_abundance_file = gene_abundance_files[y]
    cog_file = cog_files[y]
    prefix = names(KO_file)
    print(KO_file)
    tKO = data.frame(fread(KO_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    colnames(tKO)[1] = "gene_id"
    colnames(tKO)[3] = "KO"
    print(cog_file)
    tCOG = data.frame(fread(cog_file, sep="\t", header=F, fill=TRUE), check.names=FALSE)
    head(tCOG)
    print(gene_abundance_file)
    tAbundance = data.frame(fread(gene_abundance_file, sep="\t", header=T, fill=TRUE), check.names=FALSE)
    colnames(tAbundance)[1] = "gene_id"
    reca = tCOG[grepl(".*COG0468.*", tCOG$V13),]$V1
    if(length(reca) == 0){
      print("No recA genes found...")
      reca = "gene_id_0" 
      # Here, artificially populate reca genes to at least 1 count.
      reca = rbind(tAbundance, c(reca, rep(1, length(col_order))))
    }else{
      reca = tAbundance[tAbundance$gene_id %in% reca,]
    }
    row.names(reca) = reca$gene_id
    reca$gene_id = NULL
    reca = reca[,col_order]
    for(w in 1:ncol(reca)){
      reca[,w] = as.numeric(reca[,w])
    }
    reca = rbind(reca, colSums(reca))
    
    tAbundance = merge(tAbundance, tKO[,c("gene_id", "KO")], by="gene_id")
    tAbundance[tAbundance$KO != "NULL", ]
    row.names(tAbundance) = tAbundance$gene_id
    tAbundance$gene_id = NULL
    # Then, the magic happens here with dplyr.
    tAbundance = tAbundance %>% 
      dplyr::group_by(KO) %>% # means that the following summarize function will be done by the kegg_entry variable.
      dplyr::summarise_at(.vars = as.character(unique(mapping$`#SampleID`)),
                          .funs = c(sum="sum")) %>% 
      as.data.frame()
    colnames(tAbundance) = gsub("_sum", "", colnames(tAbundance))
    tAbundance = tAbundance[!is.na(tAbundance$KO),]
    tAbundance = tAbundance[,c("KO", col_order)]
    row.names(tAbundance) = tAbundance$KO
    tAbundance$KO = NULL
    tAbundance = t(t(tAbundance) / as.numeric(reca[nrow(reca),]))
    tAbundance = data.frame(tAbundance)
    
    tAbundance$Type = prefix
    
    tAbundance$KO = row.names(tAbundance)
    
    if(y == 1){
      final_df = tAbundance
    }else{
      final_df = rbind(final_df, tAbundance)
    }
  }
  return(final_df)
}

##############################################
## END - Load functions and libraries ########
##############################################