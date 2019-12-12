# Load libraries and user functions ----------------------------------------------------------

# library(devtools)
# # withr::with_libpaths("/Users/heustonef/Desktop/10XGenomicsData/", code = devtools::install_version(package='Seurat', version = package_version('2.3.0')))
# devtools::install_version(package='Seurat', version = package_version('2.3.0'), )
# library(Seurat, lib.loc = "/Users/heustonef/Desktop/10XGenomicsData/")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("dplyr", "colorRamps", "monocle", "stringr", "plyr", "dplyr")
# Load libraries
ipak(packages)
library(Seurat)
sessionInfo()
# Create basic colour palette

basic_color_palette<-c("#cb4bbe",
                       "lightskyblue",
                       "grey37",
                       "#53ba48",
                       "moccasin",
                       "#dcca44",
                       "#502d70",
                       "#afd671",
                       "#cb4471",
                       "#69da99",
                       "#d44732",
                       "#6fd3cf",
                       "#5d2539",
                       "#cfd0a0",
                       "blue",
                       "#d2883b",
                       "#6793c0",
                       "#898938",
                       "#c98ac2",
                       "yellow",
                       "#c4c0cc",
                       "#7d3d23",
                       "#00a5ff",
                       "#d68d7a",
                       "#a2c1a3")

# Adjust colour palette
adjust_palette_size <- function(object_length, basic_color_palette){
  if(length(unique(object_length)) > length(basic_color_palette)){
    new_length <- length(unique(object_length)) - length(basic_color_palette)
    my_palette <- c(basic_color_palette, primary.colors(new_length))
  } else {
    my_palette <- basic_color_palette
  }
}
my_palellet<-c(basic_color_palette, primary.colors(8))

# 20190808_New Figures ----------------------------------------------------

setwd("monocle_msAggr/")
setwd("~/Desktop/10XGenomicsData/msAggr/msAggr_dim36/monocle_msAggr/")
msdata<-readRDS("20180517_msAggrdim36res2.5_UnsupClustMonocle.rds")
colnames(pData(msdata))
msdata@phenoData@data$invPseudo<--1*msdata@phenoData@data$Pseudotime
head(msdata@phenoData@data)

plot_cell_trajectory(msdata, color_by = "invPseudo")
plot_cell_trajectory(msdata, color_by = "orig.ident") + facet_wrap(~orig.ident)
plot_cell_trajectory(msdata, color_by = "Pseudotime") + facet_wrap(~orig.ident)


msSeurat<-readRDS("~/Desktop/10XGenomicsData/20180517_msAggrdim36res2.5_tsne.rds")
msSeurat<-SetAllIdent(msSeurat, id = "res.2")

TSNEPlot(msSeurat, colors.use = c("#2777AE", "#F67F0F", "#309F2B", "#D52728"), group.by = "orig.ident")

msdata<-readRDS("../20180517_msAggrdim36res2.5_UnsupClustMonocle_rev.rds")
plot_cell_trajectory(msdata, color_by = "Pseudotime")
plot_cell_trajectory(msdata, color_by = "State")
# will run msdata_state3<-orderCells(msdata, root_state = 3) on biowulf, mem=300g, cpus-per-task=32

msdata<-readRDS("/Users/heustonef/Desktop/20180517_msAggrdim36res2.5_UnsupClustMonocle_state3.rds")
plot_cell_trajectory(msdata, color_by = "Pseudotime")
plot_cell_trajectory(msdata, color_by = "State")

setwd("monocle_dim36Subsets/")
msdata<-readRDS("LSKSeuratSubset_UnsupClustMonocle.rds")
my_plot_cell_trajectory(msdata, color_by = "res.2", point_colors = my_palellet)

msdata<-readRDS("CMPSeuratSubset_UnsupClustMonocle.rds")
my_plot_cell_trajectory(msdata, color_by = "res.2", point_colors = my_palellet)

msdata<-readRDS("MEPSeuratSubset_UnsupClustMonocle.rds")
my_plot_cell_trajectory(msdata, color_by = "res.2", point_colors = my_palellet)

msdata<-readRDS("GMPSeuratSubset_UnsupClustMonocle.rds")
my_plot_cell_trajectory(msdata, color_by = "res.2", point_colors = my_palellet)


