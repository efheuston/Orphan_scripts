
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "colorRamps", "parallel", "future", "R.utils", "monocle")

# Load libraries
suppressMessages(ipak(packages))


# setwd("~/Desktop/bmDBA/MonocleRun/")
# bmdba<-readRDS("20190830_all_ccadim20-FINAL.rds")
# TSNEPlot(bmdba)
# unique(bmdba@meta.data$orig.ident)
# length(unique(bmdba@meta.data$orig.ident))
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


color_pallet = c("gray", "gray", "gray", "lightskyblue", "gray","gray","gray","#cb4bbe", "gray","gray","gray","#53ba48", "deeppink4", "#dcca44", "#502d70")
TSNEPlot(bmdba, colors.use = color_pallet, group.by = "orig.ident")
png(filename = "20190830_all_ccadim20_tsnebysmplType.png", width = 800, height = 800)
TSNEPlot(bmdba, colors.use = color_pallet, group.by = "orig.ident")
dev.off()

setwd("/Users/heustonef/Desktop/10XGenomicsData/msAggr/msAggr_dim36/monocle_msAggr/")
msaggr<-readRDS("20180517_msAggrdim36res2.5_UnsupClustMonocle.rds")
png(filename = "Monocle_matchScanpyColors.png", width = 800, height = 800)
my_plot_cell_trajectory(msaggr, color_by = "orig.ident", point_colors = c("#4076AA", "#E68636", "#539D3E", "#C43832"))
dev.off()




# 20190925 ----------------------------------------------------------------

setwd("/Volumes/Labs/Bodine/ToDeena/20190821_CellRanger_BPsamples/20190924_SeuratMonocle/")
setwd("curatedDBA/")
dba<-readRDS("20190924_ConDBA_ccadim20-FINAL_UnsupClustMonocle.rds")


my_plot_cell_trajectory(dba, cell_size = 1, point_colors = basic_color_palette)
adjust_palette_size <- function(object_length, basic_color_palette){
  if(length(unique(object_length)) > length(basic_color_palette)){
    new_length <- length(unique(object_length)) - length(basic_color_palette)
    my_palette <- c(basic_color_palette, primary.colors(new_length))
  } else {
    my_palette <- basic_color_palette
  }
}
length(unique(dba$res.2))
my_palette <- c(basic_color_palette, primary.colors(4))
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2')



png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.2.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.2FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~res.2)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.2FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~orig.ident)
dev.off()

png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.1.5.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.1.5FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~res.1.5)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.1.5FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~orig.ident)
dev.off()

png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.1.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.1FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~res.1)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.1FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~orig.ident)
dev.off()

png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.0.6.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.0.6FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~res.0.6)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_trajectory-res.0.6FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~orig.ident)
dev.off()



png("20190924_ConDBA_ccadim20-FINAL_clstr-res.2.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.2FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~res.2)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.2FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~orig.ident)
dev.off()

png("20190924_ConDBA_ccadim20-FINAL_clstr-res.1.5.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.1.5FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~res.1.5)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.1.5FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~orig.ident)
dev.off()

png("20190924_ConDBA_ccadim20-FINAL_clstr-res.1.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.1FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~res.1)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.1FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~orig.ident)
dev.off()

png("20190924_ConDBA_ccadim20-FINAL_clstr-res.0.6.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6')
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.0.6FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~res.0.6)
dev.off()
png("20190924_ConDBA_ccadim20-FINAL_clstr-res.0.6FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~orig.ident)
dev.off()



setwd("../BP03vsBP05//")
dba<-readRDS("20190924_BP03BP05_ccadim20-FINAL_UnsupClustMonocle.rds")


png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.2.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.2FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~res.2)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.2FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~orig.ident)
dev.off()

png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.1.5.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.1.5FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~res.1.5)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.1.5FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~orig.ident)
dev.off()

png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.1.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.1FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~res.1)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.1FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~orig.ident)
dev.off()

png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.0.6.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.0.6FACET.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~res.0.6)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_trajectory-res.0.6FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_trajectory(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~orig.ident)
dev.off()



png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.2.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.2FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~res.2)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.2FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.2') + facet_wrap(~orig.ident)
dev.off()

png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.1.5.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.1.5FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~res.1.5)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.1.5FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1.5') + facet_wrap(~orig.ident)
dev.off()

png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.1.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.1FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~res.1)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.1FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.1') + facet_wrap(~orig.ident)
dev.off()

png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.0.6.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6')
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.0.6FACET.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~res.0.6)
dev.off()
png("20190924_BP03BP05_ccadim20-FINAL_clstr-res.0.6FACETorig.ident.png", height = 800, width = 800)
my_plot_cell_clusters(dba, cell_size = 1, point_colors = my_palette, color_by = 'res.0.6') + facet_wrap(~orig.ident)
dev.off()


# 20190930_Curated analysis -----------------------------------------------

setwd("/Volumes/Labs/Bodine/ToDeena/20190821_CellRanger_BPsamples/20190924_SeuratMonocle/curatedDBA/")
bmdba<-readRDS("20190924_ConDBA_ccadim20-FINAL.rds")

bmdba<-SetAllIdent(bmdba, id = "res.2")
color_pallet = c("gray", "steelblue", "gray", "dodgerblue1", "gray", "violet", "blue3", "deeppink", "firebrick")
TSNEPlot(bmdba, group.by = "orig.ident", colors.use = color_pallet)
TSNEPlot(bmdba, colors.use = c(basic_color_palette, primary.colors(4)), do.label = TRUE)

png("20190924_ConDBA_ccadim20_tsne-bySampleType.png", height = 1600, width = 1600)
TSNEPlot(bmdba, group.by = "orig.ident", colors.use = color_pallet, pt.size = 3) + ggtitle("TSNE plot: controls = Gray; RPLs = Blue, RPSs = Red \n note clusters 3, 8, parts of 11 and 24")
dev.off()
