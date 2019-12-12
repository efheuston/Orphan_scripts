# Install Packages --------------------------------------------------------

# install.packages("fossil")
# install.packages("clues")
library(Seurat)
library(dplyr)
library(Matrix)
library(clues)
library(fossil)
library(monocle)
# library(MyPlotTrajectoryPackage)

# Color Pallet ------------------------------------------------------------

color_pallet<-c("#cb4bbe",
                 "lightskyblue",
                 "#dcca44",
                 "#502d70",
                 "green4",
                 "red",
                 "springgreen",
                 "#5d2539",
                 "#cfd0a0",
                 "blue",
                 "#d2883b",
                 "maroon2",
                 "#898938",
                 "#c98ac2",
                 "yellow",
                 "#c4c0cc",
                 "#7d3d23",
                 "#00a5ff",
                 "#d68d7a",
                 "#a2c1a3")


# Compare DBA_dim14res1 vs DBA_dim17res1 ----------------------------------

setwd("DBA")
load("tsneIterations/nope/dba_cca_tsne-dim17_res1.Robj")
dba_dim17res1<-iterated_tsne
remove(iterated_tsne)
load("dba_cca_tsne-dim14_res1.Robj")
dba_dim14res1<-iterated_tsne
remove(iterated_tsne)

head(dba_dim14res1@meta.data)

dim14res1_cluster<-dba_dim14res1@meta.data$res.1
dim17res1_cluster<-dba_dim17res1@meta.data$res.1


george<-adjustedRand(as.numeric(dim17res1_cluster), as.numeric(dim14res1_cluster))
max(dim14res1_cluster)
max(as.numeric(dim14res1_cluster))
max(as.numeric(dim17res1_cluster))

dim14res1_classes<-dba_dim14res1@meta.data$group
dim17res1_classes<-dba_dim17res1@meta.data$group

adjustedRand(as.numeric(dim14res1_cluster), dim14res1_classes)

library(mclust)
george<-adjustedRandIndex(dim14res1_cluster, dim14res1_classes)
adjustedRandIndex(dim17res1_cluster, dim17res1_classes)
load("dba_cca_tsne-dim14_res1-v7_mdim7_monocle/dba_cca_tsne-dim14_res1-v7_mdim7_dpFeature.Robj")
dba_monocle<-cell_pop_dpFeature
remove(cell_pop_dpFeature)








