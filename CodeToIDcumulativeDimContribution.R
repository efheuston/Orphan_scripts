mep<-readRDS("MEPSeuratSubset_tsne.rds")
slotNames(mep)
slotNames(mep@dr$pca)

head(mep@dr$pca@sdev)
george<-mep@dr$pca@sdev
length(george)
library(Seurat)
PCElbowPlot(object = mep, num.pc = 36)
george

george<-(mep@dr$pca)


george
# Get cell embeddings for the first 2 PCs
Seurat::GetCellEmbeddings(object = object, reduction.type = "pca", dims.use = 1:2)

george<-mep

# george<-my_RunPCA(george, pc.genes = george@var.genes, pcs.compute = 36, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

george
george@dr$pca@jackstraw
JackStrawPlot(george, PCs = 1:36)
PCElbowPlot(george)



george<-RunPCA(george, pc.genes = george@var.genes, pcs.compute = 50, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCElbowPlot(george, num.pc = 50)
george<-JackStraw(george, num.pc = 50)

JackStrawPlot(george, PCs = 1:50)
slotNames(george@dr$pca@jackstraw)
head(george@dr$pca@jackstraw@overall.p.values)
