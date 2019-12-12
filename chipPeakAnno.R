# Attempt 1 ---------------------------------------------------------------


setwd("/Users/heustonef/Desktop/ATACpeaks/HomerPeaks/")
cmpe<-as.data.frame(read.table("ErySubset.bed", sep = "\t", header = FALSE, col.names = c("space", "start", "end"), blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))

head(cmpe)
rownames(cmpe)<-1:nrow(cmpe)
cmpe<-toGRanges(cmpe, format = "BED")
cmpe

cmp<-as.data.frame(read.table("CMP.bed", sep = "\t", header = FALSE, col.names = c("space", "start", "end"), blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
head(cmp)
rownames(cmp)<-1:nrow(cmp)
cmp<-toGRanges(cmp, format = "BED")

ery<-as.data.frame(read.table("ErySubset.bed", sep = "\t", header = FALSE, col.names = c("space", "start", "end"), blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
ery<-as.data.frame(read.table("ErySubset.bed", sep = "\t", header = FALSE, blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
head(ery)
ery$V1<-gsub("chr", "", ery$V1)
ery$id<-1:nrow(ery)
ery<-GRanges(seqnames = ery$V1, ranges=IRanges(start = ery$V2, end = ery$V3), names=ery$V4)
head(ery)

cmpe
cmp
ery

ol<-findOverlapsOfPeaks(cmp, cmpe, ery, connectedPeaks = "min", maxgap = 1000)
makeVennDiagram(ol, totalTest = 1e+2)



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")


library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(rtracklayer)
mmdb<-org.Mm.eg.db
# Attempt 2 ---------------------------------------------------------------



cmp<-as.data.frame(read.table("CMP.bed", sep = "\t", header = FALSE, blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
# head(cmp)
cmp$V1<-gsub("chr", "", cmp$V1)
cmp$id<-1:nrow(cmp)
cmp_avg<-mean(cmp$V3-cmp$V2)
cmp_avg
cmp<-GRanges(seqnames = cmp$V1, ranges=IRanges(start = cmp$V2, end = cmp$V3, names=paste("cmp", 1:nrow(cmp), sep = "")), seqinfo = cmp$V4)
head(cmp)



cmpe<-as.data.frame(read.table("ErySubset.bed", sep = "\t", header = FALSE, blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
# head(cmpe)
cmpe$V1<-gsub("chr", "", cmpe$V1)
cmpe$id<-1:nrow(cmpe)
cmpe_avg<-mean(cmpe$V3 - cmpe$V2)
cmpe_avg
cmpe<-GRanges(seqnames = cmpe$V1, ranges=IRanges(start = cmpe$V2, end = cmpe$V3, names=paste("cmpe", 1:nrow(cmpe), sep = "")), seqinfo = cmpe$V4)
head(cmpe)



ery<-as.data.frame(read.table("Ery.bed", sep = "\t", header = FALSE, blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
# head(ery)
ery$V1<-gsub("chr", "", ery$V1)
ery$id<-1:nrow(ery)
ery_avg<-mean(ery$V3-ery$V2)
ery_avg
ery<-GRanges(seqnames = ery$V1, ranges=IRanges(start = ery$V2, end = ery$V3, names=paste("ery", 1:nrow(ery), sep = "")), seqinfo = ery$V4)
head(ery)
max(ery_avg, cmp_avg, cmpe_avg)

# Calculate totalTest estimation
peakwidth_mean<-250
genm_size<-2.5e9
intragenic_size = 0.05 #percent of genome
upper_limit = genm_size*intragenic_size/peakwidth_mean
upper_limit
cmpe
cmp
ery

# ol<-findOverlapsOfPeaks(cmp, cmpe, ery, connectedPeaks = "min")
ol<-findOverlapsOfPeaks(cmp, cmpe, ery, connectedPeaks = "merge", maxgap = 1000)
# makeVennDiagram(ol, totalTest =10e5)
makeVennDiagram(ol, totalTest = upper_limit)
# install.packages(c("GenomicFeatures", "TxDb.Mmusculus.mm10.knownGene"))
install.packages("SuperExactTest")
library(SuperExactTest)
data("eqtls")


total = 2.5e9
num.expected.overlap = upper_limit
max_overlap = 15100 + 10700 + 500 + 8
length.gene.sets <- c(26942, 122922, 52244)
(length.gene.sets[2][1])
names(length.gene.sets) <-c("cmpe", "cmp", "ery")
length.gene.sets
p=sapply(0:max_overlap,function(i) dpsets(x = i, L = length.gene.sets, n=total, log.p = TRUE))

length

(length.gene.sets=sapply(cis.eqtls, length))
total=18196
(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
?dpsets
cmplength<-70600 + 15100 + 10700 + 22000
erylength<- 15100 + 10700 + 22000
cmpelength<-500 + 15100 + 10700
install.packages("extraDistr")
library(extraDistr)
dmvhyper(x = c('cmp' = 10700, 'ery' = 10700, 'cmpe' = 10700), n = c('cmp' =cmplength, 'ery' = erylength, 'cmpe' = cmpelength), k = (1e4))




peaks1 <- GRanges(seqnames=c("1", "2", "3", "4", "5", "6", 
                              "2", "6", "6", "6", "6", "5"),
                   ranges=IRanges(start=c(967654, 2010897, 2496704, 3075869, 
                                          3123260, 3857501, 201089, 1543200, 
                                          1557200, 1563000, 1569800, 167889600),
                                  end= c(967754, 2010997, 2496804, 3075969, 
                                         3123360, 3857601, 201089, 1555199,
                                         1560599, 1565199, 1573799, 167893599),
                                  names=paste("Site", 1:12, sep="")),
                  strand="+")

peaks2 <- GRanges(seqnames=c("1", "2", "3", "4", "5", "6", "1", "2", "3", 
                                     "4", "5", "6", "6", "6", "6", "6", "5"),
                          ranges=IRanges(start=c(967659, 2010898, 2496700, 
                                                 3075866, 3123260, 3857500, 
                                                 96765, 201089, 249670, 307586, 
                                                 312326, 385750, 1549800, 
                                                 1554400, 1565000, 1569400,
                                                 167888600), 
                                         end=c(967869, 2011108, 2496920, 
                                               3076166,3123470, 3857780, 
                                               96985, 201299, 249890, 307796, 
                                               312586, 385960, 1550599, 1560799,
                                               1565399, 1571199, 167888999), 
                                         names=paste("t", 1:17, sep="")),
                          strand=c("+", "+", "+", "+", "+", "+", "-", "-", "-", 
                                   "-", "-", "-", "+", "+", "+", "+", "+"))

ol <- findOverlapsOfPeaks(peaks1, peaks2, maxgap=1000)
peaklist <- ol$peaklist
peaks3 <- GRanges(seqnames=c("1", "2", "3", "4", "5", 
                             "6", "1", "2", "3", "4"),
                   ranges=IRanges(start=c(967859, 2010868, 2496500, 3075966,
                                          3123460, 3851500, 96865, 201189, 
                                          249600, 307386),
                                  end= c(967969, 2011908, 2496720, 3076166,
                                         3123470, 3857680, 96985, 201299, 
                                         249890, 307796),
                                  names=paste("p", 1:10, sep="")),
                  strand=c("+", "+", "+", "+", "+", 
                           "+", "-", "-", "-", "-"))

ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, maxgap=1000, 
                          connectedPeaks="min")
makeVennDiagram(ol, totalTest=1e+2)



# ChipSeeker --------------------------------------------------------------

# BiocManager::install("ChIPseeker")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("ChIPseeker")
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("/Users/heustonef/Desktop/ATACpeaks/HomerPeaks/")
# cmp<-as.data.frame(read.table("CMP.bed", sep = "\t", header = FALSE, blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
# cmplength<-nrow(cmp)
# head(cmplength)
# cmp<-GRanges(seqnames = cmp$V1, ranges=IRanges(start = cmp$V2, end = cmp$V3, names=paste("cmp", 1:nrow(cmp), sep = "")), seqinfo = cmp$V4)
# cmp$id<-paste("peak_", 1:cmplength, sep = "")
# cmp$fakescore <-1
# cmp
# 
# 
# cmpe<-as.data.frame(read.table("ErySubset.bed", sep = "\t", header = FALSE, blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
# cmpelength<-nrow(cmpe)
# cmpe$V1<-gsub("chr", "", cmpe$V1)
# cmpe<-GRanges(seqnames = cmpe$V1, ranges=IRanges(start = cmpe$V2, end = cmpe$V3, names=paste("cmpe", 1:nrow(cmpe), sep = "")), seqinfo = cmpe$V4)
# cmpe$id<-paste("peak_", 1:cmpelength, sep = "")
# cmpe$fakescore <-1
# cmpe
# 
# 
# 
# ery<-as.data.frame(read.table("Ery.bed", sep = "\t", header = FALSE, blank.lines.skip = TRUE,strip.white = TRUE, skipNul = TRUE))
# erylength<-nrow(ery)
# ery$V1<-gsub("chr", "", ery$V1)
# ery<-GRanges(seqnames = ery$V1, ranges=IRanges(start = ery$V2, end = ery$V3, names=paste("ery", 1:nrow(ery), sep = "")), seqinfo = ery$V4)
# ery$id<-paste("peak_", 1:erylength, sep = "")
# ery$fakescore<-1
# ery


files<-c("CMP.bed", "Ery.bed", "ErySubset.bed")
names(files)<-c("cmp", "ery", "cmpe")
files
cmp<-readPeakFile("CMP.bed")
ery<-readPeakFile("Ery.bed")
cmpe<-readPeakFile("ErySubset.bed")
# covplot(cmp, weightCol = "fakescore")
promoter<-getPromoters(TxDb<-txdb, upstream = 3000, downstream = 3000)
# cmpeTagMatrix<-getTagMatrix(cmpe, windows = promoter)
# tagHeatmap(cmpeTagMatrix, xlim = c(-1000, 1000), color = "black")
cmpePeakAnno<-annotatePeak(cmpe, TxDb = txdb, annoDb = "org.Mm.eg.db")
plotAnnoPie(cmpePeakAnno)
cmpPeakAnno<-annotatePeak(cmp, TxDb = txdb, annoDb = "org.Mm.eg.db")
plotAnnoPie(cmpPeakAnno)
eryPeakAnno<-annotatePeak(ery, TxDb = txdb, annoDb = "org.Mm.eg.db")
plotAnnoPie(eryPeakAnno)


tagMatrixList<-lapply(files, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim = c(-3000, 3000), TxDb = txdb)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)

cmpefiles<-c("CMP.bed", "Ery.bed", "ErySubset.bed")
names(cmpefiles)<-c("cmp", "ery", "cmpe")
cmpe_overlapEnrichment<- enrichPeakOverlap(queryPeak = files[[3]], targetPeak = unlist(files[1:2]), TxDb = txdb, pAdjustMethod = "BH", nShuffle = 1000, chainFile = NULL, verbose = TRUE)
cmpe_overlapEnrichment


cmpmkefiles<-c("CMP.bed", "Mk.bed", "MkErySub.bed")
names(cmpmkefiles)<-c("cmp", "mk", "cmpmke")
cmpmke_overlapEnrichment<- enrichPeakOverlap(queryPeak = files[[3]], targetPeak = unlist(files[1:2]), TxDb = txdb, pAdjustMethod = "BH", nShuffle = 1000, chainFile = NULL, verbose = TRUE)
cmpmke_overlapEnrichment





