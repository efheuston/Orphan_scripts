ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("ChIPseeker","TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db", "ggplot2")

## Install Libraries
ipak(packages)

# graphingSet=c("plotAnnoBar","plotAnnoPie", "vennpie", "upsetplot", "plotDistToTSS"),
# varblGraphs=c("covplot", "tagHeatmap", "plotAvgProf")

originalDIR<-getwd()

##Define R function input and 
autoChipSeeker<-function(useDir='',promoterRegion=c(-3000,3000), graphingSet=c("plotAnnoBar","plotAnnoPie", "vennpie", "upsetplot", "plotDistToTSS"),varblGraphs=c("covplot", "tagHeatmap", "plotAvgProf"))
{
  args <- commandArgs(trailingOnly=TRUE)  ## means when you type " Rscript --vanilla sillyScript.R iris.txt out.txt", it will start counting after "--" and pass the arguments to the vector "args" (separated by space)
  print("Running ChipSeekerRscript.r")
  try(setwd(useDir), silent = T)
  print(paste("Using ", getwd(), sep = ""))
  myfile<-list.files(pattern=".bed$")
  
  ##Define ChIPSeeker-required Variables`
  txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
  promoter<-getPromoters(TxDb=txdb, upstream=abs(promoterRegion[1]), downstream=abs(promoterRegion[2]))
  
  
  ##Pre-process for ChIPSeeker
  for (i in myfile){
    print(i)
    peak<-readPeakFile(i)
    tagMatrix<-getTagMatrix(peak, windows=promoter)
    histPA<-annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
    paFile<-paste(gsub(".bed", "", i), "-annoPks.txt", sep="")
    write.table(histPA@annoStat, file=paFile, sep="\t",row.names=FALSE)
    #tssFile<-paste(gsub(".bed", "", i), "-TSSdist.txt", sep="")
    #write.table(data.frame(histPA@anno)[,15], file=tssFile, sep="\t", row.names=FALSE, col.names = FALSE)
    
    # Plot functinos for ChIP profiles
    for(varbl in varblGraphs){
      namevec<-paste(gsub(".bed", "", i),"-", varbl, ".png", sep="")
      png(filename = paste(gsub(".bed", "", i),"-", varbl, ".png", sep=""), width = 1100, height = 600)
      print (varbl)
      pngFxn<-get(varbl)
      if(grepl("covplot", varbl)){
        try(x<-pngFxn(peak))
        try(plot(x))
      } else if(grepl("tagHeatmap",varbl)){
        try(pngFxn(tagMatrix, xlim=promoterRegion, color="black"))
      } else if(grepl("plotAvgProf", varbl)){
        try(x<-pngFxn(tagMatrix, xlim=promoterRegion, xlab="Genomic Region (5' -> 3')", ylab = "Read Count Frequency"))
        try(plot(x))
      }
      Sys.sleep(3)
      dev.off()
    }
    
    # Plot functions for annotated peaks
    for(fig in graphingSet){
      namevec<-paste(gsub(".bed", "", i),"-", fig, ".png", sep="")
      png(filename = paste(gsub(".bed", "", i),"-", fig, ".png", sep=""), width = 1100, height = 600)
      print (fig)
      pngFxn<-get(fig)
      if(grepl("plotAnnoPie|vennpie|upsetplot", fig)){
        pngFxn(histPA)
      } else {
        x<-pngFxn(histPA)
        plot(x)
        
      }
      Sys.sleep(3)
      dev.off()
    }
  }
}
  try(setwd(originalDIR))

# ipak function: install and load multiple R packages. Check to see if packages are installed. Install them if they are not, then load them into the R session.
#thanks to stevenworthington from GitHub (accessed 2016.08.05)
