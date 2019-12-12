
####################Load libraries####################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # if (length(new.pkg))
  #   install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, library, character.only = TRUE)
}

packages<-c("Seurat", "dplyr", "Matrix", "monocle", "reshape2", "colorRamps")

# Load libraries
ipak(packages)





my_findMarkerslist<-function(x){
  tsne_obj<-readRDS(x)
  
  my_name<-gsub(pattern = "20180419_msAggr-", replacement = "", x = x)
  my_name<-gsub(pattern = ".png", replacement = "", x = my_name)
  
  
  
  mouse.markers<-FindAllMarkers(tsne_obj, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
  
  mouse.differential_markers_top100 <- mouse.markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
  
  
  write.table(mouse.differential_markers_top100, file = paste("20180419_mouse_top100_differentiatlMarkers-",my_name,".txt", sep = ""), sep="\t", quote=F)
  mouse.differential_markers_all <- mouse.markers %>% group_by(cluster)
  write.table(mouse.differential_markers_all, file =paste("20180419_mouse_all_differentiatlMarkers-",my_name,".txt", sep=""), sep="\t", quote=F)
  
}


my_rdsList<-list.files(pattern = ".rds$")
lapply(my_rdsList, my_findMarkerslist)