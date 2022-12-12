# load GLS file locations


gls_files_all         <- list.files(path="data", pattern = ".lux$", recursive = T)
for(i in 1:length(gls_files_all)){
  gls <- as.character(read.table(paste("data",gls_files_all[i],sep="/"),skip=2,nrow=1,header = F)[1,3])
  if(gls=="DWX558") gls = "B118"
  
  gls <- data.frame(gls=gls, driftadj=T)
  gls$driftadj[!grepl("driftadj", gls_files_all[i])]<-F
  gls$file <- gls_files_all[i]
  
  if(i==1) gls.ids <- gls else gls.ids <- rbind(gls.ids, gls)
}
gls.ids <- gls.ids[order(gls.ids$driftadj, decreasing = T),]
save(gls.ids ,file="data/GLS files.RData")



deg_files_all         <- list.files(path="data", pattern = ".deg$", recursive = T)
for(i in 1:length(deg_files_all)){
  gls <- as.character(read.table(paste("data",deg_files_all[i],sep="/"),skip=2,nrow=1,header = F)[1,3])
  if(gls=="DWX558") gls = "B118"
  
  deg <- data.frame(gls=gls)
  deg$file <- deg_files_all[i]
  
  if(i==1) deg.ids <- deg else deg.ids <- rbind(deg.ids, deg)
}
save(deg.ids ,file="data/DEG files.RData")