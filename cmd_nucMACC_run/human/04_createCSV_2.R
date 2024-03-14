#!/usr/bin/env Rscript 

USER<-Sys.info()["user"]

path<-paste0("/home/",USER,"/nucMACC_Paper")
pathSTD<-paste0(path,"/QC_out/")

library(stringr)

setwd(pathSTD)
files.mono<-list.files("RUN/00_ALIGNMENT/monoNuc/")
files.sub<-list.files("RUN/00_ALIGNMENT/subNuc/")

# Get sample name
files.mono <- c("example_r1_data", "sample_r2_info")  # Example input

# Extract Sample_Name by splitting at "_r" and taking the first part
Sample_Name <- sapply(strsplit(files.mono, "_r"), `[`, 1)

# Extract replicate by splitting at "_" and taking the fifth element
replicate <- sapply(strsplit(files.mono, "_"), `[`, 5

MNase_U<-rep(15,length(Sample_Name))
MNase_U[grep("_H",Sample_Name)]<-30

df<-data.frame(Sample_Name=Sample_Name,
               replicate=replicate,
           path_mono=paste0(pathSTD,"RUN/00_ALIGNMENT/monoNuc/",files.mono),
           path_sub=paste0(pathSTD,"RUN/00_ALIGNMENT/subNuc/",files.sub),
           MNase_U)

write.csv(df, row.names = F,
          file=paste0($path,"/data/samples_WC_nucMACC.csv"), quote = F)
