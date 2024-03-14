#!/usr/bin/env Rscript 

# set working directory
setwd("QC_out/")

# save working directory as pathSTD
pathSTD<-getwd()

# list mono-nucleosome and sub-nucleosome BAM files
files.mono<-list.files("RUN/00_ALIGNMENT/monoNuc/")
files.sub<-list.files("RUN/00_ALIGNMENT/subNuc/")

# Change sample name by splitting at "_MNase_" and taking the first part
Sample_Name <- sapply(strsplit(files.mono, "_MNase_"), `[`, 1)

# Extract replicate number by splitting at "_" and taking the sixth element
replicate <- sapply(strsplit(files.mono, "_"), `[`, 6)

# Extract MNase units by splitting at "_" and taking the fourth element, then removing the last character
MNase_U <- sapply(strsplit(files.mono, "_"), function(x) substr(x[4], 1, nchar(x[4])-1))

# create data frame with sample name and corresponding file names for mono- and sub- nucleosomes
df<-data.frame(Sample_Name=Sample_Name,
               replicate=replicate,
           path_mono=paste0(pathSTD,"/RUN/00_ALIGNMENT/monoNuc/",files.mono),
           path_sub=paste0(pathSTD,"/RUN/00_ALIGNMENT/subNuc/",files.sub),
           MNase_U)

# write .csv file with sample names and paths
write.csv(df, row.names = F,
          file = paste0(pathSTD,"/data/samples_nucMACC.csv", quote = F)
