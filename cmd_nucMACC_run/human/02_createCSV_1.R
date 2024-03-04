pathSTD<-"/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/"
library(stringr)
library(tidyverse)
setwd(pathSTD)
files<-list.files("data/", recursive = T) 



fwd<-grep("1.fastq.gz",files, value = T)
rev<-grep("2.fastq.gz",files, value = T)

name<-fwd %>% str_split_i("/",2) %>% 
    str_split_i("_R1",1)


MNase_U<-rep(15,length(name))
MNase_U[grep("_H_",name)]<-30


df<-data.frame(Sample_Name=name,
           path_fwdReads=paste0(pathSTD,"/data/",fwd),
           path_revReads=paste0(pathSTD,"/data/",rev),
           MNase_U)

h2az<-df %>% filter(grepl("A_H2AZ",Sample_Name))

write.csv(h2az, row.names = F,
          file="script/addData//samples_H2AZ.csv", quote = F)

input<-df %>% filter(grepl("A_Inp",Sample_Name))

write.csv(input, row.names = F,
          file="script/addData//samples_WC.csv", quote = F)


input.shh2az<-df %>% filter(grepl("A_shH2AZ",Sample_Name))

write.csv(input.shh2az, row.names = F,
          file="script/addData//samples_WC_shH2AZ.csv", quote = F)
