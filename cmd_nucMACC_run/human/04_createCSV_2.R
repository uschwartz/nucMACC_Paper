pathSTD<-"/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/WC_QC_out/"
library(stringr)

setwd(pathSTD)
files.mono<-list.files("RUN/00_ALIGNMENT/monoNuc/")
files.sub<-list.files("RUN/00_ALIGNMENT/subNuc/")

##sample name
Sample_Name<-files.mono %>% 
    str_split_i("_r",1)

replicate<-files.mono %>% 
    str_split_i("_",4)


MNase_U<-rep(15,length(Sample_Name))
MNase_U[grep("_H",Sample_Name)]<-30

df<-data.frame(Sample_Name=Sample_Name,
               replicate=replicate,
           path_mono=paste0(pathSTD,"RUN/00_ALIGNMENT/monoNuc/",files.mono),
           path_sub=paste0(pathSTD,"RUN/00_ALIGNMENT/subNuc/",files.sub),
           MNase_U)

write.csv(df, row.names = F,
          file="../script/addData/samples_WC_nucMACC.csv", quote = F)
