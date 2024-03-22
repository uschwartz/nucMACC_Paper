setwd("QC_out/")
pathSTD<-getwd()
library(stringr)



files.mono<-list.files("RUN/00_ALIGNMENT/monoNuc/")
files.sub<-list.files("RUN/00_ALIGNMENT/subNuc/")

##sample name
Sample_Name<-files.mono %>%
    str_split_i("_MNase_",1)

replicate<-files.mono %>%
    str_split_i("_",5)


MNase_U<-files.mono %>%
    str_split_i("_",3) %>%
    str_sub(1,-2)


df<-data.frame(Sample_Name=Sample_Name,
               replicate=replicate,
               path_mono=paste0(pathSTD,"/RUN/00_ALIGNMENT/monoNuc/",files.mono),
               path_sub=paste0(pathSTD,"/RUN/00_ALIGNMENT/subNuc/",files.sub),
               MNase_U)

write.csv(df, row.names = F,
          file="../data/samples_nucMACC.csv", quote = F)
