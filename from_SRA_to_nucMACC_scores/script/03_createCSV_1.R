pathSTD<-getwd()
library(stringr)

files<-list.files("data/raw_data/")

fwd<-grep("1.fastq.gz",files, value = T)
rev<-grep("2.fastq.gz",files, value = T)

name<-fwd %>% str_sub(4,-12)

MNase_U<-name %>% str_split_i("_",3) %>%
    str_sub(1,-2)


df<-data.frame(Sample_Name=name,
           path_fwdReads=paste0(pathSTD,"/data/raw_data/",fwd),
           path_revReads=paste0(pathSTD,"/data/raw_data/",rev),
           MNase_U)

write.csv(df, row.names = F,
          file="data/samples.csv", quote = F)
