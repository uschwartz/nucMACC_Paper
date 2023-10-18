setwd("/Volumes/admin//Analysis/R001_nucMacc/manuscript_figures/Fig1/data/expression/")

## load gtf
gtf<-read.delim("dm3.refGene.gtf", header = F)

gtf.tx<-gtf[which(gtf$V3 == "transcript"),]


#gene ids
id.gtf<-sapply(strsplit(as.character(gtf.tx$V9), fixed = T, split = ";"),
               function(x) sapply(
                   strsplit(grep("transcript_id", x, value = T), split = " ", fixed = T),
                   function(y) y[3]    
               ) 
)

## take protein coding
gtf.tx.c<-gtf.tx[grep("NM_",id.gtf),]

#gene names
genes.gtf<-sapply(strsplit(as.character(gtf.tx.c$V9), fixed = T, split = ";"),
                  function(x) sapply(
                      strsplit(grep("gene_id", x, value = T), split = " ", fixed = T),
                      function(y) y[2]    
                  ) 
)

gtf.tx.c.u<-gtf.tx.c[!duplicated(genes.gtf),]

write.table(cbind(gtf.tx.c.u[,c(1,4,5)],
                  genes.gtf[!duplicated(genes.gtf)],
                  gtf.tx.c.u[,c(6,7)]), 
            file="dm3.refGene_coding.bed", 
            col.names = F, row.names = F,
            quote=F, sep="\t")

