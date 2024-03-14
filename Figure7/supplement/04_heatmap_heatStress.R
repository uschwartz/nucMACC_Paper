#!/usr/bin/env Rscript 

setwd("~/Analysis/R001_nucMacc/manuscript_figures/Fig6/heatShock_SAGA//")

library(stringr)

profiles<-list.files(pattern=".txt",recursive = T)

hyper_n<-list()
hypo_n<-list()
hyper_hs<-list()
hypo_hs<-list()

for(i in profiles){
    #i="Bdf1/Bdf1_Profile.txt"
    dat.table<-read.delim(i)
    mx<-dat.table[-1,3:102]
    mx.n<-mx[1:2,]
    row.names(mx.n) <- dat.table$X[2:3]
    mx.hs<-mx[3:4,]
    row.names(mx.hs) <- dat.table$X[4:5]
    hyper_n[[str_split_i(i,"/",1)]]<-mx.n["hyper_TSS.bed",]
    hypo_n[[str_split_i(i,"/",1)]]<-mx.n["hypo_TSS.bed",]
    hyper_hs[[str_split_i(i,"/",1)]]<-mx.hs["hyper_TSS.bed",]
    hypo_hs[[str_split_i(i,"/",1)]]<-mx.hs["hypo_TSS.bed",]
}

mx.hyper_n<-do.call(rbind,hyper_n)
mx.hypo_n<-do.call(rbind,hypo_n)
mx.hyper_hs<-do.call(rbind,hyper_hs)
mx.hypo_hs<-do.call(rbind,hypo_hs)

colnames(mx.hyper_n)<-rep("",100)
colnames(mx.hyper_n)[c(1,50,100)]<-c("-0.5kb","TSS","0.5kb")
colnames(mx.hypo_n)<-rep("",100)
colnames(mx.hypo_n)[c(1,50,100)]<-c("-0.5kb","TSS","0.5kb")

colnames(mx.hyper_hs)<-rep("",100)
colnames(mx.hyper_hs)[c(1,50,100)]<-c("-0.5kb","TSS","0.5kb")
colnames(mx.hypo_hs)<-rep("",100)
colnames(mx.hypo_hs)[c(1,50,100)]<-c("-0.5kb","TSS","0.5kb")

pdf("hypo_plusOne_hs.pdf", width = 4, height = 1.5)
pheatmap::pheatmap(mx.hypo_hs, scale = "none", cluster_rows = F, cluster_cols = F,
                    breaks = seq(1,15,length.out = 101),border_color = NA)
dev.off()

pdf("hyper_plusOne_hs.pdf", width = 4, height = 1.5)
pheatmap::pheatmap(mx.hyper_hs, scale = "none", cluster_rows = F, cluster_cols = F,
                    breaks = seq(1,15,length.out = 101),border_color = NA)
dev.off()

pdf("hypo_plusOne_n.pdf", width = 4, height = 1.5)
pheatmap::pheatmap(mx.hypo_n, scale = "none", cluster_rows = F, cluster_cols = F,
                breaks = seq(1,15,length.out = 101),border_color = NA)
dev.off()

pdf("hyper_plusOne_n.pdf", width = 4, height = 1.5)
pheatmap::pheatmap(mx.hyper_n, scale = "none", cluster_rows = F, cluster_cols = F,
                    breaks = seq(1,15,length.out = 101),border_color = NA)
dev.off()
