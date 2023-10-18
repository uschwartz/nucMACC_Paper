setwd("~/Analysis/R001_nucMacc/manuscript_figures/Fig6/Initiation/")

library(stringr)


dirs<-list.dirs()
profiles<-list.files(pattern=".txt",recursive = T)


hyper<-list()
hypo<-list()

for(i in profiles){
    #i="Bdf1/Bdf1_Profile.txt"
    dat.table<-read.delim(i)
    mx<-dat.table[-1,-(1:2)]
    row.names(mx) <- dat.table$X[2:3]
    scale.fac<-apply(mx,1,median)
    #mx.norm<-log2(mx+0.001)
    mx.norm<-mx
    hyper[[str_split_i(i,"/",1)]]<-mx.norm["hyper_TSS.bed",]
    hypo[[str_split_i(i,"/",1)]]<-mx.norm["hypo_TSS.bed",]
}


mx.hyper<-do.call(rbind,hyper)
mx.hypo<-do.call(rbind,hypo)

colnames(mx.hyper)<-rep("",100)
colnames(mx.hyper)[c(1,50,100)]<-c("-0.5kb","TSS","0.5kb")
colnames(mx.hypo)<-rep("",100)
colnames(mx.hypo)[c(1,50,100)]<-c("-0.5kb","TSS","0.5kb")

pdf("hypo_plusOne.pdf", width = 4, height = 3)
    pheatmap::pheatmap(mx.hypo, scale = "none", cluster_rows = F, cluster_cols = F,
                   breaks = seq(1,12,length.out = 101), border_color = NA)
dev.off()

pdf("hyper_plusOne.pdf", width = 4, height = 3)
    pheatmap::pheatmap(mx.hyper, scale = "none", cluster_rows = F, cluster_cols = F,
                   breaks = seq(1,12,length.out = 101), border_color = NA)
dev.off()