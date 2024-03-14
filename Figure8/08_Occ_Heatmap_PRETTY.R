#!/usr/bin/env Rscript 

library(pheatmap)

# Import Data -------------------------------------------------------------

setwd("/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/TSSplots/shH2AZ_unstable_PRETTY/")

read.delim("data/values_Heatmap_occ.txt", header = F,nrows = 1)
values <- read.delim("data/values_Heatmap_occ.txt", skip = 1 ,header = F)

# extract best enrichment of sub-nucMACC scores -----------------------------

values.unstbl.shH2AZ<-values[,7:186]

# heatmap color scale
colNucs<-colorRampPalette(c(rep("#FFFFFF",1),"#D95F02"))(100)

colnames(values.unstbl.shH2AZ)<-c("-900", rep(".",89),"TSS",
                                  rep(".",88),"+ 900")

# plot together
values.unstbl.wc<-values[,187:(ncol(values)-180)]
colnames(values.unstbl.wc)<-c("-900", rep(".",89),"TSS",
                                  rep(".",88),"+ 900")

values.unstbl.H2AZ<-values[,(ncol(values)-179):ncol(values)]
colnames(values.unstbl.H2AZ)<-c("-900", rep(".",89),"TSS",
                              rep(".",88),"+ 900")

gaps_col=c(rep(180,4),rep(360,4))

mx.merged<-cbind(values.unstbl.wc,values.unstbl.shH2AZ, values.unstbl.H2AZ)

png("Occ_heatmap3_unique.png",res = 700, width = 2700, height = 1800)
    pheatmap(mx.merged, scale="none",
         breaks=seq(0,50000,length.out = 101),
         show_colnames = T, show_rownames = F,
         color = colNucs,gaps_col=gaps_col,
         cluster_rows = F, cluster_cols = F,
         na_col =  "#FFFFFF")
dev.off()

# OCC profile

profile.wc<-apply(values.unstbl.wc,2,mean)
profile.sh<-apply(values.unstbl.shH2AZ,2,mean)
profile.H2AZ<-apply(values.unstbl.H2AZ,2,mean)

pdf("profile_Occupancy_unique.pdf", width = 4, height=3)
    plot(seq(-900,890,10), profile.H2AZ, type="l", 
         xlab="distance from TSS", ylab = "nucleosome occupancy",
         lwd=2, main="shH2A.Z occupancy TSS", col="orange")
    lines(seq(-900,890,10), profile.wc, col="darkgrey",
          lwd=2)
    lines(seq(-900,890,10), profile.sh, col="lightseagreen",
          lwd=2)
    abline(h=0, lty=2)
    legend("topright", lwd=2, col=c("darkgrey","lightseagreen", "orange"),
           bty="n", legend = c("normal", "shH2A.Z", "H2A.Z"))
dev.off()


