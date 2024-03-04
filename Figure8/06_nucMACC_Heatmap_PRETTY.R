library(pheatmap)

# Import Data -------------------------------------------------------------
setwd("/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/TSSplots/shH2AZ_unstable_PRETTY//")

## header
read.delim("data/values_Heatmap_nucMACC.txt", header = F,nrows = 1)

values <- read.delim("data/values_Heatmap_nucMACC.txt", skip = 1 ,header = F)


# extract best enrichment of sub-nucMACC scores -----------------------------------------------------
values.unstbl.shH2AZ<-values[,7:186]

#heatmap color scale
colmonoMACC<-colorRampPalette(c("#984EA3","#FFFFFF","#1B9E77"))(100)

colnames(values.unstbl.shH2AZ)<-c("-900", rep(".",89),"TSS",
                                  rep(".",88),"+ 900")

##########

#plot together
values.unstbl.wc<-values[,187:ncol(values)]
colnames(values.unstbl.wc)<-c("-900", rep(".",89),"TSS",
                                  rep(".",88),"+ 900")


gaps_col=c(rep(180,4))


mx.merged<-cbind(values.unstbl.wc,values.unstbl.shH2AZ)


png("nucMACC_heatmap_unique.png",res = 700, 
    width = 2000, height = 1800)
    
    pheatmap(mx.merged, scale="none",
         breaks=seq(-1,1,length.out = 101),
         show_colnames = T, show_rownames = F,
         color = colmonoMACC,gaps_col=gaps_col,
         cluster_rows = F, cluster_cols = F,
         na_col =  "#FFFFFF")
    
dev.off()


########### nucMACC #########

profile.wc<-apply(values.unstbl.wc,2,mean)
profile.sh<-apply(values.unstbl.shH2AZ,2,mean)


#######
pdf("profile_nucMACC_unique.pdf", width = 4, height=3)
    plot(seq(-900,890,10), profile.sh, type="l", ylim=c(-0.25,0.1),
         xlab="distance from TSS", ylab = "nucMACC score",
         lwd=2, main="shH2A.Z nucMACCS TSS", col="lightseagreen")
    lines(seq(-900,890,10), profile.wc, col="darkgrey",
          lwd=2)
    abline(h=0, lty=2)
    legend("topleft", lwd=2, col=c("darkgrey","lightseagreen"),
           bty="n", legend = c("normal", "shH2A.Z"))
dev.off()


