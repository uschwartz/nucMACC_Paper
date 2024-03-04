library(pheatmap)

# Import Data -------------------------------------------------------------
setwd("/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/TSSplots/shH2AZ_unstable_PRETTY//")

## header
#read.delim("data/values_Heatmap_sub-nucMACC.txt", header = F,nrows = 1)

values <- read.delim("data/values_Heatmap_sub-nucMACC.txt", skip = 1 ,header = F)

regions <- read.delim("data/sortedRegions_Heatmap_sub_bed.txt")

tx2gene<-read.delim("~/Annotation/GRCh38/tx2gene_INFO.txt")

# extract best enrichment of sub-nucMACC scores -----------------------------------------------------

len_sub <- 120/10

(ncol(values)-6)/2

values.unstbl.shH2AZ<-values[,7:186]


# set range -100 / +10
sum.unstbl<-apply(values.unstbl.shH2AZ[,80:90], MARGIN = 1, function(x) sum(x))

#order by unstbl nuc
ord<-order(sum.unstbl,decreasing = F)

## remove duplicated genes
genes.sort<-regions[ord,]

mx<-match(genes.sort$name, tx2gene$tx.id)
genes.sort$name<-tx2gene$gene.id[mx]

genes.unique<-genes.sort[!duplicated(genes.sort$name),]

#heatmap color scale
colsubMACC<-colorRampPalette(c("#2166AC","#FFFFFF","#B2182B"))(100)


colnames(values.unstbl.shH2AZ)<-c("-900", rep(".",89),"TSS",
                                  rep(".",88),"+ 900")


##########

#plot together
values.unstbl.wc<-values[,187:ncol(values)]
colnames(values.unstbl.wc)<-c("-900", rep(".",89),"TSS",
                                  rep(".",88),"+ 900")


gaps_col=c(rep(180,4))


mx.merged<-cbind(values.unstbl.wc,values.unstbl.shH2AZ)


png("sub-nucMACC_heatmap_unique.png",res = 700, 
    width = 2000, height = 1800)
    
    pheatmap(mx.merged[ord[!duplicated(genes.sort$name)],], scale="none",
         breaks=seq(-2,2,length.out = 101),
         show_colnames = T, show_rownames = F,
         color = colsubMACC,gaps_col=gaps_col,
         cluster_rows = F, cluster_cols = F,
         na_col =  "#FFFFFF")
    
dev.off()


############ sorted TSS file

write.table(genes.unique,file = "data/unstable_TSS_MANsort.bed", sep="\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


########### sub-nucMACC #########

profile.wc<-apply(values.unstbl.wc[ord[!duplicated(genes.sort$name)],],2,mean)
profile.sh<-apply(values.unstbl.shH2AZ[ord[!duplicated(genes.sort$name)],],2,mean)


#######
pdf("profile_sub-nucMACC_unique.pdf", width = 4, height=3)
    plot(seq(-900,890,10), profile.sh, type="l", 
         xlab="distance from TSS", ylab = "sub-nucMACC score",
         lwd=2, main="shH2A.Z sub-nucs TSS", col="lightseagreen")
    lines(seq(-900,890,10), profile.wc, col="darkgrey",
          lwd=2)
    abline(h=0, lty=2)
    legend("bottomright", lwd=2, col=c("darkgrey","lightseagreen"),
           bty="n", legend = c("normal", "shH2A.Z"))
dev.off()
