#!/usr/bin/env Rscript 

path<-paste0("/home/",USER,"/nucMACC_Paper/data/NETseq/NETseq_Profile_unstable_NDR/")
setwd(path)

# Load Packages -----------------------------------------------------------

library(ggplot2)
library(RColorBrewer)

# experiment with the multiplier to find the perfect position -------------

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
}

# Import Data -------------------------------------------------------------

dark2 <- c(RColorBrewer::brewer.pal(8, "Dark2"))

heat.val <- read.delim("values_heat_5scale.txt", skip=2)
genes.list <- read.delim("sortedRegions_Heatmap.bed")

genes.NDR <- genes.list[which(genes.list$deepTools_group=="TSS_NDR.bed"),]
genes.unstable <- genes.list[which(genes.list$deepTools_group=="TSS_unstable.bed"),]

row.names(heat.val)<-as.character(genes.list$name)
colnames(heat.val)

NDR <- heat.val[which(genes.list$deepTools_group=="TSS_NDR.bed"),-c(1,2)]
unstable <- heat.val[which(genes.list$deepTools_group=="TSS_unstable.bed"),-c(1,2)]

# convert NAs to 0
woNA.NDR<-NDR
woNA.NDR[is.na(woNA.NDR)]<-0
woNA.unstable<-unstable
woNA.unstable[is.na(woNA.unstable)]<-0

NDR.Median <- apply(woNA.NDR,MARGIN=2, function(x) median(x))
unstable.Median <- apply(woNA.unstable,MARGIN=2, function(x) median(x))

min.value <- min(NDR.Median,unstable.Median)
max.value <- max(NDR.Median,unstable.Median)

pdf("Profile_300bp.pdf", width = 6, height = 5)
plot(seq(1:length(NDR.Median)),NDR.Median, type = "l", xaxt = "n", lwd=2.5, col = dark2[1], xlab = "Position", ylab = "Median",ylim=c(min.value,max.value))
lines(seq(1:length(unstable.Median)),unstable.Median, col =dark2[2], lwd=2.5)
axis(1, at=c(0,30,60), labels=c("-300","TSS","300"))
legend(0, 135, legend=c("NDR", "unstable"),
      col=c(dark2[1], dark2[2]), lty=c(1,1), cex=0.8,box.lty=0)
title("NETseq signal")
dev.off()

# filter 150bp after TSS
NDR.flt<-NDR[,31:46]
unstable.flt<-unstable[,31:46]

# calculate sum of 150bp region
NDR.flt.sum<-apply(NDR.flt, MARGIN = 1, sum)
unstable.flt.sum<-apply(unstable.flt, MARGIN = 1, sum)

# plot boxplot of 150bp sum according NDR and unstable 
data.NDR<-data.frame(NDR.flt.sum, rep("NDR", length(NDR.flt.sum)))
data.unst<-data.frame(unstable.flt.sum, rep("unstable", length(unstable.flt.sum)))

colnames(data.NDR)<-c("Sum", "Cat")
colnames(data.unst)<-c("Sum", "Cat")

data<-rbind(data.NDR, data.unst)

data$Log<-log2(data$Sum+1)
data$gene<-rownames(data)

write.table(data, "../data/NETseq_Value.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

pdf("NETseq_Sum_150bp.pdf", width = 4, height = 6)
ggplot(data, aes(x=Cat, y=Log, col=as.factor(Cat)))+
  geom_boxplot(aes(group=Cat), outlier.shape = NA)+
  theme_bw()+
  xlab("")+
  ylim(9,14)+
  ylab("sum(NETseq signal)")+
  ggtitle("NETseq Signal 150bp after TSS")+
  scale_color_manual(values=dark2, name="")
dev.off()

shapiro.test(data$Log[data$Cat=="unstable"])
shapiro.test(data$Log[data$Cat=="NDR"])

test<-wilcox.test(data$Log[data$Cat=="unstable"],data$Log[data$Cat=="NDR"])
test$p.value

write.table(test$p.value,file="P_Value_150bp.txt")
