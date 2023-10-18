# Load Packages -----------------------------------------------------------



library(ggplot2)
library(RColorBrewer)


give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# Import Data -------------------------------------------------------------
setwd("/Users/mac-pro3/Analysis/NETseq/")

dark2 <- c(RColorBrewer::brewer.pal(8, "Dark2"))

heat.val <- read.delim("Pausing_Profile_sub/values_heat_5scale.txt", skip=2)
genes.list <- read.delim("Pausing_Profile_sub/sortedRegions_Heatmap.bed")
genes.NDR <- genes.list[which(genes.list$deepTools_group=="NDR_TSS_TTS.bed"),]
genes.unstable <- genes.list[which(genes.list$deepTools_group=="unstable_TSS_TTS.bed"),]

row.names(heat.val)<-as.character(genes.list$name)
colnames(heat.val)



NDR <- heat.val[which(genes.list$deepTools_group=="NDR_TSS_TTS.bed"),1:315]
unstable <- heat.val[which(genes.list$deepTools_group=="unstable_TSS_TTS.bed"),1:315]
#convert NAs to 0
woNA.NDR<-NDR
woNA.NDR[is.na(woNA.NDR)]<-0
woNA.unstable<-unstable
woNA.unstable[is.na(woNA.unstable)]<-0

NDR.Median <- apply(woNA.NDR,MARGIN=2, function(x) median(x))
unstable.Median <- apply(woNA.unstable,MARGIN=2, function(x) median(x))

min.value <- min(NDR.Median,unstable.Median)
max.value <- max(NDR.Median,unstable.Median)

pdf("Pausing_Profile_sub/Median_profile.pdf", width = 6, height = 5)
  plot(seq(1:length(NDR.Median)),NDR.Median, type = "l", xaxt = "n", lwd=2.5, col = dark2[1], xlab = "Position", ylab = "Median",ylim=c(min.value,max.value))
  lines(seq(1:length(unstable.Median)),unstable.Median, col =dark2[2], lwd=2.5)
  axis(1, at=c(0,50,65,265,315), labels=c("-500","TSS","unscaled","TES","500"))
  legend(220, 2.7, legend=c("NDR", "unstable"),
       col=c(dark2[1], dark2[2]), lty=c(1,1), cex=0.8,box.lty=0)
dev.off()

# pausing Index NDR
NDR.un<-NDR[,51:65]
NDR.body<-NDR[,66:265]
pausing.NDR <- apply(NDR.un,MARGIN = 1,function(x) median(x))/apply(NDR.body,MARGIN = 1,function(x) median(x))

# pausing Index NDR
unstable.un<-unstable[,51:65]
unstable.body<-unstable[,66:265]
pausing.unstable <- apply(unstable.un,MARGIN = 1,function(x) median(x))/apply(unstable.body,MARGIN = 1,function(x) median(x))

mean(pausing.NDR)
mean(pausing.unstable)

median(pausing.NDR)
median(pausing.unstable)


pausing <- c(pausing.NDR,pausing.unstable)
cat <- c(rep("NDR",length(pausing.NDR)),rep("unstable",length(pausing.unstable)))
df <- data.frame(pausing,cat)

pdf("Pausing_Profile_sub/Pausing_Violin.pdf", width = 6, height = 5)
ggplot(df, aes(x=cat, y=pausing, color=cat, fill=cat))+
  geom_violin()+
  theme_bw()+
  ylim(0,5)+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")
dev.off()

pdf("Pausing_Profile_sub/Pausing_Box.pdf", width = 4, height = 5)
ggplot(df, aes(x=cat, y=pausing, color=cat))+
  geom_boxplot()+
  scale_color_brewer(palette="Dark2")+
  #ylim(0,5)+
  ylab("pausing index")+
  theme_bw()
dev.off()

shapiro.test(df$pausing[df$cat=="NDR"])
shapiro.test(df$pausing[df$cat=="unstable"])

test<-wilcox.test(df$pausing[df$cat=="NDR"],df$pausing[df$cat=="unstable"])
test$p.value

write.table(test$p.value,file="Pausing_Profile_sub/P_Value.txt")
