# Load Packages -----------------------------------------------------------



library(ggplot2)
library(RColorBrewer)


# Import Data -------------------------------------------------------------





setwd("/Users/mac-pro3/Analysis/Drosophila/Pausing/")

heat.val <- read.delim("values_heat_5scale.txt", skip=2)


NDR <- heat.val[1589:6939,1:315]
unstable <- heat.val[1:1588,1:315]
#convert NAs to 0
woNA.NDR<-NDR
woNA.NDR[is.na(woNA.NDR)]<-0
woNA.unstable<-unstable
woNA.unstable[is.na(woNA.unstable)]<-0

NDR.Median <- apply(woNA.NDR,MARGIN=2, function(x) median(x))
unstable.Median <- apply(woNA.unstable,MARGIN=2, function(x) median(x))

min.value <- min(NDR.Median,unstable.Median)
max.value <- max(NDR.Median,unstable.Median)

plot(seq(1:length(NDR.Median)),NDR.Median, type = "l", xaxt = "n", col = "blue", xlab = "Position", ylab = "Median",ylim=c(min.value,max.value))
lines(seq(1:length(unstable.Median)),unstable.Median, col ="green")
axis(1, at=c(0,50,70,265,315), labels=c("-500","TSS","unscaled","TES","500"))
legend(220, 2.7, legend=c("NDR", "unstable"),
       col=c("blue", "green"), lty=1:2, cex=0.8,box.lty=0)

# pausing Index NDR
NDR.un<-NDR[,51:70]
NDR.body<-NDR[,71:265]
pausing.NDR <- apply(NDR.un,MARGIN = 1,function(x) mean(x)) / apply(NDR.body,MARGIN = 1,function(x) mean(x))

# pausing Index unstable
unstable.un<-unstable[,51:70]
unstable.body<-unstable[,61:265]
pausing.unstable <- apply(unstable.un,MARGIN = 1,function(x) mean(x))/apply(unstable.body,MARGIN = 1,function(x) mean(x))

mean(pausing.NDR, na.rm = TRUE)
mean(pausing.unstable)

median(pausing.NDR, na.rm = TRUE)
median(pausing.unstable)


pausing <- c(pausing.NDR,pausing.unstable)
cat <- c(rep("NDR",length(pausing.NDR)),rep("unstable",length(pausing.unstable)))
df <- data.frame(pausing,cat)

ggplot(df, aes(x=cat, y=log(pausing), color=cat))+
  geom_violin()+
  scale_color_brewer(palette="Dark2")



ggplot(df, aes(x=cat, y=pausing, color=cat))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_brewer(palette="Dark2")+
  ylim(-1,10)




