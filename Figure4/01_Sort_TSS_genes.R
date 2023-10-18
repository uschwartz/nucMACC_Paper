
# Load Packages -----------------------------------------------------------

library(ggpubr)
library(readxl)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

# Import Data -------------------------------------------------------------


setwd("/Users/mac-pro3/Analysis/Drosophila/")
TSS<-read.delim("dm3.refGene_coding.bed", header = FALSE)
values <- read.delim("TSS_profile/values_Heatmap.txt", skip = 2 ,header = TRUE)[,-301]
regions <- read.delim("TSS_profile/sortedRegions_Heatmap.txt")
matrix <- read.delim("TSS_profile/computeMatrix2txt_unstable.txt.gz", skip = 1, header=FALSE)
expr <- read.delim("RPKM_protCod_expressed_genes.txt")
val <- read.delim("TSS_profile/values_Profile_unstable.txt")

# Define unstable TSS -----------------------------------------------------


len_sub <- 120/10

values<-values[,100:180]
idx<-apply(values, MARGIN = 1, function(x) sum(x))
idx<-ifelse(idx>=len_sub,TRUE,FALSE)
idx_stable<-ifelse(idx==TRUE,FALSE,TRUE)
unstable_TSS<-regions[idx,]
stable_TSS<-regions[idx_stable,]
table(idx)

write.table(unstable_TSS$name,file = "unstable_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(stable_TSS$name,file = "stable_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Assign Gene expression to unstable TSS ----------------------------------


idx2<-NULL

unstable_TSS$name<-as.character(unstable_TSS$name)
expr$refSeq<- as.character(expr$refSeq)



for(i in 1:nrow(unstable_TSS)){
  if(is.na(unstable_TSS[i,4])){
    break
  }
  for(n in 1:nrow(expr)){
    if(is.na(expr[n,2])){
      break
    }
    if(unstable_TSS[i,4] == expr[n,2]){
      idx2<-c(idx2,n)
      break
    }
  }
}

unst.expr <- expr[idx2,]
stabl.expr <- expr[-idx2,]


# create column with gene names in the epressed genes list ----------------

expr[idx2,4]<-"unstable"
expr[-idx2,4]<-"stable"
colnames(expr)[4]<-"cat"
expr[,5]<-rownames(expr)
colnames(expr)[5]<-"Gene"
unst.expr[,4]<-rownames(unst.expr)
colnames(unst.expr)[4]<-"Gene"
stabl.expr[,4]<-rownames(stabl.expr)
colnames(stabl.expr)[4]<-"Gene"


# Violin plot of gene expression ------------------------------------------


sample_size = expr %>% group_by(cat) %>% summarize(num=n())

expr %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(cat, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=genes.logRPKMs, fill=cat)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Expression of genes according stability") +
  xlab("")+
  stat_compare_means( label.x = 1.3, label.y = 13.5)


compare_means(log2TPM ~ cat, data = expr)


shapiro.test(unst.expr$genes.logRPKMs)
shapiro.test(stabl.expr$genes.logRPKMs)
wilcox.test(unst.expr$genes.logRPKMs,stabl.expr$genes.logRPKMs)





# type cast gene names to character ---------------------------------------


stabl.expr$Gene<- as.character(stabl.expr$Gene)
unst.expr$Gene<- as.character(unst.expr$Gene)
TSS$V4<- as.character(TSS$V4)
str(TSS)


# getting idx of unstable expressed TSS ------------------------------------


idx.uns.TSS<-NULL
for(i in 1:nrow(unst.expr)){
  if(is.na(unst.expr[i,4])){
    break
  }
  for(n in 1:nrow(TSS)){
    if(is.na(TSS[n,4])){
      break
    }
    if(unst.expr[i,4] == TSS[n,4]){
      idx.uns.TSS<-c(idx.uns.TSS,n)
      break
    }
  }
}


# getting id of stable expressed TSS --------------------------------------


idx.stab.TSS<-NULL
for(i in 1:nrow(stabl.expr)){
  if(is.na(stabl.expr[i,4])){
    break
  }
  for(n in 1:nrow(TSS)){
    if(is.na(TSS[n,4])){
      break
    }
    if(stabl.expr[i,4] == TSS[n,4]){
      idx.stab.TSS<-c(idx.stab.TSS,n)
      break
    }
  }
}




# Wrong Order of the if loops --> includes non-expressed genes ------------


for(i in 1:nrow(TSS)){
  if(is.na(TSS[i,4])){
    break
  }
  for(n in 1:nrow(expr)){
    if(is.na(expr[n,5])){
      break
    }
    if(TSS[i,4] == expr[n,5]){
      TSS[i,7] <- expr[n,4]
      break
    }
  }
}


TSS_unstable <- TSS[which(TSS$V7=="unstable"),1:6]
TSS_stable <- TSS[which(TSS$V7=="stable"|is.na(TSS$V7)),1:6]



# rigth commands ----------------------------------------------------------

TSS_unstable <- TSS[idx.uns.TSS,]
TSS_stable <- TSS[idx.stab.TSS,]


write.table(TSS_unstable,file = "TSS_unstable.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TSS_stable,file = "TSS_stable.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(TSS_unstable$V4,file = "unstable_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TSS_stable$V4,file = "stable_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
