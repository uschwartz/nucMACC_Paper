#!/usr/bin/env Rscript 

USER<-Sys.info()["user"]

path<-paste0("/home/",USER,"/nucMACC_Paper/data/MNase_Yeast/")
setwd(path)

# Load Packages -----------------------------------------------------------

library(ggpubr)
library(readxl)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(stringr)

# Import Data -------------------------------------------------------------

TSS<-read.delim("annotation/Yeast_TSS.bed", header = FALSE)
values <- read.delim("H4_IP/results/RUN/11_TSS_profile/subNucs/unstable_TSS/values_Heatmap.txt", skip = 2 ,header = TRUE)[,-301]
regions <- read.delim("H4_IP/results/RUN/11_TSS_profile/subNucs/unstable_TSS/sortedRegions_Heatmap.txt")
matrix<-read.delim("H4_IP/results/RUN/11_TSS_profile/subNucs/unstable_TSS/computeMatrix2plot_unstable.txt", skip = 1, header=FALSE)
prom.cat <- read_excel("41586_2021_3314_MOESM3_ESM.xlsx", sheet = 2)
expr <- read.delim("annotation/TPM_ProteinCoding.txt")
expr$log2TPM<-log2(expr$TPM+1)
expr<-expr[which(expr$TPM>1.5),]
gtf <- read.delim("../NETseq/data/Saccharomyces_cerevisiae.R64-1-1.96.gtf", header = FALSE, skip=5)

# Remove RLP genes ----------------------------------------------------

# Table form Rossi et al. 
RPL<-prom.cat[which(prom.cat$`Feature class Level 1`=="01_RP"),]
RPL<-RPL$`Systematic ID`

table(duplicated(RPL))
idx.RPL.dup<-ifelse(duplicated(RPL),FALSE,TRUE)
RPL<-RPL[idx.RPL.dup]

TSS.woRP<-TSS[!(TSS$V4 %in% RPL),]
expr.woRP<-expr[!(expr$gene %in% RPL),]

# Remove non-matching TSS and expression genes ----------------------------

expr.match<-expr.woRP[(expr.woRP$gene %in% TSS.woRP$V4),]
TSS.match<-TSS.woRP[(TSS.woRP$V4 %in% expr.match$gene),]
regions.match<-regions[(regions$name %in% expr.match$gene),]

# Define unstable TSS -----------------------------------------------------

len_sub <- 120/10
values<-values[,100:180]
idx<-apply(values, MARGIN = 1, function(x) sum(x))
idx<-ifelse(idx>=len_sub,TRUE,FALSE)
idx_stable<-ifelse(idx==TRUE,FALSE,TRUE)
unstable_TSS<-regions[idx,]
NDR_TSS<-regions[idx_stable,]
table(idx)
table(idx_stable)

write.table(unstable_TSS$name,file = "annotation/unstable_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(NDR_TSS$name,file = "annotation/NDR_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(unstable_TSS,file = "annotation/unstable_TSS.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(NDR_TSS,file = "annotation/NDR_TSS.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Assign Gene expression TSS ----------------------------------------------

unstabl.expr <- expr.match[(expr.match$gene %in% unstable_TSS$name),]
stabl.expr <- expr.match[(expr.match$gene %in% NDR_TSS$name),]

table(unstabl.expr$gene %in% stabl.expr$gene)

# create column with gene names in the expressed genes list ----------------

expr.sort<-rbind(unstabl.expr,stabl.expr)
expr.sort[1:nrow(unstabl.expr),4]<-"unstable"
expr.sort[nrow(unstabl.expr)+1:nrow(stabl.expr),4]<-"NDR"
colnames(expr.sort)[4]<-"cat"
table(expr.sort$cat)

# Violin plot of gene expression ------------------------------------------

sample_size = expr.sort %>% group_by(cat) %>% summarize(num=n())

expr.sort %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(cat, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=log2TPM, fill=cat)) +
  geom_violin(width=1.0) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Expression of genes according stability") +
  xlab("")+
  ylab("log2(TPM+1)")

shapiro.test(expr.sort$log2TPM[expr.sort$cat=="NDR"])
shapiro.test(expr.sort$log2TPM[expr.sort$cat=="unstable"])
test<-wilcox.test(expr.sort$log2TPM[expr.sort$cat=="NDR"], expr.sort$log2TPM[expr.sort$cat=="unstable"])
p.value <- test$p.value

write.table(p.value,"../NETseq/script/Figures/p_value_gene_expression.txt")
mean(stabl.expr$log2TPM)

# getting unstable expressed TSS ------------------------------------

TSS.unstable <- TSS[(TSS$V4 %in% unstabl.expr$gene),]

# getting id of stable expressed TSS --------------------------------------

TSS.NDR<-TSS[(TSS$V4 %in% stabl.expr$gene),]

# writing files ----------------------------------------------------------

write.table(TSS.unstable,file = "../NETseq/data/TSS_unstable.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TSS.NDR,file = "../NETseq/data/TSS_NDR.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(TSS.unstable$V4,file = "../NETseq/data/unstable_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(TSS.NDR$V4,file = "../NETseq/data/NDR_TSS.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# get all genes in RPL ----------------------------------------------------

gene.raw<-sapply(strsplit(as.character(gtf$V9), split = ";", fixed=T),
                       function(x) x[grep("gene_id",x)])

gene<-sapply(strsplit(gene.raw, split = " "),function(x) x[2])

gtf<-gtf[-(which(gene %in% RPL)),]

# get only NDR TSS -------------------------------------------------------

gtf.NDR<-gtf[(which(gene %in% TSS.NDR$V4)),]

# get only unstable TSS --------------------------------------------------

gtf.unstable<-gtf[(which(gene %in% TSS.unstable$V4)),]

write.table(gtf.NDR,file="../NETseq/data/NDR.gtf", quote=F,
            sep="\t", row.names = F, col.names = F)

write.table(gtf.unstable,file="../NETseq/data/unstable.gtf", quote=F,
            sep="\t", row.names = F, col.names = F)


