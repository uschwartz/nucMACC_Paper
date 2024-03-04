

# Import Data -------------------------------------------------------------
setwd("/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/TSSplots/shH2AZ/data/")
values <- read.delim("values_Heatmap_sub_bed.txt", skip = 1 ,header = F)
regions <- read.delim("sortedRegions_Heatmap_sub_bed.txt")



# Define unstable TSS -----------------------------------------------------

len_sub <- 120/10

(ncol(values)-6)/2

values.unstbl<-values[,187:ncol(values)]

# set range -150 / +50
idx<-apply(values.unstbl[,75:95], MARGIN = 1, function(x) sum(x))
idx<-ifelse(idx>=len_sub,TRUE,FALSE)
unstable_TSS<-regions[idx,]


write.table(unstable_TSS,file = "unstable_TSS_shH2AZ.bed", sep="\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
