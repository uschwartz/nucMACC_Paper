#!/usr/bin/env Rscript 

library(ggplot2)
library(ggpubr)

# Import Data -------------------------------------------------------------

setwd("/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/TSSplots/")

df<-data.frame(numb.unstableTSS=c(61,680),
           condition=c("normal", "shH2A.Z"))

pdf("Barplot_enrich_unstable_TSS_unique.pdf", width=3, height = 4)
    ggbarplot(df, x = "condition", y = "numb.unstableTSS",
          fill = "condition", label = TRUE,color = "condition",
          palette = c("darkgrey","lightseagreen"))
dev.off()
