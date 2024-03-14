#!/usr/bin/env Rscript 

USER<-Sys.info()["user"]

OUT<-paste0("/home/",USER,"/nucMACC_Paper/data/R001_nucMacc/manuscript_figures/Fig5/pairwise")

library(corrplot)

mx<-matrix(data = c(0.9,0.76,0.73, 0.53, 0.53, 0.47), ncol=1,
        dimnames = list(
        "row"=c("1.5U+100U","1.5U+25U","6.25U+100U","1.5U+6.25U",
                "6.25U+25U","1.5U+6.25U" ),
                "col"="1.5U+6.25U+25U+100U"))

pdf(file=paste0(OUT,"/corPlot.pdf"), height = 5)
corrplot(mx,col =RColorBrewer::brewer.pal(5,"Greens"),col.lim = c(0,1),
        addCoef.col = "white",tl.col="black"
        )
dev.off()

mx_sub<-matrix(data = c(0.92,0.78,0.19), ncol=1,
        dimnames = list(
                "row"=c("1.5U+100U","1.5U+25U","1.5U+6.25U" ),
                "col"="1.5U+6.25U+25U+100U"))

pdf(file=paste0(OUT,"/corPlot_sub.pdf"), height = 4)
corrplot(mx_sub,col =RColorBrewer::brewer.pal(5,"Greens"),col.lim = c(0,1),
        addCoef.col = "white",tl.col="black")
dev.off()

