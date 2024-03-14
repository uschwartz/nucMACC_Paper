#!/usr/bin/env Rscript 

path<-paste0("/home/",USER,"/nucMACC_Paper/data/nucMACC_H2AZ/TSSplots")
setwd(path)

# mono-nucleosomes
# 1. shH2A.Z

mono.shH2AZ<-read.delim("TSS_mono-nucs_shH2AZ.tsv")
val.mono.shH2AZ<-mono.shH2AZ[2:3,3:182]
row.names(val.mono.shH2AZ)<-mono.shH2AZ[2:3,1]
pos<-seq(-890,900,10)

pdf("TSS_mono-nucs_shH2AZ_R.pdf", width = 4, height = 3)
plot(pos,val.mono.shH2AZ[1,], type="l", lwd=2, col="#A3BBD7",
    main="mono-nucleosomes shH2A.Z", xlab="distance to TSS",
    ylab="norm. MNase fragment density", ylim=c(1300,4700))
lines(pos,val.mono.shH2AZ[2,], lwd=2, col="#761D78")
dev.off()

# 2. Whole chromatin

mono.WC<-read.delim("TSS_mono-nucs_WC.tsv")

val.mono.WC<-mono.WC[2:3,3:182]
row.names(val.mono.WC)<-mono.WC[2:3,1]

pdf("TSS_mono-nucs_WC_R.pdf", width = 4, height = 3)
plot(pos,val.mono.WC[1,], type="l", lwd=2, col="#A3BBD7",
    main="mono-nucleosomes WC", xlab="distance to TSS",
    ylab="norm. MNase fragment density", ylim=c(1300,4700))
lines(pos,val.mono.WC[2,], lwd=2, col="#761D78")
legend("bottomleft", col=c("#A3BBD7","#761D78"),
    legend=c("high MNase", "low MNase"),lwd=2, bty="n")
dev.off()


# sub-nucleosomes 
# 1. shH2A.Z

sub.shH2AZ<-read.delim("TSS_sub-nucs_shH2AZ.tsv")
val.sub.shH2AZ<-sub.shH2AZ[2:3,3:182]
row.names(val.sub.shH2AZ)<-sub.shH2AZ[2:3,1]
pos<-seq(-890,900,10)

pdf("TSS_sub-nucs_shH2AZ_R.pdf", width = 4, height = 3)
plot(pos,val.sub.shH2AZ[1,], type="l", lwd=2, col="#EE924F",
    main="sub-nucleosomes shH2A.Z", xlab="distance to TSS",
    ylab="norm. MNase fragment density", ylim=c(1500,9000))
lines(pos,val.sub.shH2AZ[2,], lwd=2, col="#751428")
dev.off()

# 2. Whole chromatin

sub.WC<-read.delim("TSS_sub-nucs_WC.tsv")
val.sub.WC<-sub.WC[2:3,3:182]
row.names(val.sub.WC)<-sub.WC[2:3,1]

pdf("TSS_sub-nucs_WC_R.pdf", width = 4, height = 3)
plot(pos,val.sub.WC[1,], type="l", lwd=2, col="#EE924F",
    main="sub-nucleosomes WC", xlab="distance to TSS",
    ylab="norm. MNase fragment density", ylim=c(1500,9000))
lines(pos,val.sub.WC[2,], lwd=2, col="#751428")
legend("bottomleft", col=c("#EE924F","#751428"),
    legend=c("high MNase", "low MNase"),lwd=2, bty="n")
dev.off()
