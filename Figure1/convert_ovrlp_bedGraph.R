#!/usr/bin/env Rscript 

args<-commandArgs(TRUE)

library(rtracklayer)

input<-import.bedGraph(args[1])

disjoint<-disjoin(input)

Overl<-findOverlaps(disjoint,input)
score.idx.list<-split(Overl@to, Overl@from )

score.disjoin<-sapply(score.idx.list, function(x) mean(input$score[x]))

disjoint$score<-score.disjoin
# sort 
seqlevels(disjoint)<-sort(seqlevels(disjoint))
disjoint.sort<-sort(disjoint)
#export
export.bedGraph(disjoint.sort,args[2])
