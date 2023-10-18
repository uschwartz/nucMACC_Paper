#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outDir="/Users/admin/Analysis/R001_nucMacc/manuscript_figures/Fig6/Initiation/"
params.inputDir="/Users/admin/Analysis/R001_nucMacc/manuscript_figures/Fig6/data/Initiation/*.bw"
params.Tss="/Users/admin/Analysis/R001_nucMacc/manuscript_figures/Fig6/data/*.bed"


channel.fromPath(params.inputDir).map{tuple(it.name.split('_')[0], it.name.split('_')[1],it )}set{input}
channel.fromPath(params.Tss).set{TSS}


process computeMatrix{
input:
 tuple val(gene),val(rep),file(input),file(TSS)

output:
  tuple val(gene), file ("*_Matrix.gz"), val(rep)

script:
"""
computeMatrix reference-point -S $input --samplesLabel $rep \
 -R $TSS \
 --referencePoint TSS \
 -o $gene"_Matrix.gz" \
 -b 500 -a 500 --smartLabels -p max
"""
}


process plot_Profile{
publishDir "${params.outDir}/${name}/", mode: 'copy'

input:
tuple val(name), file(matrix2plot), val(rep)


output:
tuple file("*_Profile.txt"), val(name), val(rep)
file("*.pdf")

script:
"""
plotProfile -m $matrix2plot -out 'DefaultProfile.pdf' \
     --outFileNameData $name"_Profile.txt"

"""
}

process plot_Heatmap{

publishDir "${params.outDir}/${name}/", mode: 'copy'

input:
  tuple val(name),file(matrix2plot),val(rep)

output:
  file("*.pdf")

script:
"""
plotHeatmap -m $matrix2plot --plotTitle $name --legendLocation upper-right \
 --colorMap Greens -out $name"_heatmap.pdf"
"""
}

workflow{

computeMatrix(input.combine(TSS.collect().toList()))

plot_Profile(computeMatrix.out)
plot_Heatmap(computeMatrix.out)
}
