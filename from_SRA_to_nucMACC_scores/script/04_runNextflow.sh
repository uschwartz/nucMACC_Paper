MY_PWD=$1
cd $MY_PWD

nextflow run nucMACC \
--analysis 'MNaseQC' \
--csvInput 'data/samples.csv' \
--outDir QC_out \
--genomeIdx $MY_PWD'/data/Bowtie2Index/genome' \
--genomeSize 14372600200 \
--publishBamFlt \
--TSS $MY_PWD'/data/genes.gtf' 
