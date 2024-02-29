MY_PWD=$1
cd $MY_PWD


nextflow run ~/00_scripts/nextflow/nucMACC \
--analysis 'nucMACC' \
--csvInput 'data/samples_nucMACC.csv' \
--outDir Run_out \
--genomeIdx $MY_PWD'/data/Bowtie2Index/genome' \
--genomeSize 14372600200 \
--genome $MY_PWD'/data/genome.fa' \
--bamEntry \
--TSS $MY_PWD'/data/genes.gtf' 
