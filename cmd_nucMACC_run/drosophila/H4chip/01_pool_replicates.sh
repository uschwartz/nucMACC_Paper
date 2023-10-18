AnalysisDir=/Users/admin/Analysis/R001_nucMacc

cd $AnalysisDir/data/H4/

mkdir pooled

cat *_1.5U_1.fastq.gz >pooled/H4_ChIP_1.5U_1.fastq.gz
cat *_1.5U_2.fastq.gz >pooled/H4_ChIP_1.5U_2.fastq.gz
cat *_6.25U_1.fastq.gz >pooled/H4_ChIP_6.25U_1.fastq.gz
cat *_6.25U_2.fastq.gz >pooled/H4_ChIP_6.25U_2.fastq.gz

cat *_25U_1.fastq.gz >pooled/H4_ChIP_25U_1.fastq.gz
cat *_25U_2.fastq.gz >pooled/H4_ChIP_25U_2.fastq.gz
cat *_100U_1.fastq.gz >pooled/H4_ChIP_100U_1.fastq.gz
cat *_100U_2.fastq.gz >pooled/H4_ChIP_100U_2.fastq.gz


## check fastqc

mkdir fastqc

fastqc -o fastqc *.gz

