AnalysisDir=$1
OutDir=$AnalysisDir/data/H3

# download raw data of GEO accession GSE78984 at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78984
mkdir -p $OutDir

while IFS=$'\t' read -r col1 col2
do
    echo "get data for experiment:$col1|$col2"
    #load data from sra
    prefetch -O $OutDir"/"$col2 $col1


done < $AnalysisDir/script/H3chip/additionalFiles/SRA_index_file.tsv

cd $OutDir

for i in *_H3_*; do
echo $i

## convert format
fastq-dump --split-files -O $i $i"/"*".sra"

#gzip
gzip $i"/"*_1.fastq
gzip $i"/"*_2.fastq

#remove uneeded files
rm $OutDir"/"$i"/"*.sra

done

##change name
for i in *_H3_*; do
echo $i

#change name
mv $i"/"*1.fastq.gz $i"_1.fastq.gz"
mv $i"/"*2.fastq.gz $i"_2.fastq.gz"

done

##pool replicates
mkdir pooled

cat *_1.5U_1.fastq.gz >pooled/H3_ChIP_1.5U_1.fastq.gz
cat *_1.5U_2.fastq.gz >pooled/H3_ChIP_1.5U_2.fastq.gz
cat *_6.25U_1.fastq.gz >pooled/H3_ChIP_6.25U_1.fastq.gz
cat *_6.25U_2.fastq.gz >pooled/H3_ChIP_6.25U_2.fastq.gz

cat *_25U_1.fastq.gz >pooled/H3_ChIP_25U_1.fastq.gz
cat *_25U_2.fastq.gz >pooled/H3_ChIP_25U_2.fastq.gz
cat *_100U_1.fastq.gz >pooled/H3_ChIP_100U_1.fastq.gz
cat *_100U_2.fastq.gz >pooled/H3_ChIP_100U_2.fastq.gz
