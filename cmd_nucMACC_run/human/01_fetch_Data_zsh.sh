
AnalysisDir="/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ"
OutDir=$AnalysisDir/data/

# download raw data of GEO accession
mkdir -p $OutDir

while IFS=$'\t' read -r col1 col2
do
    echo "get data for experiment:$col1|$col2"
    #load data from sra
    prefetch -O $OutDir"/"$col2 $col1
    ## convert format
    fastq-dump --split-files -O $OutDir"/"$col2 $OutDir"/"$col2"/"*

    #merge fastqs and gzip
    cat $OutDir"/"$col2"/"*_1.fastq | gzip >$OutDir"/"$col2"/"$col2"_1.fastq.gz"
    cat $OutDir"/"$col2"/"*_2.fastq | gzip >$OutDir"/"$col2"/"$col2"_2.fastq.gz"

    #remove uneeded files
    rm $OutDir"/"$col2"/"*.fastq
    rm -r $OutDir"/"$col2"/"$col1"/"

done < $AnalysisDir/script/SRA/sra_to_name.tsv
