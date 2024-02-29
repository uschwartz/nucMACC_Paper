MY_PWD=$1

mkdir $MY_PWD"/data"
cd $MY_PWD"/data"


wget -O genome.fa.gz https://ftp.ensembl.org/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
gunzip genome.fa.gz

wget -O genes.gtf.gz https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz
gunzip genes.gtf.gz


mkdir Bowtie2Index
bowtie2-build genome.fa Bowtie2Index/genome

cd ..
git clone https://github.com/uschwartz/nucMACC.git
