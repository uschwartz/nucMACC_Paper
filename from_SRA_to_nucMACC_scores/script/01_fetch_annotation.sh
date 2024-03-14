#!/usr/bin/env bash

MY_PWD=$1

mkdir $MY_PWD"/data"
cd $MY_PWD"/data"

# Download the reference genome for D. melanogaster
wget -O genome.fa.gz https://ftp.ensembl.org/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz

# Unzip genome file
gunzip genome.fa.gz

# Download annotation file for D. melanogaster
wget -O genes.gtf.gz https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz
gunzip genes.gtf.gz

# Create bowtie2 index from the reference genome
mkdir Bowtie2Index
bowtie2-build genome.fa Bowtie2Index/genome

# Get the nucMACC pipeline
cd ..
git clone https://github.com/uschwartz/nucMACC.git
