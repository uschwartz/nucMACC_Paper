# Scripts to download all the required files and software used to run the nucMACC pipeline

## Requirements
*Docker (tested with v4.17.0)
*nextflow (tested with version 23.12.0-edge.5901)
*Bowtie2 (tested with version 2.3.4.3)
*sratoolkit (tested with version 3.0.10)
*R (tested with version 4.2.2)
*R package stringr (tested with version 1.5.0)

## Content
* 00_master_script.sh: Contains the command to subsequently execute the scripts. To start, copy the "script" folder into a new analysis folder and change the directory to the analysis folder. Once the required tools are installed the 00_master_script.sh can be executed to automate execution of the other scripts.
* 01_fetch_annotation.sh: downloads nucMACC pipeline, genome and gene annotation and generates bowtie2 index
* 02_fetch_Data_zsh.sh: downloads fastq files from sra
* 03_createCSV_1.R: creates input CSV file for MNaseQC workflow in nucMACC pipeline
* 04_runNextflow.sh: runs MNaseQC workflow in nucMACC pipeline
* 05_createCSV_2.R: creates input CSV file for nucMACC workflow
* 06_runNextflow_nucMACC.sh: runs nucMACC pipeline
