#!/usr/bin/env Rscript 

# Get working directory
pathSTD<-getwd()

# Load libraries
library(stringr)

# List fastq files
files<-list.files("data/raw_data/")

# Define forward/reverse reads
fwd<-grep("1.fastq.gz",files, value = T)
rev<-grep("2.fastq.gz",files, value = T)

# Define sample name (removing the file ending)
name<-fwd %>% str_sub(4,-12)

# Function to extract the MNase units
extract_units <- function(name) {
  split_string <- strsplit(name, "_")[[1]]
  # Assuming the units are always in the 3rd position, followed by 'U'
  units <- substr(split_string[3], 1, nchar(split_string[3])-1)
  return(units)
}

# Apply the function to each sample in the vector of strings
MNase_U <- sapply(name, extract_units)

#print(MNase_U)

# Create a dataframe with file paths and file names
df<-data.frame(Sample_Name=name,
           path_fwdReads=paste0(pathSTD,"/data/raw_data/",fwd),
           path_revReads=paste0(pathSTD,"/data/raw_data/",rev),
           MNase_U)

# Write a .csv file with fastq file names and paths
write.csv(df, row.names = F,
          file="data/samples.csv", quote = F)
