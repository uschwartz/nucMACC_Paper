#!/usr/bin/env bash

AnalysisDir="/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ"
bw=$AnalysisDir"/TSSplots/"
AnnoPath="/Users/admin/Annotation/GRCh38/"

cd $AnalysisDir

# nuccMACC scores
bwOccshH2AZ=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/shH2AZ_nucMACC_out/RUN/02_SUB-NUCS_PROFILE/
bwOccWC=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/WC_nucMACC_out/RUN/02_SUB-NUCS_PROFILE/

computeMatrix reference-point \
  -S $bwOccshH2AZ"/"*".bw" \
    -R $AnnoPath'/nextflow/protein_coding.gtf' \
    --referencePoint TSS  \
    --skipZeros \
    --missingDataAsZero \
    -o "TSSplots/computeMatrix2plot_TSS.txt.gz" \
    -b 900 -a 900 --smartLabels -p 12

plotProfile -m "TSSplots/computeMatrix2plot_TSS.txt.gz" \
               -out "TSSplots/TSS_sub-nucs_shH2AZ.pdf" \
               --perGroup --outFileNameData "TSSplots/TSS_sub-nucs_shH2AZ.tsv"

# normal

computeMatrix reference-point \
 -S $bwOccWC"/"*".bw" \
   -R $AnnoPath'/nextflow/protein_coding.gtf' \
   --referencePoint TSS  \
   --skipZeros \
   --missingDataAsZero \
   -o "TSSplots/computeMatrix2plot_TSS_WC.txt.gz" \
   -b 900 -a 900 --smartLabels -p 12

plotProfile -m "TSSplots/computeMatrix2plot_TSS_WC.txt.gz" \
              -out "TSSplots/TSS_sub-nucs_WC.pdf" \
              --perGroup --outFileNameData "TSSplots/TSS_sub-nucs_WC.tsv"

# mono Nucs 
bwOccshH2AZ=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/shH2AZ_nucMACC_out/RUN/01_MONO-NUCS_PROFILE/
bwOccWC=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/WC_nucMACC_out/RUN/01_MONO-NUCS_PROFILE/

computeMatrix reference-point \
-S $bwOccshH2AZ"/A"*".bw" \
  -R $AnnoPath'/nextflow/protein_coding.gtf' \
  --referencePoint TSS  \
  --skipZeros \
  --missingDataAsZero \
  -o "TSSplots/computeMatrix2plot_TSS_mono.txt.gz" \
  -b 900 -a 900 --smartLabels -p 12


plotProfile -m "TSSplots/computeMatrix2plot_TSS_mono.txt.gz" \
             -out "TSSplots/TSS_mono-nucs_shH2AZ.pdf" \
             --perGroup  --outFileNameData "TSSplots/TSS_mono-nucs_shH2AZ.tsv"

#### normal

computeMatrix reference-point \
-S $bwOccWC"/A"*".bw" \
 -R $AnnoPath'/nextflow/protein_coding.gtf' \
 --referencePoint TSS  \
 --skipZeros \
 --missingDataAsZero \
 -o "TSSplots/computeMatrix2plot_TSS_WC_mono.txt.gz" \
 -b 900 -a 900 --smartLabels -p 12

plotProfile -m "TSSplots/computeMatrix2plot_TSS_WC_mono.txt.gz" \
            -out "TSSplots/TSS_mono-nucs_WC.pdf" \
            --perGroup --outFileNameData "TSSplots/TSS_mono-nucs_WC.tsv"
