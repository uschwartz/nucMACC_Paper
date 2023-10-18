AnalysisDir=~/Analysis/R001_nucMacc

cd $AnalysisDir

outDir=manuscript_figures/Fig4/data/GROseq/raw

nextflow run ~/00_scripts/nextflow/RNAseq/main.nf  \
        --fastqPath $AnalysisDir"/"$outDir \
	--outPath $AnalysisDir"/manuscript_figures/Fig4/GROseq/" \
	--STARidxPath ~/Annotation/Drosophila_melanogaster_UCSC_dm3/STARidx \
	--gtfPath ~/Annotation/Drosophila_melanogaster_UCSC_dm3/nextflow \
	--gtfFile protein_coding.gtf --strandness unstranded
