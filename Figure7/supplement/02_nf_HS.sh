AnalysisDir=/Users/admin/Analysis/R001_nucMacc/
cd $AnalysisDir

nextflow run script/manuscript_figures/Fig6/02_nf_HS \
-w manuscript_figures/work_init/Fig6/ -resume
