#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper/data/R001_nucMacc
cd $AnalysisDir

nextflow run $AnalysisDir/script/manuscript_figures/Fig6/02_nf_HS \
-w $AnalysisDir/manuscript_figures/work_init/Fig6/
