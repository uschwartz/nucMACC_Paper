#!/usr/bin/env bash

AnalysisDir=/Users/admin/Analysis/R001_nucMacc/
cd $AnalysisDir

nextflow run script/manuscript_figures/Fig6/01_nf_initiation \
-w manuscript_figures/work_init/Fig6/ -resume
