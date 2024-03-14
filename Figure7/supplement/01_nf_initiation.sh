#!/usr/bin/env bash

projectDir=/home/$USER/nucMACC_Paper/data/R001_nucMacc

cd $projectDir

nextflow run $projectDir/script/manuscript_figures/Fig6/01_nf_initiation \
-w $projectDir/manuscript_figures/work_init/Fig6/ 
