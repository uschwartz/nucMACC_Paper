#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper/data/R001_nucMacc/manuscript_figures/Fig4/motifs/
cd $AnalysisDir

# find motifs (using HOMER tool) 100bp upstream and 100 bp downstream of TSS

# for unstable TSS
findMotifs.pl ../data/TSSgroups/unstable_TSS.txt fly unstableNuc \
-start -100 -end 100 -p 20 -bg ../data/TSSgroups/all_expr_TSS.txt

# for stable TSS
findMotifs.pl ../data/TSSgroups/stable_TSS.txt fly NDR \
-start -100 -end 100 -p 20 -bg ../data/TSSgroups/all_expr_TSS.txt
