AnalysisDir=~/Analysis/R001_nucMacc/manuscript_figures/Fig4/motifs/

cd $AnalysisDir

findMotifs.pl ../data/TSSgroups/unstable_TSS.txt fly unstableNuc \
-start -100 -end 100 -p 20 -bg ../data/TSSgroups/all_expr_TSS.txt


findMotifs.pl ../data/TSSgroups/stable_TSS.txt fly NDR \
-start -100 -end 100 -p 20 -bg ../data/TSSgroups/all_expr_TSS.txt
