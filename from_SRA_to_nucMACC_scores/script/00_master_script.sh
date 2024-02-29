export MY_PWD=$(pwd)

./script/01_fetch_annotation.sh $MY_PWD
./script/02_fetch_Data_zsh.sh $MY_PWD
Rscript ./script/03_createCSV_1.R
./script/04_runNextflow.sh $MY_PWD
Rscript ./script/05_createCSV_2.R
./script/06_runNextflow_nucMACC.sh $MY_PWD
