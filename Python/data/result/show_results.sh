#! /bin/bash
name="_res_5.txt"
for dataset in USAir NS PB Yeast Celegans Power Router Ecoli
do
cat $dataset$name | awk '{if (NR==1) print;}'
done
for dataset in USAir NS PB Yeast Celegans Power Router Ecoli
do
cat $dataset$name | awk '{if (NR==2) print;}'
done
for dataset in USAir NS PB Yeast Celegans Power Router Ecoli
do
cat $dataset$name | awk '{if (NR==3) print;}'
done
for dataset in USAir NS PB Yeast Celegans Power Router Ecoli
do
cat $dataset$name | awk '{if (NR==4) print;}'
done
