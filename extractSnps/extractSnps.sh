#!/bin/bash

# e.g. ./extractSnps.sh ../data/crapsnps.txt /ibscratch/wrayvisscher/imputation/twge/data/imputed/chrCHR/twge_1kg_p1v3_CHR ../data/twge_crap

snplistfile=${1}
plinkrt=${2}
outfile=${3}

touch ${outfile}_mergelist.txt
rm ${outfile}_mergelist.txt
touch ${outfile}_mergelist.txt

for i in {1..22}
do
        filename=$(sed -e "s/CHR/$i/g" <<< ${plinkrt})
        plink --noweb --bfile ${filename} --extract ${snplistfile} --make-bed --out ${outfile}_${i}
        echo "${outfile}_${i}.bed ${outfile}_${i}.bim ${outfile}_${i}.fam" >> ${outfile}_mergelist.txt
done
cat ${outfile}_mergelist.txt
sed -i 1d ${outfile}_mergelist.txt
cat ${outfile}_mergelist.txt
plink --noweb --bfile ${outfile}_1 --merge-list ${outfile}_mergelist.txt --make-bed --out ${outfile}

rm ${outfile}_*