#!/bin/bash

# e.g. ./extractSnps.sh snplist.txt 

snplistfile=${1}
keepfile=${2}
plinkrt=${3}
outfile=${4}

touch ${outfile}_mergelist.txt
rm ${outfile}_mergelist.txt
touch ${outfile}_mergelist.txt

for i in {1..22}
do
        filename=$(sed -e "s/CHR/$i/g" <<< ${plinkrt})
        plink1.90 --bfile ${filename} --keep ${keepfile} --extract ${snplistfile} --make-bed --out ${outfile}_${i}
        if [ -f "${outfile}_${i}.bim" ]; then
	        echo "${outfile}_${i}.bed ${outfile}_${i}.bim ${outfile}_${i}.fam" >> ${outfile}_mergelist.txt
	   	fi
done
cat ${outfile}_mergelist.txt
sed -i 1d ${outfile}_mergelist.txt
cat ${outfile}_mergelist.txt
plink1.90 --bfile ${outfile}_1 --merge-list ${outfile}_mergelist.txt --make-bed --out ${outfile}

rm ${outfile}_*