#!/bin/bash

set -e

snpfile=${1}
improot=${2}
keepfile=${3}
output=${4}

snplist1="${output}.snplist1"
snplist2="${output}.snplist2"
snplist3="${output}.snplist3"

extractSnps="/panfs/panasas01/sscm/gh13047/repo/randomScripts/extractSnps/extractSnps.sh"

# Extract SNPs
awk '{ print $1 }' ${snpfile} | sed 1d > ${snplist1}
${extractSnps} ${snplist1} ${keepfile} ${improot} ${output}_snps

# Clump SNPs
echo "SNP P" > ${snplist2}
awk '{ print $1, $4 }' ${snpfile} | sed 1d >> ${snplist2}
plink1.90 --bfile ${output}_snps --clump ${snplist2} --make-bed --out ${output}_clump

# Construct profile score
awk '{ print $1, $6, $5 }' ${snpfile} | sed 1d > ${output}_temp
awk -v f1="${output}_temp" -v f2="${output}_clump.clumped" ' FILENAME==f1 {arr[$1]=$0; next} 
FILENAME==f2 {print arr[$3]} ' ${output}_temp ${output}_clump.clumped > ${snplist3}

plink --noweb --bfile ${output}_snps --score ${snplist3} --out ${output}
