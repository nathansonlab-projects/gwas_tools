#!/bin/bash

# pluta 4/3/2017

# script to split plink file by chromosome, convert into vcf format, and check some parameters
# against HRC data (essentially pre phasing)

# phasing is handled by the HRC

# get user defined input
NPARAM=$#

INFILE=$1
CHR=$2
OUTDIR=$3

if [ $NPARAM -lt 3 ]
then
	echo ""
	echo "USAGE :: "
	echo "./preImputeCheck.sh INFILE CHR OUTDIR"
	echo "this program may use a lot of memory, set the -M flag in bsub if needed."
	echo "INFILE = root of the input PLINK file"
	echo "CHR = chromosome"
	echo "OUTDIR = the output directory"
	echo "currently only supports grch37"
	echo ""
	exit
fi



# ---------------- constants -------------------------------------- #
# user-defined
# directory where HRC-1000G-check-bim.pl exists
SCRIPTDIR=/project/knathanslab/TECAC/scripts

# location of HRC.r1-1.GRCh37.wgs.mac5.sites.tab reference file
REF=/project/knathanslab/REF/HHRC/grch37

# plink executable
PLLINKBIN=~/plink1.09

# vcf-sort executable
VCFSORTBIN=~/vcftools_0.1.13/bin/vcf-sort
# ----------------------------------------------------------------- #

BED=${INFILE}.bed
BIM=${INFILE}.bim
FAM=${INFILE}.fam

# do the necessary plink files exist?


if [ ! -d $OUTDIR ]
then
	mkdir $OUTDIR
fi

for i in $BED $BIM $FAM
do
	if [ ! -f ${i} ]
	then
		echo "$i not found!"
		exit 1	
	fi
done

echo "PLINK files found"
echo "extracting chr${CHR}..."

# extract chromosome
CMD="${PLINKBIN} --bfile $INFILE 
                   --chr $CHR 
                   --recode 
                   --make-bed 
                   --out ${OUTDIR}/chr${CHR}"
eval $CMD

if [ $? -ne 0 ]
then
	echo "$CMD failed! Aborting."
	exit 1
else
	echo "Done!"
fi


cd $OUTDIR

# get MAF
${PLINKBIN} --bfile chr${CHR} --freq --out maf --noweb

echo "align to HRC ..."

# compare phased data to HRC and apply corrections
# this part of the script is time and memory intensive
# the smallest chromosome, chr22, used 22244 MB of memory
# chr1 is the largest, and is about 6 times larger than chr22
CMD="perl ${SCRIPTDIR}/HRC-1000G-check-bim.pl -b chr${CHR}.bim 
                                              -f maf.frq 
                                              -r ${REF}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h"
eval $CMD

if [ $? -ne 0 ]
then
	echo "alignment failed!"
	echo "Aborting."
	exit 1
else
	echo "Done!"
fi

if [ ! -f Exclude-chr${CHR}-HRC.txt ]
then
	echo "HRC-1000G-check-bim.pl output missing."
	echo "Aborting"
	exit 1
fi


echo "Correcting PLINK files..."

# these files are all output from the above script
${PLINKBIN} --bfile chr${CHR} --exclude Exclude-chr${CHR}-HRC.txt --make-bed --out TEMP1
${PLINKBIN} --bfile TEMP1 --update-map Chromosome-chr${CHR}-HRC.txt --update-chr --make-bed --out TEMP2
${PLINKBIN} --bfile TEMP2 --update-map Position-chr${CHR}-HRC.txt --make-bed --out TEMP3
${PLINKBIN} --bfile TEMP3 --flip Strand-Flip-chr${CHR}-HRC.txt --make-bed --out TEMP4
${PLINKBIN} --bfile TEMP4 --reference-allele Force-Allele1-chr${CHR}-HRC.txt --make-bed --out chr${CHR}-fix 

echo "done!"
rm TEMP*
rm temp*

echo "converting to VCF...."
CMD="${PLINKBIN} --bfile chr${CHR}-fix --recode vcf-iid --out chr${CHR}"
eval $CMD

if [ $? -ne 0 ]
then 
	echo "$CMD failed!"
	echo "Aborting"
	exit 1
fi

$VCFSORTBIN chr${CHR}.vcf | bgzip -c > chr${CHR}.vcf.gz


