#/bin/bash
#
# This is extract_geno.sh
#
# This extracts genotypes for individuals belonging to each population.

INFILE=$1
DATADIR=$2
POP=$3

POPDIR=${DATADIR}/geno/${POP}

mkdir -p $POPDIR

SAMPS=${DATADIR}/misc/pops/in5M/${POP}.txt

OUTPREF=${POPDIR}/${POP}.5M.imputed.dose.biallelic

~/bin/vcftools --gzvcf $INFILE \
    --keep $SAMPS \
    --mac 1 \
    --recode \
    --out $OUTPREF

bgzip -f ${OUTPREF}.recode.vcf
tabix -p vcf -f ${OUTPREF}.recode.vcf.gz
