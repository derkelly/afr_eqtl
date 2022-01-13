#!/bin/bash

# master directories
export PROJDIR="/project/tishkofflab/rna_redo/afr_wbrna"
export READDIR="/project/tishkofflab/afiimani/STUDY/reads"

# code directories
export PIPEDIR=${PROJDIR}/pipeline
export ERRDIR=${PIPEDIR}/err
export LOGDIR=${PIPEDIR}/log

# data directories
export DATADIR=${PROJDIR}/data
export TEMPDIR=${DATADIR}/tmp
export SAMPDIR=${DATADIR}/samples
export MISCDIR=${DATADIR}/misc
export GENDIR=${DATADIR}/geno
export GENODIR=${DATADIR}/genomes
export RASQDIR=${DATADIR}/rasqual
export METADIR=${DATADIR}/metal
export LEAFDIR=${DATADIR}/leafcutter
export FASTDIR=${DATADIR}/fastqtl
export FSTDIR=${DATADIR}/fst
export EQTLDIR=${DATADIR}/eqtl
export CHUNKDIR=${DATADIR}/chunks

# list of samples
readarray -t SAMPLES<${MISCDIR}/samples.txt
readarray -t SAMPS<${MISCDIR}/samps.txt
readarray -t POPS<${MISCDIR}/all_pops.txt

export $SAMPLES
export $SAMPS
export $POPS
