#!/bin/bash
#
# This is map_STAR_firstpass.sh
#
# This script calls 'map_chunk_STAR_firstpass.sh' in order to map
# chunks of samples at once.

PROJDIR="/project/tishkofflab/rna_redo/afr_wbrna"
READDIR="/project/tishkofflab/afiimani/STUDY/reads"
PIPEDIR=${PROJDIR}/pipeline
DATADIR=${PROJDIR}/data
OUTDIR=${DATADIR}/samples

CHUNKS=$(cd ${DATADIR}/chunks/ ; ls chunk*)

# call map_chunk_STAR_firstpass.sh for each chunk
for CHUNK in ${CHUNKS[@]};
do

    ERR=${PIPEDIR}/err
    LOG=${PIPEDIR}/log

    bsub -n 17 \
	-M 35840 \
	-e ${ERR}/${CHUNK}_firstpass.err \
    	-o ${LOG}/${CHUNK}_firstpass.log \
    	bash map_chunk_STAR_firstpass.sh \
    	${DATADIR}/chunks/${CHUNK} \
    	$READDIR \
    	$OUTDIR

done
