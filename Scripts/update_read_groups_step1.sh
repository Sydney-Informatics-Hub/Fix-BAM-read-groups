#!/bin/bash 

#########################################################
#
# Platform: NCI Gadi HPC
# Description:  BAM read group update step 1
# see github.com/Sydney-Informatics-Hub/Fix-BAM-read-groups
#
# Author/s: Cali Willet
# cali.willet@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the scientific and technical assistance
# <or e.g. bioinformatics assistance of <PERSON>> of Sydney Informatics
# Hub and resources and services from the National Computational
# Infrastructure (NCI), which is supported by the Australian Government
# with access facilitated by the University of Sydney.
#
#########################################################

sample=`echo $1 | cut -d ',' -f 1`
bam=`echo $1 | cut -d ',' -f 2`
outdir=`echo $1 | cut -d ',' -f 3`
lib=`echo $1 | cut -d ',' -f 4`
centre=`echo $1 | cut -d ',' -f 5`
platform=`echo $1 | cut -d ',' -f 6`

prefix=$(basename $bam | sed 's/\.bam//')

# Get header: 
samtools view -@ $NCPUS -Ho ${outdir}/${prefix}.header ${bam}

# Get SAM: 
# Inlcude the header so 'reheader' tool can be used at step 3
samtools view -@ $NCPUS -ho ${outdir}/${prefix}.sam ${bam}
