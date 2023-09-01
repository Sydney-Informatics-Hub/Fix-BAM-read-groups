#!/bin/bash 

#########################################################
#
# Platform: NCI Gadi HPC
# Description:  Convert updated SAM to BAM
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
bai=`echo $1 | cut -d ',' -f 3`
outdir=`echo $1 | cut -d ',' -f 4`

# Convert updated SAM to BAM: 
samtools view -@ $NCPUS -bo ${outdir}/${sample}_newRG.bam ${outdir}/${sample}_newRG.sam

# Add the header
samtools reheader -P ${outdir}/${sample}.header ${outdir}/${sample}_newRG.bam > ${outdir}/${sample}.final.bam

# Re-index: 
samtools index -@ ${NCPUS} ${outdir}/${sample}.final.bam

# remove temp files:
rm -rf ${outdir}/${sample}.header ${outdir}/${sample}*.sam ${outdir}/${sample}_newRG.bam
