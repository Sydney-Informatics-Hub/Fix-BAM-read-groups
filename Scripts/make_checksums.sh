#!/bin/bash 

#########################################################
#
# Platform: NCI Gadi HPC
# Description:  Make checksums for new BAMs 
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
old_bam=`echo $1 | cut -d ',' -f 2`
bam_dir=`echo $1 | cut -d ',' -f 3`

prefix=$(basename $old_bam | sed 's/\.bam//')
new_bam=${bam_dir}/${prefix}.rh.bam

md5sum ${new_bam} > ${new_bam}.md5
