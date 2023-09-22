#!/bin/bash 

#########################################################
#
# Platform: NCI Gadi HPC
# Description:  Simple check for compelteness with flagstats
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


old_stats=${bam_dir}/${prefix}.OLD.flagstats
new_stats=${bam_dir}/${prefix}.rh.flagstats
diff_stats=${bam_dir}/${prefix}.diffStats

samtools flagstat -@ $PBS_NCPUS $old_bam > $old_stats
samtools flagstat -@ $PBS_NCPUS $new_bam > $new_stats

sdiff -s $old_stats $new_stats > ${diff_stats}

if ! [ -s ${diff_stats} ]
then 
	# Diff file is empty
	rm -rf ${diff_stats} ${old_stats}
fi
