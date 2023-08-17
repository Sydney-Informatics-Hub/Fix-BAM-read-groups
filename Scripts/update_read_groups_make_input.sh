#!/bin/bash 

#########################################################
#
# Platform: NCI Gadi HPC
# Description:  Create input text file for parallel BAM read group update
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


#############################
# UPDATE VARIABLES: 

# Library ID
# if no specific library ID for this sample extraction, leave as 1
lib=1 

# Sequencing centre
centre=Ramaciotti 

# Sequencing platform 
platform=illumina 

# Location of input BAMs 
indir=/g/data/qe80/bams  

# Check assumption: that all BAMs in this dir are to be used
# Update if required
bams=($(ls ${indir}/*bam))
bais=($(ls ${indir}/*bai))

# Check assumption: that the sample ID is separated from BAM suffix by '.'
# Update of required 
samples=($(ls ${indir}/*bam | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1 )) 

# Output directory for new BAMs
outdir=/scratch/er01/PIPE-4346-BW-WGS/qe80_updated_readGroup_bams 

# Output directory for job logs:
logdir=./PBS_logs 

# END UPDATE
#############################

inputs=./Inputs/

mkdir -p $outdir $logdir $inputs

inputs=./Inputs/update_read_groups.inputs
rm -rf $inputs

for (( i = 0; i < ${#bams[@]}; i++ ))
do
	bam=${bams[$i]}
	sample=${samples[$i]}
	bai=${bais[$i]}
	
	printf "${sample},${bam},${bai},${outdir},${lib},${centre},${platform}\n" >> $inputs 
	
done 

echo `wc -l < $inputs` tasks written to $inputs 


