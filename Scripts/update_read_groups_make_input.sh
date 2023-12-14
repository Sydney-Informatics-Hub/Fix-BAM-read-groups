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

project=qe80
lstorage="scratch/qe80+gdata/qe80+scratch/er01"


# Config
config=/scratch/qe80/reheader_bams.config

# Output directory for new BAMs
outdir=/scratch/er01/PIPE-4346-BW-WGS/qe80_updated_readGroup_bams_360s 

# Sequencing platform
# Note that if not illumina, you will need to check the flowcell and lane variables 
# are correctly filled depending on your read id format 
platform=illumina 


# END UPDATE
#############################


# Update project and storage: 
sed -i "s|^#PBS -P.*|#PBS -P ${project}|g" ./Scripts/*pbs
sed -i "s|^#PBS -l[ ]*storage.*|#PBS -l storage=${lstorage}|g" ./Scripts/*pbs


# Set up output dirs
inputs=./Inputs
logdir=./PBS_logs

mkdir -p $outdir $logdir $inputs


# Write the inputs file
inputs=./Inputs/update_read_groups.inputs
rm -rf $inputs
while read line
do 
	sample=$(echo $line | awk '{print $1}')
	bam=$(echo $line | awk '{print $2}')
	centre=$(echo $line | awk '{print $3}')
	lib=$(echo $line | awk '{print $4}')
	
	size=$(ls -l ${bam} | awk '{print $5}')
	
	if ! [[ $lib ]]
	then
		lib=1
	fi
	
	printf "${sample},${bam},${outdir},${lib},${centre},${platform},${size}\n" >> $inputs
	
done < ${config}


# Add scale to assist efficient resourcing 
min=$(awk -F , '{print $7}' ${inputs} | sort -n | head -1)

with_scale=${inputs}-scaled
rm -rf ${with_scale}

while read line
do 
	size=$(echo $line | awk -F , '{print $7}')
	scale=$( printf '%.2f\n' $(echo "${size}/${min}" | bc -l) )
		
	printf "${line},${scale}\n" >> ${with_scale}
done < ${inputs} 


# Size sort to improve parallel efficiency
sort -t ',' -rnk 8 ${with_scale} > ${inputs}
rm ${with_scale}

echo `wc -l < ${inputs}` tasks written to ${inputs} 


