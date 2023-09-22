#!/bin/bash


inputs=Inputs/update_read_groups.inputs
err=0

# Check PBS logs: validate
val_log=PBS_logs/validate_bams.o
val_e=PBS_logs/validate_bams.e
all_exit=$(grep "Exit" ${val_log} | awk '{print $3}')
if [[ $all_exit -eq 0 ]]
then
        printf "Parent validate job exit status OK\n"
else 
        printf "ERROR: Parent validate job exit status ${all_exit}\n\n"
fi

tasks=$(wc -l < $inputs)
success_validate=$(grep -zoP 'Tool returned:\n0' $val_log | wc -l)
if [[ $success_validate -lt $tasks ]]
then 
	printf "ERROR: $tasks total tasks yet only $success_validate successful in $val_log\n"
else
	printf "All validate tasks returned 0 from tool\n"	
fi

success_validate=$(grep 'exited with status 0' $val_e| wc -l)
if [[ $success_validate -lt $tasks ]]
then 
	printf "ERROR: $tasks total validate tasks yet only $success_validate exited with status 0\n"
else
	printf "All validate tasks exited with status 0\n"	
fi



# Check PBS logs: flagstats
flag_log=PBS_logs/compare_flagstats.o
flag_e=PBS_logs/compare_flagstats.e
all_exit=$(grep "Exit" ${flag_log} | awk '{print $3}')
if [[ $all_exit -eq 0 ]]
then
        printf "Parent flagstats job exit status OK\n"
else 
        printf "ERROR: Parent flagstats job exit status ${all_exit}\n\n"
fi

success_flagstats=$(grep 'exited with status 0' $flag_e | wc -l)
if [[ $success_flagstats -lt $tasks ]]
then 
	printf "ERROR: $tasks total flagstats tasks yet only $success_flagstats exited with status 0\n"
else
	printf "All flagstats tasks exited with status 0\n"	
fi

# Check per task outputs: 
while read line
do
	sample=`echo $line | cut -d ',' -f 1`
	old_bam=`echo $line | cut -d ',' -f 2`
	bam_dir=`echo $line | cut -d ',' -f 3`

	prefix=$(basename $old_bam | sed 's/\.bam//')
	
	validate=${bam_dir}/${prefix}.validate.log
	diff=${bam_dir}/${prefix}.diffStats
	
	# Check that there is no diff file from flagstats (empty diffs were deleted per task)
        if [[ -s $diff ]]
        then
		printf "ERROR: $diff exists and is non-empty\n" 
		((err+=1))
	fi
	
	# Check that the picard validation file reported "No errors found"
	if ! grep -q "No errors found" ${validate}
	then 
		printf "ERROR: ${validate} is missing 'No errors found' status\n"
		((err+=1))
	fi
done < ${inputs}

if [[ $err -eq 0 ]]
then
	printf "Picard ValidateSamFiles and samtools flagstats detected no issues with the new BAM files.\n"
	printf "Please proceed to make checksums.\n"
else
	printf "$err errors detected - please check logs and outputs commencing from step 1.\n"
fi



