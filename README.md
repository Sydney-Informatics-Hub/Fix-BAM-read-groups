# Fix-BAM-read-groups
Change the read group metadata within a BAM file. Operates on the header as well as the individual SAM output lines.

 
# Update variables below, then run this script which will submit
# update_read_groups.pbs over each sample in the specified BAM input directory

## When to use

If you simply need to change the header, you can use [samtools reheader](http://www.htslib.org/doc/samtools-reheader.html). 

If you need to update read groups and you are CERTAIN that all of the reads in the BAM come from the same read group, you can use [Picard AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/13832752709403-AddOrReplaceReadGroups-Picard-).

If you need to update the read groups and have no idea how many read groups are contained within the BAM and/or know the read groups but need to make a change, this set of scripts will be useful.


If you are planning on using GATK with your BAMs, the following minimum read group fields are required: ID, PU*, SM, PL, LB. PU is not required by GATK but takes precedence over ID for base recalibration if it is present. ID (or PU) and LB are important for BQSR and MarkDuplicates, see [this blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671) for discussion.   


## Overview

Takes a BAM file and updates the @RG header lines and read group IDs within a BAM file, by first extracting the headers, then converting BAM to SAM, reading the SAM file line by line, and capturing the flowcell and lane from the unique read ID. The flowcell-lane info is saved throughput, as a new SAM file is written using the updated read group information. 

After the new SAM has been written, the collection of read groups (can be from one to many) are written as separate @RG headers within the headers file, over-writing any existing @RG headers. The new headers and new SAM are concatenated, and converted to a new updated BAM file. 

## Assumptions

- The BAM is indexed
- The data is Illumina, and has flowcell and lane as fields 3 and 4 respectively within the ':' delimited read ID
- Running on NCI Gadi HPC 
- You have sufficient disk space to temporarily hold 2 SAM files and 2 new BAM file per input BAM file. Depending on compression level of the input BAM, you may need ~15 X the disk space of the input BAM (ie 1.5 TB for a 100 GB BAM). If you are limited with processing space, you will need to adjust the workflow described below to process fewer samples at a time.  

## Workflow

### Make parallel inputs file

Open `./Scripts/update_read_groups_make_input.sh` and adjust variables specific to your dataset:

- lib: Library ID. If no specific library ID for this sample extraction, leave as 1
- centre: Sequencing centre
- platform: Sequencing platform (should be illumina)  
- indir: Location of input BAMs. Will operate on all BAMs in this directory, and assumes that the sample ID of these BAMs is separated from BAM suffix by '.'
- outdir: Output directory for new BAMs
- logdir: Output directory for PBS job logs

Save the script, then run with:
```
bash ./Scripts/update_read_groups_make_input.sh
```

Output will be `./Inputs/update_read_groups.inputs` which contains the required inforamtion for each parallel task (one per input BAM) to be executed at the next step. 

### Step 1 - convert BAM to SAM

This step simply converts the BAM to a SAM file and extracts the current headers to a separate file. This step has been separated from step 2 as it multithreads where step 2 does not.
 

Open `Scripts/update_read_groups_run_parallel_step1.pbs` and edit for your project, lstorage, and number of samples. 7 Broadwell CPU with 9 GB RAM per CPU is a good trade-off between walltime, CPU efficiency and SU. Other resource configurations may be faster, but cost more and are less efficient. 

For example if you have 20 samples, applying 7 CPU per sample on the Broadwell normal queues, request a total of 7 x 20 CPUS = 1240 CPUs and 140 x 9 GB RAM = 1260 GB RAM. 

Walltime for one 89 GB BAM which uncompressed to 720 GB SAM was 24 minutes on the above resource settings. 

Submit step 1 with:

```
qsub Scripts/update_read_groups_run_parallel_step1.pbs
```

Output will be SAM and header files for each sample in inputs file, in the outdir specified in `update_read_groups_make_input.sh`. 


### Step 2 - update read groups in the SAM and RG headers in the headers file

Open `Scripts/update_read_groups_run_parallel_step2.pbs` and edit for your project, lstorage, and number of samples. This step does not use of multiple thread. If the number of samples is greater than the number of CPU on the node (48 for Cascade Lake or 28 for Broadwell) request whole nodes. 

This step is the slowest (multiple hours) so if you have a great discrepancy in the size of your input BAM files, you may have idle CPU in the job for many hours. You may choose to increase efficiency and reduce cost by for:

- Size sorting the inputs according to BAM size (largest to smallest) and submitting with less CPU than required to run all tasks at once
- Splitting the inputs into two or more similarly sized input lists and submit them as separate jobs (change .o and .e output log file names to avoid over-write) 
- Running a `for` loop over the inputs file, parsing the qsub directives on the command line, and allowing a 2 second sleep between each iteration of the loop

To submit the job in the default way (not using any of the above 3 efficiency strategies, run:
```
qsub Scripts/update_read_groups_run_parallel_step2.pbs
```

Output will be `<outdir>/<sample>_newRG.SAM` per sample. The headers files created at step 1 will be edited in-place. 


### Step 3 - convert updated SAM to BAM

Like step 1, this step multithreads so is separated from step 2. The same CPU and mem settings for step 1 (7 Broadwell CPU with 9 GB RAM per CPU) provides the optimal tradeoff between walltime, service unis and CPU efficiency. 

Open `Scripts/update_read_groups_run_parallel_step3.pbs` and edit for your project, lstorage, and number of samples.

Submit with:

```
qsub Scripts/update_read_groups_run_parallel_step3.pbs
```
Output will be `<outdir>/<sample>.final.bam` and `<outdir>/<sample>.final.bai`. The previous BAI file is simply copied from the old BAI in order to update filename and timestamp. The SAM and header files created during steps 1 and 2 are removed. 



