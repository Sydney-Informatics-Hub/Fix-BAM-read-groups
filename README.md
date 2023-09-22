# Fix-BAM-read-groups
Change the read group metadata within a BAM file. Operates on the header as well as the individual SAM output lines.


## When to use

If you simply need to change the header, you can use [samtools reheader](http://www.htslib.org/doc/samtools-reheader.html). 

If you need to update read groups and you are CERTAIN that all of the reads in the BAM come from the same read group, you can use [Picard AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/13832752709403-AddOrReplaceReadGroups-Picard-).

If you need to update the read groups and have no idea how many read groups are contained within the BAM and/or know the read groups but need to make a change, this set of scripts will be useful.


If you are planning on using GATK with your BAMs, the following minimum read group fields are required: 

- ID
- PU*
- SM
- PL
- LB

PU is not required by GATK but takes precedence over ID for base recalibration if it is present. ID (or PU) and LB are important for BQSR and MarkDuplicates, see [this blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671) for discussion.

This workflow adds these read group fields to you BAM, as well as CN (sequencing centre).    


## Overview

Takes a BAM file and updates the @RG header lines and read group IDs within a BAM file, by first extracting the headers, then converting BAM to SAM, reading the SAM file line by line, and capturing the flowcell and lane from the unique read ID in order to update both the RG headers and the read group IDs. 

The new RG headers (based on the user-specifed values and the flowcell and lane derived from the read IDs) over-write any existing @RG headers. The updated SAM file is then converted back to BAM format and new BAI index created. 

The workflow is broken up into 3 steps to separate the multi-threading steps (1 and 3) from the single-threaded step 2.

Steps 1 and 3 are simply SAM/BAM conversions, while step 2 performs the editing. 

The read group IDs (attached to every read in the BAM file) are updated to:

`<flowcell>.<lane>.<sampleID>_<lib>`

The read group headers (one to many, reflecting the total number of flowcell-lanes found in the BAM) are updated to:

`@RG\tID:<flowcell>.<lane>.<sampleID>_<lib>\tPL:<platform>\tPU:<flowcell>.<lane>\tSM:<sampleID>\tLB:<sampleID>\_<lib>\tCN:<centre>`

This verbose format is adopted to adhere to universal standards and ensure downstream compatibility with GATK and other tools that make use of read group metadata. 


## Assumptions

- The data is Illumina, and has flowcell and lane as fields 3 and 4 respectively within the ':' delimited read ID
- Running on NCI Gadi HPC 
- You have sufficient disk space to temporarily hold 2 SAM files and 2 new BAM files per input BAM file. Depending on compression level of the input BAM, you may need ~15 X the disk space of the input BAM (ie 1.5 TB for a 100 GB BAM). If you are limited with processing space, you will need to adjust the workflow described below to process fewer samples at a time.  

## Workflow

### Make a sample config file

This is a plain TSV that contains 3 columns: sample ID (to be used for the `'SM'` tag), full path to the BAM file, and sequencing centre. An optional fourth column can be added containing the library ID. If no fourth column, the default library value of 1 will be used. 

Example format:
```
sample1	/scratch/myproj/sample1_lane1.bam	KCG
sample1	/scratch/myproj/sample1_lane2.bam	KCG
sample2	/scratch/myproj/sample2_lane3.bam	Ramaciotti	2a
sample2 /scratch/myproj/sample2_lane7.bam	Ramaciotti	2b 
```

### Make parallel inputs file

Open `./Scripts/update_read_groups_make_input.sh` and adjust variables specific to your dataset:

- project: NCI project code for accounting
- lstorage: list of storage paths required for the job, in the format required for the NCI lstorage PBS directive eg "scratch/myproj1+gdata/myproj2"
- config: Config file as described above, including path
- outdir: Output directory for new BAMs (will be created if it does not exist)


Save the script, then run with:
```
bash ./Scripts/update_read_groups_make_input.sh
```

Output will be `./Inputs/update_read_groups.inputs` which contains the required information for each parallel task (one per input BAM) to be executed at the next step.

The comma-delimited inputs file will be sorted from largest BAM to smallest. It will list the size in bytes of the BAM file (column 7) and the ratio of the BAM size to the smallest BAM (column 8). This is to aid in resourcing and compute efficiency for jobs where there is a big size difference between the BAMs to be reheadered. This is important as the walltime per BAM is dependent on BAM size.   


### Step 1: Convert BAM to SAM

This step simply converts the BAM to a SAM file and extracts the current headers to a separate file. This step has been separated from step 2 as it multithreads where step 2 does not.
 

Open `Scripts/update_read_groups_run_parallel_step1.pbs` and adjust resources for your number of samples. 

7 Broadwell CPU with 9 GB RAM per CPU is a good trade-off between walltime, CPU efficiency and SU. Other resource configurations may be faster, but cost more and are less efficient. 

For example if you have 20 samples, applying 7 CPU per sample on the Broadwell normal queues, request a total of 7 x 20 CPUS = 1240 CPUs and 140 x 9 GB RAM = 1260 GB RAM to run all tasks at the same time. Note that if you have very large and very small BAMs, idle CPU will occur when the small BAMs complete long before the large BAMs. Since the inputs are size-sorted, you can request less CPU than required to run all tasks in parallel,  to increase overall CPU efficiency and decrease SU usage. Column 8 of the inputs file is designed to help you estimate the apropriate resourcing here. 

Walltime for one 89 GB BAM which uncompressed to 720 GB SAM was 24 minutes on the above resource settings. 

Submit step 1 with:

```
qsub Scripts/update_read_groups_run_parallel_step1.pbs
```

Output will be `<outdir>/<bam_prefix>.sam` and <outdir>/<bam_prefix>.header` files for each sample in `./Inputs/update_read_groups.inputs`. 


### Step 2: Update read groups in the SAM and RG headers in the headers file

Open `Scripts/update_read_groups_run_parallel_step2.pbs` and adjust resources for your number of samples. 

This step does not multi-thread. If the number of samples is greater than the number of CPU on the node (48 for Cascade Lake or 28 for Broadwell) request whole nodes. 

This step is the slowest (3.5 hours for test sample with 89 GB BAM / 720 GB SAM) so if you have a great discrepancy in the size of your input BAM files, you may have idle CPU in the job for multiple hours. You may choose to increase efficiency and reduce cost by:

- Since the inputs are size sorted according to BAM size (largest to smallest), submitting with less CPU than required to run all tasks at once will enable small tasks to cycle through while larger tasks continue to run
- Splitting the inputs into two or more similarly sized input lists and submit them as separate jobs (change .o and .e output log file names in the PBS directives to avoid over-write) 
- Running a `for` loop over the inputs file, parsing the qsub directives on the command line, and adding a 2 second `sleep` between each iteration of the loop, so that different walltimes can be requested for different tasks

To submit the job in the default way, run:
```
qsub Scripts/update_read_groups_run_parallel_step2.pbs
```

Output will be `<outdir>/<bam_prefix>_newRG.sam` per sample. The headers files created at step 1 will be edited in-place. 


### Step 3: Convert updated SAM to BAM

Like step 1, this step multithreads so is separated from step 2. The same CPU and mem settings for step 1 (7 Broadwell CPU with 9 GB RAM per CPU) provides the optimal trade-off between walltime, service units and CPU efficiency. 

Open `Scripts/update_read_groups_run_parallel_step3.pbs` and adjust resources for your number of samples.

Submit with:

```
qsub Scripts/update_read_groups_run_parallel_step3.pbs
```
Output will be `<outdir>/<bam_prefix>.rh.bam` and `<outdir>/<bam_prefix>.rh.bai`. The SAM and header files created during steps 1 and 2 are removed. 


## Checking the output

### Manually inspect the headers and read group IDs

You should check that the read group headers and read group IDs look correct on at least one updated BAM. 

Check the read group headers:

```
module load samtools
samtools view -H <updated_BAM> | grep "@RG"
```

Example expected output (this will vary for your data depending on user-specified parameter values and flowcell-lanes detected in your BAM):
```
@RG     ID:H5WGHDSX2.1.STHD_F_Daenery_1 PL:illumina     PU:H5WGHDSX2.1  SM:STHD_F_Daenery       LB:STHD_F_Daenery_1     CN:Ramaciotti
@RG     ID:H5WGHDSX2.2.STHD_F_Daenery_1 PL:illumina     PU:H5WGHDSX2.2  SM:STHD_F_Daenery       LB:STHD_F_Daenery_1     CN:Ramaciotti
@RG     ID:H5WGHDSX2.3.STHD_F_Daenery_1 PL:illumina     PU:H5WGHDSX2.3  SM:STHD_F_Daenery       LB:STHD_F_Daenery_1     CN:Ramaciotti
@RG     ID:H5WGHDSX2.4.STHD_F_Daenery_1 PL:illumina     PU:H5WGHDSX2.4  SM:STHD_F_Daenery       LB:STHD_F_Daenery_1     CN:Ramaciotti
```


Check the read group IDs:

```
module load samtools
samtools view <updated_BAM> | head -3
```

Example expected output (read group ID is field  'RG:Z:<rgid>' - note the column number of this field can vary between datasets:
```
A00152:414:H5WGHDSX2:1:1213:24352:16861 163     MSTS01000001.1  1       2       60S87M3S        =       484     634     GTACATAGTATTAGCAATATTGTGTGATGATCAACTATGAATGACTCAGTTCCCAGCAAAACAATGATCCAAGACAAGTACAAAGAACTCACAAAGGAAAATGCTATCCACGTCCAGATAAAGAACTGATGGACTTTGTAGATCAAAGTA    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF:FFF:FFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFF,FFFFFFFFFF:,FF:FFF:FFFFFFFFFFFFFFFFFFFFFFF  MC:Z:94M1D56M   PG:Z:MarkDuplicates     RG:Z:H5WGHDSX2.1.STHD_F_Daenery_1       NM:i:0  MQ:i:60 AS:i:87 XS:i:112
A00152:414:H5WGHDSX2:2:1359:30219:3192  163     MSTS01000001.1  1       32      49S87M1I13M     =       598     746     TAGCAATATTGTGTGATGATCAACTATGAATGACTCAGTTCCCAGCAAAACAATGATCCAAGACAAGTACAAAGAACTCACAAAGGAAAATGCTATCCACGTCCAGATAAAGAACTGATGGACTTTGTAGATCAAAGTACATATTTTTTA    FF:FFFFFFFF:FFFFFF:FFFF,:FFFFFFFFFF:FFFF:FFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFF,FFFFFFFFFFFF,FF:FFFFFFFFFFFFFF:FFFFFFFF:,FFFFFF:FFFFFFFF  MC:Z:76M1I73M   PG:Z:MarkDuplicates     RG:Z:H5WGHDSX2.2.STHD_F_Daenery_1       NM:i:1  MQ:i:60 AS:i:93 XS:i:106
A00152:414:H5WGHDSX2:2:2171:16034:28479 99      MSTS01000001.1  1       60      27S84M1I39M     =       270     420     ACTATGAATGACTCAGTTCCCAGCAAAACAATGATCCAAGACAAGTACAAAGAACTCACAAAGGAAAATGCTATCCACGTCCAGATAAAGAACTGATGGACTTTGTAGATCAAAATACATATTTTTTAAACTTTATTCTTTTTTGTGTGTT   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF MC:Z:151M       PG:Z:MarkDuplicates     RG:Z:H5WGHDSX2.2.STHD_F_Daenery_1       NM:i:1  MQ:i:60 XQ:i:156 AS:i:116 XS:i:89
```

### Compare flagstat output

If something went awry with the read line printing, for example missing a read, samtools flagstats would emit different values. This repository includes a quick checker script that can be used to compare the outputs of `samtools flagstats`.

Adjust resources depending on your number of samples, allowing 7 CPUs and 9 GB per CPU per sample on the normal Broadwell queue, and submit with:
```
qsub Scripts/compare_flagstats_run_parallel.pbs
```

Output will be flagstats files in the same directory that the newly created BAM files are in. 

The script runs flagstats on the new and the original BAM files, and then uses Linux sdiff to compare them. If the flagstats files are the same, the 'old bam' flagstats file will be deleted. If there is a difference, the flagstats for the original BAM will not be deleted, and a `<outdir>/<bam_prefix>.diffStats` file will be present.  

Errors at this stage can be detected by the final check script, to be run after flagstats AND validate check steps are run. 

### Validate BAM files with Picard

This step can be carried out at the same time as the `compare_flagstats` job. 

This runs PicardTools ValidateSamFiles to detect format issues with the new BAMs. It does not compare the original and new BAM files, just checks for correct formatting. 

Adjust resources depending on your number of samples, allowing 1 hugemem CPU and 31 GB per CPU per sample, and submit with:
```
qsub Scripts/validate_bams_run_parallel.pbs
```

Output will be `<outdir>/<bam_prefix>.validate.log`. This job can be checked with the checking script (next). 

### Check the validate and flagstats jobs

Hopefully you have been checking the outputs and exit statuses of the first 3 steps of this workflow as you progress. If not, this step should hopefully alert to any issues you have missed. 

After both `compare_flagstats` and `validate_bams` are finished, run this step with:

```
bash Scripts/final_check.sh
```

A successful check will output this:
```
Parent validate job exit status OK
All validate tasks returned 0 from tool
All validate tasks exited with status 0
Parent flagstats job exit status OK
All flagstats tasks exited with status 0
Picard ValidateSamFiles and samtools flagstats detected no issues with the new BAM files.
Please proceed to make checksums.
```

If errors are detected, the specific error will be printed to screen.
 

### Checksums

If your `final_check` was successful, move on to create md5 checksums for the new BAMs. Transfer these along with your BAMs and use these to verify successful transfer to wherever you need to store them. 

Adjust resources depending on your number of samples, and submit with:
```
qsub Scripts/make_checksums_run_parallel.pbs
```

Output will be `<outdir>/<bam_prefix>.rh.bam.md5` per BAM. 
