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

- The BAM is indexed
- The data is Illumina, and has flowcell and lane as fields 3 and 4 respectively within the ':' delimited read ID
- Running on NCI Gadi HPC 
- You have sufficient disk space to temporarily hold 2 SAM files and 2 new BAM files per input BAM file. Depending on compression level of the input BAM, you may need ~15 X the disk space of the input BAM (ie 1.5 TB for a 100 GB BAM). If you are limited with processing space, you will need to adjust the workflow described below to process fewer samples at a time.  

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

Output will be `./Inputs/update_read_groups.inputs` which contains the required information for each parallel task (one per input BAM) to be executed at the next step. 

### Step 1: Convert BAM to SAM

This step simply converts the BAM to a SAM file and extracts the current headers to a separate file. This step has been separated from step 2 as it multithreads where step 2 does not.
 

Open `Scripts/update_read_groups_run_parallel_step1.pbs` and edit for your project, lstorage, and number of samples. 

7 Broadwell CPU with 9 GB RAM per CPU is a good trade-off between walltime, CPU efficiency and SU. Other resource configurations may be faster, but cost more and are less efficient. 

For example if you have 20 samples, applying 7 CPU per sample on the Broadwell normal queues, request a total of 7 x 20 CPUS = 1240 CPUs and 140 x 9 GB RAM = 1260 GB RAM. 

Walltime for one 89 GB BAM which uncompressed to 720 GB SAM was 24 minutes on the above resource settings. 

Submit step 1 with:

```
qsub Scripts/update_read_groups_run_parallel_step1.pbs
```

Output will be `<outdir>/<sample>.SAM` and <outdir>/<sample>.header` files for each sample in `./Inputs/update_read_groups.inputs`. 


### Step 2: Update read groups in the SAM and RG headers in the headers file

Open `Scripts/update_read_groups_run_parallel_step2.pbs` and edit for your project, lstorage, and number of samples. 

This step does not multi-thread. If the number of samples is greater than the number of CPU on the node (48 for Cascade Lake or 28 for Broadwell) request whole nodes. 

This step is the slowest (3.5 hours for test sample with 89 GB BAM / 720 GB SAM) so if you have a great discrepancy in the size of your input BAM files, you may have idle CPU in the job for multiple hours. You may choose to increase efficiency and reduce cost by:

- Size sorting the inputs according to BAM size (largest to smallest) and submitting with less CPU than required to run all tasks at once
- Splitting the inputs into two or more similarly sized input lists and submit them as separate jobs (change .o and .e output log file names in the PBS directives to avoid over-write) 
- Running a `for` loop over the inputs file, parsing the qsub directives on the command line, and adding a 2 second `sleep` between each iteration of the loop

To submit the job in the default way (not using any of the above 3 efficiency strategies, run:
```
qsub Scripts/update_read_groups_run_parallel_step2.pbs
```

Output will be `<outdir>/<sample>_newRG.SAM` per sample. The headers files created at step 1 will be edited in-place. 


### Step 3: Convert updated SAM to BAM

Like step 1, this step multithreads so is separated from step 2. The same CPU and mem settings for step 1 (7 Broadwell CPU with 9 GB RAM per CPU) provides the optimal tradeoff between walltime, service unis and CPU efficiency. 

Open `Scripts/update_read_groups_run_parallel_step3.pbs` and edit for your project, lstorage, and number of samples.

Submit with:

```
qsub Scripts/update_read_groups_run_parallel_step3.pbs
```
Output will be `<outdir>/<sample>.final.bam` and `<outdir>/<sample>.final.bai`. The SAM and header files created during steps 1 and 2 are removed. 


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

If something went awry with the read line printing, for example missing a read or duplicating a mapping line, samtools flagstats would emit different values. This repository includes a qick checker script that can be used to compare the outputs of `samtools flagstats`.

Note that it is not written as an `nci-parallel` job like the previous 3 steps. It is written to run from a bash `for` loop, using the sample IDs from the `Inputs/updated_read_groups_make_input.sh` file. The PBS logs will be made with the default naming method, and placed int he directory where the script is submitted. If you want to tidy this behaviour, adapt this script to the `nci-parallel` method, or simly create a directory for the PBS logs, cd into it, and submit the script from there, ensuring that the path to the script and the inputs file is update din the run command below.

Before submitting, open `Scripts/sam_flagstats.pbs` and edit the variables `outdir=<outdir>`, `old_bam_dir=<old_bam_dir>`, and `new_bam_dir=<new_bam_dir>`. Update the file name format for `old_bam` to match your input data. Save, then submit with:

```
samples=($(awk -F "," '{print $1}' ./Inputs/update_read_groups.inputs))
for sample in ${samples[@]}
	do echo $sample
	qsub -v sample="$sample" ./Scripts/sam_flagstats.pbs
	sleep 2
done
```

Output will be flagstats files in `<outdir>/new_stats` and `<outdir>/old_stats`. The script runs Linux `sdiff` on these files: if there is no difference, the exit status will be zero and no 'diff' line in the PBS .o log file. If there is a difference detected, the exit status will be 1 and there will be one or more lines of difference written to the PBS .o log file. 

Simply grep for "Exit Status" from the `stats.o*` files. 

### Checksums

Once you have checked the BAMs, as always, create checksums, and use these to verify successful transfer to wherever you need to store them. 

 

