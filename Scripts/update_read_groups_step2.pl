#!/usr/bin/env perl

#########################################################
#
# Platform: NCI Gadi HPC
# Description:  Update SAM read group IDs and RG headers
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

use warnings;
use strict; 

my $inputs = $ARGV[0];
my ($sample, $bam, $outdir, $lib, $centre, $platform, @rest) = split(',', $inputs); 

my @patharray = split('/', $bam); 
my $prefix = $patharray[-1]; 
$prefix =~ s/\.bam//;  

my $sam = "$outdir\/$prefix\.sam";  
my $new_sam = "$outdir\/$prefix\_newRG.sam";
my $header = "$outdir\/$prefix\.header";
my $new_header = "$outdir\/$prefix\_new.header"; 

my $rgidhash = {}; 

open (S, $sam) || die "$! $sam\n";
open (O, ">$new_sam") || die "$! write $new_sam\n";

my $m = 0; 
 
while (my $line = <S> ) {
	chomp $line;
	
	if ($line !~m/^\@/) { # Not a header
		# Collect the read group info:
		my @cols = split('\t', $line);
		my $read_id = $cols[0]; 
		my @ids = split(':', $read_id); 
		my $flowcell = $ids[2];
		my $lane = $ids[3]; 
	
		# Replace the read group ID in the SAM line:
		my $rg_id = "$flowcell\.$lane\.$sample\_$lib";
		$line =~ s/RG\:Z\:(\S+)/RG\:Z\:$rg_id/;
		my $old_rg_id = $1; 
	
		# Print updated line to new SAM:
		print O "$line\n"; 
	
		# Store the read group metadata for later replacing in the header: 
		my $new_rg_id = "\@RG\tID\:$rg_id\tPL\:$platform\tPU\:$flowcell\.$lane\tSM\:$sample\tLB\:$sample\_$lib\tCN\:$centre"; 
		$rgidhash->{$new_rg_id}->{old_rg_id} = $old_rg_id;
		
	} 
	else { # A header - will be updated at step 3 with samtools reheader (can't update @RG headers until all read groups are known)
		print O "$line\n";	
	}
	 
} close S; close O; 

# Replace read groups in headers
# This method ensures that all read groups identified from the reads are written to headers
# Direct find and replace proved fallible when the original BAM had fewer RG headers than correct
open (H, $header) || die "$! $header\n";
open (NH, ">$new_header") || die "$! $new_header\n";
my $done_rg = 0; 

while (my $line = <H>) {
	chomp $line; 
	if ($line =~ m/^\@RG/) { 
		if ( ! $done_rg) {
			foreach my $new_rg_id (sort keys %{$rgidhash}) {
				print NH  "$new_rg_id\n"; 
			}
			$done_rg++; 
		}				
	}
	else {
		print NH "$line\n"; 
	} 
} close H; close NH; 

`mv $new_header $header`;
