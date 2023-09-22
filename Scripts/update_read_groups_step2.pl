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
		my $full_rg = "\@RG\tID\:$rg_id\tPL\:$platform\tPU\:$flowcell\.$lane\tSM\:$sample\tLB\:$sample\_$lib\tCN\:$centre"; 
		$rgidhash->{$old_rg_id}->{full_rg} = $full_rg;
	} 
	else { # A header - will be updated at step 3 with samtools reheader (can't update @RG headers until all read groups are known)
		print O "$line\n";	
	}
	 
} close S; close O; 

# Update headers file with new read group metadata: 
foreach my $old_rg_id (sort keys %{$rgidhash}) {
	my $new_id = $rgidhash->{$old_rg_id}->{full_rg}; 
	`sed -i 's/.*ID\:$old_rg_id.*/$new_id/'  $header`; 
}
