#!/usr/bin/perl
## combines the results of each "divide" step and generates the adjacency graph needed for
## cytoscape visualization
##
## argument 2: the number of chromsomes
## argument 4: the experimental resolution
##
## Kimberly MacKay July 4, 2016
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, USA.

use strict;
use warnings;

## check to ensure two arguments was passed in
die "ERROR: must pass in two arguments." if @ARGV != 2;

my $num_chr = $ARGV[0];
my $experimental_resolution = $ARGV[1];

# scan each of the results files and collect the result in this hash table
my %interactions;
my %dynamics_coefficents;

###########################################################################################
##	initialize variables
###########################################################################################
## get the ending index of each chromosome
my @stop_index;
for(my $chr = 1; $chr <= $num_chr; $chr++)
{
	print STDERR "What is the ending index of CHR".$chr."? ";
	chomp(my $input = <STDIN>);
	$stop_index[$chr] = int($input);
}

$stop_index[0] = 0;

## make an array that maps genomic bin to the corresponding chromosome
my @chrs;
my $bin = 1;
for(my $j = 1; $j <= $num_chr; $j++)
{
	for(my $i = $bin; $i <= $stop_index[$j]; $i++)
	{
			$chrs[$bin] = $j;
			$bin = $bin +1;
	}
}

my $num_variables = $stop_index[-1];

## print the header of the output file
print "source_node \t target_node \t type_of_interaction \t frequency \t node1_chr \t node2_chr \n";

###########################################################################################
## print out the linear interactions and their "distances" according to the 
## frequency value or experimental resolution (for s.pombe was 10 kb)
###########################################################################################
for(my $row = 1; $row <= $num_variables; $row++)
{
	
	if($row < $stop_index[1])
	{
		print "bin".$row."\tbin".($row+1)."\tlinear1\t".(1/$experimental_resolution)."\t".$chrs[$row]."\t".$chrs[($row+1)]."\n";
	}
	elsif($row < $stop_index[2] && $row != $stop_index[1])
	{
		print "bin".$row."\tbin".($row+1)."\tlinear2\t".(1/$experimental_resolution)."\t".$chrs[$row]."\t".$chrs[($row+1)]."\n";
	}
	elsif($row < $stop_index[3] && $row != $stop_index[1] && $row != $stop_index[2])
	{
		print "bin".$row."\tbin".($row+1)."\tlinear3\t".(1/$experimental_resolution)."\t".$chrs[$row]."\t".$chrs[($row+1)]."\n";
	}
	
	# initialize the %dynamics_coefficients hash for each bin to be = 0
	$dynamics_coefficents{"bin".$row} = 0;
}

###########################################################################################
##	Merge the results of the divide and conquer
###########################################################################################
# gut check - uncomment to see the file ordering based on glob
#print STDERR (glob("$results_dir/*.txt"));

# for each intra-interaction
for(my $chr = 0; $chr < $num_chr; $chr++)
{
	# read in the three result files
	print STDERR "cis-interaction: give the path for the freq file: ";
	chomp(my $freq_file = <STDIN>);
	open FREQ_FILE, "$freq_file" or die "ERROR: $freq_file could not be opened";
	chomp(my @clp_freq_results = <FREQ_FILE>);
	close FREQ_FILE;
	
	print STDERR "cis-interaction: give the path for the row file: ";
	chomp(my $row_file = <STDIN>);
	open ROW_FILE, "$row_file" or die "ERROR: $row_file could not be opened";
	chomp(my @clp_row_results = <ROW_FILE>);
	close ROW_FILE;
	
	## sanity check to make sure the files are the same length
	if($#clp_row_results == $#clp_freq_results)
	{
		for(my $i = 0; $i <= $#clp_row_results; $i++)
		{	
		
			## if it is a non-zero interaction
			if($clp_freq_results[$i] != 0)
			{
				## get the value of bin 1 - note it needs to be offest by the chr start index
				my $node1 = "bin".(($i+1) + $stop_index[$chr]);
				
				## get the corresponding frequency
				my $freq = $clp_freq_results[$i];
				
				## get node(s)2
				## split the line from the clp file
				my @row_results =  split /\s+/, ($clp_row_results[$i]);
		
				# if the interaction could be bound to one or more column(s)
				for(my $j = 0; $j <= $#row_results; $j++)
				{
					my $node2 = "bin".($row_results[$j] + $stop_index[$chr]);
					
					# add it to the hash
					$interactions{"cis_".$node1."_".$node2}{NODE1_BIN} = $node1;
					$interactions{"cis_".$node1."_".$node2}{NODE1_CHR} = ($chr+1);
					
					$interactions{"cis_".$node1."_".$node2}{NODE2_BIN} = $node2;
					$interactions{"cis_".$node1."_".$node2}{NODE2_CHR} = ($chr+1);
					
					$interactions{"cis_".$node1."_".$node2}{INTERACTION_FREQ} = $freq;
					$interactions{"cis_".$node1."_".$node2}{INTERACTION_TYPE} = "cis";
									
					$dynamics_coefficents{$node1} = $dynamics_coefficents{$node1} + 1;
					$dynamics_coefficents{$node2} = $dynamics_coefficents{$node2} + 1;	
				}
			}
		}
	}
	else
	{
		print "ERROR: files are not the same length and they should be.";
	}
}

# holds the stop index that needs to be added to the second chromosomes bin number
my @temp = [558, 558, 1012];

# for each inter-interaction
for(my $chr = 0; $chr < $num_chr; $chr++)
{
	# read in the two result files
	print STDERR "trans-interaction: give the path for the freq file: ";
	chomp(my $freq_file = <STDIN>);
	open FREQ_FILE, "$freq_file" or die "ERROR: $freq_file could not be opened";
	chomp(my @clp_freq_results = <FREQ_FILE>);
	close FREQ_FILE;
		
	print STDERR "trans-interaction: give the path for the row file: ";
	chomp(my $row_file = <STDIN>);
	open ROW_FILE, "$row_file" or die "ERROR: $row_file could not be opened";
	chomp(my @clp_row_results = <ROW_FILE>);
	close ROW_FILE;
	
	# get the first chromosome of the trans interaction
	print STDERR "What is the first chormosome? ";
	chomp(my $chr1 = <STDIN>);
	
	# get the second chromosome of the trans interaction
	print STDERR "What is the second chormosome? ";
	chomp(my $chr2 = <STDIN>);
	
	## sanity check to make sure the files are the same length
	if($#clp_row_results == $#clp_freq_results)
	{
		for(my $i = 0; $i <= $#clp_row_results; $i++)
		{	
			## if it is a non-zero interaction
			if($clp_freq_results[$i] != 0)
			{
				## get the value of bin 1 - note it needs to be offest by the chr start index
				my $node1 = "bin".(($i+1) + $stop_index[$chr1-1]);
				
				## get the corresponding frequency
				my $freq = $clp_freq_results[$i];
				
				## get node(s)2
				## split the line from the clp file
				my @row_results =  split /\s+/, $clp_row_results[$i];
				
				# if the interaction could be bound to one or more column(s)
				for(my $j = 0; $j <= $#row_results; $j++)
				{
					my $node2 = "bin".($row_results[$j] + $stop_index[$chr2-1]);
					
					# add it to the hash
					$interactions{"trans_".$node1."_".$node2}{NODE1_BIN} = $node1;
					$interactions{"trans_".$node1."_".$node2}{NODE1_CHR} = $chr1;
					
					$interactions{"trans_".$node1."_".$node2}{NODE2_BIN} = $node2;
					$interactions{"trans_".$node1."_".$node2}{NODE2_CHR} = $chr2;
					
					$interactions{"trans_".$node1."_".$node2}{INTERACTION_FREQ} = $freq;
					$interactions{"trans_".$node1."_".$node2}{INTERACTION_TYPE} = "trans";
									
					$dynamics_coefficents{$node1} = $dynamics_coefficents{$node1} + 1;
					$dynamics_coefficents{$node2} = $dynamics_coefficents{$node2} + 1;
				}
			}
		}
	}
	else
	{
		print "ERROR: files are not the same length and they should be.";
	}
}

###########################################################################################
## print out the cis and trans interactions and their "distances" according to the 
## frequency value and dynamics coefficient
###########################################################################################

foreach my $selected_interaction (sort keys %interactions) 
{
		# calculate the average dynamics coefficient between the two bins involved in the interactions
		my $d_coefficient = ($dynamics_coefficents{$interactions{$selected_interaction}{NODE1_BIN}} + $dynamics_coefficents{$interactions{$selected_interaction}{NODE2_BIN}})/2;

		# pint the interaction	
		print $interactions{$selected_interaction}{NODE1_BIN}."\t".$interactions{$selected_interaction}{NODE2_BIN}."\t".$interactions{$selected_interaction}{INTERACTION_TYPE}."\t".($d_coefficient/$interactions{$selected_interaction}{INTERACTION_FREQ})."\t".$interactions{$selected_interaction}{NODE1_CHR}."\t".$interactions{$selected_interaction}{NODE2_CHR}."\n";
}	
