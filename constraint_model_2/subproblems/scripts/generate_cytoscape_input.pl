#!/usr/bin/perl
## combines the results of each "divide" step and generates the adjacency graph needed for
## cytoscape visualization
##
## argument 1: the interaction matrix
## argument 2: the number of chromsomes
## argument 3: the value to scale the interaction frequency by (in order to convert it to
## an integer
## argument 4: the experimental resolution
##
## Kimberly MacKay July 4, 2016
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, USA.

use strict;
use warnings;

## check to ensure four arguments was passed in
die "ERROR: must pass in four arguments." if @ARGV != 4;

my $HiC_file = $ARGV[0];
my $num_chr = $ARGV[1];
my $scale = $ARGV[2]; # for s.pombe work this was set to 1000
my $experimental_resolution = $ARGV[3];

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
	my $input = <STDIN>;
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

###########################################################################################
##	Generate the linear interactions
##		parse the interaction matrix
###########################################################################################
## open the interaction matrix file
open HIC, "$HiC_file" or die "ERROR: $HiC_file could not be opened.";
chomp(my @matrix_file = <HIC>);
close HIC;

## $#matrix_file will be the size of the whole-genome contact map since the file also 
## contains a header line
my $num_variables = $#matrix_file;

my @frequencies;

## for each line after the header line
for(my $row = 1; $row <= $#matrix_file; $row++)
{
	## split the line
	my @matrix_line = split /\t/, $matrix_file[$row];
	
	## loop through the entire file to extract the frequencies
	## note: we only have to extract one half of the matrix since it is symmetric
	## along the diagonal
	for(my $col = $row; $col <= $num_variables; $col++)
	{
		## adjusts NA's to 0's
		if($matrix_line[$col] =~ "NA")
		{
			$frequencies[$row][$col] = 0;
		}
		else
		{
			# convert the frequency to an integer to improve speed
			$frequencies[$row][$col] = int($matrix_line[$col]*$scale);
		}
	}
}

###########################################################################################
## print out the linear interactions and their "distances" according to the 
## frequency value or experimental resolution (for s.pombe was 10 kb)
###########################################################################################
for(my $row = 1; $row <= $num_variables; $row++)
{
	my $col = $row;
	
	if($row != $stop_index[1] && $row != $stop_index[2] && $row != $stop_index[3])
	{
		if($frequencies[$row][$col] != 0)
		{
			print "bin".$row."\tbin".($row+1)."\tlinear\t".(1/$frequencies[$row][$col])."\t".$chrs[$row]."\t".$chrs[($row+1)]."\n";
		}
		else
		{
			print "bin".$row."\tbin".($row+1)."\tlinear\t".(1/$experimental_resolution)."\t".$chrs[$row]."\t".$chrs[($row+1)]."\n";
		}
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
	print STDERR "Give the path for the freq file: ";
	chomp(my $freq_file = <STDIN>);
	open FREQ_FILE, "$freq_file" or die "ERROR: $freq_file could not be opened";
	chomp(my @clp_freq_results = <FREQ_FILE>);
	close FREQ_FILE;
	
	print STDERR "Give the path for the row file: ";
	chomp(my $row_file = <STDIN>);
	open ROW_FILE, "$row_file" or die "ERROR: $row_file could not be opened";
	chomp(my @clp_row_results = <ROW_FILE>);
	close ROW_FILE;
	
	print STDERR "Give the path for the variable file: ";
	chomp(my $var_file = <STDIN>);
	open VAR_FILE, "$var_file" or die "ERROR: $var_file could not be opened";
	chomp(my @non_zero_vars = <VAR_FILE>);
	close VAR_FILE;
	
	## sanity check to make sure the files are the same length
	if($#clp_row_results == $#non_zero_vars && $#clp_row_results == $#clp_freq_results)
	{
		#print STDERR $#non_zero_vars."\n";
		
		for(my $i = 0; $i <= $#non_zero_vars; $i++)
		{	
			## get node1
			my $node1 = "bin".$non_zero_vars[$i];

			## get node(s)2
			## split the line from the clp file
			my @row_results =  split /\s+/, $clp_row_results[$i];
			
			## get the corresponding frequency
			my $freq = $clp_freq_results[$i];
	
			# if the interaction could be bound to more then one column
			for(my $j = 0; $j <= $#row_results; $j++)
			{
				if($row_results[$j] != 0)
				{
					my $node2 = "bin".$row_results[$j];

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

# for each inter-interaction
for(my $chr = 0; $chr < $num_chr; $chr++)
{
	# read in the three result files
	print STDERR "Give the path for the freq file: ";
	chomp(my $freq_file = <STDIN>);
	open FREQ_FILE, "$freq_file" or die "ERROR: $freq_file could not be opened";
	chomp(my @clp_freq_results = <FREQ_FILE>);
	close FREQ_FILE;
		
	print STDERR "Give the path for the row file: ";
	chomp(my $row_file = <STDIN>);
	open ROW_FILE, "$row_file" or die "ERROR: $row_file could not be opened";
	chomp(my @clp_row_results = <ROW_FILE>);
	close ROW_FILE;
	
	print STDERR "Give the path for the variable file: ";
	chomp(my $var_file = <STDIN>);
	open VAR_FILE, "$var_file" or die "ERROR: $var_file could not be opened";
	chomp(my @non_zero_vars = <VAR_FILE>);
	close VAR_FILE;
	
	## sanity check to make sure the files are the same length
	if($#clp_row_results == $#non_zero_vars && $#clp_row_results == $#clp_freq_results)
	{
		for(my $i = 0; $i <= $#non_zero_vars; $i++)
		{	
			## get node1
			my $node1 = "bin".$non_zero_vars[$i];

			## get the chromosome that corresponds to node 1
			my $chr1 = 0;
			
			my $found = 0;
			
			for(my $k = 1; $k <= $num_chr && !$found; $k++)
			{
				if($non_zero_vars[$i] <= $stop_index[$k])
				{
					$chr1 = $k;
					$found = 1;
				}
			}

			## get node(s)2
			## split the line from the clp file
			my @row_results =  split /\s+/, $clp_row_results[$i];
			
			## get the corresponding frequency
			my $freq = $clp_freq_results[$i];
	
			# if the interaction could be bound to more then one column
			for(my $j = 0; $j <= $#row_results; $j++)
			{
				if($row_results[$j] != 0)
				{
					my $node2 = "bin".$row_results[$j];
					
					## get the chromosome that corresponds to node 1
					my $chr2;
					$found = 0;
					for(my $l = 1; $l <= $num_chr && !$found; $l++)
					{
						if($row_results[$j] <= $stop_index[$l])
						{
							$chr2 = $l;
							$found = 1;
						}
					}

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
