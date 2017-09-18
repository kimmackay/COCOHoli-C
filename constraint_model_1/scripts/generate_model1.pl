#!/usr/bin/perl
## will generate the CLP eclispe program for constraint model 1 based on a given 
## interaction matrix
##
## argument 1: the interaction matrix
## argument 2: the value to scale the interaction frequency by (in order to convert it to
## an integer
##
## Kimberly MacKay May 10, 2015
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, USA.

use strict;
use warnings;

## check to ensure two arguments was passed in
die "ERROR: must pass in two argumnets." if @ARGV != 2;

my $HiC_file = $ARGV[0];
my $scale = $ARGV[1]; # for s.pombe work this was set to 1000

###########################################################################################
##	parse the interaction matrix
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
for(my $i = 1; $i <= $#matrix_file; $i++)
{
	## split the line
	my @matrix_line = split /\t/, $matrix_file[$i];
	
	## loop through the entire file to extract the frequencies
	for(my $j = 1; $j <= $#matrix_file; $j++)
	{
		## adjusts NA's to 0's
		if($matrix_line[$j] =~ "NA")
		{
			$frequencies[$i][$j] = 0;
		}
		else
		{
			# convert the frequency to an integer to improve speed
			$frequencies[$i][$j] = int($matrix_line[$j]*$scale);
		}
	}
}

###########################################################################################
##	printing the program
###########################################################################################

print "% Load the relevant libraries \n";
print ":- lib(gfd). \n\n";

print "% maximize(Threshold, ContactMap) is true if \n";
print "% for each detected Hi-C interaction within ContactMap \n";
print "% involving a specific genomic bin at least Threshold \n";
print "% terms are ground to a zero value and the non-zero \n";
print "% terms represent the subset with maximal sum of rounded \n"; 
print "% and scaled interaction frequencies. Threshold must be \n"; 
print "% bound to a value between 0 and N. \n";
print "maximize(Threshold, ContactMap) :- \n\n";

print "\t% Variable Declarations: \n";
print "\t% The list ContactMap has one variable for each row and \n";
print "\t% column of the whole-genome contact map \n";

# print each variable that represents the contact map
print "\tContactMap = [";
for(my $i = 1; $i <= $num_variables; $i++)
{
	for(my $j = 1; $j <= $num_variables; $j++)
	{
		# account for the last variable, it won't have a comma
		if($i == $num_variables && $j == $num_variables)
		{
			print "V".$i."_".$j."],\n\n";
		}
		else
		{
			print "V".$i."_".$j.", ";
		}
	}
}

print "\t% Each Bin<i> involves all the contact map information for \n";
print "\t% row <i> and for column <i> \n";
for(my $i = 1; $i <= $num_variables; $i++)
{
	print "\tBin".$i."=[";
	for(my $j = 1; $j <= $num_variables; $j++)
	{
			## if it is the last set
			if($j == $num_variables)
			{
				if($i != $j)
				{
					print "V".$i."_".$j.", V".$j."_".$i."],\n";
				}
				else
				{
					print "V".$i."_".$j."],\n";	
				}	
			}
			else
			{
				if($i != $j)
				{
					print "V".$i."_".$j.", V".$j."_".$i.", ";
				}
				else
				{
					print "V".$i."_".$j.", ";	
				}	
			}
	}	
}

print "\n\t% Representation of the Genome: \n";
print "\t% Each cell (V<i>_<j>) of the ContactMap can assume \n";
print "\t% the rounded and scaled integer value based on the \n";
print "\t% interaction frequency or `0' where `0'  represents an \n";
print "\t% interaction not being selected. Elements along the \n";
print "\t% diagonal of ContactMap assume a value of `0' to \n";
print "\t% prevent a region interacting with itself \n";

for(my $i = 1; $i <= $num_variables; $i++)
{
	for(my $j = 1; $j <= $num_variables; $j++)
	{
		# account for the self-self interaction
		if($i == $j || $frequencies[$i][$j] == 0)
		{
			print "\tV".$i."_".$j." :: [0],\n";
		}
		else
		{
			print "\tV".$i."_".$j." :: [0, ".$frequencies[$i][$j]."],\n";
		}
	}
}

print "\n\t% Constraints: \n";
print "\t% Each genomic bin has at least Threshold interactions that \n";
print "\t% are zero, meaning that each bin can be actively involve in \n";
print "\t% only N minus Threshold detected (i.e. nonzero) Hi-C \n";
print "\t% interactions (where N is the number of genomic bins in the \n";
print "\t% original whole-genome contact map) \n";
for(my $i = 1; $i <= $num_variables; $i++)
{
	print "\tatleast(Threshold, Bin".$i.", 0), \n";
}

print "\n\t% Optimize: the sum of the selected interaction \n";
print "\t% frequencies is maximal (it is the minimum of \n";
print "\t% the additive inverse of the sum) \n";
print "\tCost #= -sum(ContactMap), \n";
print "\tsearch(ContactMap, 0, input_order, indomain_max, bb_min(Cost), []). \n";
