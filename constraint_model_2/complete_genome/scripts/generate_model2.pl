#!/usr/bin/perl
## will generate the CLP eclispe program for constraint model 2 based on a given 
## interaction matrix
##
## argument 1: the interaction matrix
## argument 2: the value to scale the interaction frequency by (in order to convert it to
## an integer
##
## Kimberly MacKay June 8, 2016
## updated: July 28, 2017; removal of residual and redundant code
##
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, USA.

use strict;
use warnings;
use List::MoreUtils 'uniq';

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
##	printing the program
###########################################################################################
print "% license: This work is licensed under the Creative Commons Attribution-NonCommercial-\n";
print "% ShareAlike 3.0 Unported License. To view a copy of this license, visit \n";
print "% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, \n";
print "% PO Box 1866, Mountain View, CA 94042, USA.\n\n";

print "% Load the relevant libraries \n";
print ":- lib(gfd). \n";
print ":- lib(branch_and_bound).\n\n";

print "% alldifferent_except(Vars) is true if each term in list Vars is  \n";
print "% pairwise different from every other term, or has a value    \n";
print "% of zero.  \n";
print "alldifferent_except(Vars) :-  \n";
	print "\t% the list Vars has length N \n";
	print "\tlength(Vars, N), \n\n";

	print "\t% for each pair of terms in Vars, check if they are \n";
	print "\t% different or zero \n";
	print "\t( for(I,1,N), param(Vars, N) do \n\n";
	
		print "\t\t% The variable I is the list position of the element X \n";
		print "\t\t% in Vars. Note in ECLiPSE, indexing starts at 1. \n";
		print "\t\telement(I, Vars, X),\n\n";
	
		print "\t\t( for(J, 1, N), param(Vars, I, X) do\n\n";

			print "\t\t\t% The variable J is the list position of the element Y \n";
			print "\t\t\t% in Vars. Note in ECLiPSE, indexing starts at 1. \n"; 
			print "\t\t\telement(J, Vars, Y), \n\n";

			print "\t\t\t% If I and J do not correspond to the same position in \n";
			print "\t\t\t% Vars and X or Y is non-zero \n";
			print "\t\t\t( (I \\= J) -> \n\n";
				
				print "\t\t\t\t% case 1: both X and Y are non-zero \n";
				print "\t\t\t\t% case 2: X is non-zero and Y is zero \n";
				print "\t\t\t\t% case 3: X is zero and Y is non-zero \n";
				print "\t\t\t\t% case 4: both X and Y are zero \n";
				print "\t\t\t\t (X #\\= 0 and Y #\\= 0 and X #\\= Y) or\n";
				print "\t\t\t\t (X #\\= 0 and Y #= 0) or\n";
				print "\t\t\t\t (X #= 0 and Y #\\= 0) or\n";
				print "\t\t\t\t (X #= 0 and Y #= 0) \n\n";
		
			print "\t\t\t; \n";
				print "\t\t\t\ttrue \n";
			print "\t\t\t) \n";
		print "\t\t) \n";
	print "\t).\n\n";

print "% enforce_symmetry(Vars) ensures that for each \n";
print "% for each term in Vars bound to a non-zero   \n";
print "% value the corresponding element at index 'term' in  \n";
print "% Vars (ie. Var[term]) is bound to zero. This \n";
print "% ensures that 	each genomic bin can only be involved \n";
print "% in one selected interaction in the solution set.\n";
print "enforce_symmetry(Vars) :- \n";
	print "\t% for each terms in Vars, check if it is non-zero\n";
	print "\t( foreach(X, Vars), param(Vars) do \n\n";
		
		print "\t% chose a value for X from it's domain\n";
		print "\tgfd_update,\n";
		print "\tindomain(X), \n\n";

		print "\t\t% If the value bound to X is non-zero \n";
		print "\t\t( (X \\= 0) -> \n\n";
	
			print "\t\t\t% The variable X is the list position of the element K \n";
			print "\t\t\t% in Vars. Note in ECLiPSE, indexing starts at 1. \n";
			print "\t\t\telement(X, Vars, K), \n\n";
		
			print "\t\t\t% K must be bound to zero \n";
			print "\t\t\tK = 0 \n";
		print "\t\t; \n";
			print "\t\t\ttrue \n";
		print "\t\t) \n";
	print "\t). \n\n";

print "% maximize(Rows, Freqs)  is true if the elements in Row are all   \n";
print "% different or zero, the corresponding elements in Freqs   \n";
print "% are bound to zero or the associated rounded and scaled   \n";
print "% integer value (based on the interaction frequency from the \n";
print "% whole-genome contact map), and the elements in Freqs  \n";
print "% represent the subset of rounded and scaled interaction  \n";
print "% frequencies that have maximum sum \n";
print "maximize(Rows, Freqs) :- \n\n";

	print "\t% Variable Declarations: \n";
	print "\t% The list Rows has one variable for each row of the \n";
	print "\t% whole-genome contact map \n";
	print "\tRows = [";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		# account for the last variable, it won't have a comma
		if($i == $num_variables)
		{
			print "Row".$i."],\n\n";
		}
		else
		{
			print "Row".$i.", ";
		}
	}
	
	print "\t% The list Freqs has one variable for each row of the \n";
	print "\t% whole-genome contact map	\n";
	print "\tFreqs = [";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		# account for the last variable, it won't have a comma
		if($i == $num_variables)
		{
			print "Freq".$i."],\n\n";
		}
		else
		{
			print "Freq".$i.", ";
		}
	}

	print "\t% Representation of the Genome: \n";
	print "\t% Each Row term can assume a value based on interacting bin \n";
	print "\t% indices where `0' represents an interaction not being \n";
	print "\t% selected and a non-zero value (ranging from 1 to N) \n";
	print "\t% represents which genomic bin is involved in the selected \n";
	print "\t% interaction \n";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		# account for the last number
		if($i == $num_variables)
		{
			print "\tRow".$i." :: [0], \n\n";
		}
		elsif($i+1 == $num_variables)
		{
			print "\tRow".$i." :: [0, ".($i+1)."], \n";
		}
		else
		{
			print "\tRow".$i." :: [0, ".($i+1)."..".$num_variables."], \n";
		}
	}
	
	print "\t% Each frequency term can assume either the rounded \n";
	print "\t% and scaled integer value (based on the corresponding \n";
	print "\t% interaction frequency from the whole-genome contact \n";
	print "\t% map) or a value of `0' where `0'  represents an \n";
	print "\t% interaction not being selected \n";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
	
		# account for the last number
		if($i == $num_variables)
		{
			print "\tFreq".$i." :: [0],";
		}
		else
		{
			print "\tFreq".$i." :: [0, ";
	
			# note we only have to loop through half of the matrix since it is symmetric
			for(my $j = $i+1; $j <= $num_variables; $j++)
			{
				# account for the last number, it won't have a comma
				if($j == $num_variables)
				{
					print $frequencies[$i][$j]."], \n";
				}
				else
				{
					print $frequencies[$i][$j].", ";
				}
			}
		}
	}

	print "\n\n\t% Constraints: \n";
	print "\t% Each pair of corresponding (Row<i>, Freq<i>) variables \n";
	print "\t% must assume dependent values based on data from the \n";
	print "\t% whole-genome contact map; A (Row, Freq) pair ground to \n";
	print "\t% (0,0) encodes that nothing is chosen \n\n";
	for(my $i = 1; $i <= $num_variables; $i++)
	{
		# account for the last number
		if($i == $num_variables)
		{
			print "\t (Row".$i." #= 0 and Freq".$i." #= 0), \n\n";
		}
		else
		{
		
			# print the first conjuctions
			my $j = $i + 1;
			print "\t((Row".$i." #= ".$j." and Freq".$i." #= ".$frequencies[$i][$j].") or\n";
			
			# if it is the last variable, print the no case
			if($j >= $num_variables)
			{
				print "\t (Row".$i." #= 0 and Freq".$i." #= 0)), \n\n"; 
			}
		
			# note we only have to loop through half of the matrix since it is symmetric
			for($j = $i+2; $j <= $num_variables; $j++)
			{
				print "\t (Row".$i." #= ".$j." and Freq".$i." #= ".$frequencies[$i][$j].") or\n";
				
				# if it is the last variable, print the no case
				if($j >= $num_variables)
				{
					print "\t (Row".$i." #= 0 and Freq".$i." #= 0)), \n\n"; 
				}
			}

		}
	}

	print "\t% All of the values assumed by the Row<i> variables must be  \n";
	print "\t% all different or zero; multiple zeros are allowed  \n";
	print "\talldifferent_except(Rows), \n\n";

	print "\t% Optimize: the sum of the selected interaction \n";
	print "\t% frequencies is maximal (it is the minimum of \n";
	print "\t% the additive inverse of the sum). The predicate \n";
	print "\t% enforce_symmetry/1 ensures each genomic bin is \n";
	print "\t% involved in only one interaction in the solution set.\n";
	print "\tCost #= -sum(Freqs), \n";
	print "\tminimize(enforce_symmetry(Rows), Cost).\n";
