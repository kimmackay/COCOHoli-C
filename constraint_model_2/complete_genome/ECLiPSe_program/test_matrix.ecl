% license: This work is licensed under the Creative Commons Attribution-NonCommercial-
% ShareAlike 3.0 Unported License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
% PO Box 1866, Mountain View, CA 94042, USA.

% Load the relevant libraries 
:- lib(gfd). 
:- lib(gfd_search).
:- lib(branch_and_bound).

% alldifferent_except(Vars) is true if each term in list Vars is  
% pairwise different from every other term, or has a value    
% of zero.  
alldifferent_except(Vars) :-  
	% the list Vars has length N 
	length(Vars, N), 

	% for each pair of terms in Vars, check if they are 
	% different or zero 
	( for(I,1,N), param(Vars, N) do 

		% The variable I is the list position of the element X 
		% in Vars. Note in ECLiPSE, indexing starts at 1. 
		element(I, Vars, X),

		( for(J, 1, N), param(Vars, I, X) do

			% The variable J is the list position of the element Y 
			% in Vars. Note in ECLiPSE, indexing starts at 1. 
			element(J, Vars, Y), 

			% If I and J do not correspond to the same position in 
			% Vars and X or Y is non-zero 
			( (I \= J) -> 

				% case 1: both X and Y are non-zero 
				% case 2: X is non-zero and Y is zero 
				% case 3: X is zero and Y is non-zero 
				% case 4: both X and Y are zero 
				 (X #\= 0 and Y #\= 0 and X #\= Y) or
				 (X #\= 0 and Y #= 0) or
				 (X #= 0 and Y #\= 0) or
				 (X #= 0 and Y #= 0) 

			; 
				true 
			) 
		) 
	).

% enforce_symmetry(Vars) ensures that for each 
% for each term in Vars bound to a non-zero   
% value the corresponding element at index 'term' in  
% Vars (ie. Var[term]) is bound to zero. This 
% ensures that 	each genomic bin can only be involved 
% in one selected interaction in the solution set.
enforce_symmetry(Vars) :- 
	% for each terms in Vars, check if it is non-zero
	( foreach(X, Vars), param(Vars) do 

	% chose a value for X from it's domain
	gfd_update,
	indomain(X, min), 

		% If the value bound to X is non-zero 
		( (X \= 0) -> 

			% The variable X is the list position of the element K 
			% in Vars. Note in ECLiPSE, indexing starts at 1. 
			element(X, Vars, K), 

			% K must be bound to zero 
			K = 0 
		; 
			true 
		) 
	). 

% maximize(Rows, Freqs)  is true if the elements in Row are all   
% different or zero, the corresponding elements in Freqs   
% are bound to zero or the associated rounded and scaled   
% integer value (based on the interaction frequency from the 
% whole-genome contact map), and the elements in Freqs  
% represent the subset of rounded and scaled interaction  
% frequencies that have maximum sum 
maximize(Rows, Freqs) :- 

	% Variable Declarations: 
	% The list Rows has one variable for each row of the 
	% whole-genome contact map 
	Rows = [Row1, Row2, Row3, Row4, Row5, Row6],

	% The list Freqs has one variable for each row of the 
	% whole-genome contact map	
	Freqs = [Freq1, Freq2, Freq3, Freq4, Freq5, Freq6],

	% Representation of the Genome: 
	% Each Row term can assume a value based on interacting bin 
	% indices where `0' represents an interaction not being 
	% selected and a non-zero value (ranging from 1 to N) 
	% represents which genomic bin is involved in the selected 
	% interaction 
	Row1 :: [0, 2..6], 
	Row2 :: [0, 3..6], 
	Row3 :: [0, 4..6], 
	Row4 :: [0, 5..6], 
	Row5 :: [0, 6], 
	Row6 :: [0], 

	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of `0' where `0'  represents an 
	% interaction not being selected 
	Freq1 :: [0, 1, 2, 3, 4, 8], 
	Freq2 :: [0, 3, 2, 5, 1], 
	Freq3 :: [0, 5, 8, 2], 
	Freq4 :: [0, 2, 9], 
	Freq5 :: [0, 4], 
	Freq6 :: [0],

	% Constraints: 
	% Each pair of corresponding (Row<i>, Freq<i>) variables 
	% must assume dependent values based on data from the 
	% whole-genome contact map; A (Row, Freq) pair ground to 
	% (0,0) encodes that nothing is chosen 

	((Row1 #= 2 and Freq1 #= 1) or
	 (Row1 #= 3 and Freq1 #= 2) or
	 (Row1 #= 4 and Freq1 #= 3) or
	 (Row1 #= 5 and Freq1 #= 4) or
	 (Row1 #= 6 and Freq1 #= 8) or
	 (Row1 #= 0 and Freq1 #= 0)), 

	((Row2 #= 3 and Freq2 #= 3) or
	 (Row2 #= 4 and Freq2 #= 2) or
	 (Row2 #= 5 and Freq2 #= 5) or
	 (Row2 #= 6 and Freq2 #= 1) or
	 (Row2 #= 0 and Freq2 #= 0)), 

	((Row3 #= 4 and Freq3 #= 5) or
	 (Row3 #= 5 and Freq3 #= 8) or
	 (Row3 #= 6 and Freq3 #= 2) or
	 (Row3 #= 0 and Freq3 #= 0)), 

	((Row4 #= 5 and Freq4 #= 2) or
	 (Row4 #= 6 and Freq4 #= 9) or
	 (Row4 #= 0 and Freq4 #= 0)), 

	((Row5 #= 6 and Freq5 #= 4) or
	 (Row5 #= 0 and Freq5 #= 0)), 

	 (Row6 #= 0 and Freq6 #= 0), 

	% All of the values assumed by the Row<i> variables must be  
	% all different or zero; multiple zeros are allowed  
	alldifferent_except(Rows), 

	% Optimize: the sum of the selected interaction 
	% frequencies is maximal (it is the minimum of 
	% the additive inverse of the sum). The predicate 
	% enforce_symmetry/1 ensures each genomic bin is 
	% involved in only one interaction in the solution set.
	Cost #= -sum(Freqs), 
	minimize(enforce_symmetry(Rows), Cost).
