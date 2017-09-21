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
			K #= 0 
		; 
			true 
		) 
	). 

% maximize(RowFile, FreqFile, Rows) is true if the elements in Rows are all  
% different or zero, the corresponding elements in Freqs  
% are bound to zero or the associated rounded and scaled  
% integer value (based on the interaction frequency from the 
% whole-genome contact map), and the elements in Freqs  
% represent the subset of rounded and scaled interaction  
% frequencies that have maximum sum 
maximize(RowFile, FreqFile, Rows) :-

	% Variable Declarations: 
	% The list Rows has one variable for each row of the 
	% whole-genome contact map 
	Rows = [Row1, Row2, Row3, Row4, Row5, Row6, Row7, Row8, Row9, Row10, Row11, Row12, Row13, Row14, Row15, Row16, Row17, Row18, Row19, Row20, Row21, Row22, Row23, Row24, Row25],

	% The list Freqs has one variable for each row of the 
	% whole-genome contact map	
	Freqs = [Freq1, Freq2, Freq3, Freq4, Freq5, Freq6, Freq7, Freq8, Freq9, Freq10, Freq11, Freq12, Freq13, Freq14, Freq15, Freq16, Freq17, Freq18, Freq19, Freq20, Freq21, Freq22, Freq23, Freq24, Freq25],

	% Representation of the Genome: 
	% Each Row term can assume a value based on interacting bin 
	% indices where `0' represents an interaction not being 
	% selected and a non-zero value (ranging from 1 to N) 
	% represents which genomic bin is involved in the selected 
	% interaction 
	Row1 :: [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row2 :: [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row3 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row4 :: [0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row5 :: [0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row6 :: [0, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row7 :: [0, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row8 :: [0, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row9 :: [0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row10 :: [0, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row11 :: [0, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row12 :: [0, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row13 :: [0, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row14 :: [0, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row15 :: [0, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row16 :: [0, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row17 :: [0, 18, 19, 20, 21, 22, 23, 24, 25],
	Row18 :: [0, 19, 20, 21, 22, 23, 24, 25],
	Row19 :: [0, 20, 21, 22, 23, 24, 25],
	Row20 :: [0, 21, 22, 23, 24, 25],
	Row21 :: [0, 22, 23, 24, 25],
	Row22 :: [0, 23, 24, 25],
	Row23 :: [0, 24, 25],
	Row24 :: [0, 25],
	Row25 :: [0],

	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of `0' where `0'  represents an 
	% interaction not being selected 
	Freq1 :: [0, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25],
	Freq2 :: [0, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 2, 3],
	Freq3 :: [0, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 2, 3, 4, 5],
	Freq4 :: [0, 11, 13, 15, 17, 19, 21, 23, 25, 2, 3, 4, 5, 7, 9],
	Freq5 :: [0, 15, 17, 19, 21, 23, 25, 2, 3, 4, 5, 7, 9, 11, 13],
	Freq6 :: [0, 19, 21, 23, 25, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17],
	Freq7 :: [0, 22, 24, 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25],
	Freq8 :: [0, 25, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23],
	Freq9 :: [0, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25],
	Freq10 :: [0, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 2, 3],
	Freq11 :: [0, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 2, 3, 4],
	Freq12 :: [0, 9, 11, 13, 15, 17, 19, 21, 23, 25, 2, 3, 4, 5],
	Freq13 :: [0, 11, 13, 15, 17, 19, 21, 23, 25, 2, 3, 4, 5],
	Freq14 :: [0, 15, 17, 19, 21, 23, 25, 2, 3, 4, 5, 7],
	Freq15 :: [0, 17, 19, 21, 23, 25, 2, 3, 4, 5, 7],
	Freq16 :: [0, 21, 23, 25, 2, 3, 4, 5, 7, 9],
	Freq17 :: [0, 23, 25, 2, 3, 4, 5, 7, 9],
	Freq18 :: [0, 2, 3, 4, 5, 7, 9, 11],
	Freq19 :: [0, 3, 4, 5, 7, 9, 11],
	Freq20 :: [0, 5, 7, 9, 11, 13],
	Freq21 :: [0, 7, 9, 11, 13],
	Freq22 :: [0, 11, 13, 15],
	Freq23 :: [0, 13, 15],
	Freq24 :: [0, 17],
	Freq25 :: [0],

	% Constraints: 
	% Each pair of corresponding (Row<i>, Freq<i>) variables 
	% must assume dependent values based on data from the 
	% whole-genome contact map; A (Row, Freq) pair ground to 
	% (0,0) encodes that nothing is chosen 

	((Row1 #= 2 and Freq1 #= 2) or
	 (Row1 #= 3 and Freq1 #= 3) or
	 (Row1 #= 4 and Freq1 #= 4) or
	 (Row1 #= 5 and Freq1 #= 5) or
	 (Row1 #= 6 and Freq1 #= 7) or
	 (Row1 #= 7 and Freq1 #= 9) or
	 (Row1 #= 8 and Freq1 #= 11) or
	 (Row1 #= 9 and Freq1 #= 13) or
	 (Row1 #= 10 and Freq1 #= 15) or
	 (Row1 #= 11 and Freq1 #= 17) or
	 (Row1 #= 12 and Freq1 #= 19) or
	 (Row1 #= 13 and Freq1 #= 21) or
	 (Row1 #= 14 and Freq1 #= 23) or
	 (Row1 #= 15 and Freq1 #= 25) or
	 (Row1 #= 16 and Freq1 #= 2) or
	 (Row1 #= 17 and Freq1 #= 3) or
	 (Row1 #= 18 and Freq1 #= 4) or
	 (Row1 #= 19 and Freq1 #= 5) or
	 (Row1 #= 20 and Freq1 #= 7) or
	 (Row1 #= 21 and Freq1 #= 9) or
	 (Row1 #= 22 and Freq1 #= 11) or
	 (Row1 #= 23 and Freq1 #= 13) or
	 (Row1 #= 24 and Freq1 #= 15) or
	 (Row1 #= 25 and Freq1 #= 17) or
	 (Row1 #= 0 and Freq1 #= 0)), 

	((Row2 #= 3 and Freq2 #= 4) or
	 (Row2 #= 4 and Freq2 #= 5) or
	 (Row2 #= 5 and Freq2 #= 7) or
	 (Row2 #= 6 and Freq2 #= 9) or
	 (Row2 #= 7 and Freq2 #= 11) or
	 (Row2 #= 8 and Freq2 #= 13) or
	 (Row2 #= 9 and Freq2 #= 15) or
	 (Row2 #= 10 and Freq2 #= 17) or
	 (Row2 #= 11 and Freq2 #= 19) or
	 (Row2 #= 12 and Freq2 #= 21) or
	 (Row2 #= 13 and Freq2 #= 23) or
	 (Row2 #= 14 and Freq2 #= 25) or
	 (Row2 #= 15 and Freq2 #= 2) or
	 (Row2 #= 16 and Freq2 #= 3) or
	 (Row2 #= 17 and Freq2 #= 4) or
	 (Row2 #= 18 and Freq2 #= 5) or
	 (Row2 #= 19 and Freq2 #= 7) or
	 (Row2 #= 20 and Freq2 #= 9) or
	 (Row2 #= 21 and Freq2 #= 11) or
	 (Row2 #= 22 and Freq2 #= 13) or
	 (Row2 #= 23 and Freq2 #= 15) or
	 (Row2 #= 24 and Freq2 #= 17) or
	 (Row2 #= 25 and Freq2 #= 19) or
	 (Row2 #= 0 and Freq2 #= 0)), 

	((Row3 #= 4 and Freq3 #= 7) or
	 (Row3 #= 5 and Freq3 #= 9) or
	 (Row3 #= 6 and Freq3 #= 11) or
	 (Row3 #= 7 and Freq3 #= 13) or
	 (Row3 #= 8 and Freq3 #= 15) or
	 (Row3 #= 9 and Freq3 #= 17) or
	 (Row3 #= 10 and Freq3 #= 19) or
	 (Row3 #= 11 and Freq3 #= 21) or
	 (Row3 #= 12 and Freq3 #= 23) or
	 (Row3 #= 13 and Freq3 #= 25) or
	 (Row3 #= 14 and Freq3 #= 2) or
	 (Row3 #= 15 and Freq3 #= 3) or
	 (Row3 #= 16 and Freq3 #= 4) or
	 (Row3 #= 17 and Freq3 #= 5) or
	 (Row3 #= 18 and Freq3 #= 7) or
	 (Row3 #= 19 and Freq3 #= 9) or
	 (Row3 #= 20 and Freq3 #= 11) or
	 (Row3 #= 21 and Freq3 #= 13) or
	 (Row3 #= 22 and Freq3 #= 15) or
	 (Row3 #= 23 and Freq3 #= 17) or
	 (Row3 #= 24 and Freq3 #= 19) or
	 (Row3 #= 25 and Freq3 #= 21) or
	 (Row3 #= 0 and Freq3 #= 0)), 

	((Row4 #= 5 and Freq4 #= 11) or
	 (Row4 #= 6 and Freq4 #= 13) or
	 (Row4 #= 7 and Freq4 #= 15) or
	 (Row4 #= 8 and Freq4 #= 17) or
	 (Row4 #= 9 and Freq4 #= 19) or
	 (Row4 #= 10 and Freq4 #= 21) or
	 (Row4 #= 11 and Freq4 #= 23) or
	 (Row4 #= 12 and Freq4 #= 25) or
	 (Row4 #= 13 and Freq4 #= 2) or
	 (Row4 #= 14 and Freq4 #= 3) or
	 (Row4 #= 15 and Freq4 #= 4) or
	 (Row4 #= 16 and Freq4 #= 5) or
	 (Row4 #= 17 and Freq4 #= 7) or
	 (Row4 #= 18 and Freq4 #= 9) or
	 (Row4 #= 19 and Freq4 #= 11) or
	 (Row4 #= 20 and Freq4 #= 13) or
	 (Row4 #= 21 and Freq4 #= 15) or
	 (Row4 #= 22 and Freq4 #= 17) or
	 (Row4 #= 23 and Freq4 #= 19) or
	 (Row4 #= 24 and Freq4 #= 21) or
	 (Row4 #= 25 and Freq4 #= 23) or
	 (Row4 #= 0 and Freq4 #= 0)), 

	((Row5 #= 6 and Freq5 #= 15) or
	 (Row5 #= 7 and Freq5 #= 17) or
	 (Row5 #= 8 and Freq5 #= 19) or
	 (Row5 #= 9 and Freq5 #= 21) or
	 (Row5 #= 10 and Freq5 #= 23) or
	 (Row5 #= 11 and Freq5 #= 25) or
	 (Row5 #= 12 and Freq5 #= 2) or
	 (Row5 #= 13 and Freq5 #= 3) or
	 (Row5 #= 14 and Freq5 #= 4) or
	 (Row5 #= 15 and Freq5 #= 5) or
	 (Row5 #= 16 and Freq5 #= 7) or
	 (Row5 #= 17 and Freq5 #= 9) or
	 (Row5 #= 18 and Freq5 #= 11) or
	 (Row5 #= 19 and Freq5 #= 13) or
	 (Row5 #= 20 and Freq5 #= 15) or
	 (Row5 #= 21 and Freq5 #= 17) or
	 (Row5 #= 22 and Freq5 #= 19) or
	 (Row5 #= 23 and Freq5 #= 21) or
	 (Row5 #= 24 and Freq5 #= 23) or
	 (Row5 #= 25 and Freq5 #= 25) or
	 (Row5 #= 0 and Freq5 #= 0)), 

	((Row6 #= 7 and Freq6 #= 19) or
	 (Row6 #= 8 and Freq6 #= 21) or
	 (Row6 #= 9 and Freq6 #= 23) or
	 (Row6 #= 10 and Freq6 #= 25) or
	 (Row6 #= 11 and Freq6 #= 2) or
	 (Row6 #= 12 and Freq6 #= 3) or
	 (Row6 #= 13 and Freq6 #= 4) or
	 (Row6 #= 14 and Freq6 #= 5) or
	 (Row6 #= 15 and Freq6 #= 7) or
	 (Row6 #= 16 and Freq6 #= 9) or
	 (Row6 #= 17 and Freq6 #= 11) or
	 (Row6 #= 18 and Freq6 #= 13) or
	 (Row6 #= 19 and Freq6 #= 15) or
	 (Row6 #= 20 and Freq6 #= 17) or
	 (Row6 #= 21 and Freq6 #= 19) or
	 (Row6 #= 22 and Freq6 #= 21) or
	 (Row6 #= 23 and Freq6 #= 23) or
	 (Row6 #= 24 and Freq6 #= 25) or
	 (Row6 #= 25 and Freq6 #= 2) or
	 (Row6 #= 0 and Freq6 #= 0)), 

	((Row7 #= 8 and Freq7 #= 22) or
	 (Row7 #= 9 and Freq7 #= 24) or
	 (Row7 #= 10 and Freq7 #= 1) or
	 (Row7 #= 11 and Freq7 #= 2) or
	 (Row7 #= 12 and Freq7 #= 3) or
	 (Row7 #= 13 and Freq7 #= 4) or
	 (Row7 #= 14 and Freq7 #= 5) or
	 (Row7 #= 15 and Freq7 #= 7) or
	 (Row7 #= 16 and Freq7 #= 9) or
	 (Row7 #= 17 and Freq7 #= 11) or
	 (Row7 #= 18 and Freq7 #= 13) or
	 (Row7 #= 19 and Freq7 #= 15) or
	 (Row7 #= 20 and Freq7 #= 17) or
	 (Row7 #= 21 and Freq7 #= 19) or
	 (Row7 #= 22 and Freq7 #= 21) or
	 (Row7 #= 23 and Freq7 #= 23) or
	 (Row7 #= 24 and Freq7 #= 25) or
	 (Row7 #= 25 and Freq7 #= 2) or
	 (Row7 #= 0 and Freq7 #= 0)), 

	((Row8 #= 9 and Freq8 #= 25) or
	 (Row8 #= 10 and Freq8 #= 2) or
	 (Row8 #= 11 and Freq8 #= 3) or
	 (Row8 #= 12 and Freq8 #= 4) or
	 (Row8 #= 13 and Freq8 #= 5) or
	 (Row8 #= 14 and Freq8 #= 7) or
	 (Row8 #= 15 and Freq8 #= 9) or
	 (Row8 #= 16 and Freq8 #= 11) or
	 (Row8 #= 17 and Freq8 #= 13) or
	 (Row8 #= 18 and Freq8 #= 15) or
	 (Row8 #= 19 and Freq8 #= 17) or
	 (Row8 #= 20 and Freq8 #= 19) or
	 (Row8 #= 21 and Freq8 #= 21) or
	 (Row8 #= 22 and Freq8 #= 23) or
	 (Row8 #= 23 and Freq8 #= 25) or
	 (Row8 #= 24 and Freq8 #= 2) or
	 (Row8 #= 25 and Freq8 #= 3) or
	 (Row8 #= 0 and Freq8 #= 0)), 

	((Row9 #= 10 and Freq9 #= 2) or
	 (Row9 #= 11 and Freq9 #= 3) or
	 (Row9 #= 12 and Freq9 #= 4) or
	 (Row9 #= 13 and Freq9 #= 5) or
	 (Row9 #= 14 and Freq9 #= 7) or
	 (Row9 #= 15 and Freq9 #= 9) or
	 (Row9 #= 16 and Freq9 #= 11) or
	 (Row9 #= 17 and Freq9 #= 13) or
	 (Row9 #= 18 and Freq9 #= 15) or
	 (Row9 #= 19 and Freq9 #= 17) or
	 (Row9 #= 20 and Freq9 #= 19) or
	 (Row9 #= 21 and Freq9 #= 21) or
	 (Row9 #= 22 and Freq9 #= 23) or
	 (Row9 #= 23 and Freq9 #= 25) or
	 (Row9 #= 24 and Freq9 #= 2) or
	 (Row9 #= 25 and Freq9 #= 3) or
	 (Row9 #= 0 and Freq9 #= 0)), 

	((Row10 #= 11 and Freq10 #= 4) or
	 (Row10 #= 12 and Freq10 #= 5) or
	 (Row10 #= 13 and Freq10 #= 7) or
	 (Row10 #= 14 and Freq10 #= 9) or
	 (Row10 #= 15 and Freq10 #= 11) or
	 (Row10 #= 16 and Freq10 #= 13) or
	 (Row10 #= 17 and Freq10 #= 15) or
	 (Row10 #= 18 and Freq10 #= 17) or
	 (Row10 #= 19 and Freq10 #= 19) or
	 (Row10 #= 20 and Freq10 #= 21) or
	 (Row10 #= 21 and Freq10 #= 23) or
	 (Row10 #= 22 and Freq10 #= 25) or
	 (Row10 #= 23 and Freq10 #= 2) or
	 (Row10 #= 24 and Freq10 #= 3) or
	 (Row10 #= 25 and Freq10 #= 4) or
	 (Row10 #= 0 and Freq10 #= 0)), 

	((Row11 #= 12 and Freq11 #= 5) or
	 (Row11 #= 13 and Freq11 #= 7) or
	 (Row11 #= 14 and Freq11 #= 9) or
	 (Row11 #= 15 and Freq11 #= 11) or
	 (Row11 #= 16 and Freq11 #= 13) or
	 (Row11 #= 17 and Freq11 #= 15) or
	 (Row11 #= 18 and Freq11 #= 17) or
	 (Row11 #= 19 and Freq11 #= 19) or
	 (Row11 #= 20 and Freq11 #= 21) or
	 (Row11 #= 21 and Freq11 #= 23) or
	 (Row11 #= 22 and Freq11 #= 25) or
	 (Row11 #= 23 and Freq11 #= 2) or
	 (Row11 #= 24 and Freq11 #= 3) or
	 (Row11 #= 25 and Freq11 #= 4) or
	 (Row11 #= 0 and Freq11 #= 0)), 

	((Row12 #= 13 and Freq12 #= 9) or
	 (Row12 #= 14 and Freq12 #= 11) or
	 (Row12 #= 15 and Freq12 #= 13) or
	 (Row12 #= 16 and Freq12 #= 15) or
	 (Row12 #= 17 and Freq12 #= 17) or
	 (Row12 #= 18 and Freq12 #= 19) or
	 (Row12 #= 19 and Freq12 #= 21) or
	 (Row12 #= 20 and Freq12 #= 23) or
	 (Row12 #= 21 and Freq12 #= 25) or
	 (Row12 #= 22 and Freq12 #= 2) or
	 (Row12 #= 23 and Freq12 #= 3) or
	 (Row12 #= 24 and Freq12 #= 4) or
	 (Row12 #= 25 and Freq12 #= 5) or
	 (Row12 #= 0 and Freq12 #= 0)), 

	((Row13 #= 14 and Freq13 #= 11) or
	 (Row13 #= 15 and Freq13 #= 13) or
	 (Row13 #= 16 and Freq13 #= 15) or
	 (Row13 #= 17 and Freq13 #= 17) or
	 (Row13 #= 18 and Freq13 #= 19) or
	 (Row13 #= 19 and Freq13 #= 21) or
	 (Row13 #= 20 and Freq13 #= 23) or
	 (Row13 #= 21 and Freq13 #= 25) or
	 (Row13 #= 22 and Freq13 #= 2) or
	 (Row13 #= 23 and Freq13 #= 3) or
	 (Row13 #= 24 and Freq13 #= 4) or
	 (Row13 #= 25 and Freq13 #= 5) or
	 (Row13 #= 0 and Freq13 #= 0)), 

	((Row14 #= 15 and Freq14 #= 15) or
	 (Row14 #= 16 and Freq14 #= 17) or
	 (Row14 #= 17 and Freq14 #= 19) or
	 (Row14 #= 18 and Freq14 #= 21) or
	 (Row14 #= 19 and Freq14 #= 23) or
	 (Row14 #= 20 and Freq14 #= 25) or
	 (Row14 #= 21 and Freq14 #= 2) or
	 (Row14 #= 22 and Freq14 #= 3) or
	 (Row14 #= 23 and Freq14 #= 4) or
	 (Row14 #= 24 and Freq14 #= 5) or
	 (Row14 #= 25 and Freq14 #= 7) or
	 (Row14 #= 0 and Freq14 #= 0)), 

	((Row15 #= 16 and Freq15 #= 17) or
	 (Row15 #= 17 and Freq15 #= 19) or
	 (Row15 #= 18 and Freq15 #= 21) or
	 (Row15 #= 19 and Freq15 #= 23) or
	 (Row15 #= 20 and Freq15 #= 25) or
	 (Row15 #= 21 and Freq15 #= 2) or
	 (Row15 #= 22 and Freq15 #= 3) or
	 (Row15 #= 23 and Freq15 #= 4) or
	 (Row15 #= 24 and Freq15 #= 5) or
	 (Row15 #= 25 and Freq15 #= 7) or
	 (Row15 #= 0 and Freq15 #= 0)), 

	((Row16 #= 17 and Freq16 #= 21) or
	 (Row16 #= 18 and Freq16 #= 23) or
	 (Row16 #= 19 and Freq16 #= 25) or
	 (Row16 #= 20 and Freq16 #= 2) or
	 (Row16 #= 21 and Freq16 #= 3) or
	 (Row16 #= 22 and Freq16 #= 4) or
	 (Row16 #= 23 and Freq16 #= 5) or
	 (Row16 #= 24 and Freq16 #= 7) or
	 (Row16 #= 25 and Freq16 #= 9) or
	 (Row16 #= 0 and Freq16 #= 0)), 

	((Row17 #= 18 and Freq17 #= 23) or
	 (Row17 #= 19 and Freq17 #= 25) or
	 (Row17 #= 20 and Freq17 #= 2) or
	 (Row17 #= 21 and Freq17 #= 3) or
	 (Row17 #= 22 and Freq17 #= 4) or
	 (Row17 #= 23 and Freq17 #= 5) or
	 (Row17 #= 24 and Freq17 #= 7) or
	 (Row17 #= 25 and Freq17 #= 9) or
	 (Row17 #= 0 and Freq17 #= 0)), 

	((Row18 #= 19 and Freq18 #= 2) or
	 (Row18 #= 20 and Freq18 #= 3) or
	 (Row18 #= 21 and Freq18 #= 4) or
	 (Row18 #= 22 and Freq18 #= 5) or
	 (Row18 #= 23 and Freq18 #= 7) or
	 (Row18 #= 24 and Freq18 #= 9) or
	 (Row18 #= 25 and Freq18 #= 11) or
	 (Row18 #= 0 and Freq18 #= 0)), 

	((Row19 #= 20 and Freq19 #= 3) or
	 (Row19 #= 21 and Freq19 #= 4) or
	 (Row19 #= 22 and Freq19 #= 5) or
	 (Row19 #= 23 and Freq19 #= 7) or
	 (Row19 #= 24 and Freq19 #= 9) or
	 (Row19 #= 25 and Freq19 #= 11) or
	 (Row19 #= 0 and Freq19 #= 0)), 

	((Row20 #= 21 and Freq20 #= 5) or
	 (Row20 #= 22 and Freq20 #= 7) or
	 (Row20 #= 23 and Freq20 #= 9) or
	 (Row20 #= 24 and Freq20 #= 11) or
	 (Row20 #= 25 and Freq20 #= 13) or
	 (Row20 #= 0 and Freq20 #= 0)), 

	((Row21 #= 22 and Freq21 #= 7) or
	 (Row21 #= 23 and Freq21 #= 9) or
	 (Row21 #= 24 and Freq21 #= 11) or
	 (Row21 #= 25 and Freq21 #= 13) or
	 (Row21 #= 0 and Freq21 #= 0)), 

	((Row22 #= 23 and Freq22 #= 11) or
	 (Row22 #= 24 and Freq22 #= 13) or
	 (Row22 #= 25 and Freq22 #= 15) or
	 (Row22 #= 0 and Freq22 #= 0)), 

	((Row23 #= 24 and Freq23 #= 13) or
	 (Row23 #= 25 and Freq23 #= 15) or
	 (Row23 #= 0 and Freq23 #= 0)), 

	((Row24 #= 25 and Freq24 #= 17) or
	 (Row24 #= 0 and Freq24 #= 0)), 

	((Row25 #= 0 and Freq25 #= 0)), 

	% All of the values assumed by the Row<i> variables must be  
	% all different or zero; multiple zeros are allowed  
	alldifferent_except(Rows), 

	% Optimize: the sum of the selected interaction 
	% frequencies is maximal (it is the minimum of 
	% the additive inverse of the sum). The predicate 
	% enforce_symmetry/1 ensures each genomic bin is 
	% involved in only one interaction in the solution set.
	Cost #= -sum(Freqs), 
	minimize(enforce_symmetry(Rows), Cost),

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% 	Output the results
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	open(FreqFile, 'write', FREQ_OUT),
	%%list the frequencies
	(foreach(X,Freqs),
		param(FREQ_OUT) do
			get_domain_as_list(X, DomList),
				(foreach(Y,DomList),
					param(FREQ_OUT) do
						write(FREQ_OUT, Y),
						write(FREQ_OUT, ' ')
				),
			write(FREQ_OUT, "\n")
	),
	close(FREQ_OUT),

	%% list the potential rows
	open(RowFile, 'write', ROW_OUT),
	(foreach(X,Rows),
		param(ROW_OUT) do
			get_domain_as_list(X, DomList),
				(foreach(Y,DomList),
					param(ROW_OUT) do
						write(ROW_OUT, Y),
						write(ROW_OUT, ' ')
				),
				write(ROW_OUT, "\n")
	),
	close(ROW_OUT).
