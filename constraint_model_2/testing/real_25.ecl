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
	Row1 :: [0],
	Row2 :: [0],
	Row3 :: [0],
	Row4 :: [0],
	Row5 :: [0],
	Row6 :: [0, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row7 :: [0, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row8 :: [0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row9 :: [0, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row10 :: [0, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row11 :: [0, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row12 :: [0, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row13 :: [0, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row14 :: [0, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row15 :: [0, 17, 18, 19, 20, 21, 22, 23, 24, 25],
	Row16 :: [0, 18, 19, 20, 21, 22, 23, 24, 25],
	Row17 :: [0, 19, 20, 21, 22, 23, 24, 25],
	Row18 :: [0, 20, 21, 22, 23, 24, 25],
	Row19 :: [0, 21, 22, 23, 24, 25],
	Row20 :: [0, 22, 23, 24, 25],
	Row21 :: [0, 23, 24, 25],
	Row22 :: [0, 24, 25],
	Row23 :: [0, 25],
	Row24 :: [0],
	Row25 :: [0],
	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of `0' where `0'  represents an 
	% interaction not being selected 
	Freq1 :: [0],
	Freq2 :: [0],
	Freq3 :: [0],
	Freq4 :: [0],
	Freq5 :: [0],
	Freq6 :: [0, 47, 41, 33, 21, 19, 17, 14, 13, 11, 8, 7, 6, 5, 4, 3],
	Freq7 :: [0, 58, 42, 23, 20, 21, 18, 16, 14, 11, 8, 7, 6, 5, 4],
	Freq8 :: [0, 53, 23, 22, 24, 20, 17, 13, 11, 7, 6, 5, 4, 3],
	Freq9 :: [0, 32, 28, 30, 24, 21, 17, 13, 9, 7, 6, 5, 4],
	Freq10 :: [0, 41, 37, 31, 28, 22, 15, 11, 9, 6, 7, 5, 4],
	Freq11 :: [0, 40, 31, 28, 22, 16, 14, 12, 10, 9, 8, 7, 6],
	Freq12 :: [0, 41, 33, 26, 21, 16, 13, 11, 9, 8, 7, 6],
	Freq13 :: [0, 43, 32, 25, 19, 15, 11, 10, 8, 7, 6],
	Freq14 :: [0, 41, 31, 23, 19, 14, 13, 12, 9, 7],
	Freq15 :: [0, 43, 32, 23, 17, 15, 14, 12, 10, 9],
	Freq16 :: [0, 45, 33, 24, 20, 19, 15, 13, 10],
	Freq17 :: [0, 45, 30, 26, 25, 20, 14, 13],
	Freq18 :: [0, 39, 32, 29, 23, 17, 15],
	Freq19 :: [0, 44, 36, 28, 19, 16],
	Freq20 :: [0, 46, 33, 23, 19],
	Freq21 :: [0, 41, 26, 23],
	Freq22 :: [0, 35, 31],
	Freq23 :: [0, 42],
	Freq24 :: [0],
	Freq25 :: [0],

	% Constraints: 
	% Each pair of corresponding (Row<i>, Freq<i>) variables 
	% must assume dependent values based on data from the 
	% whole-genome contact map; A (Row, Freq) pair ground to 
	% (0,0) encodes that nothing is chosen 

	((Row1 #= 0 and Freq1 #= 0)), 

	((Row2 #= 0 and Freq2 #= 0)), 

	((Row3 #= 0 and Freq3 #= 0)), 

	((Row4 #= 0 and Freq4 #= 0)), 

	((Row5 #= 0 and Freq5 #= 0)), 

	((Row6 #= 8 and Freq6 #= 47) or
	 (Row6 #= 9 and Freq6 #= 41) or
	 (Row6 #= 10 and Freq6 #= 33) or
	 (Row6 #= 11 and Freq6 #= 21) or
	 (Row6 #= 12 and Freq6 #= 19) or
	 (Row6 #= 13 and Freq6 #= 19) or
	 (Row6 #= 14 and Freq6 #= 17) or
	 (Row6 #= 15 and Freq6 #= 14) or
	 (Row6 #= 16 and Freq6 #= 13) or
	 (Row6 #= 17 and Freq6 #= 11) or
	 (Row6 #= 18 and Freq6 #= 8) or
	 (Row6 #= 19 and Freq6 #= 7) or
	 (Row6 #= 20 and Freq6 #= 6) or
	 (Row6 #= 21 and Freq6 #= 5) or
	 (Row6 #= 22 and Freq6 #= 5) or
	 (Row6 #= 23 and Freq6 #= 4) or
	 (Row6 #= 24 and Freq6 #= 4) or
	 (Row6 #= 25 and Freq6 #= 3) or
	 (Row6 #= 0 and Freq6 #= 0)), 

	((Row7 #= 9 and Freq7 #= 58) or
	 (Row7 #= 10 and Freq7 #= 42) or
	 (Row7 #= 11 and Freq7 #= 23) or
	 (Row7 #= 12 and Freq7 #= 20) or
	 (Row7 #= 13 and Freq7 #= 21) or
	 (Row7 #= 14 and Freq7 #= 18) or
	 (Row7 #= 15 and Freq7 #= 16) or
	 (Row7 #= 16 and Freq7 #= 14) or
	 (Row7 #= 17 and Freq7 #= 11) or
	 (Row7 #= 18 and Freq7 #= 8) or
	 (Row7 #= 19 and Freq7 #= 7) or
	 (Row7 #= 20 and Freq7 #= 6) or
	 (Row7 #= 21 and Freq7 #= 5) or
	 (Row7 #= 22 and Freq7 #= 5) or
	 (Row7 #= 23 and Freq7 #= 4) or
	 (Row7 #= 24 and Freq7 #= 4) or
	 (Row7 #= 25 and Freq7 #= 4) or
	 (Row7 #= 0 and Freq7 #= 0)), 

	((Row8 #= 10 and Freq8 #= 53) or
	 (Row8 #= 11 and Freq8 #= 23) or
	 (Row8 #= 12 and Freq8 #= 22) or
	 (Row8 #= 13 and Freq8 #= 24) or
	 (Row8 #= 14 and Freq8 #= 20) or
	 (Row8 #= 15 and Freq8 #= 17) or
	 (Row8 #= 16 and Freq8 #= 13) or
	 (Row8 #= 17 and Freq8 #= 11) or
	 (Row8 #= 18 and Freq8 #= 7) or
	 (Row8 #= 19 and Freq8 #= 6) or
	 (Row8 #= 20 and Freq8 #= 5) or
	 (Row8 #= 21 and Freq8 #= 5) or
	 (Row8 #= 22 and Freq8 #= 4) or
	 (Row8 #= 23 and Freq8 #= 4) or
	 (Row8 #= 24 and Freq8 #= 3) or
	 (Row8 #= 25 and Freq8 #= 3) or
	 (Row8 #= 0 and Freq8 #= 0)), 

	((Row9 #= 11 and Freq9 #= 32) or
	 (Row9 #= 12 and Freq9 #= 28) or
	 (Row9 #= 13 and Freq9 #= 30) or
	 (Row9 #= 14 and Freq9 #= 24) or
	 (Row9 #= 15 and Freq9 #= 21) or
	 (Row9 #= 16 and Freq9 #= 17) or
	 (Row9 #= 17 and Freq9 #= 13) or
	 (Row9 #= 18 and Freq9 #= 9) or
	 (Row9 #= 19 and Freq9 #= 7) or
	 (Row9 #= 20 and Freq9 #= 6) or
	 (Row9 #= 21 and Freq9 #= 5) or
	 (Row9 #= 22 and Freq9 #= 4) or
	 (Row9 #= 23 and Freq9 #= 4) or
	 (Row9 #= 24 and Freq9 #= 4) or
	 (Row9 #= 25 and Freq9 #= 4) or
	 (Row9 #= 0 and Freq9 #= 0)), 

	((Row10 #= 12 and Freq10 #= 41) or
	 (Row10 #= 13 and Freq10 #= 37) or
	 (Row10 #= 14 and Freq10 #= 31) or
	 (Row10 #= 15 and Freq10 #= 28) or
	 (Row10 #= 16 and Freq10 #= 22) or
	 (Row10 #= 17 and Freq10 #= 15) or
	 (Row10 #= 18 and Freq10 #= 11) or
	 (Row10 #= 19 and Freq10 #= 9) or
	 (Row10 #= 20 and Freq10 #= 6) or
	 (Row10 #= 21 and Freq10 #= 7) or
	 (Row10 #= 22 and Freq10 #= 6) or
	 (Row10 #= 23 and Freq10 #= 5) or
	 (Row10 #= 24 and Freq10 #= 5) or
	 (Row10 #= 25 and Freq10 #= 4) or
	 (Row10 #= 0 and Freq10 #= 0)), 

	((Row11 #= 13 and Freq11 #= 40) or
	 (Row11 #= 14 and Freq11 #= 31) or
	 (Row11 #= 15 and Freq11 #= 28) or
	 (Row11 #= 16 and Freq11 #= 22) or
	 (Row11 #= 17 and Freq11 #= 16) or
	 (Row11 #= 18 and Freq11 #= 14) or
	 (Row11 #= 19 and Freq11 #= 12) or
	 (Row11 #= 20 and Freq11 #= 10) or
	 (Row11 #= 21 and Freq11 #= 9) or
	 (Row11 #= 22 and Freq11 #= 8) or
	 (Row11 #= 23 and Freq11 #= 7) or
	 (Row11 #= 24 and Freq11 #= 7) or
	 (Row11 #= 25 and Freq11 #= 6) or
	 (Row11 #= 0 and Freq11 #= 0)), 

	((Row12 #= 14 and Freq12 #= 41) or
	 (Row12 #= 15 and Freq12 #= 33) or
	 (Row12 #= 16 and Freq12 #= 26) or
	 (Row12 #= 17 and Freq12 #= 21) or
	 (Row12 #= 18 and Freq12 #= 16) or
	 (Row12 #= 19 and Freq12 #= 13) or
	 (Row12 #= 20 and Freq12 #= 11) or
	 (Row12 #= 21 and Freq12 #= 9) or
	 (Row12 #= 22 and Freq12 #= 9) or
	 (Row12 #= 23 and Freq12 #= 8) or
	 (Row12 #= 24 and Freq12 #= 7) or
	 (Row12 #= 25 and Freq12 #= 6) or
	 (Row12 #= 0 and Freq12 #= 0)), 

	((Row13 #= 15 and Freq13 #= 43) or
	 (Row13 #= 16 and Freq13 #= 32) or
	 (Row13 #= 17 and Freq13 #= 25) or
	 (Row13 #= 18 and Freq13 #= 19) or
	 (Row13 #= 19 and Freq13 #= 15) or
	 (Row13 #= 20 and Freq13 #= 11) or
	 (Row13 #= 21 and Freq13 #= 10) or
	 (Row13 #= 22 and Freq13 #= 10) or
	 (Row13 #= 23 and Freq13 #= 8) or
	 (Row13 #= 24 and Freq13 #= 7) or
	 (Row13 #= 25 and Freq13 #= 6) or
	 (Row13 #= 0 and Freq13 #= 0)), 

	((Row14 #= 16 and Freq14 #= 41) or
	 (Row14 #= 17 and Freq14 #= 31) or
	 (Row14 #= 18 and Freq14 #= 23) or
	 (Row14 #= 19 and Freq14 #= 19) or
	 (Row14 #= 20 and Freq14 #= 14) or
	 (Row14 #= 21 and Freq14 #= 13) or
	 (Row14 #= 22 and Freq14 #= 12) or
	 (Row14 #= 23 and Freq14 #= 9) or
	 (Row14 #= 24 and Freq14 #= 9) or
	 (Row14 #= 25 and Freq14 #= 7) or
	 (Row14 #= 0 and Freq14 #= 0)), 

	((Row15 #= 17 and Freq15 #= 43) or
	 (Row15 #= 18 and Freq15 #= 32) or
	 (Row15 #= 19 and Freq15 #= 23) or
	 (Row15 #= 20 and Freq15 #= 17) or
	 (Row15 #= 21 and Freq15 #= 15) or
	 (Row15 #= 22 and Freq15 #= 14) or
	 (Row15 #= 23 and Freq15 #= 12) or
	 (Row15 #= 24 and Freq15 #= 10) or
	 (Row15 #= 25 and Freq15 #= 9) or
	 (Row15 #= 0 and Freq15 #= 0)), 

	((Row16 #= 18 and Freq16 #= 45) or
	 (Row16 #= 19 and Freq16 #= 33) or
	 (Row16 #= 20 and Freq16 #= 24) or
	 (Row16 #= 21 and Freq16 #= 20) or
	 (Row16 #= 22 and Freq16 #= 19) or
	 (Row16 #= 23 and Freq16 #= 15) or
	 (Row16 #= 24 and Freq16 #= 13) or
	 (Row16 #= 25 and Freq16 #= 10) or
	 (Row16 #= 0 and Freq16 #= 0)), 

	((Row17 #= 19 and Freq17 #= 45) or
	 (Row17 #= 20 and Freq17 #= 30) or
	 (Row17 #= 21 and Freq17 #= 26) or
	 (Row17 #= 22 and Freq17 #= 25) or
	 (Row17 #= 23 and Freq17 #= 20) or
	 (Row17 #= 24 and Freq17 #= 14) or
	 (Row17 #= 25 and Freq17 #= 13) or
	 (Row17 #= 0 and Freq17 #= 0)), 

	((Row18 #= 20 and Freq18 #= 39) or
	 (Row18 #= 21 and Freq18 #= 32) or
	 (Row18 #= 22 and Freq18 #= 29) or
	 (Row18 #= 23 and Freq18 #= 23) or
	 (Row18 #= 24 and Freq18 #= 17) or
	 (Row18 #= 25 and Freq18 #= 15) or
	 (Row18 #= 0 and Freq18 #= 0)), 

	((Row19 #= 21 and Freq19 #= 44) or
	 (Row19 #= 22 and Freq19 #= 36) or
	 (Row19 #= 23 and Freq19 #= 28) or
	 (Row19 #= 24 and Freq19 #= 19) or
	 (Row19 #= 25 and Freq19 #= 16) or
	 (Row19 #= 0 and Freq19 #= 0)), 

	((Row20 #= 22 and Freq20 #= 46) or
	 (Row20 #= 23 and Freq20 #= 33) or
	 (Row20 #= 24 and Freq20 #= 23) or
	 (Row20 #= 25 and Freq20 #= 19) or
	 (Row20 #= 0 and Freq20 #= 0)), 

	((Row21 #= 23 and Freq21 #= 41) or
	 (Row21 #= 24 and Freq21 #= 26) or
	 (Row21 #= 25 and Freq21 #= 23) or
	 (Row21 #= 0 and Freq21 #= 0)), 

	((Row22 #= 24 and Freq22 #= 35) or
	 (Row22 #= 25 and Freq22 #= 31) or
	 (Row22 #= 0 and Freq22 #= 0)), 

	((Row23 #= 25 and Freq23 #= 42) or
	 (Row23 #= 0 and Freq23 #= 0)), 

	((Row24 #= 0 and Freq24 #= 0)), 

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
