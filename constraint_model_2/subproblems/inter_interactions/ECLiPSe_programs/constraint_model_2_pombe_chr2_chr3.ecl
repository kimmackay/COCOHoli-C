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
	Rows = [Row1, Row2, Row3, Row4, Row5, Row6, Row7, Row8, Row9, Row10, Row11, Row12, Row13, Row14, Row15, Row16, Row17, Row18, Row19, Row20, Row21, Row22, Row23, Row24, Row25, Row26, Row27, Row28, Row29, Row30, Row31, Row32, Row33, Row34, Row35, Row36, Row37, Row38, Row39, Row40, Row41, Row42, Row43, Row44, Row45, Row46, Row47, Row48, Row49, Row50, Row51, Row52, Row53, Row54, Row55, Row56, Row57, Row58, Row59, Row60, Row61, Row62, Row63, Row64, Row65, Row66, Row67, Row68, Row69, Row70, Row71, Row72, Row73, Row74, Row75, Row76, Row77, Row78, Row79, Row80, Row81, Row82, Row83, Row84, Row85, Row86, Row87, Row88, Row89, Row90, Row91, Row92, Row93, Row94, Row95, Row96, Row97, Row98, Row99, Row100, Row101, Row102, Row103, Row104, Row105, Row106, Row107, Row108, Row109, Row110, Row111, Row112, Row113, Row114, Row115, Row116, Row117, Row118, Row119, Row120, Row121, Row122, Row123, Row124, Row125, Row126, Row127, Row128, Row129, Row130, Row131, Row132, Row133, Row134, Row135, Row136, Row137, Row138, Row139, Row140, Row141, Row142, Row143, Row144, Row145, Row146, Row147, Row148, Row149, Row150, Row151, Row152, Row153, Row154, Row155, Row156, Row157, Row158, Row159, Row160, Row161, Row162, Row163, Row164, Row165, Row166, Row167, Row168, Row169, Row170, Row171, Row172, Row173, Row174, Row175, Row176, Row177, Row178, Row179, Row180, Row181, Row182, Row183, Row184, Row185, Row186, Row187, Row188, Row189, Row190, Row191, Row192, Row193, Row194, Row195, Row196, Row197, Row198, Row199, Row200, Row201, Row202, Row203, Row204, Row205, Row206, Row207, Row208, Row209, Row210, Row211, Row212, Row213, Row214, Row215, Row216, Row217, Row218, Row219, Row220, Row221, Row222, Row223, Row224, Row225, Row226, Row227, Row228, Row229, Row230, Row231, Row232, Row233, Row234, Row235, Row236, Row237, Row238, Row239, Row240, Row241, Row242, Row243, Row244, Row245, Row246, Row247, Row248, Row249, Row250, Row251, Row252, Row253, Row254, Row255, Row256, Row257, Row258, Row259, Row260, Row261, Row262, Row263, Row264, Row265, Row266, Row267, Row268, Row269, Row270, Row271, Row272, Row273, Row274, Row275, Row276, Row277, Row278, Row279, Row280, Row281, Row282, Row283, Row284, Row285, Row286, Row287, Row288, Row289, Row290, Row291, Row292, Row293, Row294, Row295, Row296, Row297, Row298, Row299, Row300, Row301, Row302, Row303, Row304, Row305, Row306, Row307, Row308, Row309, Row310, Row311, Row312, Row313, Row314, Row315, Row316, Row317, Row318, Row319, Row320, Row321, Row322, Row323, Row324, Row325, Row326, Row327, Row328, Row329, Row330, Row331, Row332, Row333, Row334, Row335, Row336, Row337, Row338, Row339, Row340, Row341, Row342, Row343, Row344, Row345, Row346, Row347, Row348, Row349, Row350, Row351, Row352, Row353, Row354, Row355, Row356, Row357, Row358, Row359, Row360, Row361, Row362, Row363, Row364, Row365, Row366, Row367, Row368, Row369, Row370, Row371, Row372, Row373, Row374, Row375, Row376, Row377, Row378, Row379, Row380, Row381, Row382, Row383, Row384, Row385, Row386, Row387, Row388, Row389, Row390, Row391, Row392, Row393, Row394, Row395, Row396, Row397, Row398, Row399, Row400, Row401, Row402, Row403, Row404, Row405, Row406, Row407, Row408, Row409, Row410, Row411, Row412, Row413, Row414, Row415, Row416, Row417, Row418, Row419, Row420, Row421, Row422, Row423, Row424, Row425, Row426, Row427, Row428, Row429, Row430, Row431, Row432, Row433, Row434, Row435, Row436, Row437, Row438, Row439, Row440, Row441, Row442, Row443, Row444, Row445, Row446, Row447, Row448, Row449, Row450, Row451, Row452, Row453, Row454],

	% The list Freqs has one variable for each row of the 
	% whole-genome contact map	
	Freqs = [Freq1, Freq2, Freq3, Freq4, Freq5, Freq6, Freq7, Freq8, Freq9, Freq10, Freq11, Freq12, Freq13, Freq14, Freq15, Freq16, Freq17, Freq18, Freq19, Freq20, Freq21, Freq22, Freq23, Freq24, Freq25, Freq26, Freq27, Freq28, Freq29, Freq30, Freq31, Freq32, Freq33, Freq34, Freq35, Freq36, Freq37, Freq38, Freq39, Freq40, Freq41, Freq42, Freq43, Freq44, Freq45, Freq46, Freq47, Freq48, Freq49, Freq50, Freq51, Freq52, Freq53, Freq54, Freq55, Freq56, Freq57, Freq58, Freq59, Freq60, Freq61, Freq62, Freq63, Freq64, Freq65, Freq66, Freq67, Freq68, Freq69, Freq70, Freq71, Freq72, Freq73, Freq74, Freq75, Freq76, Freq77, Freq78, Freq79, Freq80, Freq81, Freq82, Freq83, Freq84, Freq85, Freq86, Freq87, Freq88, Freq89, Freq90, Freq91, Freq92, Freq93, Freq94, Freq95, Freq96, Freq97, Freq98, Freq99, Freq100, Freq101, Freq102, Freq103, Freq104, Freq105, Freq106, Freq107, Freq108, Freq109, Freq110, Freq111, Freq112, Freq113, Freq114, Freq115, Freq116, Freq117, Freq118, Freq119, Freq120, Freq121, Freq122, Freq123, Freq124, Freq125, Freq126, Freq127, Freq128, Freq129, Freq130, Freq131, Freq132, Freq133, Freq134, Freq135, Freq136, Freq137, Freq138, Freq139, Freq140, Freq141, Freq142, Freq143, Freq144, Freq145, Freq146, Freq147, Freq148, Freq149, Freq150, Freq151, Freq152, Freq153, Freq154, Freq155, Freq156, Freq157, Freq158, Freq159, Freq160, Freq161, Freq162, Freq163, Freq164, Freq165, Freq166, Freq167, Freq168, Freq169, Freq170, Freq171, Freq172, Freq173, Freq174, Freq175, Freq176, Freq177, Freq178, Freq179, Freq180, Freq181, Freq182, Freq183, Freq184, Freq185, Freq186, Freq187, Freq188, Freq189, Freq190, Freq191, Freq192, Freq193, Freq194, Freq195, Freq196, Freq197, Freq198, Freq199, Freq200, Freq201, Freq202, Freq203, Freq204, Freq205, Freq206, Freq207, Freq208, Freq209, Freq210, Freq211, Freq212, Freq213, Freq214, Freq215, Freq216, Freq217, Freq218, Freq219, Freq220, Freq221, Freq222, Freq223, Freq224, Freq225, Freq226, Freq227, Freq228, Freq229, Freq230, Freq231, Freq232, Freq233, Freq234, Freq235, Freq236, Freq237, Freq238, Freq239, Freq240, Freq241, Freq242, Freq243, Freq244, Freq245, Freq246, Freq247, Freq248, Freq249, Freq250, Freq251, Freq252, Freq253, Freq254, Freq255, Freq256, Freq257, Freq258, Freq259, Freq260, Freq261, Freq262, Freq263, Freq264, Freq265, Freq266, Freq267, Freq268, Freq269, Freq270, Freq271, Freq272, Freq273, Freq274, Freq275, Freq276, Freq277, Freq278, Freq279, Freq280, Freq281, Freq282, Freq283, Freq284, Freq285, Freq286, Freq287, Freq288, Freq289, Freq290, Freq291, Freq292, Freq293, Freq294, Freq295, Freq296, Freq297, Freq298, Freq299, Freq300, Freq301, Freq302, Freq303, Freq304, Freq305, Freq306, Freq307, Freq308, Freq309, Freq310, Freq311, Freq312, Freq313, Freq314, Freq315, Freq316, Freq317, Freq318, Freq319, Freq320, Freq321, Freq322, Freq323, Freq324, Freq325, Freq326, Freq327, Freq328, Freq329, Freq330, Freq331, Freq332, Freq333, Freq334, Freq335, Freq336, Freq337, Freq338, Freq339, Freq340, Freq341, Freq342, Freq343, Freq344, Freq345, Freq346, Freq347, Freq348, Freq349, Freq350, Freq351, Freq352, Freq353, Freq354, Freq355, Freq356, Freq357, Freq358, Freq359, Freq360, Freq361, Freq362, Freq363, Freq364, Freq365, Freq366, Freq367, Freq368, Freq369, Freq370, Freq371, Freq372, Freq373, Freq374, Freq375, Freq376, Freq377, Freq378, Freq379, Freq380, Freq381, Freq382, Freq383, Freq384, Freq385, Freq386, Freq387, Freq388, Freq389, Freq390, Freq391, Freq392, Freq393, Freq394, Freq395, Freq396, Freq397, Freq398, Freq399, Freq400, Freq401, Freq402, Freq403, Freq404, Freq405, Freq406, Freq407, Freq408, Freq409, Freq410, Freq411, Freq412, Freq413, Freq414, Freq415, Freq416, Freq417, Freq418, Freq419, Freq420, Freq421, Freq422, Freq423, Freq424, Freq425, Freq426, Freq427, Freq428, Freq429, Freq430, Freq431, Freq432, Freq433, Freq434, Freq435, Freq436, Freq437, Freq438, Freq439, Freq440, Freq441, Freq442, Freq443, Freq444, Freq445, Freq446, Freq447, Freq448, Freq449, Freq450, Freq451, Freq452, Freq453, Freq454],

	% Representation of the Genome: 
	% Each Row term can assume a value based on interacting bin 
	% indices where `0' represents an interaction not being 
	% selected and a non-zero value (ranging from 1 to N) 
	% represents which genomic bin is involved in the selected 
	% interaction 
	Row1 :: [0],
	Row2 :: [0],
	Row3 :: [0],
	Row4 :: [0, 245],
	Row5 :: [0, 245],
	Row6 :: [0, 246],
	Row7 :: [0],
	Row8 :: [0, 246],
	Row9 :: [0],
	Row10 :: [0, 245],
	Row11 :: [0],
	Row12 :: [0],
	Row13 :: [0, 245],
	Row14 :: [0],
	Row15 :: [0],
	Row16 :: [0],
	Row17 :: [0],
	Row18 :: [0],
	Row19 :: [0],
	Row20 :: [0],
	Row21 :: [0],
	Row22 :: [0, 246],
	Row23 :: [0],
	Row24 :: [0, 245],
	Row25 :: [0],
	Row26 :: [0, 245],
	Row27 :: [0, 246],
	Row28 :: [0],
	Row29 :: [0],
	Row30 :: [0],
	Row31 :: [0],
	Row32 :: [0],
	Row33 :: [0, 245, 246],
	Row34 :: [0, 245],
	Row35 :: [0, 245],
	Row36 :: [0, 245, 246],
	Row37 :: [0],
	Row38 :: [0, 245],
	Row39 :: [0],
	Row40 :: [0],
	Row41 :: [0, 245],
	Row42 :: [0],
	Row43 :: [0, 245],
	Row44 :: [0],
	Row45 :: [0, 245],
	Row46 :: [0, 245],
	Row47 :: [0, 245],
	Row48 :: [0],
	Row49 :: [0],
	Row50 :: [0],
	Row51 :: [0],
	Row52 :: [0],
	Row53 :: [0],
	Row54 :: [0],
	Row55 :: [0],
	Row56 :: [0],
	Row57 :: [0],
	Row58 :: [0, 82],
	Row59 :: [0],
	Row60 :: [0],
	Row61 :: [0],
	Row62 :: [0],
	Row63 :: [0],
	Row64 :: [0],
	Row65 :: [0],
	Row66 :: [0],
	Row67 :: [0],
	Row68 :: [0],
	Row69 :: [0],
	Row70 :: [0],
	Row71 :: [0],
	Row72 :: [0],
	Row73 :: [0],
	Row74 :: [0],
	Row75 :: [0],
	Row76 :: [0],
	Row77 :: [0],
	Row78 :: [0],
	Row79 :: [0],
	Row80 :: [0],
	Row81 :: [0],
	Row82 :: [0],
	Row83 :: [0],
	Row84 :: [0],
	Row85 :: [0],
	Row86 :: [0],
	Row87 :: [0],
	Row88 :: [0],
	Row89 :: [0],
	Row90 :: [0],
	Row91 :: [0],
	Row92 :: [0],
	Row93 :: [0],
	Row94 :: [0],
	Row95 :: [0],
	Row96 :: [0],
	Row97 :: [0],
	Row98 :: [0],
	Row99 :: [0],
	Row100 :: [0],
	Row101 :: [0],
	Row102 :: [0],
	Row103 :: [0],
	Row104 :: [0],
	Row105 :: [0],
	Row106 :: [0],
	Row107 :: [0, 245],
	Row108 :: [0],
	Row109 :: [0],
	Row110 :: [0],
	Row111 :: [0],
	Row112 :: [0],
	Row113 :: [0],
	Row114 :: [0],
	Row115 :: [0],
	Row116 :: [0],
	Row117 :: [0],
	Row118 :: [0],
	Row119 :: [0, 152],
	Row120 :: [0, 145],
	Row121 :: [0],
	Row122 :: [0],
	Row123 :: [0, 245],
	Row124 :: [0, 144, 246],
	Row125 :: [0, 138, 143, 144, 145, 146, 148],
	Row126 :: [0, 138, 141, 142, 143, 144, 145, 147],
	Row127 :: [0, 140, 141, 142, 143, 144, 145, 146, 147, 148, 150, 153],
	Row128 :: [0, 137, 143, 144, 145, 146, 147, 148, 150],
	Row129 :: [0, 138, 143, 144, 145, 146, 147, 148],
	Row130 :: [0, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149],
	Row131 :: [0, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152],
	Row132 :: [0, 39, 130, 133, 136, 137, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 155, 245],
	Row133 :: [0, 130, 134, 135, 136, 137, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 156, 245],
	Row134 :: [0, 101, 129, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150],
	Row135 :: [0, 129, 131, 134, 135, 136, 137, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150],
	Row136 :: [0, 88, 92, 93, 94, 102, 105, 129, 131, 132, 134, 135, 136, 137, 138, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 245],
	Row137 :: [0, 97, 99, 100, 102, 105, 121, 126, 127, 129, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 154, 245],
	Row138 :: [0, 96, 97, 98, 101, 104, 105, 106, 127, 129, 131, 134, 136, 137, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 150],
	Row139 :: [0, 91, 97, 99, 100, 104, 121, 126, 128, 129, 130, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 150],
	Row140 :: [0, 87, 92, 93, 94, 96, 97, 98, 99, 100, 101, 102, 104, 123, 124, 125, 126, 127, 128, 129, 131, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148],
	Row141 :: [0, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 118, 121, 122, 123, 125, 126, 127, 128, 129, 130, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152, 153],
	Row142 :: [0, 88, 89, 90, 92, 93, 96, 97, 98, 99, 100, 101, 102, 105, 116, 117, 120, 121, 125, 126, 128, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152],
	Row143 :: [0, 87, 90, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 129, 130, 131, 132, 133, 134, 135, 136, 137, 139, 140, 141, 142, 143, 144, 145, 146, 147, 150],
	Row144 :: [0, 81, 84, 88, 89, 90, 91, 92, 93, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152],
	Row145 :: [0, 49, 77, 79, 81, 82, 84, 85, 86, 87, 90, 91, 92, 93, 94, 95, 97, 98, 99, 100, 101, 102, 104, 105, 106, 115, 116, 117, 118, 120, 121, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 149, 159],
	Row146 :: [0, 72, 80, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 150, 151, 152, 153, 155, 245],
	Row147 :: [0, 73, 80, 82, 84, 86, 87, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 153, 154],
	Row148 :: [0, 39, 69, 73, 80, 82, 83, 84, 86, 87, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 153, 154],
	Row149 :: [0, 79, 81, 83, 84, 86, 87, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153],
	Row150 :: [0, 81, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149],
	Row151 :: [0, 80, 82, 83, 84, 85, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 151, 152],
	Row152 :: [0, 81, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147],
	Row153 :: [0, 80, 82, 84, 86, 87, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 144, 145, 146, 147],
	Row154 :: [0, 80, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144],
	Row155 :: [0, 80, 82, 83, 84, 87, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 143, 146],
	Row156 :: [0, 82, 83, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 144, 146],
	Row157 :: [0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 143, 146],
	Row158 :: [0, 77, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 144],
	Row159 :: [0, 82, 83, 84, 85, 86, 87, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142],
	Row160 :: [0],
	Row161 :: [0],
	Row162 :: [0],
	Row163 :: [0],
	Row164 :: [0],
	Row165 :: [0],
	Row166 :: [0, 82, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 134, 137, 138],
	Row167 :: [0, 81, 82, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 136, 137, 138],
	Row168 :: [0, 82, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 136, 137],
	Row169 :: [0, 82, 84, 87, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 135, 136, 137, 140],
	Row170 :: [0, 90, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 138],
	Row171 :: [0, 92, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 136],
	Row172 :: [0, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 129, 130, 131, 140],
	Row173 :: [0, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 129, 130, 131],
	Row174 :: [0, 88, 89, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 124, 125, 128, 129, 130, 131],
	Row175 :: [0, 83, 88, 89, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 124, 125, 126, 130, 132],
	Row176 :: [0, 82, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 129, 130, 131, 132, 136, 140],
	Row177 :: [0, 90, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 125, 126, 128, 129, 131],
	Row178 :: [0, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 116, 117, 118, 119, 120, 121, 123, 125, 129, 132, 140],
	Row179 :: [0, 91, 92, 93, 95, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 122, 124, 125, 126],
	Row180 :: [0, 91, 92, 93, 98, 99, 100, 101, 103, 105, 106, 115, 116, 117, 120, 121, 125, 126, 130],
	Row181 :: [0, 87, 91, 92, 93, 98, 99, 100, 101, 105, 117, 120, 121, 125, 130],
	Row182 :: [0, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 117, 122, 137],
	Row183 :: [0, 85, 88, 93, 98, 100, 101],
	Row184 :: [0, 97, 104],
	Row185 :: [0, 93],
	Row186 :: [0, 96, 97, 99],
	Row187 :: [0, 98, 104, 105],
	Row188 :: [0, 97, 99],
	Row189 :: [0],
	Row190 :: [0],
	Row191 :: [0],
	Row192 :: [0, 104],
	Row193 :: [0],
	Row194 :: [0, 96],
	Row195 :: [0, 77],
	Row196 :: [0],
	Row197 :: [0, 152],
	Row198 :: [0],
	Row199 :: [0],
	Row200 :: [0],
	Row201 :: [0],
	Row202 :: [0],
	Row203 :: [0],
	Row204 :: [0],
	Row205 :: [0],
	Row206 :: [0],
	Row207 :: [0],
	Row208 :: [0],
	Row209 :: [0],
	Row210 :: [0],
	Row211 :: [0],
	Row212 :: [0],
	Row213 :: [0],
	Row214 :: [0],
	Row215 :: [0],
	Row216 :: [0],
	Row217 :: [0],
	Row218 :: [0],
	Row219 :: [0],
	Row220 :: [0],
	Row221 :: [0],
	Row222 :: [0],
	Row223 :: [0],
	Row224 :: [0],
	Row225 :: [0],
	Row226 :: [0],
	Row227 :: [0],
	Row228 :: [0],
	Row229 :: [0],
	Row230 :: [0],
	Row231 :: [0],
	Row232 :: [0],
	Row233 :: [0],
	Row234 :: [0],
	Row235 :: [0],
	Row236 :: [0],
	Row237 :: [0],
	Row238 :: [0],
	Row239 :: [0],
	Row240 :: [0],
	Row241 :: [0],
	Row242 :: [0],
	Row243 :: [0],
	Row244 :: [0],
	Row245 :: [0],
	Row246 :: [0],
	Row247 :: [0],
	Row248 :: [0],
	Row249 :: [0],
	Row250 :: [0],
	Row251 :: [0],
	Row252 :: [0],
	Row253 :: [0],
	Row254 :: [0],
	Row255 :: [0],
	Row256 :: [0],
	Row257 :: [0],
	Row258 :: [0],
	Row259 :: [0],
	Row260 :: [0],
	Row261 :: [0],
	Row262 :: [0],
	Row263 :: [0],
	Row264 :: [0],
	Row265 :: [0, 245],
	Row266 :: [0],
	Row267 :: [0],
	Row268 :: [0],
	Row269 :: [0],
	Row270 :: [0],
	Row271 :: [0],
	Row272 :: [0],
	Row273 :: [0],
	Row274 :: [0],
	Row275 :: [0],
	Row276 :: [0],
	Row277 :: [0],
	Row278 :: [0],
	Row279 :: [0],
	Row280 :: [0],
	Row281 :: [0],
	Row282 :: [0],
	Row283 :: [0],
	Row284 :: [0],
	Row285 :: [0],
	Row286 :: [0],
	Row287 :: [0],
	Row288 :: [0],
	Row289 :: [0],
	Row290 :: [0],
	Row291 :: [0],
	Row292 :: [0],
	Row293 :: [0],
	Row294 :: [0],
	Row295 :: [0],
	Row296 :: [0],
	Row297 :: [0],
	Row298 :: [0],
	Row299 :: [0],
	Row300 :: [0],
	Row301 :: [0],
	Row302 :: [0],
	Row303 :: [0],
	Row304 :: [0],
	Row305 :: [0],
	Row306 :: [0],
	Row307 :: [0],
	Row308 :: [0],
	Row309 :: [0],
	Row310 :: [0],
	Row311 :: [0],
	Row312 :: [0],
	Row313 :: [0],
	Row314 :: [0],
	Row315 :: [0],
	Row316 :: [0],
	Row317 :: [0],
	Row318 :: [0],
	Row319 :: [0, 245],
	Row320 :: [0],
	Row321 :: [0],
	Row322 :: [0],
	Row323 :: [0],
	Row324 :: [0],
	Row325 :: [0],
	Row326 :: [0],
	Row327 :: [0],
	Row328 :: [0],
	Row329 :: [0],
	Row330 :: [0],
	Row331 :: [0],
	Row332 :: [0],
	Row333 :: [0],
	Row334 :: [0],
	Row335 :: [0],
	Row336 :: [0],
	Row337 :: [0],
	Row338 :: [0],
	Row339 :: [0],
	Row340 :: [0, 245],
	Row341 :: [0],
	Row342 :: [0],
	Row343 :: [0],
	Row344 :: [0],
	Row345 :: [0],
	Row346 :: [0],
	Row347 :: [0],
	Row348 :: [0],
	Row349 :: [0],
	Row350 :: [0, 245],
	Row351 :: [0],
	Row352 :: [0],
	Row353 :: [0],
	Row354 :: [0],
	Row355 :: [0],
	Row356 :: [0],
	Row357 :: [0],
	Row358 :: [0],
	Row359 :: [0, 245],
	Row360 :: [0],
	Row361 :: [0],
	Row362 :: [0],
	Row363 :: [0],
	Row364 :: [0, 245],
	Row365 :: [0, 245],
	Row366 :: [0, 245, 246],
	Row367 :: [0, 245],
	Row368 :: [0, 245],
	Row369 :: [0, 245],
	Row370 :: [0],
	Row371 :: [0],
	Row372 :: [0],
	Row373 :: [0],
	Row374 :: [0],
	Row375 :: [0],
	Row376 :: [0],
	Row377 :: [0],
	Row378 :: [0],
	Row379 :: [0],
	Row380 :: [0],
	Row381 :: [0],
	Row382 :: [0],
	Row383 :: [0],
	Row384 :: [0],
	Row385 :: [0],
	Row386 :: [0],
	Row387 :: [0],
	Row388 :: [0],
	Row389 :: [0],
	Row390 :: [0],
	Row391 :: [0],
	Row392 :: [0],
	Row393 :: [0, 245],
	Row394 :: [0],
	Row395 :: [0],
	Row396 :: [0],
	Row397 :: [0],
	Row398 :: [0],
	Row399 :: [0],
	Row400 :: [0],
	Row401 :: [0],
	Row402 :: [0],
	Row403 :: [0],
	Row404 :: [0],
	Row405 :: [0],
	Row406 :: [0, 245],
	Row407 :: [0],
	Row408 :: [0],
	Row409 :: [0, 245],
	Row410 :: [0, 245, 246],
	Row411 :: [0, 245],
	Row412 :: [0, 246],
	Row413 :: [0, 245],
	Row414 :: [0, 246],
	Row415 :: [0],
	Row416 :: [0, 245],
	Row417 :: [0, 245],
	Row418 :: [0, 245],
	Row419 :: [0],
	Row420 :: [0, 245],
	Row421 :: [0, 245, 246],
	Row422 :: [0, 245, 246],
	Row423 :: [0, 245, 246],
	Row424 :: [0, 245, 246],
	Row425 :: [0, 245, 246],
	Row426 :: [0, 245, 246],
	Row427 :: [0, 245, 246],
	Row428 :: [0, 246],
	Row429 :: [0, 245, 246],
	Row430 :: [0, 245],
	Row431 :: [0],
	Row432 :: [0, 245],
	Row433 :: [0],
	Row434 :: [0, 245, 246],
	Row435 :: [0],
	Row436 :: [0, 245, 246],
	Row437 :: [0, 245],
	Row438 :: [0, 245],
	Row439 :: [0, 245],
	Row440 :: [0, 245, 246],
	Row441 :: [0, 245, 246],
	Row442 :: [0, 245, 246],
	Row443 :: [0],
	Row444 :: [0, 245],
	Row445 :: [0, 245, 246],
	Row446 :: [0, 245, 246],
	Row447 :: [0, 245, 246],
	Row448 :: [0, 245, 246],
	Row449 :: [0, 245, 246],
	Row450 :: [0, 245],
	Row451 :: [0],
	Row452 :: [0],
	Row453 :: [0],
	Row454 :: [0],
	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of 0 where 0  represents an 
	% interaction not being selected 
	Freq1 :: [0],
	Freq2 :: [0],
	Freq3 :: [0],
	Freq4 :: [0, 1],
	Freq5 :: [0, 1],
	Freq6 :: [0, 1],
	Freq7 :: [0],
	Freq8 :: [0, 1],
	Freq9 :: [0],
	Freq10 :: [0, 1],
	Freq11 :: [0],
	Freq12 :: [0],
	Freq13 :: [0, 1],
	Freq14 :: [0],
	Freq15 :: [0],
	Freq16 :: [0],
	Freq17 :: [0],
	Freq18 :: [0],
	Freq19 :: [0],
	Freq20 :: [0],
	Freq21 :: [0],
	Freq22 :: [0, 1],
	Freq23 :: [0],
	Freq24 :: [0, 1],
	Freq25 :: [0],
	Freq26 :: [0, 1],
	Freq27 :: [0, 1],
	Freq28 :: [0],
	Freq29 :: [0],
	Freq30 :: [0],
	Freq31 :: [0],
	Freq32 :: [0],
	Freq33 :: [0, 1],
	Freq34 :: [0, 1],
	Freq35 :: [0, 1],
	Freq36 :: [0, 1],
	Freq37 :: [0],
	Freq38 :: [0, 1],
	Freq39 :: [0],
	Freq40 :: [0],
	Freq41 :: [0, 1],
	Freq42 :: [0],
	Freq43 :: [0, 1],
	Freq44 :: [0],
	Freq45 :: [0, 1],
	Freq46 :: [0, 1],
	Freq47 :: [0, 1],
	Freq48 :: [0],
	Freq49 :: [0],
	Freq50 :: [0],
	Freq51 :: [0],
	Freq52 :: [0],
	Freq53 :: [0],
	Freq54 :: [0],
	Freq55 :: [0],
	Freq56 :: [0],
	Freq57 :: [0],
	Freq58 :: [0, 1],
	Freq59 :: [0],
	Freq60 :: [0],
	Freq61 :: [0],
	Freq62 :: [0],
	Freq63 :: [0],
	Freq64 :: [0],
	Freq65 :: [0],
	Freq66 :: [0],
	Freq67 :: [0],
	Freq68 :: [0],
	Freq69 :: [0],
	Freq70 :: [0],
	Freq71 :: [0],
	Freq72 :: [0],
	Freq73 :: [0],
	Freq74 :: [0],
	Freq75 :: [0],
	Freq76 :: [0],
	Freq77 :: [0],
	Freq78 :: [0],
	Freq79 :: [0],
	Freq80 :: [0],
	Freq81 :: [0],
	Freq82 :: [0],
	Freq83 :: [0],
	Freq84 :: [0],
	Freq85 :: [0],
	Freq86 :: [0],
	Freq87 :: [0],
	Freq88 :: [0],
	Freq89 :: [0],
	Freq90 :: [0],
	Freq91 :: [0],
	Freq92 :: [0],
	Freq93 :: [0],
	Freq94 :: [0],
	Freq95 :: [0],
	Freq96 :: [0],
	Freq97 :: [0],
	Freq98 :: [0],
	Freq99 :: [0],
	Freq100 :: [0],
	Freq101 :: [0],
	Freq102 :: [0],
	Freq103 :: [0],
	Freq104 :: [0],
	Freq105 :: [0],
	Freq106 :: [0],
	Freq107 :: [0, 1],
	Freq108 :: [0],
	Freq109 :: [0],
	Freq110 :: [0],
	Freq111 :: [0],
	Freq112 :: [0],
	Freq113 :: [0],
	Freq114 :: [0],
	Freq115 :: [0],
	Freq116 :: [0],
	Freq117 :: [0],
	Freq118 :: [0],
	Freq119 :: [0, 1],
	Freq120 :: [0, 1],
	Freq121 :: [0],
	Freq122 :: [0],
	Freq123 :: [0, 1],
	Freq124 :: [0, 1],
	Freq125 :: [0, 1],
	Freq126 :: [0, 1],
	Freq127 :: [0, 1],
	Freq128 :: [0, 1],
	Freq129 :: [0, 1],
	Freq130 :: [0, 1, 2],
	Freq131 :: [0, 1, 2],
	Freq132 :: [0, 1, 2],
	Freq133 :: [0, 1, 2],
	Freq134 :: [0, 1],
	Freq135 :: [0, 1],
	Freq136 :: [0, 1],
	Freq137 :: [0, 1],
	Freq138 :: [0, 1],
	Freq139 :: [0, 1],
	Freq140 :: [0, 1],
	Freq141 :: [0, 1],
	Freq142 :: [0, 1],
	Freq143 :: [0, 1],
	Freq144 :: [0, 1],
	Freq145 :: [0, 1],
	Freq146 :: [0, 1],
	Freq147 :: [0, 1, 2],
	Freq148 :: [0, 1, 2],
	Freq149 :: [0, 1, 2],
	Freq150 :: [0, 1, 2],
	Freq151 :: [0, 1, 2],
	Freq152 :: [0, 1, 2],
	Freq153 :: [0, 1, 2],
	Freq154 :: [0, 1, 2, 3],
	Freq155 :: [0, 1, 2, 3, 4],
	Freq156 :: [0, 1, 2, 3, 4],
	Freq157 :: [0, 1, 2, 3, 4],
	Freq158 :: [0, 1, 2, 3, 4, 5],
	Freq159 :: [0, 1, 2, 3, 4, 5, 6],
	Freq160 :: [0],
	Freq161 :: [0],
	Freq162 :: [0],
	Freq163 :: [0],
	Freq164 :: [0],
	Freq165 :: [0],
	Freq166 :: [0, 1, 2, 3, 4, 5, 6],
	Freq167 :: [0, 1, 2, 3, 4, 5, 6],
	Freq168 :: [0, 1, 2, 3, 4, 5],
	Freq169 :: [0, 1, 2, 3, 4],
	Freq170 :: [0, 1, 2, 3],
	Freq171 :: [0, 1, 2, 3],
	Freq172 :: [0, 1, 2],
	Freq173 :: [0, 1, 2],
	Freq174 :: [0, 1],
	Freq175 :: [0, 1],
	Freq176 :: [0, 1, 2],
	Freq177 :: [0, 1, 2],
	Freq178 :: [0, 1, 2],
	Freq179 :: [0, 1],
	Freq180 :: [0, 1],
	Freq181 :: [0, 1],
	Freq182 :: [0, 1],
	Freq183 :: [0, 1],
	Freq184 :: [0, 1],
	Freq185 :: [0, 1],
	Freq186 :: [0, 1],
	Freq187 :: [0, 1],
	Freq188 :: [0, 1],
	Freq189 :: [0],
	Freq190 :: [0],
	Freq191 :: [0],
	Freq192 :: [0, 1],
	Freq193 :: [0],
	Freq194 :: [0, 1],
	Freq195 :: [0, 1],
	Freq196 :: [0],
	Freq197 :: [0, 1],
	Freq198 :: [0],
	Freq199 :: [0],
	Freq200 :: [0],
	Freq201 :: [0],
	Freq202 :: [0],
	Freq203 :: [0],
	Freq204 :: [0],
	Freq205 :: [0],
	Freq206 :: [0],
	Freq207 :: [0],
	Freq208 :: [0],
	Freq209 :: [0],
	Freq210 :: [0],
	Freq211 :: [0],
	Freq212 :: [0],
	Freq213 :: [0],
	Freq214 :: [0],
	Freq215 :: [0],
	Freq216 :: [0],
	Freq217 :: [0],
	Freq218 :: [0],
	Freq219 :: [0],
	Freq220 :: [0],
	Freq221 :: [0],
	Freq222 :: [0],
	Freq223 :: [0],
	Freq224 :: [0],
	Freq225 :: [0],
	Freq226 :: [0],
	Freq227 :: [0],
	Freq228 :: [0],
	Freq229 :: [0],
	Freq230 :: [0],
	Freq231 :: [0],
	Freq232 :: [0],
	Freq233 :: [0],
	Freq234 :: [0],
	Freq235 :: [0],
	Freq236 :: [0],
	Freq237 :: [0],
	Freq238 :: [0],
	Freq239 :: [0],
	Freq240 :: [0],
	Freq241 :: [0],
	Freq242 :: [0],
	Freq243 :: [0],
	Freq244 :: [0],
	Freq245 :: [0],
	Freq246 :: [0],
	Freq247 :: [0],
	Freq248 :: [0],
	Freq249 :: [0],
	Freq250 :: [0],
	Freq251 :: [0],
	Freq252 :: [0],
	Freq253 :: [0],
	Freq254 :: [0],
	Freq255 :: [0],
	Freq256 :: [0],
	Freq257 :: [0],
	Freq258 :: [0],
	Freq259 :: [0],
	Freq260 :: [0],
	Freq261 :: [0],
	Freq262 :: [0],
	Freq263 :: [0],
	Freq264 :: [0],
	Freq265 :: [0, 1],
	Freq266 :: [0],
	Freq267 :: [0],
	Freq268 :: [0],
	Freq269 :: [0],
	Freq270 :: [0],
	Freq271 :: [0],
	Freq272 :: [0],
	Freq273 :: [0],
	Freq274 :: [0],
	Freq275 :: [0],
	Freq276 :: [0],
	Freq277 :: [0],
	Freq278 :: [0],
	Freq279 :: [0],
	Freq280 :: [0],
	Freq281 :: [0],
	Freq282 :: [0],
	Freq283 :: [0],
	Freq284 :: [0],
	Freq285 :: [0],
	Freq286 :: [0],
	Freq287 :: [0],
	Freq288 :: [0],
	Freq289 :: [0],
	Freq290 :: [0],
	Freq291 :: [0],
	Freq292 :: [0],
	Freq293 :: [0],
	Freq294 :: [0],
	Freq295 :: [0],
	Freq296 :: [0],
	Freq297 :: [0],
	Freq298 :: [0],
	Freq299 :: [0],
	Freq300 :: [0],
	Freq301 :: [0],
	Freq302 :: [0],
	Freq303 :: [0],
	Freq304 :: [0],
	Freq305 :: [0],
	Freq306 :: [0],
	Freq307 :: [0],
	Freq308 :: [0],
	Freq309 :: [0],
	Freq310 :: [0],
	Freq311 :: [0],
	Freq312 :: [0],
	Freq313 :: [0],
	Freq314 :: [0],
	Freq315 :: [0],
	Freq316 :: [0],
	Freq317 :: [0],
	Freq318 :: [0],
	Freq319 :: [0, 1],
	Freq320 :: [0],
	Freq321 :: [0],
	Freq322 :: [0],
	Freq323 :: [0],
	Freq324 :: [0],
	Freq325 :: [0],
	Freq326 :: [0],
	Freq327 :: [0],
	Freq328 :: [0],
	Freq329 :: [0],
	Freq330 :: [0],
	Freq331 :: [0],
	Freq332 :: [0],
	Freq333 :: [0],
	Freq334 :: [0],
	Freq335 :: [0],
	Freq336 :: [0],
	Freq337 :: [0],
	Freq338 :: [0],
	Freq339 :: [0],
	Freq340 :: [0, 1],
	Freq341 :: [0],
	Freq342 :: [0],
	Freq343 :: [0],
	Freq344 :: [0],
	Freq345 :: [0],
	Freq346 :: [0],
	Freq347 :: [0],
	Freq348 :: [0],
	Freq349 :: [0],
	Freq350 :: [0, 1],
	Freq351 :: [0],
	Freq352 :: [0],
	Freq353 :: [0],
	Freq354 :: [0],
	Freq355 :: [0],
	Freq356 :: [0],
	Freq357 :: [0],
	Freq358 :: [0],
	Freq359 :: [0, 1],
	Freq360 :: [0],
	Freq361 :: [0],
	Freq362 :: [0],
	Freq363 :: [0],
	Freq364 :: [0, 1],
	Freq365 :: [0, 1],
	Freq366 :: [0, 1],
	Freq367 :: [0, 1],
	Freq368 :: [0, 1],
	Freq369 :: [0, 1],
	Freq370 :: [0],
	Freq371 :: [0],
	Freq372 :: [0],
	Freq373 :: [0],
	Freq374 :: [0],
	Freq375 :: [0],
	Freq376 :: [0],
	Freq377 :: [0],
	Freq378 :: [0],
	Freq379 :: [0],
	Freq380 :: [0],
	Freq381 :: [0],
	Freq382 :: [0],
	Freq383 :: [0],
	Freq384 :: [0],
	Freq385 :: [0],
	Freq386 :: [0],
	Freq387 :: [0],
	Freq388 :: [0],
	Freq389 :: [0],
	Freq390 :: [0],
	Freq391 :: [0],
	Freq392 :: [0],
	Freq393 :: [0, 1],
	Freq394 :: [0],
	Freq395 :: [0],
	Freq396 :: [0],
	Freq397 :: [0],
	Freq398 :: [0],
	Freq399 :: [0],
	Freq400 :: [0],
	Freq401 :: [0],
	Freq402 :: [0],
	Freq403 :: [0],
	Freq404 :: [0],
	Freq405 :: [0],
	Freq406 :: [0, 1],
	Freq407 :: [0],
	Freq408 :: [0],
	Freq409 :: [0, 1],
	Freq410 :: [0, 1],
	Freq411 :: [0, 1],
	Freq412 :: [0, 1],
	Freq413 :: [0, 1],
	Freq414 :: [0, 1],
	Freq415 :: [0],
	Freq416 :: [0, 1],
	Freq417 :: [0, 1],
	Freq418 :: [0, 1],
	Freq419 :: [0],
	Freq420 :: [0, 1],
	Freq421 :: [0, 1],
	Freq422 :: [0, 1],
	Freq423 :: [0, 1],
	Freq424 :: [0, 1],
	Freq425 :: [0, 1],
	Freq426 :: [0, 2],
	Freq427 :: [0, 1],
	Freq428 :: [0, 1],
	Freq429 :: [0, 1],
	Freq430 :: [0, 1],
	Freq431 :: [0],
	Freq432 :: [0, 1],
	Freq433 :: [0],
	Freq434 :: [0, 1],
	Freq435 :: [0],
	Freq436 :: [0, 1],
	Freq437 :: [0, 1],
	Freq438 :: [0, 1],
	Freq439 :: [0, 1],
	Freq440 :: [0, 1],
	Freq441 :: [0, 1],
	Freq442 :: [0, 1],
	Freq443 :: [0],
	Freq444 :: [0, 1],
	Freq445 :: [0, 3, 2],
	Freq446 :: [0, 1],
	Freq447 :: [0, 1],
	Freq448 :: [0, 1],
	Freq449 :: [0, 1],
	Freq450 :: [0, 1],
	Freq451 :: [0],
	Freq452 :: [0],
	Freq453 :: [0],
	Freq454 :: [0],


	% Constraints: 
	% Each pair of corresponding (Row<i>, Freq<i>) variables 
	% must assume dependent values based on data from the 
	% whole-genome contact map; A (Row, Freq) pair ground to 
	% (0,0) encodes that nothing is chosen 

	((Row1 #= 0 and Freq1 #= 0)), 

	((Row2 #= 0 and Freq2 #= 0)), 

	((Row3 #= 0 and Freq3 #= 0)), 

	((Row4 #= 245 and Freq4 #= 1) or
	(Row4 #= 0 and Freq4 #= 0)), 

	((Row5 #= 245 and Freq5 #= 1) or
	(Row5 #= 0 and Freq5 #= 0)), 

	((Row6 #= 246 and Freq6 #= 1) or
	(Row6 #= 0 and Freq6 #= 0)), 

	((Row7 #= 0 and Freq7 #= 0)), 

	((Row8 #= 246 and Freq8 #= 1) or
	(Row8 #= 0 and Freq8 #= 0)), 

	((Row9 #= 0 and Freq9 #= 0)), 

	((Row10 #= 245 and Freq10 #= 1) or
	(Row10 #= 0 and Freq10 #= 0)), 

	((Row11 #= 0 and Freq11 #= 0)), 

	((Row12 #= 0 and Freq12 #= 0)), 

	((Row13 #= 245 and Freq13 #= 1) or
	(Row13 #= 0 and Freq13 #= 0)), 

	((Row14 #= 0 and Freq14 #= 0)), 

	((Row15 #= 0 and Freq15 #= 0)), 

	((Row16 #= 0 and Freq16 #= 0)), 

	((Row17 #= 0 and Freq17 #= 0)), 

	((Row18 #= 0 and Freq18 #= 0)), 

	((Row19 #= 0 and Freq19 #= 0)), 

	((Row20 #= 0 and Freq20 #= 0)), 

	((Row21 #= 0 and Freq21 #= 0)), 

	((Row22 #= 246 and Freq22 #= 1) or
	(Row22 #= 0 and Freq22 #= 0)), 

	((Row23 #= 0 and Freq23 #= 0)), 

	((Row24 #= 245 and Freq24 #= 1) or
	(Row24 #= 0 and Freq24 #= 0)), 

	((Row25 #= 0 and Freq25 #= 0)), 

	((Row26 #= 245 and Freq26 #= 1) or
	(Row26 #= 0 and Freq26 #= 0)), 

	((Row27 #= 246 and Freq27 #= 1) or
	(Row27 #= 0 and Freq27 #= 0)), 

	((Row28 #= 0 and Freq28 #= 0)), 

	((Row29 #= 0 and Freq29 #= 0)), 

	((Row30 #= 0 and Freq30 #= 0)), 

	((Row31 #= 0 and Freq31 #= 0)), 

	((Row32 #= 0 and Freq32 #= 0)), 

	((Row33 #= 245 and Freq33 #= 1) or
	(Row33 #= 246 and Freq33 #= 1) or
	(Row33 #= 0 and Freq33 #= 0)), 

	((Row34 #= 245 and Freq34 #= 1) or
	(Row34 #= 0 and Freq34 #= 0)), 

	((Row35 #= 245 and Freq35 #= 1) or
	(Row35 #= 0 and Freq35 #= 0)), 

	((Row36 #= 245 and Freq36 #= 1) or
	(Row36 #= 246 and Freq36 #= 1) or
	(Row36 #= 0 and Freq36 #= 0)), 

	((Row37 #= 0 and Freq37 #= 0)), 

	((Row38 #= 245 and Freq38 #= 1) or
	(Row38 #= 0 and Freq38 #= 0)), 

	((Row39 #= 0 and Freq39 #= 0)), 

	((Row40 #= 0 and Freq40 #= 0)), 

	((Row41 #= 245 and Freq41 #= 1) or
	(Row41 #= 0 and Freq41 #= 0)), 

	((Row42 #= 0 and Freq42 #= 0)), 

	((Row43 #= 245 and Freq43 #= 1) or
	(Row43 #= 0 and Freq43 #= 0)), 

	((Row44 #= 0 and Freq44 #= 0)), 

	((Row45 #= 245 and Freq45 #= 1) or
	(Row45 #= 0 and Freq45 #= 0)), 

	((Row46 #= 245 and Freq46 #= 1) or
	(Row46 #= 0 and Freq46 #= 0)), 

	((Row47 #= 245 and Freq47 #= 1) or
	(Row47 #= 0 and Freq47 #= 0)), 

	((Row48 #= 0 and Freq48 #= 0)), 

	((Row49 #= 0 and Freq49 #= 0)), 

	((Row50 #= 0 and Freq50 #= 0)), 

	((Row51 #= 0 and Freq51 #= 0)), 

	((Row52 #= 0 and Freq52 #= 0)), 

	((Row53 #= 0 and Freq53 #= 0)), 

	((Row54 #= 0 and Freq54 #= 0)), 

	((Row55 #= 0 and Freq55 #= 0)), 

	((Row56 #= 0 and Freq56 #= 0)), 

	((Row57 #= 0 and Freq57 #= 0)), 

	((Row58 #= 82 and Freq58 #= 1) or
	(Row58 #= 0 and Freq58 #= 0)), 

	((Row59 #= 0 and Freq59 #= 0)), 

	((Row60 #= 0 and Freq60 #= 0)), 

	((Row61 #= 0 and Freq61 #= 0)), 

	((Row62 #= 0 and Freq62 #= 0)), 

	((Row63 #= 0 and Freq63 #= 0)), 

	((Row64 #= 0 and Freq64 #= 0)), 

	((Row65 #= 0 and Freq65 #= 0)), 

	((Row66 #= 0 and Freq66 #= 0)), 

	((Row67 #= 0 and Freq67 #= 0)), 

	((Row68 #= 0 and Freq68 #= 0)), 

	((Row69 #= 0 and Freq69 #= 0)), 

	((Row70 #= 0 and Freq70 #= 0)), 

	((Row71 #= 0 and Freq71 #= 0)), 

	((Row72 #= 0 and Freq72 #= 0)), 

	((Row73 #= 0 and Freq73 #= 0)), 

	((Row74 #= 0 and Freq74 #= 0)), 

	((Row75 #= 0 and Freq75 #= 0)), 

	((Row76 #= 0 and Freq76 #= 0)), 

	((Row77 #= 0 and Freq77 #= 0)), 

	((Row78 #= 0 and Freq78 #= 0)), 

	((Row79 #= 0 and Freq79 #= 0)), 

	((Row80 #= 0 and Freq80 #= 0)), 

	((Row81 #= 0 and Freq81 #= 0)), 

	((Row82 #= 0 and Freq82 #= 0)), 

	((Row83 #= 0 and Freq83 #= 0)), 

	((Row84 #= 0 and Freq84 #= 0)), 

	((Row85 #= 0 and Freq85 #= 0)), 

	((Row86 #= 0 and Freq86 #= 0)), 

	((Row87 #= 0 and Freq87 #= 0)), 

	((Row88 #= 0 and Freq88 #= 0)), 

	((Row89 #= 0 and Freq89 #= 0)), 

	((Row90 #= 0 and Freq90 #= 0)), 

	((Row91 #= 0 and Freq91 #= 0)), 

	((Row92 #= 0 and Freq92 #= 0)), 

	((Row93 #= 0 and Freq93 #= 0)), 

	((Row94 #= 0 and Freq94 #= 0)), 

	((Row95 #= 0 and Freq95 #= 0)), 

	((Row96 #= 0 and Freq96 #= 0)), 

	((Row97 #= 0 and Freq97 #= 0)), 

	((Row98 #= 0 and Freq98 #= 0)), 

	((Row99 #= 0 and Freq99 #= 0)), 

	((Row100 #= 0 and Freq100 #= 0)), 

	((Row101 #= 0 and Freq101 #= 0)), 

	((Row102 #= 0 and Freq102 #= 0)), 

	((Row103 #= 0 and Freq103 #= 0)), 

	((Row104 #= 0 and Freq104 #= 0)), 

	((Row105 #= 0 and Freq105 #= 0)), 

	((Row106 #= 0 and Freq106 #= 0)), 

	((Row107 #= 245 and Freq107 #= 1) or
	(Row107 #= 0 and Freq107 #= 0)), 

	((Row108 #= 0 and Freq108 #= 0)), 

	((Row109 #= 0 and Freq109 #= 0)), 

	((Row110 #= 0 and Freq110 #= 0)), 

	((Row111 #= 0 and Freq111 #= 0)), 

	((Row112 #= 0 and Freq112 #= 0)), 

	((Row113 #= 0 and Freq113 #= 0)), 

	((Row114 #= 0 and Freq114 #= 0)), 

	((Row115 #= 0 and Freq115 #= 0)), 

	((Row116 #= 0 and Freq116 #= 0)), 

	((Row117 #= 0 and Freq117 #= 0)), 

	((Row118 #= 0 and Freq118 #= 0)), 

	((Row119 #= 152 and Freq119 #= 1) or
	(Row119 #= 0 and Freq119 #= 0)), 

	((Row120 #= 145 and Freq120 #= 1) or
	(Row120 #= 0 and Freq120 #= 0)), 

	((Row121 #= 0 and Freq121 #= 0)), 

	((Row122 #= 0 and Freq122 #= 0)), 

	((Row123 #= 245 and Freq123 #= 1) or
	(Row123 #= 0 and Freq123 #= 0)), 

	((Row124 #= 144 and Freq124 #= 1) or
	(Row124 #= 246 and Freq124 #= 1) or
	(Row124 #= 0 and Freq124 #= 0)), 

	((Row125 #= 138 and Freq125 #= 1) or
	(Row125 #= 143 and Freq125 #= 1) or
	(Row125 #= 144 and Freq125 #= 1) or
	(Row125 #= 145 and Freq125 #= 1) or
	(Row125 #= 146 and Freq125 #= 1) or
	(Row125 #= 148 and Freq125 #= 1) or
	(Row125 #= 0 and Freq125 #= 0)), 

	((Row126 #= 138 and Freq126 #= 1) or
	(Row126 #= 141 and Freq126 #= 1) or
	(Row126 #= 142 and Freq126 #= 1) or
	(Row126 #= 143 and Freq126 #= 1) or
	(Row126 #= 144 and Freq126 #= 1) or
	(Row126 #= 145 and Freq126 #= 1) or
	(Row126 #= 147 and Freq126 #= 1) or
	(Row126 #= 0 and Freq126 #= 0)), 

	((Row127 #= 140 and Freq127 #= 1) or
	(Row127 #= 141 and Freq127 #= 1) or
	(Row127 #= 142 and Freq127 #= 1) or
	(Row127 #= 143 and Freq127 #= 1) or
	(Row127 #= 144 and Freq127 #= 1) or
	(Row127 #= 145 and Freq127 #= 1) or
	(Row127 #= 146 and Freq127 #= 1) or
	(Row127 #= 147 and Freq127 #= 1) or
	(Row127 #= 148 and Freq127 #= 1) or
	(Row127 #= 150 and Freq127 #= 1) or
	(Row127 #= 153 and Freq127 #= 1) or
	(Row127 #= 0 and Freq127 #= 0)), 

	((Row128 #= 137 and Freq128 #= 1) or
	(Row128 #= 143 and Freq128 #= 1) or
	(Row128 #= 144 and Freq128 #= 1) or
	(Row128 #= 145 and Freq128 #= 1) or
	(Row128 #= 146 and Freq128 #= 1) or
	(Row128 #= 147 and Freq128 #= 1) or
	(Row128 #= 148 and Freq128 #= 1) or
	(Row128 #= 150 and Freq128 #= 1) or
	(Row128 #= 0 and Freq128 #= 0)), 

	((Row129 #= 138 and Freq129 #= 1) or
	(Row129 #= 143 and Freq129 #= 1) or
	(Row129 #= 144 and Freq129 #= 1) or
	(Row129 #= 145 and Freq129 #= 1) or
	(Row129 #= 146 and Freq129 #= 1) or
	(Row129 #= 147 and Freq129 #= 1) or
	(Row129 #= 148 and Freq129 #= 1) or
	(Row129 #= 0 and Freq129 #= 0)), 

	((Row130 #= 138 and Freq130 #= 1) or
	(Row130 #= 140 and Freq130 #= 1) or
	(Row130 #= 141 and Freq130 #= 1) or
	(Row130 #= 142 and Freq130 #= 1) or
	(Row130 #= 143 and Freq130 #= 1) or
	(Row130 #= 144 and Freq130 #= 1) or
	(Row130 #= 145 and Freq130 #= 1) or
	(Row130 #= 146 and Freq130 #= 2) or
	(Row130 #= 147 and Freq130 #= 1) or
	(Row130 #= 148 and Freq130 #= 1) or
	(Row130 #= 149 and Freq130 #= 1) or
	(Row130 #= 0 and Freq130 #= 0)), 

	((Row131 #= 138 and Freq131 #= 1) or
	(Row131 #= 140 and Freq131 #= 1) or
	(Row131 #= 141 and Freq131 #= 1) or
	(Row131 #= 142 and Freq131 #= 1) or
	(Row131 #= 143 and Freq131 #= 1) or
	(Row131 #= 144 and Freq131 #= 1) or
	(Row131 #= 145 and Freq131 #= 2) or
	(Row131 #= 146 and Freq131 #= 2) or
	(Row131 #= 147 and Freq131 #= 1) or
	(Row131 #= 148 and Freq131 #= 1) or
	(Row131 #= 149 and Freq131 #= 1) or
	(Row131 #= 150 and Freq131 #= 1) or
	(Row131 #= 152 and Freq131 #= 1) or
	(Row131 #= 0 and Freq131 #= 0)), 

	((Row132 #= 39 and Freq132 #= 1) or
	(Row132 #= 130 and Freq132 #= 1) or
	(Row132 #= 133 and Freq132 #= 1) or
	(Row132 #= 136 and Freq132 #= 1) or
	(Row132 #= 137 and Freq132 #= 1) or
	(Row132 #= 138 and Freq132 #= 1) or
	(Row132 #= 140 and Freq132 #= 1) or
	(Row132 #= 141 and Freq132 #= 1) or
	(Row132 #= 142 and Freq132 #= 1) or
	(Row132 #= 143 and Freq132 #= 1) or
	(Row132 #= 144 and Freq132 #= 2) or
	(Row132 #= 145 and Freq132 #= 2) or
	(Row132 #= 146 and Freq132 #= 2) or
	(Row132 #= 147 and Freq132 #= 2) or
	(Row132 #= 148 and Freq132 #= 1) or
	(Row132 #= 149 and Freq132 #= 1) or
	(Row132 #= 150 and Freq132 #= 1) or
	(Row132 #= 151 and Freq132 #= 1) or
	(Row132 #= 152 and Freq132 #= 1) or
	(Row132 #= 155 and Freq132 #= 1) or
	(Row132 #= 245 and Freq132 #= 1) or
	(Row132 #= 0 and Freq132 #= 0)), 

	((Row133 #= 130 and Freq133 #= 1) or
	(Row133 #= 134 and Freq133 #= 1) or
	(Row133 #= 135 and Freq133 #= 1) or
	(Row133 #= 136 and Freq133 #= 1) or
	(Row133 #= 137 and Freq133 #= 1) or
	(Row133 #= 138 and Freq133 #= 1) or
	(Row133 #= 140 and Freq133 #= 1) or
	(Row133 #= 141 and Freq133 #= 1) or
	(Row133 #= 142 and Freq133 #= 1) or
	(Row133 #= 143 and Freq133 #= 1) or
	(Row133 #= 144 and Freq133 #= 1) or
	(Row133 #= 145 and Freq133 #= 2) or
	(Row133 #= 146 and Freq133 #= 2) or
	(Row133 #= 147 and Freq133 #= 2) or
	(Row133 #= 148 and Freq133 #= 1) or
	(Row133 #= 149 and Freq133 #= 1) or
	(Row133 #= 150 and Freq133 #= 1) or
	(Row133 #= 151 and Freq133 #= 1) or
	(Row133 #= 156 and Freq133 #= 1) or
	(Row133 #= 245 and Freq133 #= 1) or
	(Row133 #= 0 and Freq133 #= 0)), 

	((Row134 #= 101 and Freq134 #= 1) or
	(Row134 #= 129 and Freq134 #= 1) or
	(Row134 #= 135 and Freq134 #= 1) or
	(Row134 #= 136 and Freq134 #= 1) or
	(Row134 #= 137 and Freq134 #= 1) or
	(Row134 #= 138 and Freq134 #= 1) or
	(Row134 #= 139 and Freq134 #= 1) or
	(Row134 #= 140 and Freq134 #= 1) or
	(Row134 #= 141 and Freq134 #= 1) or
	(Row134 #= 142 and Freq134 #= 1) or
	(Row134 #= 143 and Freq134 #= 1) or
	(Row134 #= 144 and Freq134 #= 1) or
	(Row134 #= 145 and Freq134 #= 1) or
	(Row134 #= 146 and Freq134 #= 1) or
	(Row134 #= 147 and Freq134 #= 1) or
	(Row134 #= 148 and Freq134 #= 1) or
	(Row134 #= 149 and Freq134 #= 1) or
	(Row134 #= 150 and Freq134 #= 1) or
	(Row134 #= 0 and Freq134 #= 0)), 

	((Row135 #= 129 and Freq135 #= 1) or
	(Row135 #= 131 and Freq135 #= 1) or
	(Row135 #= 134 and Freq135 #= 1) or
	(Row135 #= 135 and Freq135 #= 1) or
	(Row135 #= 136 and Freq135 #= 1) or
	(Row135 #= 137 and Freq135 #= 1) or
	(Row135 #= 138 and Freq135 #= 1) or
	(Row135 #= 140 and Freq135 #= 1) or
	(Row135 #= 141 and Freq135 #= 1) or
	(Row135 #= 142 and Freq135 #= 1) or
	(Row135 #= 143 and Freq135 #= 1) or
	(Row135 #= 144 and Freq135 #= 1) or
	(Row135 #= 145 and Freq135 #= 1) or
	(Row135 #= 146 and Freq135 #= 1) or
	(Row135 #= 147 and Freq135 #= 1) or
	(Row135 #= 148 and Freq135 #= 1) or
	(Row135 #= 149 and Freq135 #= 1) or
	(Row135 #= 150 and Freq135 #= 1) or
	(Row135 #= 0 and Freq135 #= 0)), 

	((Row136 #= 88 and Freq136 #= 1) or
	(Row136 #= 92 and Freq136 #= 1) or
	(Row136 #= 93 and Freq136 #= 1) or
	(Row136 #= 94 and Freq136 #= 1) or
	(Row136 #= 102 and Freq136 #= 1) or
	(Row136 #= 105 and Freq136 #= 1) or
	(Row136 #= 129 and Freq136 #= 1) or
	(Row136 #= 131 and Freq136 #= 1) or
	(Row136 #= 132 and Freq136 #= 1) or
	(Row136 #= 134 and Freq136 #= 1) or
	(Row136 #= 135 and Freq136 #= 1) or
	(Row136 #= 136 and Freq136 #= 1) or
	(Row136 #= 137 and Freq136 #= 1) or
	(Row136 #= 138 and Freq136 #= 1) or
	(Row136 #= 141 and Freq136 #= 1) or
	(Row136 #= 142 and Freq136 #= 1) or
	(Row136 #= 143 and Freq136 #= 1) or
	(Row136 #= 144 and Freq136 #= 1) or
	(Row136 #= 145 and Freq136 #= 1) or
	(Row136 #= 146 and Freq136 #= 1) or
	(Row136 #= 147 and Freq136 #= 1) or
	(Row136 #= 148 and Freq136 #= 1) or
	(Row136 #= 149 and Freq136 #= 1) or
	(Row136 #= 150 and Freq136 #= 1) or
	(Row136 #= 151 and Freq136 #= 1) or
	(Row136 #= 245 and Freq136 #= 1) or
	(Row136 #= 0 and Freq136 #= 0)), 

	((Row137 #= 97 and Freq137 #= 1) or
	(Row137 #= 99 and Freq137 #= 1) or
	(Row137 #= 100 and Freq137 #= 1) or
	(Row137 #= 102 and Freq137 #= 1) or
	(Row137 #= 105 and Freq137 #= 1) or
	(Row137 #= 121 and Freq137 #= 1) or
	(Row137 #= 126 and Freq137 #= 1) or
	(Row137 #= 127 and Freq137 #= 1) or
	(Row137 #= 129 and Freq137 #= 1) or
	(Row137 #= 132 and Freq137 #= 1) or
	(Row137 #= 133 and Freq137 #= 1) or
	(Row137 #= 134 and Freq137 #= 1) or
	(Row137 #= 135 and Freq137 #= 1) or
	(Row137 #= 136 and Freq137 #= 1) or
	(Row137 #= 137 and Freq137 #= 1) or
	(Row137 #= 138 and Freq137 #= 1) or
	(Row137 #= 140 and Freq137 #= 1) or
	(Row137 #= 141 and Freq137 #= 1) or
	(Row137 #= 142 and Freq137 #= 1) or
	(Row137 #= 143 and Freq137 #= 1) or
	(Row137 #= 144 and Freq137 #= 1) or
	(Row137 #= 145 and Freq137 #= 1) or
	(Row137 #= 146 and Freq137 #= 1) or
	(Row137 #= 147 and Freq137 #= 1) or
	(Row137 #= 148 and Freq137 #= 1) or
	(Row137 #= 154 and Freq137 #= 1) or
	(Row137 #= 245 and Freq137 #= 1) or
	(Row137 #= 0 and Freq137 #= 0)), 

	((Row138 #= 96 and Freq138 #= 1) or
	(Row138 #= 97 and Freq138 #= 1) or
	(Row138 #= 98 and Freq138 #= 1) or
	(Row138 #= 101 and Freq138 #= 1) or
	(Row138 #= 104 and Freq138 #= 1) or
	(Row138 #= 105 and Freq138 #= 1) or
	(Row138 #= 106 and Freq138 #= 1) or
	(Row138 #= 127 and Freq138 #= 1) or
	(Row138 #= 129 and Freq138 #= 1) or
	(Row138 #= 131 and Freq138 #= 1) or
	(Row138 #= 134 and Freq138 #= 1) or
	(Row138 #= 136 and Freq138 #= 1) or
	(Row138 #= 137 and Freq138 #= 1) or
	(Row138 #= 139 and Freq138 #= 1) or
	(Row138 #= 140 and Freq138 #= 1) or
	(Row138 #= 141 and Freq138 #= 1) or
	(Row138 #= 142 and Freq138 #= 1) or
	(Row138 #= 143 and Freq138 #= 1) or
	(Row138 #= 144 and Freq138 #= 1) or
	(Row138 #= 145 and Freq138 #= 1) or
	(Row138 #= 146 and Freq138 #= 1) or
	(Row138 #= 147 and Freq138 #= 1) or
	(Row138 #= 148 and Freq138 #= 1) or
	(Row138 #= 150 and Freq138 #= 1) or
	(Row138 #= 0 and Freq138 #= 0)), 

	((Row139 #= 91 and Freq139 #= 1) or
	(Row139 #= 97 and Freq139 #= 1) or
	(Row139 #= 99 and Freq139 #= 1) or
	(Row139 #= 100 and Freq139 #= 1) or
	(Row139 #= 104 and Freq139 #= 1) or
	(Row139 #= 121 and Freq139 #= 1) or
	(Row139 #= 126 and Freq139 #= 1) or
	(Row139 #= 128 and Freq139 #= 1) or
	(Row139 #= 129 and Freq139 #= 1) or
	(Row139 #= 130 and Freq139 #= 1) or
	(Row139 #= 132 and Freq139 #= 1) or
	(Row139 #= 133 and Freq139 #= 1) or
	(Row139 #= 134 and Freq139 #= 1) or
	(Row139 #= 135 and Freq139 #= 1) or
	(Row139 #= 136 and Freq139 #= 1) or
	(Row139 #= 137 and Freq139 #= 1) or
	(Row139 #= 138 and Freq139 #= 1) or
	(Row139 #= 140 and Freq139 #= 1) or
	(Row139 #= 141 and Freq139 #= 1) or
	(Row139 #= 142 and Freq139 #= 1) or
	(Row139 #= 143 and Freq139 #= 1) or
	(Row139 #= 144 and Freq139 #= 1) or
	(Row139 #= 145 and Freq139 #= 1) or
	(Row139 #= 146 and Freq139 #= 1) or
	(Row139 #= 147 and Freq139 #= 1) or
	(Row139 #= 148 and Freq139 #= 1) or
	(Row139 #= 150 and Freq139 #= 1) or
	(Row139 #= 0 and Freq139 #= 0)), 

	((Row140 #= 87 and Freq140 #= 1) or
	(Row140 #= 92 and Freq140 #= 1) or
	(Row140 #= 93 and Freq140 #= 1) or
	(Row140 #= 94 and Freq140 #= 1) or
	(Row140 #= 96 and Freq140 #= 1) or
	(Row140 #= 97 and Freq140 #= 1) or
	(Row140 #= 98 and Freq140 #= 1) or
	(Row140 #= 99 and Freq140 #= 1) or
	(Row140 #= 100 and Freq140 #= 1) or
	(Row140 #= 101 and Freq140 #= 1) or
	(Row140 #= 102 and Freq140 #= 1) or
	(Row140 #= 104 and Freq140 #= 1) or
	(Row140 #= 123 and Freq140 #= 1) or
	(Row140 #= 124 and Freq140 #= 1) or
	(Row140 #= 125 and Freq140 #= 1) or
	(Row140 #= 126 and Freq140 #= 1) or
	(Row140 #= 127 and Freq140 #= 1) or
	(Row140 #= 128 and Freq140 #= 1) or
	(Row140 #= 129 and Freq140 #= 1) or
	(Row140 #= 131 and Freq140 #= 1) or
	(Row140 #= 133 and Freq140 #= 1) or
	(Row140 #= 134 and Freq140 #= 1) or
	(Row140 #= 135 and Freq140 #= 1) or
	(Row140 #= 136 and Freq140 #= 1) or
	(Row140 #= 137 and Freq140 #= 1) or
	(Row140 #= 138 and Freq140 #= 1) or
	(Row140 #= 139 and Freq140 #= 1) or
	(Row140 #= 140 and Freq140 #= 1) or
	(Row140 #= 141 and Freq140 #= 1) or
	(Row140 #= 142 and Freq140 #= 1) or
	(Row140 #= 143 and Freq140 #= 1) or
	(Row140 #= 144 and Freq140 #= 1) or
	(Row140 #= 145 and Freq140 #= 1) or
	(Row140 #= 146 and Freq140 #= 1) or
	(Row140 #= 147 and Freq140 #= 1) or
	(Row140 #= 148 and Freq140 #= 1) or
	(Row140 #= 0 and Freq140 #= 0)), 

	((Row141 #= 86 and Freq141 #= 1) or
	(Row141 #= 87 and Freq141 #= 1) or
	(Row141 #= 88 and Freq141 #= 1) or
	(Row141 #= 89 and Freq141 #= 1) or
	(Row141 #= 90 and Freq141 #= 1) or
	(Row141 #= 91 and Freq141 #= 1) or
	(Row141 #= 92 and Freq141 #= 1) or
	(Row141 #= 93 and Freq141 #= 1) or
	(Row141 #= 94 and Freq141 #= 1) or
	(Row141 #= 95 and Freq141 #= 1) or
	(Row141 #= 96 and Freq141 #= 1) or
	(Row141 #= 97 and Freq141 #= 1) or
	(Row141 #= 98 and Freq141 #= 1) or
	(Row141 #= 99 and Freq141 #= 1) or
	(Row141 #= 100 and Freq141 #= 1) or
	(Row141 #= 101 and Freq141 #= 1) or
	(Row141 #= 102 and Freq141 #= 1) or
	(Row141 #= 103 and Freq141 #= 1) or
	(Row141 #= 104 and Freq141 #= 1) or
	(Row141 #= 105 and Freq141 #= 1) or
	(Row141 #= 106 and Freq141 #= 1) or
	(Row141 #= 118 and Freq141 #= 1) or
	(Row141 #= 121 and Freq141 #= 1) or
	(Row141 #= 122 and Freq141 #= 1) or
	(Row141 #= 123 and Freq141 #= 1) or
	(Row141 #= 125 and Freq141 #= 1) or
	(Row141 #= 126 and Freq141 #= 1) or
	(Row141 #= 127 and Freq141 #= 1) or
	(Row141 #= 128 and Freq141 #= 1) or
	(Row141 #= 129 and Freq141 #= 1) or
	(Row141 #= 130 and Freq141 #= 1) or
	(Row141 #= 132 and Freq141 #= 1) or
	(Row141 #= 133 and Freq141 #= 1) or
	(Row141 #= 134 and Freq141 #= 1) or
	(Row141 #= 135 and Freq141 #= 1) or
	(Row141 #= 136 and Freq141 #= 1) or
	(Row141 #= 137 and Freq141 #= 1) or
	(Row141 #= 138 and Freq141 #= 1) or
	(Row141 #= 139 and Freq141 #= 1) or
	(Row141 #= 140 and Freq141 #= 1) or
	(Row141 #= 141 and Freq141 #= 1) or
	(Row141 #= 142 and Freq141 #= 1) or
	(Row141 #= 143 and Freq141 #= 1) or
	(Row141 #= 144 and Freq141 #= 1) or
	(Row141 #= 145 and Freq141 #= 1) or
	(Row141 #= 146 and Freq141 #= 1) or
	(Row141 #= 147 and Freq141 #= 1) or
	(Row141 #= 148 and Freq141 #= 1) or
	(Row141 #= 149 and Freq141 #= 1) or
	(Row141 #= 150 and Freq141 #= 1) or
	(Row141 #= 152 and Freq141 #= 1) or
	(Row141 #= 153 and Freq141 #= 1) or
	(Row141 #= 0 and Freq141 #= 0)), 

	((Row142 #= 88 and Freq142 #= 1) or
	(Row142 #= 89 and Freq142 #= 1) or
	(Row142 #= 90 and Freq142 #= 1) or
	(Row142 #= 92 and Freq142 #= 1) or
	(Row142 #= 93 and Freq142 #= 1) or
	(Row142 #= 96 and Freq142 #= 1) or
	(Row142 #= 97 and Freq142 #= 1) or
	(Row142 #= 98 and Freq142 #= 1) or
	(Row142 #= 99 and Freq142 #= 1) or
	(Row142 #= 100 and Freq142 #= 1) or
	(Row142 #= 101 and Freq142 #= 1) or
	(Row142 #= 102 and Freq142 #= 1) or
	(Row142 #= 105 and Freq142 #= 1) or
	(Row142 #= 116 and Freq142 #= 1) or
	(Row142 #= 117 and Freq142 #= 1) or
	(Row142 #= 120 and Freq142 #= 1) or
	(Row142 #= 121 and Freq142 #= 1) or
	(Row142 #= 125 and Freq142 #= 1) or
	(Row142 #= 126 and Freq142 #= 1) or
	(Row142 #= 128 and Freq142 #= 1) or
	(Row142 #= 130 and Freq142 #= 1) or
	(Row142 #= 131 and Freq142 #= 1) or
	(Row142 #= 132 and Freq142 #= 1) or
	(Row142 #= 133 and Freq142 #= 1) or
	(Row142 #= 134 and Freq142 #= 1) or
	(Row142 #= 135 and Freq142 #= 1) or
	(Row142 #= 136 and Freq142 #= 1) or
	(Row142 #= 137 and Freq142 #= 1) or
	(Row142 #= 138 and Freq142 #= 1) or
	(Row142 #= 139 and Freq142 #= 1) or
	(Row142 #= 140 and Freq142 #= 1) or
	(Row142 #= 141 and Freq142 #= 1) or
	(Row142 #= 142 and Freq142 #= 1) or
	(Row142 #= 143 and Freq142 #= 1) or
	(Row142 #= 144 and Freq142 #= 1) or
	(Row142 #= 145 and Freq142 #= 1) or
	(Row142 #= 146 and Freq142 #= 1) or
	(Row142 #= 147 and Freq142 #= 1) or
	(Row142 #= 148 and Freq142 #= 1) or
	(Row142 #= 149 and Freq142 #= 1) or
	(Row142 #= 150 and Freq142 #= 1) or
	(Row142 #= 152 and Freq142 #= 1) or
	(Row142 #= 0 and Freq142 #= 0)), 

	((Row143 #= 87 and Freq143 #= 1) or
	(Row143 #= 90 and Freq143 #= 1) or
	(Row143 #= 94 and Freq143 #= 1) or
	(Row143 #= 95 and Freq143 #= 1) or
	(Row143 #= 96 and Freq143 #= 1) or
	(Row143 #= 97 and Freq143 #= 1) or
	(Row143 #= 98 and Freq143 #= 1) or
	(Row143 #= 99 and Freq143 #= 1) or
	(Row143 #= 100 and Freq143 #= 1) or
	(Row143 #= 101 and Freq143 #= 1) or
	(Row143 #= 102 and Freq143 #= 1) or
	(Row143 #= 103 and Freq143 #= 1) or
	(Row143 #= 104 and Freq143 #= 1) or
	(Row143 #= 105 and Freq143 #= 1) or
	(Row143 #= 106 and Freq143 #= 1) or
	(Row143 #= 115 and Freq143 #= 1) or
	(Row143 #= 116 and Freq143 #= 1) or
	(Row143 #= 117 and Freq143 #= 1) or
	(Row143 #= 118 and Freq143 #= 1) or
	(Row143 #= 119 and Freq143 #= 1) or
	(Row143 #= 120 and Freq143 #= 1) or
	(Row143 #= 121 and Freq143 #= 1) or
	(Row143 #= 122 and Freq143 #= 1) or
	(Row143 #= 123 and Freq143 #= 1) or
	(Row143 #= 124 and Freq143 #= 1) or
	(Row143 #= 125 and Freq143 #= 1) or
	(Row143 #= 127 and Freq143 #= 1) or
	(Row143 #= 129 and Freq143 #= 1) or
	(Row143 #= 130 and Freq143 #= 1) or
	(Row143 #= 131 and Freq143 #= 1) or
	(Row143 #= 132 and Freq143 #= 1) or
	(Row143 #= 133 and Freq143 #= 1) or
	(Row143 #= 134 and Freq143 #= 1) or
	(Row143 #= 135 and Freq143 #= 1) or
	(Row143 #= 136 and Freq143 #= 1) or
	(Row143 #= 137 and Freq143 #= 1) or
	(Row143 #= 139 and Freq143 #= 1) or
	(Row143 #= 140 and Freq143 #= 1) or
	(Row143 #= 141 and Freq143 #= 1) or
	(Row143 #= 142 and Freq143 #= 1) or
	(Row143 #= 143 and Freq143 #= 1) or
	(Row143 #= 144 and Freq143 #= 1) or
	(Row143 #= 145 and Freq143 #= 1) or
	(Row143 #= 146 and Freq143 #= 1) or
	(Row143 #= 147 and Freq143 #= 1) or
	(Row143 #= 150 and Freq143 #= 1) or
	(Row143 #= 0 and Freq143 #= 0)), 

	((Row144 #= 81 and Freq144 #= 1) or
	(Row144 #= 84 and Freq144 #= 1) or
	(Row144 #= 88 and Freq144 #= 1) or
	(Row144 #= 89 and Freq144 #= 1) or
	(Row144 #= 90 and Freq144 #= 1) or
	(Row144 #= 91 and Freq144 #= 1) or
	(Row144 #= 92 and Freq144 #= 1) or
	(Row144 #= 93 and Freq144 #= 1) or
	(Row144 #= 95 and Freq144 #= 1) or
	(Row144 #= 96 and Freq144 #= 1) or
	(Row144 #= 97 and Freq144 #= 1) or
	(Row144 #= 98 and Freq144 #= 1) or
	(Row144 #= 99 and Freq144 #= 1) or
	(Row144 #= 100 and Freq144 #= 1) or
	(Row144 #= 101 and Freq144 #= 1) or
	(Row144 #= 102 and Freq144 #= 1) or
	(Row144 #= 103 and Freq144 #= 1) or
	(Row144 #= 104 and Freq144 #= 1) or
	(Row144 #= 105 and Freq144 #= 1) or
	(Row144 #= 106 and Freq144 #= 1) or
	(Row144 #= 115 and Freq144 #= 1) or
	(Row144 #= 116 and Freq144 #= 1) or
	(Row144 #= 117 and Freq144 #= 1) or
	(Row144 #= 118 and Freq144 #= 1) or
	(Row144 #= 121 and Freq144 #= 1) or
	(Row144 #= 122 and Freq144 #= 1) or
	(Row144 #= 123 and Freq144 #= 1) or
	(Row144 #= 124 and Freq144 #= 1) or
	(Row144 #= 125 and Freq144 #= 1) or
	(Row144 #= 126 and Freq144 #= 1) or
	(Row144 #= 127 and Freq144 #= 1) or
	(Row144 #= 128 and Freq144 #= 1) or
	(Row144 #= 129 and Freq144 #= 1) or
	(Row144 #= 130 and Freq144 #= 1) or
	(Row144 #= 131 and Freq144 #= 1) or
	(Row144 #= 132 and Freq144 #= 1) or
	(Row144 #= 133 and Freq144 #= 1) or
	(Row144 #= 134 and Freq144 #= 1) or
	(Row144 #= 135 and Freq144 #= 1) or
	(Row144 #= 136 and Freq144 #= 1) or
	(Row144 #= 137 and Freq144 #= 1) or
	(Row144 #= 138 and Freq144 #= 1) or
	(Row144 #= 139 and Freq144 #= 1) or
	(Row144 #= 140 and Freq144 #= 1) or
	(Row144 #= 141 and Freq144 #= 1) or
	(Row144 #= 142 and Freq144 #= 1) or
	(Row144 #= 143 and Freq144 #= 1) or
	(Row144 #= 144 and Freq144 #= 1) or
	(Row144 #= 145 and Freq144 #= 1) or
	(Row144 #= 146 and Freq144 #= 1) or
	(Row144 #= 147 and Freq144 #= 1) or
	(Row144 #= 148 and Freq144 #= 1) or
	(Row144 #= 149 and Freq144 #= 1) or
	(Row144 #= 150 and Freq144 #= 1) or
	(Row144 #= 152 and Freq144 #= 1) or
	(Row144 #= 0 and Freq144 #= 0)), 

	((Row145 #= 49 and Freq145 #= 1) or
	(Row145 #= 77 and Freq145 #= 1) or
	(Row145 #= 79 and Freq145 #= 1) or
	(Row145 #= 81 and Freq145 #= 1) or
	(Row145 #= 82 and Freq145 #= 1) or
	(Row145 #= 84 and Freq145 #= 1) or
	(Row145 #= 85 and Freq145 #= 1) or
	(Row145 #= 86 and Freq145 #= 1) or
	(Row145 #= 87 and Freq145 #= 1) or
	(Row145 #= 90 and Freq145 #= 1) or
	(Row145 #= 91 and Freq145 #= 1) or
	(Row145 #= 92 and Freq145 #= 1) or
	(Row145 #= 93 and Freq145 #= 1) or
	(Row145 #= 94 and Freq145 #= 1) or
	(Row145 #= 95 and Freq145 #= 1) or
	(Row145 #= 97 and Freq145 #= 1) or
	(Row145 #= 98 and Freq145 #= 1) or
	(Row145 #= 99 and Freq145 #= 1) or
	(Row145 #= 100 and Freq145 #= 1) or
	(Row145 #= 101 and Freq145 #= 1) or
	(Row145 #= 102 and Freq145 #= 1) or
	(Row145 #= 104 and Freq145 #= 1) or
	(Row145 #= 105 and Freq145 #= 1) or
	(Row145 #= 106 and Freq145 #= 1) or
	(Row145 #= 115 and Freq145 #= 1) or
	(Row145 #= 116 and Freq145 #= 1) or
	(Row145 #= 117 and Freq145 #= 1) or
	(Row145 #= 118 and Freq145 #= 1) or
	(Row145 #= 120 and Freq145 #= 1) or
	(Row145 #= 121 and Freq145 #= 1) or
	(Row145 #= 124 and Freq145 #= 1) or
	(Row145 #= 125 and Freq145 #= 1) or
	(Row145 #= 126 and Freq145 #= 1) or
	(Row145 #= 127 and Freq145 #= 1) or
	(Row145 #= 128 and Freq145 #= 1) or
	(Row145 #= 129 and Freq145 #= 1) or
	(Row145 #= 130 and Freq145 #= 1) or
	(Row145 #= 131 and Freq145 #= 1) or
	(Row145 #= 132 and Freq145 #= 1) or
	(Row145 #= 133 and Freq145 #= 1) or
	(Row145 #= 134 and Freq145 #= 1) or
	(Row145 #= 135 and Freq145 #= 1) or
	(Row145 #= 136 and Freq145 #= 1) or
	(Row145 #= 137 and Freq145 #= 1) or
	(Row145 #= 138 and Freq145 #= 1) or
	(Row145 #= 139 and Freq145 #= 1) or
	(Row145 #= 140 and Freq145 #= 1) or
	(Row145 #= 141 and Freq145 #= 1) or
	(Row145 #= 142 and Freq145 #= 1) or
	(Row145 #= 143 and Freq145 #= 1) or
	(Row145 #= 144 and Freq145 #= 1) or
	(Row145 #= 145 and Freq145 #= 1) or
	(Row145 #= 146 and Freq145 #= 1) or
	(Row145 #= 147 and Freq145 #= 1) or
	(Row145 #= 149 and Freq145 #= 1) or
	(Row145 #= 159 and Freq145 #= 1) or
	(Row145 #= 0 and Freq145 #= 0)), 

	((Row146 #= 72 and Freq146 #= 1) or
	(Row146 #= 80 and Freq146 #= 1) or
	(Row146 #= 82 and Freq146 #= 1) or
	(Row146 #= 84 and Freq146 #= 1) or
	(Row146 #= 85 and Freq146 #= 1) or
	(Row146 #= 86 and Freq146 #= 1) or
	(Row146 #= 87 and Freq146 #= 1) or
	(Row146 #= 88 and Freq146 #= 1) or
	(Row146 #= 89 and Freq146 #= 1) or
	(Row146 #= 90 and Freq146 #= 1) or
	(Row146 #= 91 and Freq146 #= 1) or
	(Row146 #= 92 and Freq146 #= 1) or
	(Row146 #= 93 and Freq146 #= 1) or
	(Row146 #= 94 and Freq146 #= 1) or
	(Row146 #= 95 and Freq146 #= 1) or
	(Row146 #= 96 and Freq146 #= 1) or
	(Row146 #= 97 and Freq146 #= 1) or
	(Row146 #= 98 and Freq146 #= 1) or
	(Row146 #= 99 and Freq146 #= 1) or
	(Row146 #= 100 and Freq146 #= 1) or
	(Row146 #= 101 and Freq146 #= 1) or
	(Row146 #= 102 and Freq146 #= 1) or
	(Row146 #= 103 and Freq146 #= 1) or
	(Row146 #= 104 and Freq146 #= 1) or
	(Row146 #= 105 and Freq146 #= 1) or
	(Row146 #= 106 and Freq146 #= 1) or
	(Row146 #= 115 and Freq146 #= 1) or
	(Row146 #= 116 and Freq146 #= 1) or
	(Row146 #= 117 and Freq146 #= 1) or
	(Row146 #= 118 and Freq146 #= 1) or
	(Row146 #= 119 and Freq146 #= 1) or
	(Row146 #= 120 and Freq146 #= 1) or
	(Row146 #= 121 and Freq146 #= 1) or
	(Row146 #= 122 and Freq146 #= 1) or
	(Row146 #= 123 and Freq146 #= 1) or
	(Row146 #= 124 and Freq146 #= 1) or
	(Row146 #= 125 and Freq146 #= 1) or
	(Row146 #= 126 and Freq146 #= 1) or
	(Row146 #= 127 and Freq146 #= 1) or
	(Row146 #= 128 and Freq146 #= 1) or
	(Row146 #= 129 and Freq146 #= 1) or
	(Row146 #= 130 and Freq146 #= 1) or
	(Row146 #= 131 and Freq146 #= 1) or
	(Row146 #= 132 and Freq146 #= 1) or
	(Row146 #= 133 and Freq146 #= 1) or
	(Row146 #= 134 and Freq146 #= 1) or
	(Row146 #= 135 and Freq146 #= 1) or
	(Row146 #= 136 and Freq146 #= 1) or
	(Row146 #= 137 and Freq146 #= 1) or
	(Row146 #= 138 and Freq146 #= 1) or
	(Row146 #= 139 and Freq146 #= 1) or
	(Row146 #= 140 and Freq146 #= 1) or
	(Row146 #= 141 and Freq146 #= 1) or
	(Row146 #= 142 and Freq146 #= 1) or
	(Row146 #= 143 and Freq146 #= 1) or
	(Row146 #= 144 and Freq146 #= 1) or
	(Row146 #= 145 and Freq146 #= 1) or
	(Row146 #= 146 and Freq146 #= 1) or
	(Row146 #= 147 and Freq146 #= 1) or
	(Row146 #= 148 and Freq146 #= 1) or
	(Row146 #= 150 and Freq146 #= 1) or
	(Row146 #= 151 and Freq146 #= 1) or
	(Row146 #= 152 and Freq146 #= 1) or
	(Row146 #= 153 and Freq146 #= 1) or
	(Row146 #= 155 and Freq146 #= 1) or
	(Row146 #= 245 and Freq146 #= 1) or
	(Row146 #= 0 and Freq146 #= 0)), 

	((Row147 #= 73 and Freq147 #= 1) or
	(Row147 #= 80 and Freq147 #= 1) or
	(Row147 #= 82 and Freq147 #= 1) or
	(Row147 #= 84 and Freq147 #= 1) or
	(Row147 #= 86 and Freq147 #= 1) or
	(Row147 #= 87 and Freq147 #= 1) or
	(Row147 #= 89 and Freq147 #= 1) or
	(Row147 #= 90 and Freq147 #= 1) or
	(Row147 #= 91 and Freq147 #= 1) or
	(Row147 #= 92 and Freq147 #= 1) or
	(Row147 #= 93 and Freq147 #= 1) or
	(Row147 #= 94 and Freq147 #= 1) or
	(Row147 #= 95 and Freq147 #= 1) or
	(Row147 #= 96 and Freq147 #= 1) or
	(Row147 #= 97 and Freq147 #= 1) or
	(Row147 #= 98 and Freq147 #= 1) or
	(Row147 #= 99 and Freq147 #= 1) or
	(Row147 #= 100 and Freq147 #= 1) or
	(Row147 #= 101 and Freq147 #= 1) or
	(Row147 #= 102 and Freq147 #= 1) or
	(Row147 #= 103 and Freq147 #= 1) or
	(Row147 #= 104 and Freq147 #= 1) or
	(Row147 #= 105 and Freq147 #= 1) or
	(Row147 #= 106 and Freq147 #= 1) or
	(Row147 #= 115 and Freq147 #= 1) or
	(Row147 #= 116 and Freq147 #= 1) or
	(Row147 #= 117 and Freq147 #= 1) or
	(Row147 #= 118 and Freq147 #= 1) or
	(Row147 #= 119 and Freq147 #= 1) or
	(Row147 #= 120 and Freq147 #= 1) or
	(Row147 #= 121 and Freq147 #= 1) or
	(Row147 #= 122 and Freq147 #= 1) or
	(Row147 #= 123 and Freq147 #= 1) or
	(Row147 #= 124 and Freq147 #= 1) or
	(Row147 #= 125 and Freq147 #= 1) or
	(Row147 #= 126 and Freq147 #= 1) or
	(Row147 #= 127 and Freq147 #= 1) or
	(Row147 #= 128 and Freq147 #= 1) or
	(Row147 #= 129 and Freq147 #= 1) or
	(Row147 #= 130 and Freq147 #= 1) or
	(Row147 #= 131 and Freq147 #= 1) or
	(Row147 #= 132 and Freq147 #= 1) or
	(Row147 #= 133 and Freq147 #= 1) or
	(Row147 #= 134 and Freq147 #= 1) or
	(Row147 #= 135 and Freq147 #= 1) or
	(Row147 #= 136 and Freq147 #= 1) or
	(Row147 #= 137 and Freq147 #= 1) or
	(Row147 #= 138 and Freq147 #= 1) or
	(Row147 #= 139 and Freq147 #= 1) or
	(Row147 #= 140 and Freq147 #= 1) or
	(Row147 #= 141 and Freq147 #= 1) or
	(Row147 #= 142 and Freq147 #= 1) or
	(Row147 #= 143 and Freq147 #= 1) or
	(Row147 #= 144 and Freq147 #= 1) or
	(Row147 #= 145 and Freq147 #= 1) or
	(Row147 #= 146 and Freq147 #= 2) or
	(Row147 #= 147 and Freq147 #= 1) or
	(Row147 #= 148 and Freq147 #= 1) or
	(Row147 #= 149 and Freq147 #= 1) or
	(Row147 #= 150 and Freq147 #= 1) or
	(Row147 #= 153 and Freq147 #= 1) or
	(Row147 #= 154 and Freq147 #= 1) or
	(Row147 #= 0 and Freq147 #= 0)), 

	((Row148 #= 39 and Freq148 #= 1) or
	(Row148 #= 69 and Freq148 #= 1) or
	(Row148 #= 73 and Freq148 #= 1) or
	(Row148 #= 80 and Freq148 #= 1) or
	(Row148 #= 82 and Freq148 #= 1) or
	(Row148 #= 83 and Freq148 #= 1) or
	(Row148 #= 84 and Freq148 #= 1) or
	(Row148 #= 86 and Freq148 #= 1) or
	(Row148 #= 87 and Freq148 #= 1) or
	(Row148 #= 89 and Freq148 #= 1) or
	(Row148 #= 90 and Freq148 #= 1) or
	(Row148 #= 91 and Freq148 #= 1) or
	(Row148 #= 92 and Freq148 #= 1) or
	(Row148 #= 93 and Freq148 #= 1) or
	(Row148 #= 94 and Freq148 #= 1) or
	(Row148 #= 95 and Freq148 #= 1) or
	(Row148 #= 96 and Freq148 #= 1) or
	(Row148 #= 97 and Freq148 #= 1) or
	(Row148 #= 98 and Freq148 #= 1) or
	(Row148 #= 99 and Freq148 #= 1) or
	(Row148 #= 100 and Freq148 #= 1) or
	(Row148 #= 101 and Freq148 #= 1) or
	(Row148 #= 102 and Freq148 #= 1) or
	(Row148 #= 103 and Freq148 #= 1) or
	(Row148 #= 105 and Freq148 #= 1) or
	(Row148 #= 106 and Freq148 #= 1) or
	(Row148 #= 115 and Freq148 #= 1) or
	(Row148 #= 116 and Freq148 #= 1) or
	(Row148 #= 117 and Freq148 #= 1) or
	(Row148 #= 118 and Freq148 #= 1) or
	(Row148 #= 119 and Freq148 #= 1) or
	(Row148 #= 120 and Freq148 #= 1) or
	(Row148 #= 121 and Freq148 #= 1) or
	(Row148 #= 122 and Freq148 #= 1) or
	(Row148 #= 123 and Freq148 #= 1) or
	(Row148 #= 124 and Freq148 #= 1) or
	(Row148 #= 125 and Freq148 #= 1) or
	(Row148 #= 126 and Freq148 #= 1) or
	(Row148 #= 127 and Freq148 #= 1) or
	(Row148 #= 128 and Freq148 #= 1) or
	(Row148 #= 129 and Freq148 #= 1) or
	(Row148 #= 130 and Freq148 #= 1) or
	(Row148 #= 131 and Freq148 #= 1) or
	(Row148 #= 132 and Freq148 #= 1) or
	(Row148 #= 133 and Freq148 #= 1) or
	(Row148 #= 134 and Freq148 #= 1) or
	(Row148 #= 135 and Freq148 #= 1) or
	(Row148 #= 136 and Freq148 #= 1) or
	(Row148 #= 137 and Freq148 #= 1) or
	(Row148 #= 138 and Freq148 #= 1) or
	(Row148 #= 139 and Freq148 #= 1) or
	(Row148 #= 140 and Freq148 #= 1) or
	(Row148 #= 141 and Freq148 #= 1) or
	(Row148 #= 142 and Freq148 #= 1) or
	(Row148 #= 143 and Freq148 #= 1) or
	(Row148 #= 144 and Freq148 #= 1) or
	(Row148 #= 145 and Freq148 #= 2) or
	(Row148 #= 146 and Freq148 #= 2) or
	(Row148 #= 147 and Freq148 #= 2) or
	(Row148 #= 148 and Freq148 #= 1) or
	(Row148 #= 149 and Freq148 #= 1) or
	(Row148 #= 150 and Freq148 #= 1) or
	(Row148 #= 153 and Freq148 #= 1) or
	(Row148 #= 154 and Freq148 #= 1) or
	(Row148 #= 0 and Freq148 #= 0)), 

	((Row149 #= 79 and Freq149 #= 1) or
	(Row149 #= 81 and Freq149 #= 1) or
	(Row149 #= 83 and Freq149 #= 1) or
	(Row149 #= 84 and Freq149 #= 1) or
	(Row149 #= 86 and Freq149 #= 1) or
	(Row149 #= 87 and Freq149 #= 1) or
	(Row149 #= 89 and Freq149 #= 1) or
	(Row149 #= 90 and Freq149 #= 1) or
	(Row149 #= 91 and Freq149 #= 1) or
	(Row149 #= 92 and Freq149 #= 1) or
	(Row149 #= 93 and Freq149 #= 1) or
	(Row149 #= 94 and Freq149 #= 1) or
	(Row149 #= 95 and Freq149 #= 1) or
	(Row149 #= 96 and Freq149 #= 1) or
	(Row149 #= 97 and Freq149 #= 1) or
	(Row149 #= 98 and Freq149 #= 1) or
	(Row149 #= 99 and Freq149 #= 1) or
	(Row149 #= 100 and Freq149 #= 1) or
	(Row149 #= 101 and Freq149 #= 2) or
	(Row149 #= 102 and Freq149 #= 1) or
	(Row149 #= 103 and Freq149 #= 1) or
	(Row149 #= 104 and Freq149 #= 2) or
	(Row149 #= 105 and Freq149 #= 2) or
	(Row149 #= 106 and Freq149 #= 1) or
	(Row149 #= 115 and Freq149 #= 1) or
	(Row149 #= 116 and Freq149 #= 1) or
	(Row149 #= 117 and Freq149 #= 1) or
	(Row149 #= 118 and Freq149 #= 1) or
	(Row149 #= 119 and Freq149 #= 1) or
	(Row149 #= 120 and Freq149 #= 1) or
	(Row149 #= 121 and Freq149 #= 1) or
	(Row149 #= 122 and Freq149 #= 1) or
	(Row149 #= 123 and Freq149 #= 1) or
	(Row149 #= 124 and Freq149 #= 1) or
	(Row149 #= 125 and Freq149 #= 1) or
	(Row149 #= 126 and Freq149 #= 1) or
	(Row149 #= 127 and Freq149 #= 1) or
	(Row149 #= 128 and Freq149 #= 1) or
	(Row149 #= 129 and Freq149 #= 1) or
	(Row149 #= 130 and Freq149 #= 1) or
	(Row149 #= 131 and Freq149 #= 1) or
	(Row149 #= 132 and Freq149 #= 1) or
	(Row149 #= 133 and Freq149 #= 1) or
	(Row149 #= 134 and Freq149 #= 1) or
	(Row149 #= 135 and Freq149 #= 1) or
	(Row149 #= 136 and Freq149 #= 1) or
	(Row149 #= 137 and Freq149 #= 1) or
	(Row149 #= 138 and Freq149 #= 1) or
	(Row149 #= 139 and Freq149 #= 1) or
	(Row149 #= 140 and Freq149 #= 1) or
	(Row149 #= 141 and Freq149 #= 1) or
	(Row149 #= 142 and Freq149 #= 1) or
	(Row149 #= 143 and Freq149 #= 1) or
	(Row149 #= 144 and Freq149 #= 1) or
	(Row149 #= 145 and Freq149 #= 1) or
	(Row149 #= 146 and Freq149 #= 1) or
	(Row149 #= 147 and Freq149 #= 1) or
	(Row149 #= 148 and Freq149 #= 1) or
	(Row149 #= 149 and Freq149 #= 1) or
	(Row149 #= 150 and Freq149 #= 1) or
	(Row149 #= 151 and Freq149 #= 1) or
	(Row149 #= 152 and Freq149 #= 1) or
	(Row149 #= 153 and Freq149 #= 1) or
	(Row149 #= 0 and Freq149 #= 0)), 

	((Row150 #= 81 and Freq150 #= 1) or
	(Row150 #= 84 and Freq150 #= 1) or
	(Row150 #= 85 and Freq150 #= 1) or
	(Row150 #= 86 and Freq150 #= 1) or
	(Row150 #= 87 and Freq150 #= 1) or
	(Row150 #= 88 and Freq150 #= 1) or
	(Row150 #= 89 and Freq150 #= 1) or
	(Row150 #= 90 and Freq150 #= 1) or
	(Row150 #= 91 and Freq150 #= 1) or
	(Row150 #= 92 and Freq150 #= 1) or
	(Row150 #= 93 and Freq150 #= 1) or
	(Row150 #= 94 and Freq150 #= 1) or
	(Row150 #= 95 and Freq150 #= 1) or
	(Row150 #= 96 and Freq150 #= 1) or
	(Row150 #= 97 and Freq150 #= 1) or
	(Row150 #= 98 and Freq150 #= 1) or
	(Row150 #= 99 and Freq150 #= 1) or
	(Row150 #= 100 and Freq150 #= 1) or
	(Row150 #= 101 and Freq150 #= 1) or
	(Row150 #= 102 and Freq150 #= 1) or
	(Row150 #= 103 and Freq150 #= 1) or
	(Row150 #= 104 and Freq150 #= 2) or
	(Row150 #= 105 and Freq150 #= 1) or
	(Row150 #= 106 and Freq150 #= 1) or
	(Row150 #= 115 and Freq150 #= 1) or
	(Row150 #= 116 and Freq150 #= 2) or
	(Row150 #= 117 and Freq150 #= 2) or
	(Row150 #= 118 and Freq150 #= 1) or
	(Row150 #= 119 and Freq150 #= 1) or
	(Row150 #= 120 and Freq150 #= 1) or
	(Row150 #= 121 and Freq150 #= 1) or
	(Row150 #= 122 and Freq150 #= 1) or
	(Row150 #= 123 and Freq150 #= 1) or
	(Row150 #= 124 and Freq150 #= 1) or
	(Row150 #= 125 and Freq150 #= 1) or
	(Row150 #= 126 and Freq150 #= 1) or
	(Row150 #= 127 and Freq150 #= 1) or
	(Row150 #= 128 and Freq150 #= 1) or
	(Row150 #= 129 and Freq150 #= 1) or
	(Row150 #= 130 and Freq150 #= 1) or
	(Row150 #= 131 and Freq150 #= 1) or
	(Row150 #= 132 and Freq150 #= 1) or
	(Row150 #= 133 and Freq150 #= 1) or
	(Row150 #= 134 and Freq150 #= 1) or
	(Row150 #= 135 and Freq150 #= 1) or
	(Row150 #= 136 and Freq150 #= 1) or
	(Row150 #= 137 and Freq150 #= 1) or
	(Row150 #= 138 and Freq150 #= 1) or
	(Row150 #= 139 and Freq150 #= 1) or
	(Row150 #= 140 and Freq150 #= 1) or
	(Row150 #= 141 and Freq150 #= 1) or
	(Row150 #= 142 and Freq150 #= 1) or
	(Row150 #= 143 and Freq150 #= 1) or
	(Row150 #= 144 and Freq150 #= 1) or
	(Row150 #= 145 and Freq150 #= 1) or
	(Row150 #= 146 and Freq150 #= 1) or
	(Row150 #= 147 and Freq150 #= 1) or
	(Row150 #= 148 and Freq150 #= 1) or
	(Row150 #= 149 and Freq150 #= 1) or
	(Row150 #= 0 and Freq150 #= 0)), 

	((Row151 #= 80 and Freq151 #= 1) or
	(Row151 #= 82 and Freq151 #= 1) or
	(Row151 #= 83 and Freq151 #= 1) or
	(Row151 #= 84 and Freq151 #= 1) or
	(Row151 #= 85 and Freq151 #= 1) or
	(Row151 #= 87 and Freq151 #= 1) or
	(Row151 #= 88 and Freq151 #= 1) or
	(Row151 #= 89 and Freq151 #= 1) or
	(Row151 #= 90 and Freq151 #= 1) or
	(Row151 #= 91 and Freq151 #= 1) or
	(Row151 #= 92 and Freq151 #= 1) or
	(Row151 #= 93 and Freq151 #= 1) or
	(Row151 #= 94 and Freq151 #= 1) or
	(Row151 #= 95 and Freq151 #= 1) or
	(Row151 #= 96 and Freq151 #= 1) or
	(Row151 #= 97 and Freq151 #= 1) or
	(Row151 #= 98 and Freq151 #= 1) or
	(Row151 #= 99 and Freq151 #= 2) or
	(Row151 #= 100 and Freq151 #= 2) or
	(Row151 #= 101 and Freq151 #= 2) or
	(Row151 #= 102 and Freq151 #= 2) or
	(Row151 #= 103 and Freq151 #= 2) or
	(Row151 #= 104 and Freq151 #= 2) or
	(Row151 #= 105 and Freq151 #= 2) or
	(Row151 #= 106 and Freq151 #= 1) or
	(Row151 #= 115 and Freq151 #= 1) or
	(Row151 #= 116 and Freq151 #= 2) or
	(Row151 #= 117 and Freq151 #= 2) or
	(Row151 #= 118 and Freq151 #= 1) or
	(Row151 #= 119 and Freq151 #= 1) or
	(Row151 #= 120 and Freq151 #= 1) or
	(Row151 #= 121 and Freq151 #= 1) or
	(Row151 #= 122 and Freq151 #= 1) or
	(Row151 #= 123 and Freq151 #= 1) or
	(Row151 #= 124 and Freq151 #= 1) or
	(Row151 #= 125 and Freq151 #= 1) or
	(Row151 #= 126 and Freq151 #= 1) or
	(Row151 #= 127 and Freq151 #= 1) or
	(Row151 #= 128 and Freq151 #= 1) or
	(Row151 #= 129 and Freq151 #= 1) or
	(Row151 #= 130 and Freq151 #= 1) or
	(Row151 #= 131 and Freq151 #= 1) or
	(Row151 #= 132 and Freq151 #= 1) or
	(Row151 #= 133 and Freq151 #= 1) or
	(Row151 #= 134 and Freq151 #= 1) or
	(Row151 #= 135 and Freq151 #= 1) or
	(Row151 #= 136 and Freq151 #= 1) or
	(Row151 #= 137 and Freq151 #= 1) or
	(Row151 #= 138 and Freq151 #= 1) or
	(Row151 #= 139 and Freq151 #= 1) or
	(Row151 #= 140 and Freq151 #= 1) or
	(Row151 #= 141 and Freq151 #= 1) or
	(Row151 #= 142 and Freq151 #= 1) or
	(Row151 #= 143 and Freq151 #= 1) or
	(Row151 #= 144 and Freq151 #= 1) or
	(Row151 #= 145 and Freq151 #= 1) or
	(Row151 #= 146 and Freq151 #= 1) or
	(Row151 #= 147 and Freq151 #= 1) or
	(Row151 #= 151 and Freq151 #= 1) or
	(Row151 #= 152 and Freq151 #= 1) or
	(Row151 #= 0 and Freq151 #= 0)), 

	((Row152 #= 81 and Freq152 #= 1) or
	(Row152 #= 82 and Freq152 #= 1) or
	(Row152 #= 84 and Freq152 #= 1) or
	(Row152 #= 85 and Freq152 #= 1) or
	(Row152 #= 86 and Freq152 #= 1) or
	(Row152 #= 87 and Freq152 #= 1) or
	(Row152 #= 88 and Freq152 #= 1) or
	(Row152 #= 89 and Freq152 #= 1) or
	(Row152 #= 90 and Freq152 #= 1) or
	(Row152 #= 91 and Freq152 #= 1) or
	(Row152 #= 92 and Freq152 #= 1) or
	(Row152 #= 93 and Freq152 #= 1) or
	(Row152 #= 94 and Freq152 #= 1) or
	(Row152 #= 95 and Freq152 #= 1) or
	(Row152 #= 96 and Freq152 #= 1) or
	(Row152 #= 97 and Freq152 #= 1) or
	(Row152 #= 98 and Freq152 #= 1) or
	(Row152 #= 99 and Freq152 #= 2) or
	(Row152 #= 100 and Freq152 #= 2) or
	(Row152 #= 101 and Freq152 #= 2) or
	(Row152 #= 102 and Freq152 #= 1) or
	(Row152 #= 103 and Freq152 #= 2) or
	(Row152 #= 104 and Freq152 #= 2) or
	(Row152 #= 105 and Freq152 #= 2) or
	(Row152 #= 106 and Freq152 #= 2) or
	(Row152 #= 115 and Freq152 #= 2) or
	(Row152 #= 116 and Freq152 #= 2) or
	(Row152 #= 117 and Freq152 #= 1) or
	(Row152 #= 118 and Freq152 #= 1) or
	(Row152 #= 119 and Freq152 #= 1) or
	(Row152 #= 120 and Freq152 #= 2) or
	(Row152 #= 121 and Freq152 #= 2) or
	(Row152 #= 122 and Freq152 #= 1) or
	(Row152 #= 123 and Freq152 #= 2) or
	(Row152 #= 124 and Freq152 #= 1) or
	(Row152 #= 125 and Freq152 #= 1) or
	(Row152 #= 126 and Freq152 #= 1) or
	(Row152 #= 127 and Freq152 #= 1) or
	(Row152 #= 128 and Freq152 #= 1) or
	(Row152 #= 129 and Freq152 #= 1) or
	(Row152 #= 130 and Freq152 #= 1) or
	(Row152 #= 131 and Freq152 #= 1) or
	(Row152 #= 132 and Freq152 #= 1) or
	(Row152 #= 133 and Freq152 #= 1) or
	(Row152 #= 134 and Freq152 #= 1) or
	(Row152 #= 135 and Freq152 #= 1) or
	(Row152 #= 136 and Freq152 #= 1) or
	(Row152 #= 137 and Freq152 #= 1) or
	(Row152 #= 138 and Freq152 #= 1) or
	(Row152 #= 139 and Freq152 #= 1) or
	(Row152 #= 140 and Freq152 #= 1) or
	(Row152 #= 141 and Freq152 #= 1) or
	(Row152 #= 142 and Freq152 #= 1) or
	(Row152 #= 143 and Freq152 #= 1) or
	(Row152 #= 144 and Freq152 #= 1) or
	(Row152 #= 145 and Freq152 #= 1) or
	(Row152 #= 146 and Freq152 #= 1) or
	(Row152 #= 147 and Freq152 #= 1) or
	(Row152 #= 0 and Freq152 #= 0)), 

	((Row153 #= 80 and Freq153 #= 1) or
	(Row153 #= 82 and Freq153 #= 1) or
	(Row153 #= 84 and Freq153 #= 1) or
	(Row153 #= 86 and Freq153 #= 1) or
	(Row153 #= 87 and Freq153 #= 1) or
	(Row153 #= 89 and Freq153 #= 1) or
	(Row153 #= 90 and Freq153 #= 1) or
	(Row153 #= 91 and Freq153 #= 1) or
	(Row153 #= 92 and Freq153 #= 1) or
	(Row153 #= 93 and Freq153 #= 1) or
	(Row153 #= 94 and Freq153 #= 1) or
	(Row153 #= 95 and Freq153 #= 1) or
	(Row153 #= 96 and Freq153 #= 1) or
	(Row153 #= 97 and Freq153 #= 2) or
	(Row153 #= 98 and Freq153 #= 1) or
	(Row153 #= 99 and Freq153 #= 2) or
	(Row153 #= 100 and Freq153 #= 2) or
	(Row153 #= 101 and Freq153 #= 2) or
	(Row153 #= 102 and Freq153 #= 2) or
	(Row153 #= 103 and Freq153 #= 2) or
	(Row153 #= 104 and Freq153 #= 2) or
	(Row153 #= 105 and Freq153 #= 2) or
	(Row153 #= 106 and Freq153 #= 2) or
	(Row153 #= 115 and Freq153 #= 2) or
	(Row153 #= 116 and Freq153 #= 2) or
	(Row153 #= 117 and Freq153 #= 2) or
	(Row153 #= 118 and Freq153 #= 2) or
	(Row153 #= 119 and Freq153 #= 1) or
	(Row153 #= 120 and Freq153 #= 2) or
	(Row153 #= 121 and Freq153 #= 2) or
	(Row153 #= 122 and Freq153 #= 1) or
	(Row153 #= 123 and Freq153 #= 2) or
	(Row153 #= 124 and Freq153 #= 1) or
	(Row153 #= 125 and Freq153 #= 2) or
	(Row153 #= 126 and Freq153 #= 1) or
	(Row153 #= 127 and Freq153 #= 1) or
	(Row153 #= 128 and Freq153 #= 1) or
	(Row153 #= 129 and Freq153 #= 2) or
	(Row153 #= 130 and Freq153 #= 1) or
	(Row153 #= 131 and Freq153 #= 1) or
	(Row153 #= 132 and Freq153 #= 1) or
	(Row153 #= 133 and Freq153 #= 1) or
	(Row153 #= 134 and Freq153 #= 1) or
	(Row153 #= 135 and Freq153 #= 1) or
	(Row153 #= 136 and Freq153 #= 1) or
	(Row153 #= 137 and Freq153 #= 1) or
	(Row153 #= 138 and Freq153 #= 1) or
	(Row153 #= 140 and Freq153 #= 1) or
	(Row153 #= 141 and Freq153 #= 1) or
	(Row153 #= 142 and Freq153 #= 1) or
	(Row153 #= 144 and Freq153 #= 1) or
	(Row153 #= 145 and Freq153 #= 1) or
	(Row153 #= 146 and Freq153 #= 1) or
	(Row153 #= 147 and Freq153 #= 1) or
	(Row153 #= 0 and Freq153 #= 0)), 

	((Row154 #= 80 and Freq154 #= 1) or
	(Row154 #= 82 and Freq154 #= 1) or
	(Row154 #= 83 and Freq154 #= 1) or
	(Row154 #= 84 and Freq154 #= 1) or
	(Row154 #= 85 and Freq154 #= 1) or
	(Row154 #= 86 and Freq154 #= 1) or
	(Row154 #= 87 and Freq154 #= 1) or
	(Row154 #= 88 and Freq154 #= 1) or
	(Row154 #= 89 and Freq154 #= 1) or
	(Row154 #= 90 and Freq154 #= 1) or
	(Row154 #= 91 and Freq154 #= 1) or
	(Row154 #= 92 and Freq154 #= 1) or
	(Row154 #= 93 and Freq154 #= 1) or
	(Row154 #= 94 and Freq154 #= 1) or
	(Row154 #= 95 and Freq154 #= 1) or
	(Row154 #= 96 and Freq154 #= 2) or
	(Row154 #= 97 and Freq154 #= 2) or
	(Row154 #= 98 and Freq154 #= 1) or
	(Row154 #= 99 and Freq154 #= 2) or
	(Row154 #= 100 and Freq154 #= 2) or
	(Row154 #= 101 and Freq154 #= 2) or
	(Row154 #= 102 and Freq154 #= 2) or
	(Row154 #= 103 and Freq154 #= 2) or
	(Row154 #= 104 and Freq154 #= 2) or
	(Row154 #= 105 and Freq154 #= 3) or
	(Row154 #= 106 and Freq154 #= 2) or
	(Row154 #= 115 and Freq154 #= 2) or
	(Row154 #= 116 and Freq154 #= 2) or
	(Row154 #= 117 and Freq154 #= 2) or
	(Row154 #= 118 and Freq154 #= 2) or
	(Row154 #= 119 and Freq154 #= 2) or
	(Row154 #= 120 and Freq154 #= 1) or
	(Row154 #= 121 and Freq154 #= 2) or
	(Row154 #= 122 and Freq154 #= 2) or
	(Row154 #= 123 and Freq154 #= 1) or
	(Row154 #= 124 and Freq154 #= 2) or
	(Row154 #= 125 and Freq154 #= 1) or
	(Row154 #= 126 and Freq154 #= 1) or
	(Row154 #= 127 and Freq154 #= 1) or
	(Row154 #= 128 and Freq154 #= 1) or
	(Row154 #= 129 and Freq154 #= 1) or
	(Row154 #= 130 and Freq154 #= 1) or
	(Row154 #= 131 and Freq154 #= 1) or
	(Row154 #= 132 and Freq154 #= 1) or
	(Row154 #= 133 and Freq154 #= 1) or
	(Row154 #= 134 and Freq154 #= 1) or
	(Row154 #= 135 and Freq154 #= 1) or
	(Row154 #= 136 and Freq154 #= 1) or
	(Row154 #= 137 and Freq154 #= 1) or
	(Row154 #= 138 and Freq154 #= 1) or
	(Row154 #= 139 and Freq154 #= 1) or
	(Row154 #= 140 and Freq154 #= 1) or
	(Row154 #= 141 and Freq154 #= 1) or
	(Row154 #= 142 and Freq154 #= 1) or
	(Row154 #= 143 and Freq154 #= 1) or
	(Row154 #= 144 and Freq154 #= 1) or
	(Row154 #= 0 and Freq154 #= 0)), 

	((Row155 #= 80 and Freq155 #= 1) or
	(Row155 #= 82 and Freq155 #= 1) or
	(Row155 #= 83 and Freq155 #= 1) or
	(Row155 #= 84 and Freq155 #= 1) or
	(Row155 #= 87 and Freq155 #= 1) or
	(Row155 #= 89 and Freq155 #= 1) or
	(Row155 #= 90 and Freq155 #= 1) or
	(Row155 #= 91 and Freq155 #= 1) or
	(Row155 #= 92 and Freq155 #= 1) or
	(Row155 #= 93 and Freq155 #= 1) or
	(Row155 #= 94 and Freq155 #= 1) or
	(Row155 #= 95 and Freq155 #= 1) or
	(Row155 #= 96 and Freq155 #= 2) or
	(Row155 #= 97 and Freq155 #= 2) or
	(Row155 #= 98 and Freq155 #= 2) or
	(Row155 #= 99 and Freq155 #= 2) or
	(Row155 #= 100 and Freq155 #= 2) or
	(Row155 #= 101 and Freq155 #= 2) or
	(Row155 #= 102 and Freq155 #= 2) or
	(Row155 #= 103 and Freq155 #= 3) or
	(Row155 #= 104 and Freq155 #= 4) or
	(Row155 #= 105 and Freq155 #= 3) or
	(Row155 #= 106 and Freq155 #= 3) or
	(Row155 #= 115 and Freq155 #= 2) or
	(Row155 #= 116 and Freq155 #= 3) or
	(Row155 #= 117 and Freq155 #= 2) or
	(Row155 #= 118 and Freq155 #= 2) or
	(Row155 #= 119 and Freq155 #= 2) or
	(Row155 #= 120 and Freq155 #= 2) or
	(Row155 #= 121 and Freq155 #= 2) or
	(Row155 #= 122 and Freq155 #= 2) or
	(Row155 #= 123 and Freq155 #= 2) or
	(Row155 #= 124 and Freq155 #= 2) or
	(Row155 #= 125 and Freq155 #= 1) or
	(Row155 #= 126 and Freq155 #= 2) or
	(Row155 #= 127 and Freq155 #= 1) or
	(Row155 #= 128 and Freq155 #= 1) or
	(Row155 #= 129 and Freq155 #= 1) or
	(Row155 #= 130 and Freq155 #= 1) or
	(Row155 #= 131 and Freq155 #= 1) or
	(Row155 #= 132 and Freq155 #= 1) or
	(Row155 #= 133 and Freq155 #= 1) or
	(Row155 #= 134 and Freq155 #= 1) or
	(Row155 #= 135 and Freq155 #= 1) or
	(Row155 #= 136 and Freq155 #= 1) or
	(Row155 #= 137 and Freq155 #= 1) or
	(Row155 #= 138 and Freq155 #= 1) or
	(Row155 #= 139 and Freq155 #= 1) or
	(Row155 #= 140 and Freq155 #= 1) or
	(Row155 #= 141 and Freq155 #= 1) or
	(Row155 #= 143 and Freq155 #= 1) or
	(Row155 #= 146 and Freq155 #= 1) or
	(Row155 #= 0 and Freq155 #= 0)), 

	((Row156 #= 82 and Freq156 #= 1) or
	(Row156 #= 83 and Freq156 #= 1) or
	(Row156 #= 85 and Freq156 #= 1) or
	(Row156 #= 86 and Freq156 #= 1) or
	(Row156 #= 87 and Freq156 #= 1) or
	(Row156 #= 88 and Freq156 #= 1) or
	(Row156 #= 89 and Freq156 #= 1) or
	(Row156 #= 90 and Freq156 #= 1) or
	(Row156 #= 91 and Freq156 #= 1) or
	(Row156 #= 92 and Freq156 #= 1) or
	(Row156 #= 93 and Freq156 #= 1) or
	(Row156 #= 94 and Freq156 #= 1) or
	(Row156 #= 95 and Freq156 #= 1) or
	(Row156 #= 96 and Freq156 #= 1) or
	(Row156 #= 97 and Freq156 #= 2) or
	(Row156 #= 98 and Freq156 #= 2) or
	(Row156 #= 99 and Freq156 #= 2) or
	(Row156 #= 100 and Freq156 #= 2) or
	(Row156 #= 101 and Freq156 #= 2) or
	(Row156 #= 102 and Freq156 #= 3) or
	(Row156 #= 103 and Freq156 #= 3) or
	(Row156 #= 104 and Freq156 #= 4) or
	(Row156 #= 105 and Freq156 #= 4) or
	(Row156 #= 106 and Freq156 #= 3) or
	(Row156 #= 115 and Freq156 #= 3) or
	(Row156 #= 116 and Freq156 #= 3) or
	(Row156 #= 117 and Freq156 #= 2) or
	(Row156 #= 118 and Freq156 #= 2) or
	(Row156 #= 119 and Freq156 #= 2) or
	(Row156 #= 120 and Freq156 #= 2) or
	(Row156 #= 121 and Freq156 #= 2) or
	(Row156 #= 122 and Freq156 #= 2) or
	(Row156 #= 123 and Freq156 #= 1) or
	(Row156 #= 124 and Freq156 #= 2) or
	(Row156 #= 125 and Freq156 #= 2) or
	(Row156 #= 126 and Freq156 #= 2) or
	(Row156 #= 127 and Freq156 #= 1) or
	(Row156 #= 128 and Freq156 #= 1) or
	(Row156 #= 129 and Freq156 #= 1) or
	(Row156 #= 130 and Freq156 #= 1) or
	(Row156 #= 131 and Freq156 #= 1) or
	(Row156 #= 132 and Freq156 #= 1) or
	(Row156 #= 133 and Freq156 #= 1) or
	(Row156 #= 134 and Freq156 #= 1) or
	(Row156 #= 135 and Freq156 #= 1) or
	(Row156 #= 136 and Freq156 #= 1) or
	(Row156 #= 137 and Freq156 #= 1) or
	(Row156 #= 138 and Freq156 #= 1) or
	(Row156 #= 139 and Freq156 #= 1) or
	(Row156 #= 140 and Freq156 #= 1) or
	(Row156 #= 141 and Freq156 #= 1) or
	(Row156 #= 144 and Freq156 #= 1) or
	(Row156 #= 146 and Freq156 #= 1) or
	(Row156 #= 0 and Freq156 #= 0)), 

	((Row157 #= 81 and Freq157 #= 1) or
	(Row157 #= 82 and Freq157 #= 1) or
	(Row157 #= 83 and Freq157 #= 1) or
	(Row157 #= 84 and Freq157 #= 1) or
	(Row157 #= 85 and Freq157 #= 1) or
	(Row157 #= 86 and Freq157 #= 1) or
	(Row157 #= 87 and Freq157 #= 1) or
	(Row157 #= 88 and Freq157 #= 1) or
	(Row157 #= 89 and Freq157 #= 1) or
	(Row157 #= 90 and Freq157 #= 1) or
	(Row157 #= 91 and Freq157 #= 1) or
	(Row157 #= 92 and Freq157 #= 1) or
	(Row157 #= 93 and Freq157 #= 2) or
	(Row157 #= 94 and Freq157 #= 2) or
	(Row157 #= 95 and Freq157 #= 1) or
	(Row157 #= 96 and Freq157 #= 2) or
	(Row157 #= 97 and Freq157 #= 2) or
	(Row157 #= 98 and Freq157 #= 2) or
	(Row157 #= 99 and Freq157 #= 2) or
	(Row157 #= 100 and Freq157 #= 3) or
	(Row157 #= 101 and Freq157 #= 3) or
	(Row157 #= 102 and Freq157 #= 3) or
	(Row157 #= 103 and Freq157 #= 3) or
	(Row157 #= 104 and Freq157 #= 4) or
	(Row157 #= 105 and Freq157 #= 4) or
	(Row157 #= 106 and Freq157 #= 4) or
	(Row157 #= 115 and Freq157 #= 4) or
	(Row157 #= 116 and Freq157 #= 4) or
	(Row157 #= 117 and Freq157 #= 3) or
	(Row157 #= 118 and Freq157 #= 3) or
	(Row157 #= 119 and Freq157 #= 2) or
	(Row157 #= 120 and Freq157 #= 3) or
	(Row157 #= 121 and Freq157 #= 3) or
	(Row157 #= 122 and Freq157 #= 2) or
	(Row157 #= 123 and Freq157 #= 2) or
	(Row157 #= 124 and Freq157 #= 2) or
	(Row157 #= 125 and Freq157 #= 2) or
	(Row157 #= 126 and Freq157 #= 1) or
	(Row157 #= 127 and Freq157 #= 1) or
	(Row157 #= 128 and Freq157 #= 2) or
	(Row157 #= 129 and Freq157 #= 1) or
	(Row157 #= 130 and Freq157 #= 1) or
	(Row157 #= 131 and Freq157 #= 1) or
	(Row157 #= 132 and Freq157 #= 1) or
	(Row157 #= 133 and Freq157 #= 1) or
	(Row157 #= 134 and Freq157 #= 1) or
	(Row157 #= 135 and Freq157 #= 1) or
	(Row157 #= 136 and Freq157 #= 1) or
	(Row157 #= 137 and Freq157 #= 1) or
	(Row157 #= 138 and Freq157 #= 1) or
	(Row157 #= 140 and Freq157 #= 1) or
	(Row157 #= 141 and Freq157 #= 1) or
	(Row157 #= 142 and Freq157 #= 1) or
	(Row157 #= 143 and Freq157 #= 1) or
	(Row157 #= 146 and Freq157 #= 1) or
	(Row157 #= 0 and Freq157 #= 0)), 

	((Row158 #= 77 and Freq158 #= 1) or
	(Row158 #= 81 and Freq158 #= 1) or
	(Row158 #= 82 and Freq158 #= 1) or
	(Row158 #= 83 and Freq158 #= 1) or
	(Row158 #= 84 and Freq158 #= 1) or
	(Row158 #= 85 and Freq158 #= 1) or
	(Row158 #= 86 and Freq158 #= 1) or
	(Row158 #= 87 and Freq158 #= 1) or
	(Row158 #= 88 and Freq158 #= 1) or
	(Row158 #= 89 and Freq158 #= 1) or
	(Row158 #= 90 and Freq158 #= 1) or
	(Row158 #= 91 and Freq158 #= 1) or
	(Row158 #= 92 and Freq158 #= 1) or
	(Row158 #= 93 and Freq158 #= 1) or
	(Row158 #= 94 and Freq158 #= 2) or
	(Row158 #= 95 and Freq158 #= 1) or
	(Row158 #= 96 and Freq158 #= 2) or
	(Row158 #= 97 and Freq158 #= 2) or
	(Row158 #= 98 and Freq158 #= 2) or
	(Row158 #= 99 and Freq158 #= 3) or
	(Row158 #= 100 and Freq158 #= 2) or
	(Row158 #= 101 and Freq158 #= 3) or
	(Row158 #= 102 and Freq158 #= 3) or
	(Row158 #= 103 and Freq158 #= 4) or
	(Row158 #= 104 and Freq158 #= 5) or
	(Row158 #= 105 and Freq158 #= 5) or
	(Row158 #= 106 and Freq158 #= 5) or
	(Row158 #= 115 and Freq158 #= 4) or
	(Row158 #= 116 and Freq158 #= 4) or
	(Row158 #= 117 and Freq158 #= 4) or
	(Row158 #= 118 and Freq158 #= 3) or
	(Row158 #= 119 and Freq158 #= 3) or
	(Row158 #= 120 and Freq158 #= 3) or
	(Row158 #= 121 and Freq158 #= 3) or
	(Row158 #= 122 and Freq158 #= 2) or
	(Row158 #= 123 and Freq158 #= 2) or
	(Row158 #= 124 and Freq158 #= 2) or
	(Row158 #= 125 and Freq158 #= 2) or
	(Row158 #= 126 and Freq158 #= 2) or
	(Row158 #= 127 and Freq158 #= 2) or
	(Row158 #= 128 and Freq158 #= 2) or
	(Row158 #= 129 and Freq158 #= 2) or
	(Row158 #= 130 and Freq158 #= 1) or
	(Row158 #= 131 and Freq158 #= 2) or
	(Row158 #= 132 and Freq158 #= 1) or
	(Row158 #= 133 and Freq158 #= 1) or
	(Row158 #= 134 and Freq158 #= 1) or
	(Row158 #= 135 and Freq158 #= 1) or
	(Row158 #= 136 and Freq158 #= 1) or
	(Row158 #= 137 and Freq158 #= 1) or
	(Row158 #= 138 and Freq158 #= 1) or
	(Row158 #= 139 and Freq158 #= 1) or
	(Row158 #= 140 and Freq158 #= 1) or
	(Row158 #= 141 and Freq158 #= 1) or
	(Row158 #= 142 and Freq158 #= 1) or
	(Row158 #= 144 and Freq158 #= 1) or
	(Row158 #= 0 and Freq158 #= 0)), 

	((Row159 #= 82 and Freq159 #= 1) or
	(Row159 #= 83 and Freq159 #= 1) or
	(Row159 #= 84 and Freq159 #= 1) or
	(Row159 #= 85 and Freq159 #= 1) or
	(Row159 #= 86 and Freq159 #= 1) or
	(Row159 #= 87 and Freq159 #= 1) or
	(Row159 #= 89 and Freq159 #= 1) or
	(Row159 #= 90 and Freq159 #= 1) or
	(Row159 #= 91 and Freq159 #= 1) or
	(Row159 #= 92 and Freq159 #= 1) or
	(Row159 #= 93 and Freq159 #= 1) or
	(Row159 #= 94 and Freq159 #= 1) or
	(Row159 #= 95 and Freq159 #= 2) or
	(Row159 #= 96 and Freq159 #= 2) or
	(Row159 #= 97 and Freq159 #= 3) or
	(Row159 #= 98 and Freq159 #= 2) or
	(Row159 #= 99 and Freq159 #= 3) or
	(Row159 #= 100 and Freq159 #= 3) or
	(Row159 #= 101 and Freq159 #= 3) or
	(Row159 #= 102 and Freq159 #= 3) or
	(Row159 #= 103 and Freq159 #= 4) or
	(Row159 #= 104 and Freq159 #= 5) or
	(Row159 #= 105 and Freq159 #= 5) or
	(Row159 #= 106 and Freq159 #= 6) or
	(Row159 #= 115 and Freq159 #= 4) or
	(Row159 #= 116 and Freq159 #= 5) or
	(Row159 #= 117 and Freq159 #= 5) or
	(Row159 #= 118 and Freq159 #= 4) or
	(Row159 #= 119 and Freq159 #= 3) or
	(Row159 #= 120 and Freq159 #= 3) or
	(Row159 #= 121 and Freq159 #= 3) or
	(Row159 #= 122 and Freq159 #= 3) or
	(Row159 #= 123 and Freq159 #= 2) or
	(Row159 #= 124 and Freq159 #= 2) or
	(Row159 #= 125 and Freq159 #= 2) or
	(Row159 #= 126 and Freq159 #= 2) or
	(Row159 #= 127 and Freq159 #= 2) or
	(Row159 #= 128 and Freq159 #= 2) or
	(Row159 #= 129 and Freq159 #= 1) or
	(Row159 #= 130 and Freq159 #= 1) or
	(Row159 #= 131 and Freq159 #= 1) or
	(Row159 #= 132 and Freq159 #= 1) or
	(Row159 #= 133 and Freq159 #= 1) or
	(Row159 #= 134 and Freq159 #= 1) or
	(Row159 #= 135 and Freq159 #= 1) or
	(Row159 #= 136 and Freq159 #= 1) or
	(Row159 #= 137 and Freq159 #= 1) or
	(Row159 #= 138 and Freq159 #= 1) or
	(Row159 #= 139 and Freq159 #= 1) or
	(Row159 #= 140 and Freq159 #= 1) or
	(Row159 #= 141 and Freq159 #= 1) or
	(Row159 #= 142 and Freq159 #= 1) or
	(Row159 #= 0 and Freq159 #= 0)), 

	((Row160 #= 0 and Freq160 #= 0)), 

	((Row161 #= 0 and Freq161 #= 0)), 

	((Row162 #= 0 and Freq162 #= 0)), 

	((Row163 #= 0 and Freq163 #= 0)), 

	((Row164 #= 0 and Freq164 #= 0)), 

	((Row165 #= 0 and Freq165 #= 0)), 

	((Row166 #= 82 and Freq166 #= 1) or
	(Row166 #= 89 and Freq166 #= 1) or
	(Row166 #= 90 and Freq166 #= 1) or
	(Row166 #= 91 and Freq166 #= 1) or
	(Row166 #= 92 and Freq166 #= 1) or
	(Row166 #= 93 and Freq166 #= 1) or
	(Row166 #= 94 and Freq166 #= 1) or
	(Row166 #= 95 and Freq166 #= 1) or
	(Row166 #= 96 and Freq166 #= 2) or
	(Row166 #= 97 and Freq166 #= 2) or
	(Row166 #= 98 and Freq166 #= 2) or
	(Row166 #= 99 and Freq166 #= 2) or
	(Row166 #= 100 and Freq166 #= 3) or
	(Row166 #= 101 and Freq166 #= 3) or
	(Row166 #= 102 and Freq166 #= 3) or
	(Row166 #= 103 and Freq166 #= 3) or
	(Row166 #= 104 and Freq166 #= 4) or
	(Row166 #= 105 and Freq166 #= 5) or
	(Row166 #= 106 and Freq166 #= 5) or
	(Row166 #= 115 and Freq166 #= 6) or
	(Row166 #= 116 and Freq166 #= 5) or
	(Row166 #= 117 and Freq166 #= 5) or
	(Row166 #= 118 and Freq166 #= 4) or
	(Row166 #= 119 and Freq166 #= 3) or
	(Row166 #= 120 and Freq166 #= 3) or
	(Row166 #= 121 and Freq166 #= 3) or
	(Row166 #= 122 and Freq166 #= 2) or
	(Row166 #= 123 and Freq166 #= 2) or
	(Row166 #= 124 and Freq166 #= 1) or
	(Row166 #= 125 and Freq166 #= 2) or
	(Row166 #= 126 and Freq166 #= 1) or
	(Row166 #= 127 and Freq166 #= 1) or
	(Row166 #= 128 and Freq166 #= 1) or
	(Row166 #= 129 and Freq166 #= 1) or
	(Row166 #= 130 and Freq166 #= 1) or
	(Row166 #= 131 and Freq166 #= 1) or
	(Row166 #= 132 and Freq166 #= 1) or
	(Row166 #= 134 and Freq166 #= 1) or
	(Row166 #= 137 and Freq166 #= 1) or
	(Row166 #= 138 and Freq166 #= 1) or
	(Row166 #= 0 and Freq166 #= 0)), 

	((Row167 #= 81 and Freq167 #= 1) or
	(Row167 #= 82 and Freq167 #= 1) or
	(Row167 #= 85 and Freq167 #= 1) or
	(Row167 #= 86 and Freq167 #= 1) or
	(Row167 #= 87 and Freq167 #= 1) or
	(Row167 #= 88 and Freq167 #= 1) or
	(Row167 #= 89 and Freq167 #= 1) or
	(Row167 #= 90 and Freq167 #= 1) or
	(Row167 #= 91 and Freq167 #= 1) or
	(Row167 #= 92 and Freq167 #= 1) or
	(Row167 #= 93 and Freq167 #= 1) or
	(Row167 #= 94 and Freq167 #= 1) or
	(Row167 #= 95 and Freq167 #= 1) or
	(Row167 #= 96 and Freq167 #= 2) or
	(Row167 #= 97 and Freq167 #= 2) or
	(Row167 #= 98 and Freq167 #= 2) or
	(Row167 #= 99 and Freq167 #= 2) or
	(Row167 #= 100 and Freq167 #= 3) or
	(Row167 #= 101 and Freq167 #= 3) or
	(Row167 #= 102 and Freq167 #= 3) or
	(Row167 #= 103 and Freq167 #= 4) or
	(Row167 #= 104 and Freq167 #= 5) or
	(Row167 #= 105 and Freq167 #= 6) or
	(Row167 #= 106 and Freq167 #= 5) or
	(Row167 #= 115 and Freq167 #= 6) or
	(Row167 #= 116 and Freq167 #= 5) or
	(Row167 #= 117 and Freq167 #= 5) or
	(Row167 #= 118 and Freq167 #= 3) or
	(Row167 #= 119 and Freq167 #= 3) or
	(Row167 #= 120 and Freq167 #= 3) or
	(Row167 #= 121 and Freq167 #= 3) or
	(Row167 #= 122 and Freq167 #= 2) or
	(Row167 #= 123 and Freq167 #= 2) or
	(Row167 #= 124 and Freq167 #= 1) or
	(Row167 #= 125 and Freq167 #= 2) or
	(Row167 #= 126 and Freq167 #= 2) or
	(Row167 #= 127 and Freq167 #= 1) or
	(Row167 #= 128 and Freq167 #= 1) or
	(Row167 #= 129 and Freq167 #= 1) or
	(Row167 #= 130 and Freq167 #= 1) or
	(Row167 #= 131 and Freq167 #= 1) or
	(Row167 #= 132 and Freq167 #= 1) or
	(Row167 #= 133 and Freq167 #= 1) or
	(Row167 #= 136 and Freq167 #= 1) or
	(Row167 #= 137 and Freq167 #= 1) or
	(Row167 #= 138 and Freq167 #= 1) or
	(Row167 #= 0 and Freq167 #= 0)), 

	((Row168 #= 82 and Freq168 #= 1) or
	(Row168 #= 87 and Freq168 #= 1) or
	(Row168 #= 88 and Freq168 #= 1) or
	(Row168 #= 89 and Freq168 #= 1) or
	(Row168 #= 90 and Freq168 #= 1) or
	(Row168 #= 91 and Freq168 #= 1) or
	(Row168 #= 92 and Freq168 #= 1) or
	(Row168 #= 93 and Freq168 #= 1) or
	(Row168 #= 94 and Freq168 #= 1) or
	(Row168 #= 95 and Freq168 #= 1) or
	(Row168 #= 96 and Freq168 #= 1) or
	(Row168 #= 97 and Freq168 #= 2) or
	(Row168 #= 98 and Freq168 #= 2) or
	(Row168 #= 99 and Freq168 #= 2) or
	(Row168 #= 100 and Freq168 #= 2) or
	(Row168 #= 101 and Freq168 #= 2) or
	(Row168 #= 102 and Freq168 #= 2) or
	(Row168 #= 103 and Freq168 #= 3) or
	(Row168 #= 104 and Freq168 #= 4) or
	(Row168 #= 105 and Freq168 #= 4) or
	(Row168 #= 106 and Freq168 #= 5) or
	(Row168 #= 115 and Freq168 #= 4) or
	(Row168 #= 116 and Freq168 #= 4) or
	(Row168 #= 117 and Freq168 #= 3) or
	(Row168 #= 118 and Freq168 #= 3) or
	(Row168 #= 119 and Freq168 #= 3) or
	(Row168 #= 120 and Freq168 #= 3) or
	(Row168 #= 121 and Freq168 #= 2) or
	(Row168 #= 122 and Freq168 #= 2) or
	(Row168 #= 123 and Freq168 #= 2) or
	(Row168 #= 124 and Freq168 #= 2) or
	(Row168 #= 125 and Freq168 #= 1) or
	(Row168 #= 126 and Freq168 #= 1) or
	(Row168 #= 127 and Freq168 #= 1) or
	(Row168 #= 128 and Freq168 #= 1) or
	(Row168 #= 129 and Freq168 #= 1) or
	(Row168 #= 130 and Freq168 #= 1) or
	(Row168 #= 131 and Freq168 #= 1) or
	(Row168 #= 132 and Freq168 #= 1) or
	(Row168 #= 133 and Freq168 #= 1) or
	(Row168 #= 136 and Freq168 #= 1) or
	(Row168 #= 137 and Freq168 #= 1) or
	(Row168 #= 0 and Freq168 #= 0)), 

	((Row169 #= 82 and Freq169 #= 1) or
	(Row169 #= 84 and Freq169 #= 1) or
	(Row169 #= 87 and Freq169 #= 1) or
	(Row169 #= 92 and Freq169 #= 1) or
	(Row169 #= 93 and Freq169 #= 1) or
	(Row169 #= 94 and Freq169 #= 1) or
	(Row169 #= 95 and Freq169 #= 1) or
	(Row169 #= 96 and Freq169 #= 1) or
	(Row169 #= 97 and Freq169 #= 2) or
	(Row169 #= 98 and Freq169 #= 2) or
	(Row169 #= 99 and Freq169 #= 2) or
	(Row169 #= 100 and Freq169 #= 2) or
	(Row169 #= 101 and Freq169 #= 2) or
	(Row169 #= 102 and Freq169 #= 2) or
	(Row169 #= 103 and Freq169 #= 3) or
	(Row169 #= 104 and Freq169 #= 3) or
	(Row169 #= 105 and Freq169 #= 3) or
	(Row169 #= 106 and Freq169 #= 4) or
	(Row169 #= 115 and Freq169 #= 3) or
	(Row169 #= 116 and Freq169 #= 4) or
	(Row169 #= 117 and Freq169 #= 3) or
	(Row169 #= 118 and Freq169 #= 2) or
	(Row169 #= 119 and Freq169 #= 2) or
	(Row169 #= 120 and Freq169 #= 2) or
	(Row169 #= 121 and Freq169 #= 2) or
	(Row169 #= 122 and Freq169 #= 2) or
	(Row169 #= 123 and Freq169 #= 1) or
	(Row169 #= 124 and Freq169 #= 1) or
	(Row169 #= 125 and Freq169 #= 1) or
	(Row169 #= 126 and Freq169 #= 1) or
	(Row169 #= 127 and Freq169 #= 1) or
	(Row169 #= 128 and Freq169 #= 1) or
	(Row169 #= 129 and Freq169 #= 1) or
	(Row169 #= 130 and Freq169 #= 1) or
	(Row169 #= 131 and Freq169 #= 1) or
	(Row169 #= 132 and Freq169 #= 1) or
	(Row169 #= 135 and Freq169 #= 1) or
	(Row169 #= 136 and Freq169 #= 1) or
	(Row169 #= 137 and Freq169 #= 1) or
	(Row169 #= 140 and Freq169 #= 1) or
	(Row169 #= 0 and Freq169 #= 0)), 

	((Row170 #= 90 and Freq170 #= 1) or
	(Row170 #= 92 and Freq170 #= 1) or
	(Row170 #= 93 and Freq170 #= 1) or
	(Row170 #= 94 and Freq170 #= 1) or
	(Row170 #= 95 and Freq170 #= 1) or
	(Row170 #= 96 and Freq170 #= 1) or
	(Row170 #= 97 and Freq170 #= 1) or
	(Row170 #= 98 and Freq170 #= 1) or
	(Row170 #= 99 and Freq170 #= 1) or
	(Row170 #= 100 and Freq170 #= 2) or
	(Row170 #= 101 and Freq170 #= 1) or
	(Row170 #= 102 and Freq170 #= 2) or
	(Row170 #= 103 and Freq170 #= 2) or
	(Row170 #= 104 and Freq170 #= 2) or
	(Row170 #= 105 and Freq170 #= 3) or
	(Row170 #= 106 and Freq170 #= 3) or
	(Row170 #= 115 and Freq170 #= 3) or
	(Row170 #= 116 and Freq170 #= 3) or
	(Row170 #= 117 and Freq170 #= 2) or
	(Row170 #= 118 and Freq170 #= 2) or
	(Row170 #= 119 and Freq170 #= 1) or
	(Row170 #= 120 and Freq170 #= 2) or
	(Row170 #= 121 and Freq170 #= 2) or
	(Row170 #= 122 and Freq170 #= 1) or
	(Row170 #= 123 and Freq170 #= 1) or
	(Row170 #= 124 and Freq170 #= 1) or
	(Row170 #= 125 and Freq170 #= 1) or
	(Row170 #= 126 and Freq170 #= 1) or
	(Row170 #= 127 and Freq170 #= 1) or
	(Row170 #= 128 and Freq170 #= 1) or
	(Row170 #= 129 and Freq170 #= 1) or
	(Row170 #= 130 and Freq170 #= 1) or
	(Row170 #= 131 and Freq170 #= 1) or
	(Row170 #= 132 and Freq170 #= 1) or
	(Row170 #= 133 and Freq170 #= 1) or
	(Row170 #= 138 and Freq170 #= 1) or
	(Row170 #= 0 and Freq170 #= 0)), 

	((Row171 #= 92 and Freq171 #= 1) or
	(Row171 #= 94 and Freq171 #= 1) or
	(Row171 #= 95 and Freq171 #= 1) or
	(Row171 #= 96 and Freq171 #= 1) or
	(Row171 #= 97 and Freq171 #= 1) or
	(Row171 #= 98 and Freq171 #= 2) or
	(Row171 #= 99 and Freq171 #= 1) or
	(Row171 #= 100 and Freq171 #= 1) or
	(Row171 #= 101 and Freq171 #= 2) or
	(Row171 #= 102 and Freq171 #= 2) or
	(Row171 #= 103 and Freq171 #= 2) or
	(Row171 #= 104 and Freq171 #= 2) or
	(Row171 #= 105 and Freq171 #= 3) or
	(Row171 #= 106 and Freq171 #= 3) or
	(Row171 #= 115 and Freq171 #= 2) or
	(Row171 #= 116 and Freq171 #= 3) or
	(Row171 #= 117 and Freq171 #= 2) or
	(Row171 #= 118 and Freq171 #= 2) or
	(Row171 #= 119 and Freq171 #= 1) or
	(Row171 #= 120 and Freq171 #= 2) or
	(Row171 #= 121 and Freq171 #= 1) or
	(Row171 #= 122 and Freq171 #= 1) or
	(Row171 #= 123 and Freq171 #= 1) or
	(Row171 #= 124 and Freq171 #= 1) or
	(Row171 #= 125 and Freq171 #= 1) or
	(Row171 #= 126 and Freq171 #= 1) or
	(Row171 #= 127 and Freq171 #= 1) or
	(Row171 #= 128 and Freq171 #= 1) or
	(Row171 #= 129 and Freq171 #= 1) or
	(Row171 #= 136 and Freq171 #= 1) or
	(Row171 #= 0 and Freq171 #= 0)), 

	((Row172 #= 93 and Freq172 #= 1) or
	(Row172 #= 94 and Freq172 #= 1) or
	(Row172 #= 95 and Freq172 #= 1) or
	(Row172 #= 96 and Freq172 #= 1) or
	(Row172 #= 97 and Freq172 #= 1) or
	(Row172 #= 98 and Freq172 #= 1) or
	(Row172 #= 99 and Freq172 #= 2) or
	(Row172 #= 100 and Freq172 #= 1) or
	(Row172 #= 101 and Freq172 #= 2) or
	(Row172 #= 102 and Freq172 #= 1) or
	(Row172 #= 103 and Freq172 #= 1) or
	(Row172 #= 104 and Freq172 #= 2) or
	(Row172 #= 105 and Freq172 #= 2) or
	(Row172 #= 106 and Freq172 #= 2) or
	(Row172 #= 115 and Freq172 #= 2) or
	(Row172 #= 116 and Freq172 #= 2) or
	(Row172 #= 117 and Freq172 #= 2) or
	(Row172 #= 118 and Freq172 #= 2) or
	(Row172 #= 119 and Freq172 #= 1) or
	(Row172 #= 120 and Freq172 #= 1) or
	(Row172 #= 121 and Freq172 #= 1) or
	(Row172 #= 122 and Freq172 #= 1) or
	(Row172 #= 123 and Freq172 #= 1) or
	(Row172 #= 124 and Freq172 #= 1) or
	(Row172 #= 125 and Freq172 #= 1) or
	(Row172 #= 126 and Freq172 #= 1) or
	(Row172 #= 127 and Freq172 #= 1) or
	(Row172 #= 129 and Freq172 #= 1) or
	(Row172 #= 130 and Freq172 #= 1) or
	(Row172 #= 131 and Freq172 #= 1) or
	(Row172 #= 140 and Freq172 #= 1) or
	(Row172 #= 0 and Freq172 #= 0)), 

	((Row173 #= 92 and Freq173 #= 1) or
	(Row173 #= 93 and Freq173 #= 1) or
	(Row173 #= 94 and Freq173 #= 1) or
	(Row173 #= 95 and Freq173 #= 1) or
	(Row173 #= 96 and Freq173 #= 1) or
	(Row173 #= 97 and Freq173 #= 1) or
	(Row173 #= 98 and Freq173 #= 1) or
	(Row173 #= 99 and Freq173 #= 1) or
	(Row173 #= 100 and Freq173 #= 1) or
	(Row173 #= 101 and Freq173 #= 1) or
	(Row173 #= 102 and Freq173 #= 1) or
	(Row173 #= 103 and Freq173 #= 2) or
	(Row173 #= 104 and Freq173 #= 2) or
	(Row173 #= 105 and Freq173 #= 2) or
	(Row173 #= 106 and Freq173 #= 2) or
	(Row173 #= 115 and Freq173 #= 1) or
	(Row173 #= 116 and Freq173 #= 2) or
	(Row173 #= 117 and Freq173 #= 1) or
	(Row173 #= 118 and Freq173 #= 1) or
	(Row173 #= 119 and Freq173 #= 1) or
	(Row173 #= 120 and Freq173 #= 1) or
	(Row173 #= 121 and Freq173 #= 1) or
	(Row173 #= 122 and Freq173 #= 1) or
	(Row173 #= 123 and Freq173 #= 1) or
	(Row173 #= 124 and Freq173 #= 1) or
	(Row173 #= 125 and Freq173 #= 1) or
	(Row173 #= 126 and Freq173 #= 1) or
	(Row173 #= 127 and Freq173 #= 1) or
	(Row173 #= 129 and Freq173 #= 1) or
	(Row173 #= 130 and Freq173 #= 1) or
	(Row173 #= 131 and Freq173 #= 1) or
	(Row173 #= 0 and Freq173 #= 0)), 

	((Row174 #= 88 and Freq174 #= 1) or
	(Row174 #= 89 and Freq174 #= 1) or
	(Row174 #= 94 and Freq174 #= 1) or
	(Row174 #= 95 and Freq174 #= 1) or
	(Row174 #= 96 and Freq174 #= 1) or
	(Row174 #= 97 and Freq174 #= 1) or
	(Row174 #= 98 and Freq174 #= 1) or
	(Row174 #= 99 and Freq174 #= 1) or
	(Row174 #= 100 and Freq174 #= 1) or
	(Row174 #= 101 and Freq174 #= 1) or
	(Row174 #= 102 and Freq174 #= 1) or
	(Row174 #= 103 and Freq174 #= 1) or
	(Row174 #= 104 and Freq174 #= 1) or
	(Row174 #= 105 and Freq174 #= 1) or
	(Row174 #= 106 and Freq174 #= 1) or
	(Row174 #= 115 and Freq174 #= 1) or
	(Row174 #= 116 and Freq174 #= 1) or
	(Row174 #= 117 and Freq174 #= 1) or
	(Row174 #= 118 and Freq174 #= 1) or
	(Row174 #= 119 and Freq174 #= 1) or
	(Row174 #= 120 and Freq174 #= 1) or
	(Row174 #= 121 and Freq174 #= 1) or
	(Row174 #= 122 and Freq174 #= 1) or
	(Row174 #= 124 and Freq174 #= 1) or
	(Row174 #= 125 and Freq174 #= 1) or
	(Row174 #= 128 and Freq174 #= 1) or
	(Row174 #= 129 and Freq174 #= 1) or
	(Row174 #= 130 and Freq174 #= 1) or
	(Row174 #= 131 and Freq174 #= 1) or
	(Row174 #= 0 and Freq174 #= 0)), 

	((Row175 #= 83 and Freq175 #= 1) or
	(Row175 #= 88 and Freq175 #= 1) or
	(Row175 #= 89 and Freq175 #= 1) or
	(Row175 #= 92 and Freq175 #= 1) or
	(Row175 #= 93 and Freq175 #= 1) or
	(Row175 #= 94 and Freq175 #= 1) or
	(Row175 #= 95 and Freq175 #= 1) or
	(Row175 #= 96 and Freq175 #= 1) or
	(Row175 #= 97 and Freq175 #= 1) or
	(Row175 #= 98 and Freq175 #= 1) or
	(Row175 #= 99 and Freq175 #= 1) or
	(Row175 #= 100 and Freq175 #= 1) or
	(Row175 #= 101 and Freq175 #= 1) or
	(Row175 #= 102 and Freq175 #= 1) or
	(Row175 #= 103 and Freq175 #= 1) or
	(Row175 #= 104 and Freq175 #= 1) or
	(Row175 #= 105 and Freq175 #= 1) or
	(Row175 #= 106 and Freq175 #= 1) or
	(Row175 #= 115 and Freq175 #= 1) or
	(Row175 #= 116 and Freq175 #= 1) or
	(Row175 #= 117 and Freq175 #= 1) or
	(Row175 #= 118 and Freq175 #= 1) or
	(Row175 #= 119 and Freq175 #= 1) or
	(Row175 #= 120 and Freq175 #= 1) or
	(Row175 #= 121 and Freq175 #= 1) or
	(Row175 #= 124 and Freq175 #= 1) or
	(Row175 #= 125 and Freq175 #= 1) or
	(Row175 #= 126 and Freq175 #= 1) or
	(Row175 #= 130 and Freq175 #= 1) or
	(Row175 #= 132 and Freq175 #= 1) or
	(Row175 #= 0 and Freq175 #= 0)), 

	((Row176 #= 82 and Freq176 #= 1) or
	(Row176 #= 90 and Freq176 #= 1) or
	(Row176 #= 91 and Freq176 #= 1) or
	(Row176 #= 92 and Freq176 #= 1) or
	(Row176 #= 93 and Freq176 #= 1) or
	(Row176 #= 94 and Freq176 #= 1) or
	(Row176 #= 95 and Freq176 #= 1) or
	(Row176 #= 96 and Freq176 #= 1) or
	(Row176 #= 97 and Freq176 #= 1) or
	(Row176 #= 98 and Freq176 #= 1) or
	(Row176 #= 99 and Freq176 #= 1) or
	(Row176 #= 100 and Freq176 #= 1) or
	(Row176 #= 101 and Freq176 #= 1) or
	(Row176 #= 102 and Freq176 #= 1) or
	(Row176 #= 103 and Freq176 #= 1) or
	(Row176 #= 104 and Freq176 #= 1) or
	(Row176 #= 105 and Freq176 #= 1) or
	(Row176 #= 106 and Freq176 #= 1) or
	(Row176 #= 115 and Freq176 #= 1) or
	(Row176 #= 116 and Freq176 #= 2) or
	(Row176 #= 117 and Freq176 #= 1) or
	(Row176 #= 118 and Freq176 #= 1) or
	(Row176 #= 119 and Freq176 #= 1) or
	(Row176 #= 120 and Freq176 #= 1) or
	(Row176 #= 121 and Freq176 #= 1) or
	(Row176 #= 122 and Freq176 #= 1) or
	(Row176 #= 123 and Freq176 #= 1) or
	(Row176 #= 124 and Freq176 #= 1) or
	(Row176 #= 125 and Freq176 #= 1) or
	(Row176 #= 127 and Freq176 #= 1) or
	(Row176 #= 129 and Freq176 #= 1) or
	(Row176 #= 130 and Freq176 #= 1) or
	(Row176 #= 131 and Freq176 #= 1) or
	(Row176 #= 132 and Freq176 #= 1) or
	(Row176 #= 136 and Freq176 #= 1) or
	(Row176 #= 140 and Freq176 #= 1) or
	(Row176 #= 0 and Freq176 #= 0)), 

	((Row177 #= 90 and Freq177 #= 1) or
	(Row177 #= 92 and Freq177 #= 1) or
	(Row177 #= 93 and Freq177 #= 1) or
	(Row177 #= 94 and Freq177 #= 1) or
	(Row177 #= 95 and Freq177 #= 1) or
	(Row177 #= 96 and Freq177 #= 1) or
	(Row177 #= 97 and Freq177 #= 1) or
	(Row177 #= 98 and Freq177 #= 1) or
	(Row177 #= 99 and Freq177 #= 1) or
	(Row177 #= 100 and Freq177 #= 1) or
	(Row177 #= 101 and Freq177 #= 1) or
	(Row177 #= 102 and Freq177 #= 1) or
	(Row177 #= 103 and Freq177 #= 1) or
	(Row177 #= 104 and Freq177 #= 2) or
	(Row177 #= 105 and Freq177 #= 1) or
	(Row177 #= 106 and Freq177 #= 1) or
	(Row177 #= 115 and Freq177 #= 1) or
	(Row177 #= 116 and Freq177 #= 1) or
	(Row177 #= 117 and Freq177 #= 1) or
	(Row177 #= 118 and Freq177 #= 1) or
	(Row177 #= 119 and Freq177 #= 1) or
	(Row177 #= 120 and Freq177 #= 1) or
	(Row177 #= 121 and Freq177 #= 1) or
	(Row177 #= 122 and Freq177 #= 1) or
	(Row177 #= 123 and Freq177 #= 1) or
	(Row177 #= 125 and Freq177 #= 1) or
	(Row177 #= 126 and Freq177 #= 1) or
	(Row177 #= 128 and Freq177 #= 1) or
	(Row177 #= 129 and Freq177 #= 1) or
	(Row177 #= 131 and Freq177 #= 1) or
	(Row177 #= 0 and Freq177 #= 0)), 

	((Row178 #= 92 and Freq178 #= 1) or
	(Row178 #= 93 and Freq178 #= 1) or
	(Row178 #= 94 and Freq178 #= 1) or
	(Row178 #= 95 and Freq178 #= 1) or
	(Row178 #= 96 and Freq178 #= 1) or
	(Row178 #= 97 and Freq178 #= 1) or
	(Row178 #= 98 and Freq178 #= 1) or
	(Row178 #= 99 and Freq178 #= 1) or
	(Row178 #= 100 and Freq178 #= 1) or
	(Row178 #= 101 and Freq178 #= 1) or
	(Row178 #= 102 and Freq178 #= 1) or
	(Row178 #= 103 and Freq178 #= 1) or
	(Row178 #= 104 and Freq178 #= 2) or
	(Row178 #= 105 and Freq178 #= 1) or
	(Row178 #= 106 and Freq178 #= 1) or
	(Row178 #= 116 and Freq178 #= 1) or
	(Row178 #= 117 and Freq178 #= 1) or
	(Row178 #= 118 and Freq178 #= 1) or
	(Row178 #= 119 and Freq178 #= 1) or
	(Row178 #= 120 and Freq178 #= 1) or
	(Row178 #= 121 and Freq178 #= 1) or
	(Row178 #= 123 and Freq178 #= 1) or
	(Row178 #= 125 and Freq178 #= 1) or
	(Row178 #= 129 and Freq178 #= 1) or
	(Row178 #= 132 and Freq178 #= 1) or
	(Row178 #= 140 and Freq178 #= 1) or
	(Row178 #= 0 and Freq178 #= 0)), 

	((Row179 #= 91 and Freq179 #= 1) or
	(Row179 #= 92 and Freq179 #= 1) or
	(Row179 #= 93 and Freq179 #= 1) or
	(Row179 #= 95 and Freq179 #= 1) or
	(Row179 #= 98 and Freq179 #= 1) or
	(Row179 #= 99 and Freq179 #= 1) or
	(Row179 #= 100 and Freq179 #= 1) or
	(Row179 #= 101 and Freq179 #= 1) or
	(Row179 #= 102 and Freq179 #= 1) or
	(Row179 #= 103 and Freq179 #= 1) or
	(Row179 #= 104 and Freq179 #= 1) or
	(Row179 #= 105 and Freq179 #= 1) or
	(Row179 #= 106 and Freq179 #= 1) or
	(Row179 #= 115 and Freq179 #= 1) or
	(Row179 #= 116 and Freq179 #= 1) or
	(Row179 #= 117 and Freq179 #= 1) or
	(Row179 #= 122 and Freq179 #= 1) or
	(Row179 #= 124 and Freq179 #= 1) or
	(Row179 #= 125 and Freq179 #= 1) or
	(Row179 #= 126 and Freq179 #= 1) or
	(Row179 #= 0 and Freq179 #= 0)), 

	((Row180 #= 91 and Freq180 #= 1) or
	(Row180 #= 92 and Freq180 #= 1) or
	(Row180 #= 93 and Freq180 #= 1) or
	(Row180 #= 98 and Freq180 #= 1) or
	(Row180 #= 99 and Freq180 #= 1) or
	(Row180 #= 100 and Freq180 #= 1) or
	(Row180 #= 101 and Freq180 #= 1) or
	(Row180 #= 103 and Freq180 #= 1) or
	(Row180 #= 105 and Freq180 #= 1) or
	(Row180 #= 106 and Freq180 #= 1) or
	(Row180 #= 115 and Freq180 #= 1) or
	(Row180 #= 116 and Freq180 #= 1) or
	(Row180 #= 117 and Freq180 #= 1) or
	(Row180 #= 120 and Freq180 #= 1) or
	(Row180 #= 121 and Freq180 #= 1) or
	(Row180 #= 125 and Freq180 #= 1) or
	(Row180 #= 126 and Freq180 #= 1) or
	(Row180 #= 130 and Freq180 #= 1) or
	(Row180 #= 0 and Freq180 #= 0)), 

	((Row181 #= 87 and Freq181 #= 1) or
	(Row181 #= 91 and Freq181 #= 1) or
	(Row181 #= 92 and Freq181 #= 1) or
	(Row181 #= 93 and Freq181 #= 1) or
	(Row181 #= 98 and Freq181 #= 1) or
	(Row181 #= 99 and Freq181 #= 1) or
	(Row181 #= 100 and Freq181 #= 1) or
	(Row181 #= 101 and Freq181 #= 1) or
	(Row181 #= 105 and Freq181 #= 1) or
	(Row181 #= 117 and Freq181 #= 1) or
	(Row181 #= 120 and Freq181 #= 1) or
	(Row181 #= 121 and Freq181 #= 1) or
	(Row181 #= 125 and Freq181 #= 1) or
	(Row181 #= 130 and Freq181 #= 1) or
	(Row181 #= 0 and Freq181 #= 0)), 

	((Row182 #= 94 and Freq182 #= 1) or
	(Row182 #= 96 and Freq182 #= 1) or
	(Row182 #= 97 and Freq182 #= 1) or
	(Row182 #= 98 and Freq182 #= 1) or
	(Row182 #= 99 and Freq182 #= 1) or
	(Row182 #= 100 and Freq182 #= 1) or
	(Row182 #= 101 and Freq182 #= 1) or
	(Row182 #= 102 and Freq182 #= 1) or
	(Row182 #= 103 and Freq182 #= 1) or
	(Row182 #= 104 and Freq182 #= 1) or
	(Row182 #= 105 and Freq182 #= 1) or
	(Row182 #= 106 and Freq182 #= 1) or
	(Row182 #= 117 and Freq182 #= 1) or
	(Row182 #= 122 and Freq182 #= 1) or
	(Row182 #= 137 and Freq182 #= 1) or
	(Row182 #= 0 and Freq182 #= 0)), 

	((Row183 #= 85 and Freq183 #= 1) or
	(Row183 #= 88 and Freq183 #= 1) or
	(Row183 #= 93 and Freq183 #= 1) or
	(Row183 #= 98 and Freq183 #= 1) or
	(Row183 #= 100 and Freq183 #= 1) or
	(Row183 #= 101 and Freq183 #= 1) or
	(Row183 #= 0 and Freq183 #= 0)), 

	((Row184 #= 97 and Freq184 #= 1) or
	(Row184 #= 104 and Freq184 #= 1) or
	(Row184 #= 0 and Freq184 #= 0)), 

	((Row185 #= 93 and Freq185 #= 1) or
	(Row185 #= 0 and Freq185 #= 0)), 

	((Row186 #= 96 and Freq186 #= 1) or
	(Row186 #= 97 and Freq186 #= 1) or
	(Row186 #= 99 and Freq186 #= 1) or
	(Row186 #= 0 and Freq186 #= 0)), 

	((Row187 #= 98 and Freq187 #= 1) or
	(Row187 #= 104 and Freq187 #= 1) or
	(Row187 #= 105 and Freq187 #= 1) or
	(Row187 #= 0 and Freq187 #= 0)), 

	((Row188 #= 97 and Freq188 #= 1) or
	(Row188 #= 99 and Freq188 #= 1) or
	(Row188 #= 0 and Freq188 #= 0)), 

	((Row189 #= 0 and Freq189 #= 0)), 

	((Row190 #= 0 and Freq190 #= 0)), 

	((Row191 #= 0 and Freq191 #= 0)), 

	((Row192 #= 104 and Freq192 #= 1) or
	(Row192 #= 0 and Freq192 #= 0)), 

	((Row193 #= 0 and Freq193 #= 0)), 

	((Row194 #= 96 and Freq194 #= 1) or
	(Row194 #= 0 and Freq194 #= 0)), 

	((Row195 #= 77 and Freq195 #= 1) or
	(Row195 #= 0 and Freq195 #= 0)), 

	((Row196 #= 0 and Freq196 #= 0)), 

	((Row197 #= 152 and Freq197 #= 1) or
	(Row197 #= 0 and Freq197 #= 0)), 

	((Row198 #= 0 and Freq198 #= 0)), 

	((Row199 #= 0 and Freq199 #= 0)), 

	((Row200 #= 0 and Freq200 #= 0)), 

	((Row201 #= 0 and Freq201 #= 0)), 

	((Row202 #= 0 and Freq202 #= 0)), 

	((Row203 #= 0 and Freq203 #= 0)), 

	((Row204 #= 0 and Freq204 #= 0)), 

	((Row205 #= 0 and Freq205 #= 0)), 

	((Row206 #= 0 and Freq206 #= 0)), 

	((Row207 #= 0 and Freq207 #= 0)), 

	((Row208 #= 0 and Freq208 #= 0)), 

	((Row209 #= 0 and Freq209 #= 0)), 

	((Row210 #= 0 and Freq210 #= 0)), 

	((Row211 #= 0 and Freq211 #= 0)), 

	((Row212 #= 0 and Freq212 #= 0)), 

	((Row213 #= 0 and Freq213 #= 0)), 

	((Row214 #= 0 and Freq214 #= 0)), 

	((Row215 #= 0 and Freq215 #= 0)), 

	((Row216 #= 0 and Freq216 #= 0)), 

	((Row217 #= 0 and Freq217 #= 0)), 

	((Row218 #= 0 and Freq218 #= 0)), 

	((Row219 #= 0 and Freq219 #= 0)), 

	((Row220 #= 0 and Freq220 #= 0)), 

	((Row221 #= 0 and Freq221 #= 0)), 

	((Row222 #= 0 and Freq222 #= 0)), 

	((Row223 #= 0 and Freq223 #= 0)), 

	((Row224 #= 0 and Freq224 #= 0)), 

	((Row225 #= 0 and Freq225 #= 0)), 

	((Row226 #= 0 and Freq226 #= 0)), 

	((Row227 #= 0 and Freq227 #= 0)), 

	((Row228 #= 0 and Freq228 #= 0)), 

	((Row229 #= 0 and Freq229 #= 0)), 

	((Row230 #= 0 and Freq230 #= 0)), 

	((Row231 #= 0 and Freq231 #= 0)), 

	((Row232 #= 0 and Freq232 #= 0)), 

	((Row233 #= 0 and Freq233 #= 0)), 

	((Row234 #= 0 and Freq234 #= 0)), 

	((Row235 #= 0 and Freq235 #= 0)), 

	((Row236 #= 0 and Freq236 #= 0)), 

	((Row237 #= 0 and Freq237 #= 0)), 

	((Row238 #= 0 and Freq238 #= 0)), 

	((Row239 #= 0 and Freq239 #= 0)), 

	((Row240 #= 0 and Freq240 #= 0)), 

	((Row241 #= 0 and Freq241 #= 0)), 

	((Row242 #= 0 and Freq242 #= 0)), 

	((Row243 #= 0 and Freq243 #= 0)), 

	((Row244 #= 0 and Freq244 #= 0)), 

	((Row245 #= 0 and Freq245 #= 0)), 

	((Row246 #= 0 and Freq246 #= 0)), 

	((Row247 #= 0 and Freq247 #= 0)), 

	((Row248 #= 0 and Freq248 #= 0)), 

	((Row249 #= 0 and Freq249 #= 0)), 

	((Row250 #= 0 and Freq250 #= 0)), 

	((Row251 #= 0 and Freq251 #= 0)), 

	((Row252 #= 0 and Freq252 #= 0)), 

	((Row253 #= 0 and Freq253 #= 0)), 

	((Row254 #= 0 and Freq254 #= 0)), 

	((Row255 #= 0 and Freq255 #= 0)), 

	((Row256 #= 0 and Freq256 #= 0)), 

	((Row257 #= 0 and Freq257 #= 0)), 

	((Row258 #= 0 and Freq258 #= 0)), 

	((Row259 #= 0 and Freq259 #= 0)), 

	((Row260 #= 0 and Freq260 #= 0)), 

	((Row261 #= 0 and Freq261 #= 0)), 

	((Row262 #= 0 and Freq262 #= 0)), 

	((Row263 #= 0 and Freq263 #= 0)), 

	((Row264 #= 0 and Freq264 #= 0)), 

	((Row265 #= 245 and Freq265 #= 1) or
	(Row265 #= 0 and Freq265 #= 0)), 

	((Row266 #= 0 and Freq266 #= 0)), 

	((Row267 #= 0 and Freq267 #= 0)), 

	((Row268 #= 0 and Freq268 #= 0)), 

	((Row269 #= 0 and Freq269 #= 0)), 

	((Row270 #= 0 and Freq270 #= 0)), 

	((Row271 #= 0 and Freq271 #= 0)), 

	((Row272 #= 0 and Freq272 #= 0)), 

	((Row273 #= 0 and Freq273 #= 0)), 

	((Row274 #= 0 and Freq274 #= 0)), 

	((Row275 #= 0 and Freq275 #= 0)), 

	((Row276 #= 0 and Freq276 #= 0)), 

	((Row277 #= 0 and Freq277 #= 0)), 

	((Row278 #= 0 and Freq278 #= 0)), 

	((Row279 #= 0 and Freq279 #= 0)), 

	((Row280 #= 0 and Freq280 #= 0)), 

	((Row281 #= 0 and Freq281 #= 0)), 

	((Row282 #= 0 and Freq282 #= 0)), 

	((Row283 #= 0 and Freq283 #= 0)), 

	((Row284 #= 0 and Freq284 #= 0)), 

	((Row285 #= 0 and Freq285 #= 0)), 

	((Row286 #= 0 and Freq286 #= 0)), 

	((Row287 #= 0 and Freq287 #= 0)), 

	((Row288 #= 0 and Freq288 #= 0)), 

	((Row289 #= 0 and Freq289 #= 0)), 

	((Row290 #= 0 and Freq290 #= 0)), 

	((Row291 #= 0 and Freq291 #= 0)), 

	((Row292 #= 0 and Freq292 #= 0)), 

	((Row293 #= 0 and Freq293 #= 0)), 

	((Row294 #= 0 and Freq294 #= 0)), 

	((Row295 #= 0 and Freq295 #= 0)), 

	((Row296 #= 0 and Freq296 #= 0)), 

	((Row297 #= 0 and Freq297 #= 0)), 

	((Row298 #= 0 and Freq298 #= 0)), 

	((Row299 #= 0 and Freq299 #= 0)), 

	((Row300 #= 0 and Freq300 #= 0)), 

	((Row301 #= 0 and Freq301 #= 0)), 

	((Row302 #= 0 and Freq302 #= 0)), 

	((Row303 #= 0 and Freq303 #= 0)), 

	((Row304 #= 0 and Freq304 #= 0)), 

	((Row305 #= 0 and Freq305 #= 0)), 

	((Row306 #= 0 and Freq306 #= 0)), 

	((Row307 #= 0 and Freq307 #= 0)), 

	((Row308 #= 0 and Freq308 #= 0)), 

	((Row309 #= 0 and Freq309 #= 0)), 

	((Row310 #= 0 and Freq310 #= 0)), 

	((Row311 #= 0 and Freq311 #= 0)), 

	((Row312 #= 0 and Freq312 #= 0)), 

	((Row313 #= 0 and Freq313 #= 0)), 

	((Row314 #= 0 and Freq314 #= 0)), 

	((Row315 #= 0 and Freq315 #= 0)), 

	((Row316 #= 0 and Freq316 #= 0)), 

	((Row317 #= 0 and Freq317 #= 0)), 

	((Row318 #= 0 and Freq318 #= 0)), 

	((Row319 #= 245 and Freq319 #= 1) or
	(Row319 #= 0 and Freq319 #= 0)), 

	((Row320 #= 0 and Freq320 #= 0)), 

	((Row321 #= 0 and Freq321 #= 0)), 

	((Row322 #= 0 and Freq322 #= 0)), 

	((Row323 #= 0 and Freq323 #= 0)), 

	((Row324 #= 0 and Freq324 #= 0)), 

	((Row325 #= 0 and Freq325 #= 0)), 

	((Row326 #= 0 and Freq326 #= 0)), 

	((Row327 #= 0 and Freq327 #= 0)), 

	((Row328 #= 0 and Freq328 #= 0)), 

	((Row329 #= 0 and Freq329 #= 0)), 

	((Row330 #= 0 and Freq330 #= 0)), 

	((Row331 #= 0 and Freq331 #= 0)), 

	((Row332 #= 0 and Freq332 #= 0)), 

	((Row333 #= 0 and Freq333 #= 0)), 

	((Row334 #= 0 and Freq334 #= 0)), 

	((Row335 #= 0 and Freq335 #= 0)), 

	((Row336 #= 0 and Freq336 #= 0)), 

	((Row337 #= 0 and Freq337 #= 0)), 

	((Row338 #= 0 and Freq338 #= 0)), 

	((Row339 #= 0 and Freq339 #= 0)), 

	((Row340 #= 245 and Freq340 #= 1) or
	(Row340 #= 0 and Freq340 #= 0)), 

	((Row341 #= 0 and Freq341 #= 0)), 

	((Row342 #= 0 and Freq342 #= 0)), 

	((Row343 #= 0 and Freq343 #= 0)), 

	((Row344 #= 0 and Freq344 #= 0)), 

	((Row345 #= 0 and Freq345 #= 0)), 

	((Row346 #= 0 and Freq346 #= 0)), 

	((Row347 #= 0 and Freq347 #= 0)), 

	((Row348 #= 0 and Freq348 #= 0)), 

	((Row349 #= 0 and Freq349 #= 0)), 

	((Row350 #= 245 and Freq350 #= 1) or
	(Row350 #= 0 and Freq350 #= 0)), 

	((Row351 #= 0 and Freq351 #= 0)), 

	((Row352 #= 0 and Freq352 #= 0)), 

	((Row353 #= 0 and Freq353 #= 0)), 

	((Row354 #= 0 and Freq354 #= 0)), 

	((Row355 #= 0 and Freq355 #= 0)), 

	((Row356 #= 0 and Freq356 #= 0)), 

	((Row357 #= 0 and Freq357 #= 0)), 

	((Row358 #= 0 and Freq358 #= 0)), 

	((Row359 #= 245 and Freq359 #= 1) or
	(Row359 #= 0 and Freq359 #= 0)), 

	((Row360 #= 0 and Freq360 #= 0)), 

	((Row361 #= 0 and Freq361 #= 0)), 

	((Row362 #= 0 and Freq362 #= 0)), 

	((Row363 #= 0 and Freq363 #= 0)), 

	((Row364 #= 245 and Freq364 #= 1) or
	(Row364 #= 0 and Freq364 #= 0)), 

	((Row365 #= 245 and Freq365 #= 1) or
	(Row365 #= 0 and Freq365 #= 0)), 

	((Row366 #= 245 and Freq366 #= 1) or
	(Row366 #= 246 and Freq366 #= 1) or
	(Row366 #= 0 and Freq366 #= 0)), 

	((Row367 #= 245 and Freq367 #= 1) or
	(Row367 #= 0 and Freq367 #= 0)), 

	((Row368 #= 245 and Freq368 #= 1) or
	(Row368 #= 0 and Freq368 #= 0)), 

	((Row369 #= 245 and Freq369 #= 1) or
	(Row369 #= 0 and Freq369 #= 0)), 

	((Row370 #= 0 and Freq370 #= 0)), 

	((Row371 #= 0 and Freq371 #= 0)), 

	((Row372 #= 0 and Freq372 #= 0)), 

	((Row373 #= 0 and Freq373 #= 0)), 

	((Row374 #= 0 and Freq374 #= 0)), 

	((Row375 #= 0 and Freq375 #= 0)), 

	((Row376 #= 0 and Freq376 #= 0)), 

	((Row377 #= 0 and Freq377 #= 0)), 

	((Row378 #= 0 and Freq378 #= 0)), 

	((Row379 #= 0 and Freq379 #= 0)), 

	((Row380 #= 0 and Freq380 #= 0)), 

	((Row381 #= 0 and Freq381 #= 0)), 

	((Row382 #= 0 and Freq382 #= 0)), 

	((Row383 #= 0 and Freq383 #= 0)), 

	((Row384 #= 0 and Freq384 #= 0)), 

	((Row385 #= 0 and Freq385 #= 0)), 

	((Row386 #= 0 and Freq386 #= 0)), 

	((Row387 #= 0 and Freq387 #= 0)), 

	((Row388 #= 0 and Freq388 #= 0)), 

	((Row389 #= 0 and Freq389 #= 0)), 

	((Row390 #= 0 and Freq390 #= 0)), 

	((Row391 #= 0 and Freq391 #= 0)), 

	((Row392 #= 0 and Freq392 #= 0)), 

	((Row393 #= 245 and Freq393 #= 1) or
	(Row393 #= 0 and Freq393 #= 0)), 

	((Row394 #= 0 and Freq394 #= 0)), 

	((Row395 #= 0 and Freq395 #= 0)), 

	((Row396 #= 0 and Freq396 #= 0)), 

	((Row397 #= 0 and Freq397 #= 0)), 

	((Row398 #= 0 and Freq398 #= 0)), 

	((Row399 #= 0 and Freq399 #= 0)), 

	((Row400 #= 0 and Freq400 #= 0)), 

	((Row401 #= 0 and Freq401 #= 0)), 

	((Row402 #= 0 and Freq402 #= 0)), 

	((Row403 #= 0 and Freq403 #= 0)), 

	((Row404 #= 0 and Freq404 #= 0)), 

	((Row405 #= 0 and Freq405 #= 0)), 

	((Row406 #= 245 and Freq406 #= 1) or
	(Row406 #= 0 and Freq406 #= 0)), 

	((Row407 #= 0 and Freq407 #= 0)), 

	((Row408 #= 0 and Freq408 #= 0)), 

	((Row409 #= 245 and Freq409 #= 1) or
	(Row409 #= 0 and Freq409 #= 0)), 

	((Row410 #= 245 and Freq410 #= 1) or
	(Row410 #= 246 and Freq410 #= 1) or
	(Row410 #= 0 and Freq410 #= 0)), 

	((Row411 #= 245 and Freq411 #= 1) or
	(Row411 #= 0 and Freq411 #= 0)), 

	((Row412 #= 246 and Freq412 #= 1) or
	(Row412 #= 0 and Freq412 #= 0)), 

	((Row413 #= 245 and Freq413 #= 1) or
	(Row413 #= 0 and Freq413 #= 0)), 

	((Row414 #= 246 and Freq414 #= 1) or
	(Row414 #= 0 and Freq414 #= 0)), 

	((Row415 #= 0 and Freq415 #= 0)), 

	((Row416 #= 245 and Freq416 #= 1) or
	(Row416 #= 0 and Freq416 #= 0)), 

	((Row417 #= 245 and Freq417 #= 1) or
	(Row417 #= 0 and Freq417 #= 0)), 

	((Row418 #= 245 and Freq418 #= 1) or
	(Row418 #= 0 and Freq418 #= 0)), 

	((Row419 #= 0 and Freq419 #= 0)), 

	((Row420 #= 245 and Freq420 #= 1) or
	(Row420 #= 0 and Freq420 #= 0)), 

	((Row421 #= 245 and Freq421 #= 1) or
	(Row421 #= 246 and Freq421 #= 1) or
	(Row421 #= 0 and Freq421 #= 0)), 

	((Row422 #= 245 and Freq422 #= 1) or
	(Row422 #= 246 and Freq422 #= 1) or
	(Row422 #= 0 and Freq422 #= 0)), 

	((Row423 #= 245 and Freq423 #= 1) or
	(Row423 #= 246 and Freq423 #= 1) or
	(Row423 #= 0 and Freq423 #= 0)), 

	((Row424 #= 245 and Freq424 #= 1) or
	(Row424 #= 246 and Freq424 #= 1) or
	(Row424 #= 0 and Freq424 #= 0)), 

	((Row425 #= 245 and Freq425 #= 1) or
	(Row425 #= 246 and Freq425 #= 1) or
	(Row425 #= 0 and Freq425 #= 0)), 

	((Row426 #= 245 and Freq426 #= 2) or
	(Row426 #= 246 and Freq426 #= 2) or
	(Row426 #= 0 and Freq426 #= 0)), 

	((Row427 #= 245 and Freq427 #= 1) or
	(Row427 #= 246 and Freq427 #= 1) or
	(Row427 #= 0 and Freq427 #= 0)), 

	((Row428 #= 246 and Freq428 #= 1) or
	(Row428 #= 0 and Freq428 #= 0)), 

	((Row429 #= 245 and Freq429 #= 1) or
	(Row429 #= 246 and Freq429 #= 1) or
	(Row429 #= 0 and Freq429 #= 0)), 

	((Row430 #= 245 and Freq430 #= 1) or
	(Row430 #= 0 and Freq430 #= 0)), 

	((Row431 #= 0 and Freq431 #= 0)), 

	((Row432 #= 245 and Freq432 #= 1) or
	(Row432 #= 0 and Freq432 #= 0)), 

	((Row433 #= 0 and Freq433 #= 0)), 

	((Row434 #= 245 and Freq434 #= 1) or
	(Row434 #= 246 and Freq434 #= 1) or
	(Row434 #= 0 and Freq434 #= 0)), 

	((Row435 #= 0 and Freq435 #= 0)), 

	((Row436 #= 245 and Freq436 #= 1) or
	(Row436 #= 246 and Freq436 #= 1) or
	(Row436 #= 0 and Freq436 #= 0)), 

	((Row437 #= 245 and Freq437 #= 1) or
	(Row437 #= 0 and Freq437 #= 0)), 

	((Row438 #= 245 and Freq438 #= 1) or
	(Row438 #= 0 and Freq438 #= 0)), 

	((Row439 #= 245 and Freq439 #= 1) or
	(Row439 #= 0 and Freq439 #= 0)), 

	((Row440 #= 245 and Freq440 #= 1) or
	(Row440 #= 246 and Freq440 #= 1) or
	(Row440 #= 0 and Freq440 #= 0)), 

	((Row441 #= 245 and Freq441 #= 1) or
	(Row441 #= 246 and Freq441 #= 1) or
	(Row441 #= 0 and Freq441 #= 0)), 

	((Row442 #= 245 and Freq442 #= 1) or
	(Row442 #= 246 and Freq442 #= 1) or
	(Row442 #= 0 and Freq442 #= 0)), 

	((Row443 #= 0 and Freq443 #= 0)), 

	((Row444 #= 245 and Freq444 #= 1) or
	(Row444 #= 0 and Freq444 #= 0)), 

	((Row445 #= 245 and Freq445 #= 3) or
	(Row445 #= 246 and Freq445 #= 2) or
	(Row445 #= 0 and Freq445 #= 0)), 

	((Row446 #= 245 and Freq446 #= 1) or
	(Row446 #= 246 and Freq446 #= 1) or
	(Row446 #= 0 and Freq446 #= 0)), 

	((Row447 #= 245 and Freq447 #= 1) or
	(Row447 #= 246 and Freq447 #= 1) or
	(Row447 #= 0 and Freq447 #= 0)), 

	((Row448 #= 245 and Freq448 #= 1) or
	(Row448 #= 246 and Freq448 #= 1) or
	(Row448 #= 0 and Freq448 #= 0)), 

	((Row449 #= 245 and Freq449 #= 1) or
	(Row449 #= 246 and Freq449 #= 1) or
	(Row449 #= 0 and Freq449 #= 0)), 

	((Row450 #= 245 and Freq450 #= 1) or
	(Row450 #= 0 and Freq450 #= 0)), 

	((Row451 #= 0 and Freq451 #= 0)), 

	((Row452 #= 0 and Freq452 #= 0)), 

	((Row453 #= 0 and Freq453 #= 0)), 

	((Row454 #= 0 and Freq454 #= 0)), 

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
