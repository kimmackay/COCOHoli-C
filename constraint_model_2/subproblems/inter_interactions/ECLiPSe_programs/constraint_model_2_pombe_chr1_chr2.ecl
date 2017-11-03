% license: This work is licensed under the Creative Commons Attribution-NonCommercial-
% ShareAlike 3.0 Unported License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
% PO Box 1866, Mountain View, CA 94042, USA.

% define a structure for encoding an interaction
:- local struct( interaction(row, freq) ).

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

% % enforce_symmetry(AllVars) ensures that for each  
% for each term in AllVars bound to a non-zero   
% value the corresponding element at index 'term' in  
% AllVars (ie. AllVars[term]) is bound to zero. This 
% ensures that 	each genomic bin can only be involved 
% in one selected interaction in the solution set.
enforce_symmetry(RowFile, FreqFile, AllVars, Rows, Freqs) :- 
	% for each terms in AllVars
	( fromto(AllVars, Vars, VarsRem, []), param(Rows) do 

	% chose a value for X from it's domain
	gfd_update,
	delete(Var, Vars, VarsRem, 0, first_fail), % dynamic var-select
	arg(freq of interaction, Var, F), % get the frequency of the interaction pair
	indomain(F, max), % select a value

		% If the value bound to X is non-zero 
		( (F \= 0) -> 

			% The variable Var is the list position of the element K 
			% in AllVars. Note in ECLiPSE, indexing starts at 1. 
			arg(row of interaction, Var, R), % get the row of the interaction pair
			element(R, Rows, K), 

			% K must be bound to zero 
			K #= 0 
		; 
			true 
		) 
	),
	
	% output the results
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
	close(ROW_OUT),
	
	%% print the statistics
	statistics.

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
	Rows = [Row1, Row2, Row3, Row4, Row5, Row6, Row7, Row8, Row9, Row10, Row11, Row12, Row13, Row14, Row15, Row16, Row17, Row18, Row19, Row20, Row21, Row22, Row23, Row24, Row25, Row26, Row27, Row28, Row29, Row30, Row31, Row32, Row33, Row34, Row35, Row36, Row37, Row38, Row39, Row40, Row41, Row42, Row43, Row44, Row45, Row46, Row47, Row48, Row49, Row50, Row51, Row52, Row53, Row54, Row55, Row56, Row57, Row58, Row59, Row60, Row61, Row62, Row63, Row64, Row65, Row66, Row67, Row68, Row69, Row70, Row71, Row72, Row73, Row74, Row75, Row76, Row77, Row78, Row79, Row80, Row81, Row82, Row83, Row84, Row85, Row86, Row87, Row88, Row89, Row90, Row91, Row92, Row93, Row94, Row95, Row96, Row97, Row98, Row99, Row100, Row101, Row102, Row103, Row104, Row105, Row106, Row107, Row108, Row109, Row110, Row111, Row112, Row113, Row114, Row115, Row116, Row117, Row118, Row119, Row120, Row121, Row122, Row123, Row124, Row125, Row126, Row127, Row128, Row129, Row130, Row131, Row132, Row133, Row134, Row135, Row136, Row137, Row138, Row139, Row140, Row141, Row142, Row143, Row144, Row145, Row146, Row147, Row148, Row149, Row150, Row151, Row152, Row153, Row154, Row155, Row156, Row157, Row158, Row159, Row160, Row161, Row162, Row163, Row164, Row165, Row166, Row167, Row168, Row169, Row170, Row171, Row172, Row173, Row174, Row175, Row176, Row177, Row178, Row179, Row180, Row181, Row182, Row183, Row184, Row185, Row186, Row187, Row188, Row189, Row190, Row191, Row192, Row193, Row194, Row195, Row196, Row197, Row198, Row199, Row200, Row201, Row202, Row203, Row204, Row205, Row206, Row207, Row208, Row209, Row210, Row211, Row212, Row213, Row214, Row215, Row216, Row217, Row218, Row219, Row220, Row221, Row222, Row223, Row224, Row225, Row226, Row227, Row228, Row229, Row230, Row231, Row232, Row233, Row234, Row235, Row236, Row237, Row238, Row239, Row240, Row241, Row242, Row243, Row244, Row245, Row246, Row247, Row248, Row249, Row250, Row251, Row252, Row253, Row254, Row255, Row256, Row257, Row258, Row259, Row260, Row261, Row262, Row263, Row264, Row265, Row266, Row267, Row268, Row269, Row270, Row271, Row272, Row273, Row274, Row275, Row276, Row277, Row278, Row279, Row280, Row281, Row282, Row283, Row284, Row285, Row286, Row287, Row288, Row289, Row290, Row291, Row292, Row293, Row294, Row295, Row296, Row297, Row298, Row299, Row300, Row301, Row302, Row303, Row304, Row305, Row306, Row307, Row308, Row309, Row310, Row311, Row312, Row313, Row314, Row315, Row316, Row317, Row318, Row319, Row320, Row321, Row322, Row323, Row324, Row325, Row326, Row327, Row328, Row329, Row330, Row331, Row332, Row333, Row334, Row335, Row336, Row337, Row338, Row339, Row340, Row341, Row342, Row343, Row344, Row345, Row346, Row347, Row348, Row349, Row350, Row351, Row352, Row353, Row354, Row355, Row356, Row357, Row358, Row359, Row360, Row361, Row362, Row363, Row364, Row365, Row366, Row367, Row368, Row369, Row370, Row371, Row372, Row373, Row374, Row375, Row376, Row377, Row378, Row379, Row380, Row381, Row382, Row383, Row384, Row385, Row386, Row387, Row388, Row389, Row390, Row391, Row392, Row393, Row394, Row395, Row396, Row397, Row398, Row399, Row400, Row401, Row402, Row403, Row404, Row405, Row406, Row407, Row408, Row409, Row410, Row411, Row412, Row413, Row414, Row415, Row416, Row417, Row418, Row419, Row420, Row421, Row422, Row423, Row424, Row425, Row426, Row427, Row428, Row429, Row430, Row431, Row432, Row433, Row434, Row435, Row436, Row437, Row438, Row439, Row440, Row441, Row442, Row443, Row444, Row445, Row446, Row447, Row448, Row449, Row450, Row451, Row452, Row453, Row454, Row455, Row456, Row457, Row458, Row459, Row460, Row461, Row462, Row463, Row464, Row465, Row466, Row467, Row468, Row469, Row470, Row471, Row472, Row473, Row474, Row475, Row476, Row477, Row478, Row479, Row480, Row481, Row482, Row483, Row484, Row485, Row486, Row487, Row488, Row489, Row490, Row491, Row492, Row493, Row494, Row495, Row496, Row497, Row498, Row499, Row500, Row501, Row502, Row503, Row504, Row505, Row506, Row507, Row508, Row509, Row510, Row511, Row512, Row513, Row514, Row515, Row516, Row517, Row518, Row519, Row520, Row521, Row522, Row523, Row524, Row525, Row526, Row527, Row528, Row529, Row530, Row531, Row532, Row533, Row534, Row535, Row536, Row537, Row538, Row539, Row540, Row541, Row542, Row543, Row544, Row545, Row546, Row547, Row548, Row549, Row550, Row551, Row552, Row553, Row554, Row555, Row556, Row557, Row558],

	% The list Freqs has one variable for each row of the 
	% whole-genome contact map	
	Freqs = [Freq1, Freq2, Freq3, Freq4, Freq5, Freq6, Freq7, Freq8, Freq9, Freq10, Freq11, Freq12, Freq13, Freq14, Freq15, Freq16, Freq17, Freq18, Freq19, Freq20, Freq21, Freq22, Freq23, Freq24, Freq25, Freq26, Freq27, Freq28, Freq29, Freq30, Freq31, Freq32, Freq33, Freq34, Freq35, Freq36, Freq37, Freq38, Freq39, Freq40, Freq41, Freq42, Freq43, Freq44, Freq45, Freq46, Freq47, Freq48, Freq49, Freq50, Freq51, Freq52, Freq53, Freq54, Freq55, Freq56, Freq57, Freq58, Freq59, Freq60, Freq61, Freq62, Freq63, Freq64, Freq65, Freq66, Freq67, Freq68, Freq69, Freq70, Freq71, Freq72, Freq73, Freq74, Freq75, Freq76, Freq77, Freq78, Freq79, Freq80, Freq81, Freq82, Freq83, Freq84, Freq85, Freq86, Freq87, Freq88, Freq89, Freq90, Freq91, Freq92, Freq93, Freq94, Freq95, Freq96, Freq97, Freq98, Freq99, Freq100, Freq101, Freq102, Freq103, Freq104, Freq105, Freq106, Freq107, Freq108, Freq109, Freq110, Freq111, Freq112, Freq113, Freq114, Freq115, Freq116, Freq117, Freq118, Freq119, Freq120, Freq121, Freq122, Freq123, Freq124, Freq125, Freq126, Freq127, Freq128, Freq129, Freq130, Freq131, Freq132, Freq133, Freq134, Freq135, Freq136, Freq137, Freq138, Freq139, Freq140, Freq141, Freq142, Freq143, Freq144, Freq145, Freq146, Freq147, Freq148, Freq149, Freq150, Freq151, Freq152, Freq153, Freq154, Freq155, Freq156, Freq157, Freq158, Freq159, Freq160, Freq161, Freq162, Freq163, Freq164, Freq165, Freq166, Freq167, Freq168, Freq169, Freq170, Freq171, Freq172, Freq173, Freq174, Freq175, Freq176, Freq177, Freq178, Freq179, Freq180, Freq181, Freq182, Freq183, Freq184, Freq185, Freq186, Freq187, Freq188, Freq189, Freq190, Freq191, Freq192, Freq193, Freq194, Freq195, Freq196, Freq197, Freq198, Freq199, Freq200, Freq201, Freq202, Freq203, Freq204, Freq205, Freq206, Freq207, Freq208, Freq209, Freq210, Freq211, Freq212, Freq213, Freq214, Freq215, Freq216, Freq217, Freq218, Freq219, Freq220, Freq221, Freq222, Freq223, Freq224, Freq225, Freq226, Freq227, Freq228, Freq229, Freq230, Freq231, Freq232, Freq233, Freq234, Freq235, Freq236, Freq237, Freq238, Freq239, Freq240, Freq241, Freq242, Freq243, Freq244, Freq245, Freq246, Freq247, Freq248, Freq249, Freq250, Freq251, Freq252, Freq253, Freq254, Freq255, Freq256, Freq257, Freq258, Freq259, Freq260, Freq261, Freq262, Freq263, Freq264, Freq265, Freq266, Freq267, Freq268, Freq269, Freq270, Freq271, Freq272, Freq273, Freq274, Freq275, Freq276, Freq277, Freq278, Freq279, Freq280, Freq281, Freq282, Freq283, Freq284, Freq285, Freq286, Freq287, Freq288, Freq289, Freq290, Freq291, Freq292, Freq293, Freq294, Freq295, Freq296, Freq297, Freq298, Freq299, Freq300, Freq301, Freq302, Freq303, Freq304, Freq305, Freq306, Freq307, Freq308, Freq309, Freq310, Freq311, Freq312, Freq313, Freq314, Freq315, Freq316, Freq317, Freq318, Freq319, Freq320, Freq321, Freq322, Freq323, Freq324, Freq325, Freq326, Freq327, Freq328, Freq329, Freq330, Freq331, Freq332, Freq333, Freq334, Freq335, Freq336, Freq337, Freq338, Freq339, Freq340, Freq341, Freq342, Freq343, Freq344, Freq345, Freq346, Freq347, Freq348, Freq349, Freq350, Freq351, Freq352, Freq353, Freq354, Freq355, Freq356, Freq357, Freq358, Freq359, Freq360, Freq361, Freq362, Freq363, Freq364, Freq365, Freq366, Freq367, Freq368, Freq369, Freq370, Freq371, Freq372, Freq373, Freq374, Freq375, Freq376, Freq377, Freq378, Freq379, Freq380, Freq381, Freq382, Freq383, Freq384, Freq385, Freq386, Freq387, Freq388, Freq389, Freq390, Freq391, Freq392, Freq393, Freq394, Freq395, Freq396, Freq397, Freq398, Freq399, Freq400, Freq401, Freq402, Freq403, Freq404, Freq405, Freq406, Freq407, Freq408, Freq409, Freq410, Freq411, Freq412, Freq413, Freq414, Freq415, Freq416, Freq417, Freq418, Freq419, Freq420, Freq421, Freq422, Freq423, Freq424, Freq425, Freq426, Freq427, Freq428, Freq429, Freq430, Freq431, Freq432, Freq433, Freq434, Freq435, Freq436, Freq437, Freq438, Freq439, Freq440, Freq441, Freq442, Freq443, Freq444, Freq445, Freq446, Freq447, Freq448, Freq449, Freq450, Freq451, Freq452, Freq453, Freq454, Freq455, Freq456, Freq457, Freq458, Freq459, Freq460, Freq461, Freq462, Freq463, Freq464, Freq465, Freq466, Freq467, Freq468, Freq469, Freq470, Freq471, Freq472, Freq473, Freq474, Freq475, Freq476, Freq477, Freq478, Freq479, Freq480, Freq481, Freq482, Freq483, Freq484, Freq485, Freq486, Freq487, Freq488, Freq489, Freq490, Freq491, Freq492, Freq493, Freq494, Freq495, Freq496, Freq497, Freq498, Freq499, Freq500, Freq501, Freq502, Freq503, Freq504, Freq505, Freq506, Freq507, Freq508, Freq509, Freq510, Freq511, Freq512, Freq513, Freq514, Freq515, Freq516, Freq517, Freq518, Freq519, Freq520, Freq521, Freq522, Freq523, Freq524, Freq525, Freq526, Freq527, Freq528, Freq529, Freq530, Freq531, Freq532, Freq533, Freq534, Freq535, Freq536, Freq537, Freq538, Freq539, Freq540, Freq541, Freq542, Freq543, Freq544, Freq545, Freq546, Freq547, Freq548, Freq549, Freq550, Freq551, Freq552, Freq553, Freq554, Freq555, Freq556, Freq557, Freq558],

	% The list Interactions has one variable for each interaction 
	% in the solution set	
	Interactions = [interaction(Row1, Freq1), interaction(Row2, Freq2), interaction(Row3, Freq3), interaction(Row4, Freq4), interaction(Row5, Freq5), interaction(Row6, Freq6), interaction(Row7, Freq7), interaction(Row8, Freq8), interaction(Row9, Freq9), interaction(Row10, Freq10), interaction(Row11, Freq11), interaction(Row12, Freq12), interaction(Row13, Freq13), interaction(Row14, Freq14), interaction(Row15, Freq15), interaction(Row16, Freq16), interaction(Row17, Freq17), interaction(Row18, Freq18), interaction(Row19, Freq19), interaction(Row20, Freq20), interaction(Row21, Freq21), interaction(Row22, Freq22), interaction(Row23, Freq23), interaction(Row24, Freq24), interaction(Row25, Freq25), interaction(Row26, Freq26), interaction(Row27, Freq27), interaction(Row28, Freq28), interaction(Row29, Freq29), interaction(Row30, Freq30), interaction(Row31, Freq31), interaction(Row32, Freq32), interaction(Row33, Freq33), interaction(Row34, Freq34), interaction(Row35, Freq35), interaction(Row36, Freq36), interaction(Row37, Freq37), interaction(Row38, Freq38), interaction(Row39, Freq39), interaction(Row40, Freq40), interaction(Row41, Freq41), interaction(Row42, Freq42), interaction(Row43, Freq43), interaction(Row44, Freq44), interaction(Row45, Freq45), interaction(Row46, Freq46), interaction(Row47, Freq47), interaction(Row48, Freq48), interaction(Row49, Freq49), interaction(Row50, Freq50), interaction(Row51, Freq51), interaction(Row52, Freq52), interaction(Row53, Freq53), interaction(Row54, Freq54), interaction(Row55, Freq55), interaction(Row56, Freq56), interaction(Row57, Freq57), interaction(Row58, Freq58), interaction(Row59, Freq59), interaction(Row60, Freq60), interaction(Row61, Freq61), interaction(Row62, Freq62), interaction(Row63, Freq63), interaction(Row64, Freq64), interaction(Row65, Freq65), interaction(Row66, Freq66), interaction(Row67, Freq67), interaction(Row68, Freq68), interaction(Row69, Freq69), interaction(Row70, Freq70), interaction(Row71, Freq71), interaction(Row72, Freq72), interaction(Row73, Freq73), interaction(Row74, Freq74), interaction(Row75, Freq75), interaction(Row76, Freq76), interaction(Row77, Freq77), interaction(Row78, Freq78), interaction(Row79, Freq79), interaction(Row80, Freq80), interaction(Row81, Freq81), interaction(Row82, Freq82), interaction(Row83, Freq83), interaction(Row84, Freq84), interaction(Row85, Freq85), interaction(Row86, Freq86), interaction(Row87, Freq87), interaction(Row88, Freq88), interaction(Row89, Freq89), interaction(Row90, Freq90), interaction(Row91, Freq91), interaction(Row92, Freq92), interaction(Row93, Freq93), interaction(Row94, Freq94), interaction(Row95, Freq95), interaction(Row96, Freq96), interaction(Row97, Freq97), interaction(Row98, Freq98), interaction(Row99, Freq99), interaction(Row100, Freq100), interaction(Row101, Freq101), interaction(Row102, Freq102), interaction(Row103, Freq103), interaction(Row104, Freq104), interaction(Row105, Freq105), interaction(Row106, Freq106), interaction(Row107, Freq107), interaction(Row108, Freq108), interaction(Row109, Freq109), interaction(Row110, Freq110), interaction(Row111, Freq111), interaction(Row112, Freq112), interaction(Row113, Freq113), interaction(Row114, Freq114), interaction(Row115, Freq115), interaction(Row116, Freq116), interaction(Row117, Freq117), interaction(Row118, Freq118), interaction(Row119, Freq119), interaction(Row120, Freq120), interaction(Row121, Freq121), interaction(Row122, Freq122), interaction(Row123, Freq123), interaction(Row124, Freq124), interaction(Row125, Freq125), interaction(Row126, Freq126), interaction(Row127, Freq127), interaction(Row128, Freq128), interaction(Row129, Freq129), interaction(Row130, Freq130), interaction(Row131, Freq131), interaction(Row132, Freq132), interaction(Row133, Freq133), interaction(Row134, Freq134), interaction(Row135, Freq135), interaction(Row136, Freq136), interaction(Row137, Freq137), interaction(Row138, Freq138), interaction(Row139, Freq139), interaction(Row140, Freq140), interaction(Row141, Freq141), interaction(Row142, Freq142), interaction(Row143, Freq143), interaction(Row144, Freq144), interaction(Row145, Freq145), interaction(Row146, Freq146), interaction(Row147, Freq147), interaction(Row148, Freq148), interaction(Row149, Freq149), interaction(Row150, Freq150), interaction(Row151, Freq151), interaction(Row152, Freq152), interaction(Row153, Freq153), interaction(Row154, Freq154), interaction(Row155, Freq155), interaction(Row156, Freq156), interaction(Row157, Freq157), interaction(Row158, Freq158), interaction(Row159, Freq159), interaction(Row160, Freq160), interaction(Row161, Freq161), interaction(Row162, Freq162), interaction(Row163, Freq163), interaction(Row164, Freq164), interaction(Row165, Freq165), interaction(Row166, Freq166), interaction(Row167, Freq167), interaction(Row168, Freq168), interaction(Row169, Freq169), interaction(Row170, Freq170), interaction(Row171, Freq171), interaction(Row172, Freq172), interaction(Row173, Freq173), interaction(Row174, Freq174), interaction(Row175, Freq175), interaction(Row176, Freq176), interaction(Row177, Freq177), interaction(Row178, Freq178), interaction(Row179, Freq179), interaction(Row180, Freq180), interaction(Row181, Freq181), interaction(Row182, Freq182), interaction(Row183, Freq183), interaction(Row184, Freq184), interaction(Row185, Freq185), interaction(Row186, Freq186), interaction(Row187, Freq187), interaction(Row188, Freq188), interaction(Row189, Freq189), interaction(Row190, Freq190), interaction(Row191, Freq191), interaction(Row192, Freq192), interaction(Row193, Freq193), interaction(Row194, Freq194), interaction(Row195, Freq195), interaction(Row196, Freq196), interaction(Row197, Freq197), interaction(Row198, Freq198), interaction(Row199, Freq199), interaction(Row200, Freq200), interaction(Row201, Freq201), interaction(Row202, Freq202), interaction(Row203, Freq203), interaction(Row204, Freq204), interaction(Row205, Freq205), interaction(Row206, Freq206), interaction(Row207, Freq207), interaction(Row208, Freq208), interaction(Row209, Freq209), interaction(Row210, Freq210), interaction(Row211, Freq211), interaction(Row212, Freq212), interaction(Row213, Freq213), interaction(Row214, Freq214), interaction(Row215, Freq215), interaction(Row216, Freq216), interaction(Row217, Freq217), interaction(Row218, Freq218), interaction(Row219, Freq219), interaction(Row220, Freq220), interaction(Row221, Freq221), interaction(Row222, Freq222), interaction(Row223, Freq223), interaction(Row224, Freq224), interaction(Row225, Freq225), interaction(Row226, Freq226), interaction(Row227, Freq227), interaction(Row228, Freq228), interaction(Row229, Freq229), interaction(Row230, Freq230), interaction(Row231, Freq231), interaction(Row232, Freq232), interaction(Row233, Freq233), interaction(Row234, Freq234), interaction(Row235, Freq235), interaction(Row236, Freq236), interaction(Row237, Freq237), interaction(Row238, Freq238), interaction(Row239, Freq239), interaction(Row240, Freq240), interaction(Row241, Freq241), interaction(Row242, Freq242), interaction(Row243, Freq243), interaction(Row244, Freq244), interaction(Row245, Freq245), interaction(Row246, Freq246), interaction(Row247, Freq247), interaction(Row248, Freq248), interaction(Row249, Freq249), interaction(Row250, Freq250), interaction(Row251, Freq251), interaction(Row252, Freq252), interaction(Row253, Freq253), interaction(Row254, Freq254), interaction(Row255, Freq255), interaction(Row256, Freq256), interaction(Row257, Freq257), interaction(Row258, Freq258), interaction(Row259, Freq259), interaction(Row260, Freq260), interaction(Row261, Freq261), interaction(Row262, Freq262), interaction(Row263, Freq263), interaction(Row264, Freq264), interaction(Row265, Freq265), interaction(Row266, Freq266), interaction(Row267, Freq267), interaction(Row268, Freq268), interaction(Row269, Freq269), interaction(Row270, Freq270), interaction(Row271, Freq271), interaction(Row272, Freq272), interaction(Row273, Freq273), interaction(Row274, Freq274), interaction(Row275, Freq275), interaction(Row276, Freq276), interaction(Row277, Freq277), interaction(Row278, Freq278), interaction(Row279, Freq279), interaction(Row280, Freq280), interaction(Row281, Freq281), interaction(Row282, Freq282), interaction(Row283, Freq283), interaction(Row284, Freq284), interaction(Row285, Freq285), interaction(Row286, Freq286), interaction(Row287, Freq287), interaction(Row288, Freq288), interaction(Row289, Freq289), interaction(Row290, Freq290), interaction(Row291, Freq291), interaction(Row292, Freq292), interaction(Row293, Freq293), interaction(Row294, Freq294), interaction(Row295, Freq295), interaction(Row296, Freq296), interaction(Row297, Freq297), interaction(Row298, Freq298), interaction(Row299, Freq299), interaction(Row300, Freq300), interaction(Row301, Freq301), interaction(Row302, Freq302), interaction(Row303, Freq303), interaction(Row304, Freq304), interaction(Row305, Freq305), interaction(Row306, Freq306), interaction(Row307, Freq307), interaction(Row308, Freq308), interaction(Row309, Freq309), interaction(Row310, Freq310), interaction(Row311, Freq311), interaction(Row312, Freq312), interaction(Row313, Freq313), interaction(Row314, Freq314), interaction(Row315, Freq315), interaction(Row316, Freq316), interaction(Row317, Freq317), interaction(Row318, Freq318), interaction(Row319, Freq319), interaction(Row320, Freq320), interaction(Row321, Freq321), interaction(Row322, Freq322), interaction(Row323, Freq323), interaction(Row324, Freq324), interaction(Row325, Freq325), interaction(Row326, Freq326), interaction(Row327, Freq327), interaction(Row328, Freq328), interaction(Row329, Freq329), interaction(Row330, Freq330), interaction(Row331, Freq331), interaction(Row332, Freq332), interaction(Row333, Freq333), interaction(Row334, Freq334), interaction(Row335, Freq335), interaction(Row336, Freq336), interaction(Row337, Freq337), interaction(Row338, Freq338), interaction(Row339, Freq339), interaction(Row340, Freq340), interaction(Row341, Freq341), interaction(Row342, Freq342), interaction(Row343, Freq343), interaction(Row344, Freq344), interaction(Row345, Freq345), interaction(Row346, Freq346), interaction(Row347, Freq347), interaction(Row348, Freq348), interaction(Row349, Freq349), interaction(Row350, Freq350), interaction(Row351, Freq351), interaction(Row352, Freq352), interaction(Row353, Freq353), interaction(Row354, Freq354), interaction(Row355, Freq355), interaction(Row356, Freq356), interaction(Row357, Freq357), interaction(Row358, Freq358), interaction(Row359, Freq359), interaction(Row360, Freq360), interaction(Row361, Freq361), interaction(Row362, Freq362), interaction(Row363, Freq363), interaction(Row364, Freq364), interaction(Row365, Freq365), interaction(Row366, Freq366), interaction(Row367, Freq367), interaction(Row368, Freq368), interaction(Row369, Freq369), interaction(Row370, Freq370), interaction(Row371, Freq371), interaction(Row372, Freq372), interaction(Row373, Freq373), interaction(Row374, Freq374), interaction(Row375, Freq375), interaction(Row376, Freq376), interaction(Row377, Freq377), interaction(Row378, Freq378), interaction(Row379, Freq379), interaction(Row380, Freq380), interaction(Row381, Freq381), interaction(Row382, Freq382), interaction(Row383, Freq383), interaction(Row384, Freq384), interaction(Row385, Freq385), interaction(Row386, Freq386), interaction(Row387, Freq387), interaction(Row388, Freq388), interaction(Row389, Freq389), interaction(Row390, Freq390), interaction(Row391, Freq391), interaction(Row392, Freq392), interaction(Row393, Freq393), interaction(Row394, Freq394), interaction(Row395, Freq395), interaction(Row396, Freq396), interaction(Row397, Freq397), interaction(Row398, Freq398), interaction(Row399, Freq399), interaction(Row400, Freq400), interaction(Row401, Freq401), interaction(Row402, Freq402), interaction(Row403, Freq403), interaction(Row404, Freq404), interaction(Row405, Freq405), interaction(Row406, Freq406), interaction(Row407, Freq407), interaction(Row408, Freq408), interaction(Row409, Freq409), interaction(Row410, Freq410), interaction(Row411, Freq411), interaction(Row412, Freq412), interaction(Row413, Freq413), interaction(Row414, Freq414), interaction(Row415, Freq415), interaction(Row416, Freq416), interaction(Row417, Freq417), interaction(Row418, Freq418), interaction(Row419, Freq419), interaction(Row420, Freq420), interaction(Row421, Freq421), interaction(Row422, Freq422), interaction(Row423, Freq423), interaction(Row424, Freq424), interaction(Row425, Freq425), interaction(Row426, Freq426), interaction(Row427, Freq427), interaction(Row428, Freq428), interaction(Row429, Freq429), interaction(Row430, Freq430), interaction(Row431, Freq431), interaction(Row432, Freq432), interaction(Row433, Freq433), interaction(Row434, Freq434), interaction(Row435, Freq435), interaction(Row436, Freq436), interaction(Row437, Freq437), interaction(Row438, Freq438), interaction(Row439, Freq439), interaction(Row440, Freq440), interaction(Row441, Freq441), interaction(Row442, Freq442), interaction(Row443, Freq443), interaction(Row444, Freq444), interaction(Row445, Freq445), interaction(Row446, Freq446), interaction(Row447, Freq447), interaction(Row448, Freq448), interaction(Row449, Freq449), interaction(Row450, Freq450), interaction(Row451, Freq451), interaction(Row452, Freq452), interaction(Row453, Freq453), interaction(Row454, Freq454), interaction(Row455, Freq455), interaction(Row456, Freq456), interaction(Row457, Freq457), interaction(Row458, Freq458), interaction(Row459, Freq459), interaction(Row460, Freq460), interaction(Row461, Freq461), interaction(Row462, Freq462), interaction(Row463, Freq463), interaction(Row464, Freq464), interaction(Row465, Freq465), interaction(Row466, Freq466), interaction(Row467, Freq467), interaction(Row468, Freq468), interaction(Row469, Freq469), interaction(Row470, Freq470), interaction(Row471, Freq471), interaction(Row472, Freq472), interaction(Row473, Freq473), interaction(Row474, Freq474), interaction(Row475, Freq475), interaction(Row476, Freq476), interaction(Row477, Freq477), interaction(Row478, Freq478), interaction(Row479, Freq479), interaction(Row480, Freq480), interaction(Row481, Freq481), interaction(Row482, Freq482), interaction(Row483, Freq483), interaction(Row484, Freq484), interaction(Row485, Freq485), interaction(Row486, Freq486), interaction(Row487, Freq487), interaction(Row488, Freq488), interaction(Row489, Freq489), interaction(Row490, Freq490), interaction(Row491, Freq491), interaction(Row492, Freq492), interaction(Row493, Freq493), interaction(Row494, Freq494), interaction(Row495, Freq495), interaction(Row496, Freq496), interaction(Row497, Freq497), interaction(Row498, Freq498), interaction(Row499, Freq499), interaction(Row500, Freq500), interaction(Row501, Freq501), interaction(Row502, Freq502), interaction(Row503, Freq503), interaction(Row504, Freq504), interaction(Row505, Freq505), interaction(Row506, Freq506), interaction(Row507, Freq507), interaction(Row508, Freq508), interaction(Row509, Freq509), interaction(Row510, Freq510), interaction(Row511, Freq511), interaction(Row512, Freq512), interaction(Row513, Freq513), interaction(Row514, Freq514), interaction(Row515, Freq515), interaction(Row516, Freq516), interaction(Row517, Freq517), interaction(Row518, Freq518), interaction(Row519, Freq519), interaction(Row520, Freq520), interaction(Row521, Freq521), interaction(Row522, Freq522), interaction(Row523, Freq523), interaction(Row524, Freq524), interaction(Row525, Freq525), interaction(Row526, Freq526), interaction(Row527, Freq527), interaction(Row528, Freq528), interaction(Row529, Freq529), interaction(Row530, Freq530), interaction(Row531, Freq531), interaction(Row532, Freq532), interaction(Row533, Freq533), interaction(Row534, Freq534), interaction(Row535, Freq535), interaction(Row536, Freq536), interaction(Row537, Freq537), interaction(Row538, Freq538), interaction(Row539, Freq539), interaction(Row540, Freq540), interaction(Row541, Freq541), interaction(Row542, Freq542), interaction(Row543, Freq543), interaction(Row544, Freq544), interaction(Row545, Freq545), interaction(Row546, Freq546), interaction(Row547, Freq547), interaction(Row548, Freq548), interaction(Row549, Freq549), interaction(Row550, Freq550), interaction(Row551, Freq551), interaction(Row552, Freq552), interaction(Row553, Freq553), interaction(Row554, Freq554), interaction(Row555, Freq555), interaction(Row556, Freq556), interaction(Row557, Freq557), interaction(Row558, Freq558)],

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
	Row6 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 34, 426, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row7 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 424, 425, 426, 429, 430, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row8 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 431, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row9 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 431, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row10 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row11 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 420, 425, 426, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row12 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 28, 426, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row13 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 428, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row14 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 33, 426, 432, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row15 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row16 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 429, 430, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row17 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 430, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row18 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 425, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row19 :: [0, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 425, 428, 429, 430, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row20 :: [0, 5, 6, 7, 9, 10, 11, 12, 13, 14, 16, 17, 19, 420, 424, 425, 426, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row21 :: [0, 6, 12, 13, 14, 15, 17, 19, 425, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row22 :: [0, 11, 14, 18, 427, 429, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 449],
	Row23 :: [0, 6, 10, 13, 17, 434, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 447, 449],
	Row24 :: [0, 6, 10, 11, 13, 14, 17, 19, 403, 425, 430, 435, 436, 439, 440, 441, 442, 443, 444, 445, 447, 448, 449],
	Row25 :: [0, 13, 436, 439, 440, 441, 444, 445, 446, 450],
	Row26 :: [0, 436, 438, 440, 441, 442, 444, 446],
	Row27 :: [0, 19, 441],
	Row28 :: [0, 435],
	Row29 :: [0, 440, 442],
	Row30 :: [0],
	Row31 :: [0],
	Row32 :: [0, 442],
	Row33 :: [0],
	Row34 :: [0],
	Row35 :: [0, 441],
	Row36 :: [0],
	Row37 :: [0],
	Row38 :: [0],
	Row39 :: [0],
	Row40 :: [0],
	Row41 :: [0],
	Row42 :: [0],
	Row43 :: [0],
	Row44 :: [0],
	Row45 :: [0],
	Row46 :: [0],
	Row47 :: [0],
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
	Row58 :: [0],
	Row59 :: [0, 426],
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
	Row94 :: [0, 111],
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
	Row107 :: [0],
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
	Row119 :: [0],
	Row120 :: [0],
	Row121 :: [0],
	Row122 :: [0],
	Row123 :: [0],
	Row124 :: [0],
	Row125 :: [0],
	Row126 :: [0],
	Row127 :: [0],
	Row128 :: [0],
	Row129 :: [0],
	Row130 :: [0],
	Row131 :: [0],
	Row132 :: [0],
	Row133 :: [0],
	Row134 :: [0],
	Row135 :: [0],
	Row136 :: [0],
	Row137 :: [0],
	Row138 :: [0],
	Row139 :: [0],
	Row140 :: [0],
	Row141 :: [0],
	Row142 :: [0],
	Row143 :: [0],
	Row144 :: [0],
	Row145 :: [0],
	Row146 :: [0],
	Row147 :: [0],
	Row148 :: [0],
	Row149 :: [0],
	Row150 :: [0],
	Row151 :: [0],
	Row152 :: [0],
	Row153 :: [0],
	Row154 :: [0],
	Row155 :: [0],
	Row156 :: [0],
	Row157 :: [0],
	Row158 :: [0],
	Row159 :: [0],
	Row160 :: [0],
	Row161 :: [0],
	Row162 :: [0],
	Row163 :: [0],
	Row164 :: [0],
	Row165 :: [0],
	Row166 :: [0],
	Row167 :: [0],
	Row168 :: [0],
	Row169 :: [0],
	Row170 :: [0],
	Row171 :: [0],
	Row172 :: [0],
	Row173 :: [0],
	Row174 :: [0],
	Row175 :: [0],
	Row176 :: [0],
	Row177 :: [0],
	Row178 :: [0],
	Row179 :: [0],
	Row180 :: [0],
	Row181 :: [0],
	Row182 :: [0],
	Row183 :: [0],
	Row184 :: [0],
	Row185 :: [0],
	Row186 :: [0, 244],
	Row187 :: [0],
	Row188 :: [0],
	Row189 :: [0],
	Row190 :: [0],
	Row191 :: [0],
	Row192 :: [0],
	Row193 :: [0],
	Row194 :: [0],
	Row195 :: [0],
	Row196 :: [0],
	Row197 :: [0],
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
	Row265 :: [0],
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
	Row319 :: [0],
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
	Row340 :: [0],
	Row341 :: [0],
	Row342 :: [0],
	Row343 :: [0],
	Row344 :: [0],
	Row345 :: [0],
	Row346 :: [0],
	Row347 :: [0],
	Row348 :: [0, 176],
	Row349 :: [0, 146, 168, 182],
	Row350 :: [0, 144, 178],
	Row351 :: [0],
	Row352 :: [0, 158, 173],
	Row353 :: [0, 154, 168, 169, 176],
	Row354 :: [0, 167, 169],
	Row355 :: [0, 153, 155, 166, 176],
	Row356 :: [0, 152, 154, 155, 158, 166, 167, 168, 169, 172, 176, 177],
	Row357 :: [0, 158, 159, 166, 167, 168, 169, 177],
	Row358 :: [0, 146, 147, 151, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 175, 176, 177, 178, 179, 180, 181, 185],
	Row359 :: [0, 151, 152, 153, 154, 155, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 175, 176, 177],
	Row360 :: [0, 151, 152, 153, 155, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 178, 179, 182],
	Row361 :: [0, 149, 151, 152, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 173, 174, 175, 176, 177],
	Row362 :: [0, 149, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 182],
	Row363 :: [0, 150, 151, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 179, 184],
	Row364 :: [0, 146, 149, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 187],
	Row365 :: [0, 147, 148, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182],
	Row366 :: [0, 142, 145, 147, 148, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184],
	Row367 :: [0, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 185],
	Row368 :: [0, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 184, 185],
	Row369 :: [0, 144, 145, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183],
	Row370 :: [0, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183],
	Row371 :: [0, 146, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 188],
	Row372 :: [0, 142, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182],
	Row373 :: [0, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 185],
	Row374 :: [0, 144, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 183, 186, 187, 209],
	Row375 :: [0],
	Row376 :: [0],
	Row377 :: [0],
	Row378 :: [0],
	Row379 :: [0],
	Row380 :: [0, 148, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 184],
	Row381 :: [0, 145, 146, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 187],
	Row382 :: [0, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181],
	Row383 :: [0, 146, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 182],
	Row384 :: [0, 146, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 183, 185],
	Row385 :: [0, 144, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183],
	Row386 :: [0, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182],
	Row387 :: [0, 149, 151, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 187],
	Row388 :: [0, 142, 144, 150, 151, 152, 153, 154, 155, 156, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 182, 183],
	Row389 :: [0, 150, 152, 153, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 182, 183],
	Row390 :: [0, 152, 154, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 183, 186],
	Row391 :: [0, 145, 149, 151, 152, 153, 154, 155, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 185, 186],
	Row392 :: [0, 151, 153, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 173, 174, 176, 177, 178, 180, 181, 183],
	Row393 :: [0, 150, 155, 156, 157, 158, 166, 167, 168, 169, 170, 171, 172, 173, 175, 176, 177, 178, 179, 180, 182, 196],
	Row394 :: [0, 153, 155, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 183, 188],
	Row395 :: [0, 151, 155, 156, 157, 158, 159, 166, 167, 168, 169, 170, 171, 172, 174, 175, 176, 177, 178],
	Row396 :: [0, 152, 154, 156, 157, 158, 166, 167, 168, 169, 171, 172, 173, 175, 176, 177, 178, 179, 180, 181, 183],
	Row397 :: [0, 152, 157, 158, 167, 169, 173, 175, 177, 179, 180, 186],
	Row398 :: [0, 157, 169],
	Row399 :: [0, 157, 158, 167, 169, 172, 175, 176, 181],
	Row400 :: [0, 167, 175, 182],
	Row401 :: [0, 153, 169, 175, 176, 181, 182],
	Row402 :: [0, 176],
	Row403 :: [0, 195],
	Row404 :: [0],
	Row405 :: [0],
	Row406 :: [0],
	Row407 :: [0],
	Row408 :: [0],
	Row409 :: [0],
	Row410 :: [0],
	Row411 :: [0],
	Row412 :: [0],
	Row413 :: [0],
	Row414 :: [0],
	Row415 :: [0],
	Row416 :: [0],
	Row417 :: [0],
	Row418 :: [0],
	Row419 :: [0],
	Row420 :: [0],
	Row421 :: [0],
	Row422 :: [0],
	Row423 :: [0],
	Row424 :: [0, 151],
	Row425 :: [0],
	Row426 :: [0],
	Row427 :: [0],
	Row428 :: [0],
	Row429 :: [0],
	Row430 :: [0],
	Row431 :: [0],
	Row432 :: [0],
	Row433 :: [0],
	Row434 :: [0],
	Row435 :: [0],
	Row436 :: [0],
	Row437 :: [0],
	Row438 :: [0],
	Row439 :: [0],
	Row440 :: [0],
	Row441 :: [0],
	Row442 :: [0],
	Row443 :: [0],
	Row444 :: [0],
	Row445 :: [0],
	Row446 :: [0],
	Row447 :: [0],
	Row448 :: [0],
	Row449 :: [0],
	Row450 :: [0],
	Row451 :: [0],
	Row452 :: [0],
	Row453 :: [0],
	Row454 :: [0],
	Row455 :: [0],
	Row456 :: [0],
	Row457 :: [0],
	Row458 :: [0],
	Row459 :: [0],
	Row460 :: [0],
	Row461 :: [0],
	Row462 :: [0, 4],
	Row463 :: [0],
	Row464 :: [0],
	Row465 :: [0],
	Row466 :: [0],
	Row467 :: [0],
	Row468 :: [0],
	Row469 :: [0],
	Row470 :: [0],
	Row471 :: [0],
	Row472 :: [0],
	Row473 :: [0],
	Row474 :: [0],
	Row475 :: [0],
	Row476 :: [0],
	Row477 :: [0],
	Row478 :: [0],
	Row479 :: [0],
	Row480 :: [0],
	Row481 :: [0],
	Row482 :: [0],
	Row483 :: [0],
	Row484 :: [0],
	Row485 :: [0],
	Row486 :: [0],
	Row487 :: [0],
	Row488 :: [0],
	Row489 :: [0],
	Row490 :: [0],
	Row491 :: [0],
	Row492 :: [0],
	Row493 :: [0],
	Row494 :: [0],
	Row495 :: [0],
	Row496 :: [0],
	Row497 :: [0],
	Row498 :: [0],
	Row499 :: [0],
	Row500 :: [0],
	Row501 :: [0],
	Row502 :: [0],
	Row503 :: [0],
	Row504 :: [0],
	Row505 :: [0],
	Row506 :: [0],
	Row507 :: [0],
	Row508 :: [0],
	Row509 :: [0],
	Row510 :: [0],
	Row511 :: [0],
	Row512 :: [0],
	Row513 :: [0],
	Row514 :: [0],
	Row515 :: [0],
	Row516 :: [0],
	Row517 :: [0],
	Row518 :: [0],
	Row519 :: [0],
	Row520 :: [0],
	Row521 :: [0],
	Row522 :: [0],
	Row523 :: [0],
	Row524 :: [0],
	Row525 :: [0],
	Row526 :: [0],
	Row527 :: [0],
	Row528 :: [0],
	Row529 :: [0],
	Row530 :: [0],
	Row531 :: [0, 11],
	Row532 :: [0],
	Row533 :: [0, 11],
	Row534 :: [0],
	Row535 :: [0, 11, 12],
	Row536 :: [0, 10, 19],
	Row537 :: [0, 6, 14, 450],
	Row538 :: [0, 6, 7, 11, 14, 15, 440, 442, 445, 446, 449],
	Row539 :: [0, 8, 13, 14, 16, 18, 438, 439, 440, 441],
	Row540 :: [0, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 18, 20, 439, 440, 441, 442, 443, 444, 445, 446, 448, 449],
	Row541 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 20, 436, 439, 440, 441, 442, 443, 444, 445, 447, 449, 450],
	Row542 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 439, 440, 441, 442, 443, 444, 445, 447, 448, 449, 450],
	Row543 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 437, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row544 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 435, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row545 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 433, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row546 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row547 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row548 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row549 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row550 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 433, 434, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row551 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row552 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row553 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 21, 214, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row554 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 213, 214, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row555 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 213, 214, 215, 216, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row556 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 21, 211, 212, 213, 214, 215, 216, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row557 :: [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 23, 166, 210, 211, 212, 213, 214, 215, 216, 434, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450],
	Row558 :: [0],
	% Each frequency term can assume either the rounded 
	% and scaled integer value (based on the corresponding 
	% interaction frequency from the whole-genome contact 
	% map) or a value of 0 where 0  represents an 
	% interaction not being selected 
	Freq1 :: [0],
	Freq2 :: [0],
	Freq3 :: [0],
	Freq4 :: [0],
	Freq5 :: [0],
	Freq6 :: [0, 7, 8, 6, 4, 3, 2, 1, 5, 9, 10, 11],
	Freq7 :: [0, 9, 11, 8, 7, 5, 4, 3, 2, 1, 10, 12, 13],
	Freq8 :: [0, 12, 10, 8, 5, 4, 6, 3, 2, 1, 13, 14, 16],
	Freq9 :: [0, 10, 12, 9, 11, 8, 5, 3, 6, 4, 2, 1, 13, 15],
	Freq10 :: [0, 6, 7, 8, 5, 4, 3, 2, 1, 9, 10],
	Freq11 :: [0, 3, 4, 2, 1, 5],
	Freq12 :: [0, 2, 3, 4, 1],
	Freq13 :: [0, 3, 2, 1, 4],
	Freq14 :: [0, 2, 3, 1, 4],
	Freq15 :: [0, 2, 3, 1],
	Freq16 :: [0, 1, 2, 3],
	Freq17 :: [0, 1, 2],
	Freq18 :: [0, 1, 2],
	Freq19 :: [0, 1],
	Freq20 :: [0, 1, 2],
	Freq21 :: [0, 1],
	Freq22 :: [0, 1],
	Freq23 :: [0, 1],
	Freq24 :: [0, 1],
	Freq25 :: [0, 1],
	Freq26 :: [0, 1],
	Freq27 :: [0, 1],
	Freq28 :: [0, 1],
	Freq29 :: [0, 1],
	Freq30 :: [0],
	Freq31 :: [0],
	Freq32 :: [0, 1],
	Freq33 :: [0],
	Freq34 :: [0],
	Freq35 :: [0, 1],
	Freq36 :: [0],
	Freq37 :: [0],
	Freq38 :: [0],
	Freq39 :: [0],
	Freq40 :: [0],
	Freq41 :: [0],
	Freq42 :: [0],
	Freq43 :: [0],
	Freq44 :: [0],
	Freq45 :: [0],
	Freq46 :: [0],
	Freq47 :: [0],
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
	Freq58 :: [0],
	Freq59 :: [0, 1],
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
	Freq94 :: [0, 1],
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
	Freq107 :: [0],
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
	Freq119 :: [0],
	Freq120 :: [0],
	Freq121 :: [0],
	Freq122 :: [0],
	Freq123 :: [0],
	Freq124 :: [0],
	Freq125 :: [0],
	Freq126 :: [0],
	Freq127 :: [0],
	Freq128 :: [0],
	Freq129 :: [0],
	Freq130 :: [0],
	Freq131 :: [0],
	Freq132 :: [0],
	Freq133 :: [0],
	Freq134 :: [0],
	Freq135 :: [0],
	Freq136 :: [0],
	Freq137 :: [0],
	Freq138 :: [0],
	Freq139 :: [0],
	Freq140 :: [0],
	Freq141 :: [0],
	Freq142 :: [0],
	Freq143 :: [0],
	Freq144 :: [0],
	Freq145 :: [0],
	Freq146 :: [0],
	Freq147 :: [0],
	Freq148 :: [0],
	Freq149 :: [0],
	Freq150 :: [0],
	Freq151 :: [0],
	Freq152 :: [0],
	Freq153 :: [0],
	Freq154 :: [0],
	Freq155 :: [0],
	Freq156 :: [0],
	Freq157 :: [0],
	Freq158 :: [0],
	Freq159 :: [0],
	Freq160 :: [0],
	Freq161 :: [0],
	Freq162 :: [0],
	Freq163 :: [0],
	Freq164 :: [0],
	Freq165 :: [0],
	Freq166 :: [0],
	Freq167 :: [0],
	Freq168 :: [0],
	Freq169 :: [0],
	Freq170 :: [0],
	Freq171 :: [0],
	Freq172 :: [0],
	Freq173 :: [0],
	Freq174 :: [0],
	Freq175 :: [0],
	Freq176 :: [0],
	Freq177 :: [0],
	Freq178 :: [0],
	Freq179 :: [0],
	Freq180 :: [0],
	Freq181 :: [0],
	Freq182 :: [0],
	Freq183 :: [0],
	Freq184 :: [0],
	Freq185 :: [0],
	Freq186 :: [0, 1],
	Freq187 :: [0],
	Freq188 :: [0],
	Freq189 :: [0],
	Freq190 :: [0],
	Freq191 :: [0],
	Freq192 :: [0],
	Freq193 :: [0],
	Freq194 :: [0],
	Freq195 :: [0],
	Freq196 :: [0],
	Freq197 :: [0],
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
	Freq265 :: [0],
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
	Freq319 :: [0],
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
	Freq340 :: [0],
	Freq341 :: [0],
	Freq342 :: [0],
	Freq343 :: [0],
	Freq344 :: [0],
	Freq345 :: [0],
	Freq346 :: [0],
	Freq347 :: [0],
	Freq348 :: [0, 1],
	Freq349 :: [0, 1],
	Freq350 :: [0, 1],
	Freq351 :: [0],
	Freq352 :: [0, 1],
	Freq353 :: [0, 1],
	Freq354 :: [0, 1],
	Freq355 :: [0, 1],
	Freq356 :: [0, 1],
	Freq357 :: [0, 1],
	Freq358 :: [0, 1],
	Freq359 :: [0, 1],
	Freq360 :: [0, 1],
	Freq361 :: [0, 1],
	Freq362 :: [0, 1],
	Freq363 :: [0, 1],
	Freq364 :: [0, 1, 2],
	Freq365 :: [0, 1, 2],
	Freq366 :: [0, 1, 2],
	Freq367 :: [0, 1, 2],
	Freq368 :: [0, 1, 2, 3],
	Freq369 :: [0, 1, 2, 3],
	Freq370 :: [0, 1, 2, 3],
	Freq371 :: [0, 1, 2, 4, 3],
	Freq372 :: [0, 1, 2, 3, 4],
	Freq373 :: [0, 1, 2, 3, 5, 4],
	Freq374 :: [0, 1, 2, 3, 4, 6, 5],
	Freq375 :: [0],
	Freq376 :: [0],
	Freq377 :: [0],
	Freq378 :: [0],
	Freq379 :: [0],
	Freq380 :: [0, 1, 2, 3, 4, 6],
	Freq381 :: [0, 1, 2, 3, 4, 6, 5],
	Freq382 :: [0, 1, 2, 4, 5, 3],
	Freq383 :: [0, 1, 2, 3, 4, 5],
	Freq384 :: [0, 1, 2, 3, 4],
	Freq385 :: [0, 1, 2, 3, 4],
	Freq386 :: [0, 1, 2, 3],
	Freq387 :: [0, 1, 2, 3],
	Freq388 :: [0, 1, 2],
	Freq389 :: [0, 1, 2],
	Freq390 :: [0, 1, 2],
	Freq391 :: [0, 1, 2],
	Freq392 :: [0, 1],
	Freq393 :: [0, 1],
	Freq394 :: [0, 1],
	Freq395 :: [0, 1],
	Freq396 :: [0, 1],
	Freq397 :: [0, 1],
	Freq398 :: [0, 1],
	Freq399 :: [0, 1],
	Freq400 :: [0, 1],
	Freq401 :: [0, 1],
	Freq402 :: [0, 1],
	Freq403 :: [0, 1],
	Freq404 :: [0],
	Freq405 :: [0],
	Freq406 :: [0],
	Freq407 :: [0],
	Freq408 :: [0],
	Freq409 :: [0],
	Freq410 :: [0],
	Freq411 :: [0],
	Freq412 :: [0],
	Freq413 :: [0],
	Freq414 :: [0],
	Freq415 :: [0],
	Freq416 :: [0],
	Freq417 :: [0],
	Freq418 :: [0],
	Freq419 :: [0],
	Freq420 :: [0],
	Freq421 :: [0],
	Freq422 :: [0],
	Freq423 :: [0],
	Freq424 :: [0, 1],
	Freq425 :: [0],
	Freq426 :: [0],
	Freq427 :: [0],
	Freq428 :: [0],
	Freq429 :: [0],
	Freq430 :: [0],
	Freq431 :: [0],
	Freq432 :: [0],
	Freq433 :: [0],
	Freq434 :: [0],
	Freq435 :: [0],
	Freq436 :: [0],
	Freq437 :: [0],
	Freq438 :: [0],
	Freq439 :: [0],
	Freq440 :: [0],
	Freq441 :: [0],
	Freq442 :: [0],
	Freq443 :: [0],
	Freq444 :: [0],
	Freq445 :: [0],
	Freq446 :: [0],
	Freq447 :: [0],
	Freq448 :: [0],
	Freq449 :: [0],
	Freq450 :: [0],
	Freq451 :: [0],
	Freq452 :: [0],
	Freq453 :: [0],
	Freq454 :: [0],
	Freq455 :: [0],
	Freq456 :: [0],
	Freq457 :: [0],
	Freq458 :: [0],
	Freq459 :: [0],
	Freq460 :: [0],
	Freq461 :: [0],
	Freq462 :: [0, 7],
	Freq463 :: [0],
	Freq464 :: [0],
	Freq465 :: [0],
	Freq466 :: [0],
	Freq467 :: [0],
	Freq468 :: [0],
	Freq469 :: [0],
	Freq470 :: [0],
	Freq471 :: [0],
	Freq472 :: [0],
	Freq473 :: [0],
	Freq474 :: [0],
	Freq475 :: [0],
	Freq476 :: [0],
	Freq477 :: [0],
	Freq478 :: [0],
	Freq479 :: [0],
	Freq480 :: [0],
	Freq481 :: [0],
	Freq482 :: [0],
	Freq483 :: [0],
	Freq484 :: [0],
	Freq485 :: [0],
	Freq486 :: [0],
	Freq487 :: [0],
	Freq488 :: [0],
	Freq489 :: [0],
	Freq490 :: [0],
	Freq491 :: [0],
	Freq492 :: [0],
	Freq493 :: [0],
	Freq494 :: [0],
	Freq495 :: [0],
	Freq496 :: [0],
	Freq497 :: [0],
	Freq498 :: [0],
	Freq499 :: [0],
	Freq500 :: [0],
	Freq501 :: [0],
	Freq502 :: [0],
	Freq503 :: [0],
	Freq504 :: [0],
	Freq505 :: [0],
	Freq506 :: [0],
	Freq507 :: [0],
	Freq508 :: [0],
	Freq509 :: [0],
	Freq510 :: [0],
	Freq511 :: [0],
	Freq512 :: [0],
	Freq513 :: [0],
	Freq514 :: [0],
	Freq515 :: [0],
	Freq516 :: [0],
	Freq517 :: [0],
	Freq518 :: [0],
	Freq519 :: [0],
	Freq520 :: [0],
	Freq521 :: [0],
	Freq522 :: [0],
	Freq523 :: [0],
	Freq524 :: [0],
	Freq525 :: [0],
	Freq526 :: [0],
	Freq527 :: [0],
	Freq528 :: [0],
	Freq529 :: [0],
	Freq530 :: [0],
	Freq531 :: [0, 1],
	Freq532 :: [0],
	Freq533 :: [0, 1],
	Freq534 :: [0],
	Freq535 :: [0, 1],
	Freq536 :: [0, 1],
	Freq537 :: [0, 1],
	Freq538 :: [0, 1],
	Freq539 :: [0, 1],
	Freq540 :: [0, 1],
	Freq541 :: [0, 1],
	Freq542 :: [0, 1],
	Freq543 :: [0, 1, 2],
	Freq544 :: [0, 1, 2],
	Freq545 :: [0, 2, 3, 1],
	Freq546 :: [0, 2, 3, 1],
	Freq547 :: [0, 3, 4, 2, 1],
	Freq548 :: [0, 4, 5, 3, 2, 1],
	Freq549 :: [0, 4, 2, 3, 1, 5],
	Freq550 :: [0, 6, 7, 5, 4, 3, 2, 1],
	Freq551 :: [0, 9, 11, 8, 10, 7, 5, 3, 4, 2, 1],
	Freq552 :: [0, 11, 8, 10, 7, 4, 2, 1, 3, 9],
	Freq553 :: [0, 12, 9, 8, 6, 3, 4, 2, 1, 10, 11],
	Freq554 :: [0, 10, 7, 8, 6, 5, 2, 3, 1, 4, 9, 11],
	Freq555 :: [0, 7, 8, 6, 5, 4, 3, 2, 1],
	Freq556 :: [0, 9, 8, 5, 6, 3, 2, 1, 4, 7],
	Freq557 :: [0, 8, 6, 5, 4, 3, 2, 1, 7],
	Freq558 :: [0],


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

	((Row6 #= 4 and Freq6 #= 7) or
	(Row6 #= 5 and Freq6 #= 8) or
	(Row6 #= 6 and Freq6 #= 8) or
	(Row6 #= 7 and Freq6 #= 8) or
	(Row6 #= 8 and Freq6 #= 7) or
	(Row6 #= 9 and Freq6 #= 6) or
	(Row6 #= 10 and Freq6 #= 4) or
	(Row6 #= 11 and Freq6 #= 4) or
	(Row6 #= 12 and Freq6 #= 4) or
	(Row6 #= 13 and Freq6 #= 4) or
	(Row6 #= 14 and Freq6 #= 3) or
	(Row6 #= 15 and Freq6 #= 3) or
	(Row6 #= 16 and Freq6 #= 3) or
	(Row6 #= 17 and Freq6 #= 2) or
	(Row6 #= 18 and Freq6 #= 1) or
	(Row6 #= 19 and Freq6 #= 1) or
	(Row6 #= 20 and Freq6 #= 1) or
	(Row6 #= 21 and Freq6 #= 2) or
	(Row6 #= 22 and Freq6 #= 1) or
	(Row6 #= 23 and Freq6 #= 1) or
	(Row6 #= 24 and Freq6 #= 1) or
	(Row6 #= 26 and Freq6 #= 1) or
	(Row6 #= 34 and Freq6 #= 1) or
	(Row6 #= 426 and Freq6 #= 1) or
	(Row6 #= 429 and Freq6 #= 1) or
	(Row6 #= 430 and Freq6 #= 1) or
	(Row6 #= 431 and Freq6 #= 1) or
	(Row6 #= 432 and Freq6 #= 1) or
	(Row6 #= 433 and Freq6 #= 1) or
	(Row6 #= 434 and Freq6 #= 1) or
	(Row6 #= 435 and Freq6 #= 1) or
	(Row6 #= 436 and Freq6 #= 1) or
	(Row6 #= 437 and Freq6 #= 2) or
	(Row6 #= 438 and Freq6 #= 2) or
	(Row6 #= 439 and Freq6 #= 2) or
	(Row6 #= 440 and Freq6 #= 3) or
	(Row6 #= 441 and Freq6 #= 4) or
	(Row6 #= 442 and Freq6 #= 6) or
	(Row6 #= 443 and Freq6 #= 5) or
	(Row6 #= 444 and Freq6 #= 6) or
	(Row6 #= 445 and Freq6 #= 6) or
	(Row6 #= 446 and Freq6 #= 9) or
	(Row6 #= 447 and Freq6 #= 10) or
	(Row6 #= 448 and Freq6 #= 10) or
	(Row6 #= 449 and Freq6 #= 11) or
	(Row6 #= 450 and Freq6 #= 11) or
	(Row6 #= 0 and Freq6 #= 0)), 

	((Row7 #= 4 and Freq7 #= 9) or
	(Row7 #= 5 and Freq7 #= 11) or
	(Row7 #= 6 and Freq7 #= 9) or
	(Row7 #= 7 and Freq7 #= 9) or
	(Row7 #= 8 and Freq7 #= 8) or
	(Row7 #= 9 and Freq7 #= 7) or
	(Row7 #= 10 and Freq7 #= 5) or
	(Row7 #= 11 and Freq7 #= 4) or
	(Row7 #= 12 and Freq7 #= 5) or
	(Row7 #= 13 and Freq7 #= 5) or
	(Row7 #= 14 and Freq7 #= 4) or
	(Row7 #= 15 and Freq7 #= 3) or
	(Row7 #= 16 and Freq7 #= 3) or
	(Row7 #= 17 and Freq7 #= 2) or
	(Row7 #= 18 and Freq7 #= 1) or
	(Row7 #= 19 and Freq7 #= 1) or
	(Row7 #= 20 and Freq7 #= 1) or
	(Row7 #= 21 and Freq7 #= 1) or
	(Row7 #= 22 and Freq7 #= 1) or
	(Row7 #= 23 and Freq7 #= 1) or
	(Row7 #= 424 and Freq7 #= 1) or
	(Row7 #= 425 and Freq7 #= 1) or
	(Row7 #= 426 and Freq7 #= 1) or
	(Row7 #= 429 and Freq7 #= 1) or
	(Row7 #= 430 and Freq7 #= 1) or
	(Row7 #= 432 and Freq7 #= 1) or
	(Row7 #= 433 and Freq7 #= 1) or
	(Row7 #= 434 and Freq7 #= 1) or
	(Row7 #= 435 and Freq7 #= 1) or
	(Row7 #= 436 and Freq7 #= 1) or
	(Row7 #= 437 and Freq7 #= 2) or
	(Row7 #= 438 and Freq7 #= 2) or
	(Row7 #= 439 and Freq7 #= 2) or
	(Row7 #= 440 and Freq7 #= 2) or
	(Row7 #= 441 and Freq7 #= 4) or
	(Row7 #= 442 and Freq7 #= 5) or
	(Row7 #= 443 and Freq7 #= 5) or
	(Row7 #= 444 and Freq7 #= 5) or
	(Row7 #= 445 and Freq7 #= 5) or
	(Row7 #= 446 and Freq7 #= 10) or
	(Row7 #= 447 and Freq7 #= 11) or
	(Row7 #= 448 and Freq7 #= 11) or
	(Row7 #= 449 and Freq7 #= 12) or
	(Row7 #= 450 and Freq7 #= 13) or
	(Row7 #= 0 and Freq7 #= 0)), 

	((Row8 #= 4 and Freq8 #= 12) or
	(Row8 #= 5 and Freq8 #= 12) or
	(Row8 #= 6 and Freq8 #= 10) or
	(Row8 #= 7 and Freq8 #= 12) or
	(Row8 #= 8 and Freq8 #= 10) or
	(Row8 #= 9 and Freq8 #= 8) or
	(Row8 #= 10 and Freq8 #= 5) or
	(Row8 #= 11 and Freq8 #= 4) or
	(Row8 #= 12 and Freq8 #= 6) or
	(Row8 #= 13 and Freq8 #= 5) or
	(Row8 #= 14 and Freq8 #= 3) or
	(Row8 #= 15 and Freq8 #= 3) or
	(Row8 #= 16 and Freq8 #= 3) or
	(Row8 #= 17 and Freq8 #= 2) or
	(Row8 #= 18 and Freq8 #= 1) or
	(Row8 #= 19 and Freq8 #= 1) or
	(Row8 #= 20 and Freq8 #= 1) or
	(Row8 #= 21 and Freq8 #= 1) or
	(Row8 #= 22 and Freq8 #= 1) or
	(Row8 #= 23 and Freq8 #= 1) or
	(Row8 #= 431 and Freq8 #= 1) or
	(Row8 #= 433 and Freq8 #= 1) or
	(Row8 #= 434 and Freq8 #= 1) or
	(Row8 #= 435 and Freq8 #= 1) or
	(Row8 #= 436 and Freq8 #= 1) or
	(Row8 #= 437 and Freq8 #= 1) or
	(Row8 #= 438 and Freq8 #= 2) or
	(Row8 #= 439 and Freq8 #= 2) or
	(Row8 #= 440 and Freq8 #= 2) or
	(Row8 #= 441 and Freq8 #= 4) or
	(Row8 #= 442 and Freq8 #= 5) or
	(Row8 #= 443 and Freq8 #= 5) or
	(Row8 #= 444 and Freq8 #= 5) or
	(Row8 #= 445 and Freq8 #= 5) or
	(Row8 #= 446 and Freq8 #= 12) or
	(Row8 #= 447 and Freq8 #= 13) or
	(Row8 #= 448 and Freq8 #= 13) or
	(Row8 #= 449 and Freq8 #= 14) or
	(Row8 #= 450 and Freq8 #= 16) or
	(Row8 #= 0 and Freq8 #= 0)), 

	((Row9 #= 4 and Freq9 #= 10) or
	(Row9 #= 5 and Freq9 #= 12) or
	(Row9 #= 6 and Freq9 #= 9) or
	(Row9 #= 7 and Freq9 #= 11) or
	(Row9 #= 8 and Freq9 #= 10) or
	(Row9 #= 9 and Freq9 #= 8) or
	(Row9 #= 10 and Freq9 #= 5) or
	(Row9 #= 11 and Freq9 #= 3) or
	(Row9 #= 12 and Freq9 #= 5) or
	(Row9 #= 13 and Freq9 #= 6) or
	(Row9 #= 14 and Freq9 #= 4) or
	(Row9 #= 15 and Freq9 #= 4) or
	(Row9 #= 16 and Freq9 #= 3) or
	(Row9 #= 17 and Freq9 #= 2) or
	(Row9 #= 18 and Freq9 #= 1) or
	(Row9 #= 19 and Freq9 #= 1) or
	(Row9 #= 20 and Freq9 #= 1) or
	(Row9 #= 21 and Freq9 #= 1) or
	(Row9 #= 22 and Freq9 #= 1) or
	(Row9 #= 431 and Freq9 #= 1) or
	(Row9 #= 434 and Freq9 #= 1) or
	(Row9 #= 435 and Freq9 #= 1) or
	(Row9 #= 436 and Freq9 #= 1) or
	(Row9 #= 437 and Freq9 #= 1) or
	(Row9 #= 438 and Freq9 #= 2) or
	(Row9 #= 439 and Freq9 #= 2) or
	(Row9 #= 440 and Freq9 #= 2) or
	(Row9 #= 441 and Freq9 #= 3) or
	(Row9 #= 442 and Freq9 #= 5) or
	(Row9 #= 443 and Freq9 #= 5) or
	(Row9 #= 444 and Freq9 #= 5) or
	(Row9 #= 445 and Freq9 #= 5) or
	(Row9 #= 446 and Freq9 #= 11) or
	(Row9 #= 447 and Freq9 #= 13) or
	(Row9 #= 448 and Freq9 #= 15) or
	(Row9 #= 449 and Freq9 #= 13) or
	(Row9 #= 450 and Freq9 #= 13) or
	(Row9 #= 0 and Freq9 #= 0)), 

	((Row10 #= 4 and Freq10 #= 6) or
	(Row10 #= 5 and Freq10 #= 7) or
	(Row10 #= 6 and Freq10 #= 8) or
	(Row10 #= 7 and Freq10 #= 8) or
	(Row10 #= 8 and Freq10 #= 7) or
	(Row10 #= 9 and Freq10 #= 6) or
	(Row10 #= 10 and Freq10 #= 5) or
	(Row10 #= 11 and Freq10 #= 4) or
	(Row10 #= 12 and Freq10 #= 4) or
	(Row10 #= 13 and Freq10 #= 4) or
	(Row10 #= 14 and Freq10 #= 3) or
	(Row10 #= 15 and Freq10 #= 3) or
	(Row10 #= 16 and Freq10 #= 2) or
	(Row10 #= 17 and Freq10 #= 2) or
	(Row10 #= 18 and Freq10 #= 1) or
	(Row10 #= 19 and Freq10 #= 1) or
	(Row10 #= 20 and Freq10 #= 1) or
	(Row10 #= 21 and Freq10 #= 1) or
	(Row10 #= 22 and Freq10 #= 1) or
	(Row10 #= 23 and Freq10 #= 1) or
	(Row10 #= 432 and Freq10 #= 1) or
	(Row10 #= 433 and Freq10 #= 1) or
	(Row10 #= 434 and Freq10 #= 1) or
	(Row10 #= 435 and Freq10 #= 1) or
	(Row10 #= 436 and Freq10 #= 1) or
	(Row10 #= 437 and Freq10 #= 1) or
	(Row10 #= 438 and Freq10 #= 2) or
	(Row10 #= 439 and Freq10 #= 2) or
	(Row10 #= 440 and Freq10 #= 2) or
	(Row10 #= 441 and Freq10 #= 4) or
	(Row10 #= 442 and Freq10 #= 5) or
	(Row10 #= 443 and Freq10 #= 5) or
	(Row10 #= 444 and Freq10 #= 5) or
	(Row10 #= 445 and Freq10 #= 4) or
	(Row10 #= 446 and Freq10 #= 7) or
	(Row10 #= 447 and Freq10 #= 9) or
	(Row10 #= 448 and Freq10 #= 10) or
	(Row10 #= 449 and Freq10 #= 10) or
	(Row10 #= 450 and Freq10 #= 9) or
	(Row10 #= 0 and Freq10 #= 0)), 

	((Row11 #= 4 and Freq11 #= 3) or
	(Row11 #= 5 and Freq11 #= 3) or
	(Row11 #= 6 and Freq11 #= 4) or
	(Row11 #= 7 and Freq11 #= 3) or
	(Row11 #= 8 and Freq11 #= 3) or
	(Row11 #= 9 and Freq11 #= 3) or
	(Row11 #= 10 and Freq11 #= 4) or
	(Row11 #= 11 and Freq11 #= 3) or
	(Row11 #= 12 and Freq11 #= 4) or
	(Row11 #= 13 and Freq11 #= 3) or
	(Row11 #= 14 and Freq11 #= 2) or
	(Row11 #= 15 and Freq11 #= 3) or
	(Row11 #= 16 and Freq11 #= 2) or
	(Row11 #= 17 and Freq11 #= 2) or
	(Row11 #= 18 and Freq11 #= 2) or
	(Row11 #= 19 and Freq11 #= 1) or
	(Row11 #= 20 and Freq11 #= 1) or
	(Row11 #= 21 and Freq11 #= 1) or
	(Row11 #= 22 and Freq11 #= 1) or
	(Row11 #= 23 and Freq11 #= 1) or
	(Row11 #= 24 and Freq11 #= 1) or
	(Row11 #= 25 and Freq11 #= 1) or
	(Row11 #= 26 and Freq11 #= 1) or
	(Row11 #= 420 and Freq11 #= 1) or
	(Row11 #= 425 and Freq11 #= 1) or
	(Row11 #= 426 and Freq11 #= 1) or
	(Row11 #= 431 and Freq11 #= 1) or
	(Row11 #= 432 and Freq11 #= 1) or
	(Row11 #= 433 and Freq11 #= 1) or
	(Row11 #= 434 and Freq11 #= 1) or
	(Row11 #= 435 and Freq11 #= 1) or
	(Row11 #= 436 and Freq11 #= 1) or
	(Row11 #= 437 and Freq11 #= 2) or
	(Row11 #= 438 and Freq11 #= 2) or
	(Row11 #= 439 and Freq11 #= 2) or
	(Row11 #= 440 and Freq11 #= 2) or
	(Row11 #= 441 and Freq11 #= 3) or
	(Row11 #= 442 and Freq11 #= 4) or
	(Row11 #= 443 and Freq11 #= 3) or
	(Row11 #= 444 and Freq11 #= 4) or
	(Row11 #= 445 and Freq11 #= 4) or
	(Row11 #= 446 and Freq11 #= 4) or
	(Row11 #= 447 and Freq11 #= 4) or
	(Row11 #= 448 and Freq11 #= 4) or
	(Row11 #= 449 and Freq11 #= 5) or
	(Row11 #= 450 and Freq11 #= 4) or
	(Row11 #= 0 and Freq11 #= 0)), 

	((Row12 #= 4 and Freq12 #= 2) or
	(Row12 #= 5 and Freq12 #= 3) or
	(Row12 #= 6 and Freq12 #= 4) or
	(Row12 #= 7 and Freq12 #= 3) or
	(Row12 #= 8 and Freq12 #= 3) or
	(Row12 #= 9 and Freq12 #= 2) or
	(Row12 #= 10 and Freq12 #= 2) or
	(Row12 #= 11 and Freq12 #= 3) or
	(Row12 #= 12 and Freq12 #= 3) or
	(Row12 #= 13 and Freq12 #= 2) or
	(Row12 #= 14 and Freq12 #= 2) or
	(Row12 #= 15 and Freq12 #= 2) or
	(Row12 #= 16 and Freq12 #= 2) or
	(Row12 #= 17 and Freq12 #= 2) or
	(Row12 #= 18 and Freq12 #= 1) or
	(Row12 #= 19 and Freq12 #= 1) or
	(Row12 #= 20 and Freq12 #= 1) or
	(Row12 #= 21 and Freq12 #= 1) or
	(Row12 #= 22 and Freq12 #= 1) or
	(Row12 #= 28 and Freq12 #= 1) or
	(Row12 #= 426 and Freq12 #= 1) or
	(Row12 #= 428 and Freq12 #= 1) or
	(Row12 #= 429 and Freq12 #= 1) or
	(Row12 #= 430 and Freq12 #= 1) or
	(Row12 #= 431 and Freq12 #= 1) or
	(Row12 #= 432 and Freq12 #= 1) or
	(Row12 #= 433 and Freq12 #= 1) or
	(Row12 #= 434 and Freq12 #= 1) or
	(Row12 #= 435 and Freq12 #= 1) or
	(Row12 #= 436 and Freq12 #= 1) or
	(Row12 #= 437 and Freq12 #= 1) or
	(Row12 #= 438 and Freq12 #= 1) or
	(Row12 #= 439 and Freq12 #= 2) or
	(Row12 #= 440 and Freq12 #= 2) or
	(Row12 #= 441 and Freq12 #= 2) or
	(Row12 #= 442 and Freq12 #= 2) or
	(Row12 #= 443 and Freq12 #= 3) or
	(Row12 #= 444 and Freq12 #= 3) or
	(Row12 #= 445 and Freq12 #= 3) or
	(Row12 #= 446 and Freq12 #= 4) or
	(Row12 #= 447 and Freq12 #= 4) or
	(Row12 #= 448 and Freq12 #= 4) or
	(Row12 #= 449 and Freq12 #= 4) or
	(Row12 #= 450 and Freq12 #= 4) or
	(Row12 #= 0 and Freq12 #= 0)), 

	((Row13 #= 4 and Freq13 #= 3) or
	(Row13 #= 5 and Freq13 #= 3) or
	(Row13 #= 6 and Freq13 #= 3) or
	(Row13 #= 7 and Freq13 #= 3) or
	(Row13 #= 8 and Freq13 #= 3) or
	(Row13 #= 9 and Freq13 #= 3) or
	(Row13 #= 10 and Freq13 #= 3) or
	(Row13 #= 11 and Freq13 #= 2) or
	(Row13 #= 12 and Freq13 #= 2) or
	(Row13 #= 13 and Freq13 #= 2) or
	(Row13 #= 14 and Freq13 #= 2) or
	(Row13 #= 15 and Freq13 #= 2) or
	(Row13 #= 16 and Freq13 #= 2) or
	(Row13 #= 17 and Freq13 #= 1) or
	(Row13 #= 18 and Freq13 #= 1) or
	(Row13 #= 19 and Freq13 #= 1) or
	(Row13 #= 20 and Freq13 #= 1) or
	(Row13 #= 21 and Freq13 #= 1) or
	(Row13 #= 22 and Freq13 #= 1) or
	(Row13 #= 23 and Freq13 #= 1) or
	(Row13 #= 428 and Freq13 #= 1) or
	(Row13 #= 430 and Freq13 #= 1) or
	(Row13 #= 431 and Freq13 #= 1) or
	(Row13 #= 432 and Freq13 #= 1) or
	(Row13 #= 433 and Freq13 #= 1) or
	(Row13 #= 434 and Freq13 #= 1) or
	(Row13 #= 435 and Freq13 #= 1) or
	(Row13 #= 436 and Freq13 #= 1) or
	(Row13 #= 437 and Freq13 #= 1) or
	(Row13 #= 438 and Freq13 #= 1) or
	(Row13 #= 439 and Freq13 #= 2) or
	(Row13 #= 440 and Freq13 #= 2) or
	(Row13 #= 441 and Freq13 #= 2) or
	(Row13 #= 442 and Freq13 #= 2) or
	(Row13 #= 443 and Freq13 #= 3) or
	(Row13 #= 444 and Freq13 #= 3) or
	(Row13 #= 445 and Freq13 #= 3) or
	(Row13 #= 446 and Freq13 #= 3) or
	(Row13 #= 447 and Freq13 #= 4) or
	(Row13 #= 448 and Freq13 #= 4) or
	(Row13 #= 449 and Freq13 #= 4) or
	(Row13 #= 450 and Freq13 #= 4) or
	(Row13 #= 0 and Freq13 #= 0)), 

	((Row14 #= 4 and Freq14 #= 2) or
	(Row14 #= 5 and Freq14 #= 3) or
	(Row14 #= 6 and Freq14 #= 3) or
	(Row14 #= 7 and Freq14 #= 3) or
	(Row14 #= 8 and Freq14 #= 3) or
	(Row14 #= 9 and Freq14 #= 2) or
	(Row14 #= 10 and Freq14 #= 2) or
	(Row14 #= 11 and Freq14 #= 3) or
	(Row14 #= 12 and Freq14 #= 3) or
	(Row14 #= 13 and Freq14 #= 2) or
	(Row14 #= 14 and Freq14 #= 2) or
	(Row14 #= 15 and Freq14 #= 2) or
	(Row14 #= 16 and Freq14 #= 2) or
	(Row14 #= 17 and Freq14 #= 1) or
	(Row14 #= 18 and Freq14 #= 1) or
	(Row14 #= 19 and Freq14 #= 1) or
	(Row14 #= 20 and Freq14 #= 1) or
	(Row14 #= 21 and Freq14 #= 1) or
	(Row14 #= 22 and Freq14 #= 1) or
	(Row14 #= 33 and Freq14 #= 1) or
	(Row14 #= 426 and Freq14 #= 1) or
	(Row14 #= 432 and Freq14 #= 1) or
	(Row14 #= 434 and Freq14 #= 1) or
	(Row14 #= 435 and Freq14 #= 1) or
	(Row14 #= 436 and Freq14 #= 1) or
	(Row14 #= 437 and Freq14 #= 1) or
	(Row14 #= 438 and Freq14 #= 1) or
	(Row14 #= 439 and Freq14 #= 2) or
	(Row14 #= 440 and Freq14 #= 2) or
	(Row14 #= 441 and Freq14 #= 2) or
	(Row14 #= 442 and Freq14 #= 2) or
	(Row14 #= 443 and Freq14 #= 2) or
	(Row14 #= 444 and Freq14 #= 2) or
	(Row14 #= 445 and Freq14 #= 2) or
	(Row14 #= 446 and Freq14 #= 3) or
	(Row14 #= 447 and Freq14 #= 4) or
	(Row14 #= 448 and Freq14 #= 4) or
	(Row14 #= 449 and Freq14 #= 4) or
	(Row14 #= 450 and Freq14 #= 4) or
	(Row14 #= 0 and Freq14 #= 0)), 

	((Row15 #= 4 and Freq15 #= 2) or
	(Row15 #= 5 and Freq15 #= 2) or
	(Row15 #= 6 and Freq15 #= 3) or
	(Row15 #= 7 and Freq15 #= 2) or
	(Row15 #= 8 and Freq15 #= 2) or
	(Row15 #= 9 and Freq15 #= 1) or
	(Row15 #= 10 and Freq15 #= 2) or
	(Row15 #= 11 and Freq15 #= 1) or
	(Row15 #= 12 and Freq15 #= 2) or
	(Row15 #= 13 and Freq15 #= 2) or
	(Row15 #= 14 and Freq15 #= 2) or
	(Row15 #= 15 and Freq15 #= 2) or
	(Row15 #= 16 and Freq15 #= 1) or
	(Row15 #= 17 and Freq15 #= 1) or
	(Row15 #= 18 and Freq15 #= 1) or
	(Row15 #= 19 and Freq15 #= 1) or
	(Row15 #= 20 and Freq15 #= 1) or
	(Row15 #= 21 and Freq15 #= 1) or
	(Row15 #= 22 and Freq15 #= 1) or
	(Row15 #= 23 and Freq15 #= 1) or
	(Row15 #= 25 and Freq15 #= 1) or
	(Row15 #= 26 and Freq15 #= 1) or
	(Row15 #= 431 and Freq15 #= 1) or
	(Row15 #= 432 and Freq15 #= 1) or
	(Row15 #= 433 and Freq15 #= 1) or
	(Row15 #= 434 and Freq15 #= 1) or
	(Row15 #= 435 and Freq15 #= 1) or
	(Row15 #= 436 and Freq15 #= 1) or
	(Row15 #= 437 and Freq15 #= 1) or
	(Row15 #= 438 and Freq15 #= 1) or
	(Row15 #= 439 and Freq15 #= 1) or
	(Row15 #= 440 and Freq15 #= 2) or
	(Row15 #= 441 and Freq15 #= 2) or
	(Row15 #= 442 and Freq15 #= 3) or
	(Row15 #= 443 and Freq15 #= 2) or
	(Row15 #= 444 and Freq15 #= 2) or
	(Row15 #= 445 and Freq15 #= 2) or
	(Row15 #= 446 and Freq15 #= 2) or
	(Row15 #= 447 and Freq15 #= 3) or
	(Row15 #= 448 and Freq15 #= 3) or
	(Row15 #= 449 and Freq15 #= 3) or
	(Row15 #= 450 and Freq15 #= 3) or
	(Row15 #= 0 and Freq15 #= 0)), 

	((Row16 #= 4 and Freq16 #= 1) or
	(Row16 #= 5 and Freq16 #= 1) or
	(Row16 #= 6 and Freq16 #= 2) or
	(Row16 #= 7 and Freq16 #= 2) or
	(Row16 #= 8 and Freq16 #= 2) or
	(Row16 #= 9 and Freq16 #= 2) or
	(Row16 #= 10 and Freq16 #= 2) or
	(Row16 #= 11 and Freq16 #= 2) or
	(Row16 #= 12 and Freq16 #= 1) or
	(Row16 #= 13 and Freq16 #= 2) or
	(Row16 #= 14 and Freq16 #= 1) or
	(Row16 #= 15 and Freq16 #= 1) or
	(Row16 #= 16 and Freq16 #= 1) or
	(Row16 #= 17 and Freq16 #= 1) or
	(Row16 #= 18 and Freq16 #= 1) or
	(Row16 #= 19 and Freq16 #= 1) or
	(Row16 #= 20 and Freq16 #= 1) or
	(Row16 #= 21 and Freq16 #= 1) or
	(Row16 #= 429 and Freq16 #= 1) or
	(Row16 #= 430 and Freq16 #= 1) or
	(Row16 #= 433 and Freq16 #= 1) or
	(Row16 #= 434 and Freq16 #= 1) or
	(Row16 #= 435 and Freq16 #= 1) or
	(Row16 #= 436 and Freq16 #= 1) or
	(Row16 #= 437 and Freq16 #= 1) or
	(Row16 #= 438 and Freq16 #= 1) or
	(Row16 #= 439 and Freq16 #= 1) or
	(Row16 #= 440 and Freq16 #= 1) or
	(Row16 #= 441 and Freq16 #= 1) or
	(Row16 #= 442 and Freq16 #= 1) or
	(Row16 #= 443 and Freq16 #= 2) or
	(Row16 #= 444 and Freq16 #= 2) or
	(Row16 #= 445 and Freq16 #= 2) or
	(Row16 #= 446 and Freq16 #= 2) or
	(Row16 #= 447 and Freq16 #= 2) or
	(Row16 #= 448 and Freq16 #= 2) or
	(Row16 #= 449 and Freq16 #= 3) or
	(Row16 #= 450 and Freq16 #= 2) or
	(Row16 #= 0 and Freq16 #= 0)), 

	((Row17 #= 4 and Freq17 #= 1) or
	(Row17 #= 5 and Freq17 #= 1) or
	(Row17 #= 6 and Freq17 #= 1) or
	(Row17 #= 7 and Freq17 #= 1) or
	(Row17 #= 8 and Freq17 #= 1) or
	(Row17 #= 9 and Freq17 #= 1) or
	(Row17 #= 10 and Freq17 #= 1) or
	(Row17 #= 11 and Freq17 #= 1) or
	(Row17 #= 12 and Freq17 #= 1) or
	(Row17 #= 13 and Freq17 #= 1) or
	(Row17 #= 14 and Freq17 #= 1) or
	(Row17 #= 15 and Freq17 #= 1) or
	(Row17 #= 16 and Freq17 #= 1) or
	(Row17 #= 17 and Freq17 #= 1) or
	(Row17 #= 18 and Freq17 #= 1) or
	(Row17 #= 19 and Freq17 #= 1) or
	(Row17 #= 20 and Freq17 #= 1) or
	(Row17 #= 21 and Freq17 #= 1) or
	(Row17 #= 430 and Freq17 #= 1) or
	(Row17 #= 432 and Freq17 #= 1) or
	(Row17 #= 433 and Freq17 #= 1) or
	(Row17 #= 434 and Freq17 #= 1) or
	(Row17 #= 435 and Freq17 #= 1) or
	(Row17 #= 436 and Freq17 #= 1) or
	(Row17 #= 437 and Freq17 #= 1) or
	(Row17 #= 438 and Freq17 #= 1) or
	(Row17 #= 439 and Freq17 #= 1) or
	(Row17 #= 440 and Freq17 #= 1) or
	(Row17 #= 441 and Freq17 #= 1) or
	(Row17 #= 442 and Freq17 #= 1) or
	(Row17 #= 443 and Freq17 #= 1) or
	(Row17 #= 444 and Freq17 #= 1) or
	(Row17 #= 445 and Freq17 #= 2) or
	(Row17 #= 446 and Freq17 #= 1) or
	(Row17 #= 447 and Freq17 #= 2) or
	(Row17 #= 448 and Freq17 #= 2) or
	(Row17 #= 449 and Freq17 #= 2) or
	(Row17 #= 450 and Freq17 #= 2) or
	(Row17 #= 0 and Freq17 #= 0)), 

	((Row18 #= 4 and Freq18 #= 1) or
	(Row18 #= 5 and Freq18 #= 1) or
	(Row18 #= 6 and Freq18 #= 1) or
	(Row18 #= 7 and Freq18 #= 1) or
	(Row18 #= 8 and Freq18 #= 1) or
	(Row18 #= 9 and Freq18 #= 1) or
	(Row18 #= 10 and Freq18 #= 2) or
	(Row18 #= 11 and Freq18 #= 1) or
	(Row18 #= 12 and Freq18 #= 1) or
	(Row18 #= 13 and Freq18 #= 1) or
	(Row18 #= 14 and Freq18 #= 1) or
	(Row18 #= 15 and Freq18 #= 1) or
	(Row18 #= 16 and Freq18 #= 1) or
	(Row18 #= 17 and Freq18 #= 1) or
	(Row18 #= 18 and Freq18 #= 1) or
	(Row18 #= 19 and Freq18 #= 1) or
	(Row18 #= 20 and Freq18 #= 1) or
	(Row18 #= 21 and Freq18 #= 1) or
	(Row18 #= 425 and Freq18 #= 1) or
	(Row18 #= 433 and Freq18 #= 1) or
	(Row18 #= 434 and Freq18 #= 1) or
	(Row18 #= 435 and Freq18 #= 1) or
	(Row18 #= 436 and Freq18 #= 1) or
	(Row18 #= 437 and Freq18 #= 1) or
	(Row18 #= 438 and Freq18 #= 1) or
	(Row18 #= 439 and Freq18 #= 1) or
	(Row18 #= 440 and Freq18 #= 1) or
	(Row18 #= 441 and Freq18 #= 1) or
	(Row18 #= 442 and Freq18 #= 1) or
	(Row18 #= 443 and Freq18 #= 1) or
	(Row18 #= 444 and Freq18 #= 1) or
	(Row18 #= 445 and Freq18 #= 1) or
	(Row18 #= 446 and Freq18 #= 1) or
	(Row18 #= 447 and Freq18 #= 1) or
	(Row18 #= 448 and Freq18 #= 1) or
	(Row18 #= 449 and Freq18 #= 1) or
	(Row18 #= 450 and Freq18 #= 1) or
	(Row18 #= 0 and Freq18 #= 0)), 

	((Row19 #= 5 and Freq19 #= 1) or
	(Row19 #= 6 and Freq19 #= 1) or
	(Row19 #= 7 and Freq19 #= 1) or
	(Row19 #= 8 and Freq19 #= 1) or
	(Row19 #= 10 and Freq19 #= 1) or
	(Row19 #= 11 and Freq19 #= 1) or
	(Row19 #= 12 and Freq19 #= 1) or
	(Row19 #= 13 and Freq19 #= 1) or
	(Row19 #= 14 and Freq19 #= 1) or
	(Row19 #= 15 and Freq19 #= 1) or
	(Row19 #= 16 and Freq19 #= 1) or
	(Row19 #= 18 and Freq19 #= 1) or
	(Row19 #= 19 and Freq19 #= 1) or
	(Row19 #= 20 and Freq19 #= 1) or
	(Row19 #= 425 and Freq19 #= 1) or
	(Row19 #= 428 and Freq19 #= 1) or
	(Row19 #= 429 and Freq19 #= 1) or
	(Row19 #= 430 and Freq19 #= 1) or
	(Row19 #= 432 and Freq19 #= 1) or
	(Row19 #= 433 and Freq19 #= 1) or
	(Row19 #= 434 and Freq19 #= 1) or
	(Row19 #= 435 and Freq19 #= 1) or
	(Row19 #= 436 and Freq19 #= 1) or
	(Row19 #= 437 and Freq19 #= 1) or
	(Row19 #= 438 and Freq19 #= 1) or
	(Row19 #= 439 and Freq19 #= 1) or
	(Row19 #= 440 and Freq19 #= 1) or
	(Row19 #= 441 and Freq19 #= 1) or
	(Row19 #= 442 and Freq19 #= 1) or
	(Row19 #= 443 and Freq19 #= 1) or
	(Row19 #= 444 and Freq19 #= 1) or
	(Row19 #= 445 and Freq19 #= 1) or
	(Row19 #= 446 and Freq19 #= 1) or
	(Row19 #= 447 and Freq19 #= 1) or
	(Row19 #= 448 and Freq19 #= 1) or
	(Row19 #= 449 and Freq19 #= 1) or
	(Row19 #= 450 and Freq19 #= 1) or
	(Row19 #= 0 and Freq19 #= 0)), 

	((Row20 #= 5 and Freq20 #= 1) or
	(Row20 #= 6 and Freq20 #= 1) or
	(Row20 #= 7 and Freq20 #= 1) or
	(Row20 #= 9 and Freq20 #= 1) or
	(Row20 #= 10 and Freq20 #= 1) or
	(Row20 #= 11 and Freq20 #= 1) or
	(Row20 #= 12 and Freq20 #= 1) or
	(Row20 #= 13 and Freq20 #= 1) or
	(Row20 #= 14 and Freq20 #= 1) or
	(Row20 #= 16 and Freq20 #= 1) or
	(Row20 #= 17 and Freq20 #= 1) or
	(Row20 #= 19 and Freq20 #= 1) or
	(Row20 #= 420 and Freq20 #= 1) or
	(Row20 #= 424 and Freq20 #= 1) or
	(Row20 #= 425 and Freq20 #= 1) or
	(Row20 #= 426 and Freq20 #= 1) or
	(Row20 #= 429 and Freq20 #= 1) or
	(Row20 #= 430 and Freq20 #= 1) or
	(Row20 #= 431 and Freq20 #= 1) or
	(Row20 #= 432 and Freq20 #= 1) or
	(Row20 #= 433 and Freq20 #= 1) or
	(Row20 #= 434 and Freq20 #= 1) or
	(Row20 #= 435 and Freq20 #= 1) or
	(Row20 #= 436 and Freq20 #= 1) or
	(Row20 #= 437 and Freq20 #= 1) or
	(Row20 #= 438 and Freq20 #= 1) or
	(Row20 #= 439 and Freq20 #= 1) or
	(Row20 #= 440 and Freq20 #= 1) or
	(Row20 #= 441 and Freq20 #= 1) or
	(Row20 #= 442 and Freq20 #= 1) or
	(Row20 #= 443 and Freq20 #= 1) or
	(Row20 #= 444 and Freq20 #= 1) or
	(Row20 #= 445 and Freq20 #= 2) or
	(Row20 #= 446 and Freq20 #= 1) or
	(Row20 #= 447 and Freq20 #= 1) or
	(Row20 #= 448 and Freq20 #= 1) or
	(Row20 #= 449 and Freq20 #= 1) or
	(Row20 #= 450 and Freq20 #= 1) or
	(Row20 #= 0 and Freq20 #= 0)), 

	((Row21 #= 6 and Freq21 #= 1) or
	(Row21 #= 12 and Freq21 #= 1) or
	(Row21 #= 13 and Freq21 #= 1) or
	(Row21 #= 14 and Freq21 #= 1) or
	(Row21 #= 15 and Freq21 #= 1) or
	(Row21 #= 17 and Freq21 #= 1) or
	(Row21 #= 19 and Freq21 #= 1) or
	(Row21 #= 425 and Freq21 #= 1) or
	(Row21 #= 434 and Freq21 #= 1) or
	(Row21 #= 435 and Freq21 #= 1) or
	(Row21 #= 436 and Freq21 #= 1) or
	(Row21 #= 437 and Freq21 #= 1) or
	(Row21 #= 438 and Freq21 #= 1) or
	(Row21 #= 439 and Freq21 #= 1) or
	(Row21 #= 440 and Freq21 #= 1) or
	(Row21 #= 441 and Freq21 #= 1) or
	(Row21 #= 442 and Freq21 #= 1) or
	(Row21 #= 443 and Freq21 #= 1) or
	(Row21 #= 444 and Freq21 #= 1) or
	(Row21 #= 445 and Freq21 #= 1) or
	(Row21 #= 446 and Freq21 #= 1) or
	(Row21 #= 447 and Freq21 #= 1) or
	(Row21 #= 448 and Freq21 #= 1) or
	(Row21 #= 449 and Freq21 #= 1) or
	(Row21 #= 450 and Freq21 #= 1) or
	(Row21 #= 0 and Freq21 #= 0)), 

	((Row22 #= 11 and Freq22 #= 1) or
	(Row22 #= 14 and Freq22 #= 1) or
	(Row22 #= 18 and Freq22 #= 1) or
	(Row22 #= 427 and Freq22 #= 1) or
	(Row22 #= 429 and Freq22 #= 1) or
	(Row22 #= 435 and Freq22 #= 1) or
	(Row22 #= 436 and Freq22 #= 1) or
	(Row22 #= 437 and Freq22 #= 1) or
	(Row22 #= 438 and Freq22 #= 1) or
	(Row22 #= 439 and Freq22 #= 1) or
	(Row22 #= 440 and Freq22 #= 1) or
	(Row22 #= 441 and Freq22 #= 1) or
	(Row22 #= 442 and Freq22 #= 1) or
	(Row22 #= 443 and Freq22 #= 1) or
	(Row22 #= 444 and Freq22 #= 1) or
	(Row22 #= 445 and Freq22 #= 1) or
	(Row22 #= 446 and Freq22 #= 1) or
	(Row22 #= 449 and Freq22 #= 1) or
	(Row22 #= 0 and Freq22 #= 0)), 

	((Row23 #= 6 and Freq23 #= 1) or
	(Row23 #= 10 and Freq23 #= 1) or
	(Row23 #= 13 and Freq23 #= 1) or
	(Row23 #= 17 and Freq23 #= 1) or
	(Row23 #= 434 and Freq23 #= 1) or
	(Row23 #= 436 and Freq23 #= 1) or
	(Row23 #= 437 and Freq23 #= 1) or
	(Row23 #= 438 and Freq23 #= 1) or
	(Row23 #= 439 and Freq23 #= 1) or
	(Row23 #= 440 and Freq23 #= 1) or
	(Row23 #= 441 and Freq23 #= 1) or
	(Row23 #= 442 and Freq23 #= 1) or
	(Row23 #= 443 and Freq23 #= 1) or
	(Row23 #= 444 and Freq23 #= 1) or
	(Row23 #= 445 and Freq23 #= 1) or
	(Row23 #= 447 and Freq23 #= 1) or
	(Row23 #= 449 and Freq23 #= 1) or
	(Row23 #= 0 and Freq23 #= 0)), 

	((Row24 #= 6 and Freq24 #= 1) or
	(Row24 #= 10 and Freq24 #= 1) or
	(Row24 #= 11 and Freq24 #= 1) or
	(Row24 #= 13 and Freq24 #= 1) or
	(Row24 #= 14 and Freq24 #= 1) or
	(Row24 #= 17 and Freq24 #= 1) or
	(Row24 #= 19 and Freq24 #= 1) or
	(Row24 #= 403 and Freq24 #= 1) or
	(Row24 #= 425 and Freq24 #= 1) or
	(Row24 #= 430 and Freq24 #= 1) or
	(Row24 #= 435 and Freq24 #= 1) or
	(Row24 #= 436 and Freq24 #= 1) or
	(Row24 #= 439 and Freq24 #= 1) or
	(Row24 #= 440 and Freq24 #= 1) or
	(Row24 #= 441 and Freq24 #= 1) or
	(Row24 #= 442 and Freq24 #= 1) or
	(Row24 #= 443 and Freq24 #= 1) or
	(Row24 #= 444 and Freq24 #= 1) or
	(Row24 #= 445 and Freq24 #= 1) or
	(Row24 #= 447 and Freq24 #= 1) or
	(Row24 #= 448 and Freq24 #= 1) or
	(Row24 #= 449 and Freq24 #= 1) or
	(Row24 #= 0 and Freq24 #= 0)), 

	((Row25 #= 13 and Freq25 #= 1) or
	(Row25 #= 436 and Freq25 #= 1) or
	(Row25 #= 439 and Freq25 #= 1) or
	(Row25 #= 440 and Freq25 #= 1) or
	(Row25 #= 441 and Freq25 #= 1) or
	(Row25 #= 444 and Freq25 #= 1) or
	(Row25 #= 445 and Freq25 #= 1) or
	(Row25 #= 446 and Freq25 #= 1) or
	(Row25 #= 450 and Freq25 #= 1) or
	(Row25 #= 0 and Freq25 #= 0)), 

	((Row26 #= 436 and Freq26 #= 1) or
	(Row26 #= 438 and Freq26 #= 1) or
	(Row26 #= 440 and Freq26 #= 1) or
	(Row26 #= 441 and Freq26 #= 1) or
	(Row26 #= 442 and Freq26 #= 1) or
	(Row26 #= 444 and Freq26 #= 1) or
	(Row26 #= 446 and Freq26 #= 1) or
	(Row26 #= 0 and Freq26 #= 0)), 

	((Row27 #= 19 and Freq27 #= 1) or
	(Row27 #= 441 and Freq27 #= 1) or
	(Row27 #= 0 and Freq27 #= 0)), 

	((Row28 #= 435 and Freq28 #= 1) or
	(Row28 #= 0 and Freq28 #= 0)), 

	((Row29 #= 440 and Freq29 #= 1) or
	(Row29 #= 442 and Freq29 #= 1) or
	(Row29 #= 0 and Freq29 #= 0)), 

	((Row30 #= 0 and Freq30 #= 0)), 

	((Row31 #= 0 and Freq31 #= 0)), 

	((Row32 #= 442 and Freq32 #= 1) or
	(Row32 #= 0 and Freq32 #= 0)), 

	((Row33 #= 0 and Freq33 #= 0)), 

	((Row34 #= 0 and Freq34 #= 0)), 

	((Row35 #= 441 and Freq35 #= 1) or
	(Row35 #= 0 and Freq35 #= 0)), 

	((Row36 #= 0 and Freq36 #= 0)), 

	((Row37 #= 0 and Freq37 #= 0)), 

	((Row38 #= 0 and Freq38 #= 0)), 

	((Row39 #= 0 and Freq39 #= 0)), 

	((Row40 #= 0 and Freq40 #= 0)), 

	((Row41 #= 0 and Freq41 #= 0)), 

	((Row42 #= 0 and Freq42 #= 0)), 

	((Row43 #= 0 and Freq43 #= 0)), 

	((Row44 #= 0 and Freq44 #= 0)), 

	((Row45 #= 0 and Freq45 #= 0)), 

	((Row46 #= 0 and Freq46 #= 0)), 

	((Row47 #= 0 and Freq47 #= 0)), 

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

	((Row58 #= 0 and Freq58 #= 0)), 

	((Row59 #= 426 and Freq59 #= 1) or
	(Row59 #= 0 and Freq59 #= 0)), 

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

	((Row94 #= 111 and Freq94 #= 1) or
	(Row94 #= 0 and Freq94 #= 0)), 

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

	((Row107 #= 0 and Freq107 #= 0)), 

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

	((Row119 #= 0 and Freq119 #= 0)), 

	((Row120 #= 0 and Freq120 #= 0)), 

	((Row121 #= 0 and Freq121 #= 0)), 

	((Row122 #= 0 and Freq122 #= 0)), 

	((Row123 #= 0 and Freq123 #= 0)), 

	((Row124 #= 0 and Freq124 #= 0)), 

	((Row125 #= 0 and Freq125 #= 0)), 

	((Row126 #= 0 and Freq126 #= 0)), 

	((Row127 #= 0 and Freq127 #= 0)), 

	((Row128 #= 0 and Freq128 #= 0)), 

	((Row129 #= 0 and Freq129 #= 0)), 

	((Row130 #= 0 and Freq130 #= 0)), 

	((Row131 #= 0 and Freq131 #= 0)), 

	((Row132 #= 0 and Freq132 #= 0)), 

	((Row133 #= 0 and Freq133 #= 0)), 

	((Row134 #= 0 and Freq134 #= 0)), 

	((Row135 #= 0 and Freq135 #= 0)), 

	((Row136 #= 0 and Freq136 #= 0)), 

	((Row137 #= 0 and Freq137 #= 0)), 

	((Row138 #= 0 and Freq138 #= 0)), 

	((Row139 #= 0 and Freq139 #= 0)), 

	((Row140 #= 0 and Freq140 #= 0)), 

	((Row141 #= 0 and Freq141 #= 0)), 

	((Row142 #= 0 and Freq142 #= 0)), 

	((Row143 #= 0 and Freq143 #= 0)), 

	((Row144 #= 0 and Freq144 #= 0)), 

	((Row145 #= 0 and Freq145 #= 0)), 

	((Row146 #= 0 and Freq146 #= 0)), 

	((Row147 #= 0 and Freq147 #= 0)), 

	((Row148 #= 0 and Freq148 #= 0)), 

	((Row149 #= 0 and Freq149 #= 0)), 

	((Row150 #= 0 and Freq150 #= 0)), 

	((Row151 #= 0 and Freq151 #= 0)), 

	((Row152 #= 0 and Freq152 #= 0)), 

	((Row153 #= 0 and Freq153 #= 0)), 

	((Row154 #= 0 and Freq154 #= 0)), 

	((Row155 #= 0 and Freq155 #= 0)), 

	((Row156 #= 0 and Freq156 #= 0)), 

	((Row157 #= 0 and Freq157 #= 0)), 

	((Row158 #= 0 and Freq158 #= 0)), 

	((Row159 #= 0 and Freq159 #= 0)), 

	((Row160 #= 0 and Freq160 #= 0)), 

	((Row161 #= 0 and Freq161 #= 0)), 

	((Row162 #= 0 and Freq162 #= 0)), 

	((Row163 #= 0 and Freq163 #= 0)), 

	((Row164 #= 0 and Freq164 #= 0)), 

	((Row165 #= 0 and Freq165 #= 0)), 

	((Row166 #= 0 and Freq166 #= 0)), 

	((Row167 #= 0 and Freq167 #= 0)), 

	((Row168 #= 0 and Freq168 #= 0)), 

	((Row169 #= 0 and Freq169 #= 0)), 

	((Row170 #= 0 and Freq170 #= 0)), 

	((Row171 #= 0 and Freq171 #= 0)), 

	((Row172 #= 0 and Freq172 #= 0)), 

	((Row173 #= 0 and Freq173 #= 0)), 

	((Row174 #= 0 and Freq174 #= 0)), 

	((Row175 #= 0 and Freq175 #= 0)), 

	((Row176 #= 0 and Freq176 #= 0)), 

	((Row177 #= 0 and Freq177 #= 0)), 

	((Row178 #= 0 and Freq178 #= 0)), 

	((Row179 #= 0 and Freq179 #= 0)), 

	((Row180 #= 0 and Freq180 #= 0)), 

	((Row181 #= 0 and Freq181 #= 0)), 

	((Row182 #= 0 and Freq182 #= 0)), 

	((Row183 #= 0 and Freq183 #= 0)), 

	((Row184 #= 0 and Freq184 #= 0)), 

	((Row185 #= 0 and Freq185 #= 0)), 

	((Row186 #= 244 and Freq186 #= 1) or
	(Row186 #= 0 and Freq186 #= 0)), 

	((Row187 #= 0 and Freq187 #= 0)), 

	((Row188 #= 0 and Freq188 #= 0)), 

	((Row189 #= 0 and Freq189 #= 0)), 

	((Row190 #= 0 and Freq190 #= 0)), 

	((Row191 #= 0 and Freq191 #= 0)), 

	((Row192 #= 0 and Freq192 #= 0)), 

	((Row193 #= 0 and Freq193 #= 0)), 

	((Row194 #= 0 and Freq194 #= 0)), 

	((Row195 #= 0 and Freq195 #= 0)), 

	((Row196 #= 0 and Freq196 #= 0)), 

	((Row197 #= 0 and Freq197 #= 0)), 

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

	((Row265 #= 0 and Freq265 #= 0)), 

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

	((Row319 #= 0 and Freq319 #= 0)), 

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

	((Row340 #= 0 and Freq340 #= 0)), 

	((Row341 #= 0 and Freq341 #= 0)), 

	((Row342 #= 0 and Freq342 #= 0)), 

	((Row343 #= 0 and Freq343 #= 0)), 

	((Row344 #= 0 and Freq344 #= 0)), 

	((Row345 #= 0 and Freq345 #= 0)), 

	((Row346 #= 0 and Freq346 #= 0)), 

	((Row347 #= 0 and Freq347 #= 0)), 

	((Row348 #= 176 and Freq348 #= 1) or
	(Row348 #= 0 and Freq348 #= 0)), 

	((Row349 #= 146 and Freq349 #= 1) or
	(Row349 #= 168 and Freq349 #= 1) or
	(Row349 #= 182 and Freq349 #= 1) or
	(Row349 #= 0 and Freq349 #= 0)), 

	((Row350 #= 144 and Freq350 #= 1) or
	(Row350 #= 178 and Freq350 #= 1) or
	(Row350 #= 0 and Freq350 #= 0)), 

	((Row351 #= 0 and Freq351 #= 0)), 

	((Row352 #= 158 and Freq352 #= 1) or
	(Row352 #= 173 and Freq352 #= 1) or
	(Row352 #= 0 and Freq352 #= 0)), 

	((Row353 #= 154 and Freq353 #= 1) or
	(Row353 #= 168 and Freq353 #= 1) or
	(Row353 #= 169 and Freq353 #= 1) or
	(Row353 #= 176 and Freq353 #= 1) or
	(Row353 #= 0 and Freq353 #= 0)), 

	((Row354 #= 167 and Freq354 #= 1) or
	(Row354 #= 169 and Freq354 #= 1) or
	(Row354 #= 0 and Freq354 #= 0)), 

	((Row355 #= 153 and Freq355 #= 1) or
	(Row355 #= 155 and Freq355 #= 1) or
	(Row355 #= 166 and Freq355 #= 1) or
	(Row355 #= 176 and Freq355 #= 1) or
	(Row355 #= 0 and Freq355 #= 0)), 

	((Row356 #= 152 and Freq356 #= 1) or
	(Row356 #= 154 and Freq356 #= 1) or
	(Row356 #= 155 and Freq356 #= 1) or
	(Row356 #= 158 and Freq356 #= 1) or
	(Row356 #= 166 and Freq356 #= 1) or
	(Row356 #= 167 and Freq356 #= 1) or
	(Row356 #= 168 and Freq356 #= 1) or
	(Row356 #= 169 and Freq356 #= 1) or
	(Row356 #= 172 and Freq356 #= 1) or
	(Row356 #= 176 and Freq356 #= 1) or
	(Row356 #= 177 and Freq356 #= 1) or
	(Row356 #= 0 and Freq356 #= 0)), 

	((Row357 #= 158 and Freq357 #= 1) or
	(Row357 #= 159 and Freq357 #= 1) or
	(Row357 #= 166 and Freq357 #= 1) or
	(Row357 #= 167 and Freq357 #= 1) or
	(Row357 #= 168 and Freq357 #= 1) or
	(Row357 #= 169 and Freq357 #= 1) or
	(Row357 #= 177 and Freq357 #= 1) or
	(Row357 #= 0 and Freq357 #= 0)), 

	((Row358 #= 146 and Freq358 #= 1) or
	(Row358 #= 147 and Freq358 #= 1) or
	(Row358 #= 151 and Freq358 #= 1) or
	(Row358 #= 153 and Freq358 #= 1) or
	(Row358 #= 154 and Freq358 #= 1) or
	(Row358 #= 155 and Freq358 #= 1) or
	(Row358 #= 156 and Freq358 #= 1) or
	(Row358 #= 157 and Freq358 #= 1) or
	(Row358 #= 158 and Freq358 #= 1) or
	(Row358 #= 159 and Freq358 #= 1) or
	(Row358 #= 166 and Freq358 #= 1) or
	(Row358 #= 167 and Freq358 #= 1) or
	(Row358 #= 168 and Freq358 #= 1) or
	(Row358 #= 169 and Freq358 #= 1) or
	(Row358 #= 170 and Freq358 #= 1) or
	(Row358 #= 171 and Freq358 #= 1) or
	(Row358 #= 172 and Freq358 #= 1) or
	(Row358 #= 173 and Freq358 #= 1) or
	(Row358 #= 175 and Freq358 #= 1) or
	(Row358 #= 176 and Freq358 #= 1) or
	(Row358 #= 177 and Freq358 #= 1) or
	(Row358 #= 178 and Freq358 #= 1) or
	(Row358 #= 179 and Freq358 #= 1) or
	(Row358 #= 180 and Freq358 #= 1) or
	(Row358 #= 181 and Freq358 #= 1) or
	(Row358 #= 185 and Freq358 #= 1) or
	(Row358 #= 0 and Freq358 #= 0)), 

	((Row359 #= 151 and Freq359 #= 1) or
	(Row359 #= 152 and Freq359 #= 1) or
	(Row359 #= 153 and Freq359 #= 1) or
	(Row359 #= 154 and Freq359 #= 1) or
	(Row359 #= 155 and Freq359 #= 1) or
	(Row359 #= 157 and Freq359 #= 1) or
	(Row359 #= 158 and Freq359 #= 1) or
	(Row359 #= 159 and Freq359 #= 1) or
	(Row359 #= 166 and Freq359 #= 1) or
	(Row359 #= 167 and Freq359 #= 1) or
	(Row359 #= 168 and Freq359 #= 1) or
	(Row359 #= 169 and Freq359 #= 1) or
	(Row359 #= 170 and Freq359 #= 1) or
	(Row359 #= 171 and Freq359 #= 1) or
	(Row359 #= 172 and Freq359 #= 1) or
	(Row359 #= 173 and Freq359 #= 1) or
	(Row359 #= 175 and Freq359 #= 1) or
	(Row359 #= 176 and Freq359 #= 1) or
	(Row359 #= 177 and Freq359 #= 1) or
	(Row359 #= 0 and Freq359 #= 0)), 

	((Row360 #= 151 and Freq360 #= 1) or
	(Row360 #= 152 and Freq360 #= 1) or
	(Row360 #= 153 and Freq360 #= 1) or
	(Row360 #= 155 and Freq360 #= 1) or
	(Row360 #= 157 and Freq360 #= 1) or
	(Row360 #= 158 and Freq360 #= 1) or
	(Row360 #= 159 and Freq360 #= 1) or
	(Row360 #= 166 and Freq360 #= 1) or
	(Row360 #= 167 and Freq360 #= 1) or
	(Row360 #= 168 and Freq360 #= 1) or
	(Row360 #= 169 and Freq360 #= 1) or
	(Row360 #= 170 and Freq360 #= 1) or
	(Row360 #= 171 and Freq360 #= 1) or
	(Row360 #= 172 and Freq360 #= 1) or
	(Row360 #= 173 and Freq360 #= 1) or
	(Row360 #= 174 and Freq360 #= 1) or
	(Row360 #= 175 and Freq360 #= 1) or
	(Row360 #= 176 and Freq360 #= 1) or
	(Row360 #= 178 and Freq360 #= 1) or
	(Row360 #= 179 and Freq360 #= 1) or
	(Row360 #= 182 and Freq360 #= 1) or
	(Row360 #= 0 and Freq360 #= 0)), 

	((Row361 #= 149 and Freq361 #= 1) or
	(Row361 #= 151 and Freq361 #= 1) or
	(Row361 #= 152 and Freq361 #= 1) or
	(Row361 #= 154 and Freq361 #= 1) or
	(Row361 #= 155 and Freq361 #= 1) or
	(Row361 #= 156 and Freq361 #= 1) or
	(Row361 #= 157 and Freq361 #= 1) or
	(Row361 #= 158 and Freq361 #= 1) or
	(Row361 #= 159 and Freq361 #= 1) or
	(Row361 #= 166 and Freq361 #= 1) or
	(Row361 #= 167 and Freq361 #= 1) or
	(Row361 #= 168 and Freq361 #= 1) or
	(Row361 #= 169 and Freq361 #= 1) or
	(Row361 #= 170 and Freq361 #= 1) or
	(Row361 #= 171 and Freq361 #= 1) or
	(Row361 #= 173 and Freq361 #= 1) or
	(Row361 #= 174 and Freq361 #= 1) or
	(Row361 #= 175 and Freq361 #= 1) or
	(Row361 #= 176 and Freq361 #= 1) or
	(Row361 #= 177 and Freq361 #= 1) or
	(Row361 #= 0 and Freq361 #= 0)), 

	((Row362 #= 149 and Freq362 #= 1) or
	(Row362 #= 151 and Freq362 #= 1) or
	(Row362 #= 152 and Freq362 #= 1) or
	(Row362 #= 153 and Freq362 #= 1) or
	(Row362 #= 154 and Freq362 #= 1) or
	(Row362 #= 155 and Freq362 #= 1) or
	(Row362 #= 156 and Freq362 #= 1) or
	(Row362 #= 157 and Freq362 #= 1) or
	(Row362 #= 158 and Freq362 #= 1) or
	(Row362 #= 159 and Freq362 #= 1) or
	(Row362 #= 166 and Freq362 #= 1) or
	(Row362 #= 167 and Freq362 #= 1) or
	(Row362 #= 168 and Freq362 #= 1) or
	(Row362 #= 169 and Freq362 #= 1) or
	(Row362 #= 170 and Freq362 #= 1) or
	(Row362 #= 171 and Freq362 #= 1) or
	(Row362 #= 172 and Freq362 #= 1) or
	(Row362 #= 173 and Freq362 #= 1) or
	(Row362 #= 174 and Freq362 #= 1) or
	(Row362 #= 175 and Freq362 #= 1) or
	(Row362 #= 176 and Freq362 #= 1) or
	(Row362 #= 177 and Freq362 #= 1) or
	(Row362 #= 178 and Freq362 #= 1) or
	(Row362 #= 179 and Freq362 #= 1) or
	(Row362 #= 182 and Freq362 #= 1) or
	(Row362 #= 0 and Freq362 #= 0)), 

	((Row363 #= 150 and Freq363 #= 1) or
	(Row363 #= 151 and Freq363 #= 1) or
	(Row363 #= 154 and Freq363 #= 1) or
	(Row363 #= 155 and Freq363 #= 1) or
	(Row363 #= 156 and Freq363 #= 1) or
	(Row363 #= 157 and Freq363 #= 1) or
	(Row363 #= 158 and Freq363 #= 1) or
	(Row363 #= 159 and Freq363 #= 1) or
	(Row363 #= 166 and Freq363 #= 1) or
	(Row363 #= 167 and Freq363 #= 1) or
	(Row363 #= 168 and Freq363 #= 1) or
	(Row363 #= 169 and Freq363 #= 1) or
	(Row363 #= 170 and Freq363 #= 1) or
	(Row363 #= 171 and Freq363 #= 1) or
	(Row363 #= 172 and Freq363 #= 1) or
	(Row363 #= 173 and Freq363 #= 1) or
	(Row363 #= 174 and Freq363 #= 1) or
	(Row363 #= 175 and Freq363 #= 1) or
	(Row363 #= 176 and Freq363 #= 1) or
	(Row363 #= 177 and Freq363 #= 1) or
	(Row363 #= 179 and Freq363 #= 1) or
	(Row363 #= 184 and Freq363 #= 1) or
	(Row363 #= 0 and Freq363 #= 0)), 

	((Row364 #= 146 and Freq364 #= 1) or
	(Row364 #= 149 and Freq364 #= 1) or
	(Row364 #= 151 and Freq364 #= 1) or
	(Row364 #= 152 and Freq364 #= 1) or
	(Row364 #= 153 and Freq364 #= 1) or
	(Row364 #= 154 and Freq364 #= 1) or
	(Row364 #= 155 and Freq364 #= 1) or
	(Row364 #= 156 and Freq364 #= 1) or
	(Row364 #= 157 and Freq364 #= 1) or
	(Row364 #= 158 and Freq364 #= 1) or
	(Row364 #= 159 and Freq364 #= 1) or
	(Row364 #= 166 and Freq364 #= 2) or
	(Row364 #= 167 and Freq364 #= 1) or
	(Row364 #= 168 and Freq364 #= 1) or
	(Row364 #= 169 and Freq364 #= 1) or
	(Row364 #= 170 and Freq364 #= 1) or
	(Row364 #= 171 and Freq364 #= 1) or
	(Row364 #= 172 and Freq364 #= 1) or
	(Row364 #= 173 and Freq364 #= 1) or
	(Row364 #= 174 and Freq364 #= 1) or
	(Row364 #= 175 and Freq364 #= 1) or
	(Row364 #= 176 and Freq364 #= 1) or
	(Row364 #= 177 and Freq364 #= 1) or
	(Row364 #= 178 and Freq364 #= 1) or
	(Row364 #= 179 and Freq364 #= 1) or
	(Row364 #= 180 and Freq364 #= 1) or
	(Row364 #= 181 and Freq364 #= 1) or
	(Row364 #= 187 and Freq364 #= 1) or
	(Row364 #= 0 and Freq364 #= 0)), 

	((Row365 #= 147 and Freq365 #= 1) or
	(Row365 #= 148 and Freq365 #= 1) or
	(Row365 #= 151 and Freq365 #= 1) or
	(Row365 #= 152 and Freq365 #= 1) or
	(Row365 #= 153 and Freq365 #= 1) or
	(Row365 #= 154 and Freq365 #= 1) or
	(Row365 #= 155 and Freq365 #= 1) or
	(Row365 #= 156 and Freq365 #= 1) or
	(Row365 #= 157 and Freq365 #= 1) or
	(Row365 #= 158 and Freq365 #= 1) or
	(Row365 #= 159 and Freq365 #= 1) or
	(Row365 #= 166 and Freq365 #= 2) or
	(Row365 #= 167 and Freq365 #= 2) or
	(Row365 #= 168 and Freq365 #= 1) or
	(Row365 #= 169 and Freq365 #= 2) or
	(Row365 #= 170 and Freq365 #= 1) or
	(Row365 #= 171 and Freq365 #= 1) or
	(Row365 #= 172 and Freq365 #= 1) or
	(Row365 #= 173 and Freq365 #= 1) or
	(Row365 #= 174 and Freq365 #= 1) or
	(Row365 #= 175 and Freq365 #= 1) or
	(Row365 #= 176 and Freq365 #= 1) or
	(Row365 #= 177 and Freq365 #= 1) or
	(Row365 #= 178 and Freq365 #= 1) or
	(Row365 #= 179 and Freq365 #= 1) or
	(Row365 #= 180 and Freq365 #= 1) or
	(Row365 #= 181 and Freq365 #= 1) or
	(Row365 #= 182 and Freq365 #= 1) or
	(Row365 #= 0 and Freq365 #= 0)), 

	((Row366 #= 142 and Freq366 #= 1) or
	(Row366 #= 145 and Freq366 #= 1) or
	(Row366 #= 147 and Freq366 #= 1) or
	(Row366 #= 148 and Freq366 #= 1) or
	(Row366 #= 150 and Freq366 #= 1) or
	(Row366 #= 151 and Freq366 #= 1) or
	(Row366 #= 152 and Freq366 #= 1) or
	(Row366 #= 153 and Freq366 #= 1) or
	(Row366 #= 154 and Freq366 #= 1) or
	(Row366 #= 155 and Freq366 #= 1) or
	(Row366 #= 156 and Freq366 #= 1) or
	(Row366 #= 157 and Freq366 #= 1) or
	(Row366 #= 158 and Freq366 #= 1) or
	(Row366 #= 159 and Freq366 #= 1) or
	(Row366 #= 166 and Freq366 #= 1) or
	(Row366 #= 167 and Freq366 #= 2) or
	(Row366 #= 168 and Freq366 #= 2) or
	(Row366 #= 169 and Freq366 #= 2) or
	(Row366 #= 170 and Freq366 #= 1) or
	(Row366 #= 171 and Freq366 #= 2) or
	(Row366 #= 172 and Freq366 #= 1) or
	(Row366 #= 173 and Freq366 #= 1) or
	(Row366 #= 174 and Freq366 #= 1) or
	(Row366 #= 175 and Freq366 #= 1) or
	(Row366 #= 176 and Freq366 #= 1) or
	(Row366 #= 177 and Freq366 #= 1) or
	(Row366 #= 178 and Freq366 #= 1) or
	(Row366 #= 179 and Freq366 #= 1) or
	(Row366 #= 180 and Freq366 #= 1) or
	(Row366 #= 181 and Freq366 #= 1) or
	(Row366 #= 182 and Freq366 #= 1) or
	(Row366 #= 183 and Freq366 #= 1) or
	(Row366 #= 184 and Freq366 #= 1) or
	(Row366 #= 0 and Freq366 #= 0)), 

	((Row367 #= 145 and Freq367 #= 1) or
	(Row367 #= 146 and Freq367 #= 1) or
	(Row367 #= 147 and Freq367 #= 1) or
	(Row367 #= 148 and Freq367 #= 1) or
	(Row367 #= 149 and Freq367 #= 1) or
	(Row367 #= 150 and Freq367 #= 1) or
	(Row367 #= 151 and Freq367 #= 1) or
	(Row367 #= 152 and Freq367 #= 1) or
	(Row367 #= 153 and Freq367 #= 1) or
	(Row367 #= 154 and Freq367 #= 1) or
	(Row367 #= 155 and Freq367 #= 1) or
	(Row367 #= 156 and Freq367 #= 1) or
	(Row367 #= 157 and Freq367 #= 1) or
	(Row367 #= 158 and Freq367 #= 2) or
	(Row367 #= 159 and Freq367 #= 1) or
	(Row367 #= 166 and Freq367 #= 2) or
	(Row367 #= 167 and Freq367 #= 2) or
	(Row367 #= 168 and Freq367 #= 2) or
	(Row367 #= 169 and Freq367 #= 2) or
	(Row367 #= 170 and Freq367 #= 2) or
	(Row367 #= 171 and Freq367 #= 2) or
	(Row367 #= 172 and Freq367 #= 1) or
	(Row367 #= 173 and Freq367 #= 1) or
	(Row367 #= 174 and Freq367 #= 1) or
	(Row367 #= 175 and Freq367 #= 1) or
	(Row367 #= 176 and Freq367 #= 1) or
	(Row367 #= 177 and Freq367 #= 1) or
	(Row367 #= 178 and Freq367 #= 1) or
	(Row367 #= 179 and Freq367 #= 1) or
	(Row367 #= 180 and Freq367 #= 1) or
	(Row367 #= 181 and Freq367 #= 1) or
	(Row367 #= 185 and Freq367 #= 1) or
	(Row367 #= 0 and Freq367 #= 0)), 

	((Row368 #= 145 and Freq368 #= 1) or
	(Row368 #= 146 and Freq368 #= 1) or
	(Row368 #= 147 and Freq368 #= 1) or
	(Row368 #= 148 and Freq368 #= 1) or
	(Row368 #= 149 and Freq368 #= 1) or
	(Row368 #= 150 and Freq368 #= 1) or
	(Row368 #= 151 and Freq368 #= 1) or
	(Row368 #= 152 and Freq368 #= 1) or
	(Row368 #= 153 and Freq368 #= 1) or
	(Row368 #= 154 and Freq368 #= 1) or
	(Row368 #= 155 and Freq368 #= 1) or
	(Row368 #= 156 and Freq368 #= 1) or
	(Row368 #= 157 and Freq368 #= 2) or
	(Row368 #= 158 and Freq368 #= 2) or
	(Row368 #= 159 and Freq368 #= 1) or
	(Row368 #= 166 and Freq368 #= 2) or
	(Row368 #= 167 and Freq368 #= 3) or
	(Row368 #= 168 and Freq368 #= 2) or
	(Row368 #= 169 and Freq368 #= 2) or
	(Row368 #= 170 and Freq368 #= 1) or
	(Row368 #= 171 and Freq368 #= 2) or
	(Row368 #= 172 and Freq368 #= 1) or
	(Row368 #= 173 and Freq368 #= 1) or
	(Row368 #= 174 and Freq368 #= 1) or
	(Row368 #= 175 and Freq368 #= 1) or
	(Row368 #= 176 and Freq368 #= 1) or
	(Row368 #= 177 and Freq368 #= 1) or
	(Row368 #= 178 and Freq368 #= 1) or
	(Row368 #= 179 and Freq368 #= 1) or
	(Row368 #= 180 and Freq368 #= 1) or
	(Row368 #= 181 and Freq368 #= 1) or
	(Row368 #= 184 and Freq368 #= 1) or
	(Row368 #= 185 and Freq368 #= 1) or
	(Row368 #= 0 and Freq368 #= 0)), 

	((Row369 #= 144 and Freq369 #= 1) or
	(Row369 #= 145 and Freq369 #= 1) or
	(Row369 #= 149 and Freq369 #= 1) or
	(Row369 #= 150 and Freq369 #= 1) or
	(Row369 #= 151 and Freq369 #= 1) or
	(Row369 #= 152 and Freq369 #= 1) or
	(Row369 #= 153 and Freq369 #= 1) or
	(Row369 #= 154 and Freq369 #= 1) or
	(Row369 #= 155 and Freq369 #= 1) or
	(Row369 #= 156 and Freq369 #= 1) or
	(Row369 #= 157 and Freq369 #= 2) or
	(Row369 #= 158 and Freq369 #= 2) or
	(Row369 #= 159 and Freq369 #= 2) or
	(Row369 #= 166 and Freq369 #= 3) or
	(Row369 #= 167 and Freq369 #= 2) or
	(Row369 #= 168 and Freq369 #= 2) or
	(Row369 #= 169 and Freq369 #= 2) or
	(Row369 #= 170 and Freq369 #= 2) or
	(Row369 #= 171 and Freq369 #= 2) or
	(Row369 #= 172 and Freq369 #= 1) or
	(Row369 #= 173 and Freq369 #= 1) or
	(Row369 #= 174 and Freq369 #= 1) or
	(Row369 #= 175 and Freq369 #= 1) or
	(Row369 #= 176 and Freq369 #= 1) or
	(Row369 #= 177 and Freq369 #= 1) or
	(Row369 #= 178 and Freq369 #= 1) or
	(Row369 #= 179 and Freq369 #= 1) or
	(Row369 #= 180 and Freq369 #= 1) or
	(Row369 #= 181 and Freq369 #= 1) or
	(Row369 #= 182 and Freq369 #= 1) or
	(Row369 #= 183 and Freq369 #= 1) or
	(Row369 #= 0 and Freq369 #= 0)), 

	((Row370 #= 148 and Freq370 #= 1) or
	(Row370 #= 149 and Freq370 #= 1) or
	(Row370 #= 150 and Freq370 #= 1) or
	(Row370 #= 151 and Freq370 #= 1) or
	(Row370 #= 152 and Freq370 #= 1) or
	(Row370 #= 153 and Freq370 #= 2) or
	(Row370 #= 154 and Freq370 #= 1) or
	(Row370 #= 155 and Freq370 #= 1) or
	(Row370 #= 156 and Freq370 #= 1) or
	(Row370 #= 157 and Freq370 #= 2) or
	(Row370 #= 158 and Freq370 #= 2) or
	(Row370 #= 159 and Freq370 #= 2) or
	(Row370 #= 166 and Freq370 #= 3) or
	(Row370 #= 167 and Freq370 #= 3) or
	(Row370 #= 168 and Freq370 #= 3) or
	(Row370 #= 169 and Freq370 #= 3) or
	(Row370 #= 170 and Freq370 #= 2) or
	(Row370 #= 171 and Freq370 #= 2) or
	(Row370 #= 172 and Freq370 #= 1) or
	(Row370 #= 173 and Freq370 #= 1) or
	(Row370 #= 174 and Freq370 #= 1) or
	(Row370 #= 175 and Freq370 #= 1) or
	(Row370 #= 176 and Freq370 #= 1) or
	(Row370 #= 177 and Freq370 #= 1) or
	(Row370 #= 178 and Freq370 #= 1) or
	(Row370 #= 179 and Freq370 #= 1) or
	(Row370 #= 180 and Freq370 #= 1) or
	(Row370 #= 181 and Freq370 #= 1) or
	(Row370 #= 182 and Freq370 #= 1) or
	(Row370 #= 183 and Freq370 #= 1) or
	(Row370 #= 0 and Freq370 #= 0)), 

	((Row371 #= 146 and Freq371 #= 1) or
	(Row371 #= 149 and Freq371 #= 1) or
	(Row371 #= 150 and Freq371 #= 1) or
	(Row371 #= 151 and Freq371 #= 1) or
	(Row371 #= 152 and Freq371 #= 1) or
	(Row371 #= 153 and Freq371 #= 1) or
	(Row371 #= 154 and Freq371 #= 1) or
	(Row371 #= 155 and Freq371 #= 1) or
	(Row371 #= 156 and Freq371 #= 1) or
	(Row371 #= 157 and Freq371 #= 2) or
	(Row371 #= 158 and Freq371 #= 2) or
	(Row371 #= 159 and Freq371 #= 2) or
	(Row371 #= 166 and Freq371 #= 4) or
	(Row371 #= 167 and Freq371 #= 4) or
	(Row371 #= 168 and Freq371 #= 3) or
	(Row371 #= 169 and Freq371 #= 3) or
	(Row371 #= 170 and Freq371 #= 3) or
	(Row371 #= 171 and Freq371 #= 2) or
	(Row371 #= 172 and Freq371 #= 2) or
	(Row371 #= 173 and Freq371 #= 1) or
	(Row371 #= 174 and Freq371 #= 1) or
	(Row371 #= 175 and Freq371 #= 1) or
	(Row371 #= 176 and Freq371 #= 1) or
	(Row371 #= 177 and Freq371 #= 1) or
	(Row371 #= 178 and Freq371 #= 1) or
	(Row371 #= 179 and Freq371 #= 1) or
	(Row371 #= 180 and Freq371 #= 1) or
	(Row371 #= 181 and Freq371 #= 1) or
	(Row371 #= 182 and Freq371 #= 1) or
	(Row371 #= 188 and Freq371 #= 1) or
	(Row371 #= 0 and Freq371 #= 0)), 

	((Row372 #= 142 and Freq372 #= 1) or
	(Row372 #= 146 and Freq372 #= 1) or
	(Row372 #= 147 and Freq372 #= 1) or
	(Row372 #= 148 and Freq372 #= 1) or
	(Row372 #= 149 and Freq372 #= 1) or
	(Row372 #= 150 and Freq372 #= 1) or
	(Row372 #= 151 and Freq372 #= 1) or
	(Row372 #= 152 and Freq372 #= 1) or
	(Row372 #= 153 and Freq372 #= 1) or
	(Row372 #= 154 and Freq372 #= 1) or
	(Row372 #= 155 and Freq372 #= 1) or
	(Row372 #= 156 and Freq372 #= 2) or
	(Row372 #= 157 and Freq372 #= 2) or
	(Row372 #= 158 and Freq372 #= 3) or
	(Row372 #= 159 and Freq372 #= 3) or
	(Row372 #= 166 and Freq372 #= 4) or
	(Row372 #= 167 and Freq372 #= 4) or
	(Row372 #= 168 and Freq372 #= 4) or
	(Row372 #= 169 and Freq372 #= 3) or
	(Row372 #= 170 and Freq372 #= 2) or
	(Row372 #= 171 and Freq372 #= 2) or
	(Row372 #= 172 and Freq372 #= 2) or
	(Row372 #= 173 and Freq372 #= 2) or
	(Row372 #= 174 and Freq372 #= 2) or
	(Row372 #= 175 and Freq372 #= 1) or
	(Row372 #= 176 and Freq372 #= 1) or
	(Row372 #= 177 and Freq372 #= 1) or
	(Row372 #= 178 and Freq372 #= 1) or
	(Row372 #= 179 and Freq372 #= 1) or
	(Row372 #= 180 and Freq372 #= 1) or
	(Row372 #= 181 and Freq372 #= 1) or
	(Row372 #= 182 and Freq372 #= 1) or
	(Row372 #= 0 and Freq372 #= 0)), 

	((Row373 #= 146 and Freq373 #= 1) or
	(Row373 #= 147 and Freq373 #= 1) or
	(Row373 #= 148 and Freq373 #= 1) or
	(Row373 #= 149 and Freq373 #= 1) or
	(Row373 #= 150 and Freq373 #= 1) or
	(Row373 #= 151 and Freq373 #= 1) or
	(Row373 #= 152 and Freq373 #= 1) or
	(Row373 #= 153 and Freq373 #= 1) or
	(Row373 #= 154 and Freq373 #= 2) or
	(Row373 #= 155 and Freq373 #= 2) or
	(Row373 #= 156 and Freq373 #= 2) or
	(Row373 #= 157 and Freq373 #= 2) or
	(Row373 #= 158 and Freq373 #= 3) or
	(Row373 #= 159 and Freq373 #= 3) or
	(Row373 #= 166 and Freq373 #= 5) or
	(Row373 #= 167 and Freq373 #= 5) or
	(Row373 #= 168 and Freq373 #= 5) or
	(Row373 #= 169 and Freq373 #= 4) or
	(Row373 #= 170 and Freq373 #= 3) or
	(Row373 #= 171 and Freq373 #= 3) or
	(Row373 #= 172 and Freq373 #= 2) or
	(Row373 #= 173 and Freq373 #= 2) or
	(Row373 #= 174 and Freq373 #= 2) or
	(Row373 #= 175 and Freq373 #= 2) or
	(Row373 #= 176 and Freq373 #= 1) or
	(Row373 #= 177 and Freq373 #= 2) or
	(Row373 #= 178 and Freq373 #= 2) or
	(Row373 #= 179 and Freq373 #= 1) or
	(Row373 #= 180 and Freq373 #= 1) or
	(Row373 #= 181 and Freq373 #= 1) or
	(Row373 #= 182 and Freq373 #= 1) or
	(Row373 #= 183 and Freq373 #= 1) or
	(Row373 #= 185 and Freq373 #= 1) or
	(Row373 #= 0 and Freq373 #= 0)), 

	((Row374 #= 144 and Freq374 #= 1) or
	(Row374 #= 146 and Freq374 #= 1) or
	(Row374 #= 147 and Freq374 #= 1) or
	(Row374 #= 148 and Freq374 #= 1) or
	(Row374 #= 149 and Freq374 #= 1) or
	(Row374 #= 150 and Freq374 #= 1) or
	(Row374 #= 151 and Freq374 #= 1) or
	(Row374 #= 152 and Freq374 #= 1) or
	(Row374 #= 153 and Freq374 #= 2) or
	(Row374 #= 154 and Freq374 #= 2) or
	(Row374 #= 155 and Freq374 #= 2) or
	(Row374 #= 156 and Freq374 #= 2) or
	(Row374 #= 157 and Freq374 #= 2) or
	(Row374 #= 158 and Freq374 #= 3) or
	(Row374 #= 159 and Freq374 #= 4) or
	(Row374 #= 166 and Freq374 #= 6) or
	(Row374 #= 167 and Freq374 #= 6) or
	(Row374 #= 168 and Freq374 #= 5) or
	(Row374 #= 169 and Freq374 #= 4) or
	(Row374 #= 170 and Freq374 #= 3) or
	(Row374 #= 171 and Freq374 #= 3) or
	(Row374 #= 172 and Freq374 #= 2) or
	(Row374 #= 173 and Freq374 #= 3) or
	(Row374 #= 174 and Freq374 #= 2) or
	(Row374 #= 175 and Freq374 #= 1) or
	(Row374 #= 176 and Freq374 #= 1) or
	(Row374 #= 177 and Freq374 #= 1) or
	(Row374 #= 178 and Freq374 #= 1) or
	(Row374 #= 179 and Freq374 #= 1) or
	(Row374 #= 183 and Freq374 #= 1) or
	(Row374 #= 186 and Freq374 #= 1) or
	(Row374 #= 187 and Freq374 #= 1) or
	(Row374 #= 209 and Freq374 #= 1) or
	(Row374 #= 0 and Freq374 #= 0)), 

	((Row375 #= 0 and Freq375 #= 0)), 

	((Row376 #= 0 and Freq376 #= 0)), 

	((Row377 #= 0 and Freq377 #= 0)), 

	((Row378 #= 0 and Freq378 #= 0)), 

	((Row379 #= 0 and Freq379 #= 0)), 

	((Row380 #= 148 and Freq380 #= 1) or
	(Row380 #= 150 and Freq380 #= 1) or
	(Row380 #= 151 and Freq380 #= 1) or
	(Row380 #= 152 and Freq380 #= 1) or
	(Row380 #= 153 and Freq380 #= 1) or
	(Row380 #= 154 and Freq380 #= 2) or
	(Row380 #= 155 and Freq380 #= 1) or
	(Row380 #= 156 and Freq380 #= 2) or
	(Row380 #= 157 and Freq380 #= 2) or
	(Row380 #= 158 and Freq380 #= 3) or
	(Row380 #= 159 and Freq380 #= 4) or
	(Row380 #= 166 and Freq380 #= 6) or
	(Row380 #= 167 and Freq380 #= 6) or
	(Row380 #= 168 and Freq380 #= 6) or
	(Row380 #= 169 and Freq380 #= 4) or
	(Row380 #= 170 and Freq380 #= 3) or
	(Row380 #= 171 and Freq380 #= 3) or
	(Row380 #= 172 and Freq380 #= 2) or
	(Row380 #= 173 and Freq380 #= 2) or
	(Row380 #= 174 and Freq380 #= 1) or
	(Row380 #= 175 and Freq380 #= 1) or
	(Row380 #= 176 and Freq380 #= 1) or
	(Row380 #= 177 and Freq380 #= 1) or
	(Row380 #= 178 and Freq380 #= 1) or
	(Row380 #= 184 and Freq380 #= 1) or
	(Row380 #= 0 and Freq380 #= 0)), 

	((Row381 #= 145 and Freq381 #= 1) or
	(Row381 #= 146 and Freq381 #= 1) or
	(Row381 #= 149 and Freq381 #= 1) or
	(Row381 #= 150 and Freq381 #= 1) or
	(Row381 #= 151 and Freq381 #= 1) or
	(Row381 #= 152 and Freq381 #= 1) or
	(Row381 #= 153 and Freq381 #= 1) or
	(Row381 #= 154 and Freq381 #= 2) or
	(Row381 #= 155 and Freq381 #= 2) or
	(Row381 #= 156 and Freq381 #= 3) or
	(Row381 #= 157 and Freq381 #= 3) or
	(Row381 #= 158 and Freq381 #= 4) or
	(Row381 #= 159 and Freq381 #= 4) or
	(Row381 #= 166 and Freq381 #= 6) or
	(Row381 #= 167 and Freq381 #= 6) or
	(Row381 #= 168 and Freq381 #= 5) or
	(Row381 #= 169 and Freq381 #= 5) or
	(Row381 #= 170 and Freq381 #= 3) or
	(Row381 #= 171 and Freq381 #= 4) or
	(Row381 #= 172 and Freq381 #= 2) or
	(Row381 #= 173 and Freq381 #= 2) or
	(Row381 #= 174 and Freq381 #= 2) or
	(Row381 #= 175 and Freq381 #= 2) or
	(Row381 #= 176 and Freq381 #= 2) or
	(Row381 #= 177 and Freq381 #= 2) or
	(Row381 #= 178 and Freq381 #= 1) or
	(Row381 #= 179 and Freq381 #= 1) or
	(Row381 #= 180 and Freq381 #= 1) or
	(Row381 #= 181 and Freq381 #= 1) or
	(Row381 #= 182 and Freq381 #= 1) or
	(Row381 #= 183 and Freq381 #= 1) or
	(Row381 #= 187 and Freq381 #= 1) or
	(Row381 #= 0 and Freq381 #= 0)), 

	((Row382 #= 147 and Freq382 #= 1) or
	(Row382 #= 148 and Freq382 #= 1) or
	(Row382 #= 149 and Freq382 #= 1) or
	(Row382 #= 150 and Freq382 #= 1) or
	(Row382 #= 151 and Freq382 #= 1) or
	(Row382 #= 152 and Freq382 #= 1) or
	(Row382 #= 153 and Freq382 #= 1) or
	(Row382 #= 154 and Freq382 #= 1) or
	(Row382 #= 155 and Freq382 #= 2) or
	(Row382 #= 156 and Freq382 #= 1) or
	(Row382 #= 157 and Freq382 #= 2) or
	(Row382 #= 158 and Freq382 #= 2) or
	(Row382 #= 159 and Freq382 #= 4) or
	(Row382 #= 166 and Freq382 #= 5) or
	(Row382 #= 167 and Freq382 #= 5) or
	(Row382 #= 168 and Freq382 #= 5) or
	(Row382 #= 169 and Freq382 #= 4) or
	(Row382 #= 170 and Freq382 #= 3) or
	(Row382 #= 171 and Freq382 #= 3) or
	(Row382 #= 172 and Freq382 #= 2) or
	(Row382 #= 173 and Freq382 #= 2) or
	(Row382 #= 174 and Freq382 #= 2) or
	(Row382 #= 175 and Freq382 #= 2) or
	(Row382 #= 176 and Freq382 #= 1) or
	(Row382 #= 177 and Freq382 #= 1) or
	(Row382 #= 178 and Freq382 #= 1) or
	(Row382 #= 179 and Freq382 #= 1) or
	(Row382 #= 180 and Freq382 #= 1) or
	(Row382 #= 181 and Freq382 #= 1) or
	(Row382 #= 0 and Freq382 #= 0)), 

	((Row383 #= 146 and Freq383 #= 1) or
	(Row383 #= 149 and Freq383 #= 1) or
	(Row383 #= 150 and Freq383 #= 1) or
	(Row383 #= 151 and Freq383 #= 1) or
	(Row383 #= 152 and Freq383 #= 1) or
	(Row383 #= 153 and Freq383 #= 1) or
	(Row383 #= 154 and Freq383 #= 1) or
	(Row383 #= 155 and Freq383 #= 2) or
	(Row383 #= 156 and Freq383 #= 1) or
	(Row383 #= 157 and Freq383 #= 2) or
	(Row383 #= 158 and Freq383 #= 2) or
	(Row383 #= 159 and Freq383 #= 3) or
	(Row383 #= 166 and Freq383 #= 4) or
	(Row383 #= 167 and Freq383 #= 5) or
	(Row383 #= 168 and Freq383 #= 4) or
	(Row383 #= 169 and Freq383 #= 3) or
	(Row383 #= 170 and Freq383 #= 3) or
	(Row383 #= 171 and Freq383 #= 2) or
	(Row383 #= 172 and Freq383 #= 2) or
	(Row383 #= 173 and Freq383 #= 2) or
	(Row383 #= 174 and Freq383 #= 1) or
	(Row383 #= 175 and Freq383 #= 1) or
	(Row383 #= 176 and Freq383 #= 1) or
	(Row383 #= 177 and Freq383 #= 1) or
	(Row383 #= 178 and Freq383 #= 1) or
	(Row383 #= 179 and Freq383 #= 1) or
	(Row383 #= 180 and Freq383 #= 1) or
	(Row383 #= 182 and Freq383 #= 1) or
	(Row383 #= 0 and Freq383 #= 0)), 

	((Row384 #= 146 and Freq384 #= 1) or
	(Row384 #= 150 and Freq384 #= 1) or
	(Row384 #= 151 and Freq384 #= 1) or
	(Row384 #= 152 and Freq384 #= 1) or
	(Row384 #= 153 and Freq384 #= 1) or
	(Row384 #= 154 and Freq384 #= 1) or
	(Row384 #= 155 and Freq384 #= 2) or
	(Row384 #= 156 and Freq384 #= 2) or
	(Row384 #= 157 and Freq384 #= 1) or
	(Row384 #= 158 and Freq384 #= 2) or
	(Row384 #= 159 and Freq384 #= 3) or
	(Row384 #= 166 and Freq384 #= 4) or
	(Row384 #= 167 and Freq384 #= 4) or
	(Row384 #= 168 and Freq384 #= 3) or
	(Row384 #= 169 and Freq384 #= 3) or
	(Row384 #= 170 and Freq384 #= 2) or
	(Row384 #= 171 and Freq384 #= 2) or
	(Row384 #= 172 and Freq384 #= 1) or
	(Row384 #= 173 and Freq384 #= 1) or
	(Row384 #= 174 and Freq384 #= 1) or
	(Row384 #= 175 and Freq384 #= 1) or
	(Row384 #= 176 and Freq384 #= 1) or
	(Row384 #= 177 and Freq384 #= 1) or
	(Row384 #= 178 and Freq384 #= 1) or
	(Row384 #= 179 and Freq384 #= 1) or
	(Row384 #= 183 and Freq384 #= 1) or
	(Row384 #= 185 and Freq384 #= 1) or
	(Row384 #= 0 and Freq384 #= 0)), 

	((Row385 #= 144 and Freq385 #= 1) or
	(Row385 #= 146 and Freq385 #= 1) or
	(Row385 #= 147 and Freq385 #= 1) or
	(Row385 #= 148 and Freq385 #= 1) or
	(Row385 #= 149 and Freq385 #= 1) or
	(Row385 #= 150 and Freq385 #= 1) or
	(Row385 #= 151 and Freq385 #= 1) or
	(Row385 #= 152 and Freq385 #= 1) or
	(Row385 #= 153 and Freq385 #= 1) or
	(Row385 #= 154 and Freq385 #= 1) or
	(Row385 #= 155 and Freq385 #= 1) or
	(Row385 #= 156 and Freq385 #= 2) or
	(Row385 #= 157 and Freq385 #= 2) or
	(Row385 #= 158 and Freq385 #= 2) or
	(Row385 #= 159 and Freq385 #= 2) or
	(Row385 #= 166 and Freq385 #= 3) or
	(Row385 #= 167 and Freq385 #= 4) or
	(Row385 #= 168 and Freq385 #= 4) or
	(Row385 #= 169 and Freq385 #= 3) or
	(Row385 #= 170 and Freq385 #= 2) or
	(Row385 #= 171 and Freq385 #= 2) or
	(Row385 #= 172 and Freq385 #= 2) or
	(Row385 #= 173 and Freq385 #= 2) or
	(Row385 #= 174 and Freq385 #= 1) or
	(Row385 #= 175 and Freq385 #= 2) or
	(Row385 #= 176 and Freq385 #= 1) or
	(Row385 #= 177 and Freq385 #= 1) or
	(Row385 #= 178 and Freq385 #= 1) or
	(Row385 #= 179 and Freq385 #= 1) or
	(Row385 #= 180 and Freq385 #= 1) or
	(Row385 #= 181 and Freq385 #= 1) or
	(Row385 #= 182 and Freq385 #= 1) or
	(Row385 #= 183 and Freq385 #= 1) or
	(Row385 #= 0 and Freq385 #= 0)), 

	((Row386 #= 149 and Freq386 #= 1) or
	(Row386 #= 150 and Freq386 #= 1) or
	(Row386 #= 151 and Freq386 #= 1) or
	(Row386 #= 152 and Freq386 #= 1) or
	(Row386 #= 153 and Freq386 #= 1) or
	(Row386 #= 154 and Freq386 #= 1) or
	(Row386 #= 155 and Freq386 #= 1) or
	(Row386 #= 156 and Freq386 #= 2) or
	(Row386 #= 157 and Freq386 #= 2) or
	(Row386 #= 158 and Freq386 #= 2) or
	(Row386 #= 159 and Freq386 #= 2) or
	(Row386 #= 166 and Freq386 #= 3) or
	(Row386 #= 167 and Freq386 #= 3) or
	(Row386 #= 168 and Freq386 #= 2) or
	(Row386 #= 169 and Freq386 #= 3) or
	(Row386 #= 170 and Freq386 #= 2) or
	(Row386 #= 171 and Freq386 #= 2) or
	(Row386 #= 172 and Freq386 #= 2) or
	(Row386 #= 173 and Freq386 #= 2) or
	(Row386 #= 174 and Freq386 #= 1) or
	(Row386 #= 175 and Freq386 #= 1) or
	(Row386 #= 176 and Freq386 #= 1) or
	(Row386 #= 177 and Freq386 #= 1) or
	(Row386 #= 178 and Freq386 #= 1) or
	(Row386 #= 179 and Freq386 #= 1) or
	(Row386 #= 180 and Freq386 #= 1) or
	(Row386 #= 181 and Freq386 #= 1) or
	(Row386 #= 182 and Freq386 #= 1) or
	(Row386 #= 0 and Freq386 #= 0)), 

	((Row387 #= 149 and Freq387 #= 1) or
	(Row387 #= 151 and Freq387 #= 1) or
	(Row387 #= 152 and Freq387 #= 1) or
	(Row387 #= 153 and Freq387 #= 1) or
	(Row387 #= 154 and Freq387 #= 1) or
	(Row387 #= 155 and Freq387 #= 1) or
	(Row387 #= 156 and Freq387 #= 1) or
	(Row387 #= 157 and Freq387 #= 1) or
	(Row387 #= 158 and Freq387 #= 2) or
	(Row387 #= 159 and Freq387 #= 1) or
	(Row387 #= 166 and Freq387 #= 3) or
	(Row387 #= 167 and Freq387 #= 3) or
	(Row387 #= 168 and Freq387 #= 2) or
	(Row387 #= 169 and Freq387 #= 2) or
	(Row387 #= 170 and Freq387 #= 2) or
	(Row387 #= 171 and Freq387 #= 2) or
	(Row387 #= 172 and Freq387 #= 2) or
	(Row387 #= 173 and Freq387 #= 1) or
	(Row387 #= 174 and Freq387 #= 1) or
	(Row387 #= 175 and Freq387 #= 1) or
	(Row387 #= 176 and Freq387 #= 2) or
	(Row387 #= 177 and Freq387 #= 1) or
	(Row387 #= 178 and Freq387 #= 1) or
	(Row387 #= 179 and Freq387 #= 1) or
	(Row387 #= 180 and Freq387 #= 1) or
	(Row387 #= 187 and Freq387 #= 1) or
	(Row387 #= 0 and Freq387 #= 0)), 

	((Row388 #= 142 and Freq388 #= 1) or
	(Row388 #= 144 and Freq388 #= 1) or
	(Row388 #= 150 and Freq388 #= 1) or
	(Row388 #= 151 and Freq388 #= 1) or
	(Row388 #= 152 and Freq388 #= 1) or
	(Row388 #= 153 and Freq388 #= 1) or
	(Row388 #= 154 and Freq388 #= 1) or
	(Row388 #= 155 and Freq388 #= 1) or
	(Row388 #= 156 and Freq388 #= 1) or
	(Row388 #= 158 and Freq388 #= 1) or
	(Row388 #= 159 and Freq388 #= 1) or
	(Row388 #= 166 and Freq388 #= 2) or
	(Row388 #= 167 and Freq388 #= 2) or
	(Row388 #= 168 and Freq388 #= 2) or
	(Row388 #= 169 and Freq388 #= 2) or
	(Row388 #= 170 and Freq388 #= 1) or
	(Row388 #= 171 and Freq388 #= 1) or
	(Row388 #= 172 and Freq388 #= 1) or
	(Row388 #= 173 and Freq388 #= 1) or
	(Row388 #= 174 and Freq388 #= 1) or
	(Row388 #= 175 and Freq388 #= 1) or
	(Row388 #= 176 and Freq388 #= 2) or
	(Row388 #= 177 and Freq388 #= 1) or
	(Row388 #= 178 and Freq388 #= 1) or
	(Row388 #= 182 and Freq388 #= 1) or
	(Row388 #= 183 and Freq388 #= 1) or
	(Row388 #= 0 and Freq388 #= 0)), 

	((Row389 #= 150 and Freq389 #= 1) or
	(Row389 #= 152 and Freq389 #= 1) or
	(Row389 #= 153 and Freq389 #= 1) or
	(Row389 #= 154 and Freq389 #= 1) or
	(Row389 #= 155 and Freq389 #= 1) or
	(Row389 #= 156 and Freq389 #= 1) or
	(Row389 #= 157 and Freq389 #= 1) or
	(Row389 #= 158 and Freq389 #= 1) or
	(Row389 #= 159 and Freq389 #= 1) or
	(Row389 #= 166 and Freq389 #= 2) or
	(Row389 #= 167 and Freq389 #= 2) or
	(Row389 #= 168 and Freq389 #= 2) or
	(Row389 #= 169 and Freq389 #= 2) or
	(Row389 #= 170 and Freq389 #= 1) or
	(Row389 #= 171 and Freq389 #= 1) or
	(Row389 #= 172 and Freq389 #= 1) or
	(Row389 #= 173 and Freq389 #= 1) or
	(Row389 #= 174 and Freq389 #= 1) or
	(Row389 #= 175 and Freq389 #= 1) or
	(Row389 #= 176 and Freq389 #= 1) or
	(Row389 #= 177 and Freq389 #= 1) or
	(Row389 #= 178 and Freq389 #= 1) or
	(Row389 #= 182 and Freq389 #= 1) or
	(Row389 #= 183 and Freq389 #= 1) or
	(Row389 #= 0 and Freq389 #= 0)), 

	((Row390 #= 152 and Freq390 #= 1) or
	(Row390 #= 154 and Freq390 #= 1) or
	(Row390 #= 155 and Freq390 #= 1) or
	(Row390 #= 156 and Freq390 #= 1) or
	(Row390 #= 157 and Freq390 #= 1) or
	(Row390 #= 158 and Freq390 #= 1) or
	(Row390 #= 159 and Freq390 #= 1) or
	(Row390 #= 166 and Freq390 #= 2) or
	(Row390 #= 167 and Freq390 #= 2) or
	(Row390 #= 168 and Freq390 #= 2) or
	(Row390 #= 169 and Freq390 #= 2) or
	(Row390 #= 170 and Freq390 #= 1) or
	(Row390 #= 171 and Freq390 #= 1) or
	(Row390 #= 172 and Freq390 #= 1) or
	(Row390 #= 173 and Freq390 #= 1) or
	(Row390 #= 174 and Freq390 #= 1) or
	(Row390 #= 175 and Freq390 #= 1) or
	(Row390 #= 176 and Freq390 #= 1) or
	(Row390 #= 177 and Freq390 #= 1) or
	(Row390 #= 179 and Freq390 #= 1) or
	(Row390 #= 182 and Freq390 #= 1) or
	(Row390 #= 183 and Freq390 #= 1) or
	(Row390 #= 186 and Freq390 #= 1) or
	(Row390 #= 0 and Freq390 #= 0)), 

	((Row391 #= 145 and Freq391 #= 1) or
	(Row391 #= 149 and Freq391 #= 1) or
	(Row391 #= 151 and Freq391 #= 1) or
	(Row391 #= 152 and Freq391 #= 1) or
	(Row391 #= 153 and Freq391 #= 1) or
	(Row391 #= 154 and Freq391 #= 1) or
	(Row391 #= 155 and Freq391 #= 1) or
	(Row391 #= 157 and Freq391 #= 1) or
	(Row391 #= 158 and Freq391 #= 1) or
	(Row391 #= 159 and Freq391 #= 1) or
	(Row391 #= 166 and Freq391 #= 1) or
	(Row391 #= 167 and Freq391 #= 1) or
	(Row391 #= 168 and Freq391 #= 1) or
	(Row391 #= 169 and Freq391 #= 2) or
	(Row391 #= 170 and Freq391 #= 1) or
	(Row391 #= 171 and Freq391 #= 1) or
	(Row391 #= 172 and Freq391 #= 1) or
	(Row391 #= 173 and Freq391 #= 1) or
	(Row391 #= 174 and Freq391 #= 1) or
	(Row391 #= 175 and Freq391 #= 1) or
	(Row391 #= 176 and Freq391 #= 1) or
	(Row391 #= 177 and Freq391 #= 1) or
	(Row391 #= 178 and Freq391 #= 1) or
	(Row391 #= 179 and Freq391 #= 1) or
	(Row391 #= 180 and Freq391 #= 1) or
	(Row391 #= 181 and Freq391 #= 1) or
	(Row391 #= 182 and Freq391 #= 1) or
	(Row391 #= 183 and Freq391 #= 1) or
	(Row391 #= 185 and Freq391 #= 1) or
	(Row391 #= 186 and Freq391 #= 1) or
	(Row391 #= 0 and Freq391 #= 0)), 

	((Row392 #= 151 and Freq392 #= 1) or
	(Row392 #= 153 and Freq392 #= 1) or
	(Row392 #= 155 and Freq392 #= 1) or
	(Row392 #= 156 and Freq392 #= 1) or
	(Row392 #= 157 and Freq392 #= 1) or
	(Row392 #= 158 and Freq392 #= 1) or
	(Row392 #= 159 and Freq392 #= 1) or
	(Row392 #= 166 and Freq392 #= 1) or
	(Row392 #= 167 and Freq392 #= 1) or
	(Row392 #= 168 and Freq392 #= 1) or
	(Row392 #= 169 and Freq392 #= 1) or
	(Row392 #= 170 and Freq392 #= 1) or
	(Row392 #= 171 and Freq392 #= 1) or
	(Row392 #= 173 and Freq392 #= 1) or
	(Row392 #= 174 and Freq392 #= 1) or
	(Row392 #= 176 and Freq392 #= 1) or
	(Row392 #= 177 and Freq392 #= 1) or
	(Row392 #= 178 and Freq392 #= 1) or
	(Row392 #= 180 and Freq392 #= 1) or
	(Row392 #= 181 and Freq392 #= 1) or
	(Row392 #= 183 and Freq392 #= 1) or
	(Row392 #= 0 and Freq392 #= 0)), 

	((Row393 #= 150 and Freq393 #= 1) or
	(Row393 #= 155 and Freq393 #= 1) or
	(Row393 #= 156 and Freq393 #= 1) or
	(Row393 #= 157 and Freq393 #= 1) or
	(Row393 #= 158 and Freq393 #= 1) or
	(Row393 #= 166 and Freq393 #= 1) or
	(Row393 #= 167 and Freq393 #= 1) or
	(Row393 #= 168 and Freq393 #= 1) or
	(Row393 #= 169 and Freq393 #= 1) or
	(Row393 #= 170 and Freq393 #= 1) or
	(Row393 #= 171 and Freq393 #= 1) or
	(Row393 #= 172 and Freq393 #= 1) or
	(Row393 #= 173 and Freq393 #= 1) or
	(Row393 #= 175 and Freq393 #= 1) or
	(Row393 #= 176 and Freq393 #= 1) or
	(Row393 #= 177 and Freq393 #= 1) or
	(Row393 #= 178 and Freq393 #= 1) or
	(Row393 #= 179 and Freq393 #= 1) or
	(Row393 #= 180 and Freq393 #= 1) or
	(Row393 #= 182 and Freq393 #= 1) or
	(Row393 #= 196 and Freq393 #= 1) or
	(Row393 #= 0 and Freq393 #= 0)), 

	((Row394 #= 153 and Freq394 #= 1) or
	(Row394 #= 155 and Freq394 #= 1) or
	(Row394 #= 157 and Freq394 #= 1) or
	(Row394 #= 158 and Freq394 #= 1) or
	(Row394 #= 159 and Freq394 #= 1) or
	(Row394 #= 166 and Freq394 #= 1) or
	(Row394 #= 167 and Freq394 #= 1) or
	(Row394 #= 168 and Freq394 #= 1) or
	(Row394 #= 169 and Freq394 #= 1) or
	(Row394 #= 170 and Freq394 #= 1) or
	(Row394 #= 171 and Freq394 #= 1) or
	(Row394 #= 172 and Freq394 #= 1) or
	(Row394 #= 173 and Freq394 #= 1) or
	(Row394 #= 174 and Freq394 #= 1) or
	(Row394 #= 175 and Freq394 #= 1) or
	(Row394 #= 176 and Freq394 #= 1) or
	(Row394 #= 177 and Freq394 #= 1) or
	(Row394 #= 183 and Freq394 #= 1) or
	(Row394 #= 188 and Freq394 #= 1) or
	(Row394 #= 0 and Freq394 #= 0)), 

	((Row395 #= 151 and Freq395 #= 1) or
	(Row395 #= 155 and Freq395 #= 1) or
	(Row395 #= 156 and Freq395 #= 1) or
	(Row395 #= 157 and Freq395 #= 1) or
	(Row395 #= 158 and Freq395 #= 1) or
	(Row395 #= 159 and Freq395 #= 1) or
	(Row395 #= 166 and Freq395 #= 1) or
	(Row395 #= 167 and Freq395 #= 1) or
	(Row395 #= 168 and Freq395 #= 1) or
	(Row395 #= 169 and Freq395 #= 1) or
	(Row395 #= 170 and Freq395 #= 1) or
	(Row395 #= 171 and Freq395 #= 1) or
	(Row395 #= 172 and Freq395 #= 1) or
	(Row395 #= 174 and Freq395 #= 1) or
	(Row395 #= 175 and Freq395 #= 1) or
	(Row395 #= 176 and Freq395 #= 1) or
	(Row395 #= 177 and Freq395 #= 1) or
	(Row395 #= 178 and Freq395 #= 1) or
	(Row395 #= 0 and Freq395 #= 0)), 

	((Row396 #= 152 and Freq396 #= 1) or
	(Row396 #= 154 and Freq396 #= 1) or
	(Row396 #= 156 and Freq396 #= 1) or
	(Row396 #= 157 and Freq396 #= 1) or
	(Row396 #= 158 and Freq396 #= 1) or
	(Row396 #= 166 and Freq396 #= 1) or
	(Row396 #= 167 and Freq396 #= 1) or
	(Row396 #= 168 and Freq396 #= 1) or
	(Row396 #= 169 and Freq396 #= 1) or
	(Row396 #= 171 and Freq396 #= 1) or
	(Row396 #= 172 and Freq396 #= 1) or
	(Row396 #= 173 and Freq396 #= 1) or
	(Row396 #= 175 and Freq396 #= 1) or
	(Row396 #= 176 and Freq396 #= 1) or
	(Row396 #= 177 and Freq396 #= 1) or
	(Row396 #= 178 and Freq396 #= 1) or
	(Row396 #= 179 and Freq396 #= 1) or
	(Row396 #= 180 and Freq396 #= 1) or
	(Row396 #= 181 and Freq396 #= 1) or
	(Row396 #= 183 and Freq396 #= 1) or
	(Row396 #= 0 and Freq396 #= 0)), 

	((Row397 #= 152 and Freq397 #= 1) or
	(Row397 #= 157 and Freq397 #= 1) or
	(Row397 #= 158 and Freq397 #= 1) or
	(Row397 #= 167 and Freq397 #= 1) or
	(Row397 #= 169 and Freq397 #= 1) or
	(Row397 #= 173 and Freq397 #= 1) or
	(Row397 #= 175 and Freq397 #= 1) or
	(Row397 #= 177 and Freq397 #= 1) or
	(Row397 #= 179 and Freq397 #= 1) or
	(Row397 #= 180 and Freq397 #= 1) or
	(Row397 #= 186 and Freq397 #= 1) or
	(Row397 #= 0 and Freq397 #= 0)), 

	((Row398 #= 157 and Freq398 #= 1) or
	(Row398 #= 169 and Freq398 #= 1) or
	(Row398 #= 0 and Freq398 #= 0)), 

	((Row399 #= 157 and Freq399 #= 1) or
	(Row399 #= 158 and Freq399 #= 1) or
	(Row399 #= 167 and Freq399 #= 1) or
	(Row399 #= 169 and Freq399 #= 1) or
	(Row399 #= 172 and Freq399 #= 1) or
	(Row399 #= 175 and Freq399 #= 1) or
	(Row399 #= 176 and Freq399 #= 1) or
	(Row399 #= 181 and Freq399 #= 1) or
	(Row399 #= 0 and Freq399 #= 0)), 

	((Row400 #= 167 and Freq400 #= 1) or
	(Row400 #= 175 and Freq400 #= 1) or
	(Row400 #= 182 and Freq400 #= 1) or
	(Row400 #= 0 and Freq400 #= 0)), 

	((Row401 #= 153 and Freq401 #= 1) or
	(Row401 #= 169 and Freq401 #= 1) or
	(Row401 #= 175 and Freq401 #= 1) or
	(Row401 #= 176 and Freq401 #= 1) or
	(Row401 #= 181 and Freq401 #= 1) or
	(Row401 #= 182 and Freq401 #= 1) or
	(Row401 #= 0 and Freq401 #= 0)), 

	((Row402 #= 176 and Freq402 #= 1) or
	(Row402 #= 0 and Freq402 #= 0)), 

	((Row403 #= 195 and Freq403 #= 1) or
	(Row403 #= 0 and Freq403 #= 0)), 

	((Row404 #= 0 and Freq404 #= 0)), 

	((Row405 #= 0 and Freq405 #= 0)), 

	((Row406 #= 0 and Freq406 #= 0)), 

	((Row407 #= 0 and Freq407 #= 0)), 

	((Row408 #= 0 and Freq408 #= 0)), 

	((Row409 #= 0 and Freq409 #= 0)), 

	((Row410 #= 0 and Freq410 #= 0)), 

	((Row411 #= 0 and Freq411 #= 0)), 

	((Row412 #= 0 and Freq412 #= 0)), 

	((Row413 #= 0 and Freq413 #= 0)), 

	((Row414 #= 0 and Freq414 #= 0)), 

	((Row415 #= 0 and Freq415 #= 0)), 

	((Row416 #= 0 and Freq416 #= 0)), 

	((Row417 #= 0 and Freq417 #= 0)), 

	((Row418 #= 0 and Freq418 #= 0)), 

	((Row419 #= 0 and Freq419 #= 0)), 

	((Row420 #= 0 and Freq420 #= 0)), 

	((Row421 #= 0 and Freq421 #= 0)), 

	((Row422 #= 0 and Freq422 #= 0)), 

	((Row423 #= 0 and Freq423 #= 0)), 

	((Row424 #= 151 and Freq424 #= 1) or
	(Row424 #= 0 and Freq424 #= 0)), 

	((Row425 #= 0 and Freq425 #= 0)), 

	((Row426 #= 0 and Freq426 #= 0)), 

	((Row427 #= 0 and Freq427 #= 0)), 

	((Row428 #= 0 and Freq428 #= 0)), 

	((Row429 #= 0 and Freq429 #= 0)), 

	((Row430 #= 0 and Freq430 #= 0)), 

	((Row431 #= 0 and Freq431 #= 0)), 

	((Row432 #= 0 and Freq432 #= 0)), 

	((Row433 #= 0 and Freq433 #= 0)), 

	((Row434 #= 0 and Freq434 #= 0)), 

	((Row435 #= 0 and Freq435 #= 0)), 

	((Row436 #= 0 and Freq436 #= 0)), 

	((Row437 #= 0 and Freq437 #= 0)), 

	((Row438 #= 0 and Freq438 #= 0)), 

	((Row439 #= 0 and Freq439 #= 0)), 

	((Row440 #= 0 and Freq440 #= 0)), 

	((Row441 #= 0 and Freq441 #= 0)), 

	((Row442 #= 0 and Freq442 #= 0)), 

	((Row443 #= 0 and Freq443 #= 0)), 

	((Row444 #= 0 and Freq444 #= 0)), 

	((Row445 #= 0 and Freq445 #= 0)), 

	((Row446 #= 0 and Freq446 #= 0)), 

	((Row447 #= 0 and Freq447 #= 0)), 

	((Row448 #= 0 and Freq448 #= 0)), 

	((Row449 #= 0 and Freq449 #= 0)), 

	((Row450 #= 0 and Freq450 #= 0)), 

	((Row451 #= 0 and Freq451 #= 0)), 

	((Row452 #= 0 and Freq452 #= 0)), 

	((Row453 #= 0 and Freq453 #= 0)), 

	((Row454 #= 0 and Freq454 #= 0)), 

	((Row455 #= 0 and Freq455 #= 0)), 

	((Row456 #= 0 and Freq456 #= 0)), 

	((Row457 #= 0 and Freq457 #= 0)), 

	((Row458 #= 0 and Freq458 #= 0)), 

	((Row459 #= 0 and Freq459 #= 0)), 

	((Row460 #= 0 and Freq460 #= 0)), 

	((Row461 #= 0 and Freq461 #= 0)), 

	((Row462 #= 4 and Freq462 #= 7) or
	(Row462 #= 0 and Freq462 #= 0)), 

	((Row463 #= 0 and Freq463 #= 0)), 

	((Row464 #= 0 and Freq464 #= 0)), 

	((Row465 #= 0 and Freq465 #= 0)), 

	((Row466 #= 0 and Freq466 #= 0)), 

	((Row467 #= 0 and Freq467 #= 0)), 

	((Row468 #= 0 and Freq468 #= 0)), 

	((Row469 #= 0 and Freq469 #= 0)), 

	((Row470 #= 0 and Freq470 #= 0)), 

	((Row471 #= 0 and Freq471 #= 0)), 

	((Row472 #= 0 and Freq472 #= 0)), 

	((Row473 #= 0 and Freq473 #= 0)), 

	((Row474 #= 0 and Freq474 #= 0)), 

	((Row475 #= 0 and Freq475 #= 0)), 

	((Row476 #= 0 and Freq476 #= 0)), 

	((Row477 #= 0 and Freq477 #= 0)), 

	((Row478 #= 0 and Freq478 #= 0)), 

	((Row479 #= 0 and Freq479 #= 0)), 

	((Row480 #= 0 and Freq480 #= 0)), 

	((Row481 #= 0 and Freq481 #= 0)), 

	((Row482 #= 0 and Freq482 #= 0)), 

	((Row483 #= 0 and Freq483 #= 0)), 

	((Row484 #= 0 and Freq484 #= 0)), 

	((Row485 #= 0 and Freq485 #= 0)), 

	((Row486 #= 0 and Freq486 #= 0)), 

	((Row487 #= 0 and Freq487 #= 0)), 

	((Row488 #= 0 and Freq488 #= 0)), 

	((Row489 #= 0 and Freq489 #= 0)), 

	((Row490 #= 0 and Freq490 #= 0)), 

	((Row491 #= 0 and Freq491 #= 0)), 

	((Row492 #= 0 and Freq492 #= 0)), 

	((Row493 #= 0 and Freq493 #= 0)), 

	((Row494 #= 0 and Freq494 #= 0)), 

	((Row495 #= 0 and Freq495 #= 0)), 

	((Row496 #= 0 and Freq496 #= 0)), 

	((Row497 #= 0 and Freq497 #= 0)), 

	((Row498 #= 0 and Freq498 #= 0)), 

	((Row499 #= 0 and Freq499 #= 0)), 

	((Row500 #= 0 and Freq500 #= 0)), 

	((Row501 #= 0 and Freq501 #= 0)), 

	((Row502 #= 0 and Freq502 #= 0)), 

	((Row503 #= 0 and Freq503 #= 0)), 

	((Row504 #= 0 and Freq504 #= 0)), 

	((Row505 #= 0 and Freq505 #= 0)), 

	((Row506 #= 0 and Freq506 #= 0)), 

	((Row507 #= 0 and Freq507 #= 0)), 

	((Row508 #= 0 and Freq508 #= 0)), 

	((Row509 #= 0 and Freq509 #= 0)), 

	((Row510 #= 0 and Freq510 #= 0)), 

	((Row511 #= 0 and Freq511 #= 0)), 

	((Row512 #= 0 and Freq512 #= 0)), 

	((Row513 #= 0 and Freq513 #= 0)), 

	((Row514 #= 0 and Freq514 #= 0)), 

	((Row515 #= 0 and Freq515 #= 0)), 

	((Row516 #= 0 and Freq516 #= 0)), 

	((Row517 #= 0 and Freq517 #= 0)), 

	((Row518 #= 0 and Freq518 #= 0)), 

	((Row519 #= 0 and Freq519 #= 0)), 

	((Row520 #= 0 and Freq520 #= 0)), 

	((Row521 #= 0 and Freq521 #= 0)), 

	((Row522 #= 0 and Freq522 #= 0)), 

	((Row523 #= 0 and Freq523 #= 0)), 

	((Row524 #= 0 and Freq524 #= 0)), 

	((Row525 #= 0 and Freq525 #= 0)), 

	((Row526 #= 0 and Freq526 #= 0)), 

	((Row527 #= 0 and Freq527 #= 0)), 

	((Row528 #= 0 and Freq528 #= 0)), 

	((Row529 #= 0 and Freq529 #= 0)), 

	((Row530 #= 0 and Freq530 #= 0)), 

	((Row531 #= 11 and Freq531 #= 1) or
	(Row531 #= 0 and Freq531 #= 0)), 

	((Row532 #= 0 and Freq532 #= 0)), 

	((Row533 #= 11 and Freq533 #= 1) or
	(Row533 #= 0 and Freq533 #= 0)), 

	((Row534 #= 0 and Freq534 #= 0)), 

	((Row535 #= 11 and Freq535 #= 1) or
	(Row535 #= 12 and Freq535 #= 1) or
	(Row535 #= 0 and Freq535 #= 0)), 

	((Row536 #= 10 and Freq536 #= 1) or
	(Row536 #= 19 and Freq536 #= 1) or
	(Row536 #= 0 and Freq536 #= 0)), 

	((Row537 #= 6 and Freq537 #= 1) or
	(Row537 #= 14 and Freq537 #= 1) or
	(Row537 #= 450 and Freq537 #= 1) or
	(Row537 #= 0 and Freq537 #= 0)), 

	((Row538 #= 6 and Freq538 #= 1) or
	(Row538 #= 7 and Freq538 #= 1) or
	(Row538 #= 11 and Freq538 #= 1) or
	(Row538 #= 14 and Freq538 #= 1) or
	(Row538 #= 15 and Freq538 #= 1) or
	(Row538 #= 440 and Freq538 #= 1) or
	(Row538 #= 442 and Freq538 #= 1) or
	(Row538 #= 445 and Freq538 #= 1) or
	(Row538 #= 446 and Freq538 #= 1) or
	(Row538 #= 449 and Freq538 #= 1) or
	(Row538 #= 0 and Freq538 #= 0)), 

	((Row539 #= 8 and Freq539 #= 1) or
	(Row539 #= 13 and Freq539 #= 1) or
	(Row539 #= 14 and Freq539 #= 1) or
	(Row539 #= 16 and Freq539 #= 1) or
	(Row539 #= 18 and Freq539 #= 1) or
	(Row539 #= 438 and Freq539 #= 1) or
	(Row539 #= 439 and Freq539 #= 1) or
	(Row539 #= 440 and Freq539 #= 1) or
	(Row539 #= 441 and Freq539 #= 1) or
	(Row539 #= 0 and Freq539 #= 0)), 

	((Row540 #= 4 and Freq540 #= 1) or
	(Row540 #= 5 and Freq540 #= 1) or
	(Row540 #= 6 and Freq540 #= 1) or
	(Row540 #= 7 and Freq540 #= 1) or
	(Row540 #= 10 and Freq540 #= 1) or
	(Row540 #= 11 and Freq540 #= 1) or
	(Row540 #= 12 and Freq540 #= 1) or
	(Row540 #= 13 and Freq540 #= 1) or
	(Row540 #= 14 and Freq540 #= 1) or
	(Row540 #= 15 and Freq540 #= 1) or
	(Row540 #= 16 and Freq540 #= 1) or
	(Row540 #= 18 and Freq540 #= 1) or
	(Row540 #= 20 and Freq540 #= 1) or
	(Row540 #= 439 and Freq540 #= 1) or
	(Row540 #= 440 and Freq540 #= 1) or
	(Row540 #= 441 and Freq540 #= 1) or
	(Row540 #= 442 and Freq540 #= 1) or
	(Row540 #= 443 and Freq540 #= 1) or
	(Row540 #= 444 and Freq540 #= 1) or
	(Row540 #= 445 and Freq540 #= 1) or
	(Row540 #= 446 and Freq540 #= 1) or
	(Row540 #= 448 and Freq540 #= 1) or
	(Row540 #= 449 and Freq540 #= 1) or
	(Row540 #= 0 and Freq540 #= 0)), 

	((Row541 #= 4 and Freq541 #= 1) or
	(Row541 #= 5 and Freq541 #= 1) or
	(Row541 #= 6 and Freq541 #= 1) or
	(Row541 #= 7 and Freq541 #= 1) or
	(Row541 #= 8 and Freq541 #= 1) or
	(Row541 #= 9 and Freq541 #= 1) or
	(Row541 #= 10 and Freq541 #= 1) or
	(Row541 #= 11 and Freq541 #= 1) or
	(Row541 #= 12 and Freq541 #= 1) or
	(Row541 #= 13 and Freq541 #= 1) or
	(Row541 #= 14 and Freq541 #= 1) or
	(Row541 #= 15 and Freq541 #= 1) or
	(Row541 #= 16 and Freq541 #= 1) or
	(Row541 #= 19 and Freq541 #= 1) or
	(Row541 #= 20 and Freq541 #= 1) or
	(Row541 #= 436 and Freq541 #= 1) or
	(Row541 #= 439 and Freq541 #= 1) or
	(Row541 #= 440 and Freq541 #= 1) or
	(Row541 #= 441 and Freq541 #= 1) or
	(Row541 #= 442 and Freq541 #= 1) or
	(Row541 #= 443 and Freq541 #= 1) or
	(Row541 #= 444 and Freq541 #= 1) or
	(Row541 #= 445 and Freq541 #= 1) or
	(Row541 #= 447 and Freq541 #= 1) or
	(Row541 #= 449 and Freq541 #= 1) or
	(Row541 #= 450 and Freq541 #= 1) or
	(Row541 #= 0 and Freq541 #= 0)), 

	((Row542 #= 4 and Freq542 #= 1) or
	(Row542 #= 5 and Freq542 #= 1) or
	(Row542 #= 6 and Freq542 #= 1) or
	(Row542 #= 7 and Freq542 #= 1) or
	(Row542 #= 8 and Freq542 #= 1) or
	(Row542 #= 9 and Freq542 #= 1) or
	(Row542 #= 10 and Freq542 #= 1) or
	(Row542 #= 11 and Freq542 #= 1) or
	(Row542 #= 12 and Freq542 #= 1) or
	(Row542 #= 13 and Freq542 #= 1) or
	(Row542 #= 14 and Freq542 #= 1) or
	(Row542 #= 15 and Freq542 #= 1) or
	(Row542 #= 17 and Freq542 #= 1) or
	(Row542 #= 18 and Freq542 #= 1) or
	(Row542 #= 19 and Freq542 #= 1) or
	(Row542 #= 20 and Freq542 #= 1) or
	(Row542 #= 439 and Freq542 #= 1) or
	(Row542 #= 440 and Freq542 #= 1) or
	(Row542 #= 441 and Freq542 #= 1) or
	(Row542 #= 442 and Freq542 #= 1) or
	(Row542 #= 443 and Freq542 #= 1) or
	(Row542 #= 444 and Freq542 #= 1) or
	(Row542 #= 445 and Freq542 #= 1) or
	(Row542 #= 447 and Freq542 #= 1) or
	(Row542 #= 448 and Freq542 #= 1) or
	(Row542 #= 449 and Freq542 #= 1) or
	(Row542 #= 450 and Freq542 #= 1) or
	(Row542 #= 0 and Freq542 #= 0)), 

	((Row543 #= 4 and Freq543 #= 1) or
	(Row543 #= 5 and Freq543 #= 2) or
	(Row543 #= 6 and Freq543 #= 1) or
	(Row543 #= 7 and Freq543 #= 1) or
	(Row543 #= 8 and Freq543 #= 1) or
	(Row543 #= 9 and Freq543 #= 1) or
	(Row543 #= 10 and Freq543 #= 1) or
	(Row543 #= 11 and Freq543 #= 2) or
	(Row543 #= 12 and Freq543 #= 1) or
	(Row543 #= 13 and Freq543 #= 1) or
	(Row543 #= 14 and Freq543 #= 1) or
	(Row543 #= 15 and Freq543 #= 1) or
	(Row543 #= 16 and Freq543 #= 1) or
	(Row543 #= 17 and Freq543 #= 1) or
	(Row543 #= 18 and Freq543 #= 1) or
	(Row543 #= 19 and Freq543 #= 1) or
	(Row543 #= 20 and Freq543 #= 1) or
	(Row543 #= 21 and Freq543 #= 1) or
	(Row543 #= 22 and Freq543 #= 1) or
	(Row543 #= 437 and Freq543 #= 1) or
	(Row543 #= 439 and Freq543 #= 1) or
	(Row543 #= 440 and Freq543 #= 1) or
	(Row543 #= 441 and Freq543 #= 1) or
	(Row543 #= 442 and Freq543 #= 1) or
	(Row543 #= 443 and Freq543 #= 1) or
	(Row543 #= 444 and Freq543 #= 1) or
	(Row543 #= 445 and Freq543 #= 1) or
	(Row543 #= 446 and Freq543 #= 2) or
	(Row543 #= 447 and Freq543 #= 1) or
	(Row543 #= 448 and Freq543 #= 1) or
	(Row543 #= 449 and Freq543 #= 1) or
	(Row543 #= 450 and Freq543 #= 1) or
	(Row543 #= 0 and Freq543 #= 0)), 

	((Row544 #= 4 and Freq544 #= 1) or
	(Row544 #= 5 and Freq544 #= 2) or
	(Row544 #= 6 and Freq544 #= 2) or
	(Row544 #= 7 and Freq544 #= 2) or
	(Row544 #= 8 and Freq544 #= 1) or
	(Row544 #= 9 and Freq544 #= 1) or
	(Row544 #= 10 and Freq544 #= 2) or
	(Row544 #= 11 and Freq544 #= 1) or
	(Row544 #= 12 and Freq544 #= 2) or
	(Row544 #= 13 and Freq544 #= 2) or
	(Row544 #= 14 and Freq544 #= 1) or
	(Row544 #= 15 and Freq544 #= 1) or
	(Row544 #= 16 and Freq544 #= 1) or
	(Row544 #= 17 and Freq544 #= 1) or
	(Row544 #= 18 and Freq544 #= 1) or
	(Row544 #= 19 and Freq544 #= 1) or
	(Row544 #= 20 and Freq544 #= 1) or
	(Row544 #= 435 and Freq544 #= 1) or
	(Row544 #= 437 and Freq544 #= 1) or
	(Row544 #= 438 and Freq544 #= 1) or
	(Row544 #= 439 and Freq544 #= 1) or
	(Row544 #= 440 and Freq544 #= 1) or
	(Row544 #= 441 and Freq544 #= 1) or
	(Row544 #= 442 and Freq544 #= 2) or
	(Row544 #= 443 and Freq544 #= 1) or
	(Row544 #= 444 and Freq544 #= 1) or
	(Row544 #= 445 and Freq544 #= 1) or
	(Row544 #= 446 and Freq544 #= 1) or
	(Row544 #= 447 and Freq544 #= 2) or
	(Row544 #= 448 and Freq544 #= 1) or
	(Row544 #= 449 and Freq544 #= 2) or
	(Row544 #= 450 and Freq544 #= 1) or
	(Row544 #= 0 and Freq544 #= 0)), 

	((Row545 #= 4 and Freq545 #= 2) or
	(Row545 #= 5 and Freq545 #= 3) or
	(Row545 #= 6 and Freq545 #= 3) or
	(Row545 #= 7 and Freq545 #= 2) or
	(Row545 #= 8 and Freq545 #= 2) or
	(Row545 #= 9 and Freq545 #= 2) or
	(Row545 #= 10 and Freq545 #= 2) or
	(Row545 #= 11 and Freq545 #= 2) or
	(Row545 #= 12 and Freq545 #= 2) or
	(Row545 #= 13 and Freq545 #= 2) or
	(Row545 #= 14 and Freq545 #= 2) or
	(Row545 #= 15 and Freq545 #= 1) or
	(Row545 #= 16 and Freq545 #= 1) or
	(Row545 #= 17 and Freq545 #= 1) or
	(Row545 #= 18 and Freq545 #= 1) or
	(Row545 #= 19 and Freq545 #= 1) or
	(Row545 #= 20 and Freq545 #= 1) or
	(Row545 #= 21 and Freq545 #= 1) or
	(Row545 #= 22 and Freq545 #= 1) or
	(Row545 #= 433 and Freq545 #= 1) or
	(Row545 #= 438 and Freq545 #= 1) or
	(Row545 #= 439 and Freq545 #= 1) or
	(Row545 #= 440 and Freq545 #= 1) or
	(Row545 #= 441 and Freq545 #= 1) or
	(Row545 #= 442 and Freq545 #= 2) or
	(Row545 #= 443 and Freq545 #= 2) or
	(Row545 #= 444 and Freq545 #= 1) or
	(Row545 #= 445 and Freq545 #= 2) or
	(Row545 #= 446 and Freq545 #= 2) or
	(Row545 #= 447 and Freq545 #= 2) or
	(Row545 #= 448 and Freq545 #= 2) or
	(Row545 #= 449 and Freq545 #= 3) or
	(Row545 #= 450 and Freq545 #= 2) or
	(Row545 #= 0 and Freq545 #= 0)), 

	((Row546 #= 4 and Freq546 #= 2) or
	(Row546 #= 5 and Freq546 #= 3) or
	(Row546 #= 6 and Freq546 #= 3) or
	(Row546 #= 7 and Freq546 #= 3) or
	(Row546 #= 8 and Freq546 #= 3) or
	(Row546 #= 9 and Freq546 #= 2) or
	(Row546 #= 10 and Freq546 #= 2) or
	(Row546 #= 11 and Freq546 #= 2) or
	(Row546 #= 12 and Freq546 #= 2) or
	(Row546 #= 13 and Freq546 #= 2) or
	(Row546 #= 14 and Freq546 #= 2) or
	(Row546 #= 15 and Freq546 #= 2) or
	(Row546 #= 16 and Freq546 #= 1) or
	(Row546 #= 17 and Freq546 #= 1) or
	(Row546 #= 18 and Freq546 #= 1) or
	(Row546 #= 19 and Freq546 #= 1) or
	(Row546 #= 20 and Freq546 #= 1) or
	(Row546 #= 21 and Freq546 #= 1) or
	(Row546 #= 22 and Freq546 #= 1) or
	(Row546 #= 23 and Freq546 #= 1) or
	(Row546 #= 437 and Freq546 #= 1) or
	(Row546 #= 438 and Freq546 #= 1) or
	(Row546 #= 439 and Freq546 #= 1) or
	(Row546 #= 440 and Freq546 #= 1) or
	(Row546 #= 441 and Freq546 #= 2) or
	(Row546 #= 442 and Freq546 #= 2) or
	(Row546 #= 443 and Freq546 #= 1) or
	(Row546 #= 444 and Freq546 #= 2) or
	(Row546 #= 445 and Freq546 #= 2) or
	(Row546 #= 446 and Freq546 #= 2) or
	(Row546 #= 447 and Freq546 #= 3) or
	(Row546 #= 448 and Freq546 #= 3) or
	(Row546 #= 449 and Freq546 #= 3) or
	(Row546 #= 450 and Freq546 #= 3) or
	(Row546 #= 0 and Freq546 #= 0)), 

	((Row547 #= 4 and Freq547 #= 3) or
	(Row547 #= 5 and Freq547 #= 3) or
	(Row547 #= 6 and Freq547 #= 4) or
	(Row547 #= 7 and Freq547 #= 4) or
	(Row547 #= 8 and Freq547 #= 3) or
	(Row547 #= 9 and Freq547 #= 3) or
	(Row547 #= 10 and Freq547 #= 2) or
	(Row547 #= 11 and Freq547 #= 2) or
	(Row547 #= 12 and Freq547 #= 3) or
	(Row547 #= 13 and Freq547 #= 3) or
	(Row547 #= 14 and Freq547 #= 2) or
	(Row547 #= 15 and Freq547 #= 2) or
	(Row547 #= 16 and Freq547 #= 1) or
	(Row547 #= 17 and Freq547 #= 1) or
	(Row547 #= 18 and Freq547 #= 1) or
	(Row547 #= 19 and Freq547 #= 1) or
	(Row547 #= 20 and Freq547 #= 1) or
	(Row547 #= 21 and Freq547 #= 1) or
	(Row547 #= 22 and Freq547 #= 1) or
	(Row547 #= 435 and Freq547 #= 1) or
	(Row547 #= 436 and Freq547 #= 1) or
	(Row547 #= 437 and Freq547 #= 1) or
	(Row547 #= 438 and Freq547 #= 1) or
	(Row547 #= 439 and Freq547 #= 1) or
	(Row547 #= 440 and Freq547 #= 1) or
	(Row547 #= 441 and Freq547 #= 2) or
	(Row547 #= 442 and Freq547 #= 2) or
	(Row547 #= 443 and Freq547 #= 2) or
	(Row547 #= 444 and Freq547 #= 2) or
	(Row547 #= 445 and Freq547 #= 2) or
	(Row547 #= 446 and Freq547 #= 2) or
	(Row547 #= 447 and Freq547 #= 3) or
	(Row547 #= 448 and Freq547 #= 3) or
	(Row547 #= 449 and Freq547 #= 4) or
	(Row547 #= 450 and Freq547 #= 4) or
	(Row547 #= 0 and Freq547 #= 0)), 

	((Row548 #= 4 and Freq548 #= 4) or
	(Row548 #= 5 and Freq548 #= 4) or
	(Row548 #= 6 and Freq548 #= 5) or
	(Row548 #= 7 and Freq548 #= 4) or
	(Row548 #= 8 and Freq548 #= 4) or
	(Row548 #= 9 and Freq548 #= 3) or
	(Row548 #= 10 and Freq548 #= 3) or
	(Row548 #= 11 and Freq548 #= 3) or
	(Row548 #= 12 and Freq548 #= 3) or
	(Row548 #= 13 and Freq548 #= 3) or
	(Row548 #= 14 and Freq548 #= 2) or
	(Row548 #= 15 and Freq548 #= 2) or
	(Row548 #= 16 and Freq548 #= 1) or
	(Row548 #= 17 and Freq548 #= 2) or
	(Row548 #= 18 and Freq548 #= 1) or
	(Row548 #= 19 and Freq548 #= 1) or
	(Row548 #= 20 and Freq548 #= 1) or
	(Row548 #= 21 and Freq548 #= 1) or
	(Row548 #= 22 and Freq548 #= 1) or
	(Row548 #= 434 and Freq548 #= 1) or
	(Row548 #= 435 and Freq548 #= 1) or
	(Row548 #= 436 and Freq548 #= 1) or
	(Row548 #= 437 and Freq548 #= 1) or
	(Row548 #= 438 and Freq548 #= 1) or
	(Row548 #= 439 and Freq548 #= 1) or
	(Row548 #= 440 and Freq548 #= 1) or
	(Row548 #= 441 and Freq548 #= 2) or
	(Row548 #= 442 and Freq548 #= 2) or
	(Row548 #= 443 and Freq548 #= 3) or
	(Row548 #= 444 and Freq548 #= 3) or
	(Row548 #= 445 and Freq548 #= 3) or
	(Row548 #= 446 and Freq548 #= 3) or
	(Row548 #= 447 and Freq548 #= 4) or
	(Row548 #= 448 and Freq548 #= 4) or
	(Row548 #= 449 and Freq548 #= 5) or
	(Row548 #= 450 and Freq548 #= 4) or
	(Row548 #= 0 and Freq548 #= 0)), 

	((Row549 #= 4 and Freq549 #= 4) or
	(Row549 #= 5 and Freq549 #= 4) or
	(Row549 #= 6 and Freq549 #= 4) or
	(Row549 #= 7 and Freq549 #= 4) or
	(Row549 #= 8 and Freq549 #= 4) or
	(Row549 #= 9 and Freq549 #= 4) or
	(Row549 #= 10 and Freq549 #= 2) or
	(Row549 #= 11 and Freq549 #= 2) or
	(Row549 #= 12 and Freq549 #= 3) or
	(Row549 #= 13 and Freq549 #= 3) or
	(Row549 #= 14 and Freq549 #= 3) or
	(Row549 #= 15 and Freq549 #= 2) or
	(Row549 #= 16 and Freq549 #= 2) or
	(Row549 #= 17 and Freq549 #= 1) or
	(Row549 #= 18 and Freq549 #= 1) or
	(Row549 #= 19 and Freq549 #= 1) or
	(Row549 #= 20 and Freq549 #= 1) or
	(Row549 #= 21 and Freq549 #= 1) or
	(Row549 #= 22 and Freq549 #= 1) or
	(Row549 #= 23 and Freq549 #= 1) or
	(Row549 #= 434 and Freq549 #= 1) or
	(Row549 #= 435 and Freq549 #= 1) or
	(Row549 #= 436 and Freq549 #= 1) or
	(Row549 #= 437 and Freq549 #= 1) or
	(Row549 #= 438 and Freq549 #= 1) or
	(Row549 #= 439 and Freq549 #= 1) or
	(Row549 #= 440 and Freq549 #= 2) or
	(Row549 #= 441 and Freq549 #= 2) or
	(Row549 #= 442 and Freq549 #= 2) or
	(Row549 #= 443 and Freq549 #= 2) or
	(Row549 #= 444 and Freq549 #= 2) or
	(Row549 #= 445 and Freq549 #= 3) or
	(Row549 #= 446 and Freq549 #= 4) or
	(Row549 #= 447 and Freq549 #= 4) or
	(Row549 #= 448 and Freq549 #= 3) or
	(Row549 #= 449 and Freq549 #= 5) or
	(Row549 #= 450 and Freq549 #= 5) or
	(Row549 #= 0 and Freq549 #= 0)), 

	((Row550 #= 4 and Freq550 #= 6) or
	(Row550 #= 5 and Freq550 #= 7) or
	(Row550 #= 6 and Freq550 #= 7) or
	(Row550 #= 7 and Freq550 #= 7) or
	(Row550 #= 8 and Freq550 #= 6) or
	(Row550 #= 9 and Freq550 #= 5) or
	(Row550 #= 10 and Freq550 #= 4) or
	(Row550 #= 11 and Freq550 #= 4) or
	(Row550 #= 12 and Freq550 #= 3) or
	(Row550 #= 13 and Freq550 #= 3) or
	(Row550 #= 14 and Freq550 #= 3) or
	(Row550 #= 15 and Freq550 #= 2) or
	(Row550 #= 16 and Freq550 #= 2) or
	(Row550 #= 17 and Freq550 #= 2) or
	(Row550 #= 18 and Freq550 #= 1) or
	(Row550 #= 19 and Freq550 #= 1) or
	(Row550 #= 20 and Freq550 #= 1) or
	(Row550 #= 21 and Freq550 #= 1) or
	(Row550 #= 433 and Freq550 #= 1) or
	(Row550 #= 434 and Freq550 #= 1) or
	(Row550 #= 436 and Freq550 #= 1) or
	(Row550 #= 437 and Freq550 #= 1) or
	(Row550 #= 438 and Freq550 #= 1) or
	(Row550 #= 439 and Freq550 #= 1) or
	(Row550 #= 440 and Freq550 #= 2) or
	(Row550 #= 441 and Freq550 #= 2) or
	(Row550 #= 442 and Freq550 #= 4) or
	(Row550 #= 443 and Freq550 #= 4) or
	(Row550 #= 444 and Freq550 #= 3) or
	(Row550 #= 445 and Freq550 #= 3) or
	(Row550 #= 446 and Freq550 #= 5) or
	(Row550 #= 447 and Freq550 #= 5) or
	(Row550 #= 448 and Freq550 #= 7) or
	(Row550 #= 449 and Freq550 #= 6) or
	(Row550 #= 450 and Freq550 #= 7) or
	(Row550 #= 0 and Freq550 #= 0)), 

	((Row551 #= 4 and Freq551 #= 9) or
	(Row551 #= 5 and Freq551 #= 11) or
	(Row551 #= 6 and Freq551 #= 8) or
	(Row551 #= 7 and Freq551 #= 10) or
	(Row551 #= 8 and Freq551 #= 8) or
	(Row551 #= 9 and Freq551 #= 7) or
	(Row551 #= 10 and Freq551 #= 5) or
	(Row551 #= 11 and Freq551 #= 3) or
	(Row551 #= 12 and Freq551 #= 4) or
	(Row551 #= 13 and Freq551 #= 4) or
	(Row551 #= 14 and Freq551 #= 3) or
	(Row551 #= 15 and Freq551 #= 2) or
	(Row551 #= 16 and Freq551 #= 2) or
	(Row551 #= 17 and Freq551 #= 1) or
	(Row551 #= 18 and Freq551 #= 1) or
	(Row551 #= 19 and Freq551 #= 1) or
	(Row551 #= 20 and Freq551 #= 1) or
	(Row551 #= 21 and Freq551 #= 1) or
	(Row551 #= 22 and Freq551 #= 1) or
	(Row551 #= 437 and Freq551 #= 1) or
	(Row551 #= 438 and Freq551 #= 1) or
	(Row551 #= 439 and Freq551 #= 1) or
	(Row551 #= 440 and Freq551 #= 1) or
	(Row551 #= 441 and Freq551 #= 2) or
	(Row551 #= 442 and Freq551 #= 3) or
	(Row551 #= 443 and Freq551 #= 4) or
	(Row551 #= 444 and Freq551 #= 3) or
	(Row551 #= 445 and Freq551 #= 3) or
	(Row551 #= 446 and Freq551 #= 8) or
	(Row551 #= 447 and Freq551 #= 9) or
	(Row551 #= 448 and Freq551 #= 9) or
	(Row551 #= 449 and Freq551 #= 9) or
	(Row551 #= 450 and Freq551 #= 10) or
	(Row551 #= 0 and Freq551 #= 0)), 

	((Row552 #= 4 and Freq552 #= 11) or
	(Row552 #= 5 and Freq552 #= 11) or
	(Row552 #= 6 and Freq552 #= 8) or
	(Row552 #= 7 and Freq552 #= 10) or
	(Row552 #= 8 and Freq552 #= 8) or
	(Row552 #= 9 and Freq552 #= 7) or
	(Row552 #= 10 and Freq552 #= 4) or
	(Row552 #= 11 and Freq552 #= 4) or
	(Row552 #= 12 and Freq552 #= 4) or
	(Row552 #= 13 and Freq552 #= 4) or
	(Row552 #= 14 and Freq552 #= 2) or
	(Row552 #= 15 and Freq552 #= 2) or
	(Row552 #= 16 and Freq552 #= 2) or
	(Row552 #= 17 and Freq552 #= 1) or
	(Row552 #= 18 and Freq552 #= 1) or
	(Row552 #= 19 and Freq552 #= 1) or
	(Row552 #= 20 and Freq552 #= 1) or
	(Row552 #= 21 and Freq552 #= 1) or
	(Row552 #= 22 and Freq552 #= 1) or
	(Row552 #= 435 and Freq552 #= 1) or
	(Row552 #= 436 and Freq552 #= 1) or
	(Row552 #= 437 and Freq552 #= 1) or
	(Row552 #= 438 and Freq552 #= 1) or
	(Row552 #= 439 and Freq552 #= 1) or
	(Row552 #= 440 and Freq552 #= 1) or
	(Row552 #= 441 and Freq552 #= 3) or
	(Row552 #= 442 and Freq552 #= 3) or
	(Row552 #= 443 and Freq552 #= 4) or
	(Row552 #= 444 and Freq552 #= 4) or
	(Row552 #= 445 and Freq552 #= 4) or
	(Row552 #= 446 and Freq552 #= 8) or
	(Row552 #= 447 and Freq552 #= 8) or
	(Row552 #= 448 and Freq552 #= 9) or
	(Row552 #= 449 and Freq552 #= 9) or
	(Row552 #= 450 and Freq552 #= 10) or
	(Row552 #= 0 and Freq552 #= 0)), 

	((Row553 #= 4 and Freq553 #= 12) or
	(Row553 #= 5 and Freq553 #= 12) or
	(Row553 #= 6 and Freq553 #= 9) or
	(Row553 #= 7 and Freq553 #= 9) or
	(Row553 #= 8 and Freq553 #= 8) or
	(Row553 #= 9 and Freq553 #= 6) or
	(Row553 #= 10 and Freq553 #= 3) or
	(Row553 #= 11 and Freq553 #= 3) or
	(Row553 #= 12 and Freq553 #= 4) or
	(Row553 #= 13 and Freq553 #= 3) or
	(Row553 #= 14 and Freq553 #= 2) or
	(Row553 #= 15 and Freq553 #= 2) or
	(Row553 #= 16 and Freq553 #= 2) or
	(Row553 #= 17 and Freq553 #= 1) or
	(Row553 #= 21 and Freq553 #= 1) or
	(Row553 #= 214 and Freq553 #= 1) or
	(Row553 #= 436 and Freq553 #= 1) or
	(Row553 #= 437 and Freq553 #= 1) or
	(Row553 #= 438 and Freq553 #= 1) or
	(Row553 #= 439 and Freq553 #= 1) or
	(Row553 #= 440 and Freq553 #= 1) or
	(Row553 #= 441 and Freq553 #= 2) or
	(Row553 #= 442 and Freq553 #= 3) or
	(Row553 #= 443 and Freq553 #= 4) or
	(Row553 #= 444 and Freq553 #= 3) or
	(Row553 #= 445 and Freq553 #= 4) or
	(Row553 #= 446 and Freq553 #= 8) or
	(Row553 #= 447 and Freq553 #= 9) or
	(Row553 #= 448 and Freq553 #= 10) or
	(Row553 #= 449 and Freq553 #= 10) or
	(Row553 #= 450 and Freq553 #= 11) or
	(Row553 #= 0 and Freq553 #= 0)), 

	((Row554 #= 4 and Freq554 #= 10) or
	(Row554 #= 5 and Freq554 #= 10) or
	(Row554 #= 6 and Freq554 #= 7) or
	(Row554 #= 7 and Freq554 #= 8) or
	(Row554 #= 8 and Freq554 #= 6) or
	(Row554 #= 9 and Freq554 #= 5) or
	(Row554 #= 10 and Freq554 #= 2) or
	(Row554 #= 11 and Freq554 #= 2) or
	(Row554 #= 12 and Freq554 #= 3) or
	(Row554 #= 13 and Freq554 #= 3) or
	(Row554 #= 14 and Freq554 #= 2) or
	(Row554 #= 15 and Freq554 #= 1) or
	(Row554 #= 16 and Freq554 #= 1) or
	(Row554 #= 17 and Freq554 #= 1) or
	(Row554 #= 20 and Freq554 #= 1) or
	(Row554 #= 213 and Freq554 #= 1) or
	(Row554 #= 214 and Freq554 #= 1) or
	(Row554 #= 438 and Freq554 #= 1) or
	(Row554 #= 439 and Freq554 #= 1) or
	(Row554 #= 440 and Freq554 #= 1) or
	(Row554 #= 441 and Freq554 #= 2) or
	(Row554 #= 442 and Freq554 #= 3) or
	(Row554 #= 443 and Freq554 #= 4) or
	(Row554 #= 444 and Freq554 #= 3) or
	(Row554 #= 445 and Freq554 #= 3) or
	(Row554 #= 446 and Freq554 #= 7) or
	(Row554 #= 447 and Freq554 #= 8) or
	(Row554 #= 448 and Freq554 #= 8) or
	(Row554 #= 449 and Freq554 #= 9) or
	(Row554 #= 450 and Freq554 #= 11) or
	(Row554 #= 0 and Freq554 #= 0)), 

	((Row555 #= 4 and Freq555 #= 7) or
	(Row555 #= 5 and Freq555 #= 8) or
	(Row555 #= 6 and Freq555 #= 6) or
	(Row555 #= 7 and Freq555 #= 6) or
	(Row555 #= 8 and Freq555 #= 5) or
	(Row555 #= 9 and Freq555 #= 4) or
	(Row555 #= 10 and Freq555 #= 3) or
	(Row555 #= 11 and Freq555 #= 3) or
	(Row555 #= 12 and Freq555 #= 2) or
	(Row555 #= 13 and Freq555 #= 2) or
	(Row555 #= 14 and Freq555 #= 1) or
	(Row555 #= 15 and Freq555 #= 2) or
	(Row555 #= 16 and Freq555 #= 1) or
	(Row555 #= 17 and Freq555 #= 1) or
	(Row555 #= 20 and Freq555 #= 1) or
	(Row555 #= 213 and Freq555 #= 1) or
	(Row555 #= 214 and Freq555 #= 2) or
	(Row555 #= 215 and Freq555 #= 1) or
	(Row555 #= 216 and Freq555 #= 1) or
	(Row555 #= 436 and Freq555 #= 1) or
	(Row555 #= 437 and Freq555 #= 1) or
	(Row555 #= 438 and Freq555 #= 1) or
	(Row555 #= 439 and Freq555 #= 1) or
	(Row555 #= 440 and Freq555 #= 1) or
	(Row555 #= 441 and Freq555 #= 2) or
	(Row555 #= 442 and Freq555 #= 4) or
	(Row555 #= 443 and Freq555 #= 3) or
	(Row555 #= 444 and Freq555 #= 3) or
	(Row555 #= 445 and Freq555 #= 4) or
	(Row555 #= 446 and Freq555 #= 7) or
	(Row555 #= 447 and Freq555 #= 7) or
	(Row555 #= 448 and Freq555 #= 8) or
	(Row555 #= 449 and Freq555 #= 8) or
	(Row555 #= 450 and Freq555 #= 8) or
	(Row555 #= 0 and Freq555 #= 0)), 

	((Row556 #= 4 and Freq556 #= 9) or
	(Row556 #= 5 and Freq556 #= 8) or
	(Row556 #= 6 and Freq556 #= 5) or
	(Row556 #= 7 and Freq556 #= 6) or
	(Row556 #= 8 and Freq556 #= 5) or
	(Row556 #= 9 and Freq556 #= 3) or
	(Row556 #= 10 and Freq556 #= 2) or
	(Row556 #= 11 and Freq556 #= 2) or
	(Row556 #= 12 and Freq556 #= 2) or
	(Row556 #= 13 and Freq556 #= 1) or
	(Row556 #= 14 and Freq556 #= 1) or
	(Row556 #= 15 and Freq556 #= 2) or
	(Row556 #= 16 and Freq556 #= 1) or
	(Row556 #= 17 and Freq556 #= 1) or
	(Row556 #= 21 and Freq556 #= 1) or
	(Row556 #= 211 and Freq556 #= 1) or
	(Row556 #= 212 and Freq556 #= 1) or
	(Row556 #= 213 and Freq556 #= 1) or
	(Row556 #= 214 and Freq556 #= 3) or
	(Row556 #= 215 and Freq556 #= 1) or
	(Row556 #= 216 and Freq556 #= 1) or
	(Row556 #= 436 and Freq556 #= 1) or
	(Row556 #= 437 and Freq556 #= 1) or
	(Row556 #= 438 and Freq556 #= 1) or
	(Row556 #= 439 and Freq556 #= 1) or
	(Row556 #= 440 and Freq556 #= 2) or
	(Row556 #= 441 and Freq556 #= 2) or
	(Row556 #= 442 and Freq556 #= 2) or
	(Row556 #= 443 and Freq556 #= 3) or
	(Row556 #= 444 and Freq556 #= 4) or
	(Row556 #= 445 and Freq556 #= 4) or
	(Row556 #= 446 and Freq556 #= 6) or
	(Row556 #= 447 and Freq556 #= 7) or
	(Row556 #= 448 and Freq556 #= 7) or
	(Row556 #= 449 and Freq556 #= 8) or
	(Row556 #= 450 and Freq556 #= 9) or
	(Row556 #= 0 and Freq556 #= 0)), 

	((Row557 #= 4 and Freq557 #= 8) or
	(Row557 #= 5 and Freq557 #= 6) or
	(Row557 #= 6 and Freq557 #= 6) or
	(Row557 #= 7 and Freq557 #= 5) or
	(Row557 #= 8 and Freq557 #= 4) or
	(Row557 #= 9 and Freq557 #= 3) or
	(Row557 #= 10 and Freq557 #= 2) or
	(Row557 #= 11 and Freq557 #= 2) or
	(Row557 #= 12 and Freq557 #= 2) or
	(Row557 #= 13 and Freq557 #= 2) or
	(Row557 #= 14 and Freq557 #= 2) or
	(Row557 #= 15 and Freq557 #= 2) or
	(Row557 #= 16 and Freq557 #= 1) or
	(Row557 #= 18 and Freq557 #= 1) or
	(Row557 #= 19 and Freq557 #= 1) or
	(Row557 #= 23 and Freq557 #= 1) or
	(Row557 #= 166 and Freq557 #= 1) or
	(Row557 #= 210 and Freq557 #= 1) or
	(Row557 #= 211 and Freq557 #= 1) or
	(Row557 #= 212 and Freq557 #= 1) or
	(Row557 #= 213 and Freq557 #= 1) or
	(Row557 #= 214 and Freq557 #= 3) or
	(Row557 #= 215 and Freq557 #= 1) or
	(Row557 #= 216 and Freq557 #= 1) or
	(Row557 #= 434 and Freq557 #= 1) or
	(Row557 #= 436 and Freq557 #= 1) or
	(Row557 #= 437 and Freq557 #= 2) or
	(Row557 #= 438 and Freq557 #= 1) or
	(Row557 #= 439 and Freq557 #= 1) or
	(Row557 #= 440 and Freq557 #= 2) or
	(Row557 #= 441 and Freq557 #= 3) or
	(Row557 #= 442 and Freq557 #= 2) or
	(Row557 #= 443 and Freq557 #= 3) or
	(Row557 #= 444 and Freq557 #= 4) or
	(Row557 #= 445 and Freq557 #= 5) or
	(Row557 #= 446 and Freq557 #= 5) or
	(Row557 #= 447 and Freq557 #= 7) or
	(Row557 #= 448 and Freq557 #= 7) or
	(Row557 #= 449 and Freq557 #= 8) or
	(Row557 #= 450 and Freq557 #= 8) or
	(Row557 #= 0 and Freq557 #= 0)), 

	((Row558 #= 0 and Freq558 #= 0)), 

	% All of the values assumed by the Row<i> variables must be  
	% all different or zero; multiple zeros are allowed  
	alldifferent_except(Rows), 

	% Optimize: the sum of the selected interaction 
	% frequencies is maximal (it is the minimum of 
	% the additive inverse of the sum). The predicate 
	% enforce_symmetry/1 ensures each genomic bin is 
	% involved in only one interaction in the solution set.
	Cost #= -sum(Freqs), 
	minimize(enforce_symmetry(RowFile, FreqFile, Interactions, Rows, Freqs), Cost).