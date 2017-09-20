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
	Rows = [Row1, Row2, Row3, Row4, Row5, Row6, Row7, Row8, Row9, Row10, Row11, Row12, Row13, Row14, Row15, Row16, Row17, Row18, Row19, Row20, Row21, Row22, Row23, Row24, Row25, Row26, Row27, Row28, Row29, Row30, Row31, Row32, Row33, Row34, Row35, Row36, Row37, Row38, Row39, Row40, Row41, Row42, Row43, Row44, Row45, Row46, Row47, Row48, Row49, Row50, Row51, Row52, Row53, Row54, Row55, Row56, Row57, Row58, Row59, Row60, Row61, Row62, Row63, Row64, Row65, Row66, Row67, Row68, Row69, Row70, Row71, Row72, Row73, Row74, Row75, Row76, Row77, Row78, Row79, Row80, Row81, Row82, Row83, Row84, Row85, Row86, Row87, Row88, Row89, Row90, Row91, Row92, Row93, Row94, Row95, Row96, Row97, Row98, Row99, Row100, Row101, Row102, Row103, Row104, Row105, Row106, Row107, Row108, Row109, Row110, Row111, Row112, Row113, Row114, Row115, Row116, Row117, Row118, Row119, Row120, Row121, Row122, Row123, Row124, Row125, Row126, Row127, Row128, Row129, Row130, Row131, Row132, Row133, Row134, Row135, Row136, Row137, Row138, Row139, Row140, Row141, Row142, Row143, Row144, Row145, Row146, Row147, Row148, Row149, Row150, Row151, Row152, Row153, Row154, Row155, Row156, Row157, Row158, Row159, Row160, Row161, Row162, Row163, Row164, Row165, Row166, Row167, Row168, Row169, Row170, Row171, Row172, Row173, Row174, Row175, Row176, Row177, Row178, Row179, Row180, Row181, Row182, Row183, Row184, Row185, Row186, Row187, Row188, Row189, Row190, Row191, Row192, Row193, Row194, Row195, Row196, Row197, Row198, Row199, Row200, Row201, Row202, Row203, Row204, Row205, Row206, Row207, Row208, Row209, Row210, Row211, Row212, Row213, Row214, Row215, Row216, Row217, Row218, Row219, Row220, Row221, Row222, Row223, Row224, Row225, Row226, Row227, Row228, Row229, Row230, Row231, Row232, Row233, Row234, Row235, Row236, Row237, Row238, Row239, Row240, Row241, Row242, Row243, Row244, Row245, Row246, Row247, Row248, Row249, Row250, Row251, Row252, Row253, Row254, Row255, Row256, Row257, Row258, Row259, Row260, Row261, Row262, Row263, Row264, Row265, Row266, Row267, Row268, Row269, Row270, Row271, Row272, Row273, Row274, Row275, Row276, Row277, Row278, Row279, Row280, Row281, Row282, Row283, Row284, Row285, Row286, Row287, Row288, Row289, Row290, Row291, Row292, Row293, Row294, Row295, Row296, Row297, Row298, Row299, Row300, Row301, Row302, Row303, Row304, Row305, Row306, Row307, Row308, Row309, Row310, Row311, Row312, Row313, Row314, Row315, Row316, Row317, Row318, Row319, Row320, Row321, Row322, Row323, Row324, Row325, Row326, Row327, Row328, Row329, Row330, Row331, Row332, Row333, Row334, Row335, Row336, Row337, Row338, Row339, Row340, Row341, Row342, Row343, Row344, Row345, Row346, Row347, Row348, Row349, Row350, Row351, Row352, Row353, Row354, Row355, Row356, Row357, Row358, Row359, Row360, Row361, Row362, Row363, Row364, Row365, Row366, Row367, Row368, Row369, Row370, Row371, Row372, Row373, Row374, Row375, Row376, Row377, Row378, Row379, Row380, Row381, Row382, Row383, Row384, Row385, Row386, Row387, Row388, Row389, Row390, Row391, Row392, Row393, Row394, Row395, Row396, Row397, Row398, Row399, Row400, Row401, Row402, Row403, Row404, Row405, Row406, Row407, Row408, Row409, Row410, Row411, Row412, Row413, Row414, Row415, Row416, Row417, Row418, Row419, Row420, Row421, Row422, Row423, Row424, Row425, Row426, Row427, Row428, Row429, Row430, Row431, Row432, Row433, Row434, Row435, Row436, Row437, Row438, Row439, Row440, Row441, Row442, Row443, Row444, Row445, Row446, Row447, Row448, Row449, Row450, Row451, Row452, Row453, Row454, Row455, Row456, Row457, Row458, Row459, Row460, Row461, Row462, Row463, Row464, Row465, Row466, Row467, Row468, Row469, Row470, Row471, Row472, Row473, Row474, Row475, Row476, Row477, Row478, Row479, Row480, Row481, Row482, Row483, Row484, Row485, Row486, Row487, Row488, Row489, Row490, Row491, Row492, Row493, Row494, Row495, Row496, Row497, Row498, Row499, Row500, Row501, Row502, Row503, Row504, Row505, Row506, Row507, Row508, Row509, Row510, Row511, Row512, Row513, Row514, Row515, Row516, Row517, Row518, Row519, Row520, Row521, Row522, Row523, Row524, Row525, Row526, Row527, Row528, Row529, Row530, Row531, Row532, Row533, Row534, Row535, Row536, Row537, Row538, Row539, Row540, Row541, Row542, Row543, Row544, Row545, Row546, Row547, Row548, Row549, Row550, Row551, Row552, Row553, Row554, Row555, Row556, Row557, Row558],

	% The list Freqs has one variable for each row of the 
	% whole-genome contact map	
	Freqs = [Freq1, Freq2, Freq3, Freq4, Freq5, Freq6, Freq7, Freq8, Freq9, Freq10, Freq11, Freq12, Freq13, Freq14, Freq15, Freq16, Freq17, Freq18, Freq19, Freq20, Freq21, Freq22, Freq23, Freq24, Freq25, Freq26, Freq27, Freq28, Freq29, Freq30, Freq31, Freq32, Freq33, Freq34, Freq35, Freq36, Freq37, Freq38, Freq39, Freq40, Freq41, Freq42, Freq43, Freq44, Freq45, Freq46, Freq47, Freq48, Freq49, Freq50, Freq51, Freq52, Freq53, Freq54, Freq55, Freq56, Freq57, Freq58, Freq59, Freq60, Freq61, Freq62, Freq63, Freq64, Freq65, Freq66, Freq67, Freq68, Freq69, Freq70, Freq71, Freq72, Freq73, Freq74, Freq75, Freq76, Freq77, Freq78, Freq79, Freq80, Freq81, Freq82, Freq83, Freq84, Freq85, Freq86, Freq87, Freq88, Freq89, Freq90, Freq91, Freq92, Freq93, Freq94, Freq95, Freq96, Freq97, Freq98, Freq99, Freq100, Freq101, Freq102, Freq103, Freq104, Freq105, Freq106, Freq107, Freq108, Freq109, Freq110, Freq111, Freq112, Freq113, Freq114, Freq115, Freq116, Freq117, Freq118, Freq119, Freq120, Freq121, Freq122, Freq123, Freq124, Freq125, Freq126, Freq127, Freq128, Freq129, Freq130, Freq131, Freq132, Freq133, Freq134, Freq135, Freq136, Freq137, Freq138, Freq139, Freq140, Freq141, Freq142, Freq143, Freq144, Freq145, Freq146, Freq147, Freq148, Freq149, Freq150, Freq151, Freq152, Freq153, Freq154, Freq155, Freq156, Freq157, Freq158, Freq159, Freq160, Freq161, Freq162, Freq163, Freq164, Freq165, Freq166, Freq167, Freq168, Freq169, Freq170, Freq171, Freq172, Freq173, Freq174, Freq175, Freq176, Freq177, Freq178, Freq179, Freq180, Freq181, Freq182, Freq183, Freq184, Freq185, Freq186, Freq187, Freq188, Freq189, Freq190, Freq191, Freq192, Freq193, Freq194, Freq195, Freq196, Freq197, Freq198, Freq199, Freq200, Freq201, Freq202, Freq203, Freq204, Freq205, Freq206, Freq207, Freq208, Freq209, Freq210, Freq211, Freq212, Freq213, Freq214, Freq215, Freq216, Freq217, Freq218, Freq219, Freq220, Freq221, Freq222, Freq223, Freq224, Freq225, Freq226, Freq227, Freq228, Freq229, Freq230, Freq231, Freq232, Freq233, Freq234, Freq235, Freq236, Freq237, Freq238, Freq239, Freq240, Freq241, Freq242, Freq243, Freq244, Freq245, Freq246, Freq247, Freq248, Freq249, Freq250, Freq251, Freq252, Freq253, Freq254, Freq255, Freq256, Freq257, Freq258, Freq259, Freq260, Freq261, Freq262, Freq263, Freq264, Freq265, Freq266, Freq267, Freq268, Freq269, Freq270, Freq271, Freq272, Freq273, Freq274, Freq275, Freq276, Freq277, Freq278, Freq279, Freq280, Freq281, Freq282, Freq283, Freq284, Freq285, Freq286, Freq287, Freq288, Freq289, Freq290, Freq291, Freq292, Freq293, Freq294, Freq295, Freq296, Freq297, Freq298, Freq299, Freq300, Freq301, Freq302, Freq303, Freq304, Freq305, Freq306, Freq307, Freq308, Freq309, Freq310, Freq311, Freq312, Freq313, Freq314, Freq315, Freq316, Freq317, Freq318, Freq319, Freq320, Freq321, Freq322, Freq323, Freq324, Freq325, Freq326, Freq327, Freq328, Freq329, Freq330, Freq331, Freq332, Freq333, Freq334, Freq335, Freq336, Freq337, Freq338, Freq339, Freq340, Freq341, Freq342, Freq343, Freq344, Freq345, Freq346, Freq347, Freq348, Freq349, Freq350, Freq351, Freq352, Freq353, Freq354, Freq355, Freq356, Freq357, Freq358, Freq359, Freq360, Freq361, Freq362, Freq363, Freq364, Freq365, Freq366, Freq367, Freq368, Freq369, Freq370, Freq371, Freq372, Freq373, Freq374, Freq375, Freq376, Freq377, Freq378, Freq379, Freq380, Freq381, Freq382, Freq383, Freq384, Freq385, Freq386, Freq387, Freq388, Freq389, Freq390, Freq391, Freq392, Freq393, Freq394, Freq395, Freq396, Freq397, Freq398, Freq399, Freq400, Freq401, Freq402, Freq403, Freq404, Freq405, Freq406, Freq407, Freq408, Freq409, Freq410, Freq411, Freq412, Freq413, Freq414, Freq415, Freq416, Freq417, Freq418, Freq419, Freq420, Freq421, Freq422, Freq423, Freq424, Freq425, Freq426, Freq427, Freq428, Freq429, Freq430, Freq431, Freq432, Freq433, Freq434, Freq435, Freq436, Freq437, Freq438, Freq439, Freq440, Freq441, Freq442, Freq443, Freq444, Freq445, Freq446, Freq447, Freq448, Freq449, Freq450, Freq451, Freq452, Freq453, Freq454, Freq455, Freq456, Freq457, Freq458, Freq459, Freq460, Freq461, Freq462, Freq463, Freq464, Freq465, Freq466, Freq467, Freq468, Freq469, Freq470, Freq471, Freq472, Freq473, Freq474, Freq475, Freq476, Freq477, Freq478, Freq479, Freq480, Freq481, Freq482, Freq483, Freq484, Freq485, Freq486, Freq487, Freq488, Freq489, Freq490, Freq491, Freq492, Freq493, Freq494, Freq495, Freq496, Freq497, Freq498, Freq499, Freq500, Freq501, Freq502, Freq503, Freq504, Freq505, Freq506, Freq507, Freq508, Freq509, Freq510, Freq511, Freq512, Freq513, Freq514, Freq515, Freq516, Freq517, Freq518, Freq519, Freq520, Freq521, Freq522, Freq523, Freq524, Freq525, Freq526, Freq527, Freq528, Freq529, Freq530, Freq531, Freq532, Freq533, Freq534, Freq535, Freq536, Freq537, Freq538, Freq539, Freq540, Freq541, Freq542, Freq543, Freq544, Freq545, Freq546, Freq547, Freq548, Freq549, Freq550, Freq551, Freq552, Freq553, Freq554, Freq555, Freq556, Freq557, Freq558],

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
	Row6 :: [0, 245],
	Row7 :: [0, 245],
	Row8 :: [0, 245, 246],
	Row9 :: [0],
	Row10 :: [0, 246],
	Row11 :: [0, 245],
	Row12 :: [0, 245, 246],
	Row13 :: [0, 245],
	Row14 :: [0],
	Row15 :: [0],
	Row16 :: [0],
	Row17 :: [0],
	Row18 :: [0],
	Row19 :: [0],
	Row20 :: [0, 245, 246],
	Row21 :: [0, 245],
	Row22 :: [0, 245],
	Row23 :: [0, 245],
	Row24 :: [0, 245],
	Row25 :: [0, 245],
	Row26 :: [0],
	Row27 :: [0],
	Row28 :: [0],
	Row29 :: [0],
	Row30 :: [0, 246],
	Row31 :: [0],
	Row32 :: [0],
	Row33 :: [0, 246],
	Row34 :: [0, 245],
	Row35 :: [0],
	Row36 :: [0, 246],
	Row37 :: [0],
	Row38 :: [0],
	Row39 :: [0, 245],
	Row40 :: [0],
	Row41 :: [0],
	Row42 :: [0],
	Row43 :: [0],
	Row44 :: [0],
	Row45 :: [0],
	Row46 :: [0],
	Row47 :: [0, 245],
	Row48 :: [0, 245],
	Row49 :: [0, 245],
	Row50 :: [0],
	Row51 :: [0, 245],
	Row52 :: [0],
	Row53 :: [0, 245],
	Row54 :: [0],
	Row55 :: [0],
	Row56 :: [0, 246],
	Row57 :: [0],
	Row58 :: [0, 245, 246],
	Row59 :: [0, 245],
	Row60 :: [0],
	Row61 :: [0],
	Row62 :: [0],
	Row63 :: [0],
	Row64 :: [0],
	Row65 :: [0, 245],
	Row66 :: [0],
	Row67 :: [0],
	Row68 :: [0],
	Row69 :: [0, 245],
	Row70 :: [0],
	Row71 :: [0, 245],
	Row72 :: [0],
	Row73 :: [0, 245],
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
	Row88 :: [0, 245],
	Row89 :: [0],
	Row90 :: [0],
	Row91 :: [0],
	Row92 :: [0],
	Row93 :: [0],
	Row94 :: [0, 245],
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
	Row165 :: [0, 82],
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
	Row176 :: [0, 245],
	Row177 :: [0],
	Row178 :: [0],
	Row179 :: [0],
	Row180 :: [0],
	Row181 :: [0],
	Row182 :: [0],
	Row183 :: [0],
	Row184 :: [0],
	Row185 :: [0],
	Row186 :: [0],
	Row187 :: [0],
	Row188 :: [0],
	Row189 :: [0],
	Row190 :: [0],
	Row191 :: [0],
	Row192 :: [0],
	Row193 :: [0],
	Row194 :: [0],
	Row195 :: [0],
	Row196 :: [0, 99],
	Row197 :: [0],
	Row198 :: [0],
	Row199 :: [0],
	Row200 :: [0],
	Row201 :: [0],
	Row202 :: [0],
	Row203 :: [0],
	Row204 :: [0, 245],
	Row205 :: [0],
	Row206 :: [0],
	Row207 :: [0],
	Row208 :: [0],
	Row209 :: [0],
	Row210 :: [0],
	Row211 :: [0, 245],
	Row212 :: [0],
	Row213 :: [0],
	Row214 :: [0],
	Row215 :: [0, 246],
	Row216 :: [0],
	Row217 :: [0],
	Row218 :: [0],
	Row219 :: [0],
	Row220 :: [0],
	Row221 :: [0],
	Row222 :: [0],
	Row223 :: [0],
	Row224 :: [0],
	Row225 :: [0, 82, 245, 246],
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
	Row244 :: [0, 245],
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
	Row338 :: [0, 160],
	Row339 :: [0],
	Row340 :: [0],
	Row341 :: [0],
	Row342 :: [0, 101],
	Row343 :: [0],
	Row344 :: [0],
	Row345 :: [0, 130],
	Row346 :: [0],
	Row347 :: [0],
	Row348 :: [0, 88, 96],
	Row349 :: [0, 94, 99, 124],
	Row350 :: [0, 104],
	Row351 :: [0, 130],
	Row352 :: [0, 124, 126],
	Row353 :: [0, 94, 100],
	Row354 :: [0, 99, 100, 104, 130],
	Row355 :: [0, 104, 105, 127],
	Row356 :: [0, 83, 87, 91, 93, 94, 97, 98, 99, 105, 116, 120, 122, 123, 124, 125, 126, 127, 130],
	Row357 :: [0, 91, 98, 116, 120, 121, 125],
	Row358 :: [0, 91, 92, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 129, 132],
	Row359 :: [0, 96, 97, 99, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 124, 126, 128, 129],
	Row360 :: [0, 96, 98, 100, 101, 103, 104, 105, 106, 115, 116, 117, 118, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130],
	Row361 :: [0, 91, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 115, 116, 117, 118, 119, 121, 122, 124, 125, 128],
	Row362 :: [0, 93, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 125, 126, 127, 128, 129, 130],
	Row363 :: [0, 92, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129],
	Row364 :: [0, 87, 92, 93, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 132, 137],
	Row365 :: [0, 88, 89, 91, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 133, 136, 137, 138],
	Row366 :: [0, 88, 89, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 131, 138],
	Row367 :: [0, 87, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 133, 134, 138, 140],
	Row368 :: [0, 87, 91, 92, 93, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134],
	Row369 :: [0, 86, 88, 89, 90, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 135, 136, 138],
	Row370 :: [0, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 135],
	Row371 :: [0, 12, 87, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 133],
	Row372 :: [0, 12, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136],
	Row373 :: [0, 12, 90, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 136, 137],
	Row374 :: [0, 12, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 136],
	Row375 :: [0],
	Row376 :: [0],
	Row377 :: [0],
	Row378 :: [0],
	Row379 :: [0],
	Row380 :: [0, 12, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 136],
	Row381 :: [0, 12, 84, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 134],
	Row382 :: [0, 12, 87, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 134, 137, 138],
	Row383 :: [0, 12, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132],
	Row384 :: [0, 12, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129],
	Row385 :: [0, 12, 84, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132],
	Row386 :: [0, 12, 90, 91, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 129, 130, 131, 135],
	Row387 :: [0, 12, 91, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131],
	Row388 :: [0, 12, 83, 92, 95, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 130, 136],
	Row389 :: [0, 12, 92, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130],
	Row390 :: [0, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 123, 125, 127, 129, 131],
	Row391 :: [0, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 103, 104, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 124, 125, 126, 127, 128, 129],
	Row392 :: [0, 93, 94, 98, 99, 100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 119, 120, 121, 122, 123, 125, 127, 129, 130],
	Row393 :: [0, 93, 94, 95, 97, 98, 99, 100, 101, 103, 104, 105, 116, 117, 119, 120, 121, 122, 124, 126, 129, 130, 131],
	Row394 :: [0, 92, 97, 99, 100, 101, 105, 106, 115, 116, 117, 118, 119, 120, 121, 122, 125, 128, 131],
	Row395 :: [0, 94, 97, 99, 103, 105, 115, 116, 117, 118, 119, 120, 121, 122, 125, 127],
	Row396 :: [0, 84, 96, 97, 98, 99, 100, 101, 103, 104, 105, 115, 116, 117, 121, 125, 126, 127, 129, 131],
	Row397 :: [0, 97, 98, 104, 125],
	Row398 :: [0, 104],
	Row399 :: [0, 98, 99, 100, 101, 104, 106, 124, 125],
	Row400 :: [0, 97, 99, 125],
	Row401 :: [0, 87, 98, 99, 116, 130, 131],
	Row402 :: [0],
	Row403 :: [0, 120],
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
	Row424 :: [0],
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
	Row462 :: [0],
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
	Row479 :: [0, 245],
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
	Row531 :: [0],
	Row532 :: [0],
	Row533 :: [0],
	Row534 :: [0],
	Row535 :: [0],
	Row536 :: [0],
	Row537 :: [0, 245, 246],
	Row538 :: [0],
	Row539 :: [0, 245],
	Row540 :: [0],
	Row541 :: [0, 246],
	Row542 :: [0, 245, 246],
	Row543 :: [0, 245, 246],
	Row544 :: [0],
	Row545 :: [0],
	Row546 :: [0, 245, 246],
	Row547 :: [0, 245],
	Row548 :: [0, 245],
	Row549 :: [0, 245],
	Row550 :: [0, 246],
	Row551 :: [0],
	Row552 :: [0],
	Row553 :: [0, 245, 246],
	Row554 :: [0, 245, 246],
	Row555 :: [0, 245, 246],
	Row556 :: [0, 245, 246],
	Row557 :: [0, 244, 245, 246],
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
	Freq6 :: [0, 1],
	Freq7 :: [0, 1],
	Freq8 :: [0, 1],
	Freq9 :: [0],
	Freq10 :: [0, 1],
	Freq11 :: [0, 1],
	Freq12 :: [0, 1],
	Freq13 :: [0, 1],
	Freq14 :: [0],
	Freq15 :: [0],
	Freq16 :: [0],
	Freq17 :: [0],
	Freq18 :: [0],
	Freq19 :: [0],
	Freq20 :: [0, 1],
	Freq21 :: [0, 1],
	Freq22 :: [0, 1],
	Freq23 :: [0, 1],
	Freq24 :: [0, 1],
	Freq25 :: [0, 1],
	Freq26 :: [0],
	Freq27 :: [0],
	Freq28 :: [0],
	Freq29 :: [0],
	Freq30 :: [0, 1],
	Freq31 :: [0],
	Freq32 :: [0],
	Freq33 :: [0, 1],
	Freq34 :: [0, 1],
	Freq35 :: [0],
	Freq36 :: [0, 1],
	Freq37 :: [0],
	Freq38 :: [0],
	Freq39 :: [0, 1],
	Freq40 :: [0],
	Freq41 :: [0],
	Freq42 :: [0],
	Freq43 :: [0],
	Freq44 :: [0],
	Freq45 :: [0],
	Freq46 :: [0],
	Freq47 :: [0, 1],
	Freq48 :: [0, 1],
	Freq49 :: [0, 1],
	Freq50 :: [0],
	Freq51 :: [0, 1],
	Freq52 :: [0],
	Freq53 :: [0, 1],
	Freq54 :: [0],
	Freq55 :: [0],
	Freq56 :: [0, 1],
	Freq57 :: [0],
	Freq58 :: [0, 1],
	Freq59 :: [0, 1],
	Freq60 :: [0],
	Freq61 :: [0],
	Freq62 :: [0],
	Freq63 :: [0],
	Freq64 :: [0],
	Freq65 :: [0, 1],
	Freq66 :: [0],
	Freq67 :: [0],
	Freq68 :: [0],
	Freq69 :: [0, 1],
	Freq70 :: [0],
	Freq71 :: [0, 1],
	Freq72 :: [0],
	Freq73 :: [0, 1],
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
	Freq88 :: [0, 1],
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
	Freq165 :: [0, 1],
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
	Freq176 :: [0, 1],
	Freq177 :: [0],
	Freq178 :: [0],
	Freq179 :: [0],
	Freq180 :: [0],
	Freq181 :: [0],
	Freq182 :: [0],
	Freq183 :: [0],
	Freq184 :: [0],
	Freq185 :: [0],
	Freq186 :: [0],
	Freq187 :: [0],
	Freq188 :: [0],
	Freq189 :: [0],
	Freq190 :: [0],
	Freq191 :: [0],
	Freq192 :: [0],
	Freq193 :: [0],
	Freq194 :: [0],
	Freq195 :: [0],
	Freq196 :: [0, 1],
	Freq197 :: [0],
	Freq198 :: [0],
	Freq199 :: [0],
	Freq200 :: [0],
	Freq201 :: [0],
	Freq202 :: [0],
	Freq203 :: [0],
	Freq204 :: [0, 1],
	Freq205 :: [0],
	Freq206 :: [0],
	Freq207 :: [0],
	Freq208 :: [0],
	Freq209 :: [0],
	Freq210 :: [0],
	Freq211 :: [0, 1],
	Freq212 :: [0],
	Freq213 :: [0],
	Freq214 :: [0],
	Freq215 :: [0, 1],
	Freq216 :: [0],
	Freq217 :: [0],
	Freq218 :: [0],
	Freq219 :: [0],
	Freq220 :: [0],
	Freq221 :: [0],
	Freq222 :: [0],
	Freq223 :: [0],
	Freq224 :: [0],
	Freq225 :: [0, 1],
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
	Freq244 :: [0, 1],
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
	Freq338 :: [0, 1],
	Freq339 :: [0],
	Freq340 :: [0],
	Freq341 :: [0],
	Freq342 :: [0, 1],
	Freq343 :: [0],
	Freq344 :: [0],
	Freq345 :: [0, 1],
	Freq346 :: [0],
	Freq347 :: [0],
	Freq348 :: [0, 1],
	Freq349 :: [0, 1],
	Freq350 :: [0, 1],
	Freq351 :: [0, 1],
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
	Freq364 :: [0, 1],
	Freq365 :: [0, 1, 2],
	Freq366 :: [0, 1, 2],
	Freq367 :: [0, 1, 2],
	Freq368 :: [0, 1, 2],
	Freq369 :: [0, 1, 2, 3],
	Freq370 :: [0, 1, 2, 3],
	Freq371 :: [0, 1, 2, 3],
	Freq372 :: [0, 1, 2, 3, 4],
	Freq373 :: [0, 2, 1, 3, 4, 5, 6],
	Freq374 :: [0, 3, 1, 2, 4, 5, 6],
	Freq375 :: [0],
	Freq376 :: [0],
	Freq377 :: [0],
	Freq378 :: [0],
	Freq379 :: [0],
	Freq380 :: [0, 7, 1, 2, 3, 4, 5, 6],
	Freq381 :: [0, 4, 1, 2, 3, 5, 6, 7],
	Freq382 :: [0, 3, 1, 2, 4, 6, 5],
	Freq383 :: [0, 2, 1, 3, 4],
	Freq384 :: [0, 1, 2, 3, 4],
	Freq385 :: [0, 1, 2, 3],
	Freq386 :: [0, 1, 2, 3],
	Freq387 :: [0, 1, 2],
	Freq388 :: [0, 1, 2, 3],
	Freq389 :: [0, 1, 2],
	Freq390 :: [0, 1],
	Freq391 :: [0, 1],
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
	Freq402 :: [0],
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
	Freq424 :: [0],
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
	Freq462 :: [0],
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
	Freq479 :: [0, 1],
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
	Freq531 :: [0],
	Freq532 :: [0],
	Freq533 :: [0],
	Freq534 :: [0],
	Freq535 :: [0],
	Freq536 :: [0],
	Freq537 :: [0, 1],
	Freq538 :: [0],
	Freq539 :: [0, 1],
	Freq540 :: [0],
	Freq541 :: [0, 1],
	Freq542 :: [0, 1],
	Freq543 :: [0, 1],
	Freq544 :: [0],
	Freq545 :: [0],
	Freq546 :: [0, 1],
	Freq547 :: [0, 1],
	Freq548 :: [0, 1],
	Freq549 :: [0, 1],
	Freq550 :: [0, 1],
	Freq551 :: [0],
	Freq552 :: [0],
	Freq553 :: [0, 1],
	Freq554 :: [0, 1, 2],
	Freq555 :: [0, 2, 3],
	Freq556 :: [0, 4],
	Freq557 :: [0, 1, 4, 6],
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

	((Row6 #= 245 and Freq6 #= 1) or
	(Row6 #= 0 and Freq6 #= 0)), 

	((Row7 #= 245 and Freq7 #= 1) or
	(Row7 #= 0 and Freq7 #= 0)), 

	((Row8 #= 245 and Freq8 #= 1) or
	(Row8 #= 246 and Freq8 #= 1) or
	(Row8 #= 0 and Freq8 #= 0)), 

	((Row9 #= 0 and Freq9 #= 0)), 

	((Row10 #= 246 and Freq10 #= 1) or
	(Row10 #= 0 and Freq10 #= 0)), 

	((Row11 #= 245 and Freq11 #= 1) or
	(Row11 #= 0 and Freq11 #= 0)), 

	((Row12 #= 245 and Freq12 #= 1) or
	(Row12 #= 246 and Freq12 #= 1) or
	(Row12 #= 0 and Freq12 #= 0)), 

	((Row13 #= 245 and Freq13 #= 1) or
	(Row13 #= 0 and Freq13 #= 0)), 

	((Row14 #= 0 and Freq14 #= 0)), 

	((Row15 #= 0 and Freq15 #= 0)), 

	((Row16 #= 0 and Freq16 #= 0)), 

	((Row17 #= 0 and Freq17 #= 0)), 

	((Row18 #= 0 and Freq18 #= 0)), 

	((Row19 #= 0 and Freq19 #= 0)), 

	((Row20 #= 245 and Freq20 #= 1) or
	(Row20 #= 246 and Freq20 #= 1) or
	(Row20 #= 0 and Freq20 #= 0)), 

	((Row21 #= 245 and Freq21 #= 1) or
	(Row21 #= 0 and Freq21 #= 0)), 

	((Row22 #= 245 and Freq22 #= 1) or
	(Row22 #= 0 and Freq22 #= 0)), 

	((Row23 #= 245 and Freq23 #= 1) or
	(Row23 #= 0 and Freq23 #= 0)), 

	((Row24 #= 245 and Freq24 #= 1) or
	(Row24 #= 0 and Freq24 #= 0)), 

	((Row25 #= 245 and Freq25 #= 1) or
	(Row25 #= 0 and Freq25 #= 0)), 

	((Row26 #= 0 and Freq26 #= 0)), 

	((Row27 #= 0 and Freq27 #= 0)), 

	((Row28 #= 0 and Freq28 #= 0)), 

	((Row29 #= 0 and Freq29 #= 0)), 

	((Row30 #= 246 and Freq30 #= 1) or
	(Row30 #= 0 and Freq30 #= 0)), 

	((Row31 #= 0 and Freq31 #= 0)), 

	((Row32 #= 0 and Freq32 #= 0)), 

	((Row33 #= 246 and Freq33 #= 1) or
	(Row33 #= 0 and Freq33 #= 0)), 

	((Row34 #= 245 and Freq34 #= 1) or
	(Row34 #= 0 and Freq34 #= 0)), 

	((Row35 #= 0 and Freq35 #= 0)), 

	((Row36 #= 246 and Freq36 #= 1) or
	(Row36 #= 0 and Freq36 #= 0)), 

	((Row37 #= 0 and Freq37 #= 0)), 

	((Row38 #= 0 and Freq38 #= 0)), 

	((Row39 #= 245 and Freq39 #= 1) or
	(Row39 #= 0 and Freq39 #= 0)), 

	((Row40 #= 0 and Freq40 #= 0)), 

	((Row41 #= 0 and Freq41 #= 0)), 

	((Row42 #= 0 and Freq42 #= 0)), 

	((Row43 #= 0 and Freq43 #= 0)), 

	((Row44 #= 0 and Freq44 #= 0)), 

	((Row45 #= 0 and Freq45 #= 0)), 

	((Row46 #= 0 and Freq46 #= 0)), 

	((Row47 #= 245 and Freq47 #= 1) or
	(Row47 #= 0 and Freq47 #= 0)), 

	((Row48 #= 245 and Freq48 #= 1) or
	(Row48 #= 0 and Freq48 #= 0)), 

	((Row49 #= 245 and Freq49 #= 1) or
	(Row49 #= 0 and Freq49 #= 0)), 

	((Row50 #= 0 and Freq50 #= 0)), 

	((Row51 #= 245 and Freq51 #= 1) or
	(Row51 #= 0 and Freq51 #= 0)), 

	((Row52 #= 0 and Freq52 #= 0)), 

	((Row53 #= 245 and Freq53 #= 1) or
	(Row53 #= 0 and Freq53 #= 0)), 

	((Row54 #= 0 and Freq54 #= 0)), 

	((Row55 #= 0 and Freq55 #= 0)), 

	((Row56 #= 246 and Freq56 #= 1) or
	(Row56 #= 0 and Freq56 #= 0)), 

	((Row57 #= 0 and Freq57 #= 0)), 

	((Row58 #= 245 and Freq58 #= 1) or
	(Row58 #= 246 and Freq58 #= 1) or
	(Row58 #= 0 and Freq58 #= 0)), 

	((Row59 #= 245 and Freq59 #= 1) or
	(Row59 #= 0 and Freq59 #= 0)), 

	((Row60 #= 0 and Freq60 #= 0)), 

	((Row61 #= 0 and Freq61 #= 0)), 

	((Row62 #= 0 and Freq62 #= 0)), 

	((Row63 #= 0 and Freq63 #= 0)), 

	((Row64 #= 0 and Freq64 #= 0)), 

	((Row65 #= 245 and Freq65 #= 1) or
	(Row65 #= 0 and Freq65 #= 0)), 

	((Row66 #= 0 and Freq66 #= 0)), 

	((Row67 #= 0 and Freq67 #= 0)), 

	((Row68 #= 0 and Freq68 #= 0)), 

	((Row69 #= 245 and Freq69 #= 1) or
	(Row69 #= 0 and Freq69 #= 0)), 

	((Row70 #= 0 and Freq70 #= 0)), 

	((Row71 #= 245 and Freq71 #= 1) or
	(Row71 #= 0 and Freq71 #= 0)), 

	((Row72 #= 0 and Freq72 #= 0)), 

	((Row73 #= 245 and Freq73 #= 1) or
	(Row73 #= 0 and Freq73 #= 0)), 

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

	((Row88 #= 245 and Freq88 #= 1) or
	(Row88 #= 0 and Freq88 #= 0)), 

	((Row89 #= 0 and Freq89 #= 0)), 

	((Row90 #= 0 and Freq90 #= 0)), 

	((Row91 #= 0 and Freq91 #= 0)), 

	((Row92 #= 0 and Freq92 #= 0)), 

	((Row93 #= 0 and Freq93 #= 0)), 

	((Row94 #= 245 and Freq94 #= 1) or
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

	((Row165 #= 82 and Freq165 #= 1) or
	(Row165 #= 0 and Freq165 #= 0)), 

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

	((Row176 #= 245 and Freq176 #= 1) or
	(Row176 #= 0 and Freq176 #= 0)), 

	((Row177 #= 0 and Freq177 #= 0)), 

	((Row178 #= 0 and Freq178 #= 0)), 

	((Row179 #= 0 and Freq179 #= 0)), 

	((Row180 #= 0 and Freq180 #= 0)), 

	((Row181 #= 0 and Freq181 #= 0)), 

	((Row182 #= 0 and Freq182 #= 0)), 

	((Row183 #= 0 and Freq183 #= 0)), 

	((Row184 #= 0 and Freq184 #= 0)), 

	((Row185 #= 0 and Freq185 #= 0)), 

	((Row186 #= 0 and Freq186 #= 0)), 

	((Row187 #= 0 and Freq187 #= 0)), 

	((Row188 #= 0 and Freq188 #= 0)), 

	((Row189 #= 0 and Freq189 #= 0)), 

	((Row190 #= 0 and Freq190 #= 0)), 

	((Row191 #= 0 and Freq191 #= 0)), 

	((Row192 #= 0 and Freq192 #= 0)), 

	((Row193 #= 0 and Freq193 #= 0)), 

	((Row194 #= 0 and Freq194 #= 0)), 

	((Row195 #= 0 and Freq195 #= 0)), 

	((Row196 #= 99 and Freq196 #= 1) or
	(Row196 #= 0 and Freq196 #= 0)), 

	((Row197 #= 0 and Freq197 #= 0)), 

	((Row198 #= 0 and Freq198 #= 0)), 

	((Row199 #= 0 and Freq199 #= 0)), 

	((Row200 #= 0 and Freq200 #= 0)), 

	((Row201 #= 0 and Freq201 #= 0)), 

	((Row202 #= 0 and Freq202 #= 0)), 

	((Row203 #= 0 and Freq203 #= 0)), 

	((Row204 #= 245 and Freq204 #= 1) or
	(Row204 #= 0 and Freq204 #= 0)), 

	((Row205 #= 0 and Freq205 #= 0)), 

	((Row206 #= 0 and Freq206 #= 0)), 

	((Row207 #= 0 and Freq207 #= 0)), 

	((Row208 #= 0 and Freq208 #= 0)), 

	((Row209 #= 0 and Freq209 #= 0)), 

	((Row210 #= 0 and Freq210 #= 0)), 

	((Row211 #= 245 and Freq211 #= 1) or
	(Row211 #= 0 and Freq211 #= 0)), 

	((Row212 #= 0 and Freq212 #= 0)), 

	((Row213 #= 0 and Freq213 #= 0)), 

	((Row214 #= 0 and Freq214 #= 0)), 

	((Row215 #= 246 and Freq215 #= 1) or
	(Row215 #= 0 and Freq215 #= 0)), 

	((Row216 #= 0 and Freq216 #= 0)), 

	((Row217 #= 0 and Freq217 #= 0)), 

	((Row218 #= 0 and Freq218 #= 0)), 

	((Row219 #= 0 and Freq219 #= 0)), 

	((Row220 #= 0 and Freq220 #= 0)), 

	((Row221 #= 0 and Freq221 #= 0)), 

	((Row222 #= 0 and Freq222 #= 0)), 

	((Row223 #= 0 and Freq223 #= 0)), 

	((Row224 #= 0 and Freq224 #= 0)), 

	((Row225 #= 82 and Freq225 #= 1) or
	(Row225 #= 245 and Freq225 #= 1) or
	(Row225 #= 246 and Freq225 #= 1) or
	(Row225 #= 0 and Freq225 #= 0)), 

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

	((Row244 #= 245 and Freq244 #= 1) or
	(Row244 #= 0 and Freq244 #= 0)), 

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

	((Row338 #= 160 and Freq338 #= 1) or
	(Row338 #= 0 and Freq338 #= 0)), 

	((Row339 #= 0 and Freq339 #= 0)), 

	((Row340 #= 0 and Freq340 #= 0)), 

	((Row341 #= 0 and Freq341 #= 0)), 

	((Row342 #= 101 and Freq342 #= 1) or
	(Row342 #= 0 and Freq342 #= 0)), 

	((Row343 #= 0 and Freq343 #= 0)), 

	((Row344 #= 0 and Freq344 #= 0)), 

	((Row345 #= 130 and Freq345 #= 1) or
	(Row345 #= 0 and Freq345 #= 0)), 

	((Row346 #= 0 and Freq346 #= 0)), 

	((Row347 #= 0 and Freq347 #= 0)), 

	((Row348 #= 88 and Freq348 #= 1) or
	(Row348 #= 96 and Freq348 #= 1) or
	(Row348 #= 0 and Freq348 #= 0)), 

	((Row349 #= 94 and Freq349 #= 1) or
	(Row349 #= 99 and Freq349 #= 1) or
	(Row349 #= 124 and Freq349 #= 1) or
	(Row349 #= 0 and Freq349 #= 0)), 

	((Row350 #= 104 and Freq350 #= 1) or
	(Row350 #= 0 and Freq350 #= 0)), 

	((Row351 #= 130 and Freq351 #= 1) or
	(Row351 #= 0 and Freq351 #= 0)), 

	((Row352 #= 124 and Freq352 #= 1) or
	(Row352 #= 126 and Freq352 #= 1) or
	(Row352 #= 0 and Freq352 #= 0)), 

	((Row353 #= 94 and Freq353 #= 1) or
	(Row353 #= 100 and Freq353 #= 1) or
	(Row353 #= 0 and Freq353 #= 0)), 

	((Row354 #= 99 and Freq354 #= 1) or
	(Row354 #= 100 and Freq354 #= 1) or
	(Row354 #= 104 and Freq354 #= 1) or
	(Row354 #= 130 and Freq354 #= 1) or
	(Row354 #= 0 and Freq354 #= 0)), 

	((Row355 #= 104 and Freq355 #= 1) or
	(Row355 #= 105 and Freq355 #= 1) or
	(Row355 #= 127 and Freq355 #= 1) or
	(Row355 #= 0 and Freq355 #= 0)), 

	((Row356 #= 83 and Freq356 #= 1) or
	(Row356 #= 87 and Freq356 #= 1) or
	(Row356 #= 91 and Freq356 #= 1) or
	(Row356 #= 93 and Freq356 #= 1) or
	(Row356 #= 94 and Freq356 #= 1) or
	(Row356 #= 97 and Freq356 #= 1) or
	(Row356 #= 98 and Freq356 #= 1) or
	(Row356 #= 99 and Freq356 #= 1) or
	(Row356 #= 105 and Freq356 #= 1) or
	(Row356 #= 116 and Freq356 #= 1) or
	(Row356 #= 120 and Freq356 #= 1) or
	(Row356 #= 122 and Freq356 #= 1) or
	(Row356 #= 123 and Freq356 #= 1) or
	(Row356 #= 124 and Freq356 #= 1) or
	(Row356 #= 125 and Freq356 #= 1) or
	(Row356 #= 126 and Freq356 #= 1) or
	(Row356 #= 127 and Freq356 #= 1) or
	(Row356 #= 130 and Freq356 #= 1) or
	(Row356 #= 0 and Freq356 #= 0)), 

	((Row357 #= 91 and Freq357 #= 1) or
	(Row357 #= 98 and Freq357 #= 1) or
	(Row357 #= 116 and Freq357 #= 1) or
	(Row357 #= 120 and Freq357 #= 1) or
	(Row357 #= 121 and Freq357 #= 1) or
	(Row357 #= 125 and Freq357 #= 1) or
	(Row357 #= 0 and Freq357 #= 0)), 

	((Row358 #= 91 and Freq358 #= 1) or
	(Row358 #= 92 and Freq358 #= 1) or
	(Row358 #= 94 and Freq358 #= 1) or
	(Row358 #= 95 and Freq358 #= 1) or
	(Row358 #= 96 and Freq358 #= 1) or
	(Row358 #= 97 and Freq358 #= 1) or
	(Row358 #= 98 and Freq358 #= 1) or
	(Row358 #= 99 and Freq358 #= 1) or
	(Row358 #= 100 and Freq358 #= 1) or
	(Row358 #= 101 and Freq358 #= 1) or
	(Row358 #= 102 and Freq358 #= 1) or
	(Row358 #= 103 and Freq358 #= 1) or
	(Row358 #= 104 and Freq358 #= 1) or
	(Row358 #= 105 and Freq358 #= 1) or
	(Row358 #= 106 and Freq358 #= 1) or
	(Row358 #= 115 and Freq358 #= 1) or
	(Row358 #= 116 and Freq358 #= 1) or
	(Row358 #= 117 and Freq358 #= 1) or
	(Row358 #= 118 and Freq358 #= 1) or
	(Row358 #= 119 and Freq358 #= 1) or
	(Row358 #= 120 and Freq358 #= 1) or
	(Row358 #= 121 and Freq358 #= 1) or
	(Row358 #= 122 and Freq358 #= 1) or
	(Row358 #= 123 and Freq358 #= 1) or
	(Row358 #= 124 and Freq358 #= 1) or
	(Row358 #= 125 and Freq358 #= 1) or
	(Row358 #= 126 and Freq358 #= 1) or
	(Row358 #= 129 and Freq358 #= 1) or
	(Row358 #= 132 and Freq358 #= 1) or
	(Row358 #= 0 and Freq358 #= 0)), 

	((Row359 #= 96 and Freq359 #= 1) or
	(Row359 #= 97 and Freq359 #= 1) or
	(Row359 #= 99 and Freq359 #= 1) or
	(Row359 #= 104 and Freq359 #= 1) or
	(Row359 #= 105 and Freq359 #= 1) or
	(Row359 #= 106 and Freq359 #= 1) or
	(Row359 #= 115 and Freq359 #= 1) or
	(Row359 #= 116 and Freq359 #= 1) or
	(Row359 #= 117 and Freq359 #= 1) or
	(Row359 #= 118 and Freq359 #= 1) or
	(Row359 #= 119 and Freq359 #= 1) or
	(Row359 #= 120 and Freq359 #= 1) or
	(Row359 #= 121 and Freq359 #= 1) or
	(Row359 #= 122 and Freq359 #= 1) or
	(Row359 #= 124 and Freq359 #= 1) or
	(Row359 #= 126 and Freq359 #= 1) or
	(Row359 #= 128 and Freq359 #= 1) or
	(Row359 #= 129 and Freq359 #= 1) or
	(Row359 #= 0 and Freq359 #= 0)), 

	((Row360 #= 96 and Freq360 #= 1) or
	(Row360 #= 98 and Freq360 #= 1) or
	(Row360 #= 100 and Freq360 #= 1) or
	(Row360 #= 101 and Freq360 #= 1) or
	(Row360 #= 103 and Freq360 #= 1) or
	(Row360 #= 104 and Freq360 #= 1) or
	(Row360 #= 105 and Freq360 #= 1) or
	(Row360 #= 106 and Freq360 #= 1) or
	(Row360 #= 115 and Freq360 #= 1) or
	(Row360 #= 116 and Freq360 #= 1) or
	(Row360 #= 117 and Freq360 #= 1) or
	(Row360 #= 118 and Freq360 #= 1) or
	(Row360 #= 120 and Freq360 #= 1) or
	(Row360 #= 121 and Freq360 #= 1) or
	(Row360 #= 122 and Freq360 #= 1) or
	(Row360 #= 123 and Freq360 #= 1) or
	(Row360 #= 124 and Freq360 #= 1) or
	(Row360 #= 125 and Freq360 #= 1) or
	(Row360 #= 126 and Freq360 #= 1) or
	(Row360 #= 127 and Freq360 #= 1) or
	(Row360 #= 128 and Freq360 #= 1) or
	(Row360 #= 129 and Freq360 #= 1) or
	(Row360 #= 130 and Freq360 #= 1) or
	(Row360 #= 0 and Freq360 #= 0)), 

	((Row361 #= 91 and Freq361 #= 1) or
	(Row361 #= 96 and Freq361 #= 1) or
	(Row361 #= 97 and Freq361 #= 1) or
	(Row361 #= 98 and Freq361 #= 1) or
	(Row361 #= 99 and Freq361 #= 1) or
	(Row361 #= 100 and Freq361 #= 1) or
	(Row361 #= 101 and Freq361 #= 1) or
	(Row361 #= 102 and Freq361 #= 1) or
	(Row361 #= 103 and Freq361 #= 1) or
	(Row361 #= 104 and Freq361 #= 1) or
	(Row361 #= 105 and Freq361 #= 1) or
	(Row361 #= 115 and Freq361 #= 1) or
	(Row361 #= 116 and Freq361 #= 1) or
	(Row361 #= 117 and Freq361 #= 1) or
	(Row361 #= 118 and Freq361 #= 1) or
	(Row361 #= 119 and Freq361 #= 1) or
	(Row361 #= 121 and Freq361 #= 1) or
	(Row361 #= 122 and Freq361 #= 1) or
	(Row361 #= 124 and Freq361 #= 1) or
	(Row361 #= 125 and Freq361 #= 1) or
	(Row361 #= 128 and Freq361 #= 1) or
	(Row361 #= 0 and Freq361 #= 0)), 

	((Row362 #= 93 and Freq362 #= 1) or
	(Row362 #= 95 and Freq362 #= 1) or
	(Row362 #= 96 and Freq362 #= 1) or
	(Row362 #= 97 and Freq362 #= 1) or
	(Row362 #= 98 and Freq362 #= 1) or
	(Row362 #= 99 and Freq362 #= 1) or
	(Row362 #= 100 and Freq362 #= 1) or
	(Row362 #= 101 and Freq362 #= 1) or
	(Row362 #= 102 and Freq362 #= 1) or
	(Row362 #= 103 and Freq362 #= 1) or
	(Row362 #= 104 and Freq362 #= 1) or
	(Row362 #= 105 and Freq362 #= 1) or
	(Row362 #= 106 and Freq362 #= 1) or
	(Row362 #= 115 and Freq362 #= 1) or
	(Row362 #= 116 and Freq362 #= 1) or
	(Row362 #= 117 and Freq362 #= 1) or
	(Row362 #= 118 and Freq362 #= 1) or
	(Row362 #= 119 and Freq362 #= 1) or
	(Row362 #= 120 and Freq362 #= 1) or
	(Row362 #= 121 and Freq362 #= 1) or
	(Row362 #= 122 and Freq362 #= 1) or
	(Row362 #= 125 and Freq362 #= 1) or
	(Row362 #= 126 and Freq362 #= 1) or
	(Row362 #= 127 and Freq362 #= 1) or
	(Row362 #= 128 and Freq362 #= 1) or
	(Row362 #= 129 and Freq362 #= 1) or
	(Row362 #= 130 and Freq362 #= 1) or
	(Row362 #= 0 and Freq362 #= 0)), 

	((Row363 #= 92 and Freq363 #= 1) or
	(Row363 #= 94 and Freq363 #= 1) or
	(Row363 #= 95 and Freq363 #= 1) or
	(Row363 #= 96 and Freq363 #= 1) or
	(Row363 #= 97 and Freq363 #= 1) or
	(Row363 #= 98 and Freq363 #= 1) or
	(Row363 #= 99 and Freq363 #= 1) or
	(Row363 #= 100 and Freq363 #= 1) or
	(Row363 #= 101 and Freq363 #= 1) or
	(Row363 #= 102 and Freq363 #= 1) or
	(Row363 #= 103 and Freq363 #= 1) or
	(Row363 #= 104 and Freq363 #= 1) or
	(Row363 #= 105 and Freq363 #= 1) or
	(Row363 #= 106 and Freq363 #= 1) or
	(Row363 #= 115 and Freq363 #= 1) or
	(Row363 #= 116 and Freq363 #= 1) or
	(Row363 #= 117 and Freq363 #= 1) or
	(Row363 #= 118 and Freq363 #= 1) or
	(Row363 #= 119 and Freq363 #= 1) or
	(Row363 #= 120 and Freq363 #= 1) or
	(Row363 #= 121 and Freq363 #= 1) or
	(Row363 #= 122 and Freq363 #= 1) or
	(Row363 #= 123 and Freq363 #= 1) or
	(Row363 #= 124 and Freq363 #= 1) or
	(Row363 #= 125 and Freq363 #= 1) or
	(Row363 #= 126 and Freq363 #= 1) or
	(Row363 #= 127 and Freq363 #= 1) or
	(Row363 #= 128 and Freq363 #= 1) or
	(Row363 #= 129 and Freq363 #= 1) or
	(Row363 #= 0 and Freq363 #= 0)), 

	((Row364 #= 87 and Freq364 #= 1) or
	(Row364 #= 92 and Freq364 #= 1) or
	(Row364 #= 93 and Freq364 #= 1) or
	(Row364 #= 94 and Freq364 #= 1) or
	(Row364 #= 96 and Freq364 #= 1) or
	(Row364 #= 97 and Freq364 #= 1) or
	(Row364 #= 98 and Freq364 #= 1) or
	(Row364 #= 99 and Freq364 #= 1) or
	(Row364 #= 100 and Freq364 #= 1) or
	(Row364 #= 101 and Freq364 #= 1) or
	(Row364 #= 102 and Freq364 #= 1) or
	(Row364 #= 103 and Freq364 #= 1) or
	(Row364 #= 104 and Freq364 #= 1) or
	(Row364 #= 105 and Freq364 #= 1) or
	(Row364 #= 106 and Freq364 #= 1) or
	(Row364 #= 115 and Freq364 #= 1) or
	(Row364 #= 116 and Freq364 #= 1) or
	(Row364 #= 117 and Freq364 #= 1) or
	(Row364 #= 118 and Freq364 #= 1) or
	(Row364 #= 119 and Freq364 #= 1) or
	(Row364 #= 120 and Freq364 #= 1) or
	(Row364 #= 121 and Freq364 #= 1) or
	(Row364 #= 122 and Freq364 #= 1) or
	(Row364 #= 123 and Freq364 #= 1) or
	(Row364 #= 124 and Freq364 #= 1) or
	(Row364 #= 125 and Freq364 #= 1) or
	(Row364 #= 126 and Freq364 #= 1) or
	(Row364 #= 127 and Freq364 #= 1) or
	(Row364 #= 128 and Freq364 #= 1) or
	(Row364 #= 129 and Freq364 #= 1) or
	(Row364 #= 130 and Freq364 #= 1) or
	(Row364 #= 132 and Freq364 #= 1) or
	(Row364 #= 137 and Freq364 #= 1) or
	(Row364 #= 0 and Freq364 #= 0)), 

	((Row365 #= 88 and Freq365 #= 1) or
	(Row365 #= 89 and Freq365 #= 1) or
	(Row365 #= 91 and Freq365 #= 1) or
	(Row365 #= 93 and Freq365 #= 1) or
	(Row365 #= 94 and Freq365 #= 1) or
	(Row365 #= 95 and Freq365 #= 1) or
	(Row365 #= 96 and Freq365 #= 1) or
	(Row365 #= 97 and Freq365 #= 1) or
	(Row365 #= 98 and Freq365 #= 1) or
	(Row365 #= 99 and Freq365 #= 1) or
	(Row365 #= 100 and Freq365 #= 1) or
	(Row365 #= 101 and Freq365 #= 1) or
	(Row365 #= 102 and Freq365 #= 1) or
	(Row365 #= 103 and Freq365 #= 1) or
	(Row365 #= 104 and Freq365 #= 1) or
	(Row365 #= 105 and Freq365 #= 1) or
	(Row365 #= 106 and Freq365 #= 1) or
	(Row365 #= 115 and Freq365 #= 2) or
	(Row365 #= 116 and Freq365 #= 1) or
	(Row365 #= 117 and Freq365 #= 1) or
	(Row365 #= 118 and Freq365 #= 1) or
	(Row365 #= 119 and Freq365 #= 1) or
	(Row365 #= 120 and Freq365 #= 1) or
	(Row365 #= 121 and Freq365 #= 1) or
	(Row365 #= 122 and Freq365 #= 1) or
	(Row365 #= 123 and Freq365 #= 1) or
	(Row365 #= 124 and Freq365 #= 1) or
	(Row365 #= 125 and Freq365 #= 1) or
	(Row365 #= 126 and Freq365 #= 1) or
	(Row365 #= 127 and Freq365 #= 1) or
	(Row365 #= 128 and Freq365 #= 1) or
	(Row365 #= 129 and Freq365 #= 1) or
	(Row365 #= 130 and Freq365 #= 1) or
	(Row365 #= 131 and Freq365 #= 1) or
	(Row365 #= 133 and Freq365 #= 1) or
	(Row365 #= 136 and Freq365 #= 1) or
	(Row365 #= 137 and Freq365 #= 1) or
	(Row365 #= 138 and Freq365 #= 1) or
	(Row365 #= 0 and Freq365 #= 0)), 

	((Row366 #= 88 and Freq366 #= 1) or
	(Row366 #= 89 and Freq366 #= 1) or
	(Row366 #= 92 and Freq366 #= 1) or
	(Row366 #= 93 and Freq366 #= 1) or
	(Row366 #= 94 and Freq366 #= 1) or
	(Row366 #= 95 and Freq366 #= 1) or
	(Row366 #= 96 and Freq366 #= 1) or
	(Row366 #= 97 and Freq366 #= 1) or
	(Row366 #= 98 and Freq366 #= 1) or
	(Row366 #= 99 and Freq366 #= 1) or
	(Row366 #= 100 and Freq366 #= 1) or
	(Row366 #= 101 and Freq366 #= 1) or
	(Row366 #= 102 and Freq366 #= 1) or
	(Row366 #= 103 and Freq366 #= 1) or
	(Row366 #= 104 and Freq366 #= 2) or
	(Row366 #= 105 and Freq366 #= 1) or
	(Row366 #= 106 and Freq366 #= 1) or
	(Row366 #= 115 and Freq366 #= 1) or
	(Row366 #= 116 and Freq366 #= 2) or
	(Row366 #= 117 and Freq366 #= 2) or
	(Row366 #= 118 and Freq366 #= 2) or
	(Row366 #= 119 and Freq366 #= 1) or
	(Row366 #= 120 and Freq366 #= 1) or
	(Row366 #= 121 and Freq366 #= 1) or
	(Row366 #= 122 and Freq366 #= 1) or
	(Row366 #= 123 and Freq366 #= 1) or
	(Row366 #= 124 and Freq366 #= 1) or
	(Row366 #= 125 and Freq366 #= 1) or
	(Row366 #= 126 and Freq366 #= 1) or
	(Row366 #= 127 and Freq366 #= 1) or
	(Row366 #= 128 and Freq366 #= 1) or
	(Row366 #= 129 and Freq366 #= 1) or
	(Row366 #= 131 and Freq366 #= 1) or
	(Row366 #= 138 and Freq366 #= 1) or
	(Row366 #= 0 and Freq366 #= 0)), 

	((Row367 #= 87 and Freq367 #= 1) or
	(Row367 #= 91 and Freq367 #= 1) or
	(Row367 #= 92 and Freq367 #= 1) or
	(Row367 #= 93 and Freq367 #= 1) or
	(Row367 #= 94 and Freq367 #= 1) or
	(Row367 #= 95 and Freq367 #= 1) or
	(Row367 #= 96 and Freq367 #= 1) or
	(Row367 #= 97 and Freq367 #= 1) or
	(Row367 #= 98 and Freq367 #= 1) or
	(Row367 #= 99 and Freq367 #= 1) or
	(Row367 #= 100 and Freq367 #= 1) or
	(Row367 #= 101 and Freq367 #= 1) or
	(Row367 #= 102 and Freq367 #= 1) or
	(Row367 #= 103 and Freq367 #= 1) or
	(Row367 #= 104 and Freq367 #= 2) or
	(Row367 #= 105 and Freq367 #= 2) or
	(Row367 #= 106 and Freq367 #= 1) or
	(Row367 #= 115 and Freq367 #= 1) or
	(Row367 #= 116 and Freq367 #= 2) or
	(Row367 #= 117 and Freq367 #= 1) or
	(Row367 #= 118 and Freq367 #= 1) or
	(Row367 #= 119 and Freq367 #= 1) or
	(Row367 #= 120 and Freq367 #= 2) or
	(Row367 #= 121 and Freq367 #= 1) or
	(Row367 #= 122 and Freq367 #= 1) or
	(Row367 #= 123 and Freq367 #= 1) or
	(Row367 #= 124 and Freq367 #= 1) or
	(Row367 #= 125 and Freq367 #= 1) or
	(Row367 #= 126 and Freq367 #= 1) or
	(Row367 #= 127 and Freq367 #= 1) or
	(Row367 #= 128 and Freq367 #= 1) or
	(Row367 #= 129 and Freq367 #= 1) or
	(Row367 #= 130 and Freq367 #= 1) or
	(Row367 #= 131 and Freq367 #= 1) or
	(Row367 #= 133 and Freq367 #= 1) or
	(Row367 #= 134 and Freq367 #= 1) or
	(Row367 #= 138 and Freq367 #= 1) or
	(Row367 #= 140 and Freq367 #= 1) or
	(Row367 #= 0 and Freq367 #= 0)), 

	((Row368 #= 87 and Freq368 #= 1) or
	(Row368 #= 91 and Freq368 #= 1) or
	(Row368 #= 92 and Freq368 #= 1) or
	(Row368 #= 93 and Freq368 #= 1) or
	(Row368 #= 95 and Freq368 #= 1) or
	(Row368 #= 96 and Freq368 #= 1) or
	(Row368 #= 97 and Freq368 #= 1) or
	(Row368 #= 98 and Freq368 #= 1) or
	(Row368 #= 99 and Freq368 #= 1) or
	(Row368 #= 100 and Freq368 #= 1) or
	(Row368 #= 101 and Freq368 #= 2) or
	(Row368 #= 102 and Freq368 #= 1) or
	(Row368 #= 103 and Freq368 #= 2) or
	(Row368 #= 104 and Freq368 #= 2) or
	(Row368 #= 105 and Freq368 #= 2) or
	(Row368 #= 106 and Freq368 #= 1) or
	(Row368 #= 115 and Freq368 #= 2) or
	(Row368 #= 116 and Freq368 #= 2) or
	(Row368 #= 117 and Freq368 #= 2) or
	(Row368 #= 118 and Freq368 #= 1) or
	(Row368 #= 119 and Freq368 #= 1) or
	(Row368 #= 120 and Freq368 #= 1) or
	(Row368 #= 121 and Freq368 #= 1) or
	(Row368 #= 122 and Freq368 #= 1) or
	(Row368 #= 123 and Freq368 #= 1) or
	(Row368 #= 124 and Freq368 #= 1) or
	(Row368 #= 125 and Freq368 #= 1) or
	(Row368 #= 126 and Freq368 #= 1) or
	(Row368 #= 127 and Freq368 #= 1) or
	(Row368 #= 128 and Freq368 #= 1) or
	(Row368 #= 129 and Freq368 #= 1) or
	(Row368 #= 130 and Freq368 #= 1) or
	(Row368 #= 131 and Freq368 #= 1) or
	(Row368 #= 132 and Freq368 #= 1) or
	(Row368 #= 133 and Freq368 #= 1) or
	(Row368 #= 134 and Freq368 #= 1) or
	(Row368 #= 0 and Freq368 #= 0)), 

	((Row369 #= 86 and Freq369 #= 1) or
	(Row369 #= 88 and Freq369 #= 1) or
	(Row369 #= 89 and Freq369 #= 1) or
	(Row369 #= 90 and Freq369 #= 1) or
	(Row369 #= 92 and Freq369 #= 1) or
	(Row369 #= 93 and Freq369 #= 1) or
	(Row369 #= 94 and Freq369 #= 1) or
	(Row369 #= 95 and Freq369 #= 1) or
	(Row369 #= 96 and Freq369 #= 1) or
	(Row369 #= 97 and Freq369 #= 1) or
	(Row369 #= 98 and Freq369 #= 1) or
	(Row369 #= 99 and Freq369 #= 1) or
	(Row369 #= 100 and Freq369 #= 1) or
	(Row369 #= 101 and Freq369 #= 1) or
	(Row369 #= 102 and Freq369 #= 1) or
	(Row369 #= 103 and Freq369 #= 2) or
	(Row369 #= 104 and Freq369 #= 2) or
	(Row369 #= 105 and Freq369 #= 2) or
	(Row369 #= 106 and Freq369 #= 2) or
	(Row369 #= 115 and Freq369 #= 2) or
	(Row369 #= 116 and Freq369 #= 3) or
	(Row369 #= 117 and Freq369 #= 2) or
	(Row369 #= 118 and Freq369 #= 2) or
	(Row369 #= 119 and Freq369 #= 2) or
	(Row369 #= 120 and Freq369 #= 2) or
	(Row369 #= 121 and Freq369 #= 2) or
	(Row369 #= 122 and Freq369 #= 2) or
	(Row369 #= 123 and Freq369 #= 1) or
	(Row369 #= 124 and Freq369 #= 1) or
	(Row369 #= 125 and Freq369 #= 1) or
	(Row369 #= 126 and Freq369 #= 1) or
	(Row369 #= 127 and Freq369 #= 1) or
	(Row369 #= 128 and Freq369 #= 1) or
	(Row369 #= 129 and Freq369 #= 1) or
	(Row369 #= 130 and Freq369 #= 1) or
	(Row369 #= 131 and Freq369 #= 1) or
	(Row369 #= 135 and Freq369 #= 1) or
	(Row369 #= 136 and Freq369 #= 1) or
	(Row369 #= 138 and Freq369 #= 1) or
	(Row369 #= 0 and Freq369 #= 0)), 

	((Row370 #= 87 and Freq370 #= 1) or
	(Row370 #= 88 and Freq370 #= 1) or
	(Row370 #= 89 and Freq370 #= 1) or
	(Row370 #= 90 and Freq370 #= 1) or
	(Row370 #= 91 and Freq370 #= 1) or
	(Row370 #= 92 and Freq370 #= 1) or
	(Row370 #= 93 and Freq370 #= 1) or
	(Row370 #= 94 and Freq370 #= 1) or
	(Row370 #= 95 and Freq370 #= 1) or
	(Row370 #= 96 and Freq370 #= 1) or
	(Row370 #= 97 and Freq370 #= 1) or
	(Row370 #= 98 and Freq370 #= 1) or
	(Row370 #= 99 and Freq370 #= 1) or
	(Row370 #= 100 and Freq370 #= 1) or
	(Row370 #= 101 and Freq370 #= 1) or
	(Row370 #= 102 and Freq370 #= 1) or
	(Row370 #= 103 and Freq370 #= 2) or
	(Row370 #= 104 and Freq370 #= 2) or
	(Row370 #= 105 and Freq370 #= 3) or
	(Row370 #= 106 and Freq370 #= 2) or
	(Row370 #= 115 and Freq370 #= 3) or
	(Row370 #= 116 and Freq370 #= 3) or
	(Row370 #= 117 and Freq370 #= 3) or
	(Row370 #= 118 and Freq370 #= 2) or
	(Row370 #= 119 and Freq370 #= 2) or
	(Row370 #= 120 and Freq370 #= 2) or
	(Row370 #= 121 and Freq370 #= 2) or
	(Row370 #= 122 and Freq370 #= 2) or
	(Row370 #= 123 and Freq370 #= 1) or
	(Row370 #= 124 and Freq370 #= 1) or
	(Row370 #= 125 and Freq370 #= 2) or
	(Row370 #= 126 and Freq370 #= 1) or
	(Row370 #= 127 and Freq370 #= 1) or
	(Row370 #= 128 and Freq370 #= 1) or
	(Row370 #= 129 and Freq370 #= 1) or
	(Row370 #= 130 and Freq370 #= 1) or
	(Row370 #= 131 and Freq370 #= 1) or
	(Row370 #= 132 and Freq370 #= 1) or
	(Row370 #= 135 and Freq370 #= 1) or
	(Row370 #= 0 and Freq370 #= 0)), 

	((Row371 #= 12 and Freq371 #= 1) or
	(Row371 #= 87 and Freq371 #= 1) or
	(Row371 #= 90 and Freq371 #= 1) or
	(Row371 #= 91 and Freq371 #= 1) or
	(Row371 #= 92 and Freq371 #= 1) or
	(Row371 #= 93 and Freq371 #= 1) or
	(Row371 #= 94 and Freq371 #= 1) or
	(Row371 #= 95 and Freq371 #= 1) or
	(Row371 #= 96 and Freq371 #= 1) or
	(Row371 #= 97 and Freq371 #= 1) or
	(Row371 #= 98 and Freq371 #= 1) or
	(Row371 #= 99 and Freq371 #= 1) or
	(Row371 #= 100 and Freq371 #= 2) or
	(Row371 #= 101 and Freq371 #= 2) or
	(Row371 #= 102 and Freq371 #= 2) or
	(Row371 #= 103 and Freq371 #= 2) or
	(Row371 #= 104 and Freq371 #= 3) or
	(Row371 #= 105 and Freq371 #= 3) or
	(Row371 #= 106 and Freq371 #= 2) or
	(Row371 #= 115 and Freq371 #= 2) or
	(Row371 #= 116 and Freq371 #= 3) or
	(Row371 #= 117 and Freq371 #= 3) or
	(Row371 #= 118 and Freq371 #= 2) or
	(Row371 #= 119 and Freq371 #= 2) or
	(Row371 #= 120 and Freq371 #= 2) or
	(Row371 #= 121 and Freq371 #= 2) or
	(Row371 #= 122 and Freq371 #= 2) or
	(Row371 #= 123 and Freq371 #= 2) or
	(Row371 #= 124 and Freq371 #= 1) or
	(Row371 #= 125 and Freq371 #= 2) or
	(Row371 #= 126 and Freq371 #= 1) or
	(Row371 #= 127 and Freq371 #= 1) or
	(Row371 #= 128 and Freq371 #= 1) or
	(Row371 #= 129 and Freq371 #= 1) or
	(Row371 #= 130 and Freq371 #= 1) or
	(Row371 #= 131 and Freq371 #= 1) or
	(Row371 #= 133 and Freq371 #= 1) or
	(Row371 #= 0 and Freq371 #= 0)), 

	((Row372 #= 12 and Freq372 #= 1) or
	(Row372 #= 91 and Freq372 #= 1) or
	(Row372 #= 92 and Freq372 #= 1) or
	(Row372 #= 93 and Freq372 #= 1) or
	(Row372 #= 94 and Freq372 #= 1) or
	(Row372 #= 95 and Freq372 #= 1) or
	(Row372 #= 96 and Freq372 #= 1) or
	(Row372 #= 97 and Freq372 #= 1) or
	(Row372 #= 98 and Freq372 #= 2) or
	(Row372 #= 99 and Freq372 #= 1) or
	(Row372 #= 100 and Freq372 #= 2) or
	(Row372 #= 101 and Freq372 #= 2) or
	(Row372 #= 102 and Freq372 #= 2) or
	(Row372 #= 103 and Freq372 #= 2) or
	(Row372 #= 104 and Freq372 #= 3) or
	(Row372 #= 105 and Freq372 #= 3) or
	(Row372 #= 106 and Freq372 #= 3) or
	(Row372 #= 115 and Freq372 #= 4) or
	(Row372 #= 116 and Freq372 #= 4) or
	(Row372 #= 117 and Freq372 #= 3) or
	(Row372 #= 118 and Freq372 #= 3) or
	(Row372 #= 119 and Freq372 #= 2) or
	(Row372 #= 120 and Freq372 #= 3) or
	(Row372 #= 121 and Freq372 #= 2) or
	(Row372 #= 122 and Freq372 #= 2) or
	(Row372 #= 123 and Freq372 #= 2) or
	(Row372 #= 124 and Freq372 #= 2) or
	(Row372 #= 125 and Freq372 #= 1) or
	(Row372 #= 126 and Freq372 #= 2) or
	(Row372 #= 127 and Freq372 #= 1) or
	(Row372 #= 128 and Freq372 #= 1) or
	(Row372 #= 129 and Freq372 #= 1) or
	(Row372 #= 130 and Freq372 #= 1) or
	(Row372 #= 131 and Freq372 #= 1) or
	(Row372 #= 132 and Freq372 #= 1) or
	(Row372 #= 133 and Freq372 #= 1) or
	(Row372 #= 134 and Freq372 #= 1) or
	(Row372 #= 135 and Freq372 #= 1) or
	(Row372 #= 136 and Freq372 #= 1) or
	(Row372 #= 0 and Freq372 #= 0)), 

	((Row373 #= 12 and Freq373 #= 2) or
	(Row373 #= 90 and Freq373 #= 1) or
	(Row373 #= 92 and Freq373 #= 1) or
	(Row373 #= 93 and Freq373 #= 1) or
	(Row373 #= 94 and Freq373 #= 1) or
	(Row373 #= 95 and Freq373 #= 1) or
	(Row373 #= 96 and Freq373 #= 2) or
	(Row373 #= 97 and Freq373 #= 1) or
	(Row373 #= 98 and Freq373 #= 1) or
	(Row373 #= 99 and Freq373 #= 2) or
	(Row373 #= 100 and Freq373 #= 2) or
	(Row373 #= 101 and Freq373 #= 2) or
	(Row373 #= 102 and Freq373 #= 2) or
	(Row373 #= 103 and Freq373 #= 3) or
	(Row373 #= 104 and Freq373 #= 4) or
	(Row373 #= 105 and Freq373 #= 4) or
	(Row373 #= 106 and Freq373 #= 4) or
	(Row373 #= 115 and Freq373 #= 5) or
	(Row373 #= 116 and Freq373 #= 6) or
	(Row373 #= 117 and Freq373 #= 5) or
	(Row373 #= 118 and Freq373 #= 3) or
	(Row373 #= 119 and Freq373 #= 3) or
	(Row373 #= 120 and Freq373 #= 3) or
	(Row373 #= 121 and Freq373 #= 2) or
	(Row373 #= 122 and Freq373 #= 3) or
	(Row373 #= 123 and Freq373 #= 2) or
	(Row373 #= 124 and Freq373 #= 1) or
	(Row373 #= 125 and Freq373 #= 2) or
	(Row373 #= 126 and Freq373 #= 1) or
	(Row373 #= 127 and Freq373 #= 1) or
	(Row373 #= 128 and Freq373 #= 1) or
	(Row373 #= 129 and Freq373 #= 1) or
	(Row373 #= 130 and Freq373 #= 1) or
	(Row373 #= 131 and Freq373 #= 1) or
	(Row373 #= 132 and Freq373 #= 1) or
	(Row373 #= 133 and Freq373 #= 1) or
	(Row373 #= 134 and Freq373 #= 1) or
	(Row373 #= 136 and Freq373 #= 1) or
	(Row373 #= 137 and Freq373 #= 1) or
	(Row373 #= 0 and Freq373 #= 0)), 

	((Row374 #= 12 and Freq374 #= 3) or
	(Row374 #= 90 and Freq374 #= 1) or
	(Row374 #= 91 and Freq374 #= 1) or
	(Row374 #= 92 and Freq374 #= 1) or
	(Row374 #= 93 and Freq374 #= 1) or
	(Row374 #= 94 and Freq374 #= 1) or
	(Row374 #= 95 and Freq374 #= 1) or
	(Row374 #= 96 and Freq374 #= 1) or
	(Row374 #= 97 and Freq374 #= 1) or
	(Row374 #= 98 and Freq374 #= 2) or
	(Row374 #= 99 and Freq374 #= 2) or
	(Row374 #= 100 and Freq374 #= 2) or
	(Row374 #= 101 and Freq374 #= 2) or
	(Row374 #= 102 and Freq374 #= 2) or
	(Row374 #= 103 and Freq374 #= 3) or
	(Row374 #= 104 and Freq374 #= 3) or
	(Row374 #= 105 and Freq374 #= 4) or
	(Row374 #= 106 and Freq374 #= 5) or
	(Row374 #= 115 and Freq374 #= 6) or
	(Row374 #= 116 and Freq374 #= 6) or
	(Row374 #= 117 and Freq374 #= 5) or
	(Row374 #= 118 and Freq374 #= 4) or
	(Row374 #= 119 and Freq374 #= 4) or
	(Row374 #= 120 and Freq374 #= 3) or
	(Row374 #= 121 and Freq374 #= 3) or
	(Row374 #= 122 and Freq374 #= 3) or
	(Row374 #= 123 and Freq374 #= 2) or
	(Row374 #= 124 and Freq374 #= 2) or
	(Row374 #= 125 and Freq374 #= 2) or
	(Row374 #= 126 and Freq374 #= 1) or
	(Row374 #= 127 and Freq374 #= 1) or
	(Row374 #= 128 and Freq374 #= 1) or
	(Row374 #= 129 and Freq374 #= 1) or
	(Row374 #= 130 and Freq374 #= 1) or
	(Row374 #= 131 and Freq374 #= 1) or
	(Row374 #= 132 and Freq374 #= 1) or
	(Row374 #= 133 and Freq374 #= 1) or
	(Row374 #= 134 and Freq374 #= 1) or
	(Row374 #= 136 and Freq374 #= 1) or
	(Row374 #= 0 and Freq374 #= 0)), 

	((Row375 #= 0 and Freq375 #= 0)), 

	((Row376 #= 0 and Freq376 #= 0)), 

	((Row377 #= 0 and Freq377 #= 0)), 

	((Row378 #= 0 and Freq378 #= 0)), 

	((Row379 #= 0 and Freq379 #= 0)), 

	((Row380 #= 12 and Freq380 #= 7) or
	(Row380 #= 93 and Freq380 #= 1) or
	(Row380 #= 94 and Freq380 #= 1) or
	(Row380 #= 95 and Freq380 #= 1) or
	(Row380 #= 96 and Freq380 #= 1) or
	(Row380 #= 97 and Freq380 #= 1) or
	(Row380 #= 98 and Freq380 #= 1) or
	(Row380 #= 99 and Freq380 #= 2) or
	(Row380 #= 100 and Freq380 #= 2) or
	(Row380 #= 101 and Freq380 #= 2) or
	(Row380 #= 102 and Freq380 #= 2) or
	(Row380 #= 103 and Freq380 #= 3) or
	(Row380 #= 104 and Freq380 #= 4) or
	(Row380 #= 105 and Freq380 #= 5) or
	(Row380 #= 106 and Freq380 #= 5) or
	(Row380 #= 115 and Freq380 #= 7) or
	(Row380 #= 116 and Freq380 #= 6) or
	(Row380 #= 117 and Freq380 #= 5) or
	(Row380 #= 118 and Freq380 #= 4) or
	(Row380 #= 119 and Freq380 #= 3) or
	(Row380 #= 120 and Freq380 #= 3) or
	(Row380 #= 121 and Freq380 #= 3) or
	(Row380 #= 122 and Freq380 #= 2) or
	(Row380 #= 123 and Freq380 #= 2) or
	(Row380 #= 124 and Freq380 #= 1) or
	(Row380 #= 125 and Freq380 #= 1) or
	(Row380 #= 126 and Freq380 #= 1) or
	(Row380 #= 127 and Freq380 #= 1) or
	(Row380 #= 128 and Freq380 #= 1) or
	(Row380 #= 129 and Freq380 #= 1) or
	(Row380 #= 130 and Freq380 #= 1) or
	(Row380 #= 131 and Freq380 #= 1) or
	(Row380 #= 132 and Freq380 #= 1) or
	(Row380 #= 133 and Freq380 #= 1) or
	(Row380 #= 136 and Freq380 #= 1) or
	(Row380 #= 0 and Freq380 #= 0)), 

	((Row381 #= 12 and Freq381 #= 4) or
	(Row381 #= 84 and Freq381 #= 1) or
	(Row381 #= 88 and Freq381 #= 1) or
	(Row381 #= 89 and Freq381 #= 1) or
	(Row381 #= 90 and Freq381 #= 1) or
	(Row381 #= 91 and Freq381 #= 1) or
	(Row381 #= 92 and Freq381 #= 1) or
	(Row381 #= 93 and Freq381 #= 1) or
	(Row381 #= 94 and Freq381 #= 1) or
	(Row381 #= 95 and Freq381 #= 1) or
	(Row381 #= 96 and Freq381 #= 2) or
	(Row381 #= 97 and Freq381 #= 1) or
	(Row381 #= 98 and Freq381 #= 2) or
	(Row381 #= 99 and Freq381 #= 2) or
	(Row381 #= 100 and Freq381 #= 2) or
	(Row381 #= 101 and Freq381 #= 2) or
	(Row381 #= 102 and Freq381 #= 3) or
	(Row381 #= 103 and Freq381 #= 3) or
	(Row381 #= 104 and Freq381 #= 4) or
	(Row381 #= 105 and Freq381 #= 4) or
	(Row381 #= 106 and Freq381 #= 5) or
	(Row381 #= 115 and Freq381 #= 6) or
	(Row381 #= 116 and Freq381 #= 7) or
	(Row381 #= 117 and Freq381 #= 5) or
	(Row381 #= 118 and Freq381 #= 4) or
	(Row381 #= 119 and Freq381 #= 3) or
	(Row381 #= 120 and Freq381 #= 4) or
	(Row381 #= 121 and Freq381 #= 3) or
	(Row381 #= 122 and Freq381 #= 3) or
	(Row381 #= 123 and Freq381 #= 2) or
	(Row381 #= 124 and Freq381 #= 2) or
	(Row381 #= 125 and Freq381 #= 2) or
	(Row381 #= 126 and Freq381 #= 2) or
	(Row381 #= 127 and Freq381 #= 1) or
	(Row381 #= 128 and Freq381 #= 1) or
	(Row381 #= 129 and Freq381 #= 1) or
	(Row381 #= 130 and Freq381 #= 1) or
	(Row381 #= 131 and Freq381 #= 1) or
	(Row381 #= 132 and Freq381 #= 1) or
	(Row381 #= 134 and Freq381 #= 1) or
	(Row381 #= 0 and Freq381 #= 0)), 

	((Row382 #= 12 and Freq382 #= 3) or
	(Row382 #= 87 and Freq382 #= 1) or
	(Row382 #= 92 and Freq382 #= 1) or
	(Row382 #= 93 and Freq382 #= 1) or
	(Row382 #= 94 and Freq382 #= 1) or
	(Row382 #= 95 and Freq382 #= 1) or
	(Row382 #= 96 and Freq382 #= 1) or
	(Row382 #= 97 and Freq382 #= 2) or
	(Row382 #= 98 and Freq382 #= 2) or
	(Row382 #= 99 and Freq382 #= 2) or
	(Row382 #= 100 and Freq382 #= 2) or
	(Row382 #= 101 and Freq382 #= 2) or
	(Row382 #= 102 and Freq382 #= 2) or
	(Row382 #= 103 and Freq382 #= 3) or
	(Row382 #= 104 and Freq382 #= 3) or
	(Row382 #= 105 and Freq382 #= 4) or
	(Row382 #= 106 and Freq382 #= 3) or
	(Row382 #= 115 and Freq382 #= 6) or
	(Row382 #= 116 and Freq382 #= 5) or
	(Row382 #= 117 and Freq382 #= 5) or
	(Row382 #= 118 and Freq382 #= 4) or
	(Row382 #= 119 and Freq382 #= 3) or
	(Row382 #= 120 and Freq382 #= 3) or
	(Row382 #= 121 and Freq382 #= 2) or
	(Row382 #= 122 and Freq382 #= 2) or
	(Row382 #= 123 and Freq382 #= 2) or
	(Row382 #= 124 and Freq382 #= 2) or
	(Row382 #= 125 and Freq382 #= 2) or
	(Row382 #= 126 and Freq382 #= 1) or
	(Row382 #= 127 and Freq382 #= 1) or
	(Row382 #= 128 and Freq382 #= 1) or
	(Row382 #= 129 and Freq382 #= 1) or
	(Row382 #= 130 and Freq382 #= 1) or
	(Row382 #= 131 and Freq382 #= 1) or
	(Row382 #= 134 and Freq382 #= 1) or
	(Row382 #= 137 and Freq382 #= 1) or
	(Row382 #= 138 and Freq382 #= 1) or
	(Row382 #= 0 and Freq382 #= 0)), 

	((Row383 #= 12 and Freq383 #= 2) or
	(Row383 #= 90 and Freq383 #= 1) or
	(Row383 #= 91 and Freq383 #= 1) or
	(Row383 #= 92 and Freq383 #= 1) or
	(Row383 #= 93 and Freq383 #= 1) or
	(Row383 #= 94 and Freq383 #= 1) or
	(Row383 #= 95 and Freq383 #= 1) or
	(Row383 #= 96 and Freq383 #= 1) or
	(Row383 #= 97 and Freq383 #= 1) or
	(Row383 #= 98 and Freq383 #= 1) or
	(Row383 #= 99 and Freq383 #= 2) or
	(Row383 #= 100 and Freq383 #= 2) or
	(Row383 #= 101 and Freq383 #= 2) or
	(Row383 #= 102 and Freq383 #= 2) or
	(Row383 #= 103 and Freq383 #= 2) or
	(Row383 #= 104 and Freq383 #= 3) or
	(Row383 #= 105 and Freq383 #= 3) or
	(Row383 #= 106 and Freq383 #= 4) or
	(Row383 #= 115 and Freq383 #= 4) or
	(Row383 #= 116 and Freq383 #= 4) or
	(Row383 #= 117 and Freq383 #= 3) or
	(Row383 #= 118 and Freq383 #= 3) or
	(Row383 #= 119 and Freq383 #= 3) or
	(Row383 #= 120 and Freq383 #= 2) or
	(Row383 #= 121 and Freq383 #= 3) or
	(Row383 #= 122 and Freq383 #= 2) or
	(Row383 #= 123 and Freq383 #= 2) or
	(Row383 #= 124 and Freq383 #= 2) or
	(Row383 #= 125 and Freq383 #= 2) or
	(Row383 #= 126 and Freq383 #= 1) or
	(Row383 #= 127 and Freq383 #= 1) or
	(Row383 #= 128 and Freq383 #= 1) or
	(Row383 #= 129 and Freq383 #= 1) or
	(Row383 #= 130 and Freq383 #= 1) or
	(Row383 #= 131 and Freq383 #= 1) or
	(Row383 #= 132 and Freq383 #= 1) or
	(Row383 #= 0 and Freq383 #= 0)), 

	((Row384 #= 12 and Freq384 #= 1) or
	(Row384 #= 93 and Freq384 #= 1) or
	(Row384 #= 94 and Freq384 #= 1) or
	(Row384 #= 95 and Freq384 #= 1) or
	(Row384 #= 96 and Freq384 #= 1) or
	(Row384 #= 97 and Freq384 #= 1) or
	(Row384 #= 98 and Freq384 #= 1) or
	(Row384 #= 99 and Freq384 #= 2) or
	(Row384 #= 100 and Freq384 #= 2) or
	(Row384 #= 101 and Freq384 #= 2) or
	(Row384 #= 102 and Freq384 #= 1) or
	(Row384 #= 103 and Freq384 #= 2) or
	(Row384 #= 104 and Freq384 #= 2) or
	(Row384 #= 105 and Freq384 #= 2) or
	(Row384 #= 106 and Freq384 #= 3) or
	(Row384 #= 115 and Freq384 #= 4) or
	(Row384 #= 116 and Freq384 #= 3) or
	(Row384 #= 117 and Freq384 #= 3) or
	(Row384 #= 118 and Freq384 #= 2) or
	(Row384 #= 119 and Freq384 #= 2) or
	(Row384 #= 120 and Freq384 #= 3) or
	(Row384 #= 121 and Freq384 #= 2) or
	(Row384 #= 122 and Freq384 #= 2) or
	(Row384 #= 123 and Freq384 #= 1) or
	(Row384 #= 124 and Freq384 #= 1) or
	(Row384 #= 125 and Freq384 #= 1) or
	(Row384 #= 126 and Freq384 #= 1) or
	(Row384 #= 127 and Freq384 #= 1) or
	(Row384 #= 128 and Freq384 #= 1) or
	(Row384 #= 129 and Freq384 #= 1) or
	(Row384 #= 0 and Freq384 #= 0)), 

	((Row385 #= 12 and Freq385 #= 1) or
	(Row385 #= 84 and Freq385 #= 1) or
	(Row385 #= 92 and Freq385 #= 1) or
	(Row385 #= 93 and Freq385 #= 1) or
	(Row385 #= 94 and Freq385 #= 1) or
	(Row385 #= 95 and Freq385 #= 1) or
	(Row385 #= 96 and Freq385 #= 1) or
	(Row385 #= 97 and Freq385 #= 1) or
	(Row385 #= 98 and Freq385 #= 1) or
	(Row385 #= 99 and Freq385 #= 2) or
	(Row385 #= 100 and Freq385 #= 2) or
	(Row385 #= 101 and Freq385 #= 2) or
	(Row385 #= 102 and Freq385 #= 1) or
	(Row385 #= 103 and Freq385 #= 2) or
	(Row385 #= 104 and Freq385 #= 2) or
	(Row385 #= 105 and Freq385 #= 3) or
	(Row385 #= 106 and Freq385 #= 3) or
	(Row385 #= 115 and Freq385 #= 3) or
	(Row385 #= 116 and Freq385 #= 3) or
	(Row385 #= 117 and Freq385 #= 2) or
	(Row385 #= 118 and Freq385 #= 2) or
	(Row385 #= 119 and Freq385 #= 2) or
	(Row385 #= 120 and Freq385 #= 2) or
	(Row385 #= 121 and Freq385 #= 2) or
	(Row385 #= 122 and Freq385 #= 2) or
	(Row385 #= 123 and Freq385 #= 1) or
	(Row385 #= 124 and Freq385 #= 1) or
	(Row385 #= 125 and Freq385 #= 1) or
	(Row385 #= 126 and Freq385 #= 1) or
	(Row385 #= 127 and Freq385 #= 1) or
	(Row385 #= 128 and Freq385 #= 1) or
	(Row385 #= 129 and Freq385 #= 1) or
	(Row385 #= 130 and Freq385 #= 1) or
	(Row385 #= 131 and Freq385 #= 1) or
	(Row385 #= 132 and Freq385 #= 1) or
	(Row385 #= 0 and Freq385 #= 0)), 

	((Row386 #= 12 and Freq386 #= 1) or
	(Row386 #= 90 and Freq386 #= 1) or
	(Row386 #= 91 and Freq386 #= 1) or
	(Row386 #= 93 and Freq386 #= 1) or
	(Row386 #= 94 and Freq386 #= 1) or
	(Row386 #= 95 and Freq386 #= 1) or
	(Row386 #= 96 and Freq386 #= 1) or
	(Row386 #= 97 and Freq386 #= 1) or
	(Row386 #= 98 and Freq386 #= 1) or
	(Row386 #= 99 and Freq386 #= 1) or
	(Row386 #= 100 and Freq386 #= 1) or
	(Row386 #= 101 and Freq386 #= 2) or
	(Row386 #= 102 and Freq386 #= 1) or
	(Row386 #= 103 and Freq386 #= 1) or
	(Row386 #= 104 and Freq386 #= 2) or
	(Row386 #= 105 and Freq386 #= 2) or
	(Row386 #= 106 and Freq386 #= 2) or
	(Row386 #= 115 and Freq386 #= 2) or
	(Row386 #= 116 and Freq386 #= 3) or
	(Row386 #= 117 and Freq386 #= 3) or
	(Row386 #= 118 and Freq386 #= 2) or
	(Row386 #= 119 and Freq386 #= 1) or
	(Row386 #= 120 and Freq386 #= 1) or
	(Row386 #= 121 and Freq386 #= 2) or
	(Row386 #= 122 and Freq386 #= 2) or
	(Row386 #= 123 and Freq386 #= 1) or
	(Row386 #= 124 and Freq386 #= 1) or
	(Row386 #= 125 and Freq386 #= 1) or
	(Row386 #= 126 and Freq386 #= 1) or
	(Row386 #= 127 and Freq386 #= 1) or
	(Row386 #= 129 and Freq386 #= 1) or
	(Row386 #= 130 and Freq386 #= 1) or
	(Row386 #= 131 and Freq386 #= 1) or
	(Row386 #= 135 and Freq386 #= 1) or
	(Row386 #= 0 and Freq386 #= 0)), 

	((Row387 #= 12 and Freq387 #= 1) or
	(Row387 #= 91 and Freq387 #= 1) or
	(Row387 #= 93 and Freq387 #= 1) or
	(Row387 #= 94 and Freq387 #= 1) or
	(Row387 #= 95 and Freq387 #= 1) or
	(Row387 #= 96 and Freq387 #= 1) or
	(Row387 #= 97 and Freq387 #= 1) or
	(Row387 #= 98 and Freq387 #= 1) or
	(Row387 #= 99 and Freq387 #= 1) or
	(Row387 #= 100 and Freq387 #= 1) or
	(Row387 #= 101 and Freq387 #= 1) or
	(Row387 #= 102 and Freq387 #= 1) or
	(Row387 #= 103 and Freq387 #= 2) or
	(Row387 #= 104 and Freq387 #= 2) or
	(Row387 #= 105 and Freq387 #= 2) or
	(Row387 #= 106 and Freq387 #= 2) or
	(Row387 #= 115 and Freq387 #= 2) or
	(Row387 #= 116 and Freq387 #= 2) or
	(Row387 #= 117 and Freq387 #= 2) or
	(Row387 #= 118 and Freq387 #= 2) or
	(Row387 #= 119 and Freq387 #= 1) or
	(Row387 #= 120 and Freq387 #= 1) or
	(Row387 #= 121 and Freq387 #= 1) or
	(Row387 #= 122 and Freq387 #= 1) or
	(Row387 #= 123 and Freq387 #= 1) or
	(Row387 #= 124 and Freq387 #= 1) or
	(Row387 #= 125 and Freq387 #= 1) or
	(Row387 #= 126 and Freq387 #= 1) or
	(Row387 #= 127 and Freq387 #= 1) or
	(Row387 #= 128 and Freq387 #= 1) or
	(Row387 #= 129 and Freq387 #= 1) or
	(Row387 #= 130 and Freq387 #= 1) or
	(Row387 #= 131 and Freq387 #= 1) or
	(Row387 #= 0 and Freq387 #= 0)), 

	((Row388 #= 12 and Freq388 #= 1) or
	(Row388 #= 83 and Freq388 #= 1) or
	(Row388 #= 92 and Freq388 #= 1) or
	(Row388 #= 95 and Freq388 #= 1) or
	(Row388 #= 97 and Freq388 #= 1) or
	(Row388 #= 98 and Freq388 #= 1) or
	(Row388 #= 99 and Freq388 #= 1) or
	(Row388 #= 100 and Freq388 #= 1) or
	(Row388 #= 101 and Freq388 #= 1) or
	(Row388 #= 102 and Freq388 #= 1) or
	(Row388 #= 103 and Freq388 #= 1) or
	(Row388 #= 104 and Freq388 #= 1) or
	(Row388 #= 105 and Freq388 #= 2) or
	(Row388 #= 106 and Freq388 #= 1) or
	(Row388 #= 115 and Freq388 #= 2) or
	(Row388 #= 116 and Freq388 #= 3) or
	(Row388 #= 117 and Freq388 #= 2) or
	(Row388 #= 118 and Freq388 #= 1) or
	(Row388 #= 119 and Freq388 #= 2) or
	(Row388 #= 120 and Freq388 #= 1) or
	(Row388 #= 121 and Freq388 #= 1) or
	(Row388 #= 122 and Freq388 #= 1) or
	(Row388 #= 123 and Freq388 #= 1) or
	(Row388 #= 124 and Freq388 #= 1) or
	(Row388 #= 125 and Freq388 #= 1) or
	(Row388 #= 126 and Freq388 #= 1) or
	(Row388 #= 127 and Freq388 #= 1) or
	(Row388 #= 128 and Freq388 #= 1) or
	(Row388 #= 130 and Freq388 #= 1) or
	(Row388 #= 136 and Freq388 #= 1) or
	(Row388 #= 0 and Freq388 #= 0)), 

	((Row389 #= 12 and Freq389 #= 1) or
	(Row389 #= 92 and Freq389 #= 1) or
	(Row389 #= 94 and Freq389 #= 1) or
	(Row389 #= 95 and Freq389 #= 1) or
	(Row389 #= 96 and Freq389 #= 1) or
	(Row389 #= 97 and Freq389 #= 1) or
	(Row389 #= 98 and Freq389 #= 1) or
	(Row389 #= 99 and Freq389 #= 1) or
	(Row389 #= 100 and Freq389 #= 1) or
	(Row389 #= 101 and Freq389 #= 1) or
	(Row389 #= 102 and Freq389 #= 1) or
	(Row389 #= 103 and Freq389 #= 1) or
	(Row389 #= 104 and Freq389 #= 1) or
	(Row389 #= 105 and Freq389 #= 1) or
	(Row389 #= 106 and Freq389 #= 1) or
	(Row389 #= 115 and Freq389 #= 1) or
	(Row389 #= 116 and Freq389 #= 2) or
	(Row389 #= 117 and Freq389 #= 1) or
	(Row389 #= 118 and Freq389 #= 1) or
	(Row389 #= 119 and Freq389 #= 1) or
	(Row389 #= 120 and Freq389 #= 1) or
	(Row389 #= 121 and Freq389 #= 1) or
	(Row389 #= 122 and Freq389 #= 1) or
	(Row389 #= 123 and Freq389 #= 1) or
	(Row389 #= 124 and Freq389 #= 1) or
	(Row389 #= 125 and Freq389 #= 1) or
	(Row389 #= 126 and Freq389 #= 1) or
	(Row389 #= 127 and Freq389 #= 1) or
	(Row389 #= 128 and Freq389 #= 1) or
	(Row389 #= 129 and Freq389 #= 1) or
	(Row389 #= 130 and Freq389 #= 1) or
	(Row389 #= 0 and Freq389 #= 0)), 

	((Row390 #= 94 and Freq390 #= 1) or
	(Row390 #= 95 and Freq390 #= 1) or
	(Row390 #= 96 and Freq390 #= 1) or
	(Row390 #= 97 and Freq390 #= 1) or
	(Row390 #= 98 and Freq390 #= 1) or
	(Row390 #= 99 and Freq390 #= 1) or
	(Row390 #= 100 and Freq390 #= 1) or
	(Row390 #= 101 and Freq390 #= 1) or
	(Row390 #= 102 and Freq390 #= 1) or
	(Row390 #= 103 and Freq390 #= 1) or
	(Row390 #= 104 and Freq390 #= 1) or
	(Row390 #= 105 and Freq390 #= 1) or
	(Row390 #= 106 and Freq390 #= 1) or
	(Row390 #= 115 and Freq390 #= 1) or
	(Row390 #= 116 and Freq390 #= 1) or
	(Row390 #= 117 and Freq390 #= 1) or
	(Row390 #= 118 and Freq390 #= 1) or
	(Row390 #= 119 and Freq390 #= 1) or
	(Row390 #= 120 and Freq390 #= 1) or
	(Row390 #= 121 and Freq390 #= 1) or
	(Row390 #= 122 and Freq390 #= 1) or
	(Row390 #= 123 and Freq390 #= 1) or
	(Row390 #= 125 and Freq390 #= 1) or
	(Row390 #= 127 and Freq390 #= 1) or
	(Row390 #= 129 and Freq390 #= 1) or
	(Row390 #= 131 and Freq390 #= 1) or
	(Row390 #= 0 and Freq390 #= 0)), 

	((Row391 #= 90 and Freq391 #= 1) or
	(Row391 #= 91 and Freq391 #= 1) or
	(Row391 #= 92 and Freq391 #= 1) or
	(Row391 #= 93 and Freq391 #= 1) or
	(Row391 #= 94 and Freq391 #= 1) or
	(Row391 #= 95 and Freq391 #= 1) or
	(Row391 #= 96 and Freq391 #= 1) or
	(Row391 #= 97 and Freq391 #= 1) or
	(Row391 #= 98 and Freq391 #= 1) or
	(Row391 #= 99 and Freq391 #= 1) or
	(Row391 #= 100 and Freq391 #= 1) or
	(Row391 #= 101 and Freq391 #= 1) or
	(Row391 #= 103 and Freq391 #= 1) or
	(Row391 #= 104 and Freq391 #= 1) or
	(Row391 #= 105 and Freq391 #= 1) or
	(Row391 #= 106 and Freq391 #= 1) or
	(Row391 #= 115 and Freq391 #= 1) or
	(Row391 #= 116 and Freq391 #= 1) or
	(Row391 #= 117 and Freq391 #= 1) or
	(Row391 #= 118 and Freq391 #= 1) or
	(Row391 #= 119 and Freq391 #= 1) or
	(Row391 #= 120 and Freq391 #= 1) or
	(Row391 #= 121 and Freq391 #= 1) or
	(Row391 #= 122 and Freq391 #= 1) or
	(Row391 #= 124 and Freq391 #= 1) or
	(Row391 #= 125 and Freq391 #= 1) or
	(Row391 #= 126 and Freq391 #= 1) or
	(Row391 #= 127 and Freq391 #= 1) or
	(Row391 #= 128 and Freq391 #= 1) or
	(Row391 #= 129 and Freq391 #= 1) or
	(Row391 #= 0 and Freq391 #= 0)), 

	((Row392 #= 93 and Freq392 #= 1) or
	(Row392 #= 94 and Freq392 #= 1) or
	(Row392 #= 98 and Freq392 #= 1) or
	(Row392 #= 99 and Freq392 #= 1) or
	(Row392 #= 100 and Freq392 #= 1) or
	(Row392 #= 101 and Freq392 #= 1) or
	(Row392 #= 102 and Freq392 #= 1) or
	(Row392 #= 103 and Freq392 #= 1) or
	(Row392 #= 104 and Freq392 #= 1) or
	(Row392 #= 105 and Freq392 #= 1) or
	(Row392 #= 106 and Freq392 #= 1) or
	(Row392 #= 115 and Freq392 #= 1) or
	(Row392 #= 116 and Freq392 #= 1) or
	(Row392 #= 117 and Freq392 #= 1) or
	(Row392 #= 119 and Freq392 #= 1) or
	(Row392 #= 120 and Freq392 #= 1) or
	(Row392 #= 121 and Freq392 #= 1) or
	(Row392 #= 122 and Freq392 #= 1) or
	(Row392 #= 123 and Freq392 #= 1) or
	(Row392 #= 125 and Freq392 #= 1) or
	(Row392 #= 127 and Freq392 #= 1) or
	(Row392 #= 129 and Freq392 #= 1) or
	(Row392 #= 130 and Freq392 #= 1) or
	(Row392 #= 0 and Freq392 #= 0)), 

	((Row393 #= 93 and Freq393 #= 1) or
	(Row393 #= 94 and Freq393 #= 1) or
	(Row393 #= 95 and Freq393 #= 1) or
	(Row393 #= 97 and Freq393 #= 1) or
	(Row393 #= 98 and Freq393 #= 1) or
	(Row393 #= 99 and Freq393 #= 1) or
	(Row393 #= 100 and Freq393 #= 1) or
	(Row393 #= 101 and Freq393 #= 1) or
	(Row393 #= 103 and Freq393 #= 1) or
	(Row393 #= 104 and Freq393 #= 1) or
	(Row393 #= 105 and Freq393 #= 1) or
	(Row393 #= 116 and Freq393 #= 1) or
	(Row393 #= 117 and Freq393 #= 1) or
	(Row393 #= 119 and Freq393 #= 1) or
	(Row393 #= 120 and Freq393 #= 1) or
	(Row393 #= 121 and Freq393 #= 1) or
	(Row393 #= 122 and Freq393 #= 1) or
	(Row393 #= 124 and Freq393 #= 1) or
	(Row393 #= 126 and Freq393 #= 1) or
	(Row393 #= 129 and Freq393 #= 1) or
	(Row393 #= 130 and Freq393 #= 1) or
	(Row393 #= 131 and Freq393 #= 1) or
	(Row393 #= 0 and Freq393 #= 0)), 

	((Row394 #= 92 and Freq394 #= 1) or
	(Row394 #= 97 and Freq394 #= 1) or
	(Row394 #= 99 and Freq394 #= 1) or
	(Row394 #= 100 and Freq394 #= 1) or
	(Row394 #= 101 and Freq394 #= 1) or
	(Row394 #= 105 and Freq394 #= 1) or
	(Row394 #= 106 and Freq394 #= 1) or
	(Row394 #= 115 and Freq394 #= 1) or
	(Row394 #= 116 and Freq394 #= 1) or
	(Row394 #= 117 and Freq394 #= 1) or
	(Row394 #= 118 and Freq394 #= 1) or
	(Row394 #= 119 and Freq394 #= 1) or
	(Row394 #= 120 and Freq394 #= 1) or
	(Row394 #= 121 and Freq394 #= 1) or
	(Row394 #= 122 and Freq394 #= 1) or
	(Row394 #= 125 and Freq394 #= 1) or
	(Row394 #= 128 and Freq394 #= 1) or
	(Row394 #= 131 and Freq394 #= 1) or
	(Row394 #= 0 and Freq394 #= 0)), 

	((Row395 #= 94 and Freq395 #= 1) or
	(Row395 #= 97 and Freq395 #= 1) or
	(Row395 #= 99 and Freq395 #= 1) or
	(Row395 #= 103 and Freq395 #= 1) or
	(Row395 #= 105 and Freq395 #= 1) or
	(Row395 #= 115 and Freq395 #= 1) or
	(Row395 #= 116 and Freq395 #= 1) or
	(Row395 #= 117 and Freq395 #= 1) or
	(Row395 #= 118 and Freq395 #= 1) or
	(Row395 #= 119 and Freq395 #= 1) or
	(Row395 #= 120 and Freq395 #= 1) or
	(Row395 #= 121 and Freq395 #= 1) or
	(Row395 #= 122 and Freq395 #= 1) or
	(Row395 #= 125 and Freq395 #= 1) or
	(Row395 #= 127 and Freq395 #= 1) or
	(Row395 #= 0 and Freq395 #= 0)), 

	((Row396 #= 84 and Freq396 #= 1) or
	(Row396 #= 96 and Freq396 #= 1) or
	(Row396 #= 97 and Freq396 #= 1) or
	(Row396 #= 98 and Freq396 #= 1) or
	(Row396 #= 99 and Freq396 #= 1) or
	(Row396 #= 100 and Freq396 #= 1) or
	(Row396 #= 101 and Freq396 #= 1) or
	(Row396 #= 103 and Freq396 #= 1) or
	(Row396 #= 104 and Freq396 #= 1) or
	(Row396 #= 105 and Freq396 #= 1) or
	(Row396 #= 115 and Freq396 #= 1) or
	(Row396 #= 116 and Freq396 #= 1) or
	(Row396 #= 117 and Freq396 #= 1) or
	(Row396 #= 121 and Freq396 #= 1) or
	(Row396 #= 125 and Freq396 #= 1) or
	(Row396 #= 126 and Freq396 #= 1) or
	(Row396 #= 127 and Freq396 #= 1) or
	(Row396 #= 129 and Freq396 #= 1) or
	(Row396 #= 131 and Freq396 #= 1) or
	(Row396 #= 0 and Freq396 #= 0)), 

	((Row397 #= 97 and Freq397 #= 1) or
	(Row397 #= 98 and Freq397 #= 1) or
	(Row397 #= 104 and Freq397 #= 1) or
	(Row397 #= 125 and Freq397 #= 1) or
	(Row397 #= 0 and Freq397 #= 0)), 

	((Row398 #= 104 and Freq398 #= 1) or
	(Row398 #= 0 and Freq398 #= 0)), 

	((Row399 #= 98 and Freq399 #= 1) or
	(Row399 #= 99 and Freq399 #= 1) or
	(Row399 #= 100 and Freq399 #= 1) or
	(Row399 #= 101 and Freq399 #= 1) or
	(Row399 #= 104 and Freq399 #= 1) or
	(Row399 #= 106 and Freq399 #= 1) or
	(Row399 #= 124 and Freq399 #= 1) or
	(Row399 #= 125 and Freq399 #= 1) or
	(Row399 #= 0 and Freq399 #= 0)), 

	((Row400 #= 97 and Freq400 #= 1) or
	(Row400 #= 99 and Freq400 #= 1) or
	(Row400 #= 125 and Freq400 #= 1) or
	(Row400 #= 0 and Freq400 #= 0)), 

	((Row401 #= 87 and Freq401 #= 1) or
	(Row401 #= 98 and Freq401 #= 1) or
	(Row401 #= 99 and Freq401 #= 1) or
	(Row401 #= 116 and Freq401 #= 1) or
	(Row401 #= 130 and Freq401 #= 1) or
	(Row401 #= 131 and Freq401 #= 1) or
	(Row401 #= 0 and Freq401 #= 0)), 

	((Row402 #= 0 and Freq402 #= 0)), 

	((Row403 #= 120 and Freq403 #= 1) or
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

	((Row424 #= 0 and Freq424 #= 0)), 

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

	((Row462 #= 0 and Freq462 #= 0)), 

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

	((Row479 #= 245 and Freq479 #= 1) or
	(Row479 #= 0 and Freq479 #= 0)), 

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

	((Row531 #= 0 and Freq531 #= 0)), 

	((Row532 #= 0 and Freq532 #= 0)), 

	((Row533 #= 0 and Freq533 #= 0)), 

	((Row534 #= 0 and Freq534 #= 0)), 

	((Row535 #= 0 and Freq535 #= 0)), 

	((Row536 #= 0 and Freq536 #= 0)), 

	((Row537 #= 245 and Freq537 #= 1) or
	(Row537 #= 246 and Freq537 #= 1) or
	(Row537 #= 0 and Freq537 #= 0)), 

	((Row538 #= 0 and Freq538 #= 0)), 

	((Row539 #= 245 and Freq539 #= 1) or
	(Row539 #= 0 and Freq539 #= 0)), 

	((Row540 #= 0 and Freq540 #= 0)), 

	((Row541 #= 246 and Freq541 #= 1) or
	(Row541 #= 0 and Freq541 #= 0)), 

	((Row542 #= 245 and Freq542 #= 1) or
	(Row542 #= 246 and Freq542 #= 1) or
	(Row542 #= 0 and Freq542 #= 0)), 

	((Row543 #= 245 and Freq543 #= 1) or
	(Row543 #= 246 and Freq543 #= 1) or
	(Row543 #= 0 and Freq543 #= 0)), 

	((Row544 #= 0 and Freq544 #= 0)), 

	((Row545 #= 0 and Freq545 #= 0)), 

	((Row546 #= 245 and Freq546 #= 1) or
	(Row546 #= 246 and Freq546 #= 1) or
	(Row546 #= 0 and Freq546 #= 0)), 

	((Row547 #= 245 and Freq547 #= 1) or
	(Row547 #= 0 and Freq547 #= 0)), 

	((Row548 #= 245 and Freq548 #= 1) or
	(Row548 #= 0 and Freq548 #= 0)), 

	((Row549 #= 245 and Freq549 #= 1) or
	(Row549 #= 0 and Freq549 #= 0)), 

	((Row550 #= 246 and Freq550 #= 1) or
	(Row550 #= 0 and Freq550 #= 0)), 

	((Row551 #= 0 and Freq551 #= 0)), 

	((Row552 #= 0 and Freq552 #= 0)), 

	((Row553 #= 245 and Freq553 #= 1) or
	(Row553 #= 246 and Freq553 #= 1) or
	(Row553 #= 0 and Freq553 #= 0)), 

	((Row554 #= 245 and Freq554 #= 1) or
	(Row554 #= 246 and Freq554 #= 2) or
	(Row554 #= 0 and Freq554 #= 0)), 

	((Row555 #= 245 and Freq555 #= 2) or
	(Row555 #= 246 and Freq555 #= 3) or
	(Row555 #= 0 and Freq555 #= 0)), 

	((Row556 #= 245 and Freq556 #= 4) or
	(Row556 #= 246 and Freq556 #= 4) or
	(Row556 #= 0 and Freq556 #= 0)), 

	((Row557 #= 244 and Freq557 #= 1) or
	(Row557 #= 245 and Freq557 #= 4) or
	(Row557 #= 246 and Freq557 #= 6) or
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
