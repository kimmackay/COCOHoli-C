test_matrix.ecl represents the matrix from figure 3 and figure 4 of the paper (./data/test_matrix.tsv)
(generated using: ./generate_model2_intra_2.pl test_matrix.tsv 10 1 6 > ../test_matrix.ecl)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("test_matrix_row.txt", "test_matrix_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -4
Found a solution with cost -9
Found a solution with cost -17
Found a solution with cost -18
Found no solution with cost -34.0 .. -19.0

R = [2, 0, 5, 6, 0, 0]


Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(53)))
Yes (0.00s cpu)

The resultant output files are in ./output and labeled test_matrix_row.txt and test_matrix_freq.txt

===========================================================================

test_5.ecl is the ECLiPSe file for a 5 X 5 fictional matrix (./data/test_5.tsv)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("test_5_row.txt", "test_5_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -4
Found a solution with cost -8
Found a solution with cost -9
Found a solution with cost -10
Found no solution with cost -17.0 .. -11.0

R = [5, 4, 0, 0, 0]


Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(39)))
Yes (0.00s cpu)

The resultant output files are in ./output and labeled test_5_row.txt and test_5_freq.txt

===========================================================================

test_10.ecl is the ECLiPSe file for a 10 X 10 fictional matrix (./data/test_10.tsv)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("test_10_row.txt", "test_10_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -30
Found a solution with cost -90
Found a solution with cost -120
Found a solution with cost -180
Found a solution with cost -270
Found a solution with cost -340
Found a solution with cost -360
Found no solution with cost -660.0 .. -361.0

R = [2, 0, 4, 0, 10, 9, 8, 0, 0, 0]


Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(129)))
Yes (0.09s cpu)

The resultant output files are in ./output and labeled test_10_row.txt and test_10_freq.txt

===========================================================================

test_15.ecl is the ECLiPSe file for a 15 X 15 fictional matrix (./data/test_15.tsv)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("test_15_row.txt", "test_15_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -7
Found a solution with cost -10
Found a solution with cost -15
Found a solution with cost -22
Found a solution with cost -23
Found a solution with cost -24
Found a solution with cost -31
Found a solution with cost -33
Found a solution with cost -37
Found a solution with cost -39
Found a solution with cost -46
Found a solution with cost -48
Found a solution with cost -49
Found a solution with cost -54
Found a solution with cost -56
Found a solution with cost -57
Found a solution with cost -58
Found a solution with cost -61
Found a solution with cost -63
Found a solution with cost -64
Found a solution with cost -69
Found a solution with cost -71
Found a solution with cost -73
Found a solution with cost -82
Found a solution with cost -84
Found a solution with cost -86
Found a solution with cost -88
Found a solution with cost -90
Found a solution with cost -91
Found a solution with cost -92
Found a solution with cost -93
Found a solution with cost -94
Found a solution with cost -103
Found a solution with cost -105
Found no solution with cost -182.0 .. -106.0

R = [10, 0, 8, 7, 15, 14, 0, 0, 13, 0, 12, 0, 0, 0, 0]


Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(269)))
Yes (45.16s cpu)

The resultant output files are in ./output and labeled test_15_row.txt and test_15_freq.txt

===========================================================================

test_17.ecl is the ECLiPSe file for a 17 X 17 fictional matrix (./data/test_17.tsv)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("test_17_row.txt", "test_17_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -9
Found a solution with cost -13
Found a solution with cost -17
Found a solution with cost -26
Found a solution with cost -30
Found a solution with cost -34
Found a solution with cost -43
Found a solution with cost -44
Found a solution with cost -45
Found a solution with cost -50
Found a solution with cost -52
Found a solution with cost -60
Found a solution with cost -67
Found a solution with cost -69
Found a solution with cost -71
Found a solution with cost -72
Found a solution with cost -75
Found a solution with cost -77
Found a solution with cost -78
Found a solution with cost -79
Found a solution with cost -82
Found a solution with cost -84
Found a solution with cost -86
Found a solution with cost -87
Found a solution with cost -90
Found a solution with cost -92
Found a solution with cost -94
Found a solution with cost -97
Found a solution with cost -99
Found a solution with cost -105
Found a solution with cost -107
Found a solution with cost -109
Found a solution with cost -113
Found a solution with cost -115
Found a solution with cost -117
Found a solution with cost -118
Found a solution with cost -119
Found a solution with cost -122
Found a solution with cost -130
Found a solution with cost -132
Found a solution with cost -134
Found no solution with cost -232.0 .. -135.0

R = [0, 10, 9, 7, 17, 16, 0, 15, 0, 0, 14, 13, 0, 0, 0, 0, 0]


Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(339)))
Yes (245.26s cpu)

The resultant output files are in ./output and labeled test_17_row.txt and test_17_freq.txt

===========================================================================

test_18.ecl is the ECLiPSe file for a 18 X 18 fictional matrix (./data/test_18.tsv)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("test_18_row.txt", "test_18_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -9
Found a solution with cost -13
Found a solution with cost -14
Found a solution with cost -15
Found a solution with cost -24
Found a solution with cost -28
Found a solution with cost -30
Found a solution with cost -32
Found a solution with cost -41
Found a solution with cost -45
Found a solution with cost -46
Found a solution with cost -49
Found a solution with cost -58
Found a solution with cost -59
Found a solution with cost -62
Found a solution with cost -63
Found a solution with cost -65
Found a solution with cost -66
Found a solution with cost -75
Found a solution with cost -79
Found a solution with cost -80
Found a solution with cost -82
Found a solution with cost -83
Found a solution with cost -84
Found a solution with cost -85
Found a solution with cost -86
Found a solution with cost -87
Found a solution with cost -88
Found a solution with cost -90
Found a solution with cost -94
Found a solution with cost -95
Found a solution with cost -97
Found a solution with cost -100
Found a solution with cost -101
Found a solution with cost -102
Found a solution with cost -103
Found a solution with cost -104
Found a solution with cost -105
Found a solution with cost -107
Found a solution with cost -115
Found a solution with cost -117
Found a solution with cost -122
Found a solution with cost -132
Found a solution with cost -134
Found a solution with cost -136
Found a solution with cost -137
Found a solution with cost -138
Found a solution with cost -139
Found a solution with cost -140
Found no solution with cost -237.0 .. -141.0

R = [11, 10, 13, 8, 18, 17, 16, 0, 15, 0, 0, 14, 0, 0, 0, 0, 0, 0]

Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(377)))
Yes (584.63s cpu)

The resultant output files are in ./output and labeled test_18_row.txt and test_18_freq.txt

===========================================================================

test_19.ecl is the ECLiPSe file for a 19 X 19 fictional matrix (./data/test_19.tsv)

invoking it in ECLiPSe produces the following output: 

The resultant output files are in ./output and labeled test_19_row.txt and test_19_freq.txt

===========================================================================

test_20.ecl is the ECLiPSe file for a 20 X 20 fictional matrix (./data/test_20.tsv)

invoking it in ECLiPSe produces the following output: 

The resultant output files are in ./output and labeled test_20_row.txt and test_20_freq.txt

===========================================================================

test_25.ecl is the ECLiPSe file for a 25 X 25 fictional matrix (./data/test_25.tsv)

invoking it in ECLiPSe produces the following output: 


The resultant output files are in ./output and labeled row_25.txt and freq_25.txt

===========================================================================

real_10.ecl is the ECLiPSe file for a 10 X 10 submatrix from the s. pombe whole-genome contact map
(generating using: ./generate_model2_intra.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 1 10 > ../real_10.ecl)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("real_10_row.txt", "real_10_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -53
Found a solution with cost -58
Found a solution with cost -111
Found no solution with cost -158.0 .. -112.0

R = [0, 0, 0, 0, 0, 0, 9, 10, 0, 0]


Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(123)))
Yes (0.00s cpu)

The resultant output files are in ./output and labeled real_10_row.txt and real_10_freq.txt

===========================================================================

real_25.ecl is the ECLiPSe file for a 25 X 25 submatrix from the s. pombe whole-genome contact map
(generating using: ./generate_model2_intra.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 1 25 > ../real_25.ecl)

invoking it in ECLiPSe produces the following output: 

[eclipse 2]: maximize("real_25_row.txt", "real_25_freq.txt", R).
Found a solution with cost 0
Found a solution with cost -42
Found a solution with cost -77
Found a solution with cost -88
Found a solution with cost -114
Found a solution with cost -121
Found a solution with cost -132
Found a solution with cost -160
Found a solution with cost -161
Found a solution with cost -165
Found a solution with cost -166
Found a solution with cost -177
Found a solution with cost -178
Found a solution with cost -204
Found a solution with cost -209
Found a solution with cost -220
Found a solution with cost -244
Found a solution with cost -247
Found a solution with cost -249
Found a solution with cost -251
Found a solution with cost -261
Found a solution with cost -262
Found a solution with cost -288
Found a solution with cost -290
Found a solution with cost -301
Found a solution with cost -325
Found a solution with cost -330
Found a solution with cost -341
Found a solution with cost -343
Found a solution with cost -354
Found a solution with cost -373
Found a solution with cost -383
Found a solution with cost -399
Found a solution with cost -401
Found a solution with cost -412
Found a solution with cost -430
Found no solution with cost -776.0 .. -431.0

R = [0, 0, 0, 0, 0, 8, 9, 0, 0, 12, 13, 0, 0, 16, 17, 0, 0, 20, 21, ...]


Delayed goals:
	gfd : gfd_do_propagate(gfd_prob(nvars(687)))
Yes (1885.89s cpu)

The resultant output files are in ./output and labeled real_25_row.txt and real_25_freq.txt

===========================================================================

real_100.ecl is the ECLiPSe file for a 100 X 100 submatrix from the s. pombe whole-genome contact map
(generating using: ./generate_model2_intra.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 1 100 > ../real_100.ecl)

invoking it in ECLiPSe produces the following output: 

The resultant output files are in ./output and labeled real_100_row.txt and real_100_freq.txt

===========================================================================