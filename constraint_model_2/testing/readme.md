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