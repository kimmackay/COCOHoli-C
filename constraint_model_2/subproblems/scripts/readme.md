subproblems for the S. pombe whole-genome contact map were generated using the following commands: 

INTRA-INTERACTIONS:

./generate_model2_intra.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 1 558 > ./constraint_model_2_pombe_chr1.ecl

./generate_model2_intra.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 559 1012 > ./constraint_model_2_pombe_chr2.ecl

./generate_model2_intra.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 1013 1258 > ./constraint_model_2_pombe_chr3.ecl

INTER-INTERACTIONS:

./generate_model2_inter.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 1 558 559 1012 > ./constraint_model_2_pombe_chr1_chr2.ecl 

./generate_model2_inter.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 1 558 1013 1258 > ./constraint_model_2_pombe_chr1_chr3.ecl 

./generate_model2_inter.pl GSM1379427_wt_999a-corrected-matrix_hic.tsv 1000 559 1012 1013 1258 > ./constraint_model_2_pombe_chr2_chr3.ecl 