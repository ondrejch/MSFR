#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-07-28
# note - email ALex Dephilis and Vlad
# GNU/GPL

'''Find critical size of smallest MCFR with eutectic NaCl-UCl3'''

import msfr

#radia = [160,165,170,174]
#radia = [150,153,155,157]
radia = [140,143,145]
relf_thickness = 400.0
my_path = "/home/ondrejch/APump/final_run/small_core_natural/"

for r in radia:
    refl = r + relf_thickness
    mycore = msfr.MSFR(r, refl, 0.1975, "66.66%NaCl+33.34%UCl3")
    mycore.queue     = 'fill'
    mycore.ompcores  = 64
    mycore.histories = 50000
    mycore.s.set_chlorine_37Cl_fraction(0.24)
    mycore.qsub_file = my_path + "/run.sh"
    if r == radia[0]:
        mycore.save_qsub_file()
    mycore.deck_path = my_path + "/%06.1f" % r
    mycore.save_deck()
    print (mycore.deck_path)
    mycore.run_deck()

if False: ''' RESULTS
 $ grep ANA_KEFF */*res.m
0140.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.77722E-01 0.00038  9.70477E-01 0.00037  7.34925E-03 0.00516 ];
0143.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.84753E-01 0.00037  9.77306E-01 0.00036  7.38676E-03 0.00518 ];
0145.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.89578E-01 0.00036  9.82116E-01 0.00036  7.41180E-03 0.00524 ];
0150.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00042E+00 0.00035  9.92909E-01 0.00035  7.52187E-03 0.00542 ];
0153.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00638E+00 0.00037  9.99079E-01 0.00036  7.55009E-03 0.00470 ];
0155.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01141E+00 0.00037  1.00383E+00 0.00036  7.53737E-03 0.00520 ];
0157.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01542E+00 0.00040  1.00777E+00 0.00037  7.61806E-03 0.00503 ];
0160.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.02254E+00 0.00037  1.01469E+00 0.00037  7.64013E-03 0.00516 ];
0165.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.03123E+00 0.00036  1.02367E+00 0.00035  7.69739E-03 0.00487 ];
0170.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.04043E+00 0.00037  1.03276E+00 0.00037  7.80146E-03 0.00477 ];
0174.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.04773E+00 0.00036  1.03988E+00 0.00036  7.82924E-03 0.00477 ];
'''
