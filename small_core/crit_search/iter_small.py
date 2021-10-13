#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-13
# GNU/GPL

'''Find critical size of smallest MCFR with eutectic NaCl-UCl3'''

import msfr

#radia = [100,120,140,160,180]
#radia = [125,128,130,132]
radia = [122,123,125,121]
#radia = [128]
relf_thickness = 400.0
my_path = "/home/ondrejch/APump/final_run/small_core/"

for r in radia:
    refl = r + relf_thickness
    mycore = msfr.MSFR(r, refl, 0.1975, "66.66%NaCl+33.34%UCl3")
    mycore.qsub_file = my_path + "/run.sh"
    mycore.queue     = 'fill'
    mycore.ompcores  = 64
    mycore.histories = 50000
    if r == radia[0]:
        mycore.save_qsub_file()
    mycore.deck_path = my_path + "/%06.1f" % r
    mycore.save_deck()
    print (mycore.deck_path)
    mycore.run_deck()

if False: ''' RESULTS
 $ grep ANA_KEFF */*res.m
0100.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.14613E-01 0.00038  9.07870E-01 0.00038  6.61892E-03 0.00579 ];
0120.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.92937E-01 0.00037  9.85942E-01 0.00036  7.07255E-03 0.00527 ];
0121.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.96862E-01 0.00038  9.89566E-01 0.00037  7.18263E-03 0.00544 ];
0122.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00067E+00 0.00036  9.93625E-01 0.00035  7.09650E-03 0.00523 ];
0123.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00329E+00 0.00035  9.96020E-01 0.00034  7.17840E-03 0.00552 ];
0125.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01092E+00 0.00035  1.00384E+00 0.00034  7.18190E-03 0.00469 ];
0128.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01944E+00 0.00037  1.01210E+00 0.00037  7.30856E-03 0.00524 ];
0130.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.02547E+00 0.00039  1.01814E+00 0.00038  7.34978E-03 0.00489 ];
0132.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.03228E+00 0.00034  1.02487E+00 0.00035  7.39090E-03 0.00509 ];
0140.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.05449E+00 0.00035  1.04683E+00 0.00035  7.56278E-03 0.00545 ];
0160.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.10448E+00 0.00034  1.09669E+00 0.00035  7.81430E-03 0.00483 ];
0180.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.14459E+00 0.00031  1.13661E+00 0.00030  8.08961E-03 0.00461 ];
'''

