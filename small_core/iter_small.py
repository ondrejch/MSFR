#!/usr/bin/python3
# 
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-13
# GNU/GPL

'''Find critical size of smallest MCFR with eutectic NaCl-UCl3'''

import msfr 

#radia = [100,140,180,220,260]
#radia = [110,115,120,125,130]
radia = [126,127,128,129]
relf_thickness = 200.0
my_path = "/home/ondrejch/APump/MCFR/small"

for r in radia:
    refl = r + relf_thickness
    mycore = msfr.MSFR(r, refl, 0.1975, "66.66%NaCl+33.34%UCl3")
    mycore.qsub_file = my_path + "/run.sh"
    if r == radia[0]:
        mycore.save_qsub_file()
    mycore.deck_path = my_path + "/%06.1f" % r
    mycore.save_deck()
    print (mycore.deck_path)
    mycore.run_deck()

if False: ''' RESULTS
 ~/APump/MCFR/small $ grep ANA_KEFF */*res.m
0100.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  8.92240E-01 0.00084  8.86068E-01 0.00084  6.40470E-03 0.01280 ];
0110.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.38855E-01 0.00094  9.32372E-01 0.00090  6.68576E-03 0.01303 ];
0115.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.59543E-01 0.00091  9.52694E-01 0.00088  6.88797E-03 0.01146 ];
0120.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.80060E-01 0.00095  9.72981E-01 0.00094  6.83616E-03 0.01167 ];
0125.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.98685E-01 0.00082  9.91769E-01 0.00080  6.99619E-03 0.01199 ];
0126.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.99963E-01 0.00084  9.93185E-01 0.00082  6.91912E-03 0.01230 ];
0127.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00500E+00 0.00080  9.98261E-01 0.00077  7.03605E-03 0.01318 ];
0128.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00805E+00 0.00086  1.00077E+00 0.00082  7.20752E-03 0.01200 ];
0129.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01301E+00 0.00073  1.00575E+00 0.00074  7.12221E-03 0.01212 ];
0130.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01523E+00 0.00077  1.00784E+00 0.00074  7.12564E-03 0.01139 ];
0140.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.04719E+00 0.00078  1.04005E+00 0.00076  7.31533E-03 0.01148 ];
0180.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.14315E+00 0.00068  1.13520E+00 0.00066  7.89905E-03 0.01194 ];
0220.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.20503E+00 0.00070  1.19665E+00 0.00069  8.33137E-03 0.01154 ];
0260.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.24764E+00 0.00071  1.23864E+00 0.00070  8.80210E-03 0.01121 ];
'''

