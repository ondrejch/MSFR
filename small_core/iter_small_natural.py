#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-07-28
# GNU/GPL

'''Find critical size of smallest MCFR with eutectic NaCl-UCl3 with natural Cl'''

import msfr

#radia = [140,180,220,260,300]
radia = [170,172,174,176,178]
relf_thickness = 200.0
my_path = "/home/ondrejch/APump/natural"

for r in radia:
    refl = r + relf_thickness
    mycore = msfr.MSFR(r, refl, 0.1975, "66.66%NaCl+33.34%UCl3")
    mycore.queue='fill'
    mycore.s.set_chlorine_37Cl_fraction(0.24)
    mycore.qsub_file = my_path + "/run.sh"
    if r == radia[0]:
        mycore.save_qsub_file()
    mycore.deck_path = my_path + "/%06.1f" % r
    mycore.save_deck()
    print (mycore.deck_path)
    mycore.run_deck()

if False: ''' RESULTS
 ~/APump/MCFR/small $ grep ANA_KEFF */*res.m
0140.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.36033E-01 0.00083  9.28989E-01 0.00080  7.05821E-03 0.01256 ];
0170.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  9.97502E-01 0.00081  9.89912E-01 0.00084  7.47696E-03 0.01217 ];
0172.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00185E+00 0.00077  9.94217E-01 0.00075  7.53817E-03 0.01215 ];
0174.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00612E+00 0.00081  9.98278E-01 0.00079  7.45476E-03 0.01130 ];
0176.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.00718E+00 0.00082  9.99882E-01 0.00080  7.60563E-03 0.01196 ];
0178.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01170E+00 0.00082  1.00414E+00 0.00081  7.68732E-03 0.01091 ];
0180.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.01481E+00 0.00080  1.00764E+00 0.00079  7.57501E-03 0.01165 ];
0220.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.06496E+00 0.00081  1.05736E+00 0.00079  8.08197E-03 0.01099 ];
0260.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.09900E+00 0.00075  1.09069E+00 0.00073  8.11930E-03 0.01090 ];
0300.0/msfr_res.m:ANA_KEFF                  (idx, [1:   6]) = [  1.12334E+00 0.00075  1.11491E+00 0.00075  8.33834E-03 0.01087 ];

'''
