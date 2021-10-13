#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-13
# GNU/GPL

'''Deplete smallest MCFR with eutectic NaCl-UCl3 and silver in reflector'''

import msfr
import numpy as np

r = 122.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness = 400.0  # 4m reflector
refl      = r + relf_thickness
ag_rs_old = []#np.arange(130,refl, 60)  # Positions of silver shell
ag_rs     = np.arange(130,refl, 10)  # Positions of silver shell

my_path = "/home/ondrejch/APump/final_run/deplete_small"

for ag_r in ag_rs:
    if ag_r in ag_rs_old:
        continue
    mycore = msfr.MSFR(r, refl, 0.1975, "66.66%NaCl+33.34%UCl3", ag_r)
    mycore.power = 1e9      # 1 GWth
    mycore.deplete = 10     # 10 years
    mycore.refuel_flow = 2.824e-10
    mycore.queue = 'fill'
    mycore.ompcores = 64
    mycore.histories = 200000
    mycore.qsub_file = my_path + "/run.sh"
    mycore.save_qsub_file()
    mycore.deck_path = my_path + "/ag_r-" + str(ag_r)
    mycore.save_deck()
    print (mycore.deck_path)
    mycore.run_deck()
