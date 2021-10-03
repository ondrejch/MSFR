#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-03-13
# GNU/GPL

'''Refuel rate search'''

import msfr
import numpy as np

r:float = 122.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness:float = 400.0  # 4m reflector
refl:float = r + relf_thickness
ag_r:float = 300.0
#refuel_rates = np.geomspace(1e-12,1e-8, 16)
#refuel_rates = np.linspace(2e-10,4e-10, 8)
refuel_rates = np.linspace(2.8e-10,2.9e-10, 11)

my_path:str = "/home/ondrejch/APump/final_run/refuel_search"

for refuel in refuel_rates:
    mycore = msfr.MSFR(r, refl, 0.1975, "66.66%NaCl+33.34%UCl3", ag_r)
    mycore.power = 1e9      # 1 GWth
    mycore.deplete = 10     # 10 years
    mycore.refuel_flow = refuel
    mycore.queue = 'fill'
    mycore.ompcores = 32
    mycore.histories = 10000
    mycore.qsub_file = my_path + "/run.sh"
    mycore.save_qsub_file()
    mycore.deck_path = my_path + "/refuel-" + str(refuel)
    mycore.save_deck()
    print (mycore.deck_path)
    mycore.run_deck()
