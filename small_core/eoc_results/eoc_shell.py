#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-07-10
# GNU/GPL

'''Analysis script for silver depletion due to core neutrons within the reflector'''

import agmsfr
import numpy as np

r = 122.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness = 400.0  # 4m reflector
refl = r + relf_thickness
ag_rs = np.arange(130,refl, 30)             # Positions of silver shell

my_path = "/home/ondrejch/APump/final_run/deplete_small"

a = {}      # Hash with results
# Open all files
for ag_r in ag_rs:
    print("Ag_r = ", ag_r)
    a[ag_r] = agmsfr.AgMSFRAnalyzer(my_path + "/ag_r-" + str(ag_r) +  "/msfr")
    a[ag_r].calc_agfrac()
    a[ag_r].calc_topisos()

for ag_r in ag_rs:
     print (ag_r, a[ag_r].get_EOCfrac('Ag'), a[ag_r].get_EOCfrac('Pd'), a[ag_r].get_EOCfrac('Cd'))
#    a[ag_r].plot_topisos(my_path + "/plots/ag_" + str(ag_r)+".png", "Ag_r= "+str(ag_r)+" cm")
#    a[ag_r].plot_multi(my_path + "/plots/fig_ag_" + str(ag_r)+".png", "Ag_r= "+str(ag_r)+" cm")


