#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-13
# GNU/GPL

'''Analysis script for silver depletion due to core neutrons within the reflector'''

import agmsfr
import numpy as np

r = 128.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness = 400.0  # 4m reflector
refl = r + relf_thickness
ag_rs = np.arange(130,refl, 30)             # Positions of silver shell
#ag_d = 0.05             # 0.5 mm thick silver shell - default

my_path = "/home/ondrejch/APump/MCFR/ag/small_jeff33"

a = {}      # Hash with results
# Open all files
#for ag_r in ag_rs[0:1]:
for ag_r in ag_rs:
    print("Ag_r = ", ag_r)
    a[ag_r] = agmsfr.AgMSFRAnalyzer(my_path + "/ag_r-" + str(ag_r) +  "/msfr")
    a[ag_r].calc_agfrac()
    a[ag_r].calc_topisos()

for ag_r in ag_rs:
     print (ag_r, a[ag_r].get_EOCfrac('Ag'), a[ag_r].get_EOCfrac('Pd'), a[ag_r].get_EOCfrac('Cd'))
#    a[ag_r].plot_topisos(my_path + "/plots/ag_" + str(ag_r)+".png", "Ag_r= "+str(ag_r)+" cm")
#    a[ag_r].plot_multi(my_path + "/plots/fig_ag_" + str(ag_r)+".png", "Ag_r= "+str(ag_r)+" cm")


