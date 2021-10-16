#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-07-10
# GNU/GPL

'''Analysis script for silver depletion due to core neutrons within the reflector'''

import agmsfr
import numpy as np
import subprocess


r = 122.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness = 400.0  # 4m reflector
refl = r + relf_thickness
ag_rs = np.arange(130,refl, 10)             # Positions of silver shell

my_path = "/home/ondrejch/APump/final_run/deplete_small"

a = {}      # Hash with results
# Open all files
for ag_r in ag_rs:
    print("Ag_r = ", ag_r)
    a[ag_r] = agmsfr.AgMSFRAnalyzer(my_path + "/ag_r-" + str(ag_r) +  "/msfr")
    a[ag_r].calc_agfrac()
    a[ag_r].calc_topisos()

f = open('shellEOC.dat','w')
for ag_r in ag_rs:
     print(ag_r-r, a[ag_r].get_EOCfrac('Ag'), a[ag_r].get_EOCfrac('Pd'), a[ag_r].get_EOCfrac('Cd'))
     f.write('\t'.join([str(ag_r-r), str(a[ag_r].get_EOCfrac('Ag')), str(a[ag_r].get_EOCfrac('Pd')), str(a[ag_r].get_EOCfrac('Cd')),'\n']))
#    a[ag_r].plot_topisos(my_path + "/plots/ag_" + str(ag_r)+".png", "Ag_r= "+str(ag_r)+" cm")
#    a[ag_r].plot_multi(my_path + "/plots/fig_ag_" + str(ag_r)+".png", "Ag_r= "+str(ag_r)+" cm")
f.close()

# Make plot
gnuplot_script=b'''
set terminal pngcairo size 800,600 enhanced font 'Verdana,14'
set out 'shellEOC.png'
set title 'Silver depletion due to core neutrons at EOC'
set xlabel 'Reflector thickness [cm]'
set ylabel 'Silver fraction remaining in the shell [%]'
unset log y
set grid
plot [0:300]'shellEOC.dat' u 1:($2*100) w p notit ls 1 pt 7 ps 0.9 lc rgb '#757575', ''  u 1:($2*100) w line notit lc rgb '#787878'
set out
quit
'''
p = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE)
p.stdin.write(gnuplot_script)

