#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-08-16
# GNU/GPL

'''Combined shell+wire resistivity at EOC'''

import agmsfr
import numpy as np
import subprocess

rrr      = agmsfr.Resistivity()
tempC    = 700.0        # Temperature of operation
rho      = {}           # Elemental resistivities at 700C [miloOhm cm]
ele_list = ['Ag', 'Pd', 'Cd']
fw       = {}           # Elemental fractions of EOC wire
fr       = {}           # Elemental fractions of EOC shell
for ele in ele_list:    # Get elemental resistivity
    rho[ele] = rrr.get_rho(ele, tempC)

# Read in EOC wire composition
my_path = "/home/ondrejch/APump/final_run/wire_small/520"
print('Reading wire composition from ',my_path)
aw = agmsfr.AgWireAnalyzer(my_path+'/msfr')
aw.wdeck_path = my_path
aw.read_wires()
for ele in ele_list:
    fw[ele] = aw.get_EOCfrac(ele)

# Read silver shell EOC
r = 122.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness = 400.0  # 4m reflector
refl  = r + relf_thickness
ag_rs = np.arange(130,refl, 10) # Positions of silver shell
#ag_rs = [130.0,400.0]   # For testing
ar = {}                 # Depleted core object dict
my_path = "/home/ondrejch/APump/final_run/deplete_small"

print('Reading shell compositions from ',my_path)
for ag_r in ag_rs:
    print("Ag_r = ", ag_r)
    ar[ag_r] = agmsfr.AgMSFRAnalyzer(my_path + "/ag_r-" + str(ag_r) +  "/msfr")
    ar[ag_r].calc_agfrac()
    ar[ag_r].calc_topisos()
    for ele in ele_list:
        fr[ag_r, ele] = ar[ag_r].get_EOCfrac(ele)

# Combine wire and shell
f          = {}         # Combined elemental fractions from wire and shell
rho_metals = {}         # Resistivity of metals [ele] only [miloOhm cm]
rho_comb   = {}         # Resistivity of metals [ele] only [miloOhm cm]
for ag_r in ag_rs:
    f_metals = 0.0
    for ele in ele_list:
        if ele == 'Ag':     # Silver is depleting
            f[ag_r, ele] = fr[ag_r, ele] * fw[ele]
        else:               # others are building up
            f[ag_r, ele] = fr[ag_r, ele] + fw[ele]
        f_metals += f[ag_r, ele]
    rho_metals[ag_r] = rrr.get_rho_AgPdCd(f[ag_r,'Ag'], f[ag_r,'Pd'], f[ag_r,'Cd'], tempC)
    rho_comb[ag_r]   = rrr.combine_rho(rho_metals[ag_r], 10e6, 1.0 - f_metals)

# Print results
f1 = open('res.dat','w')
f2 = open('res.tex','w')
print('\n*** RESULTS ***')
print('# refl[cm]  rho met/all [mOhm cm]    rho-relative to rho_Ag     rho increase [%]')
f1.write('# refl[cm]  rho met/all [mOhm cm]    rho-relative to rho_Ag     rho increase [%]\n')
for ag_r in ag_rs:
    refl_thick  = ag_r - r
    rho_rat_met = rho_metals[ag_r] / rho['Ag']
    rho_rat_tot = rho_comb[ag_r]   / rho['Ag']
    rho_pct_increase = (rho_rat_tot - 1.0) * 100
    print(f'    {refl_thick:5.1f}   {rho_metals[ag_r]:10.7f} {rho_comb[ag_r]:10.7f}  {rho_rat_met:12.10f}  {rho_rat_tot:12.10f}  {rho_pct_increase:17.12f}')
    f1.write(f'    {refl_thick:5.1f}   {rho_metals[ag_r]:10.7f} {rho_comb[ag_r]:10.7f}  {rho_rat_met:12.10f}  {rho_rat_tot:12.10f}  {rho_pct_increase:17.12f}\n')
    f2.write(f'{refl_thick:5.1f} & {rho_metals[ag_r]:10.7f} & {rho_comb[ag_r]:10.7f} & {rho_rat_met:12.10f} & {rho_rat_tot:12.10f} & {rho_pct_increase:17.12f} \\\\\n')
f1.close()
f2.close()

# ---------------------------------------------------------------------- #

# Make plot
gnuplot_script=b'''
set terminal pngcairo size 800,600 enhanced font 'Verdana,14'
set out 'res.png'
set title 'Silver wire resistivity increase [%]'
set xlabel 'Reflector thickness [cm]'
set ylabel 'Relative resistivity increase [%]'
set log y
set grid
plot [0:300]'res.dat' u 1:6 w p notit ls 1 pt 7 ps 0.9, '' u 1:6 w l ls 1 lw 0.4 smooth acs notit
set out
quit
'''
p = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE)
p.stdin.write(gnuplot_script)

