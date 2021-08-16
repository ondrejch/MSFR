#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-08-16
# GNU/GPL

'''Combined shell+wire resistivity at EOC'''

import agmsfr
import numpy as np

rrr      = agmsfr.Resistivity()
tempC    = 700.0        # Temperature of operation
rho      = {}           # Elemental resistivities at 700C [miloOhm cm]
ele_list = ['Ag', 'Pd', 'Cd']
fw       = {}           # Elemental fractions of EOC wire
fr       = {}           # Elemental fractions of EOC shell
for ele in ele_list:    # Get elemental resistivity
    rho[ele] = rrr.get_rho(ele, tempC)

# Read in EOC wire composition
my_path = "/home/ondrejch/APump/wire_small_jeff33/520"
print('Reading wire composition from ',my_path)
aw = agmsfr.AgWireAnalyzer(my_path+'/msfr')
aw.wdeck_path = my_path
aw.read_wires()
for ele in ele_list:
    fw[ele] = aw.get_EOCfrac(ele)

# Read silver shell EOC
r = 128.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness = 400.0  # 4m reflector
refl  = r + relf_thickness
ag_rs = np.arange(130,refl, 30) # Positions of silver shell
#ag_rs = [130.0,400.0]   # For testing
ar = {}                 # Depleted core object dict
my_path = "/home/ondrejch/APump/MCFR/ag/small_jeff33"

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
print('\n*** RESULTS ***')
print('# rAg[cm]  rho met/all [mOhm cm]    rho-relative to rho_Ag')
for ag_r in ag_rs:
    rho_rat_met = rho_metals[ag_r] / rho['Ag']
    rho_rat_tot = rho_comb[ag_r]   / rho['Ag']
    print(f'  {ag_r:5.1f}    {rho_metals[ag_r]:10.7f} {rho_comb[ag_r]:10.7f}   {rho_rat_tot:12.10f}  {rho_rat_met:12.10f}')

# ---------------------------------------------------------------------- #

results='''
# rAg[cm]  rho met/all [mOhm cm]    rho-relative to rho_Ag
  130.0    14.3461524 14.3462604   2.5017892795  2.5017704365
  160.0    13.8909742 13.8909463   2.4223887878  2.4223936647
  190.0     9.5208384  9.5208190   1.6602990739  1.6603024571
  220.0     7.3057658  7.3057709   1.2740253403  1.2740244505
  250.0     6.3391927  6.3391953   1.1054679219  1.1054674710
  280.0     5.9575691  5.9575771   1.0389190012  1.0389176002
  310.0     5.8150609  5.8150689   1.0140675408  1.0140661413
  340.0     5.7635963  5.7636114   1.0050940607  1.0050914374
  370.0     5.7450326  5.7450320   1.0018540812  1.0018541826
  400.0     5.7383887  5.7384024   1.0006979673  1.0006955658
  430.0     5.7359857  5.7359993   1.0002789038  1.0002765273
  460.0     5.7351023  5.7351120   1.0001241603  1.0001224678
  490.0     5.7347797  5.7347756   1.0000655057  1.0000662109
  520.0     5.7346288  5.7346337   1.0000407612  1.0000399006
'''

gnuplot='''
set terminal pngcairo size 800,600 enhanced font 'Verdana,14'
set out 'res.png'
set title 'Silver wire resistivity increase [%]'
set xlabel 'Reflector thickness [cm]
set ylabel 'Relative resistivity increase [%]
set log y
set grid
plot []'res.dat' u ($1-128):($5-1)*100 w lp notit
set out
'''

