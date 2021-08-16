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
print('# refl[cm]  rho met/all [mOhm cm]    rho-relative to rho_Ag     rho increase [%]')
for ag_r in ag_rs:
    refl_thick  = ag_r - r
    rho_rat_met = rho_metals[ag_r] / rho['Ag']
    rho_rat_tot = rho_comb[ag_r]   / rho['Ag']
    rho_pct_increase = (rho_rat_tot - 1.0) * 100
    print(f'    {refl_thick:5.1f}   {rho_metals[ag_r]:10.7f} {rho_comb[ag_r]:10.7f}   {rho_rat_tot:12.10f}  {rho_rat_met:12.10f} {rho_pct_increase:17.12f}')

# ---------------------------------------------------------------------- #

results='''
# refl[cm]  rho met/all [mOhm cm]    rho-relative to rho_Ag     rho increase [%]
      2.0   14.3461524 14.3462604   2.5017892795  2.5017704365  150.178927954087
     32.0   13.8909742 13.8909463   2.4223887878  2.4223936647  142.238878783911
     62.0    9.5208384  9.5208190   1.6602990739  1.6603024571   66.029907385358
     92.0    7.3057658  7.3057709   1.2740253403  1.2740244505   27.402534033065
    122.0    6.3391927  6.3391953   1.1054679219  1.1054674710   10.546792185840
    152.0    5.9575691  5.9575771   1.0389190012  1.0389176002    3.891900124550
    182.0    5.8150609  5.8150689   1.0140675408  1.0140661413    1.406754079231
    212.0    5.7635963  5.7636114   1.0050940607  1.0050914374    0.509406073834
    242.0    5.7450326  5.7450320   1.0018540812  1.0018541826    0.185408115551
    272.0    5.7383887  5.7384024   1.0006979673  1.0006955658    0.069796726060
    302.0    5.7359857  5.7359993   1.0002789038  1.0002765273    0.027890376420
    332.0    5.7351023  5.7351120   1.0001241603  1.0001224678    0.012416025252
    362.0    5.7347797  5.7347756   1.0000655057  1.0000662109    0.006550567650
    392.0    5.7346288  5.7346337   1.0000407612  1.0000399006    0.004076120720


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
plot []'res.dat' u 1:6 w p notit ls 1 pt 7 ps 0.9, '' u 1:6 w l ls 1 lw 0.4 smooth acs notit
set out
'''

