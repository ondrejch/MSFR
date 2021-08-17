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
ag_rs = np.arange(130,refl, 10) # Positions of silver shell
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
    print(f'    {refl_thick:5.1f}   {rho_metals[ag_r]:10.7f} {rho_comb[ag_r]:10.7f}  {rho_rat_met:12.10f}  {rho_rat_tot:12.10f}  {rho_pct_increase:17.12f}')
# Latex    print(f'{refl_thick:5.1f} & {rho_metals[ag_r]:10.7f} & {rho_comb[ag_r]:10.7f} & {rho_rat_met:12.10f} & {rho_rat_tot:12.10f} & {rho_pct_increase:17.12f} \\\\')

# ---------------------------------------------------------------------- #

results='''
# refl[cm]  rho met/all [mOhm cm]    rho-relative to rho_Ag     rho increase [%]
      2.0   14.3461524 14.3462604  2.5017704365  2.5017892795   150.178927954087
     12.0   16.4860415 16.4860355  2.8749374817  2.8749364440   187.493644400423
     22.0   15.6516863 15.6516211  2.7294374886  2.7294261077   172.942610767480
     32.0   13.8909742 13.8909463  2.4223936647  2.4223887878   142.238878783911
     42.0   12.1567444 12.1567004  2.1199679805  2.1199603100   111.996031003225
     52.0   10.6973064 10.6972810  1.8654621870  1.8654577701    86.545777008896
     62.0    9.5208384  9.5208190  1.6603024571  1.6602990739    66.029907385358
     72.0    8.5905946  8.5905832  1.4980808172  1.4980788165    49.807881645216
     82.0    7.8635853  7.8635739  1.3713004500  1.3712984676    37.129846757436
     92.0    7.3057658  7.3057709  1.2740244505  1.2740253403    27.402534033065
    102.0    6.8854257  6.8854188  1.2007229526  1.2007217518    20.072175178105
    112.0    6.5718124  6.5718234  1.1460331257  1.1460350571    14.603505706049
    122.0    6.3391927  6.3391953  1.1054674710  1.1054679219    10.546792185840
    132.0    6.1699063  6.1699186  1.0759462631  1.0759484132     7.594841319189
    142.0    6.0464655  6.0464774  1.0544198996  1.0544219739     5.442197394822
    152.0    5.9575691  5.9575771  1.0389176002  1.0389190012     3.891900124550
    162.0    5.8937357  5.8937348  1.0277859474  1.0277857795     2.778577947710
    172.0    5.8481103  5.8481249  1.0198294955  1.0198320480     1.983204801021
    182.0    5.8150609  5.8150689  1.0140661413  1.0140675408     1.406754079231
    192.0    5.7920283  5.7920295  1.0100495827  1.0100497836     1.004978355600
    202.0    5.7753323  5.7753328  1.0071380325  1.0071381060     0.713810597041
    212.0    5.7635963  5.7636114  1.0050914374  1.0050940607     0.509406073834
    222.0    5.7551852  5.7551810  1.0036246577  1.0036239265     0.362392653641
    232.0    5.7492998  5.7493150  1.0025983179  1.0026009677     0.260096774749
    242.0    5.7450326  5.7450320  1.0018541826  1.0018540812     0.185408115551
    252.0    5.7420501  5.7420524  1.0013340639  1.0013344719     0.133447194649
    262.0    5.7398987  5.7399115  1.0009589043  1.0009611257     0.096112574369
    272.0    5.7383887  5.7384024  1.0006955658  1.0006979673     0.069796726060
    282.0    5.7372798  5.7372836  1.0005022056  1.0005028655     0.050286554068
    292.0    5.7365198  5.7365321  1.0003696586  1.0003718007     0.037180071340
    302.0    5.7359857  5.7359993  1.0002765273  1.0002789038     0.027890376420
    312.0    5.7355733  5.7355886  1.0002046129  1.0002072751     0.020727507134
    322.0    5.7352937  5.7353073  1.0001558498  1.0001582178     0.015821778749
    332.0    5.7351023  5.7351120  1.0001224678  1.0001241603     0.012416025252
    342.0    5.7349538  5.7349717  1.0000965696  1.0000996995     0.009969950998
    352.0    5.7348484  5.7348611  1.0000781987  1.0000804137     0.008041368432
    362.0    5.7347797  5.7347756  1.0000662109  1.0000655057     0.006550567650
    372.0    5.7347083  5.7347266  1.0000537701  1.0000569592     0.005695921702
    382.0    5.7346609  5.7346813  1.0000455000  1.0000490561     0.004905613396
    392.0    5.7346288  5.7346337  1.0000399006  1.0000407612     0.004076120720
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

