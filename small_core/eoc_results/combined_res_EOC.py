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
      8.0   19.0098015 19.0098692  3.3150463037  3.3150581042   231.505810418648
     18.0   18.9228517 18.9228297  3.2998834600  3.2998796137   229.987961371875
     28.0   16.0395193 16.0394457  2.7970701912  2.7970573508   179.705735076809
     38.0   12.5805354 12.5804878  2.1938712762  2.1938629630   119.386296297335
     48.0    9.7831967  9.7831835  1.7060541056  1.7060518060    70.605180602455
     58.0    7.9396179  7.9396052  1.3845594788  1.3845572733    38.455727331849
     68.0    6.8653389  6.8653446  1.1972200880  1.1972210862    19.722108624509
     78.0    6.2889063  6.2889232  1.0966982195  1.0967011681     9.670116811636
     88.0    5.9967973  5.9968011  1.0457584524  1.0457591217     4.575912173434
     98.0    5.8550496  5.8550628  1.0210396134  1.0210419284     2.104192841725
    108.0    5.7887960  5.7888018  1.0094859110  1.0094869232     0.948692319930
    118.0    5.7584076  5.7584205  1.0041865984  1.0041888512     0.418885120826
    128.0    5.7448580  5.7448763  1.0018237221  1.0018269166     0.182691657455
    138.0    5.7390060  5.7390188  1.0008032191  1.0008054476     0.080544760394
    148.0    5.7365064  5.7365106  1.0003673275  1.0003680514     0.036805135856
    158.0    5.7354063  5.7354103  1.0001754858  1.0001761836     0.017618364760
    168.0    5.7349441  5.7349525  1.0000948893  1.0000963397     0.009633970405
    178.0    5.7347508  5.7347686  1.0000611688  1.0000642724     0.006427238581
    188.0    5.7346867  5.7346863  1.0000499936  1.0000499196     0.004991957046
    198.0    5.7346487  5.7346602  1.0000433738  1.0000453733     0.004537330692
    208.0    5.7346308  5.7346479  1.0000402426  1.0000432220     0.004322198603
    218.0    5.7346353  5.7346382  1.0000410375  1.0000415431     0.004154312977
    228.0    5.7346320  5.7346359  1.0000404552  1.0000411430     0.004114297899
    238.0    5.7346310  5.7346353  1.0000402906  1.0000410298     0.004102983394
    248.0    5.7346306  5.7346350  1.0000402155  1.0000409783     0.004097829197
    258.0    5.7346304  5.7346348  1.0000401793  1.0000409534     0.004095344279
    268.0    5.7346304  5.7346349  1.0000401816  1.0000409550     0.004095501005
    278.0    5.7346304  5.7346348  1.0000401804  1.0000409542     0.004095416525
    288.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    298.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    308.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    318.0    5.7346304  5.7346348  1.0000401787  1.0000409530     0.004095298729
    328.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    338.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    348.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    358.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    368.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    378.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    388.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
    398.0    5.7346304  5.7346348  1.0000401753  1.0000409507     0.004095070053
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

