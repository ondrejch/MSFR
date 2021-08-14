#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-08-13
# GNU/GPL

'''DEvelopment playground'''

import agmsfr
import numpy as np

my_path = "/home/ondrejch/APump/wire_small_jeff33/520"
#my_paths = ['/home/ondrejch/APump/wire_small_jeff33/520', '/home/ondrejch/APump/wire_small_jeff33/520/hs',
#'/home/ondrejch/APump/wire_small_jeff33/130', '/home/ondrejch/APump/wire_small_jeff33/130/hs']

a = agmsfr.AgWireAnalyzer(my_path+'/msfr')
a.wdeck_path = my_path
a.read_wires()

ele_list = ['Ag', 'Pd', 'Cd']

r = agmsfr.Resistivity()
f = {}      # fractions
rho = {}    # resistivities
f_metals = 0.0
for ele in ele_list:
    f[ele] = a.get_EOCfrac(ele)
    f_metals += f[ele]
    rho[ele] = r.get_rho(ele, 700)
    print( ele, f[ele], rho[ele] )
print('total metal fraction ', f_metals)
if f_metals >= 1.0:
    raise ValueError('Fraction of metals is >= 1: ', f_metals)

print('Ratio of rhos: ', rho['Pd']/rho['Ag'], rho['Cd']/rho['Ag'] )

rho_metals  = r.get_rho_AgPdCd(f['Ag'], f['Pd'], f['Cd'], 700)
resistivity = r.combine_rho(rho_metals, 10e6,  1.0-f_metals)

print('Resistivities: ', resistivity, rho_metals)

rho_rat_met = rho_metals  / rho['Ag']
rho_rat_tot = resistivity / rho['Ag']

print('Resistivity ratio to silver: ', rho_rat_tot, rho_rat_met)

#a = {}

#for my_path in my_paths:
#    a[my_path] = agmsfr.AgWireAnalyzer(my_path+'/msfr')
#    a[my_path].wdeck_path = my_path
#    a[my_path].read_wires()
#    print(my_path, a[my_path].get_EOCfrac('Ag'), a[my_path].get_EOCfrac('Pd'), a[my_path].get_EOCfrac('Cd') )

