#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-03-13
# GNU/GPL

'''Refuel rate plots'''

import serpentTools
import numpy as np

r:float = 122.0               # Smallest core has fuel salt radius of 128 cm
relf_thickness:float = 400.0  # 4m reflector
refl:float = r + relf_thickness
ag_r:float = 300.0
#refuel_rates = np.geomspace(1e-12,1e-8, 16)
#refuel_rates = np.linspace(2e-10,4e-10, 8)
refuel_rates = np.linspace(2.8e-10,2.9e-10, 11)

my_path:str = "/home/ondrejch/APump/final_run/refuel_search"

fr = {}
fd = {}
s  = {}

for refuel in refuel_rates:
#    if ag_r in ag_rs_old:
#        continue
    MYPATH = my_path + "/refuel-" + str(refuel)

    fr[refuel] = serpentTools.read(MYPATH + "/msfr_res.m")
    fd[refuel] = serpentTools.read(MYPATH + "/msfr_dep.m")
    s[refuel]=fd[refuel].materials['fuelsalt']
    s[refuel].data['burnup'] = fd[refuel].metadata['burnup']

    fig=fr[refuel].plot('burnup',['absKeff'])
    fig.set_title('Refuel rate ' + str(refuel))
    #fig.set_xlim([0,4])
    #fig.set_ylim([0.98,1.07])
    fig.get_figure().savefig(MYPATH + "absKeff.png", bbox_inches='tight')
    fig.get_figure().clf()
