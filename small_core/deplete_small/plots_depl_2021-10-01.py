#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 10:21:43 2020

@author: o
"""

import numpy as np
import serpentTools


MYPATH:str = '/home/ondrejch/APump/final_run/deplete_small/ag_r-520.0/'
fr=serpentTools.read(MYPATH + "msfr_res.m")
fd=serpentTools.read(MYPATH + "msfr_dep.m")
s=fd.materials['fuelsalt']
bdays = np.interp(20, fd.metadata['burnup'], fd.metadata['days'])
# 20 GWd/t is about 3611 days


#fr.plot('absKeff')
#fr.plot('conversionRatio')
#fr.plot('nubar')
#fr.plot('burnup',['totActivity','actinideActivity','fissionProductActivity'])

# fix missing BU data
s.data['burnup'] = fd.metadata['burnup']
#s.plot('burnup','adens',names=['U235','Pu239','Pu240','Pu241'])
fig=fd.plot('burnup','adens', names=['U235','Pu239','Pu240','Pu241'], materials=['fuelsalt'])
fig.set_xlim([0,20])
fig.set_ylim([0,8.5e-4])
fig.get_figure().savefig(MYPATH + "adens.png", bbox_inches='tight')
fig.get_figure().clf()


fig=fr.plot('burnup',['absKeff'])
fig.set_xlim([0,20])
#fig.set_ylim([0.98,1.07])
fig.get_figure().savefig(MYPATH + "absKeff.png", bbox_inches='tight')
fig.get_figure().clf()


fig=fr.plot('burnup',['conversionRatio'])
fig.set_xlim([0,20])
#fig.set_ylim([0.865,0.925])
fig.get_figure().savefig(MYPATH + "cr.png", bbox_inches='tight')
fig.get_figure().clf()


fig=fr.plot('burnup',['nubar'])
fig.set_xlim([0,20])
#fig.set_ylim([2.48,2.9])
fig.get_figure().savefig(MYPATH + "nubar.png", bbox_inches='tight')
fig.get_figure().clf()


#fig=s.plot('burnup','a',names=['total'])
fig=fd.plot('burnup','a', names=['total'], materials=['fuelsalt'])
fig.set_xlim([-0,20])
fig.get_figure().savefig(MYPATH + "salt_activity.png", bbox_inches='tight')
fig.get_figure().clf()

