#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2021-08-13
# GNU/GPL

'''Analysis script for silver depletion due to delayed neutrons'''

import agmsfr
import numpy as np

#my_path = "/home/ondrejch/APump/final_run/wire_small/520"
my_paths = [
'/home/ondrejch/APump/final_run/wire_small/520',
'/home/ondrejch/APump/final_run/wire_small/520/hs',
'/home/ondrejch/APump/final_run/wire_small/130',
'/home/ondrejch/APump/final_run/wire_small/130/hs']

a = {}

for my_path in my_paths:
    a[my_path] = agmsfr.AgWireAnalyzer(my_path+'/msfr')
    a[my_path].wdeck_path = my_path
    a[my_path].read_wires()
    print(my_path, a[my_path].get_EOCfrac('Ag'), a[my_path].get_EOCfrac('Pd'), a[my_path].get_EOCfrac('Cd') )

