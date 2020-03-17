#!/usr/bin/python3
# 
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-13
# GNU/GPL

'''Deplete smallest MCFR with eutectic NaCl-UCl3'''

import msfr 

r = 128.0               # Smallest core has fuel salt radius of 128 cm 
relf_thickness = 200.0  # 2m reflector
my_path = "/home/ondrejch/APump/MCFR/small"

refl = r + relf_thickness
mycore = msfr.MSFR(r, refl, 0.1975, "66.66%NaCl+33.34%UCl3")
mycore.power = 2e9      # 2 GWth
mycore.deplete = 50     # 50 years
mycore.qsub_file = my_path + "/run.sh"
mycore.save_qsub_file()
mycore.deck_path = my_path + "/%06.1f" % r
mycore.save_deck()
print (mycore.deck_path)
mycore.run_deck()

