#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-13
# GNU/GPL

import os
import math
from scipy import interpolate
from salts import Salt

do_plots = False
my_debug = False


class MSFRbase(object):
    '''Common base class for the MSFR project'''
    def __init__(self):
        self.tempK:float   = 900.0      # Salt temperature [K]
        self.silver_T:float= self.tempK - 50   # Temperature of the silver wire [K]
        self.power:float   = 3e9        # Core thermal power [W]
        self.deplete:float = 0          # Depletion flag, hacky, see code!
                                        # 0 - no depletion, then in years
                                        # Can be 0, 1, 10, 20, 30, 40 ... n/10 years
        self.lib:str       = '09c'      # CE xsection temp selection
        self.lib_ag:str    = '06c'      # CE xsection temp selection for silver
        self.queue:str     = 'gen6'     # NEcluster torque queue
        self.histories:int = 10000      # Neutron histories per cycle
        self.ompcores:int  = 16         # OMP core count
        self.deck_name:str = 'msfr'     # Serpent input file name
        self.deck_path:str = '/tmp'     # Where to run the lattice deck
        self.qsub_file:str = os.path.expanduser('~/') + '/run.sh'  # qsub script path

    def save_deck(self):
        'Saves Serpent deck into an input file'
        try:
            os.makedirs(self.deck_path, exist_ok = True)
            fh = open(self.deck_path + '/' + self.deck_name, 'w')
            fh.write(self.get_deck())
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to deck file: ",
                  self.deck_path + '/' + self.deck_name)
            print(e)

    def save_qsub_file(self):
        'Writes run file for TORQUE.'
        qsub_content = '''#!/bin/bash
#PBS -V
#PBS -N MSFR_S2
#PBS -q {self.queue}
#PBS -l nodes=1:ppn={self.ompcores}

hostname
rm -f done.dat
cd ${{PBS_O_WORKDIR}}
module load mpi
module load serpent

sss2 -omp {self.ompcores} {self.deck_name} > myout.out
awk 'BEGIN{{ORS="\\t"}} /ANA_KEFF/ || /CONVERSION/ {{print $7" "$8;}}' {self.deck_name}_res.m > done.out
#rm {self.deck_name}.out
'''.format(**locals())
        try:                # Write the deck
            f = open(self.qsub_file, 'w')
            f.write(qsub_content)
            f.close()
        except IOError as e:
            print("Unable to write to qsub file", f)
            print(e)

    def run_deck(self):
        'Runs the deck using qsub_file script'
        if self.queue == 'local':    # Run the deck locally
            os.chdir(self.deck_path)
            os.system(self.qsub_file)
        else:               # Submit the job on the cluster
            os.system('cd ' + self.deck_path + ' && qsub ' + self.qsub_file)



class MSFR(MSFRbase):
    '''Molten Spherical chloride salt Fast Reactor'''
    def __init__(self, r:float=300.0, refl:float=500.0, e:float=0.1083, salt="58%NaCl+42%UCl3", Ag_r:float=-1.0):
        if r<10.0 or refl<r or e>1.0 or e<0.0:  # Reject bad input
            raise ValueError("Bad parameters: ", r, refl, e)

        MSFRbase.__init__(self)         # in Python one has to initilize parent class explicitly
        # core parameters
        self.r:float       = r          # Core radius [cm]
        self.refl:float    = refl       # Outer reflector thickness [cm]
        self.salt_formula:str  = salt   # Salt formula
        self.silver_at_r:float = Ag_r   # Where to put silver semi-shpere [cm]
        self.silver_d:float    = 0.05   # Thickness of silver semi-sphere [cm]
        self.s             = Salt(self.salt_formula, e) # Salt used
        self.s.set_chlorine_37Cl_fraction(0.99999)      # Enriched chlorine-37

    def salt_volume(self) -> float:
        '''Get salt volume, twice the fuel sphere volume'''
        V = (4.0/3.0) * math.pi * self.r**3
        return 2.0 * V

    def get_cells(self) -> str:
        'Cell cards for Serpent input deck'
        cells = '''
%______________cell definitions_____________________________________
cell 11  0  fuelsalt  -1      % fuel salt
cell 31  0  refl       1 -2   % reflector'''
        if self.silver_at_r <= self.r:
            cells += '''
cell 99  0  outside    2      % graveyard
'''
        else:
            cells += '''
cell 20  0  silver     2 -3   % silver
cell 32  0  refl       3 -4   % reflector
cell 99  0  outside    4      % graveyard
'''
        return cells.format(**locals())

    def get_surfaces(self) -> str:
        'Surface cards for Serpent input deck'
        surfaces = '''
%______________surface definitions__________________________________
surf 1   sph  0.0 0.0 0.0 {self.r}      % fuel salt radius'''
        if self.silver_at_r <= self.r:
            surfaces += '''
surf 2   sph  0.0 0.0 0.0 {self.refl}   % reflector
'''
        else:
            silver_r_max = self.silver_at_r + self.silver_d
            surfaces += '''
surf 2   sph  0.0 0.0 0.0 {self.silver_at_r}   % reflector
surf 3   sph  0.0 0.0 0.0 {silver_r_max}       % silver
surf 4   sph  0.0 0.0 0.0 {self.refl}   % reflector
'''
        return surfaces.format(**locals())

    def get_materials(self) -> str:
        'Material definitions, non-salt'
        refl_lib  = self.lib
        materials = '''
% Iron reflector [density 7.874/((1+680*12e-6)^3)]
mat refl   -7.68435 tmp 900 rgb 128 128 178
26054.{refl_lib}  -0.058450   %  Fe
26056.{refl_lib}  -0.917540   %  Fe
26057.{refl_lib}  -0.021190   %  Fe
26058.{refl_lib}  -0.002820   %  Fe
'''
# https://www.sciencedirect.com/science/article/abs/pii/0022190262801882
        silver_density = 10.465 - 9.967e-4*self.silver_T # [g/cm^3]

        if self.silver_at_r > 0.0 and self.silver_at_r <= self.r:
            raise ValueError('Silver shell inside fuel ', self.silver_at_r, self.r)
        if self.silver_at_r >= self.refl:
            raise ValueError('Silver shell outside reflector ', self.silver_at_r, self.refl)
        if self.silver_at_r > self.r and self.silver_at_r < self.refl:
            materials += '''
% Silver
mat silver -{silver_density} tmp {self.silver_T} rgb 10 10 10 burn 1
47107.{self.lib_ag}  -0.51839    % Ag
47109.{self.lib_ag}  -0.48161    % Ag
'''
        return materials.format(**locals())

    def get_data_cards(self) -> str:
        'Data cards for the reactor'
        fs_volume = self.salt_volume()
        data_cards = '''
% Fuel salt volume
set mvol fuelsalt 0 {fs_volume}

% Power in thermal W
set power {self.power}

% Boundary condition
set bc 1

% Analog reaction rate
% set arr 2

% Neutron population and criticality cycles
set pop {self.histories} 240 40

% Turning off group constant generation hastens the calculation
set gcu -1
'''
        data_cards += '''
% Data Libraries
set acelib "sss_endfb7u.sssdir"
set declib "sss_endfb7.dec"
set nfylib "sss_endfb7.nfy"
'''
        if do_plots:
            data_cards += '''
% Plots
plot 3 1500 1500
% mesh 3 1500 1500
'''
        return data_cards.format(**locals())

    def get_repr_cards(self) -> str:
        'Reprocessing setup'
        refuel_lib = self.lib
        repr_cards = '''
%___________Reprocessing___________
% First we need some extra materials to do depletion with reprocessing correctly.

% stockpile of extra U1
mat U_stock -3.5096 burn 1 vol 1e8
17037.{refuel_lib} -0.3707532563
11023.{refuel_lib} -0.0633565017
92235.{refuel_lib} -0.05553085588
92238.{refuel_lib} -0.49977770291999996

% tanks for offgases
mat offgastankcore 0.0007 burn 1 vol 1e6
2004.{refuel_lib} 1

% overflow tank
mat overflow 0.0007 burn 1 vol 1e8
2004.{refuel_lib} 1

% mass flow definitions
mflow U_in
all 0.0

mflow offgasratecore
Ne 1e-2
Ar 1e-2
He 1e-2
Kr 1e-2
Xe 1e-2
Rn 1e-2

% need to account for the increase in vloume with refueling
mflow over
all 0.0

% predictor-corrector must be turned off to use depletion
set pcc 0
% dumps depletion matrices if needed. should be one per burnt material.
% set depmtx 1

%syntax:
% rc <from_mat> <to_mat> <mflow> <setting> where setting is either 0, 1 or 2.

rep source_rep
rc U_stock fuelsalt U_in 0
rc fuelsalt offgastankcore offgasratecore 1
rc fuelsalt overflow over 1
'''
        return repr_cards.format(**locals())

    def get_depl_cards(self) -> str:
        'Depletion data setup'
        depl_cards = '''
% Depletion cards
set inventory all
dep
pro source_rep
daystep
'''
        return depl_cards.format(**locals())

    def get_depl_1st_year(self) -> str:
        '1st year of depletion in daysteps'
        depl_1st_year ='''\
0.0208 0.0208 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
30 30 30 30 30 30 30 % 1 year
'''
        return depl_1st_year.format(**locals())

    def get_depl_add_9years(self) -> str:
        'Add 9 years in daysteps'
        depl_add_9years = '''\
30 30 30 30 30 30 30 30 30 30 30 30
30 30 30 30 30 30 30 30 30 30 30 30
30 30 30 30 30 30 30 30 30 30 30 30
60 60 60 60 60 60 60 60 60 60 60 60
60 60 60 60 60 60 60 60 60 60 60 60
60 60 60 60 60 60 60 60 60 60 60 60
'''
        return depl_add_9years.format(**locals())

    def get_depl_add_10years(self) -> str:
        'Add 10 years in daysteps'
        depl_add_10years = '''\
120 120 120 120 120 120 120 120 120 120 120 120
120 120 120 120 120 120 120 120 120 120 120 120
120 120 120 120 120 120
'''
        return depl_add_10years.format(**locals())

    def get_deck(self) -> str:
        'Serpent deck for the lattice'
        deck = '''\
set title "sphMCFR radius {self.r}, reflector {self.refl}"
'''
        deck += self.get_surfaces()
        deck += self.get_cells()
        deck += "\n"
        deck += self.s.serpent_mat(self.tempK)
        deck += self.get_materials()
        deck += self.get_data_cards()
        if self.deplete > 0.0:
            deck += self.get_repr_cards()
            deck += self.get_depl_cards()
        if self.deplete > 0.5:          # Hacky, but will do :)
            deck += self.get_depl_1st_year()
        if self.deplete > 9:
            deck += self.get_depl_add_9years()
        if self.deplete > 19:
            for i in range((self.deplete - 10) // 10):
                deck += self.get_depl_add_10years()
        return deck.format(**locals())

def liquidusNaClUCl3(xUCl3:float) -> float:
    '''Returns liquidus temperature of NaCl-LiCl3 based on
    https://doi.org/10.1016/j.jnucmat.2015.07.050
    using https://automeris.io/WebPlotDigitizer/
    xUCl3 is UCl3 molar fraction [0:1]
    returns temperature in Celsius'''

    x = [0.00, 0.022808207269860187, 0.05567526666146555, 0.0855318371043732,
         0.11387291763560031, 0.1407171777187291, 0.16456696209224897,
         0.1884138246188261, 0.2107710222408382, 0.23162585386851375,
         0.2524871254853691, 0.2733369482325713, 0.2956844262412841,
         0.32067373092247503, 0.3332737833099288, 0.3525508723644844,
         0.376626682469361, 0.40220160875921956, 0.43226157010619365,
         0.4653201812298379, 0.49837469417854924, 0.5314182786607735,
         0.5667049308894863, 0.5929777611919551, 0.633510191172695,
         0.6665269098414717, 0.6995413517463969, 0.7325526061819301,
         0.7655634052646927, 0.7985687401142121, 0.8315717981998796,
         0.8645721241689257, 0.8975683519630391, 0.9305641244043819,
         0.963555343318022, 0.9875390469412866, 1.0007021363224675]

    T = [798.4941460942108, 790.8183271004245, 774.4729081656745,
         755.0703266232713, 732.0772695140352, 709.2509663838887,
         686.5901650127134, 663.3413419443107, 641.6993259910652,
         619.1080228364567, 597.8127675043092, 574.2134271544534,
         550.6153390247157, 531.2315977031835, 522.0155051642334,
         541.7165932813905, 564.5434703457579, 587.6656104788574,
         609.1954483798368, 631.3998155178843, 652.7794249507294,
         671.9596805030346, 690.4545162978486, 703.3289043893245,
         722.6299707234994, 736.4034813194771, 749.7187931903425,
         762.3926268460502, 774.9748207567357, 786.4573377271508,
         797.4816559724537, 807.9561357476216, 817.6058578175869,
         827.1639401425298, 835.8056250172476, 840.0960726571891,
         841.9715489376551]
    Tliquidus = interpolate.interp1d(x, T, kind='cubic')
    return float(Tliquidus(xUCl3))


# ------------------------------------------------------------
if __name__ == '__main__':
    print("This module handles a simple lattice.")
    input("Press Ctrl+C to quit, or enter else to test it.")
    mycore = MSFR()
    mycore.deplete = 100
    print("***** Serpent deck: \n" + mycore.get_deck() + "\n***** ")
    mycore.deck_path = os.path.expanduser('~/tmp/msfr')
    print(mycore.deck_path + ' / ' + mycore.deck_name)
    mycore.qsub_file = mycore.deck_path + "/run.sh"
#    mycore.ompcores = 8
#    mycore.queue = 'fill'
    mycore.save_deck()
    mycore.save_qsub_file()
    mycore.run_deck()
