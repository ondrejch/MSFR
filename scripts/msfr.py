#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-13
# GNU/GPL

'''
This file contains classes that generate imput decks for MSFR Serpent calcaltions.

MSFR(MSFRbase) is a basic spherical MSFR, which can have a silver shell embedded in a reflector

AgWire(MSFRbase) describes trasmutation of a silver wire in a a depleted MSFR salt. The key issue is that the irraditing salt (the radiation source) has to be built manually for each depletion step. This is done at the end of method AgWire.wired_deck().

liquidusNaClUCl3 is a helper function to obtain liquidus temperature of NaCl-UCl3 salt.
'''

import os
import math
import re
from scipy import interpolate
from textwrap import dedent
from salts import Salt
import serpentTools

do_plots = True
my_debug = False

NUCLEAR_LIBRARIES = ['endf7','jeff33','endf8']


class MSFRbase(object):
    '''Common base class for the MSFR project'''
    def __init__(self):
        self.tempK:float   = 900.0      # Salt temperature [K]
        self.silver_T:float= self.tempK + 10   # Temperature of the silver wire [K]
        # The wire is assumed hotter as the salt is the coolant
        self.power:float   = 3e9        # Core thermal power [W]
        self.deplete:float = 0          # Depletion flag, hacky, see code!
                                        # 0 - no depletion, then in years
                                        # Can be 0, 1, 10, 20, 30, 40 ... n/10 years
        self.lib:str       = '09c'      # CE xsection temp selection
        self.lib_ag:str    = '09c'      # CE xsection temp selection for silver
        self.nfg:str       = None       # Energy group structure to use is any
        self.queue:str     = 'gen6'     # NEcluster torque queue
        self.histories:int = 50000      # Neutron histories per cycle
        self.ompcores:int  = 16         # OMP core count
        self.deck_name:str = 'mcfr_input'     # Serpent input file name
        self.deck_path:str = '/tmp'     # Where to run the Serpent deck
        self.nuc_libs:str  = 'jeff33'   # Nuclear data libraries
        self.qsub_file:str = os.path.expanduser('~/') + '/run.sh'  # qsub script path

    def rho_silver(self) -> float:
        # https://www.sciencedirect.com/science/article/abs/pii/0022190262801882
        return 10.465 - 9.967e-4*self.silver_T # [g/cm^3]

    def matdeck_silver(self, mat_name:str='silver', burn:int=1) -> str:
        '''Returns material cards for silver'''
        rgb:str ="110 110 110"
        if burn:
            rgb:str ="210 210 210"
        return f'''
% Silver
mat {mat_name} -{self.rho_silver()} tmp {self.silver_T} rgb {rgb} burn {burn}
47107.{self.lib_ag}  -0.51839    % Ag
47109.{self.lib_ag}  -0.48161    % Ag
'''

    def lib_deck(self) -> str:
        '''Returns cards for nuclear data libraries'''
        if self.nuc_libs in NUCLEAR_LIBRARIES:
            pass
        else:
            raise ValueError("ERROR: Nuclear data library ",nuclib," is unknown.")
        if self.nuc_libs == 'jeff33':
            return '''
% Data Libraries
set acelib "/opt/JEFF-3.3/sss_jeff33.xsdir"
set declib "/opt/JEFF-3.3/jeff33.dec"
set nfylib "/opt/JEFF-3.3/jeff33.nfy"
'''
        if self.nuc_libs == 'endf7':
            return '''
% Data Libraries
set acelib "sss_endfb7u.sssdir"
set declib "sss_endfb7.dec"
set nfylib "sss_endfb7.nfy"
'''
        if self.nuc_libs == 'endf8':
            return '''
% Data Libraries
set acelib "/opt/ENDFB-8.0/endfb80.xsdir"
set declib "/opt/ENDFB-8.0/sss_endfb80.dec"
set nfylib "/opt/ENDFB-8.0/sss_endfb80.nfy"
'''

    def run_deck(self):
        'Runs the deck using qsub_file script'
        if self.queue == 'local':    # Run the deck locally
            os.chdir(self.deck_path)
            os.system(self.qsub_file)
        else:               # Submit the job on the cluster
            os.system('cd ' + self.deck_path + ' && qsub ' + self.qsub_file)


AGWIRE_CASES = ['fully-submerged', 'half-submerged']

class AgWire(MSFRbase):
    '''Silver wire depleted in MSFR fuel salt. Usage:
# Example class usage:
import msfr
w = msfr.AgWire(0.2, 'half-submerged')
w.deck_path='/home/ondrejch/APump/wire_small_jeff33/130/hs/'
w.load_data()
w.save_decks()
w.save_qsub_file()     '''
    def __init__(self, wr:float = 0.2, case:str='fully-submerged'):
        if case in AGWIRE_CASES:
            self.case = case
        else:
            raise ValueError('Wrong case' + case)
        self.wr:float = wr      # wire radius [cm]
        self.fh:float = 100.0   # half-length of salt cylinder [cm]
        self.fr:float = 2.0     # radius of the salt cylinder [cm]
        if(self.fr < 2.0*self.wr):
            self.fr = 2.0*self.wr   # salt cylinder needs to be at least 2x the wire one
        MSFRbase.__init__(self) # in Python the parent class needs explicit initialization
        self.dep  = None        # depletion object from the baseline depletion case
        self.fuel = None        # fuel object from the baseline depletion
        self.wdeck_name:str = 'wire_step' # name of input deck for wire depletion steps
        self.qsub_file:str = os.path.expanduser('~/') + '/runwire.sh' # qsub script path

    def load_data(self):
        '''Open the depletion file and load the fuel data. Make sure that data path and names are
        set correctly in the parent class'''
        self.dep = serpentTools.read(self.deck_path + '/' + self.deck_name + '_dep.m')
        self.fuel = self.dep.materials['fuelsalt']

    def volume_wire(self) -> float:
        '''Calculates the wire volume'''
        return math.pi * self.wr**2 * 2.0*self.fh

    def volume_fuel(self) -> float:
        '''Calculates the fuel salt volume'''
        V = math.pi * 2.0*self.fh * (self.fr**2 - self.wr**2)
        if self.case == 'fully-submerged':
            return V
        if self.case == 'half-submerged':
            return V / 2.0

    def wire_deck(self, step:int=1) -> str:
        '''Returns wire-in-salt Serpent input deck for a particular burnup step calculation'''
        if(step < 1):
            return 'Error: step has to be >= 1, value passed: ' + str(step)
        prevstep = step - 1
        day      = self.fuel.days[step]
        prevday  = self.fuel.days[prevstep]
        output = 'set title "Activated wire in decaying fuel"\n'
        if self.case == 'fully-submerged':
            output += f'''
% --- surfaces ---
surf 1   cylx  0.0 0.0 {self.wr} -{self.fh} {self.fh}    % inner wire
surf 2   cylx  0.0 0.0 {self.fr} -{self.fh} {self.fh}    % fuel cylinder

% --- cells ---
cell 10  0  silver  -1      % wire
cell 11  0  fuel     1 -2   % fuel salt
cell 99  0  outside  2      % graveyard
'''
        if self.case == 'half-submerged':
            output += f'''
% --- surfaces ---
surf 1   cylx  0.0 0.0 {self.wr} -{self.fh} {self.fh}    % inner wire
surf 2   cylx  0.0 0.0 {self.fr} -{self.fh} {self.fh}    % fuel cylinder
surf 3   pz    0

% --- cells ---
cell 10  0  silver  -1        % wire
cell 11  0  fuel     1 -2  3  % fuel salt
cell 12  0  r-silver 1 -2 -3  % reflector silver, nondepleting
cell 99  0  outside  2        % graveyard
'''
        output += self.matdeck_silver()         # silver wire
        if self.case == 'half-submerged':       # surrounding silver
            output += self.matdeck_silver('r-silver',0)
        if self.case == 'fully-submerged':
            output += f'''
% Volumes
set mvol fuel   0  {self.volume_fuel()}
set mvol silver 0  {self.volume_wire()}\n'''
        if self.case == 'half-submerged':
            output += f'''
% Volumes
set mvol fuel     0  {self.volume_fuel()}
set mvol silver   0  {self.volume_wire()}
set mvol r-silver 0  {self.volume_fuel()}
'''
        output += self.lib_deck()

        if self.nfg is not None:
            output += f'''
% Use group structure for group constant generation
set micro {self.nfg}
set nfg {self.nfg}
'''
        else:
            output += f'''
% Turn off group constant generation
set gcu -1
'''

        output += f'''
% Depletion
set inventory all
dep daytot {day}

% Flux spectrum
det flux de fluxgrid dm silver
ene fluxgrid 3 500 1e-11 2e1

% Read binary restart file
set rfw 1'''
        if step > 1:
            output += f'''
set rfr -{prevday} "wire_step-{prevstep:03d}.wrk"'''

        output += f'''

% Radioactive decay source:
src 1 n sg fuel 1

% Options:
set nps 100000000

% --- materials ---
mat fuel sum fix "{self.lib}" {self.tempK} rgb 50 210 50
'''
        #
        # Write material composition for the burned salt fuel
        # (this acts as a neutron source for the simulation)
        #
        iso_has_xs:bool = True
        prevzai:int = 0
        m_offset:int = 400  # isomer offset for ZA.id
        for zai in self.fuel.zai:
            if zai == 0:    # total
                continue
            if zai == 666:  # lost
                continue
            if zai < prevzai:   # once the ZADs stop increasing, nuclides without cross sections follow
                iso_has_xs = False
            atomdensity = float(self.fuel.getValues('days', 'adens', [day], zai=zai))
            prevzai = zai
            if iso_has_xs:  # isotopes with xs data: <ZZAAA with isome offset> . library
                isoID = str(zai//10 + m_offset*(zai%10)) + "." + self.lib
            else:           # isotopes without xs data: ZAI
                isoID = str(zai)
            if atomdensity: # skip 0 atom densities
                output += f'{isoID}    {atomdensity}\n'
        return output

    def save_decks(self):
        '''Writes input wire depletion to respective files'''
        for step in range(1, len(self.fuel.days)):
            fname = f'{self.deck_path}/{self.wdeck_name}-{step:03d}'
            try:                # Write the deck
                f = open(fname, 'w')
                f.write(self.wire_deck(step))
                f.close()
            except IOError as e:
                print("Unable to write to file", fname)
                print(e)

    def save_qsub_file(self):
        '''Writes a qsub job submission file to run all steps.
        The depletion steps have to be run consecutively.'''
        try:                # Write the script
            frun = open(self.qsub_file, 'w')
            frun.write(f'''#!/bin/bash
#PBS -V
#PBS -N S2-wire
#PBS -q {self.queue}
#PBS -l nodes=1:ppn={self.ompcores}

hostname
rm -f donewire.dat
cd ${{PBS_O_WORKDIR}}
module load mpi
module load serpent
''')
        except IOError as e:
            print("Unable to write to file", fname)
            print(e)
        for step in range(1, len(self.fuel.days)):
            frun.write(f'''
sss2 -omp {self.ompcores} {self.wdeck_name}-{step:03d} > myout_{step:03d}.out''')
        frun.write('\n')
        frun.close()


class MSFR(MSFRbase):
    '''Molten Spherical chloride salt Fast Reactor
# Example class usage:
import msfr
mycore = msfr.MSFR(200, 300, 0.1975, "66.66%NaCl+33.34%UCl3")
mycore.power   = 1e9 # power 1GWth
mycore.deplete = 10  # 10 years
mycore.deck_path = '/tmp/'
mycore.save_deck()
    '''
    def __init__(self, r:float=300.0, refl:float=500.0, e:float=0.1083, salt="58%NaCl+42%UCl3", Ag_r:float=-1.0):
        if r<10.0 or refl<r or e>1.0 or e<0.0:  # Reject bad input
            raise ValueError("Bad parameters: ", r, refl, e)

        MSFRbase.__init__(self)         # in Python the parent class needs explicit initialization
        # core parameters
        self.r:float           = r      # Core radius [cm]
        self.refl:float        = refl   # Outer reflector thickness [cm]
        self.refl_lib          = '06c'  # Reflector nuclear data library
        self.refl_tempK        = 873.0  # Reflector temperature [K]
        self.salt_formula:str  = salt   # Salt formula
        self.refuel_flow:float = 0.0    # wt_fraction/s refuel flow
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
        refl_lib  = self.refl_lib
        materials = '''
% Cast iron reflector
mat refl   -7.034 tmp {self.refl_tempK} rgb 128 128 178
  6000.{refl_lib} -0.034000
%14000.{refl_lib} -0.026000
 14028.{refl_lib}  -2.38853E-02
 14029.{refl_lib}  -1.25674E-03
 14030.{refl_lib}  -8.57970E-04
 15031.{refl_lib} -0.003000
%16000.{refl_lib} -0.001000
 16032.{refl_lib}  -9.47153E-04
 16033.{refl_lib}  -7.71207E-06
 16034.{refl_lib}  -4.50224E-05
 16036.{refl_lib}  -1.12170E-07
 25055.{refl_lib} -0.006500
%26000.{refl_lib} -0.929500
 26054.{refl_lib}  -5.24755E-02
 26056.{refl_lib}  -8.54225E-01
 26057.{refl_lib}  -2.00806E-02
 26058.{refl_lib}  -2.71920E-03
%
%pure iron
%26054.{refl_lib}  -0.058450   %  Fe
%26056.{refl_lib}  -0.917540   %  Fe
%26057.{refl_lib}  -0.021190   %  Fe
%26058.{refl_lib}  -0.002820   %  Fe
'''

        if self.silver_at_r > 0.0 and self.silver_at_r <= self.r:
            raise ValueError('Silver shell inside fuel ', self.silver_at_r, self.r)
        if self.silver_at_r >= self.refl:
            raise ValueError('Silver shell outside reflector ', self.silver_at_r, self.refl)
        if self.silver_at_r > self.r and self.silver_at_r < self.refl:
            materials += self.matdeck_silver()

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

'''
        if self.silver_at_r > self.r and self.silver_at_r < self.refl:
            data_cards += '''
% Flux in silver shell
det silverflux de fluxgrid dm silver
ene fluxgrid 3 500 1e-11 2e1

'''
        if self.nfg is not None:
            data_cards += f'''
% Use group structure for group constant generation
set micro {self.nfg}
set nfg {self.nfg}
'''
        else:
            data_cards += f'''
% Turning off group constant generation hastens the calculation
set gcu -1
'''
        data_cards += self.lib_deck()

        if do_plots:
            data_cards += '''
% Plots
plot 3 1500 1500
% mesh 3 1500 1500
'''
        return data_cards.format(**locals())

    def get_repr_cards(self) -> str:
        'Reprocessing setup'
        refuel_lib = self.lib   # First, build refuel stream using the same material as fuel
        refuel_tmp   = self.s.serpent_mat(self.tempK) # Same as fresh fuel
        refuel_split = refuel_tmp.split('\n')
        refuel_rho   = re.search(r'-[0-9.]+', refuel_split[1]).group()
        if float(refuel_rho) > -1.0:     # Sanity check
            raise ValueError('Refuel density problem, ',refuel_rho)
        refuel_volume = self.salt_volume()
        refuel  = 'mat U_stock {refuel_rho} burn 1 vol {refuel_volume} tmp {self.tempK}\n'.format(**locals())
        refuel += '\n'.join(refuel_split[2:])   # Add isotopic density list
        repr_cards = '''
%___________Reprocessing___________
% First we need some extra materials to do depletion with reprocessing correctly.

% stockpile of extra refuel
{refuel}

% tanks for offgases
mat offgastankcore 0.0007 burn 1 vol 1e6 tmp {self.tempK}
2004.{refuel_lib} 1

% overflow tank
mat overflow 0.0007 burn 1 vol 1e8 tmp {self.tempK}
2004.{refuel_lib} 1

% mass flow definitions
mflow U_in
all {self.refuel_flow}

mflow offgasratecore
Ne 1e-2
Ar 1e-2
He 1e-2
Kr 1e-2
Xe 1e-2
Rn 1e-2

% need to account for the increase in volume with refueling
mflow over
all {self.refuel_flow}

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
0.05 0.15 0.3 0.5   % 1 day
1 2 3               % 1 week
7 7 7 14 14 14 14 28 28 28 28 42 42 42 44  % 1 year, 366 days
'''
        return depl_1st_year.format(**locals())

    def get_depl_add_9years(self) -> str:
        'Add 9 years in daysteps'
        depl_add_9years = '''\
52 52 52 52 52 52 53    % 365
52 52 52 52 52 52 53    % 365
52 52 52 52 52 52 54    % 366
52 52 52 52 52 52 53    % 365
52 52 52 52 52 52 53    % 365
52 52 52 52 52 52 54    % 366
52 52 52 52 52 52 53    % 365
52 52 52 52 52 52 53    % 365
52 52 52 52 52 52 53    % 365
'''
        return depl_add_9years.format(**locals())

    def get_depl_add_10years(self) -> str:
        'Add 10 years in daysteps'
        depl_add_10years = '''\
120 120 126 120 120 125 120 120 125 120 120 125
120 120 126 120 120 125 120 120 125 120 120 125
120 120 126 120 120 125
'''
        return depl_add_10years.format(**locals())

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


class MCRE(MSFRbase):
    '''Molten cylindrical Chloride salt fast Reactor Experiment
    MCRE usage:
        import msfr
        mycore = msfr.MCRE(20, 90, 35, 0.9, "66.66%NaCl+33.34%UCl3", 'MCRE')
        mycore.power   = 3e5 # power 300 kWth
        mycore.deck_path = '/tmp/'
        mycore.save_deck()
    MCFR usage:
        import msfr
        mycore = msfr.MCRE(200, 200, 200, 0.05, "66.66%NaCl+33.34%UCl3", 'MCFR')
        mycore.power   = 1.8e9 # power 1.8 GWth
        mycore.deck_path = '/tmp/'
        mycore.save_deck()
    '''
    def __init__(self, r:float=20.0, h:float=90, refl:float=35.0, e:float=0.9, salt="66.66%NaCl+33.34%UCl3", design='MCRE'):
        if r<10.0 or refl<r or e>1.0 or e<0.0:  # Reject bad input
            raise ValueError("Bad parameters: ", r, e)
        elif design!='MCRE' and design!='MCFR':
            raise ValueError("Bad parameter: ", design)

        MSFRbase.__init__(self)
        # core parameters
        self.r:float           = r      # Core radius [cm]
        self.h:float           = h      # Core height [cm]
        self.refl:float        = refl-r # Reflector thickness [cm]
        self.refl_lib          = '06c'  # Reflector nuclear data library
        self.refl_tempK        = 873.0  # Reflector temperature [K]
        self.salt_formula:str  = salt   # Salt formula
        self.refuel_flow:float = 0.0    # wt_fraction/s refuel flow
        self.s             = Salt(self.salt_formula, e) # Salt used
        if design == 'MCRE':
            self.s.set_chlorine_37Cl_fraction(0.24)      # Natural chlorine-37
        elif design == 'MCFR':
            self.s.set_chlorine_37Cl_fraction(0.90)      # Enriched chlorine-37
        self.design            = design  # MCRE or MCFR

    def salt_volume(self) -> float:
        '''Get salt volume, twice the fuel cylinder volume'''
        if self.design == 'MCRE':
            b_cone = math.pi * (self.r/4)**2 * (2/5*self.refl + 3) / 3
        elif self.design == 'MCFR':
            b_cone = math.pi * (self.r/2)**2 * (3/5*self.refl + 3) / 3
        cylinder = math.pi * self.r**2 * self.h
        t_cone = math.pi * (self.r/2)**2 * self.h/20 / 3
        V = cylinder - t_cone - b_cone
        return 2.0 * V

    def get_cells(self) -> str:
        'Cell cards for Serpent input deck'
        cells = dedent('''
            %______________cell definitions_____________________________________
            cell 30  0  refl       1 -2 -3  4        % radial reflector
            cell 31  0  refl      -2  3 -5           % upper reflector
            cell 32  0  refl      -8 -3              % upper reflector cone
            cell 33  0  refl      -2 -4  6           % lower reflector
            cell 34  0  refl      -7 4 -9            % lower reflector cone
            cell 50  0  fuelsalt  -1 -3  4 #32 #34   % fuel salt
            cell 97  0  outside    2                 % outside
            cell 98  0  outside    5                 % outside
            cell 99  0  outside    -6                % outside
            ''')
        return cells.format(**locals())

    def get_surfaces(self) -> str:
        'Surface cards for Serpent input deck'
        refl_r = self.r + self.refl  # reflector radius
        refl_top, refl_bottom = self.h+self.refl, 0.0-self.refl
        tcone_r = self.r/2
        if self.design == 'MCRE':
            cutoff = 2/5 * self.refl
            bcone_r = self.r/4
        elif self.design == 'MCFR':
            cutoff = 3/5 * self.refl
            bcone_r = self.r/2
        tcone_h, bcone_h = -self.h/20, cutoff+3
        surfaces = dedent('''
            %______________surface definitions__________________________________
            surf 1  cylz  0.0 0.0 {self.r}       % fuel salt
            surf 2  cylz  0.0 0.0 {refl_r}       % radial reflector
            surf 3  pz    {self.h}              % fuel top
            surf 4  pz    0                  % fuel bottom
            surf 5  pz    {refl_top}              % refl top
            surf 6  pz    {refl_bottom}              % refl bottom
            surf 7 cone   0 0 0 {bcone_r} {bcone_h}     % bottom refl cone, x y z r h
            surf 8 cone   0 0 {self.h} {tcone_r} {tcone_h}      % top refl cone
            surf 9 pz     {cutoff}                    % bottom cone refl cutoff
            ''')
        return surfaces.format(**locals())

    def get_materials(self) -> str:
        'Material definitions, non-salt'
        refl_lib  = self.refl_lib
        if self.design == 'MCRE':
            materials = dedent('''
                % MgO reflector
                mat refl -3.5 tmp 873.0 rgb 75 75 75
                 12024.06c 1.0
                 8016.06c 1.0
                ''')
        elif self.design == 'MCFR':
            materials = dedent('''
                % Lead reflector
                mat refl -10.4 tmp 873.0 rgb 75 75 75
                 82204.{refl_lib} 0.014
                 82206.{refl_lib} 0.241
                 82207.{refl_lib} 0.221
                 82208.{refl_lib} 0.524
                ''')
        return materials.format(**locals())

    def get_data_cards(self) -> str:
        'Data cards for the reactor'
        fs_volume = self.salt_volume()
        data_cards = dedent('''
            set mvol fuelsalt 0 {fs_volume}  % Fuel salt volume

            set bc 1  % Boundary condition, vacuum

            % set arr 2  % Analog reaction rate

            set pop {self.histories} 240 40  % N pop and criticality cycles
            ''')
        if self.design == 'MCRE':
            data_cards += dedent(f'''
                set power 300000.0  % Power, 300 thermal kW
                ''')
        elif self.design == 'MCFR':
            data_cards += dedent(f'''
                set power 1800000000.0  % Power, 1.8 thermal GW
                ''')
        if self.nfg is not None:
            data_cards += dedent(f'''
                % Use group structure for group constant generation
                set micro {self.nfg}
                set nfg {self.nfg}
                ''')
        else:
            data_cards += dedent(f'''
                set gcu -1  % Turning off group constant generation hastens the calculation
                ''')
        data_cards += self.lib_deck()

        if do_plots:
            data_cards += dedent('''
                % Plots
                plot 3 1500 1500
                plot 2 1500 1500
                ''')
        return data_cards.format(**locals())

    def get_repr_cards(self) -> str:
        'Reprocessing setup'
        refuel_lib = self.lib   # First, build refuel stream using the same material as fuel
        refuel_tmp   = self.s.serpent_mat(self.tempK) # Same as fresh fuel
        refuel_split = refuel_tmp.split('\n')
        refuel_rho   = re.search(r'-[0-9.]+', refuel_split[1]).group()
        if float(refuel_rho) > -1.0:     # Sanity check
            raise ValueError('Refuel density problem, ',refuel_rho)
        refuel  = 'mat U_stock {refuel_rho} burn 1 vol 1e8 tmp {self.tempK}\n'.format(**locals())
        refuel += '\n'.join(refuel_split[2:])   # Add isotopic density list
        repr_cards = dedent('''
            %___________Reprocessing___________
            % First we need some extra materials to do depletion with reprocessing correctly.

            {refuel}  % stockpile of extra refuel

            % tanks for offgases
            mat offgastankcore 0.0007 burn 1 vol 1e6 tmp {self.tempK}
            2004.{refuel_lib} 1

            % overflow tank
            mat overflow 0.0007 burn 1 vol 1e8 tmp {self.tempK}
            2004.{refuel_lib} 1

            % mass flow definitions
            mflow U_in
            all {self.refuel_flow}

            mflow offgasratecore
            Ne 1e-2
            Ar 1e-2
            He 1e-2
            Kr 1e-2
            Xe 1e-2
            Rn 1e-2

            % need to account for the increase in volume with refueling
            mflow over
            all {self.refuel_flow}

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
            ''')
        return repr_cards.format(**locals())

    def get_depl_cards(self) -> str:
        'Depletion data setup'
        depl_cards = dedent('''
            % Depletion cards
            set inventory all
            dep
            pro source_rep
            daystep
            ''')
        return depl_cards.format(**locals())

    def get_depl_1st_year(self) -> str:
        '1st year of depletion in daysteps'
        depl_1st_year = dedent('''\
            0.05 0.15 0.3 0.5   % 1 day
            1 2 3               % 1 week
            7 7 7 14 14 14 14 28 28 28 28 42 42 42 44  % 1 year, 366 days''')
        return depl_1st_year.format(**locals())

    def get_depl_add_9years(self) -> str:
        'Add 9 years in daysteps'
        depl_add_9years = dedent('''\
            52 52 52 52 52 52 53    % 365
            52 52 52 52 52 52 53    % 365
            52 52 52 52 52 52 54    % 366
            52 52 52 52 52 52 53    % 365
            52 52 52 52 52 52 53    % 365
            52 52 52 52 52 52 54    % 366
            52 52 52 52 52 52 53    % 365
            52 52 52 52 52 52 53    % 365
            52 52 52 52 52 52 53    % 365''')
        return depl_add_9years.format(**locals())

    def get_depl_add_10years(self) -> str:
        'Add 10 years in daysteps'
        depl_add_10years = dedent('''\
            120 120 126 120 120 125 120 120 125 120 120 125
            120 120 126 120 120 125 120 120 125 120 120 125
            120 120 126 120 120 125''')
        return depl_add_10years.format(**locals())

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
        qsub_content = dedent('''#!/bin/bash
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
            ''').format(**locals())
        try:                # Write the deck
            f = open(self.qsub_file, 'w')
            f.write(qsub_content)
            f.close()
        except IOError as e:
            print("Unable to write to qsub file", f)
            print(e)

    def get_deck(self) -> str:
        'Serpent deck for the lattice'
        deck = dedent('''\
            set title "cylMCFR radius {self.r}, height {self.h}, reflector {self.refl}" ''')
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
    print("This module handles a MSFR deck generation class.")
    input("Press Ctrl+C to quit, or enter else to test it.")
    mycore = MCRE(r=200.0, h=450.0, refl=250.0, e=0.13, salt="66.66%NaCl+33.34%UCl3", design='MCFR')  # MSFR()
    mycore.deplete = 0  # 100
    print("***** Serpent deck: \n" + mycore.get_deck() + "\n***** ")
    mycore.deck_path = os.path.expanduser('/Users/klawso28/Documents')  # ('~/tmp/msfr')
    print(mycore.deck_path + '/' + mycore.deck_name)
    mycore.qsub_file = mycore.deck_path + "/run.sh"
#    mycore.ompcores = 8
#    mycore.queue = 'fill'
    mycore.save_deck()
    mycore.save_qsub_file()
#   mycore.run_deck()

'''
# Example code usage:

import msfr
w = msfr.AgWire()
w.deck_path='/home/ondrejch/APump/MCFR/ag/jeff33-wire/run0'
w.load_data()
w.save_decks()
w.save_qsub_file()


import msfr
w = msfr.AgWire(0.2, 'half-submerged')
w.deck_path='/home/ondrejch/APump/MCFR/ag/jeff33-wire/run0/hs'
w.load_data()
w.save_decks()
w.save_qsub_file()


import msfr
mycore = msfr.MCRE(20, 90, 35, 0.9, "66.66%NaCl+33.34%UCl3", 'MCRE')
mycore.power   = 3e5 # power 300 kWth
mycore.deck_path = '/tmp/'
mycore.save_deck()

import msfr
mycore = msfr.MCRE(200, 200, 200, 0.05, "66.66%NaCl+33.34%UCl3", 'MCFR')
mycore.power   = 3e5 # power 300 kWth
mycore.deck_path = '/tmp/'
mycore.save_deck()

'''
