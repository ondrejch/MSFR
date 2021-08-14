#!/usr/bin/python3

import math
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import serpentTools
serpentTools.settings.rc['verbosity'] = 'error'

class Resistivity(object):
    '''Class relating to resistivity calcualtions [miloOhm cm]
    http://www.nessengr.com/techdata/metalresis.html'''
    def __init__(self):
        self.ele       = {}
        #    element      rho_0  delta_rho
        self.ele['Ag'] = (  1.6, 0.0038 )   # Silver
        self.ele['Pd'] = ( 10.8, 0.0035 )   # Palladium
        self.ele['Cd'] = (  7.4, 0.0040 )   # Cadmium

    def get_rho(self, element:str, temp:float) -> float:
        'Resistivity [miloOhm cm] as a function of temperature [degC]'
        (rho0, alpha) =self.ele[element]
        return rho0 * (1.0 + alpha*(temp - 20.0))

    def combine_rho(self, rho_c:float, rho_d:float, frac_d:float) -> float:
        '''Calculate resistivity of mixtures
        http://elektroarsenal.net/resistivity-of-mixtures-and-porous-materials.html
        http://anis.buet.ac.bd/Teaching/properties_of_materials/Lecture7.pdf
            rho_c   resistivity of continuos phase
            rho_d   resistivity of dispersed phase
            frac_d  fraction of dispersed phase'''
        if rho_d > 10.0*rho_c:
            return rho_c * (1.0 + 0.5*frac_d) / (1.0 - frac_d)
        if rho_d < 0.1*rho_c:
            return rho_c * (1.0 - frac_d)     / (1.0 + 2.0*frac_d)
        return rho_d*frac_d + rho_c*(1.0-frac_d)

    def get_rho_AgPdCd(self, fAg:float, fPd:float, fCd:float, temp:float) -> float:
        '''Calculates resistivity of Ag-Pd-Cd mixture using the series rule, if
        the series rule applies.'''
        rhoAg = self.get_rho('Ag', temp)
        rhoPd = self.get_rho('Pd', temp)
        rhoCd = self.get_rho('Cd', temp)
        if  rhoAg > 10.0*rhoPd or rhoAg < 0.1*rhoPd or \
            rhoAg > 10.0*rhoCd or rhoAg < 0.1*rhoCd:
            raise ValueError('Resistivities are too different to use mixture rule: Ag ', \
                rhoAg, ' Pd ', rhoPd, ' Cd ' , rhoCd)
        else:
            #print ( fAg, rhoAg,  fPd, rhoPd , fCd, rhoCd)
            return fAg*rhoAg + fPd*rhoPd + fCd*rhoCd



class AgWireAnalyzer(object):
    'Silver wire in depleted salt analysis class'
    def __init__(self, _deckname:str = '/full/path/to/msfr'):
        self.d0   = serpentTools.read(_deckname + '_dep.m')
        self.fuel = self.d0.materials['fuelsalt']
        self.wires = []
        self.wdeps = []
        self.wdeck_path:str = '/tmp'        # Path to wire depletion decks
        self.wdeck_name:str = 'wire_step'   # Wire depletion steps base name
        self.Ntopisos:int   = 10    # How many isotopes to plot
        self.topisos        = []    # List of the top EOC isotopes
        self.agtot          = []    # Total Ag adens
        self.agfrac         = []    # Fraction Ag adens
        self.adata = np.zeros(len(self.fuel.days)*self.Ntopisos). \
            reshape(len(self.fuel.days), self.Ntopisos)

    def read_wires(self):
        'Read all wire input decks'
        for step in range(1,len(self.fuel.days)):
            fname = f'{self.wdeck_path}/{self.wdeck_name}-{step:03d}_dep.m'
            #print(fname)
            d = serpentTools.read(fname)
            w = d.materials['silver']
            self.wdeps.append(d)
            self.wires.append(w)

        'Get silver fraction with depletion'
        Ag_isotopes = []
        for iso in w.names:
            if 'Ag' in iso:  # this is a silver isotope
                Ag_isotopes.append(iso)
        # Initial concentrations
        agsum:float = 0.0
        adens:float = 0.0
        for iso in self.wires[0].names:
            if iso == 'total':
                adens  = self.wires[0].getValues('days', 'adens', [self.wires[0].days[0]], [iso])[0,0]
            if 'Ag' in iso:     # Sum Ag isotopes
                agsum += self.wires[0].getValues('days', 'adens', [self.wires[0].days[0]], [iso])[0,0]
        self.agtot.append(agsum)
        self.agfrac.append(agsum/adens)
        # Rest of burnup
        for step in range(len(self.wires)):
            agsum:float = 0.0
            adens:float = 0.0
            for iso in self.wires[step].names:
                if iso == 'total':
                    adens  = self.wires[step].getValues('days', 'adens', [self.wires[step].days[-1]], [iso])[0,0]
                if 'Ag' in iso:     # Sum Ag isotopes
                    agsum += self.wires[step].getValues('days', 'adens', [self.wires[step].days[-1]], [iso])[0,0]
            self.agtot.append(agsum)
            self.agfrac.append(agsum/adens)

        'Find most abundant Ntopiso isotopes and form adata array for plotting'
        EOCiso = {}
        for iso in w.names:
            if iso == 'total' or iso == 'lost':
                continue
            EOCiso[iso] = w.getValues('days', 'adens', [w.days[-1]], [iso])[0,0]
        # Sort by concentration
        sortedEOCiso = sorted(EOCiso.items(), key=lambda x:x[1], reverse=True)
        # print (sortedEOCiso)
        # First N entries
        self.topisos = [ x[0] for x in sortedEOCiso[0:self.Ntopisos] ]
        # Initial concentrations
        for iso in self.topisos:
            atomdensity = self.wires[0].getValues('days', 'adens', [self.wires[0].days[0]], [iso])[0,0]
            # print(iso, atomdensity)
            self.adata[0, self.topisos.index(iso)] = atomdensity
        # Rest of burnup
        for step in range(len(self.wires)):
            # print('step: '+str(step))
            for iso in self.topisos:
                atomdensity = self.wires[step].getValues('days', 'adens', [self.wires[step].days[-1]], [iso])[0,0]
                # print(step, iso, atomdensity)
                self.adata[step+1, self.topisos.index(iso)] = atomdensity

    def burnup2days(self,x:float) -> float:
        return np.interp(x, self.d0.burnup, self.fuel.days)

    def days2burnup(self,x:float) -> float:
        return np.interp(x, self.fuel.days, self.d0.burnup)

    def plot_agfrac(self, plot_file:str='./plot_wire-Ag-frac.pdf', plot_title = ''):
        'Make plot of silver material evolution with burnup'
        fig    = plt.figure()
        myplot = fig.add_subplot(111)
        myplot.plot(self.d0.burnup, self.agfrac) # TODO
        myplot.set_xlabel('Burnup [MWd/kgHM]')
        myplot.set_ylabel("Silver fraction in the wire")
        #myplot.set_yscale('log')
        #secax = myplot.secondary_xaxis('top', functions=(self.burnup2days, self.days2burnup))
        #secax.set_xlabel('EFPD [days]')
        #myplot.legend(loc="best", fontsize="medium", title="Isotopes in silver wire")
        #(ymin, ymax) = myplot.get_ylim()
        #ymin = 1e-12
        #plt.ylim(ymin, ymax)
        if plot_title != '':
            plt.title(plot_title)
        if plot_file == None:
            plt.show()
        else:
            if not os.path.exists(os.path.dirname(plot_file)):
                os.makedirs(os.path.dirname(plot_file))
            plt.savefig(plot_file, bbox_inches='tight')
            print(f"Saved plot file: {plot_file}")
        plt.close()

    def get_EOCfrac(self, frac_ele='Pt') -> float:
        'Get elemental fraction in wire at EOC'
        frac_elesum:float = 0.0
        for iso in self.wires[-1].names:
            if frac_ele in iso:     # Sum Ag isotopes
                frac_elesum += self.wires[-1].getValues('days','adens',[self.wires[-1].days[-1]],[iso])[0,0]
#                print("ISO :",iso, self.wires[-1].getValues('days','adens',[self.wires[-1].days[-1]],[iso])[0,0]
        frac_elesum /= self.wires[-1].getValues('days','adens',[self.wires[-1].days[-1]],['total'])[0,0]
        return frac_elesum

    def plot_topisos(self, plot_file:str='./plot_wire-Ag-iso.pdf', plot_title = ''):
        'Make plot of isotopic evolution with burnup'
        fig = plt.figure()
        myplot = serpentTools.plot.plot(self.d0.burnup, self.adata, labels=self.topisos)
        myplot.set_xlabel('Burnup [MWd/kgHM]')
        myplot.set_ylabel("Atom density [10$^{24}$/cm$^{3}$]")
        myplot.set_yscale('log')
        secax = myplot.secondary_xaxis('top', functions=(self.burnup2days, self.days2burnup))
        secax.set_xlabel('EFPD [days]')
        myplot.legend(loc="best", fontsize="medium", title="Isotopes in silver wire")
        (ymin, ymax) = myplot.get_ylim()
        ymin = 1e-12
        plt.ylim(ymin, ymax)
        if plot_title != '':
            plt.title(plot_title)
        if plot_file == None:
            plt.show()
        else:
            if not os.path.exists(os.path.dirname(plot_file)):
                os.makedirs(os.path.dirname(plot_file))
            plt.savefig(plot_file, bbox_inches='tight')
            print(f"Saved plot file: {plot_file}")
        plt.close()

'''
import agmsfr
aa = agmsfr.AgWireAnalyzer("/home/o/UTK/research/ARPA-E_Pump_ORNL/silver/msfr_dev/msfr")
aa.wdeck_path="/home/o/UTK/research/ARPA-E_Pump_ORNL/silver/msfr_dev/"
aa.read_wires()
aa.plot_topisos()
'''

class AgMSFRAnalyzer(object):
    'Silver in MSFR shell depletion analysis class'
    def __init__(self, _deckname:str = "msfr"):
        'Path based constructor'
        self.deck_name = _deckname
        self.fr = serpentTools.read(self.deck_name + "_res.m")
        self.fd = serpentTools.read(self.deck_name + "_dep.m")
        self.s  = self.fd.materials['fuelsalt']
        self.s.data['burnup']  = self.fd.metadata['burnup']
        self.ag = self.fd.materials['silver']
        self.ag.data['burnup'] = self.fd.metadata['burnup']

        self.plot_path = "."    # File for plots
        self.agtot     = {}     # Total Ag adens
        self.agfrac    = {}     # Fraction Ag adens

        self.Ntopisos  = 10     # How many isotopes to plot
        self.topisos   = []     # List of the top EOC isotopes

    def get_EOCfrac(self, frac_ele='Pt') -> float:
        'Get elemental fraction at EOC'
        frac_elesum:float = 0.0
        for iso in self.ag.names:
            if frac_ele in iso:     # Sum Ag isotopes
                frac_elesum += self.ag.getValues('days','adens',[self.ag.days[-1]],[iso])[0,0]
#                print("ISO :",iso,  self.ag.getValues('days','adens',[self.ag.days[-1]],[iso])[0,0])
        frac_elesum /= self.ag.getValues('days','adens',[self.ag.days[-1]],['total'])[0,0]
        return frac_elesum

    def calc_agfrac(self):
        'Get silver fraction with depletion'
        for d in self.ag.days:           # For each depletion step
            agsum = 0.0
            for iso in self.ag.names:
                if 'Ag' in iso:     # Sum Ag isotopes
                    agsum += self.ag.getValues('days','adens',[d],[iso])[0,0]
            self.agtot[d]  = agsum
            self.agfrac[d] = agsum / self.ag.getValues('days','adens',[d],['total'])[0,0]

    def plot_agfrac(self, plot_file:str='./plot_Ag-frac.pdf', plot_title = ''):
        'Make plot of Ag fraction remaining in Ag shell with bunrup'
        if len(self.agfrac) < 2:
            self.calc_agfrac()
        fig = plt.plot(self.fd.days, self.agfrac.values(), color='silver', linewidth=2)
        plt.legend(loc="best", fontsize="medium", title="Fraction of silver remaining in Ag shell")
        plt.xlabel("time [d]")
        if plot_title != '':
            plt.title(plot_title)
        if plot_file == None:
            plt.show()
        else:
            if not os.path.exists(os.path.dirname(plot_file)):
                os.makedirs(os.path.dirname(plot_file))
            plt.savefig(plot_file, bbox_inches='tight')
            print(f"Saved plot file: {plot_file}")
        plt.close()

    def calc_topisos(self):
        'Find top N isotopes at EOC for plotting purposes'
        #ag.getValues('days', 'adens', [ag.days[-1]] , ['Xe135'])[0,0]
        # Geth hash of (isotope, concentratio)n at EOC, ag.days[-1]
        EOCiso = {}
        for iso in self.ag.names:
            EOCiso[iso] = self.ag.getValues('days', 'adens', [self.ag.days[-1]], [iso])[0,0]
            #    if 'Ag' in iso:
            #        print(iso)
        # Sort by concentration, highest is 'total'
        sortedEOCiso = sorted(EOCiso.items(), key=lambda x:x[1], reverse=True)
        # First N entries
        self.topisos = [ x[0] for x in sortedEOCiso[1:1+self.Ntopisos] ]

    def plot_topisos(self, plot_file:str='./plot_Ag-iso.pdf', plot_title = ''):
        'Make plot of isotopic evolution with burnup'
        if len(self.topisos) < 2:
            self.calc_topisos()           # Get top isotopes at EOC
        fig = self.fd.plot('burnup','adens', materials=['silver'], names=self.topisos)
        plt.grid(True,which="both")
        if plot_title != '':
            plt.title(plot_title)
        plt.legend(loc="best", fontsize="medium", title="Isotopes in silver")
        plt.yscale('log')
        (ymin, ymax) = fig.get_ylim()
        ymin = 1e-7
        plt.ylim(ymin, ymax)
        if plot_file == None:
            plt.show()
        else:
            if not os.path.exists(os.path.dirname(plot_file)):
                os.makedirs(os.path.dirname(plot_file))
            plt.savefig(plot_file, bbox_inches='tight')
            print(f"Saved plot file: {plot_file}")
        plt.close()

    def plot_multi(self, plot_file:str='./plot.pdf', plot_title = ''):
        'Multiplot'
        if len(self.topisos) < 2:
            self.calc_topisos()
        if len(self.agfrac) < 2:
            self.calc_agfrac()
        plt.figure(figsize=(7,10))
        plt.subplot(211)
        fig = plt.plot(self.fd.days, self.agfrac.values(), color='silver', linewidth=2)
        plt.xlabel("time [d]")
        plt.legend(loc="best", fontsize="medium", title="Fraction of silver remaining in Ag shell")
        if plot_title != '':
            plt.title(plot_title)
        plt.subplot(212)
        fig = self.fd.plot('burnup','adens', materials=['silver'], names=self.topisos)
        plt.grid(True,which="both")
        plt.legend(loc="best", fontsize="medium", title="Isotopes in silver")
        plt.yscale('log')
        (ymin, ymax) = fig.get_ylim()
        ymin = 1e-7
        plt.ylim(ymin, ymax)

        if plot_file == None:
            plt.show()
        else:
            if not os.path.exists(os.path.dirname(plot_file)):
                os.makedirs(os.path.dirname(plot_file))
            plt.savefig(plot_file, bbox_inches='tight')
            print(f"Saved plot file: {plot_file}")
        plt.close()

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    print("This module analyzes MSFR with silver.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    a = AgMSFR("/home/o/tmp/ag/small/run0/msfr")
    a.calc_agfrac()
    a.calc_topisos()
    a.plot_topisos()
