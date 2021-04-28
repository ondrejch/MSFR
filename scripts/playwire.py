#!/usr/bin/python3

import math
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import serpentTools


class PlayWire(object):
    'aluminum depletion analysis class'
    def __init__(self, _deckname:str = "wire"):
        'Path based constructor'
        self.deck_name = _deckname
        self.ag_mat = 'testwire'
        self.fr = serpentTools.read(self.deck_name + "_res.m")
        self.fd = serpentTools.read(self.deck_name + "_dep.m")
        self.s  = self.fd.materials['fuel']
        self.s.data['burnup']  = self.fd.metadata['burnup']
        self.ag = self.fd.materials[self.ag_mat]
        self.ag.data['burnup'] = self.fd.metadata['burnup']

        self.plot_path = "."    # File for plots
        self.agtot     = {}     # Total Al adens
        self.agfrac    = {}     # Fraction Al adens

        self.Ntopisos  = 10     # How many isotopes to plot
        self.topisos   = []     # List of the top EOC isotopes

    def calc_agfrac(self):
        'Get aluminum fraction with depletion'
        for d in self.ag.days:           # For each depletion step
            agsum = 0.0
            for iso in self.ag.names:
                if 'Al' in iso:     # Sum Al isotopes
                    agsum += self.ag.getValues('days','adens',[d],[iso])[0,0]
                self.agtot[d]  = agsum
                self.agfrac[d] = agsum / self.ag.getValues('days','adens',[d],['total'])[0,0]

    def plot_agfrac(self, plot_file:str='./plot_Al-frac.pdf', plot_title = ''):
        'Make plot of Al fraction remaining in the wire'
        if len(self.agfrac) < 2:
            self.calc_agfrac()
        fig = plt.plot(self.fd.days, list(self.agfrac.values()), color='black', linewidth=2)
        plt.legend(loc="best", fontsize="medium", title="Fraction of aluminum remaining in the wire")
        plt.xlabel("time [d]")
        if plot_title != '':
            plt.title(plot_title)
        if plot_file == None:
            plt.show()
        else:
            if not os.path.exists(os.path.dirname(plot_file)):
                os.makedirs(os.path.dirname(plot_file))
            plt.savefig(plot_file, bbox_inches='tight')
        plt.close()

    def calc_topisos(self):
        'Find top N isotopes at EOC for plotting purposes'
        #ag.getValues('days', 'adens', [ag.days[-1]] , ['Xe135'])[0,0]
        # Geth hash of (isotope, concentratio)n at EOC, ag.days[-1]
        EOCiso = {}
        for iso in self.ag.names:
            if iso != 'total':
                EOCiso[iso] = self.ag.getValues('days', 'adens', [self.ag.days[-1]], [iso])[0,0]
            #    if 'Al' in iso:
            #        print(iso)
        # Sort by concentration, highest is 'total'
        sortedEOCiso = sorted(EOCiso.items(), key=lambda x:x[1], reverse=True)
        # First N entries
        self.topisos = [ x[0] for x in sortedEOCiso[0:self.Ntopisos-1] ]

    def plot_topisos(self, plot_file:str='./plot_Al-iso.pdf', plot_title = ''):
        'Make plot of isotopic evolution with burnup'
        if len(self.topisos) < 2:
            self.calc_topisos()           # Get top isotopes at EOC
        #fig = self.fd.plot('burnup','adens', materials=[self.ag_mat], names=self.topisos)
        fig = self.fd.plot('days','adens', materials=[self.ag_mat], names=self.topisos)
        plt.grid(True,which="both")
        if plot_title != '':
            plt.title(plot_title)
        plt.legend(loc="best", fontsize="medium", title="Isotopes in aluminum")
        plt.yscale('log')
        (ymin, ymax) = fig.get_ylim()
        ymin = 1e-14
        plt.ylim(ymin, ymax)
        if plot_file == None:
            plt.show()
        else:
            if not os.path.exists(os.path.dirname(plot_file)):
                os.makedirs(os.path.dirname(plot_file))
            plt.savefig(plot_file, bbox_inches='tight')
        plt.close()

    def plot_multi(self, plot_file:str='./plot.pdf', plot_title = ''):
        'Multiplot'
        if len(self.topisos) < 2:
            self.calc_topisos()
        if len(self.agfrac) < 2:
            self.calc_agfrac()
        plt.figure(figsize=(7,10))
        plt.subplot(211)
        fig = plt.plot(self.fd.days, list(self.agfrac.values()), color='aluminum', linewidth=2)
        plt.xlabel("time [d]")
        plt.legend(loc="best", fontsize="medium", title="Fraction of aluminum remaining in the wire")
        if plot_title != '':
            plt.title(plot_title)
        plt.subplot(212)
        fig = self.fd.plot('burnup','adens', materials=[self.ag_mat], names=self.topisos)
        plt.grid(True,which="both")
        plt.legend(loc="best", fontsize="medium", title="Isotopes in aluminum")
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
        plt.close()

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    print("This module analyzes aluminum depletion from decay source.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    a = PlayWire()
    a.calc_agfrac()
    a.calc_topisos()
    a.plot_topisos()
    a.plot_agfrac()



