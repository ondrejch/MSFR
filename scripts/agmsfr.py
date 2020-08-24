#!/usr/bin/python3

import math
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import serpentTools


class AgMSFR(object):
    'Silver depletion analysis class'
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

    def calc_agfrac(self):
        'Get silver fraction with depletion'
        for d in self.ag.days:           # For each depletion step
            agsum = 0.0
            for iso in self.ag.names:
                if 'Ag' in iso:     # Sum Ag isotopes
                    agsum += self.ag.getValues('days','adens',[d],[iso])[0,0]
                self.agtot[d]  = agsum
                self.agfrac[d] = agsum / self.ag.getValues('days','adens',[d],['total'])[0,0]

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

    def plot_topisos(self,plot_file:str='plot.pdf', plot_title = ''):
        'Make plot ot isotopoic evolution with burnup'
        if len(self.topisos) < 2:
            self.topiso()           # Get top isotopes at EOC
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
            my_file = self.plot_path + '/' + plot_file
            if not os.path.exists(os.path.dirname(my_file)):
                os.makedirs(os.path.dirname(my_file))
            plt.savefig(my_file, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    print("This module analyzes MSFR with silver.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    a = AgMSFR("/home/o/tmp/ag/small/run0/msfr")
    a.calc_agfrac()
    a.calc_topisos()
    a.plot_topisos()



