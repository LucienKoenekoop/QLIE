import argparse
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl

import functions as f
import settings as s
import IO

class Run(object):
    def __init__(self, LIEdir, *args, **kwargs):
        self.LIEdir  = LIEdir
        self.vdw     = []
        self.el      = []
        self.avg_vdw = []
        self.avg_el  = []
        self.fl_vdw  = []
        self.fl_el   = []

    def fl_en_rep(self, en):
        a = np.cumsum(en, axis=1)
        b = np.indices(a.shape)
        b = b[1] + 1
        c = a / b
        return c

    def fl_en_avg(self, en):
        a = np.cumsum(en)
        b = np.indices(a.shape)
        b = b + 1
        c = a / b
        return c
        
    def read_LIE(self):
        LIEs     = {}
        self.MDs  = 0
        self.Reps = glob.glob(self.LIEdir + '/FEP1/*/*')
        for rep in self.Reps:
            LIEs[rep] = sorted(glob.glob(rep + '/md_LIE_*.log'))
            if len(LIEs[rep]) > self.MDs:
                self.MDs = len(LIEs[rep])
    
        for x,rep in enumerate(LIEs):
            for i,LIE_MD in enumerate(LIEs[rep]):
            
                with open(LIE_MD) as infile:
                    vdw_n   = []
                    el_n    = []                
                    for line in infile:
                        line = line.split()
                        if len(line) > 3:
                            if line[0] == 'Q-surr.' and line[1] == '1':
                                vdw_n.append(float(line[4]))
                                el_n.append(float(line[3]))

                    if i == 0: # if ../md_LIE_01.log file, append new list to energy array's
                        self.vdw.append(vdw_n)
                        self.el.append(el_n)
                    else: # if non ../md_LIE_01.log file, join list to current replicate list
                        self.vdw[x] = [*self.vdw[x], *vdw_n]
                        self.el[x] = [*self.el[x], *el_n]

        
        self.vdw = [x for x in self.vdw if x != []] # discards empty lists of failed reps
        self.el  = [x for x in self.el if x != []]
        
        for rep in self.vdw: # adds nan to equalize all list lengths, in order to define arrays
            if len(rep) < (self.MDs*10001):
                x = (self.MDs*10001) - len(rep)
                for i in range(x):
                    rep.append(np.nan)
        for rep in self.el:
            if len(rep) < (self.MDs*10001):
                x = (self.MDs*10001) - len(rep)
                for i in range(x):
                    rep.append(np.nan)
        
        self.vdw = np.array(self.vdw)
        self.el = np.array(self.el)
        
        self.x_axis = range(0, len(self.vdw[0]))
        
        self.avg_vdw = np.nanmean(self.vdw, axis=0)
        self.avg_el  = np.nanmean(self.el, axis=0)
        
        self.avg_fl_vdw = self.fl_en_avg(self.avg_vdw)
        self.avg_fl_el = self.fl_en_avg(self.avg_el)
        
        self.fl_vdw = self.fl_en_rep(self.vdw)
        self.fl_el  = self.fl_en_rep(self.el)
        
        avg_vdw = []
        avg_el  = []
        
        for rep in self.vdw:
            avg_vdw.append(np.nanmean(rep)) # list of average vdw per rep
        vdw = np.nanmean(avg_vdw) # total average vdw
        self.vdw_sem = np.nanstd(avg_vdw)/np.sqrt(len(avg_vdw))
        
        for rep in self.el:
            avg_el.append(np.nanmean(rep)) # list of average el per rep
        el = np.nanmean(avg_el) # total average el
        self.el_sem = np.nanstd(avg_el)/np.sqrt(len(avg_el)) 

    def new_plot(self): # plots the floating mean of all individual reps and the total average floating mean
        mpl.style.use('seaborn')
        for i in range(len(self.fl_vdw)):
            plt.figure('QLIE floating mean all reps', figsize=(20,20))
            plt.plot(self.x_axis, self.fl_vdw[i], color='tab:blue', alpha=0.7, linewidth=1)
            plt.plot(self.x_axis, self.fl_el[i], color='tab:orange', alpha=0.7, linewidth=1)
        plt.errorbar(self.x_axis, self.avg_fl_vdw[0], yerr=self.vdw_sem, color='tab:blue', label='vdw', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.errorbar(self.x_axis, self.avg_fl_el[0], yerr=self.el_sem, color='tab:orange', label='el', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.xlabel('time (fs)')
        plt.ylabel('Q-surr (kcal/mol)')
        plt.title('Average Energy Development QLIE and floating mean '+self.LIEdir.split('/')[-1])
        plt.legend()
        
    def average_plot(self): # plots the total average energy and it's floating mean
        mpl.style.use('seaborn')
        plt.figure('QLIE average energy and floating mean', figsize=(20,20))
        plt.plot(self.x_axis, self.avg_vdw, color='tab:blue', alpha=0.7, linewidth=1)
        plt.plot(self.x_axis, self.avg_el, color='tab:orange', alpha=0.7, linewidth=1)
        plt.errorbar(self.x_axis, self.avg_fl_vdw[0], yerr=self.vdw_sem, color='tab:blue', label='vdw', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.errorbar(self.x_axis, self.avg_fl_el[0], yerr=self.el_sem, color='tab:orange', label='el', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.xlabel('time (fs)')
        plt.ylabel('Q-surr (kcal/mol)')
        plt.title('Energy Development QLIE floating mean all reps '+self.LIEdir.split('/')[-1])
        plt.legend()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='plotLIE reps',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Plot the LIE energies of individual reps == ')
    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    args = parser.parse_args()
    run = Run(LIEdir = args.LIEdir)
    run.read_LIE()
    run.new_plot()
    run.average_plot()
    plt.show()

