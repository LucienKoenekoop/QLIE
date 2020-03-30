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
        self.y_vdw   = {}
        self.y_fl_vdw = {}
        self.y_el    = {}
        self.y_fl_el = {}

    def setup_plot(self, en, rep, i):
        en[i] = []
        for energy in rep:
            en[i].append(energy)
        return en[i]
            
    def fl_en(self, flen, rep, i):
        flen[i] = []
        for j in range(1, len(rep)+1):
            fl_en = np.average(rep[0:j])
            flen[i].append(fl_en)
        flen[i] = np.array(flen[i])
        return flen[i]
        
    def read_LIE(self):
        LIEs = sorted(glob.glob(self.LIEdir + '/*/*/*/md_LIE_*.log'))
        self.Reps = len(glob.glob(self.LIEdir + '/*/*/*'))
        self.MDs = int(len(LIEs)/self.Reps)
        for i,LIE in enumerate(LIEs):
            x = i // int(len(LIEs)/self.Reps)
            y = float(i) / int(len(LIEs)/self.Reps)

            with open(LIE) as infile:
                vdw_n   = []
                el_n    = []
                total_n = []
                for line in infile:
                    line = line.split()
                    if len(line) > 3:
                        if line[0] == 'Q-surr.' and line[1] == '1':
                            vdw_n.append(float(line[4]))
                            el_n.append(float(line[3]))
                            total_n.append(float(line[4]) + float(line[3]))
                
                if x == y:
                    self.vdw.append(vdw_n)
                    self.el.append(el_n)
                else:
                    self.vdw[x] = [*self.vdw[x], *vdw_n]
                    self.el[x] = [*self.el[x], *el_n]

        self.vdw = [x for x in self.vdw if x != []] 
        self.el  = [x for x in self.el if x != []]
        
        for i,rep in enumerate(self.vdw): # discards lists of prematurely crashed reps
            if len(rep) < (self.MDs*10001):
                self.vdw.pop(i)
        for i,rep in enumerate(self.el):
            if len(rep) < (self.MDs*10001):
                self.el.pop(i)
        
        self.vdw   = np.array(self.vdw)
        self.el    = np.array(self.el)
        
        self.x_axis = range(0, len(self.vdw[0]))
        
        self.avg_vdw = np.nanmean(self.vdw, axis=0)
        self.avg_el  = np.nanmean(self.el, axis=0)

        for i in range(1, len(self.avg_vdw)+1):
            fl_vdw = np.average(self.avg_vdw[0:i])
            fl_el  = np.average(self.avg_el[0:i])
            self.fl_vdw.append(fl_vdw)
            self.fl_el.append(fl_el)
        self.fl_vdw = np.array(self.fl_vdw)
        self.fl_el  = np.array(self.fl_el)
        
        vdw          = np.nanmean(self.avg_vdw)
        self.vdw_std = np.nanstd(self.avg_vdw)
        self.fl_vdw_std = np.nanstd(self.fl_vdw, ddof=1)
        
        el           = np.nanmean(self.avg_el)
        self.el_std  = np.nanstd(self.avg_el)
        self.fl_el_std = np.nanstd(self.fl_el, ddof=1)
        
    def plot_vdw_el(self):
        y_axis = {'vdw': [], 'fl_vdw': [], 'el': [], 'fl_el': []}
        x_axis = range(0, len(self.avg_vdw))
        for energy in self.avg_vdw:
            y_axis['vdw'].append(energy)
        for energy in self.fl_vdw:
            y_axis['fl_vdw'].append(energy)
        for energy in self.avg_el:
            y_axis['el'].append(energy)
        for energy in self.fl_el:
            y_axis['fl_el'].append(energy)
                
        mpl.style.use('seaborn')
        plt.figure('QLIE', figsize=(20,20))
        plt.plot(x_axis,y_axis['vdw'], color='0.3', linewidth=1)
        plt.plot(x_axis,y_axis['el'], color='0.7', linewidth=1)
        plt.errorbar(x_axis,y_axis['fl_vdw'],yerr=self.fl_vdw_std, color='0.1', label='vdw', linewidth=1, errorevery=100000, zorder=3)
        plt.errorbar(x_axis,y_axis['fl_el'],yerr=self.fl_el_std, color='0.5', label='el', linewidth=1, errorevery=100000, zorder=3)
        plt.xlabel('time (fs)')
        plt.ylabel('Q-surr (kcal/mol)')
        plt.title('Energy Development '+self.LIEdir.split('/')[-1])
        plt.legend()
        plt.show()
        
    def plot_reps(self): # plot the floating mean of all individal reps and their combined average floating mean
        y_axis = {'vdw': [], 'fl_vdw': [], 'el': [], 'fl_el': []}
        x_axis = range(0, len(self.avg_vdw))
        for energy in self.avg_vdw:
            y_axis['vdw'].append(energy)
        for energy in self.fl_vdw:
            y_axis['fl_vdw'].append(energy)
        for energy in self.avg_el:
            y_axis['el'].append(energy)
        for energy in self.fl_el:
            y_axis['fl_el'].append(energy)
        
        mpl.style.use('seaborn')
        for y in self.y_fl_vdw:
            plt.figure('QLIE floating mean all reps', figsize=(20,20))
            plt.plot(self.x_axis,self.y_fl_vdw[y], color='0.3', linewidth=1)
            plt.plot(self.x_axis,self.y_fl_el[y], color='0.7', linewidth=1)
        plt.errorbar(x_axis,y_axis['fl_vdw'],yerr=self.fl_vdw_std, color='0.1', label='vdw', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.errorbar(x_axis,y_axis['fl_el'],yerr=self.fl_el_std, color='0.5', label='el', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.xlabel('time (fs)')
        plt.ylabel('Q-surr (kcal/mol)')
        plt.title('Energy Development '+self.LIEdir.split('/')[-1])
        plt.legend()
        
    def plot(self, y): # plot each individual rep vdw and el with floating mean
        mpl.style.use('seaborn')
        plt.figure('QLIE rep '+str(y+1), figsize=(20,20))
        plt.plot(self.x_axis,self.y_vdw[y], color='0.3', linewidth=1)
        plt.plot(self.x_axis,self.y_el[y], color='0.7', linewidth=1)
        plt.plot(self.x_axis,self.y_fl_vdw[y], color='0.1', label='vdw', linewidth=1)
        plt.plot(self.x_axis,self.y_fl_el[y], color='0.5', label='el', linewidth=1)
        plt.xlabel('time (fs)')
        plt.ylabel('Q-surr (kcal/mol)')
        plt.title('Energy Development '+self.LIEdir.split('/')[-1]+' rep '+str(y+1))
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

    for i,rep in enumerate(run.vdw):
        run.setup_plot(run.y_vdw, rep, i)
        run.fl_en(run.y_fl_vdw, rep, i)
    for i,rep in enumerate(run.el):
        run.setup_plot(run.y_el, rep, i)
        run.fl_en(run.y_fl_el, rep, i)
    for y in run.y_vdw:
        run.plot(y)
    run.plot_reps()
    plt.show()
