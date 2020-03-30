import argparse
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl

import functions as f
import settings as s
import beta as b
import IO

class Run(object):
    """
    """
    def __init__(self, LIEdir, *args, **kwargs):
        self.LIEdir = LIEdir
        self.path   = os.getcwd()
        self.vdw    = []
        self.el     = []
        self.fl_vdw = []
        self.fl_el  = []
        self.convergence = 0
        
        
    def read_LIE(self):
        LIEs = sorted(glob.glob(self.LIEdir + '/*/*/*/md_LIE_*.log'))
        Reps = len(glob.glob(self.LIEdir + '/*/*/*'))
        self.MDs = int(len(LIEs)/Reps)
        for i,LIE in enumerate(LIEs):
            x = i // int(len(LIEs)/Reps) # counter for the current replicate
            y = float(i) / int(len(LIEs)/Reps) # switch value to recognize new replicate
            
            with open(LIE) as infile:
                vdw_n   = []
                el_n    = []                
                for line in infile:
                    line = line.split()
                    if len(line) > 3:
                        if line[0] == 'Q-surr.' and line[1] == '1':
                            vdw_n.append(float(line[4]))
                            el_n.append(float(line[3]))
                            
                if x == y: # if ../md_LIE_01.log file, append new list to energy array's
                    self.vdw.append(vdw_n)
                    self.el.append(el_n)
                else: # if non ../md_LIE_01.log file, join list to current replicate list
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

        self.vdw = np.array(self.vdw)
        self.el = np.array(self.el)
        
    def calc_LIE_avg(self):
        avg_vdw = np.nanmean(self.vdw, axis=0) # array of vdW mean's for each timestep
        avg_el = np.nanmean(self.el, axis=0) # array of el mean's for each timestep
        
        for i in range(1, len(avg_vdw)+1):
            fl_vdw = np.average(avg_vdw[0:i])
            fl_el  = np.average(avg_el[0:i])
            self.fl_vdw.append(fl_vdw)
            self.fl_el.append(fl_el)
        self.fl_vdw = np.array(self.fl_vdw)
        self.fl_el  = np.array(self.fl_el)

        vdw = np.nanmean(avg_vdw) # total average vdW
        vdw_std = np.nanstd(avg_vdw) # standard deviation over the average vdW values
        fl_vdw_std = np.nanstd(self.fl_vdw, ddof=1)
        
        el = np.nanmean(avg_el) # total average el
        el_std = np.nanstd(avg_el) # standard deviation over the average el values
        fl_el_std = np.nanstd(self.fl_el, ddof=1)
        
        vdw_gradient = np.gradient(self.fl_vdw)
        el_gradient = np.gradient(self.fl_el)

    def calc_convergence(self):
        vdw_gradient = []
        for rep in self.vdw:
            fl_vdw = []
            for i in range(1, len(rep)+1):
                fl_vdw.append(np.average(rep[0:i]))
            vdw_gradient.append(np.gradient(np.array(fl_vdw)))
        self.vdw_gradient = np.array(vdw_gradient)
        
        el_gradient = []
        for rep in self.el:
            fl_el = []
            for i in range(1, len(rep)+1):
                fl_el.append(np.average(rep[0:i]))
            el_gradient.append(np.gradient(np.array(fl_el)))
        self.el_gradient = np.array(el_gradient)
        
        d0 = abs(0.001)
        d1 = abs(0.0001)
        d2 = abs(0.00001)
        d3 = abs(0.000001)
        limit = 1000
        r = (100 / 70000)
        switch_vdw = 0
        switch_el = 0
        self.vdw_r = []
        self.el_r = []
        
        for k,rep in enumerate(self.vdw_gradient):
            curset = []
            bin1, bin2, bin3 = 0, 0, 0
            for i,grad in enumerate(rep):
                curset.append(grad)
                if i < 1000:
                    continue
                elif i == 1000:
                    for j in curset:
                        if j == curset[0]:
                            continue
                        if abs(j) < d0 and abs(j) > d1:
                            bin1 += 1
                        elif abs(j) < d1 and abs(j) > d2:
                            bin2 += 1
                        elif abs(j) < d2 and abs(j) > d3:
                            bin3 += 1
                    curset.pop(0)
                else:
                    if abs(grad) < d0 and abs(grad) > d1:
                        bin1 += 1
                    elif abs(grad) < d1 and abs(grad) > d2:
                        bin2 += 1
                    elif abs(grad) < d2 and abs(grad) > d3:
                        bin3 += 1
                        
                    ratio = bin1 / (bin2 * bin3 + 1)
                    self.vdw_r.append(ratio)
                    if ratio <= r:
                        switch_vdw += 1
#                        print(i, 'vdw_done')
                        break
#                    if i == 20000 or i == 40000 or i == 60000 or i == 80000 or i == 99000:
#                        print(bin1, bin2, bin3, ratio)
    
                    if abs(curset[0]) < d0 and abs(curset[0]) > d1:
                        bin1 -= 1
                    elif abs(curset[0]) < d1 and abs(curset[0]) > d2:
                        bin2 -= 1
                    elif abs(curset[0]) < d2 and abs(curset[0]) > d3:
                        bin3 -= 1
                    curset.pop(0)
                    
        for k,rep in enumerate(self.el_gradient):
            curset = []
            bin1, bin2, bin3 = 0, 0, 0
            for i,grad in enumerate(rep):
                curset.append(grad)
                if i < 1000:
                    continue
                elif i == 1000:
                    for j in curset:
                        if j == curset[0]:
                            continue
                        if abs(j) < d0 and abs(j) > d1:
                            bin1 += 1
                        elif abs(j) < d1 and abs(j) > d2:
                            bin2 += 1
                        elif abs(j) < d2 and abs(j) > d3:
                            bin3 += 1
                    curset.pop(0)
                else:
                    if abs(grad) < d0 and abs(grad) > d1:
                        bin1 += 1
                    elif abs(grad) < d1 and abs(grad) > d2:
                        bin2 += 1
                    elif abs(grad) < d2 and abs(grad) > d3:
                        bin3 += 1
                        
                    ratio = bin1 / (bin2 * bin3 + 1)
                    self.el_r.append(ratio)
                    if ratio <= r:
                        switch_el += 1
#                        print(i, 'el_done')
                        break
#                    if i == 20000 or i == 40000 or i == 60000 or i == 80000 or i == 99000:
#                        print(bin1, bin2, bin3, ratio)
    
                    if abs(curset[0]) < d0 and abs(curset[0]) > d1:
                        bin1 -= 1
                    elif abs(curset[0]) < d1 and abs(curset[0]) > d2:
                        bin2 -= 1
                    elif abs(curset[0]) < d2 and abs(curset[0]) > d3:
                        bin3 -= 1
                    curset.pop(0)
                    
        if switch_vdw == len(self.vdw_gradient) and switch_el == len(self.el_gradient):
            self.convergence = 1
        else:
            pass
        print(self.convergence)
        return self.convergence
    
    def plot_gradient(self):
        x_axis = range(0, len(self.vdw_2gradient[0]))
        y_vdw = {}
        for i, rep in enumerate(self.vdw_2gradient):
            y_vdw[i] = []
            for energy in rep:
                y_vdw[i].append(energy)
        y_el = {}
        for i, rep in enumerate(self.el_2gradient):
            y_el[i] = []
            for energy in rep:
                y_el[i].append(energy)
        
        mpl.style.use('seaborn')
        for y in y_vdw:
            plt.figure('QLIE vdw and el floating mean gradient all reps '+str(y), figsize=(20,20))
            plt.plot(x_axis,y_vdw[y], color='0.3', label='vdw', linewidth=1)
            plt.plot(x_axis,y_el[y], color='0.7', label='el', linewidth=1)
            plt.xlabel('time (fs)')
            plt.ylabel('dQ-surr/dt (kcal/(mol*s))')
            plt.title('Energy Development '+self.LIEdir.split('/')[-1])
            plt.ylim(bottom=-0.01, top=0.01)
            plt.legend()
        plt.show()
    
    def plot_r(self):
        x_axis = range(0, len(self.vdw_r))
        mpl.style.use('seaborn')
        plt.figure('gradient ratio', figsize=(20,20))
        plt.plot(x_axis, self.vdw_r, label='vdw', color='0.3', linewidth=1)
        plt.plot(x_axis, self.el_r, label='el', color='0.7', linewidth=1)
        plt.xlabel('time (fs)')
        plt.ylabel('r')
        plt.title('2YDO')
        plt.ylim(bottom=0, top=0.1)
        plt.legend()
        plt.show()
            
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='stopLIE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Stop LIE when converged == ')

    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    args = parser.parse_args()
    run = Run(LIEdir = args.LIEdir)
    run.read_LIE()
    run.calc_convergence()