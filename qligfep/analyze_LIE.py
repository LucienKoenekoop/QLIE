import argparse
import glob
import numpy as np
import os

import functions as f
import settings as s
import beta as b
import IO

class Run(object):
    """
    """
    def __init__(self, LIEdir, *args, **kwargs):
        self.LIEdir = LIEdir
        self.vdw    = []
        self.el     = []
        self.fl_vdw = []
        self.fl_el  = []
        self.a      = 0.18
        self.b0     = 0.33
        self.beta   = 0
        self.g      = 0
        self.path   = os.getcwd()
        
    def get_beta(self):
        
        ligfile = glob.glob('*sdf')[0]
        lig     = self.path+'/'+ligfile
        beta = b.Run(ligand = lig)
        for FG in beta.FG_Smarts:
            beta.Get_FG(FG)
        for OG in beta.OG_Smarts:
            beta.Get_OG(OG)
        self.beta = beta.beta()
        return self.beta
        
    def read_LIE(self):
        LIEs = sorted(glob.glob(self.LIEdir + '/*/*/*/md_LIE_*.log'))
        Reps = len(glob.glob(self.LIEdir + '/FEP1/*/*'))
        self.MDs = int(len(LIEs)/Reps)
        for i,LIE in enumerate(LIEs):
            print(i)
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

        self.vdw = [x for x in self.vdw if x != []] # discards empty lists of failed reps
        self.el  = [x for x in self.el if x != []]
        
        for i,rep in enumerate(self.vdw): # discards lists of prematurely crashed reps
            if len(rep) < (self.MDs*10001):
                self.vdw.pop(i)
        for i,rep in enumerate(self.el):
            if len(rep) < (self.MDs*10001):
                self.el.pop(i)
                
        self.vdw = np.array(self.vdw)
        self.el = np.array(self.el)
        
    def calc_LIE(self):
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
        
        print('vdw:', vdw, '+/-', fl_vdw_std, '| el:', el, '+/-', fl_el_std, 'beta', self.beta, ' ')
        dG_LIE = self.a * vdw + self.beta * el + self.g
        print('dG:', dG_LIE)
        self.energy = {'vdw': vdw, 'el': el}
        return self.energy
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='analyzeLIE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse LIE == ')
    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    args = parser.parse_args()
    run = Run(LIEdir = args.LIEdir)
    run.get_beta()
    run.read_LIE()
    run.calc_LIE()
