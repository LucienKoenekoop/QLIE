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
        self.LIEdir = LIEdir.strip('/')
        self.vdw    = []
        self.el     = []
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
        LIEs = {}
        self.MDs = 0
        self.Reps = sorted(glob.glob(self.LIEdir + '/FEP1/*/*'))
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
                        self.el[x]  = [*self.el[x], *el_n]

        self.vdw = [x for x in self.vdw if x != []] # discards empty lists of failed reps
        self.el  = [x for x in self.el if x != []]
    
    def Q_test(self, data):
        Q = {3: 0.941, 4: 0.765, 5: 0.642, 6: 0.560, 7: 0.507, 8: 0.437, 9: 0.437, 10: 0.412}
        Q1 = abs(data[0]-data[1])/(data[-1]-data[0])
        Q2 = abs(data[-1]-data[-2])/(data[-1]-data[0])
        if Q1 > Q[len(data)]:
            data.pop(0)
        if Q2 > Q[len(data)]:
            data.pop(-1)
        return data
        
    def calc_LIE(self):
        avg_vdw = []
        avg_el  = []
        
        for rep in self.vdw:
            avg_vdw.append(np.nanmean(rep)) # list of average vdw per rep
            
#        avg_vdw.sort() # Remove outliers
#        while True:
#            a = avg_vdw[:]
#            self.Q_test(avg_vdw)
#            if a == avg_vdw:
#                break
        
        vdw = np.nanmean(avg_vdw) # total average vdw
#        vdw = np.nanmedian(avg_vdw) # total median vdw
        vdw_sem = np.nanstd(avg_vdw)/np.sqrt(len(avg_vdw))
        
        for rep in self.el:
            avg_el.append(np.nanmean(rep)) # list of average el per rep
        
#        avg_el.sort() # Remove outliers
#        while True:
#            a = avg_el[:]
#            self.Q_test(avg_el)
#            if a == avg_el:
#                break
        
        el = np.nanmean(avg_el) # total average el
#        el = np.nanmedian(avg_el) # total median el
        el_sem = np.nanstd(avg_el)/np.sqrt(len(avg_el)) 
        
        print('vdw:', vdw, '+/-', vdw_sem, '| el:', el, '+/-', el_sem, 'beta', self.beta, ' ')
        dG_LIE = self.a * vdw + self.beta * el + self.g
        dG_sem = self.a * vdw_sem + self.beta * el_sem
        print('dG:', dG_LIE, '+/-', dG_sem)
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