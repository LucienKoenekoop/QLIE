import argparse
import glob
import numpy as np
import os

import functions as f
import settings as s
import IO

class Run(object):
    """
    """
    def __init__(self, LIEdir, rep, *args, **kwargs):
        self.LIEdir = LIEdir
        self.rep    = rep
        self.vdw    = []
        self.el     = []
        self.convergence = 0
        
    def read_LIE(self):
        LIEs = sorted(glob.glob(self.LIEdir + '/*/*/'+self.rep+'/md_LIE_*.log'))
        for i,LIE in enumerate(LIEs):
            with open(LIE) as infile:
                vdw_n   = []
                el_n    = []                
                for line in infile:
                    line = line.split()
                    if len(line) > 3:
                        if line[0] == 'Q-surr.' and line[1] == '1':
                            vdw_n.append(float(line[4]))
                            el_n.append(float(line[3]))
                            
                if i == 0:
                    self.vdw = vdw_n
                    self.el = el_n
                else:
                    self.vdw = [*self.vdw, *vdw_n]
                    self.el  = [*self.el, *el_n]

        self.vdw = np.array(self.vdw)
        self.el = np.array(self.el)

    def calc_convergence(self):
        fl_vdw = []
        for i in range(1, len(self.vdw)+1):
            fl_vdw.append(np.average(self.vdw[0:i]))
        vdw_gradient = np.gradient(np.array(fl_vdw))
        self.vdw_gradient = np.array(vdw_gradient)
        
        el_gradient = []
        fl_el = []
        for i in range(1, len(self.el)+1):
            fl_el.append(np.average(self.el[0:i]))
        el_gradient = np.gradient(np.array(fl_el))
        self.el_gradient = np.array(el_gradient)
        
        d0 = abs(0.001)
        d1 = abs(0.0001)
        d2 = abs(0.00001)
        d3 = abs(0.000001)
        limit = 1000
        r = (100 / 70000)
        switch_vdw = False
        switch_el = False
        self.vdw_r = []
        self.el_r = []
        
        curset = []
        bin1, bin2, bin3 = 0, 0, 0
        for i,grad in enumerate(self.vdw_gradient):
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
                    switch_vdw = True
                    break

                if abs(curset[0]) < d0 and abs(curset[0]) > d1:
                    bin1 -= 1
                elif abs(curset[0]) < d1 and abs(curset[0]) > d2:
                    bin2 -= 1
                elif abs(curset[0]) < d2 and abs(curset[0]) > d3:
                    bin3 -= 1
                curset.pop(0)
                    
        curset = []
        bin1, bin2, bin3 = 0, 0, 0
        for i,grad in enumerate(self.el_gradient):
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
                    switch_el = True
                    break

                if abs(curset[0]) < d0 and abs(curset[0]) > d1:
                    bin1 -= 1
                elif abs(curset[0]) < d1 and abs(curset[0]) > d2:
                    bin2 -= 1
                elif abs(curset[0]) < d2 and abs(curset[0]) > d3:
                    bin3 -= 1
                curset.pop(0)
                
        if switch_vdw == True and switch_el == True:
            self.convergence = 1
        else:
            pass
        print(self.convergence)
        return self.convergence
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='stopLIE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Stop LIE when converged == ')
    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    parser.add_argument('-r', '--rep',
                        dest = "rep",
                        required = True,
                        help = "current rep")
    
    args = parser.parse_args()
    run = Run(LIEdir = args.LIEdir,
              rep    = args.rep)
    run.read_LIE()
    run.calc_convergence()