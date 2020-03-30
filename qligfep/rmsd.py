import os
import argparse
import sys
import glob
import numpy as np

import IO

class rmsd(object):
    def __init__(self, LIEdir, rep, option, *args, **kwargs):
        self.LIEdir   = LIEdir.strip('/')
        self.rep      = rep
        self.option   = option
        self.top      = []
        self.protein  = []
        self.ligand   = []
        self.backbone = []
        self.residues = []
        self.dualtop  = self.LIEdir+'/FEP1/298/'+rep+'/dualtop.top'
        self.MD       = glob.glob(self.LIEdir+'/FEP1/298/'+rep+'/md_LIE_*.dcd')[-1]
        self.inp      = self.LIEdir+'/inputfiles/rmsd_'+option+'.inp'
        self.out      = self.LIEdir+'/FEP1/298/'+rep+'/rmsd_'+option+'.txt'
        
        top = LIEdir+'/inputfiles/top_p.pdb'
        
        with open(top, 'r') as infile:
            for line in infile:
                self.top.append(IO.pdb_parse_in(line))
        
    def General(self):
        for atom in self.top:
            if len(atom[0]) != 6:
                continue
            try:
                IO.AA(atom[4])
                self.protein.append(atom[1])
            except:
                if atom[4] == 'LIG':
                    self.ligand.append(atom[1])
                else:
                    continue
        self.commands = [self.dualtop, '1', str(self.protein[0])+' '+str(self.protein[-1]),
                         str(self.ligand[0])+' '+str(self.ligand[-1]), 'end', 'go',
                         self.MD, '.']
        with open(self.inp, 'w') as outfile:
            for line in self.commands:
                outfile.write(line+'\n')
#       
    def Backbone(self):
        for atom in self.top:
            if atom[2] == 'CA':
                self.backbone.append(atom[1])
        self.commands = [self.dualtop, '1']
        for CA in self.backbone:
            self.commands.append(str(CA))
        self.commands.append('end')
        self.commands.append('go')
        self.commands.append(self.MD)
        self.commands.append('.')
        with open(self.inp, 'w') as outfile:
            for line in self.commands:
                outfile.write(line+'\n')
#        
    def Residues(self, option):
        if option == 2:    
            resnum = 1
            residue = []
            for atom in self.top:
                if atom[6] == resnum:
                    residue.append(atom[1])
                else:
                    self.residues.append(str(residue[0])+' '+str(residue[-1]))
                    residue.clear()
                    residue.append(atom[1])
                    resnum += 1

                if atom[2] == 'OXT':
                    self.residues.append(str(residue[0])+' '+str(residue[-1]))
                    break
                    
        elif option == 1:
            for atom in self.top:
                if atom[2] == 'CB':
                    self.residues.append(str(atom[1]))
        self.commands = []
        
        x = -1
        for i,resn in enumerate(self.residues, 1):
            if i//10 != x:
                self.commands.append([self.dualtop])
            self.commands[x].append('1')
            self.commands[x].append(resn)
            self.commands[x].append('end')
            if i//10 != x:
                if i == 1:
                    x += 1
                    continue
                self.commands[x].append('go')
                self.commands[x].append(self.MD)
                self.commands[x].append('.')
                x += 1
        self.commands[x].append('go')
        self.commands[x].append(self.MD)
        self.commands[x].append('.')
        
        for i,j in enumerate(self.commands, 1):
            with open(self.LIEdir+'/inputfiles/rmsd_'+self.option+str(i)+'.inp', 'w') as outfile:
                for line in j:
                    outfile.write(line+'\n')
        
    def Ligand(self):
        for atom in self.top:
            if atom[4] == 'LIG':
                self.ligand.append(atom[1])
                
        self.commands = [self.dualtop, '1', str(self.ligand[0])+' '+str(self.ligand[-1]),
                         'end', 'go', self.MD, '.']
        with open(self.inp, 'w') as outfile:
            for line in self.commands:
                outfile.write(line+'\n')                

    def qcalc(self, infile, outfile):
        qcalc = '/home/koenekoop/software/q6/bin/qcalc'
        inp = ' < '+infile
        out = ' > '+outfile
        IO.run_command(qcalc, inp + out, string=True)
        
    def rmsd(self):
        frames = []
        with open(self.out) as infile:
            for line in infile:
                if line[0:3] == 'LIE':
                    frames.append(float(line[26:34]))
        rmsd = np.average(frames)
        return str(rmsd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='RMSD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Calculates the RMSD of the option requested == ')
    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    parser.add_argument('-r', '--rep',
                        dest = "rep",
                        required = True,
                        help = "rep to evaluate")
    
    parser.add_argument('-o', '--option',
                        dest = "option",
                        required = True,
                        choices = ['general', 'backbone', 'residues', 'ligand'],
                        help = "RMSA option to calculate")
    
    args = parser.parse_args()
    run = rmsd(LIEdir = args.LIEdir,
               rep    = args.rep,
               option = args.option)
    
    if run.option == 'general':
        run.General()
        run.qcalc(run.inp, run.out)
        print(run.option+' rmsd: '+run.rmsd())

    elif run.option == 'backbone':
        run.Backbone()
        run.qcalc(run.inp, run.out)
        print(run.option+' rmsd: '+run.rmsd())

    elif run.option == 'residues':
        while True:
                option = int(input('Side chain CBs (1) or All residue atoms (2)? '))
                if option == 1 or option == 2:
                    break
                else:
                    print("That's not a valid option!")
        
        run.Residues(option)
        inputs = glob.glob(run.LIEdir+'/inputfiles/rmsd_'+run.option+'*.inp')
        for i,infile in enumerate(inputs):
            run.qcalc(infile, run.LIEdir+'/FEP1/298/'+run.rep+'/rmsd_'+run.option+str(i)+'.txt')
        outfiles = glob.glob(run.LIEdir+'/FEP1/298/'+run.rep+'/rmsd_'+run.option+'*.txt')
        with open(run.LIEdir+'/FEP1/298/'+run.rep+'/rmsd_'+run.option+'.txt', 'w') as outfile:
            for i,file in enumerate(outfiles):
                with open(file) as infile:
                    outfile.write('Residues '+str(i*10-9)+'-'+str(i*10)+'\n')
                    outfile.write('----------------------------- Calculation results ----------------------------\n')
                    outfile.write('file                 frame  1: RMSD(A)  2: RMSD(A)  3: RMSD(A)  4: RMSD(A)  5: RMSD(A)  6: RMSD(A)  7: RMSD(A)  8: RMSD(A)  9: RMSD(A) 10: RMSD(A)\n')
                    for line in infile:
                        if line[0:3] == 'LIE':
                            outfile.write(line)
                    outfile.write('\n')
        if outfiles[0] == run.LIEdir+'/FEP1/298/'+run.rep+'/rmsd_'+run.option+'.txt':
            outfiles.pop(0)
        else:
            pass
        for file in outfiles:
            os.remove(file)

    elif run.option == 'ligand':
        run.Ligand()
        run.qcalc(run.inp, run.out)
        print(run.option+' rmsd: '+run.rmsd())
        
    else:
        print('Error: invalid argument --option')
        
    print('outputfile: '+os.getcwd()+'/'+run.out)