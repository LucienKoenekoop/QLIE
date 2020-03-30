import argparse
import glob
import numpy as np
import matplotlib
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import Draw

import functions as f
import settings as s
import IO

class Run(object):
    
    def __init__(self, ligand):
        self.ligand    = ligand
        self.mol       = Chem.rdmolfiles.MolFromMolFile(self.ligand, removeHs=False) 
        self.FG_Smarts = {                              # non-zero functional groups
            'AL':   '[#6][OX2H]',                               # Alcohols
            'AD1':  '[NX3;H2][CX3](=[OX1])[#6]',                # Primary Amide
            'AN12': '[NX3;H2,H1;!$(NC=O)]',                     # Primary of Secondary Amine
            'CA':   '[CX3](=O)[OX2H1]'                         # Carboxylic Acid
        }
        self.OG_Smarts = {                             # other (zero) groups
            'AD2':  '[NX3;H1][CX3](=[OX1])[#6]',                # Secondary Amide
            'AD3':  '[NX3;H0][CX3](=[OX1])[#6]',                # Tertiary Amide
            'AN3':  '[NX3;H0;!$(NC=[!#6]);!$(NC#[!#6])]',        # Tertiary Amines
            'KT':   '[#6][CX3](=O)[#6]',                        # Ketone
            'ADH':  '[CX3H1](=O)[#6]',                          # Aldehyde
            'TH':   '[#16X2H]',                                 # Thiol
            'ET':   '[OD2]([#6])[#6]',                          # Ether
            'ES':   '[#6][CX3](=O)[OX2H0][#6]',                 # Ester
            'NL':   '[NX1]#[CX2]',                              # Nitrile
            'NO':   '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]',  # Nitro
            'SU':   '[#16X2H0]'                                # Sulfide
        }
        self.FG = {}
        self.OG = {}
        
        self.charge = Chem.rdmolops.GetFormalCharge(self.mol)
        if self.charge < 0:
            self.FG['ANI']  = abs(self.charge)
            self.FG['CAT']  = 0
        elif self.charge > 0:
            self.FG['ANI']  = 0
            self.FG['CAT']  = abs(self.charge)
        else:
            self.FG['ANI']  = 0
            self.FG['CAT']  = 0

    
    def debug(self, FG):
        X = Chem.MolFromSmarts(self.OG_Smarts[FG])
        print(FG, Chem.MolToSmiles(X))
        match = self.mol.GetSubstructMatches(X)
        print(match, len(match))
        
    def Get_FG(self, X):
        FG = Chem.MolFromSmarts(self.FG_Smarts[X])
        match = self.mol.GetSubstructMatches(FG)
        self.FG[X] = len(match)
        return self.FG[X]
    
    def Get_OG(self, X):
        OG = Chem.MolFromSmarts(self.OG_Smarts[X])
        match = self.mol.GetSubstructMatches(OG)
        self.OG[X] = len(match)
        return self.OG[X]

    def beta(self):
        b0  = 0.43
        w_n = 1
        w_c = 11
        w_t = sum(self.OG.values()) + self.FG['CA'] + self.FG['AN12'] + self.FG['AD1'] + self.FG['AL'] + 11*(self.FG['ANI'] + self.FG['CAT'])
        db  = {
            'AL':   [-0.06, 1], 
            'AN12': [-0.04, 1],
            'AD1':  [-0.02, 1],
            'CA':   [-0.03, 1],
            'ANI':  [0.02, 11],
            'CAT':  [0.09, 11]
        }
        w_b = 0
        for i in db:
            w_b += db[i][0]*db[i][1]*self.FG[i]
        try:
            beta = b0 + (w_b/w_t)
        except:
            beta = b0
        return beta
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Beta',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Calculate LIE beta parameter of ligand sdf/mol file == ')
    
    parser.add_argument('-l', '--ligand',
                        dest = "ligand",
                        required = True,
                        help = "name of ligand")
    args = parser.parse_args()
    run = Run(ligand = args.ligand)
    
    for FG in run.FG_Smarts:
        run.Get_FG(FG)
    for OG in run.OG_Smarts:
        run.Get_OG(OG)
    run.beta()