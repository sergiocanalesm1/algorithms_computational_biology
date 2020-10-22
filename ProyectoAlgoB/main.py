#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:45:45 2020

@author: David Ricardo Figueora Blanco
@email: dr.figueroa10@uniandes.edu.co


   PROGRAM:   program name
   AUTHOR:    your name


   FUNCTION:  a short paragraph stating the purpose of the
              program.

   INPUT:     location of input, i.e.  standard input, a file on
              disk

   OUTPUT:    location and type of output, i.e.  a report
              containing a detail record for each city processed
              containing city id, Celsius temperature, Fahrenheit
              temperature and wind chill temperature.

   NOTES:     any relevant information that would be of
              additional help to someone looking at the program

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import time
#from rdkit import Chem
#from rdkit.Chem import ChemicalFeatures
#from rdkit import RDConfig
from Graph import createGraph, calcDist
#%%

with open("./1UBQ-CG.pdb","r") as pdbfile:
    content = pdbfile.readlines()


first_graph = createGraph( content, 16 )

'''
for k,v in overlap.items():
    print(k.get_Name(),k.get_resName(),k.get_resNum())
    for inter in v:
        print(inter[0].get_Name())
    print()


with open("aminoacids.raw","r") as aadata:
    content = aadata.readlines()


#fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
#factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

#m=Chem.MolFromSmiles(content[0].split("\t")[4])
#feat = factory.GetFeaturesForMol(m)
'''