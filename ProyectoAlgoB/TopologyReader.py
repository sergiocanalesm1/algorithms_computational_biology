#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 21:45:29 2020

@author: David Ricardo Figueora Blanco 
@email: dr.figueroa10@uniandes.edu.co 

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
## modules
import numpy as np 
import matplotlib.pyplot as plt
import os 
import sys
import re
import subprocess
'''
Read martini.itp to identify character and creat matrix of LJ parameters

'''
with open("./informartiniff.itp", "r") as file_handle:
    file_contents = file_handle.readlines()
martff = []
for line in file_contents:
    martff.append(re.split(r"(\s+)", line))


atomtypes = dict()
count = 0
for i in martff:
    # print(i[0])
    if (i[0] != ";"):
        if(len(i) < 14 and len(i) > 10):
            atomtypes[i[0]] = [i[0][0], count]
            count += 1

with open("./LJ_interac.txt", "r") as file_handle:
    file_contents = file_handle.readlines()
LJ_level = []
for line in file_contents:
    LJ_level.append(line.split(":"))

LJ_info = dict()
for i in LJ_level:
    LJ_info[i[0]] = [float(i[1]), float(i[2][:-1])]

matrixLJ = [[None]*39 for _ in range(39)]

for i in martff:
    if (i[0] != ";"):
        if(len(i) > 14):
            idx1 = atomtypes[i[2]][1]
            idx2 = atomtypes[i[4]][1]
            print(atomtypes[i[2]],atomtypes[i[4]])
            inttype = ("".join(i[14:-2])).split(",")[0]
            if(inttype[0:2] == "75"):
                # print(LJ_info[inttype[2:]][0]*0.75,LJ_info[inttype[2:]][1]-0.04)
                matrixLJ[idx1][idx2] = [LJ_info[inttype[2:]]
                                        [0]*0.75, LJ_info[inttype[2:]][1]-0.04]
            else:
                # print(inttype)
                # print(LJ_info[inttype])
                matrixLJ[idx1][idx2] = LJ_info[inttype]
                


'''
read protein itp  and extract bonds and atoms information. 
return dictionary with connectivity and atom names based on atom index
'''
                
with open("./Protein_A.itp","r")  as Prot_itpfile:
        file_contents = Prot_itpfile.readlines()
# get lines of bonds, angle and atoms creat dictionaries.   
atominit = subprocess.getoutput("grep  -n '\[ atoms \]' Protein_A.itp").split(":")[0]
bondsinit = subprocess.getoutput("grep  -n '\[ bonds \]' Protein_A.itp").split(":")[0]
anglinit = subprocess.getoutput("grep  -n '\[ angles \]' Protein_A.itp").split(":")[0]   

atoms_def=dict()
bonds_def=dict()
for i  in range(len(file_contents)):
    if(i>=int(atominit) and i< int(bondsinit)-2 ):
        line = re.split(r"(\s+)",file_contents[i])
        atoms_def[line[12]]=[line[4],line[14]]
    if(i>int(bondsinit) and i < int(anglinit)-2 ):
        line = re.split(r"(\s+)",file_contents[i])
        if(len(line)>8):
            bonds_def[line[2]]=line[4]
    i += 1

