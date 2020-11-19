#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:23:59 2020

@author: David Ricardo Figueora 
@email: dr.figueroa10

DESCRIPTION 
         A textual description of the functioning of the command or function. "
EXAMPLES 
         Some examples of common usage. 
SEE ALSO 
         A list of related commands or functions. "
BUGS 
         List known bugs. 
"""

import numpy as np
import matplotlib.pyplot as plt 
import sys 
from scipy.stats.stats import pearsonr
import time
import pandas as pd
#file_name = sys.argv[1]
file_name = "./expressionExampleGSE121594.tsv"


data_pd = pd.read_csv("./expressionExampleGSE121594.tsv",skiprows=0)

with open(file_name,"r") as exp_file:
    content = exp_file.readlines()
'''
data = []
for i in content:
    line = (i.split())
    #print(line[1:])
    if(line[1:] == ['0', '0', '0', '0', '0', '0', '0', '0', '0'] ):
        continue
    else:
        data.append(line)
data.pop(0)    
def create_coexp_mat(data):
    coexp_mat = np.eye((len(data)))
    for info_1 in range(len(data)):
        vec_1 = list(map(float,data[info_1][1:]))
        for info_2 in range(info_1 + 1, len(data)):
            vec_2 = list(map(float,data[info_2][1:]))
            #coexp_mat[info_1][info_2] = np.dot(vec_1,vec_2)
            #coexp_mat[info_1][info_2] = pearsonr(np.array(vec_1),np.array(vec_2))
            coexp_mat[info_1][info_2] = np.corrcoef(vec_1,vec_2)[0,1]  
            coexp_mat[info_2][info_1] = np.corrcoef(vec_1,vec_2)[0,1]  
    return coexp_mat
t0 = time.time()
mat = create_coexp_mat(data)
print("coexp",time.time()-t0)
def create_graph_mat(mat,p_min):
    graph_mat = np.zeros((len(data),len(data)))
    for i in range(len(mat)):
        for j in range(i+1,len(mat)):
            print(i,j)
            if(mat[i][j]> p_min):
                graph_mat[i][j]=1
                graph_mat[j][i]=1
            else:
                graph_mat[i][j]=0
                graph_mat[j][i]=0
t0 = time.time()
res = create_graph_mat(mat,p_min=0.8)
print("graph ",time.time()-t0)
'''