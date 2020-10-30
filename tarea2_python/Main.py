#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Sep  9 18:42:00 2020

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

# initial modules and plots parameters 
import numpy as np 
import matplotlib.pyplot as plt
import time

Fontsize=16
AX=16
TICKS=14
LEGEND=16
TITLES=18

plt.rc('lines',linewidth=0.8)
plt.rc('font', size=Fontsize)          # controls default text sizes
plt.rc('axes', titlesize=AX)     # fontsize of the axes title
plt.rc('axes', labelsize=AX)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=TICKS)    # fontsize of the tick labels
plt.rc('ytick', labelsize=TICKS)    # fontsize of the tick labels
plt.rc('legend', fontsize=LEGEND)    # legend fontsize
plt.rc('figure', titlesize=TITLES)
plt.rc('xtick.major', pad=7) 
#######################################################################

from KmerTable import *
from SimulateReads import SimulateReads
from OverlapGraph import *
import networkx as nx



def fastaReader2(filename):
    '''
    return a list with all lines of fasta file
    '''
    with open("./{}".format(filename),"r") as test1:
        text_pure =test1.readlines()  
    reads=[]
    for i in range(1,len(text_pure),4):
        reads.append(text_pure[i][:-1])
    return reads



files_fastas = ["test100x.fastq","test10x.fastq","test20x.fastq","test50x.fastq"]
#files_fastas = [test100x.fastq.gz,test10x.fastq,test20x.fastq,test50x.fastq.gz]
kmers_sizes_vales=[5,10,15,20,50,75]
colors=["b","g","r","c","k","a"]



   


'''
for name in files_fastas:
    print(name)
    for kmersize in kmers_sizes_vales:
        print(kmersize)
        lineword=fastaReader2(name)
        kmertab= create_KmersTable(lineword,kmersize)
        #distribution=kmers_distri(kmertab)
        a=len(list_Kmers(kmertab))
        print("number of kmers = ",a)
        print("mean abundance =",kmer_mean_abundance(kmertab))
    print()        
    '''


