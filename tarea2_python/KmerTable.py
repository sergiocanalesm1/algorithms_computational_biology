#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Sep  9 18:36:09 2020

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

#######################################################################
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

'''

KmersTable functions. 


'''


def create_KmersTable(data,kmersize):
    '''
    Create a dictionary with key = kmer and value = abundance of kmer
    
    '''
    start_time = time.time()
    kmers = dict()
    for line in data:
        for i in range(len(line)-kmersize+1):
            if (line[i:i+kmersize]) in kmers:
                kmers[line[i:i+kmersize]] += 1
            else: 
                kmers[line[i:i+kmersize]] = 1
    print("creation of kmersTable (ms) : {} ".format(time.time() - start_time))

    return kmers


def list_Kmers(kmerstable):
    '''
    return kmerstable keys = Kmers in a list
    '''
    #start_time = time.time()
    #print("--- %s seconds ---" % (time.time() - start_time))
    return list(kmerstable.keys())

def kmer_abundace(kmer,kmerstable):
    '''
    return kmer abundance for a specfick key.
    
    '''
    return (kmerstable.get(kmer))


def kmer_mean_abundance(kmerstable):
    return (np.mean(kmers_distri(kmerstable)))


def kmers_distri(kmertable):
    '''
    return an array with distribution of kmertable
    '''
    maxval_kmertable=max(kmertable.values())
    distri=np.zeros(maxval_kmertable+1)
    for i in kmertable.values():
        distri[i] += 1
    return distri

def GenerateDistriPlot(filename,kmersize):
    '''
    Generate plot of distribution based on a kmersize and a file
    file should be in the working directory
    '''
    #plt.figure()
    allword_1line = fastaReader2("{}".format(filename))
    KmerTable= create_KmersTable(allword_1line,kmersize)
    distribut = kmers_distri(KmerTable)
    #print("Mean abundance")
    #print(np.mean(distribut))
    plt.plot(np.arange(0,len(distribut)),distribut,label="{}".format(kmersize))
    #plt.axhline(np.mean(distribut),linestyle="-")
    plt.title("# Kmer distribution with kmersize {}".format(kmersize))
    plt.xlabel("Number of Kmers")
    plt.ylabel("Abundance ")
    #plt.yscale("log")
    plt.legend(loc=(1.05,0.35))
        


