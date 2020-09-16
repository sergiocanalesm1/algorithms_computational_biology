#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 19:22:29 2020

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
def fastaReader2(filename):
    '''
    return a list with all lines of fasta file
    '''
    with open("/home/david/Documents/BionIF/Algortimos/algorithms_computational_biology/tarea2_python/{}".format(filename),"r") as test1:
        text_pure =test1.readlines()  
    reads=[]
    for i in range(1,len(text_pure),4):
        reads.append(text_pure[i][:-1])
    return reads

def find_overlap(source ,destination ):
    i = 0;
    j = 0;
    if(source == destination):
        return 0
    while ( i < len(source) ):
        if ( source[ i ] == destination[ j ] ):
            j += 1
            i += 1
        else :
            i = i - j + 1
            j = 0
    return j


a1=' Muy pronto ha de sobrarnos oro para empedrar la casa replico su marido. Durante varios meses se emp'
b='os oro para empedrar la casa replico su marido. Durante varios meses se empenho en demostrar el acie'
c="dirlo. Muy pronto ha de sobrarnos oro para empedrar la casa replico su marido. Durante varios meses"
num=find_overlap(c,a1)
print(num)

reads = fastaReader2("test100x.fastq")
'''
for i in reads:
    if("sobrarnos" in i[:40] ):
        print(i)
'''

#%%
#def ProccessRead(minOverlap=10):
minOverlap=10
readCounts=dict()
father=dict()
overlaps=dict()
for line in (reads):
    if("Muchos" in line):
        print(line)
    if(line not in readCounts): # seq no existe
        readCounts[line] = 1
        father[line] = True
    else: # seq ya existe
        readCounts[line] += 1
        
    
    #paso 2. Encontrar los hijos de line que ya esten guardados

for line in reads:   
    Prefixes=[]
    for seq_inreadC in readCounts.keys():
        overlap_new_saved=find_overlap(line,seq_inreadC)
        if(overlap_new_saved>minOverlap): # si el overlap es mayor a min
            Prefixes.append((seq_inreadC,overlap_new_saved))# agregue
    overlaps[line]=Prefixes
    
    
    for seq_inOverlap in overlaps.keys():
        overlap_save_new=find_overlap(seq_inOverlap,line)
        if(overlap_save_new> minOverlap): # line es hijo de alguien
            overlaps[seq_inOverlap].append((line,overlap_save_new))
            if(" Muy pronto ha" in line[:14]):
                print(line)
            if(father[line]==False):
                pass
            else:
                father[line]=False
    
        
#ProccessRead()
    
        
    





#%%


MaxSuc = 0 
distriovelap = np.zeros(200)
for i in overlaps.values():
    distriovelap[len(i)] += 1
print(distriovelap)






def proccess_read_to_overlap(kmertable,minOverlap):
    '''
    Generate 2 hash maps that containns readCounts and overlap graph 

    '''
    sequence = list(kmertable.keys())
    readCounts = dict()
    overlaps = dict()
    for seq in sequence:

        '''
        add new sequence in the readcount
        if seq existence
        '''
        if (seq not in readCounts):
            readCounts[seq] = 1 
        else:
            readCounts[seq] += 1
            continue
    '''
    create list with (sq_source,sq_dest,val) wiht sq_source as new seq
    load the dictionary
    '''
    for seq_rc in readCounts.keys():
        overlaplist=[]
        for i_rc in readCounts.keys():
            if(i_rc==seq_rc):
                pass
            else:       
                val=find_overlap_lengt(seq_rc,i_rc)
                if(val > minOverlap):
                    overlaplist.append((seq_rc,i_rc,val))
        overlaps[seq_rc]=(overlaplist)
    '''
    complete overlap dict() with cases with sqe_dest as new seq

    '''
    for i in overlaps.keys():
        inv_overlap_lis=[]
        for k in overlaps.keys():
            if(i==k):
                pass
            else: 
                val= find_overlap_lengt(k,i)
                if(val>minOverlap):
                    inv_overlap_lis.append((k,i,val))
        if(len(inv_overlap_lis)>minOverlap):
            for invovs in inv_overlap_lis:
                overlaps[i].append(invovs)

    return readCounts,overlaps


def overlap_distribution(overlaps):
    '''
    return an array with distribution of kernel 
    '''
    values=np.zeros(len(overlaps.values()))
    for i in overlaps.values():
        values[len(i)] += 1 
    return values