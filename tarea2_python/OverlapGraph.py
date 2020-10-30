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
import networkx as nx
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
    '''
    
    return overlap
    '''
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


#%%
def ProccessReadOverlap(reads,minOverlap=10):
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
    
       
        Prefixes=[]
        for seq_inreadC in readCounts.keys():
            overlap_new_saved=find_overlap(line,seq_inreadC)
            if(overlap_new_saved>minOverlap): # si el overlap es mayor a min
                Prefixes.append((line,seq_inreadC,overlap_new_saved))# agregue
                father[seq_inreadC]=False
        overlaps[line]=Prefixes
        
        #paso 3. Encontrar si line es hijo de alguien guardado
        for seq_inOverlap in overlaps.keys():
            overlap_save_new=find_overlap(seq_inOverlap,line)
            if(overlap_save_new> minOverlap): # line es hijo de alguien
                overlaps[seq_inOverlap].append((seq_inOverlap,line,overlap_save_new))
                if(father[line]==False):
                    pass
                else:
                    father[line]=False
    return readCounts,overlaps,father
            
def getDistintcSequeces(readcounts):
    return readcounts.keys()

def overlap_distribution(overlaps):
    '''
    return an array with distribution of kernel 
    '''
    values=np.zeros(len(overlaps.values()))
    for i in overlaps.values():
        values[len(i)] += 1 
    return values

def getSource(father):
    for key,val in father.items():
        if(val==True):
            return key
def getAssembly(overlap,father):
    visitedseq=[]
    initial = getSource(father)
    complete=initial
    InitialOverlap=overlap[initial]
    Control=True
    while(Control==True):
        maxOverlap=0
        for i in (InitialOverlap):
            tempMaxOv=i[2]
            if(tempMaxOv>maxOverlap):
                maxOverlap=tempMaxOv
                Next=i
        if(len(Next)==0):
            break
        newSeq=Next[1]
        overlNex=Next[2]
        finalSeqnoOv=newSeq[overlNex:]
        if(newSeq in visitedseq):
            break
        visitedseq.append(newSeq)
        complete=complete+finalSeqnoOv
        InitialOverlap = overlap[newSeq]
    return complete    


def PlotGraph(overlap):
    G =nx.DiGraph()
    for i in overlap.values():
        if(len(i)==0):
            continue
        if(len(i)>1):
            for j in i:
                G.add_edge(j[0],j[1],weight=j[2])
        else:
            G.add_edge(i[0][0],i[0][1],weight=i[0][2])
            
    pos = nx.circular_layout(G)
    plt.figure(figsize=(10,5))
    nx.draw_networkx_nodes(G, pos, node_size=500)
    nx.draw_networkx_edges(G, pos, width=1,arrowsize=20)
    nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif")
    labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    plt.axis("off")
    plt.show()

#data10x=fastaReader2("./test10x.fastq")
#readC,overlap,father=ProccessReadOverlap(data10x,40)  
#print(getAssembly(overlap,father))
    
    
    
    
    




