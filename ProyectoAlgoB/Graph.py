#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:10:48 2020

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
import time 
import numpy as np
from Atom import Atom
import TopologyReader as topr

class Graph():
    def __init__(self,nodes,edges):
        self.nodes = nodes
        self.edges = edges
        
def calcDist(coord1,coord2):
    '''
    Calculate distance between two coordinates 
    coord1, coord2 = tuples (x,y,z)
    
    '''
    xs=(float(coord2[0])-float(coord1[0]))
    ys=(float(coord2[1])-float(coord1[1]))
    zs=(float(coord2[2])-float(coord1[2]))
    return( np.sqrt(xs**2 + ys**2 + zs**2)  )


def createGraph(content,cutoff):
    t0 = time.time()
    nodes=dict()
    overlaps=dict()
    for lines in content:
        '''
        Create graph with {Atom : [[Atom2,r][Atom3,r]]}
        values are list with pairs Atom2,r
        '''
        aaline1=lines.split()
        if(aaline1[0]!="ATOM" ):
            #print("remarkline")
            # if line are REMARK then pass
            continue
        interactions=[]
        A1=Atom(aaline1[2],aaline1[3],aaline1[4],aaline1[5:8])

        for otherlines in content:
            aaline2=otherlines.split()
            if(aaline1 == aaline2):
                # if line are the same then pass
                continue
            if(aaline1[0]!="ATOM" or aaline2[0]!="ATOM"  ):
                #print("remarkline")
                # if line are REMARK then pass
                continue
            c1=aaline1[6:9]
            c2=aaline2[6:9]
            dist_aa1_aa2=calcDist(c1,c2)
            # crear atom con (atomname,atom_id,res_name,res_num,coords)
            A2=Atom(aaline2[2],aaline2[3],aaline2[5],aaline2[6:9])
            if(dist_aa1_aa2<cutoff): ## nodos conectados 
                interactions.append([A2,dist_aa1_aa2])
                #overlaps[A1]=[A2,dist_aa1_aa2]
        overlaps[A1]=interactions
    return overlaps
        
    print(f"Time of graph building : {time.time()-t0:.4f} s")
        
            