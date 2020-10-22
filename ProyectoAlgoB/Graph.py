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
import TopologyReader  as topr

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
    '''
    Create graph with {Atom : [[Atom2,r][Atom3,r]]}
    values are list with pairs Atom2,r
    '''
    nodes = _createNodes( content )
    edges = _createEdges( nodes, cutoff )
    return Graph( nodes, edges )
        
    #print(f"Time of graph building : {time.time()-t0:.4f} s")

def _createNodes( content ):
    '''
    nodes are not nodes per se, they are a dictionary mapping an atom id with its respective atom
    '''
    atoms,bonds  = topr.create_dicts(protein_itp="./Protein_A.itp")
    nodes = {}
    for lines in content:

        aaline1=lines.split()
        if(aaline1[0]!="ATOM" ):
            #print("remarkline")
            # if line are REMARK then pass
            continue
        atom_type_temp=atoms[aaline1[1]][0] ## map atomtype from topology
        A1 = Atom( aaline1[1],aaline1[2],atom_type_temp,aaline1[3],aaline1[5],aaline1[6:9] )
        nodes[ A1.atom_id ] = { A1.atom_id : A1 }
    return nodes

def _createEdges( atoms, cutoff ):
    '''

    '''
    edges = {}
    ids = list( atoms.keys() )
    for i in range( len( ids ) ):
        for j in range( i + 1, len( ids ) ):
            if i != j:
                distance = calcDist( atoms[ ids[i] ][ ids[i] ].coords, atoms[ ids[j] ][ ids[j] ].coords )
                if distance <= cutoff:
                    '''
                    inserts tuple : (id of connected node, distance ) in both nodes
                    '''
                    if ids[i] in edges:
                        edges[ ids[i] ].append( ( ids[j], distance ) )
                    else:
                        edges[ ids[i] ] = [ ( ids[j], distance ) ]

                    if ids[j] in edges:
                        edges[ ids[j] ].append( ( ids[i], distance ) )
                    else:
                        edges[ ids[j] ] = [ ( ids[i], distance ) ]
    return  edges