#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 17:32:31 2020

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
import os 
import sys 
import re

def nodes_degree(graph):
    # graficar distribucion de vertices. 
    # parecido a los kmers
    # return array index as degree value at index i is the total number of nodes.

    maxval=max((len(v), k) for k,v in graph.items())
    degree=np.zeros(maxval[0]+1)
    #degrees = np.array()
    for k,v in graph.items():
        print(k,len(v))
        degree[len(v)] += 1
    return degree
    
    
def global_density(graph):
    """
    Parameters
    ----------
    graph : graph
        LJ parameters.

    Returns
    -------
    global_density of graph.

    """
    total_nodes=0
    total_edges=0
    
    for k,v in graph.items():
        total_nodes += 1
        total_edges += len(v)
        
    global_density = 2*total_nodes/(total_edges*(total_edges-1))
    return global_density

def local_density(graph,node):
    """    

    Parameters
    ----------
    graph : TYPE
        DESCRIPTION.
    nodes : TYPE
        DESCRIPTION.

    Returns
    -------
    local_density of give node.

    """
    edges_temp = graph[str(node)]
    count=0
    conexions = [i[0] for i in edges_temp] # nodes conected to node
    for edge in conexions:
        neigh = graph[str(edge)]
        conexion_neig = [i[0] for i in neigh] # nodes conected to node i (conected to node:parameter)
        for neigh_edge in conexion_neig:
            if neigh_edge in conexions:
                count += 1
                
    posible_edges = len(edges_temp)*(len(edges_temp)-1)/2
    return count/posible_edges
        
    
    
    
    
    
    
