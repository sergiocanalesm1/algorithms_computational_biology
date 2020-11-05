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

#import numpy as np
import matplotlib.pyplot as plt
#import os
#import time
#import TopologyReader as topr
import Graph as G
import Analysis_Networks as antx
#import networkx as nx

#from rdkit import Chem
#from rdkit.Chem import ChemicalFeatures
#from rdkit import RDConfig
#%%

if __name__ == "__main__":
    filename = "1UBQ-CG"
    with open("./{}.pdb".format(filename),"r") as pdbfile:
        content = pdbfile.readlines()



    G_class = G.Graph( content )
    distance_graph = G_class.createEdges ( G_class.cartesian, 16 )
    lj_graph = G_class.createEdges( G_class.LJ, 15 )

    degrees = antx.nodes_degree( lj_graph, filename )
    global_density = antx.global_density( lj_graph )
    local_densities = [ antx.local_density( lj_graph, node ) for node in lj_graph ]
    #antx.grouping_spectrum( degrees, local_densities )#no son del mismo tamaño, no se va a graficar


'''    
def RepresentGrap(graph_data):
    G =nx.Graph()
    for k,v in graph_data.items():
        for bonds in v:
            G.add_edge(k, bonds[0], weight=bonds[1])
    pos = nx.fruchterman_reingold_layout(G)
    plt.figure(figsize=(10,5))
    nx.draw_networkx_nodes(G, pos, node_size=50)
    nx.draw_networkx_edges(G, pos, width=0.1,arrowsize=10)
    #nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif")
    #labels = nx.get_edge_attributes(G,'weight')
    #nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    plt.axis("off")
    plt.show()

#RepresentGrap(distance_graph)
#RepresentGrap(lj_graph)
'''
