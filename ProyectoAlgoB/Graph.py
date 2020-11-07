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
#import time
import numpy as np
from Atom import Atom
import TopologyReader as topr


class Graph():
    '''
    The class Graph has no attribute edges, instead, call the function createEdges to obtain the connected graph you wish for
    '''

    def __init__( self, content ):
        self.LJ = "LJ"
        self.cartesian = 'CARTESIANA'
        self.water = "WATER"
        self.content = content
        self.atoms, self.bonds = topr.create_dicts( protein_itp="./Protein_A.itp" )
        self.atoms["W"] = 1
        self.nodes = self._createNodes()
        self.LJmatrix, self.atomtypes = topr.getLJmatrix()


    def _createNodes( self ):
        '''
        nodes are not nodes per se, they are a dictionary mapping an atom id with its respective atom
        '''

        nodes = {}
        for lines in self.content:
            aaline1 = lines.split()
            if (aaline1[0] != "ATOM"):
                # print("remarkline")
                # if line are REMARK then pass
                continue
            if aaline1[2] == 'W':
                atom_type_temp = "P4"
                #deberiamos agregar esto al map? de self.atoms?
            else:
                atom_type_temp = self.atoms[aaline1[1]][0]  ## map atomtype from topology
            A1 = Atom(aaline1[1], aaline1[2], atom_type_temp, aaline1[3], aaline1[5], aaline1[6:9])
            nodes[A1.atom_id] = {A1.atom_id: A1}
        return nodes

    def createEdges( self, g_type, cutoff ):

        '''
        creates edges based on the graph type
        '''
        edges = {}
        weight = 0
        ids = list( self.nodes.keys() )
        # if cartesian, then create bonds based on bons dictionary from topology
        if g_type == self.cartesian:
            for k,v in self.bonds.items():
                #print(k,v)
                for con_a in v:
                     #k=int(k)
                     #con_a=int(con_a)
                     weight = self._calcDist(self.nodes[k][k].coords, self.nodes[con_a][con_a].coords)
                     if k in edges:
                         edges[k].append((con_a, weight))
                     else:
                         edges[k] = [(con_a, weight)]
    
                     if con_a in edges:
                        edges[con_a].append((k, weight))
                     else:
                        edges[con_a] = [(k, weight)]
        if g_type in [ self.LJ, self.water ]:
            for i in range(len(ids)):
                if self.nodes[ ids[ i ] ][ ids[ i ] ].atom_type == "W":
                    break
                for j in range(i + 1, len(ids)):
                    if ids[i] == 163:
                        print("hola")
                    if i != j:
                        dist = self._calcDist(self.nodes[ids[i]][ids[i]].coords, self.nodes[ids[j]][ids[j]].coords)

                        water =  self.nodes[ ids[ j ] ][ ids[ j ] ].atom_type == "W"
                        if not water:
                            lj_w = self._calcLJ(self.nodes[ids[i]][ids[i]], self.nodes[ids[j]][ids[j]])
                            q_w=self._calcCoulb(self.atoms[ids[i]][1],self.atoms[ids[i]][1],dist)
                            energy = (lj_w,q_w)
                        else:
                            #TODO calculate equivalent of energy in water
                            energy = ()
                        if dist <= cutoff :
                            '''
                            inserts tuple : ( id of connected node, distance ) in source and destination node
                            if water, just from source to destination
                            '''
                            if ids[i] in edges:
                                edges[ids[i]].append((ids[j], energy))
                            else:
                                edges[ids[i]] = [(ids[j], energy)]
                            if not water:
                                if ids[j] in edges:
                                    edges[ids[j]].append((ids[i], energy))
                                else:
                                    edges[ids[j]] = [(ids[i], energy)]

        return edges

    def _calcDist(self, coord1, coord2):
        '''
        Calculate distance between two coordinates
        coord1, coord2 = tuples (x,y,z)

        '''
        xs = (float(coord2[0]) - float(coord1[0]))
        ys = (float(coord2[1]) - float(coord1[1]))
        zs = (float(coord2[2]) - float(coord1[2]))
        return  np.sqrt(xs ** 2 + ys ** 2 + zs ** 2)

    # print(f"Time of graph building : {time.time()-t0:.4f} s")

    def _calcLJ( self, atom1 , atom2 ):
        #TODO check order of sigma and epsilon
        sigmaij, epsilonij = self.LJmatrix[ self.atomtypes[atom1.atom_type][1] ][ self.atomtypes[atom2.atom_type][1] ]
        rij = self._calcDist( atom1.coords, atom2.coords )
        return 4 * epsilonij * ( ( sigmaij / rij)**12 - ( sigmaij / rij)**6 )
    
    def _calcCoulb(self , q1 , q2 , dist):
        return float(q1)*float(q2)/(4*np.pi*15*dist)