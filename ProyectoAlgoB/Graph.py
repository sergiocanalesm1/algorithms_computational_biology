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
from dist import calcDist_cython
from dist import calcLJ_cython


class Graph():
    '''
    The class Graph has no attribute edges, instead, call the function createEdges to obtain the connected graph you wish for
    '''

    def __init__(self, content, protein_itp):
        self.LJ = "LJ"
        self.cartesian = 'CARTESIANA'
        self.water = "WATER"
        self.content = content
        self.protein_itp = protein_itp
        self.atoms, self.bonds = topr.create_dicts(self.protein_itp)
        self.nodes, self.prot_end = self._createNodes()
        self.LJmatrix, self.atomtypes = topr.getLJmatrix()

    def _createNodes(self):
        '''
        nodes are not nodes per se, they are a dictionary mapping an atom id with its respective atom
        '''

        nodes = {}
        convertion = {}
        first_water = False
        prot_end = 0
        for lines in self.content:
            aaline1 = lines.split()
            if (aaline1[0] != "ATOM"):
                # print("remarkline")
                # if line are REMARK then pass
                continue
            if (aaline1[2] == "W"):
                if(first_water == False):
                    prot_end = aaline1[1]
                    first_water = True
                atom_type_temp = "P4"
            else:
                # map atomtype from topology
                atom_type_temp = self.atoms[aaline1[1]][0]
            A1 = Atom(aaline1[1], aaline1[2], atom_type_temp,
                      aaline1[3], aaline1[4], aaline1[6:9])

            nodes[A1.atom_id] = {A1.atom_id: A1}
        return nodes, int(prot_end)

    def createEdges(self, g_type, cutoff):
        '''
        creates edges based on the graph type
        '''
        edges = {}
        weight = 0
        ids = list(self.nodes.keys())
        # if cartesian, then create bonds based on bons dictionary from
        # topology
        if g_type == self.cartesian:
            for k, v in self.bonds.items():
                # print(k,v)
                for con_a in v:
                    # k=int(k)
                    # con_a=int(con_a)
                    #nodo_1 = self.nodes[ids[k]][ids[k]]
                    #nodo_2 = self.nodes[ids[con_a]][ids[con_a]]

                    weight = self._calcDist(
                        self.nodes[k][k].coords,
                        self.nodes[con_a][con_a].coords)
                    if k in edges:
                        edges[k].append((con_a, weight))
                    else:
                        edges[k] = [(con_a, weight)]

                    if con_a in edges:
                        edges[con_a].append((k, weight))
                    else:
                        edges[con_a] = [(k, weight)]
        water = False
        if g_type in [self.LJ, self.water]:
            for i in range(self.prot_end + 1):
                #water = self.nodes[ ids[ i ] ][ ids[ i ] ].atom_type == "W"
                # if not water:
                for j in range(i + 1, len(ids)):
                    if i != j:
                        nodo_1 = self.nodes[ids[i]][ids[i]]
                        nodo_2 = self.nodes[ids[j]][ids[j]]
                        if(nodo_1.atom_name == "W" and nodo_2.atom_name == "W"):
                            continue
                        x_n1 = list(map(float, nodo_1.coords))
                        x_n2 = list(map(float, nodo_2.coords))
                        dist = calcDist_cython(
                            x_n1[0], x_n1[1], x_n1[2], x_n2[0], x_n2[1], x_n2[2])
                        #dist = self._calcDist(self.nodes[ids[i]][ids[i]].coords, self.nodes[ids[j]][ids[j]].coords)
                        if 3 <= dist <= cutoff:
                            atom1 = self.nodes[ids[i]][ids[i]]
                            atom2 = self.nodes[ids[j]][ids[j]]
                            sigmaij, epsilonij = self.LJmatrix[self.atomtypes[atom1.atom_type]
                                                               [1]][self.atomtypes[atom2.atom_type][1]]
                            #lj_w=self._calcLJ( self.nodes[ids[i]][ids[i]], self.nodes[ids[j]][ids[j]] ,dist)
                            lj_w = calcLJ_cython(epsilonij, sigmaij, dist)
                            if(nodo_1.atom_type == "Qd" and nodo_2.atom_type == "Qd"):
                                q_w = self._calcCoulb(
                                    self.atoms[ids[i]][1], self.atoms[ids[j]][1], dist)
                            else:
                                q_w = 0
                            energy = (lj_w, q_w)
                            '''
                            inserts tuple : ( id of connected node, distance ) in source and destination node
                            if water, just from source to destination
                            '''
                            if ids[i] in edges:
                                edges[ids[i]].append((ids[j], energy))
                            else:
                                edges[ids[i]] = [(ids[j], energy)]
                            if nodo_2.atom_name != "W":
                                if ids[j] in edges:
                                    edges[ids[j]].append((ids[i], energy))
                                else:
                                    edges[ids[j]] = [(ids[i], energy)]

                        '''
                        else:
                            #TODO calculate equivalent of energy in water
                            lj_w=self._calcLJ( self.nodes[ids[i]][ids[i]], self.nodes[ids[j]][ids[j]] )
                            q_w=self._calcCoulb(self.atoms[ids[i]][1],self.atoms[ids[i]][1],dist)
                            energy = (lj_w,q_w)
                        '''
                        # if dist <= cutoff :
                        '''
                        inserts tuple : ( id of connected node, distance ) in source and destination node
                        if water, just from source to destination
                        '''
                        #    if ids[i] in edges:
                        #        edges[ids[i]].append((ids[j], energy))
                        #    else:
                        #        edges[ids[i]] = [(ids[j], energy)]
                        #    if not water:
                        #        if ids[j] in edges:
                        #            edges[ids[j]].append((ids[i], energy))
                        #        else:
                        #            edges[ids[j]] = [(ids[i], energy)]

        return edges

    def _calcDist(self, coord1, coord2):
        '''
        Calculate distance between two coordinates
        coord1, coord2 = tuples (x,y,z)

        '''
        xs = (float(coord2[0]) - float(coord1[0]))
        ys = (float(coord2[1]) - float(coord1[1]))
        zs = (float(coord2[2]) - float(coord1[2]))
        return np.sqrt(xs ** 2 + ys ** 2 + zs ** 2)

    # print(f"Time of graph building : {time.time()-t0:.4f} s")

    def _calcLJ(self, atom1, atom2, dist):
        # TODO check order of sigma and epsilon
        sigmaij, epsilonij = self.LJmatrix[self.atomtypes[atom1.atom_type][
            1]][self.atomtypes[atom2.atom_type][1]]
        #rij = self._calcDist( atom1.coords, atom2.coords )
        rij = dist
        return 4 * epsilonij * ((sigmaij / rij)**12 - (sigmaij / rij)**6)

    def _calcCoulb(self, q1, q2, dist):
        return float(q1) * float(q2) / (4 * np.pi * 15 * dist)
