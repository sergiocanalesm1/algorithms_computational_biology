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

import numpy as np
import matplotlib.pyplot as plt
import os
#import time
#import TopologyReader as topr
import Graph as G
import Analysis_Networks as antx
import networkx as nx

#from rdkit import Chem
#from rdkit.Chem import ChemicalFeatures
#from rdkit import RDConfig


if __name__ == "__main__":
    #filename = "water_prot"
    
    path="/home/david/Documents/BionIF/Algortimos/Proyecto/MD/MD_dataset/"
    #path="/hpcfs/home/bcom4006/estudiantes/DRFB/Proyecto/MD/MD_Dataset/"
    pdbs = ["1B8E","1A43"]
    
    #file_test='/home/david/Documents/BionIF/Algortimos/Proyecto/MD/MD_dataset/1A43/CA218S/'
    
    with open("results.txt","w") as results_file:
        for  i in pdbs:
            if(len(i)<5):
                results_file.write(i+"\n")
                mutants = os.listdir(path+i)
                for mutval in mutants:
                    #reads the initial pdb to create a dict[resnmae]=atom index
                    # used for locate mutation position in nodes                            
                    if(mutval[-2:]=="WT"): 
                        path_pdb_int=path+i+"/"+mutval+"/"+i+"_CG.pdb"
                    else:
                        path_pdb_int=path+i+"/"+mutval+"/"+i+"_"+mutval+"_CG.pdb"
                    with open(path_pdb_int,"r") as orig_pdb:
                        initial_pdb= orig_pdb.readlines()
                    convert_index={}
                    for line_orig in initial_pdb:
                        line_pdb = line_orig.split()
                        if line_pdb[0] != "ATOM":
                            continue
                        else:
                            if line_pdb[5] in convert_index:
                                convert_index[line_pdb[5]].append(line_pdb[1])
                            else:
                                convert_index[line_pdb[5]]=[line_pdb[1]]
                        
                    files_per_mut = os.listdir(path+i+"/"+mutval)
                    if("Protein_A.itp" and "Protein_B.itp" in files_per_mut):
                        print("avoiding")
                        continue
                   
                    
                    data_energy=np.loadtxt(path+i+"/"+mutval+"/energy.xvg",comments=["@","#"])
                    if(mutval[-2:]=="WT"): # analyze wildtype
                        for no_wt in mutants: # extract pos mut for local density calc.
                            if(no_wt[-2:]!="WT"):
                                mut_coord = no_wt[2]
                                new_aa = no_wt[0]
                                chain = no_wt[1]
                                pos_mut = no_wt[2:-1]
                                old_aa = no_wt[-1]            
                                c_temp=0
                                for files_dyn in files_per_mut:
                                    file_nm = files_dyn.split("_")
                                    if(".pdb" in files_dyn and len(file_nm)>3):
                                        c_temp += 1
                                        if(c_temp ==2):
                                            break
                                        with open(path+i+"/"+mutval+"/"+files_dyn,"r") as pdbfile:
                                            content = pdbfile.readlines()
                                        frame=int(file_nm[1])/3    
                                        path_itp=path+i+"/"+mutval+"/Protein_A.itp" 
                                        G_class = G.Graph( content,path_itp )
                                        distance_graph = G_class.createEdges ( G_class.cartesian, 16 )
                                        lj_graph = G_class.createEdges( G_class.LJ, 15 )
                                        global_density = antx.global_density( lj_graph )
                                        print(global_density)
                                        node_mut = convert_index[pos_mut]
                                        for nodes_to_mut in node_mut:
                                            loca_density = antx.local_density(lj_graph,nodes_to_mut)
                                            results_file.write("WT"+new_aa+pos_mut+old_aa+"\t\t"+str(global_density)+"\t"+str(loca_density)+"\t"+str(data_energy[int(frame)][0])+"\t"+ 
                                                       str(data_energy[int(frame)][1:][0])+"\t"+str(data_energy[int(frame)][1:][1])+"\t"+str(data_energy[int(frame)][1:][2])+"\n")
#                                        loca_density = antx.local_density(lj_graph,node_mut)
                                        print(loca_density)
                                        print("final")
                                        print(data_energy[int(frame)][0],data_energy[int(frame)][1:])
                                        #results_file.write("WT-"+str(new_aa)+str(pos_mut)+str(old_aa)+"\t\t"+str(global_density)+"\t"+str(loca_density)+"\t"+str(data_energy[int(frame)][0])+"\t"+
                                        #           str(data_energy[int(frame)][1:][0])+"\t"+str(data_energy[int(frame)][1:][1])+"\t"+str(data_energy[int(frame)][1:][2])+"\n")
                    else:        #analyzin mutants       
                        c_temp=0
                        for files_dyn in files_per_mut:
                            file_nm = files_dyn.split("_")
                            if(".pdb" in files_dyn and len(file_nm)>3):
                                c_temp += 1
                                if(c_temp ==2):
                                    break
                                frame=int(file_nm[1])/3
                                mut_coord = file_nm[2]
                                new_aa = mut_coord[0]
                                chain = mut_coord[1]
                                pos_mut = mut_coord[2:-1]
                                old_aa = mut_coord[-1]
                                with open(path+i+"/"+mutval+"/"+files_dyn,"r") as pdbfile:
                                    content = pdbfile.readlines()
                                path_itp=path+i+"/"+mutval+"/Protein_A.itp" 
                                G_class = G.Graph( content,path_itp )
                                distance_graph = G_class.createEdges ( G_class.cartesian, 16 )
                                lj_graph = G_class.createEdges( G_class.LJ, 15 )
                                global_density = antx.global_density( lj_graph )
                                print(global_density)
                                node_mut = convert_index[pos_mut]
                                for nodes_to_mut in node_mut:
                                    loca_density = antx.local_density(lj_graph,nodes_to_mut)
                                    results_file.write(mutval+"\t\t"+str(global_density)+"\t"+str(loca_density)+"\t"+str(data_energy[int(frame)][0])+"\t"+ 
                                                       str(data_energy[int(frame)][1:][0])+"\t"+str(data_energy[int(frame)][1:][1])+"\t"+str(data_energy[int(frame)][1:][2])+"\n")
                                print(loca_density)
                                print("final")
                                print(data_energy[int(frame)][0],data_energy[int(frame)][1:])
                             
                                               
                                
                            #results.write(global_density,"\t",loca_density+"\t"+data_energy[int(frame)][0]+"\t"+
                            #              data_energy[int(frame)][1:][0]+"\t"data_energy[int(frame)][1:][1]+"\t"+
                            #              data_energy[int(frame)][1:][2]+"\n")
                    
                #except:
                #    print("no data")
                        
                        
                        #results.append(data_energy[int(frame)][0],data_energy[int(frame)][1:])
    #filename = "back_wter_pro"
    
    #with open("./{}.pdb".format(filename),"r") as pdbfile:
    #    content = pdbfile.readlines()
    #G_class = G.Graph( content )
    #distance_graph = G_class.createEdges ( G_class.cartesian, 16 )
    #lj_graph = G_class.createEdges( G_class.LJ, 15 )
    #water_graph = G_class.createEdges( G_class.water, 15 )
    print()


def analyze( filename, graph ):
    #degrees = antx.nodes_degree( graph, filename )
    global_density = antx.global_density( graph )
    local_densities = [ antx.local_density( graph, node ) for node in graph ]
    # antx.grouping_spectrum( degrees, local_densities )#no son del mismo tamaño, no se va a graficar

#analyze(filename, lj_graph)

def RepresentGrap(graph_data,type_g):
      
    G =nx.Graph()
    for k,v in graph_data.items():
        for bonds in v:
            if(type_g=="bonds"):
                G.add_edge(k, bonds[0], weight=bonds[1])
            else:
                G.add_edge(k, bonds[0], weight=bonds[1][0]+bonds[1][1])
    pos = nx.fruchterman_reingold_layout(G)
    plt.figure(figsize=(10,5))
    nx.draw_networkx_nodes(G, pos, node_size=50)
    nx.draw_networkx_edges(G, pos, width=0.6,arrowsize=10)
    #nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif")
    #labels = nx.get_edge_attributes(G,'weight')
    #nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    plt.axis("off")
    plt.show()
##Holaalskjdh
#RepresentGrap(distance_graph,"bonds")
#RepresentGrap(lj_graph,"Electric")

