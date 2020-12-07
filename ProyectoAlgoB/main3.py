#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 15:30:12 2020

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
import os
import time
#import TopologyReader as topr
import Graph as G
import Analysis_Networks as antx
import networkx as nx

if __name__ == "__main__":
    #path = "/home/david/Documents/BionIF/Algortimos/Proyecto/MD/MD_dataset/"
    path="/hpcfs/home/bcom4006/estudiantes/DRFB/Proyecto/MD/MD_Dataset/"
    print("begin of the script")
    t_ini = time.time()

    #pdbs= os.listdir(path)
    #pdbs=np.loadtxt("../../MD_Dataset/Set_Pro_t1")#,dtype=String)
    with open("../../MD_Dataset/Set_Pro_t1","r") as prots:
        pdbs_t=prots.readlines()
    pdbs=[]
    for name in pdbs_t:
    #print(name[:-1])
        pdbs.append(name[:-1])
    print(pdbs)
    #pdbs = ["1A43"]
    stop_frame = 20
    with open("results_20frame_closenes.txt", "w") as results_file:
        results_file.write("Struct(WT-C-Num-new-Nd) global centrality        " +
                           "local centrality        closeness             frame       E_LJ  " +
                           "      E_C            ET  \n")
        for i in pdbs:
            if(len(i) < 5):
                results_file.write(i + "\n")
                print("analyzing " + i)
                mutants = os.listdir(path + i)
                for mutval in mutants:
                    # reads the initial pdb to create a dict[resnmae]=atom index
                    # used for locate mutation position in nodes
                    if(mutval[-2:] == "WT"):
                        path_pdb_int = path + i + "/" + mutval + "/" + i + "_CG.pdb"
                    else:
                        path_pdb_int = path + i + "/" + mutval + "/" + i + "_" + mutval + "_CG.pdb"
                    with open(path_pdb_int, "r") as orig_pdb:
                        initial_pdb = orig_pdb.readlines()
                    convert_index = {}
                    for line_orig in initial_pdb:
                        line_pdb = line_orig.split()
                        if line_pdb[0] != "ATOM":
                            continue
                        else:
                            if line_pdb[5] in convert_index:
                                convert_index[line_pdb[5]].append(line_pdb[1])
                            else:
                                convert_index[line_pdb[5]] = [line_pdb[1]]
                    files_per_mut = os.listdir(path + i + "/" + mutval)
                    if("Protein_A.itp" and "Protein_B.itp" in files_per_mut):
                        # if  directory contains 2 ipt files continue to other directory.
                        print("avoiding {} for multimeric protein".format(i))
                        continue
                    data_energy = np.loadtxt(
                        path + i + "/" + mutval + "/energy.xvg",
                        comments=["@", "#"])
                    if(mutval[-2:] == "WT"):  # analyze wildtype
                        print("analyzing WT of "+i)
                        mutants_to_eval = []
                        for no_wt in mutants:  # extract pos mut for local density calc.
                            if(no_wt[-2:] != "WT"):
                                mut_coord = no_wt[2]
                                new_aa = no_wt[0]
                                chain = no_wt[1]
                                pos_mut = no_wt[2:-1]
                                mutants_to_eval.append(pos_mut)
                                old_aa = no_wt[-1]
                                c_temp = False
                                lim_frames = 0
                                #print(pos_mut)
                        for frame_indx in range(0, 301, 3):
                            files_dyn = "f_" + \
                                str(frame_indx)+"_"+mutval+"_"+i+".pdb"
                            if(lim_frames == stop_frame):
                                break
                            lim_frames += 1
                            file_nm = files_dyn.split("_")
                            if(".pdb" in files_dyn and len(file_nm) > 3):
                                with open(path + i + "/" + mutval + "/" + files_dyn, "r") as pdbfile:
                                    content = pdbfile.readlines()
                                frame = int(file_nm[1]) / 3
                                path_itp = path + i + "/" + mutval + "/Protein_A.itp"
                                t_g = time.time()
                                G_class = G.Graph(content, path_itp)
                                if(c_temp == False):
                                    print(
                                        f"size of {i} = {len(G_class.content)-2}")
                                    print(
                                        f"Graph (s) | Centrality (all mutan ts {len(mutants_to_eval)})(s)")
                                    c_temp = True
                                prot_end = G_class.prot_end
                               #distance_graph = G_class.createEdges(
                               #    G_class.cartesian, 16)
                                lj_graph = G_class.createEdges(
                                    G_class.LJ, 15)
                                global_density = antx.global_density(
                                    lj_graph)
                                print("  {:.4f}  ".format(
                                    time.time()-t_g), end="|")
                                t_c = time.time()
                                mutation_nodes = []
                                for pos_mutation in mutants_to_eval:
                                    node_mut = convert_index[pos_mutation]
                                    for node_ind in node_mut:
                                        mutation_nodes.append(node_ind)
                                g = nx.DiGraph()
                                g.add_nodes_from(lj_graph.keys())
                                for k, v in lj_graph.items():
                                    for int_list in v:
                                        g.add_edges_from(
                                            [(k, int_list[0])], weight=int_list[1])
                                clust = nx.clustering(g, mutation_nodes)
                                #harm_centr = nx.harmonic_centrality(g,mutation_nodes)
                                for node, clust_coeff in clust.items():
                                    results_file.write("WT" +
                                                       #new_aa +
                                                       pos_mutation +
                                                       #old_aa +
                                                       "-"+node +
                                                       "\t\t" +
                                                       str(global_density) +
                                                       "\t" +
                                                       str(clust_coeff) +
                                                       "\t" +
                                                       str(nx.closeness_centrality(g, node)) +
                                                       "\t" +
                                                       str(data_energy[int(frame)][0]) +
                                                       "\t" +
                                                       str(data_energy[int(frame)][1:][0]) +
                                                       "\t" +
                                                       str(data_energy[int(frame)][1:][1]) +
                                                       "\t" +
                                                       str(data_energy[int(frame)][1:][2]) +
                                                       "\n")
                                    #betwenes = nx.betweenness_subset(g,mutation_nodes,mutation_nodes)
                                    #print(betwenes)
                                print(f"  {time.time()-t_c:.4f} ")

                    else:  # analyzin mutants

                        print("analzying {}".format(mutval))
                        c_temp = False
                        lim_frames = 0
                        for frame_indx in range(0, 301, 3):
                            files_dyn = "f_" + \
                                str(frame_indx)+"_"+mutval+"_"+i+".pdb"
                            #print(file_dyn)
                        #for files_dyn in files_per_mut:
                            if(lim_frames == stop_frame):
                                break
                            lim_frames += 1
                        #for files_dyn in files_per_mut:
                            file_nm = files_dyn.split("_")
                            if(".pdb" in files_dyn and len(file_nm) > 3):

                                frame = int(file_nm[1]) / 3
                                mut_coord = file_nm[2]
                                new_aa = mut_coord[0]
                                chain = mut_coord[1]
                                pos_mut = mut_coord[2:-1]
                                old_aa = mut_coord[-1]
                                with open(path + i + "/" + mutval + "/" + files_dyn, "r") as pdbfile:
                                    content = pdbfile.readlines()
                                path_itp = path + i + "/" + mutval + "/Protein_A.itp"
                                t_g = time.time()
                                G_class = G.Graph(content, path_itp)
                                if(c_temp == False):
                                    print(
                                        f"size of {i} = {len(G_class.content)-2}")
                                    print("Graph (s) |  Centrality(s)")
                                    c_temp = True
                                prot_end = G_class.prot_end
                                lj_graph = G_class.createEdges(G_class.LJ, 15)

                                print("  {:.4f}  ".format(
                                    time.time()-t_g), end="|")

                                t_c = time.time()
                                g = nx.DiGraph()
                                g.add_nodes_from(lj_graph.keys())
                                for k, v in lj_graph.items():
                                    for int_list in v:

                                        g.add_edges_from(
                                            [(k, int_list[0])], weight=int_list[1])

                                global_density = antx.global_density(lj_graph)
                                
                                #harm_centr = nx.harmonic_centrality(g,mutation_nodes)
                                node_mut = convert_index[pos_mut]
                                clust = nx.clustering(g, node_mut)
                                #print(node_mut[0])
                                results_file.write(mutval+"-"+node_mut[0] +
                                                   "\t\t" +
                                                   str(global_density) +
                                                   "\t" +
                                                   str(clust[node_mut[0]]) +
                                                   "\t" +
                                                   str(nx.closeness_centrality(g, node_mut[0])) +
                                                   "\t" +
                                                   str(data_energy[int(frame)][0]) +
                                                   "\t" +
                                                   str(data_energy[int(frame)][1:][0]) +
                                                   "\t" +
                                                   str(data_energy[int(frame)][1:][1]) +
                                                   "\t" +
                                                   str(data_energy[int(frame)][1:][2]) +
                                                   "\n")
                            print(f"  {time.time()-t_c:.4f} ")


