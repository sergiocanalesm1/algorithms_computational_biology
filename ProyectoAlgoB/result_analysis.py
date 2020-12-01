import numpy as np
import re
import matplotlib.pyplot as plt


def read_real_results( path ):
    '''
    result example:
    {
    '1A43': {
             '156': -2.4,
             '159': -4.55,
             '167': -4.55,
             '184': -0.7,
             '218': -3.7
             }
    }
    '''
    protein_stabilities = {}
    with open( path ) as file_data:
        lines = file_data.readlines()
    for line in lines:
        line = line.split()
        # que hacer con el archivo?
        protein = line[ 0 ].split( "." )[ 0 ]  # example: from '1A43.pdb' to 1A43
        mutant_name = line[ 2 ]  # example GA156A
        stability = float( line[ 3 ] )  # example -2.40
        if protein not in protein_stabilities:
            protein_stabilities[ protein ] = { re.sub( "[^0-9]", "", mutant_name ) : stability } # example: re.sub converts from WT156 to 156
        else:
            protein_stabilities[ protein ][ re.sub( "[^0-9]", "", mutant_name ) ] = stability
    return protein_stabilities


def read_foldx_results(path):
    '''
    result example:
    {
    '1A43': {
             '156': 0.262397,
             '159': 0.557025,
             '167': 0.077106,
             '184': 0.43004,
             '218': 3.6956
             }
    }
    '''
    protein_stabilities = {}
    with open( path ) as file_data:
        lines = file_data.readlines()
    # TODO desde qué linea leer!!!!!!
    for line in lines:
        line = line.split()
        protein = line[ 1 ].split( "_" )[ 0 ]  # example: from 1A43_1_0.pdb to 1A43.pdb  to 1A43
        mutant_name = line[ 0 ]  # example GA156A
        stability = float( line[ 2 ] )  # example 0.262397
        if protein not in protein_stabilities:
            protein_stabilities[ protein ] = { re.sub( "[^0-9]", "", mutant_name ) : stability } # example: re.sub converts from WT156 to 156
        else:
            protein_stabilities[ protein ][ re.sub( "[^0-9]", "", mutant_name ) ] = stability
    return protein_stabilities

def setup_proteins( path ):

    wt_string = "wild_types"
    mutant_string = "mutants"

    with open( path ) as file:
        lines = file.readlines()
    proteins = {}
    for line in lines[ 1: ]:
        if len( line ) < 10:
            protein_name = line.replace("\n", "")
            proteins[ protein_name ] = { wt_string : {}, mutant_string : {} }
            continue
        line = line.split()
        node = line[ 0 ]
        local_centrality = float( line[ 2 ] )
        frame = float( line[ 3 ] )
        energy = float( line[ 6 ] )

        type = mutant_string
        if node.startswith( "WT" ):
            type = wt_string
        prot_name = node.split( "-" )[ 0 ] #example: from  GA156A-18 to GA156A
        if prot_name not in proteins[ protein_name ][ type ]:
            proteins[ protein_name ][ type ][ prot_name ] = {
                "max_energy": energy,
                "min_energy": energy,
                "max_frame": frame,
                node : {
                    "local_centrality" : local_centrality,
                    "frame" : frame,
                    "energy" : energy
                }
            }
        #TODO : hacer que solo obtenga el primer frame de cada nodo para compara frame 0 con frame 57.
        else:
            proteins[ protein_name ][ type ][ prot_name ][ node ] = {
                "local_centrality": local_centrality,
                "frame": frame,
                "energy": energy
            }
            current_max_energy = proteins[ protein_name ][ type ][ prot_name ][ "max_energy" ]
            current_min_energy = proteins[ protein_name ][ type ][ prot_name ][ "min_energy" ]
            current_max_frame = proteins[ protein_name ][ type ][ prot_name ][ "max_frame" ]
            if frame > current_max_frame :
                proteins[ protein_name ][ type ][ prot_name ][ "max_frame"] = frame
            if energy < current_max_energy : #la que se convierte en 100%
                proteins[ protein_name ][ type ][ prot_name ][ "max_energy"] = energy
            elif energy > current_min_energy : #la que se convierte en 60%
                proteins[ protein_name ][ type ][ prot_name ][ "min_energy" ] = energy
    return proteins

def setup_proteins_all_frames( path ):

    wt_string = "wild_types"
    mutant_string = "mutants"

    with open( path ) as file:
        lines = file.readlines()
    proteins = {}
    for line in lines[ 1: ]:
        if len( line ) < 10:
            protein_name = line.replace("\n", "")
            proteins[ protein_name ] = { wt_string : {}, mutant_string : {} }
            continue
        line = line.split()
        node = line[ 0 ]
        local_centrality = float( line[ 2 ] )
        frame = float( line[ 3 ] )
        energy = float( line[ 6 ] )

        type = mutant_string
        if node.startswith( "WT" ):
            type = wt_string
        prot_name = node.split( "-" )[ 0 ] #example: from  GA156A-18 to GA156A
        node_index = node.split( "-" )[ 1 ]
        if prot_name not in proteins[ protein_name ][ type ]:
            general_info = {
                "max_energy": energy,
                "min_energy": energy,
                "max_frame": frame,
                }
            node = {
                    "local_centrality" : local_centrality,
                    "frame" : frame,
                    "energy" : energy,
                    "node_index" : node_index
                }
            proteins[ protein_name ][ type ][ prot_name ] = [
                general_info,node 
            ]
        else:
            node = {
                    "local_centrality" : local_centrality,
                    "frame" : frame,
                    "energy" : energy,
                     "node_index" : node_index
                }
            proteins[ protein_name ][ type ][ prot_name ].append(node)
            current_max_energy = proteins[ protein_name ][ type ][ prot_name ][0][ "max_energy" ]
            current_min_energy = proteins[ protein_name ][ type ][ prot_name ][0][ "min_energy" ]
            current_max_frame = proteins[ protein_name ][ type ][ prot_name ][0][ "max_frame" ]
            if frame > current_max_frame :
                proteins[ protein_name ][ type ][ prot_name ][0][ "max_frame"] = frame
            if energy < current_max_energy : #la que se convierte en 100%
                proteins[ protein_name ][ type ][ prot_name ][0][ "max_energy"] = energy
            elif energy > current_min_energy : #la que se convierte en 60%
                proteins[ protein_name ][ type ][ prot_name ][0][ "min_energy" ] = energy
    return proteins


def normalize_nodes( dict_with_nodes ):
    '''
    example result:
    {
      '156': 0.3939139505890136,
      '184': 0.45708672450167287,
      '218': 0.44125935162094765,
      '167': 0.3960422035548457,
      '159': 0.5209813999353262
    }
    '''
    proteins = {}
    prots = list( dict_with_nodes.keys() )
    for protein in prots :
        #max = dict_with_nodes[ protein ][ "max_energy" ]
        #min = dict_with_nodes[ protein ][ "min_energy" ]
        acumulated_centrality = 0
        total_nodes = 0
        for node in dict_with_nodes[ protein ] :
            if node not in [ "max_energy", "min_energy", "max_frame" ]:
                acumulated_centrality += dict_with_nodes[ protein ][ node ][ "local_centrality" ] #* local_centrality_percentage( max, min, dict_with_nodes[ protein ][ node ][ "energy" ] )
                total_nodes += 1
        proteins[ re.sub( "[^0-9]", "", protein ) ] =  acumulated_centrality / total_nodes # example: re.sub converts from WT156 to 156
    return proteins


def normalize_nodes_all_frames( dict_with_nodes ):
    '''
    this is a test to normalize all nodes
    example result:
    {
      '156': 0.3939139505890136,
      '184': 0.45708672450167287,
      '218': 0.44125935162094765,
      '167': 0.3960422035548457,
      '159': 0.5209813999353262
    }
    '''
    proteins = {}
    prots = list( dict_with_nodes.keys() )
    list_info=[ "max_energy", "min_energy", "max_frame" ]
    total_frame = 0
    for protein in prots :
        max = dict_with_nodes[ protein ][0][ "max_energy" ]
        min = dict_with_nodes[ protein ][0][ "min_energy" ]
        acumulated_centrality = 0
        total_nodes = 0
        
        acumulated_centrality_dict = {}
        for node in dict_with_nodes[ protein ] :
            #if node.co not in [ "max_energy", "min_energy", "max_frame" ]:
            if not any([i in node for i in list_info]):
                #acumulated_centrality +=  node[ "local_centrality" ] * local_centrality_percentage( max, min, node[ "energy" ] )
                if node['node_index'] not in acumulated_centrality_dict:
                    acumulated_centrality_dict[node["node_index"]] = node[ "local_centrality" ] * local_centrality_percentage( max, min, node[ "energy" ] )
                    total_nodes += 1
                    total_frame += 1
                else:
                    acumulated_centrality_dict[node["node_index"]] += node[ "local_centrality" ] * local_centrality_percentage( max, min, node[ "energy" ] )
                    total_frame += 1
        total_frame = total_frame/total_nodes 
 
        acumulated_centrality = 0
        for k in acumulated_centrality_dict.values():
            acumulated_centrality += k
        acumulated_centrality = acumulated_centrality/total_frame
        proteins[ re.sub( "[^0-9]", "", protein ) ] =  acumulated_centrality / total_nodes # example: re.sub converts from WT156 to 156
    return proteins

def read_graph_results( path ):
    graph_stability = {}
    proteins = setup_proteins( path )
    wt_string = "wild_types"
    mutant_string = "mutants"
    for protein_name, types in proteins.items():
        graph_stability[ protein_name ] = {}
        wild_types_dict_with_nodes = types[ wt_string ]
        mutants_dict_with_nodes = types[ mutant_string ]
        wild_types = normalize_nodes( wild_types_dict_with_nodes )
        mutants = normalize_nodes( mutants_dict_with_nodes )
        for protein_number, mutant_centrality in mutants.items():
            if protein_number in wild_types :
                graph_stability[ protein_name ][ protein_number ] = mutant_centrality - wild_types[ protein_number ]
            else:
                print( "graph error: la mutante no se encuentra en el WT" )
    return graph_stability



def read_graph_results_all_frames( path ):
    graph_stability = {}
    proteins = setup_proteins_all_frames( path )
    wt_string = "wild_types"
    mutant_string = "mutants"
    for protein_name, types in proteins.items():
        graph_stability[ protein_name ] = {}
        wild_types_dict_with_nodes = types[ wt_string ]
        mutants_dict_with_nodes = types[ mutant_string ]
        wild_types = normalize_nodes_all_frames( wild_types_dict_with_nodes )
        mutants = normalize_nodes_all_frames( mutants_dict_with_nodes )
        for protein_number, mutant_centrality in mutants.items():
            if protein_number in wild_types :
                graph_stability[ protein_name ][ protein_number ] = mutant_centrality - wild_types[ protein_number ]
            else:
                print( "graph error: la mutante no se encuentra en el WT" )
    return graph_stability

def local_centrality_percentage( max_energy, min_energy, energy, min_percentage=0.6 ) :
    max = np.abs( max_energy )
    min = np.abs( min_energy )
    x = np.abs( energy )
    return ( ( x - min ) / ( max - min ) ) * ( 1 - min_percentage ) + min_percentage

def create_confusion_matrix( real_data_results, data_results ) :
    true_positives = 0
    true_negatives = 0
    false_negatives = 0
    false_positives = 0

    for protein, v in data_results.items() :
        if protein in real_data_results :
            mutants = list( real_data_results[ protein ].keys() )
            #print(protein)
            for mutant in mutants:
                #print(mutant)
                if mutant in real_data_results[ protein ]:
                    #print("real",real_data_results[ protein ][ mutant ],"data",data_results[ protein ][ mutant ])
                    if real_data_results[ protein ][ mutant ] > 0 and data_results[ protein ][ mutant ] > 0:
                        true_positives += 1
                        #print("TP")
                    elif real_data_results[ protein ][ mutant ] < 0 and data_results[ protein ][ mutant ] < 0:
                        true_negatives += 1
                        #print("TN")
                    elif real_data_results[ protein ][ mutant ] > 0 and data_results[ protein ][ mutant ] < 0:
                        false_negatives += 1
                        #print("FN")
                    elif real_data_results[ protein ][ mutant ] < 0 and data_results[ protein ][ mutant ] > 0:
                        false_positives += 1
                        #print("FP")
                    else:
                        print( "al menos una estabilidad es 0, mutante: ", mutant)
        else:
            print( "create confusion matrix error: protein not in real data set" )
    return true_positives, false_positives, false_negatives, true_negatives

def print_matrix( name, true_positives, false_positives, false_negatives, true_negatives ):
    print( "{} Confusion Matrix \nTrue Positives: {}, False Negatives: {} \nFalse Positives: {}, True Negatives: {}".format( name, true_positives, false_positives, false_negatives, true_negatives ) )

def bar_graph( graph_results, foldx_results , title ):
    labels = [ "True\nPositives", "False\nNegatives", "False\nPositives", "True\nNegatives" ]
    x = np.arange( len( labels ) )
    width = 0.4
    fig, ax = plt.subplots()
    ax.bar( x - width / 2, graph_results, width, label = 'Graph' )
    ax.bar( x + width / 2, foldx_results, width, label = 'FoldX' )
    ax.set_ylabel( "# of results" )
    ax.set_title(f"FoldX and Graph Result Comparison \n {title}" )
    ax.set_xticks( x )
    ax.set_xticklabels( labels , rotation=0 , fontsize=14 )
    ax.legend()
    fig.tight_layout()
    plt.show()

def bar_graph_3( graph_results, foldx_results, i_mut , title ):
    labels = [ "True\nPositives", "False\nNegatives", "False\nPositives", "True\nNegatives" ]
    x = np.arange( len( labels ) )
    width = 0.3
    fig, ax = plt.subplots(figsize=(9,5))
    ax.bar( x - width , graph_results, width, label = 'Graph' )
    ax.bar( x, foldx_results, width, label = 'FoldX' )
    ax.bar( x + width , i_mut, width, label = 'i-mut' )
    ax.set_ylabel( "# of results" )
    ax.set_title(f"Result Comparison  {title}" )
    ax.set_xticks( x )
    ax.set_xticklabels( labels , rotation=0 , fontsize=14 )
    ax.legend(loc=(1.05,0.5))

    #fig.tight_layout()
    plt.show()


def bar_stats_graph( graph_stats, foldx_stats , title ):
    labels = [ "sensitivity", "specificity", "precision" ]
    x = np.arange( len( labels ) )
    width = 0.4
    fig, ax = plt.subplots()
    ax.bar( x - width / 2, graph_stats, width, label='Graph' )
    ax.bar( x + width / 2, foldx_stats, width, label='FoldX' )
    ax.set_title(f"Statistics Comparison {title}" )
    ax.set_xticks( x )
    ax.set_xticklabels( labels)
    ax.set_ylim( [ 0, 1 ] )
    ax.legend()
    
    fig.tight_layout()
    plt.show()

def bar_stats_graph_3( graph_stats, foldx_stats,i_mut_stats , title ):
    labels = [ "Sensitivity", "Specificity", "Precision" ]
    x = np.arange( len( labels ) )
    width = 0.3
    fig, ax = plt.subplots(figsize=(9,5))
    ax.bar( x - width , graph_stats, width, label='Graph' )
    ax.bar( x, foldx_stats, width, label='FoldX' )
    ax.bar( x + width , i_mut_stats, width, label='i-Mut' )
    ax.set_title(f"Statistics Comparison {title}" )
    ax.set_xticks( x )
    ax.set_xticklabels( labels)
    ax.set_ylim( [ 0, 1 ] )
    #ax.legend(loc=(1.05,0.5))
    
    fig.tight_layout()
    plt.show()


def statistics( true_positives, false_positives, false_negatives, true_negatives ) :
    if true_positives + false_negatives == 0:
        sensitivity = -1
    else:
        sensitivity = true_positives / (true_positives + false_negatives)  # also called recall
    if false_positives + true_negatives == 0:
        specificity = -1
    else:
        specificity = true_negatives / ( false_positives + true_negatives )
    if true_positives + false_positives == 0 :
        precision = -1
    else:
        precision = true_positives / ( true_positives + false_positives )
    return sensitivity, specificity, precision

if __name__ == "__main__":
    # rutas en cluster
    real_data_path = "./Results/RealData.txt"#"../../../FoldX/RealData.txt"
    foldX_data_path = "./Results/ddg_protein_analysis_FOLDX_RESULTS.txt"#"../../../FoldX/ddg_protein_analysis_FOLDX_RESULTS.txt"
    graph_data_path = "./Results/graph_results.txt"#"results.txt"
    i_mut_data_path = "./Results/i_Mut_final_data.txt"
    
    
    real_data_results = read_real_results( real_data_path )
    foldX_data_results = read_foldx_results( foldX_data_path )
    i_mut_data_results = read_foldx_results(i_mut_data_path)
    
    print("Only last frame")
    graph_data_results = read_graph_results( graph_data_path )
    true_positives_lf, false_positives_lf, false_negatives_lf, true_negatives_lf = create_confusion_matrix( real_data_results, graph_data_results )
    g_sensitivity_lf, g_specificity_lf, g_precision_lf = statistics( true_positives_lf, false_positives_lf, false_negatives_lf, true_negatives_lf )
    print_matrix( "Graph", true_positives_lf, false_positives_lf, false_negatives_lf, true_negatives_lf )   
    
    print("Sensitivity" , " Specificity", "Precision")
    print("   {:.4f},    {:.4f} ,   {:.4f}".format(g_sensitivity_lf, g_specificity_lf, g_precision_lf))
    
    print()
   
    
    #statistics
    
      
    
    
    print("All frames")
    graph_data_results = read_graph_results_all_frames( graph_data_path )
    true_positives, false_positives, false_negatives, true_negatives = create_confusion_matrix( real_data_results, graph_data_results )
    g_sensitivity, g_specificity, g_precision = statistics( true_positives, false_positives, false_negatives, true_negatives )
    print_matrix( "Graph", true_positives, false_positives, false_negatives, true_negatives )
    print()
    #foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives = create_confusion_matrix( real_data_results, foldX_data_results )
    #print_matrix( "FoldX", foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives )

    print("Sensitivity" , " Specificity", "Precision")
    print("   {:.4f},    {:.4f} ,   {:.4f}".format(g_sensitivity, g_specificity, g_precision))
    print()
    
    
    foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives = create_confusion_matrix( real_data_results, foldX_data_results )
    foldx_sensitivity, foldx_specificity, foldx_precision = statistics( foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives )
    print_matrix( "FoldX", foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives )   
    #bar_stats_graph( [ g_sensitivity, g_specificity, g_precision ], [ foldx_sensitivity, foldx_specificity, foldx_precision ] , "All frames")  
    #bar_graph( [ true_positives, false_positives, false_negatives, true_negatives ], [ foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives ] , "All frames")
        
    
    print("Sensitivity" , " Specificity", "Precision")
    print("   {:.4f},    {:.4f} ,   {:.4f}".format(foldx_sensitivity, foldx_specificity, foldx_precision))
    
    print()  

    #bar_stats_graph( [ g_sensitivity, g_specificity, g_precision ], [ foldx_sensitivity, foldx_specificity, foldx_precision ] , "Only last frame")
    #bar_graph( [ true_positives, false_positives, false_negatives, true_negatives ], [ foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives ] , "Only last frame" )

    i_mut_true_positives, i_mut_false_positives, i_mut_false_negatives, i_mut_true_negatives = create_confusion_matrix( real_data_results, i_mut_data_results )
    i_mut_sensitivity, i_mut_specificity, i_mut_precision = statistics( i_mut_true_positives, i_mut_false_positives, i_mut_false_negatives, i_mut_true_negatives )
    print_matrix( "i_mut", i_mut_true_positives, i_mut_false_positives, i_mut_false_negatives, i_mut_true_negatives )
    print("Sensitivity" , " Specificity", "Precision")
    print("   {:.4f},    {:.4f} ,   {:.4f}".format(i_mut_sensitivity, i_mut_specificity, i_mut_precision))
    
    
    
    bar_stats_graph_3( [ g_sensitivity, g_specificity, g_precision ], [ foldx_sensitivity, foldx_specificity, foldx_precision ], [ i_mut_sensitivity, i_mut_specificity, i_mut_precision ] , "All frames")
    bar_graph_3( [ true_positives, false_positives, false_negatives, true_negatives ], [ foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives ] , [ i_mut_true_positives, i_mut_false_positives, i_mut_false_negatives, i_mut_true_negatives ] ,  "All frames" )
   
    
    bar_stats_graph_3( [ g_sensitivity_lf, g_specificity_lf, g_precision_lf ], [ foldx_sensitivity, foldx_specificity, foldx_precision ], [ i_mut_sensitivity, i_mut_specificity, i_mut_precision ] , "Last Frame")
    bar_graph_3( [ true_positives_lf, false_positives_lf, false_negatives_lf, true_negatives_lf ], [ foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives ] , [ i_mut_true_positives, i_mut_false_positives, i_mut_false_negatives, i_mut_true_negatives ] ,  "Last Frames" )   
    
    data_perf=np.genfromtxt("./Results/Size_timeG_time_Cent")
    plt.scatter(data_perf[:,0],data_perf[:,1])
    plt.title("Time graph construction vs N atoms")
    plt.ylabel("Time (s)")
    plt.xlabel("# Atoms")
    
    

    '''
    with open("./Results/Size_proteins","r") as size_file:
        size = size_file.readlines()
    sizes = {}
    for i in size :
        line = i.split()
        if(line[0] not in sizes ):
            sizes[line[0]] = line[1]
    with open("./Results/TIME_FOLDX","r") as foldx_time:
        fold_perf = foldx_time.readlines()    
    plt.figure()
    data_per_foldx=[]
    for i in fold_perf:
        line = i.split("." or "," or " ")
        plt.scatter(sizes[line[0]],float(line[-1]),c="orange")
        plt.ylabel("Time (s)")
    plt.xlabel("# Atoms")
    
    #plt.xticks(np.arange(2400,8000,1000),np.arange(2400,8000,1000))
    #plt.scatter(data_perf[:,0],data_perf[:,2])
    #plt.title("Time local density vs N atoms")
    #plt.ylabel("Time (s)")
    #plt.xlabel("# Atoms")
    '''
    
