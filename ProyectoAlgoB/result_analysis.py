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
    ''' example result

-------------------WILDTYPES---------------------{
  'WT156': {
    'max_energy': -68445.21875,
    'min_energy': -68385.65625,
    'max_frame': 3.0,
    'WT156-18': {
      'local_centrality': 0.3939139505890136,
      'frame': 3.0,
      'energy': -68445.21875
    }
  },
  'WT184': {
    'max_energy': -68445.21875,
    'min_energy': -68385.65625,
    'max_frame': 3.0,
    'WT184-85': {
      'local_centrality': 0.4609281365357879,
      'frame': 3.0,
      'energy': -68445.21875
    },
    'WT184-86': {
      'local_centrality': 0.44337408015568935,
      'frame': 3.0,
      'energy': -68445.21875
    },
    'WT184-87': {
      'local_centrality': 0.42975928944322395,
      'frame': 3.0,
      'energy': -68445.21875
    },
    'WT184-88': {
      'local_centrality': 0.4740139700499994,
      'frame': 3.0,
      'energy': -68445.21875
    },
    'WT184-89': {
      'local_centrality': 0.4773581463236636,
      'frame': 3.0,
      'energy': -68445.21875
    }
  },
}-------------------MUTANTS---------------------{
  'GA156A': {
    'max_energy': -68220.5625,
    'min_energy': -67884.367188,
    'max_frame': 3.0,
    'GA156A-18': {
      'local_centrality': 0.41689550513079926,
      'frame': 3.0,
      'energy': -67884.367188
    }
  },
  'WA184A': {
    'max_energy': -66279.210938,
    'min_energy': -65931.523438,
    'max_frame': 3.0,
    'WA184A-85': {
      'local_centrality': 0.3780765311108319,
      'frame': 3.0,
      'energy': -65931.523438
    }
  },
  'CA218S': {
    'max_energy': -68322.921875,
    'min_energy': -68128.914062,
    'max_frame': 3.0,
    'CA218S-152': {
      'local_centrality': 0.4293942029666894,
      'frame': 3.0,
      'energy': -68128.914062
    },
    'CA218S-153': {
      'local_centrality': 0.4408480074142725,
      'frame': 3.0,
      'energy': -68128.914062
    }
  },
}
    '''
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
        max = dict_with_nodes[ protein ][ "max_energy" ]
        min = dict_with_nodes[ protein ][ "min_energy" ]
        acumulated_centrality = 0
        total_nodes = 0
        for node in dict_with_nodes[ protein ] :
            if node not in [ "max_energy", "min_energy", "max_frame" ]:
                acumulated_centrality += dict_with_nodes[ protein ][ node ][ "local_centrality" ] * local_centrality_percentage( max, min, dict_with_nodes[ protein ][ node ][ "energy" ] )
                total_nodes += 1
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
            for mutant in mutants:
                if mutant in real_data_results[ protein ]:
                    if real_data_results[ protein ][ mutant ] > 0 and data_results[ protein ][ mutant ] > 0:
                        true_positives += 1
                    elif real_data_results[ protein ][ mutant ] < 0 and data_results[ protein ][ mutant ] < 0:
                        true_negatives += 1
                    elif real_data_results[ protein ][ mutant ] > 0 and data_results[ protein ][ mutant ] < 0:
                        false_negatives += 1
                    elif real_data_results[ protein ][ mutant ] < 0 and data_results[ protein ][ mutant ] > 0:
                        false_positives += 1
                    else:
                        print( "al menos una estabilidad es 0, mutante: ", mutant)
        else:
            print( "create confusion matrix error: protein not in real data set" )
    return true_positives, false_positives, false_negatives, true_negatives

def print_matrix( name, true_positives, false_positives, false_negatives, true_negatives ):
    print( "{} Confusion Matrix \nTrue Positives: {}, False Negatives: {} \nFalse Positives: {}, True Negatives: {}".format( name, true_positives, false_positives, false_negatives, true_negatives ) )

def bar_graph( graph_results, foldx_results ):
    labels = [ "True\nPositives", "False\nNegatives", "False\nPositives", "True\nNegatives" ]
    x = np.arange( len( labels ) )
    width = 0.4
    fig, ax = plt.subplots()
    ax.bar( x - width / 2, graph_results, width, label = 'Graph' )
    ax.bar( x + width / 2, foldx_results, width, label = 'FoldX' )
    ax.set_ylabel( "# of results" )
    ax.set_title( "FoldX and Graph Result Comparison" )
    ax.set_xticks( x )
    ax.set_xticklabels( labels , rotation=0 , fontsize=14 )
    ax.legend()
    fig.tight_layout()
    plt.show()

def bar_stats_graph( graph_stats, foldx_stats ):
    labels = [ "sensitivity", "specificity", "precision" ]
    x = np.arange( len( labels ) )
    width = 0.4
    fig, ax = plt.subplots()
    ax.bar( x - width / 2, graph_stats, width, label='Graph' )
    ax.bar( x + width / 2, foldx_stats, width, label='FoldX' )
    ax.set_title( "Statistics Comparison" )
    ax.set_xticks( x )
    ax.set_xticklabels( labels)
    ax.set_ylim( [ 0, 1 ] )
    ax.legend()
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

    real_data_results = read_real_results( real_data_path )
    foldX_data_results = read_foldx_results( foldX_data_path )


    graph_data_results = read_graph_results( graph_data_path )
    true_positives, false_positives, false_negatives, true_negatives = create_confusion_matrix( real_data_results, graph_data_results )
    print_matrix( "Graph", true_positives, false_positives, false_negatives, true_negatives )
    print()
    foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives = create_confusion_matrix( real_data_results, foldX_data_results )
    print_matrix( "FoldX", foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives )

    bar_graph( [ true_positives, false_positives, false_negatives, true_negatives ], [ foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives ] )
    #statistics
    g_sensitivity, g_specificity, g_precision = statistics( true_positives, false_positives, false_negatives, true_negatives )
    foldx_sensitivity, foldx_specificity, foldx_precision = statistics( foldx_true_positives, foldx_false_positives, foldx_false_negatives, foldx_true_negatives )
    bar_stats_graph( [ g_sensitivity, g_specificity, g_precision ], [ foldx_sensitivity, foldx_specificity, foldx_precision ] )
