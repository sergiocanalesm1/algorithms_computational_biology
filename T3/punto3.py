import sys

def fastaReader(filename):
    '''
    return a list with all lines of fasta file
    '''
    with open("./{}".format(filename),"r") as test1:
        text_pure =test1.readlines() 
    seq=""
    for i in range(1,len(text_pure)):
        seq = seq + text_pure[i][:-1]
    return seq

def fastQReader(filename):
    '''
    return a list with all lines of fasta file
    '''
    with open("./{}".format(filename),"r") as test1:
        text_pure =test1.readlines()
    seq=""
    for i in range(1,len(text_pure),4):
        seq=seq+text_pure[i][:-1]
    return seq

def arreglo_sufijos( secuencia ):
    secuence_list = list(secuencia)
    sufix =  [ "".join(secuence_list[ i: ]) + "$" for i in range(len(secuencia))]
    sufix.append("$")
    sufix.sort()
    return sufix

def busqueda_lecturas( lista_lecturas, arreglo_sufijos ):
    posiciones = []
    for lectura in lista_lecturas:
        posicion = binary_search( lectura , arreglo_sufijos )
        posiciones.append( posicion )
    return posiciones

def binary_search( buscado, arreglo ):
    min = - 1
    max = len( arreglo )
    encontrado = - 1
    while min + 1 != max:
        pos = min + int( (max - min) / 2 )
        if arreglo[ pos ] == buscado + "$":
            encontrado = pos
            break
        if arreglo[ pos ] > buscado: #buscar abajo
            max = pos
        elif arreglo[ pos ] < buscado:
            min = pos
    return  encontrado

def BWT( secuencia_original, sufix_array ):
    arreglo_sufijos = []
    secuencia_original = secuencia_original
    for i in range( 0,len( secuencia ) + 1 ):
        splitter = sufix_array[ i ].find( "$" )
        arreglo_sufijos.append( sufix_array[ i ] + secuencia_original[ :-splitter ] )
    arreglo_sufijos[ 0 ] += secuencia_original
    transformada = "".join( [ arreglo_sufijos[ i ][ -1 ] for i in range(len(arreglo_sufijos)) ] )
    original = []
    for sufix in arreglo_sufijos:
        m = sufix.index("$")
        original.append(len(sufix)-m)
    return transformada, original

def matrixBWT( bwt_secuence ):
    ranking_map = {}
    for letter in bwt_secuence:
        if letter in ranking_map:
            ranking_map[letter] += 1
        else:
            ranking_map[letter] = 1
    ranking = { k: v for k, v in sorted(ranking_map.items(), key=lambda item: item) }
    matrix = [ list("-"+bwt_secuence) ] + [ [i] for i in list(ranking.keys()) ]
    #populate first row of indexes
    first = matrix[0][1]
    for i in range( 1, len( matrix ) ):
        if first == matrix[i][0]:
            matrix[i].append(1)
        else:
            matrix[i].append(0)
    #populate rest of the rows
    for i in range(2,len(matrix[0])):
        current = matrix[0][i]
        for j in range(1,len(matrix)):
            column = matrix[j]
            if current == column[0]:
                matrix[j].append(column[i-1]+1)
            else:
                matrix[j].append(column[i-1])
    return matrix

def lf_mapping( sub_secuence, matrix, original_array ):
    ranking = [ a[-1] for a in matrix[1:] ]
    cumulative_search = ""
    #ubicar filas
    suma = 0
    index_to_search = []
    for i in range(len(ranking)):
        if matrix[i + 1][0] == sub_secuence[-1]:
            index_to_search  = [ suma + i + 1 for i in range(ranking[i]) ]
            cumulative_search = sub_secuence[-1]
            break
        suma += ranking[i]
    #ubicar columnas
    for char in range(len(sub_secuence)-1,-1,-1):
        index_to_search_temp = []
        if cumulative_search == sub_secuence:
            return original[ index_to_search[0] - 1 ]
        char_before = sub_secuence[char - 1]
        for i in index_to_search:
            bwt_char = matrix[0][i]
            if bwt_char == char_before:
                #cuales bwt_char estan dentras de char_before?
                for j in range(1,len(matrix)):
                    if matrix[j][0] == bwt_char:
                        index_to_search_temp.append( matrix[j][i] )
        #buscar esos bwt_char en el ranking global
        suma = 0
        if len(index_to_search_temp) > 0:
            cumulative_search = sub_secuence[ char - 1 ] + cumulative_search
            for j in range( len(ranking) ):
                if matrix[ j + 1 ][0] < char_before:
                    suma += ranking[j]
            index_to_search = [ pos + suma for pos in index_to_search_temp ]
        else:
            return -1
    return -1

option = sys.argv[1]
file_fasta = sys.argv[2]
secuencia = fastaReader( file_fasta )
sufix = arreglo_sufijos( secuencia )
if option == "SUFFIX":
    file_fastq = sys.argv[3]
    lecturas = fastQReader(file_fastq)
    posiciones = busqueda_lecturas( lecturas, sufix )
    print("Each read was found on the following positions (respectively) \n {}".format(posiciones))
elif option == "BWT":
    lookfor = sys.argv[3]
    transformada, original = BWT( secuencia, sufix )
    mat = matrixBWT( transformada )
    index = lf_mapping(lookfor, mat, original)
    if index != -1 :
        print(" {} was found in the suffix array at the position: {}".format(lookfor, index))
    else:
        print(" {}  was not found".format(lookfor))



