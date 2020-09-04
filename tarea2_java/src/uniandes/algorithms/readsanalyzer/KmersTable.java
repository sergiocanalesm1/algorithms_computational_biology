package uniandes.algorithms.readsanalyzer;

import java.util.Set;
import java.util.HashMap;
import java.util.Arraylist;

import ngsep.sequences.RawRead;
/**
 * Stores abundances information on a list of subsequences of a fixed length k (k-mers)
 * @author Jorge Duitama
 */
public class KmersTable implements RawReadProcessor {

	private HashMap<String,int> kmer_table;
	private int kmer_size;
	/**
	 * Creates a new table with the given k-mer size
	 * @param kmerSize length of k-mers stored in this table
	 */
	public KmersTable(int kmerSize) {
		// TODO: Implementar metodo
		this.kmer_table = new HashMap<>();
		this.kmer_size = kmer_size;
	}

	/**
	 * Identifies k-mers in the given read
	 * @param read object to extract new k-mers
	 */
	public void processRead(RawRead read) {
		String sequence = read.getSequenceString();
		// TODO Implementar metodo. Calcular todos los k-mers del tamanho dado en la constructora y actualizar la abundancia de cada k-mer

		int i = 0; //recorre el string. desde este indice se calculan los kmers
		while(i <= sequence.length()-kmer_size){//se le resta el kmer size porque no quiero tener index out of bounds exception
			String kmer = sequence.substring(i,i+kmer_size); //+1 porque el metodo substing de la clase String funciona [parametro1,parametro2)
			if (kmer_table.containsKey(kmer)){
				int kmer_count = kmer_table.get(kmer); //guarda la cantidad de kmers que hay repetidos hasta el momento
				kmer_table.put(kmer,kmer_count+1); //le suma uno a la cantidad
			}
			else{
				kmer_table.put(kmer,1); //comienza con un solo kmer
			}
			i++;
		}
	}
	
	/**
	 * List with the different k-mers found up to this point
	 * @return Set<String> set of k-mers
	 */
	public Set<String> getDistinctKmers() {
		// TODO Implementar metodo
		return kmer_table.keySet();
	}
	
	/**
	 * Calculates the current abundance of the given k-mer 
	 * @param kmer sequence of length k
	 * @return int times that the given k-mer have been extracted from given reads
	 */
	public int getAbundance(String kmer) {
		// TODO Implementar metodo
		return kmer_table.get(kmer);
	}
	
	/**
	 * Calculates the distribution of abundances
	 * @return int [] array where the indexes are abundances and the values are the number of k-mers
	 * observed as many times as the corresponding array index. Position zero should be equal to zero
	 */
    public HashMap<Integer,Integer> calculateAbundancesDistribution() {
        // TODO Implementar metodo
        HashMap<Integer,Integer> abundances = new HashMap<>();
        kmer_table.forEach((key,value) -> {
            if(abundances.containsKey(value)){
                int count = abundances.get(value);
                abundances.put(value,count+1);
            }
            else{
                abundances.put(value,1);
            }
        });
		return abundances;
    }
}
