import java.util.HashMap;
import java.util.Set;
public class scratch{
    static HashMap<String,Integer> kmer_table = new HashMap<>();

    public static HashMap<Integer,Integer> calculateAbundancesDistribution() {
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
    public static void processRead() {
        String sequence = "ATAGTCTAGG";
        int kmer_size = 3;
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

    public static void main(String[] args) {
        
        // TODO Implementar metodo
        processRead();
        kmer_table.forEach((key,value) -> System.out.println(key + "," + value));
        HashMap<Integer,Integer> distri = calculateAbundancesDistribution();
        distri.forEach((key,value) -> System.out.println(key + "," + value));
    
    }
}