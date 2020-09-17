package uniandes.algorithms.readsanalyzer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import javax.swing.text.MaskFormatter;
import javax.swing.text.StyledEditorKit.BoldAction;
import javax.swing.text.html.MinimalHTMLWriter;

import htsjdk.samtools.util.RuntimeEOFException;
import ngsep.math.Distribution;
import ngsep.sequences.RawRead;
import java.util.stream.*;
/**
 * Represents an overlap graph for a set of reads taken from a sequence to assemble
 * @author Jorge Duitama
 *
 */
public class OverlapGraph implements RawReadProcessor {

	private int minOverlap;
	private Map<String,Integer> readCounts = new HashMap<>();
	private Map<String,ArrayList<ReadOverlap>> overlaps = new HashMap<>();
	private Map<String,Boolean> fathers = new HashMap<>();
	/**
	 * Creates a new overlap graph with the given minimum overlap
	 * @param minOverlap Minimum overlap
	 */
	public OverlapGraph(int minOverlap) {
		this.minOverlap = minOverlap;
	}

	/**
	 * Adds a new read to the overlap graph
	 * @param read object with the new read
	 */
	public void processRead(RawRead read) {
		String sequence = read.getSequenceString();
		//TODO: Paso 1. Agregar la secuencia al mapa de conteos si no existe.
		//Si ya existe, solo se le suma 1 a su conteo correspondiente y no se deben ejecutar los pasos 2 y 3 
		if(readCounts.containsKey(sequence)) { // la seq existe 
			int count_seq = readCounts.get(sequence); 
			readCounts.put(sequence,count_seq+1);
			return;
		} else { // seq no existe en readCounts
			readCounts.put(sequence,1);
			fathers.put(sequence,true);
		}
		
		//TODO: Paso 2. Actualizar el mapa de sobrelapes con los sobrelapes en los que la secuencia nueva sea predecesora de una secuencia existente
		//2.1 Crear un ArrayList para guardar las secuencias que tengan como prefijo un sufijo de la nueva secuencia
		//2.2 Recorrer las secuencias existentes para llenar este ArrayList creando los nuevos sobrelapes que se encuentren.
		//2.3 Después del recorrido para llenar la lista, agregar la nueva secuencia con su lista de sucesores al mapa de sobrelapes 		
		ArrayList<ReadOverlap> PrefixList = new ArrayList<>();
		for (String possiblePrefix : readCounts.keySet()) {
			int overlapLenght =getOverlapLength(sequence, possiblePrefix); 
			if(overlapLenght>minOverlap) { // si overlaplength es mayor al minimo overlap entonces agrege
				ReadOverlap OverLap = new ReadOverlap(sequence,possiblePrefix,overlapLenght);
				PrefixList.add(OverLap);	
				fathers.put(possiblePrefix,false);
							}
				}
		//if(PrefixList.size()>0) {
		overlaps.put(sequence,PrefixList);
		//}value
		//TODO: Paso 3. Actualizar el mapa de sobrelapes con los sobrelapes en los que la secuencia nueva sea sucesora de una secuencia existente
		// Recorrer el mapa de sobrelapes. Para cada secuencia existente que tenga como sufijo un prefijo de la nueva secuencia
		//se agrega un nuevo sobrelape a la lista de sobrelapes de la secuencia existente
		int MaxSucc=0;
		for (String seq_inOverlap : overlaps.keySet()) {
			int overlapLenght =getOverlapLength(seq_inOverlap,sequence);
			if(overlapLenght>minOverlap) { // si overlaplength es mayor al minimo overlap entonces agrege
				fathers.put(sequence,false);
				ReadOverlap Overlap_suf= new ReadOverlap(seq_inOverlap,sequence,overlapLenght);
				ArrayList<ReadOverlap> arrayofSeq = overlaps.get(seq_inOverlap);
				int tempMaxSuc=arrayofSeq.size();
				if(tempMaxSuc>MaxSucc) { // count the max number of succesosr for overlaps distribution.
					MaxSucc=tempMaxSuc;
				}
				arrayofSeq.add(Overlap_suf);
				overlaps.put(seq_inOverlap, arrayofSeq);
			}
		}
	}
	
	
	
	/**
	 * Returns the length of the maximum overlap between a suffix of sequence 1 and a prefix of sequence 2
	 * @param sequence1 Sequence to evaluate suffixes
	 * @param sequence2 Sequence to evaluate prefixes
	 * @return int Maximum overlap between a prefix of sequence2 and a suffix of sequence 1
	 */
	private int getOverlapLength(String sequence1, String sequence2) {
		// TODO Implementar metodo
		int i = 0;
		int j =0;
		int l1 = sequence1.length();
		if (sequence1.equals(sequence2)) { // sequences iguales entonces retornar 0
			return 0;
		}
		while (i < l1) {
			if (sequence1.charAt(i) == sequence2.charAt(j)) {
				j++;
				i++;
			}else {
			i = i-j+1;
			j=0;	
			}
		}
		return j;
	}
				

	/**
	 * Returns a set of the sequences that have been added to this graph 
	 * @return Set<String> of the different sequences
	 */
	public Set<String> getDistinctSequences() {
		//TODO: Implementar metodo

		return readCounts.keySet();
	}

	/**
	 * Calculates the abundance of the given sequence
	 * @param sequence to search
	 * @return int Times that the given sequence has been added to this graph
	 */
	public int getSequenceAbundance(String sequence) {
		//TODO: Implementar metodo
		return readCounts.get(sequence);
	}
	
	/**
	 * Calculates the distribution of abundances
	 * @return int [] array where the indexes are abundances and the values are the number of sequences
	 * observed as many times as the corresponding array index. Position zero should be equal to zero
	 */
	public int[] calculateAbundancesDistribution() {
		int maxValueReadCounts = Collections.max(readCounts.values());
		int distributions[] =new int[maxValueReadCounts+1];
		for (Integer abud : readCounts.values() ) {
			distributions[abud] = distributions[abud] +1 ;
		}
		
		
		double sum = IntStream.of(distributions).sum();
		System.out.println("mean abundace = "+sum/distributions.length);
		
		
		return distributions;
	}

	/**
	 * Calculates the distribution of number of successors
	 * @return int [] array where the indexes are number of successors and the values are the number of 
	 * sequences having as many successors as the corresponding array index.
	 */
	public int[] calculateOverlapDistribution() {

		
		int maxSuccNumber = 0;
		for (ArrayList<ReadOverlap> readarray : overlaps.values()) {
			int tempmaxSucc = readarray.size();
			if(tempmaxSucc> maxSuccNumber) {
				maxSuccNumber=tempmaxSucc;
			}
		}

		int distributionsOverlap[] =new int[maxSuccNumber+1];
		

		
		for (ArrayList<ReadOverlap> readarray : overlaps.values()) {
			int Nsuc= readarray.size();
			distributionsOverlap[Nsuc] = distributionsOverlap[Nsuc]+1; 
		}
		
		double sum = IntStream.of(distributionsOverlap).sum();
		System.out.println("Max Succesor number "+maxSuccNumber);
		System.out.println("mean succesors = "+sum/distributionsOverlap.length);
		return distributionsOverlap;
	}
	/**
	 * Predicts the leftmost sequence of the final assembly for this overlap graph
	 * @return String Source sequence for the layout path that will be the left most subsequence in the assembly
	 */
	public String getSourceSequence () {
		// TODO Implementar metodo recorriendo las secuencias existentes y buscando una secuencia que no tenga predecesores
		//String answer = fathers.get(true);
		String answer = null;
		
		for (String seq : fathers.keySet()) {
			Boolean ans = fathers.get(seq);
			if(ans) {
				return seq;
			}
			
		}

		return answer;
	}
	
	/**
	 * Calculates a layout path for this overlap graph
	 * @return ArrayList<ReadOverlap> List of adjacent overlaps. The destination sequence of the overlap in 
	 * position i must be the source sequence of the overlap in position i+1. 
	 */
	public String getAssembly() {
		
		//ArrayList<ReadOverlap> layout = new ArrayList<>();
		HashSet<String> visitedSequences = new HashSet<>();
		ReadOverlap Next = null;
		//ArrayList<String> visitedSequences = new ArrayList<>();
		String Initial = getSourceSequence();
		
		String Complete =Initial;
		//visitedSequences.add(Initial);
		
		ArrayList<ReadOverlap> OverlapsInitial = overlaps.get(Initial);

		while(true) {
		
		int maxOverlap= 0;
		for (int i = 0; i < OverlapsInitial.size(); i++) {
			Boolean Control =true;
			ReadOverlap OverlapInterno = OverlapsInitial.get(i);
			int tempMaxOverlpa = OverlapInterno.getOverlap();
			if(tempMaxOverlpa>maxOverlap ) { // revisar visited seqcuences.
				maxOverlap=tempMaxOverlpa;
				Next=OverlapInterno;
								}
						}
			if(Next==null) {
				break;
			}
				
			String NewSeq = Next.getDestSequence();
			int overlNexSeq = Next.getOverlap();
			String NewSeqNoOverlap= NewSeq.substring(overlNexSeq);
			if( visitedSequences.contains(NewSeq)) {
				//System.out.println(Complete);
				break;	
			}
			Complete = Complete+NewSeqNoOverlap;
			visitedSequences.add(NewSeq);
			
			OverlapsInitial = overlaps.get(NewSeq);
			
		}
		
		// TODO Implementar metodo. Comenzar por la secuencia fuente que calcula el método anterior
		// Luego, hacer un ciclo en el que en cada paso se busca la secuencia no visitada que tenga mayor sobrelape con la secuencia actual.
		// Agregar el sobrelape a la lista de respuesta y la secuencia destino al conjunto de secuencias visitadas. Parar cuando no se encuentre una secuencia nueva
		
		return Complete;
	}


	
}
