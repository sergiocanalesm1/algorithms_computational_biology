package uniandes.algorithms.readsanalyzer;

import java.util.*;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

/**
 * Simple script that simulates error free reads from a text in fasta format
 * @author Jorge Duitama
 *
 */
public class SimpleReadsSimulator {
	/**
	 * Main class that executes the program
	 * @param args Array of arguments:
	 * args[0]: Source sequence in fasta format. If many sequences are present, it only takes the first sequence
	 * args[1]: Length of the reads to simulate
	 * args[2]: Number of reads to simulate
	 * args[3]: Path to the output file
	 * args[4]: Error rate.
	 * @throws Exception If the fasta file can not be loaded
	 */
	public static void main(String[] args) throws Exception {
		String filename = args[0];
		int readLength = Integer.parseInt(args[1]);
		int numReads = Integer.parseInt(args[2]);
		String outFile = args[3];
		double errorRate= Integer.parseInt(args[4]);
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.setSequenceType(StringBuilder.class);
		QualifiedSequenceList sequences = handler.loadSequences(filename);
		if(sequences.size()==0) throw new Exception("No sequences found in file: "+filename);
		QualifiedSequence seq = sequences.get(0);
		String sequence = seq.getCharacters().toString();
		int seqLength = sequence.length();
		System.out.println("Length of the sequence to simulate reads: "+seqLength);
		double averageRD = ((double)numReads*readLength)/seqLength;
		System.out.println("Expected average RD: "+averageRD);
		char [] fixedQS = new char [readLength];
		Arrays.fill(fixedQS, '5');
		String fixedQSStr = new String(fixedQS);
		Random random = new Random();
		
		try (PrintStream out = new PrintStream(outFile)){
			//System.out.println(sequence);
			//System.out.println(seqLength);
			for (int i = 0; i < numReads; i++) {
				int randomCharPos=random.nextInt(seqLength-readLength); // avoid out of bounds
				String SubString = sequence.substring(randomCharPos,randomCharPos+readLength);
				String Mutated="";
				for (int j = 0; j < SubString.length(); j++) {
					double randomMutation = Math.random();
					if(randomMutation<errorRate/100.0) { // change base.
						// perform mutation 
						if(SubString.charAt(j) == 'A') {
					        List<String> list = Arrays.asList("G", "T", "C");
				            Collections.shuffle(list);  
				            Mutated=Mutated+list.get(0);
							}
						if(SubString.charAt(j)== 'C') {
					        List<String> list = Arrays.asList("G", "T", "A");
				            Collections.shuffle(list);  
				            Mutated=Mutated+list.get(0);
						}
						if(SubString.charAt(j)=='T') {
					        List<String> bases = Arrays.asList("G", "C", "A");
				            Collections.shuffle(bases);  
				            Mutated=Mutated+bases.get(0);
							}
						if(SubString.charAt(j)=='G') {
					        List<String> bases = Arrays.asList("C", "T", "A");
				            Collections.shuffle(bases);  
				            Mutated=Mutated+bases.get(0);
						}
					}
					
					else {
						Mutated=Mutated+SubString.charAt(j);
					}
					
					
				}
				//System.out.println(Mutated);
				String Num= Integer.toString(i);
				out.append("@"+Num+"\n");
				out.append(Mutated+"\n");
				out.append("+"+"\n");
				out.append(fixedQSStr+"\n");
				
			}


			// TODO: Generar lecturas aleatorias. Utilizar el objeto random para generar una posicion aleatoria de inicio
			// en la cadena sequence. Extraer la lectura de tamanho readLength e imprimirla en formato fastq.
			// Utilizar la cadena fixedQSStr para generar calidades fijas para el formato
			
			
		}
	}
}
