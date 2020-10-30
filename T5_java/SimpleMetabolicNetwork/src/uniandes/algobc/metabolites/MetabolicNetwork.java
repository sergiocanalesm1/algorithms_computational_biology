package uniandes.algobc.metabolites;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.io.FileWriter;   // Import the FileWriter class
//import com.sun.org.apache.bcel.internal.classfile.Node;

//import com.sun.org.apache.bcel.internal.classfile.Node;

//import jdk.jfr.internal.PrivateAccess;
/**
 * Represents a metabolic network of reactions on metabolites
 * @author Jorge Duitama
 */
public class MetabolicNetwork {

	private Map<String,Enzyme> enzymes = new TreeMap<String,Enzyme>(); 
	private Map<String,Metabolite> metabolites = new TreeMap<String,Metabolite>();
	private Set<String> compartments = new TreeSet<String>();
	private Map<String,Reaction> reactions = new TreeMap<String,Reaction>();
	/**
	 * Adds a new gene product that can catalyze reactions
	 * @param product New gene product
	 */
	public void addEnzyme(Enzyme enzyme) {
		enzymes.put(enzyme.getId(), enzyme);
	}
	/**
	 * Adds a new metabolite. If a metabolite with the given name is already added, it 
	 * @param metabolite New metabolite
	 */
	public void addMetabolite(Metabolite metabolite) {
		metabolites.put(metabolite.getId(), metabolite);
		compartments.add(metabolite.getCompartment());
	}
	/**
	 * Adds a new reaction
	 * @param r New reaction between metabolites
	 */
	public void addReaction(Reaction r) {
		reactions.put(r.getId(),r);
		
	}
	/**
	 * Returns the gene product with the given id
	 * @param id of the product to search
	 * @return GeneProduct with the given id
	 */
	public Enzyme getEnzyme (String id) {
		return enzymes.get(id);
	}
	/**
	 * Returns the metabolite with the given id
	 * @param id of the metabolite to search
	 * @return Metabolite with the given id
	 */
	public Metabolite getMetabolite (String id) {
		return metabolites.get(id);
	}
	/**
	 * @return List of metabolites in the network
	 */
	public List<Metabolite> getMetabolitesList() {
		return new ArrayList<Metabolite>(metabolites.values());
	}
	/**
	 * @return List of reactions in the network
	 */
	public List<Reaction> getReactionsList () {
		return new ArrayList<Reaction>(reactions.values());
	}
	
	public ArrayList<Metabolite> getOnlySubstrates() {
		
		ArrayList<Metabolite> Substrate  = new ArrayList<>();
		List<Reaction> RxnList = getReactionsList();
		List<Metabolite> Metabolitos = getMetabolitesList();
		for (Metabolite metabol : Metabolitos) {
			boolean var=false;
			for (Reaction reaction : RxnList) {
				if(reaction.getProducts().contains(metabol)) {
					var=true;
					break;
						}
				}
			if(var==false) {
				Substrate.add(metabol);
			}
		}
		return Substrate ;
	}
	
	public ArrayList<Metabolite> getOnlyProducs() {
		
		ArrayList<Metabolite> Substrate  = new ArrayList<>();
		List<Reaction> RxnList = getReactionsList();
		List<Metabolite> Metabolitos = getMetabolitesList();
		for (Metabolite metabol : Metabolitos) {
			boolean var=false;
			for (Reaction reaction : RxnList) {
				if(reaction.getReactants().contains(metabol)) {
					var=true;
					break;
						}
				}
			if(var==false) {
				Substrate.add(metabol);
			}
			
		}
	
		return Substrate ;
	}
	public Map<Metabolite[], Integer> createGraph() {
		//Map<Metabolite,Integer> Nodes = new HashMap<>();
		Map<Metabolite[], Integer> rxn_counts = new HashMap<>();
		
		List<Reaction> Rxns = getReactionsList();
		for (Reaction rxn : Rxns) {
			List<ReactionComponent> substrates = rxn.getReactants();
			List<ReactionComponent> producs = rxn.getProducts();
			Metabolite[] sub_prod = new Metabolite[2]; 
			for (ReactionComponent rxnComp : substrates) {
					sub_prod[0]=(rxnComp.getMetabolite());
					for (ReactionComponent rxnCompProd : producs) { //  buscar todos los productos asociados al metabolito i 
			
						sub_prod[1]=(rxnCompProd.getMetabolite());	
		
						if (!rxn_counts.containsKey(sub_prod)) { // si rxn_counts no tienen el arreglo [reac,prod]
							rxn_counts.put(sub_prod, 1); // agregelo con 1
							//System.out.println(sub_prod[0].getName()+"+"+sub_prod[1].getName());

						}else { // si rnx_counts contiene al arreglo 
							int numrxn = rxn_counts.get(sub_prod); 
							rxn_counts.put(sub_prod, numrxn+1); // adicionar uno al valor 
							
						}
						
		
					}
	
				
			}
		}
		
		return rxn_counts;	
		
		
	}
	
	public static void main(String[] args) throws IOException {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network = loader.loadNetwork("data/e_coli_core.xml");
		System.out.println("Enzymes");
		for(Enzyme enzyme:network.enzymes.values()) {
			System.out.println(enzyme.getId()+" "+enzyme.getName());
		}
		System.out.println();
		
		List<Metabolite> metabolitesList = network.getMetabolitesList();
		System.out.println("Loaded "+metabolitesList.size()+" metabolites: ");
		for(Metabolite m:metabolitesList) {
			System.out.println(m.getId()+" "+m.getName()+" "+m.getCompartment()+" "+m.getChemicalFormula());
		}
		System.out.println();
		List<Reaction> reactions = network.getReactionsList();
		System.out.println("Loaded "+reactions.size()+" reactions");
		for(Reaction r:reactions) {
			System.out.println(r.getId()+" "+r.getName()+" "+r.getReactants().size()+" "+r.getProducts().size()+" "+r.getEnzymes().size()+" "+r.getLowerBoundFlux()+" "+r.getUpperBoundFlux());
		}
		System.out.println("file information \n");
		
		Map<Metabolite[], Integer> graph = network.createGraph();
		FileWriter myWriter = new FileWriter("MetabolicNetwork_ecoli.txt");
		for (Metabolite[]  pair : graph.keySet()) {
			Metabolite reac = pair[0];
			Metabolite produc = pair[1];
			System.out.println(reac.getName()+"\t"+produc.getName()+"\t"+graph.get(pair));
			myWriter.write(reac.getName()+"\t"+produc.getName()+"\t"+graph.get(pair)+"\n");
		}
		myWriter.close();
	}
}
