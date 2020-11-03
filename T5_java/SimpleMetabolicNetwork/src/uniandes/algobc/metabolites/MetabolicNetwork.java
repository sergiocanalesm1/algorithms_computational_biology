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
	private Set<String> compartments = new TreeSet<String>();
	private Map<String,Reaction> reactions = new TreeMap<String,Reaction>();
	/*
	This is the graph
	 */
	private Map<String,Metabolite> metabolites = new TreeMap<String,Metabolite>();
	private HashMap<String,List<HashMap<String,Integer>>> edges = new HashMap<>();
	/**
	 * Adds a new gene product that can catalyze reactions
	 * @param enzyme New gene product
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
		ArrayList<Metabolite> onlySubstrate  = new ArrayList<>();
		List<Reaction> RxnList = getReactionsList();
		List<Metabolite> Metabolitos = getMetabolitesList();
		for (Metabolite metabol : Metabolitos) {
			boolean var=false;
			for (Reaction reaction : RxnList) {
				if(reaction.getProductsMap().containsKey( metabol.getId() )) {
					var=true;
					break;
						}
				}
			if(!var) {
				onlySubstrate.add(metabol);
			}
		}
		return onlySubstrate ;
	}

	public ArrayList<Metabolite> getOnlyProducs() {
		ArrayList<Metabolite> onlyProducts  = new ArrayList<>();
		List<Reaction> RxnList = getReactionsList();
		List<Metabolite> Metabolitos = getMetabolitesList();
		for (Metabolite metabol : Metabolitos) {
			boolean var=false;
			for (Reaction reaction : RxnList) {
				if(reaction.getReactantsMap().containsKey( metabol.getId() )) {
					var=true;
					break;
						}
				}
			if(!var) {
				onlyProducts.add(metabol);
			}
		}
		return onlyProducts ;
	}
	public Map<String, HashMap<String,Integer>> createGraph() {
		HashMap<String, HashMap<String,Integer>> edges = new HashMap<>();
		List<Reaction> Rxns = getReactionsList();
		Metabolite sustrate;
		Metabolite product;
		List<ReactionComponent> substrates;
		List<ReactionComponent> products;
		HashMap<String,Integer> tempEdge;
		for (Reaction rxn : Rxns) {
			 substrates = rxn.getReactants();
			 products = rxn.getProducts();
			for (ReactionComponent rxnComp : substrates) {
				sustrate = rxnComp.getMetabolite();
				for (ReactionComponent rxnCompProd : products) { //  buscar todos los productos asociados al metabolito i
					product = rxnCompProd.getMetabolite();
					if ( !edges.containsKey( sustrate.getId() )) {
						tempEdge = new HashMap<>();
						tempEdge.put(product.getId(), 1);
						edges.put( sustrate.getId(), tempEdge );
					}
					else if ( !edges.get(sustrate.getId()).containsKey(product.getId())){ //si el producto no existe dentro del sustrato
						tempEdge = edges.get( sustrate.getId() );
						tempEdge.put( product.getId(), 1);
						edges.put( sustrate.getId(), tempEdge );
					}
					else {
						tempEdge = edges.get( sustrate.getId() );
						tempEdge.put( product.getId(), tempEdge.get( product.getId() ) + 1);
						edges.put( sustrate.getId(), tempEdge ); // adicionar uno al valor
					}
				}
			}
		}
		return edges;
	}

	public void fileWriter(Map<String, HashMap<String,Integer>> graph) throws IOException {
		FileWriter myWriter = new FileWriter("docs/MetabolicNetwork_ecoli.txt");
		for (String sustrateId : graph.keySet()) {
			Metabolite sustrate = this.metabolites.get(sustrateId);
			for(String productId : graph.get( sustrateId ).keySet()){
				Metabolite product = this.metabolites.get(productId);
				System.out.println(sustrate.getName()+"\t"+product.getName()+"\t"+graph.get(sustrateId).get(productId)+"\n");
				myWriter.write(sustrate.getName()+"\t"+product.getName()+"\t"+graph.get(sustrateId).get(productId)+"\n");
			}

		}
		myWriter.close();
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
		/**
		System.out.println("ALL REACT VS PRODS");
		for(Reaction r:reactions) {
			//System.out.println(r.getId()+" "+r.getName()+" "+r.getReactants().size()+" "+r.getProducts().size()+" "+r.getEnzymes().size()+" "+r.getLowerBoundFlux()+" "+r.getUpperBoundFlux());
			//System.out.println(r.getReactantsMap().size() + "------------------" + r.getProductsMap().size());
			List<ReactionComponent> reactants = r.getReactants();
			List<ReactionComponent> prods = r.getProducts();
			for (ReactionComponent rC : reactants) {
				//System.out.println(rC.getMetabolite().getName()+"\t");
				for (ReactionComponent pC : prods) {
					System.out.println(rC.getMetabolite().getName()+"\t"+pC.getMetabolite().getName());
					
				}
				
			}
			
		}
		*/
		System.out.println("\nThese are the metabolites that are only substrates");
		ArrayList<Metabolite> substrates = network.getOnlySubstrates();
		for(Metabolite metabolite: substrates){
			System.out.println(metabolite.getId());
		}
		System.out.println("\nThese are the metabolites that are only products");
		ArrayList<Metabolite> products = network.getOnlyProducs();
		for(Metabolite metabolite: products){
			System.out.println(metabolite.getId());
		}
		System.out.println("\nFILE INFORMATION \n");

		Map<String, HashMap<String,Integer>> graph = network.createGraph();
		network.fileWriter( graph );
	}
}
