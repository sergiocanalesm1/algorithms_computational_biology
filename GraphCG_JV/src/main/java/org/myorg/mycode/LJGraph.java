package org.myorg.mycode;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import org.jgrapht.graph.DefaultWeightedEdge;


public class LJGraph {
	private double cutOff;


	private Map<Atom, Map<Atom, Double>> LJ_Graph = new HashMap<>();
	public LJGraph() {
		this.cutOff = cutOff;
	}
	

	

	public HashMap<Integer, String> Read_Itp(){
		HashMap<Integer, String> Atom_dict = new HashMap<>();
		try{
			// Open the file that is the first 
			// command line parameter
			FileInputStream fstream = new FileInputStream("/home/david/Documents/BionIF/Algortimos/Proyecto/GraphCG/data/Protein_A.itp");
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			Boolean start= false;
			//Read File Line By Line
			while ((strLine = br.readLine()) != null)   {
				// Print the content on the console
				if (strLine.startsWith("[ bonds ]") ) {
					break;
				}
				if (strLine.startsWith("[ atoms ]")) {
					start=true;
					continue;
				}
				if(start ) {
					String[] arrOfStr = strLine.split("\\s+");
					if(arrOfStr.length>1) {
						String At_indx = arrOfStr[1];
						String At_type = arrOfStr[2];
						Atom_dict.put(Integer.parseInt(At_indx),At_type);
					}
				}

				//[ATOM, 1, BB, THR, 1, 32.800, 45.200, 38.300, 1.00, 0.00]
			}

			//Close the input stream
			in.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

		return Atom_dict;

	}
	
	





}
