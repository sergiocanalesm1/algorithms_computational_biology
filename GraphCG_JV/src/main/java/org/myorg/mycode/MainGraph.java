package org.myorg.mycode;

import java.util.List;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;


public class MainGraph {

	public static Double[][] Sigma_Mat;
	public static Double[][] E_Mat;
	public static HashMap<String, Integer> atom_NumList;

	public static HashMap<String, Integer> createAtomNumList(){

		HashMap<String, Integer> Atom_NumList = new HashMap<>(); 		
		try{
			// Open the file that is the first 
			// command line parameter
			FileInputStream fstream = new FileInputStream("data/AtomType.itp");
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			//Read File Line By Line
			int count = 0;
			while ((strLine = br.readLine()) != null)   {
				// Print the content on the console	 
				Atom_NumList.put(strLine,count);
				count ++;
			}		  
			//Close the input stream
			in.close();
		}
		catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
		return Atom_NumList;
	}

	public static double calcDist(double x1,double y1,double z1,double x2,double y2, double z2) {
		double Sum_sqrs =Math.pow((x1-x2), 2) + Math.pow((y1-y2), 2) + Math.pow((z1-z2), 2);
		return Math.sqrt(Sum_sqrs);
	}


	public static Double calcLJint(Atom A1,Atom A2,Double rij ) {

		int idx1 =  atom_NumList.get(A1.getAtomType());
		int idx2 =  atom_NumList.get(A2.getAtomType());

		double Sig = Sigma_Mat[idx1][idx2];
		double Ep = E_Mat[idx1][idx2];
		return 4*Ep*(  Math.pow((Sig/rij), 12) -Math.pow((Sig/rij), 6));
	}

	public static void main(String[] args)
	{
		HashMap<Integer, String> Atom_dict = LJGraph.Read_Itp();
		atom_NumList = createAtomNumList();
		/**
		 * create the Sigma and Epsilon mat to map LJ parameters
		 * 
		 */
		Sigma_Mat= new Double[39][39];
		E_Mat= new Double[39][39];
		try{
			
			// Open the file that is the first 
			// command line parameter
			FileInputStream fstream = new FileInputStream("data/LJMatrix.itp");
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			//Read File Line By Line
			while ((strLine = br.readLine()) != null) {
				// Print the content on the console	 
				String[] interact = strLine.split(";|,");
				if(interact.length>4) {
					//System.out.println(interact[0]+" "+interact[1]+" "+interact[2]+" "+interact[3]);
					String at_1 = interact[0].replace("\t","").strip();
					String at_2 = interact[1].replace("\t","").strip();
					int idx1 = atom_NumList.get(at_1);
					int idx2 = atom_NumList.get(at_2);
					E_Mat[idx1][idx2]=Double.valueOf(interact[2]);
					Sigma_Mat[idx1][idx2]=0.43;
				}
				else {
					String at_1 = interact[0].replace("\t","").strip();
					String at_2 = interact[1].replace("\t","").strip();
					int idx1 = atom_NumList.get(at_1);
					int idx2 = atom_NumList.get(at_2);
					E_Mat[idx1][idx2]=Double.valueOf(interact[2]);
					Sigma_Mat[idx1][idx2]=Double.valueOf(interact[3]);
					//System.out.println(interact[0]+" "+interact[1]+" "+interact[2]+" "+interact[3]);
				}
			}		  
			//Close the input stream
			in.close();
		}
		catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

		long time = System.currentTimeMillis();
		SimpleWeightedGraph<Atom,DefaultWeightedEdge> LJGraph = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
		double cutOFF = 15.0;
		try{
			List<String> PDB_file = Files.readAllLines(Paths.get("data/Test_NetworkX.pdb"));

			for (int i = 0; i < PDB_file.size(); i++) {
				String strLine_1 = PDB_file.get(i);
				if(!strLine_1.startsWith("ATOM")) {
					continue;
				}
				for (int j = i+1; j < PDB_file.size(); j++) {
					/**
					 * Creation of first atom
					 */
					String[] arrOfStr1 = strLine_1.split("\\s+");
					String strLine_2 = PDB_file.get(j);
					if (!strLine_2.startsWith("ATOM")) {
						continue;
					}
					String[] arrOfStr2 = strLine_2.split("\\s+");
					double x_c1,y_c1,z_c1;
					x_c1 = Double.parseDouble(arrOfStr1[5]);
					y_c1 = Double.parseDouble(arrOfStr1[6]);
					z_c1 = Double.parseDouble(arrOfStr1[7]);

					//System.out.println(arrOfStr[0]+" "+arrOfStr[3]+" "+arrOfStr[4]+" "+arrOfStr[5]+" ");
					int atomId1 = Integer.parseInt(arrOfStr1[1]);
					int resNum1 = Integer.parseInt(arrOfStr1[4]);
					String atomName1 = arrOfStr1[2];
					String resName1 = arrOfStr1[3];
					String atomType1;
					if(atomName1.equals("W")) {
						atomType1 = "P4";
					}
					else {
						atomType1 = Atom_dict.get(atomId1);
					}

					Atom A1 = new Atom(atomName1, atomType1, resName1, resNum1, x_c1, y_c1, z_c1, atomId1);
					LJGraph.addVertex(A1);

					double x_c2,y_c2,z_c2;

					x_c2 = Double.parseDouble(arrOfStr2[5]);
					y_c2 = Double.parseDouble(arrOfStr2[6]);
					z_c2 = Double.parseDouble(arrOfStr2[7]);

					double Dist = calcDist(x_c1, y_c1, z_c1, x_c2, y_c2, z_c2);

					if(cutOFF > Dist && Dist > 4.0) {
						int atomId2 = Integer.parseInt(arrOfStr2[1]);
						int resNum2 = Integer.parseInt(arrOfStr2[4]);
						String atomName2 = arrOfStr2[2];
						String resName2 = arrOfStr2[3];
						String atomType2;
						if(atomName2.equals("W")) {
							atomType2 = "P4";
						}
						else {
							atomType2 = Atom_dict.get(atomId2);
						}

						Atom A2 = new Atom(atomName2, atomType2, resName2, resNum2, x_c2, y_c2, z_c2, atomId1);
						LJGraph.addVertex(A2);

						double LJ_int = calcLJint(A1, A2, Dist);
						DefaultWeightedEdge e = LJGraph.addEdge(A1, A2);
						LJGraph.setEdgeWeight(e, LJ_int);
					}
				}
			}
			//			String line = Files.readAllLines(Paths.get("/home/david/Documents/BionIF/Algortimos/Proyecto/GraphCG/data/Test_NetworkX.pdb")).get(1);
			//System.out.println(line);
		} 
		catch(IOException e){
			System.out.println(e.getMessage());
		}
		time = System.currentTimeMillis()-time;
		System.out.println("Time building overlap graph(s): "+time/1000);
		
		/**

		try{
			// Open the file that is the first 
			// command line parameter
			FileInputStream fstream = new FileInputStream("/home/david/Documents/BionIF/Algortimos/Proyecto/GraphCG/data/Test_NetworkX.pdb");
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			//Read File Line By Line
			while ((strLine = br.readLine()) != null)   {
				// Print the content on the console
				if (strLine.startsWith("ATOM")) {

					String[] arrOfStr = strLine.split("\\s+");
					int at_idx,res_num ;
					String at_type,res_name,at_name;
					Double x_c,y_c,z_c;
					//System.out.println(arrOfStr[0]+" "+arrOfStr[3]+" "+arrOfStr[4]+" "+arrOfStr[5]+" ");
					at_idx = Integer.parseInt(arrOfStr[1]);
					res_num = Integer.parseInt(arrOfStr[4]);
					at_name = arrOfStr[2];
					res_name = arrOfStr[3];
					x_c = Double.parseDouble(arrOfStr[5]);
					y_c = Double.parseDouble(arrOfStr[6]);
					z_c = Double.parseDouble(arrOfStr[7]);
					at_type = Atom_dict.get(at_idx);
					//Atom A1 = new Atom(at_type, res_name, res_num, x_c, y_c, z_c, at_idx);
					//[ATOM, 1, BB, THR, 1, 32.800, 45.200, 38.300, 1.00, 0.00]
				}

			}
			//Close the input stream
			in.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
		 */


	}

}


