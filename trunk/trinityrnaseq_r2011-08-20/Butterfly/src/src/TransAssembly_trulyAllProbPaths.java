
//import java.util.Set;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistanceWoVer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import gnu.getopt.Getopt;


import java.io.BufferedReader;
//import java.io.FileNotFoundException;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import gnu.getopt.*;

public class TransAssembly_trulyAllProbPaths {

	private static final boolean DEBUG = true;
	private static final SeqVertex ROOT = new SeqVertex(-1, "S",Integer.MAX_VALUE);
	private static final SeqVertex T_VERTEX = new SeqVertex(-2, "E");
	private static final int LINE_LEN = 60;

	private static int LAST_ID = -1;
	private static int LAST_REAL_ID = -1;
	private static int MAX_DEPTH = 0;

	private static final double EDGE_THR = 0.05;
	private static final int COMP_AVG_COV_THR = 1;
	private static final int INITIAL_EDGE_ABS_THR = 0;

	private static final double FLOW_THR = 0.02;
	private static final int MIN_TRIPLET_SUPPORT_THR = 2;
	private static int MIN_OUTPUT_SEQ;

	private static int VERBOSE_LEVEL = 10;
	private static int MAX_PAIR_DISTANCE = 0;
	//	private static int MAX_PAIR_FOR_LONG_MIDDLE;
	//	private static final int MIN_SINGLE_NODE = 70;
	private static final int MAX_PATH_SIZE = 150000;
	private static final int MAX_MM_ALLOWED = 6;
	private static final int EXTREME_EDGE_FLOW_FACTOR = 200;
	private static final boolean USE_TRIPLETS = true;
//	private static final int LONGEST_MM_STRETCH = 2;
	private static int K = 0;
	private static boolean COLLAPSE_SNPs = false;

	private static boolean USE_DEGENERATE_CODE = false;
	private static String[] LETTERS = new String[]{"A","C","G","T"};
//	private static String[] LETTERS_withDeg = new String[]{"A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N"};

	private static PrintStream ERR_STREAM;

	private static Map<String, String> DEGENERATE_CODE = new HashMap<String, String>() {
		private static final long serialVersionUID = 1L;

		{ 
			put("AG","R");
			put("CT","Y");
			put("CG","S");
			put("AT","W");
			put("GT","K");
			put("AC","M");
			put("CGT","B");
			put("AGT","D");
			put("ACT","H");
			put("ACG","V");
			put("ACGT","N");
		}
	};

	private static Map<String, String> DEGENERATE_CODE_REV = new HashMap<String, String>() {
		private static final long serialVersionUID = 1L;

		{ 
			put("R","AG");
			put("Y","CT");
			put("S","GC");
			put("W","AT");
			put("K","GT");
			put("M","AC");
			put("B","CGT");
			put("D","AGT");
			put("H","ACT");
			put("V","ACG");
			put("N","ACGT");
		}
	};



	public static void main(String[] args) throws Exception 
	{
		int totalNumReads = 0;

		String file = "";
		boolean printUsage = false;
		LongOpt[] longopts = new LongOpt[3];
		longopts[0] = new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h');

		Getopt g = new Getopt("TransAssembly", args, "L:F:N:C:V:Sh",longopts);
		int c;

		while ((c = g.getopt()) != -1)
		{
			switch(c)
			{
			case 'h':
				printUsage = true;
				break;
			case 'L':
				MIN_OUTPUT_SEQ = Integer.parseInt(g.getOptarg());
				break;

			case 'F':
				MAX_PAIR_DISTANCE = Integer.parseInt(g.getOptarg());
				break;

			case 'N':
				totalNumReads = Integer.parseInt(g.getOptarg());
				break;
			case 'V':
				VERBOSE_LEVEL = Integer.parseInt(g.getOptarg());
				break;
			case 'C':
				file = g.getOptarg();
				break;
			case '?':
				printUsage = true;
				break; 
				//
			default:
				printUsage = true;
			}
		}

		printUsage = printUsage || file.equals("") || totalNumReads==0 || MAX_PAIR_DISTANCE == 0;

		if (printUsage)
		{
			System.err.println("");
			System.err.println("########################################################################################");
			System.err.println("#");
			System.err.println("# Required:");
			System.err.println("#  -N  <int>     total number of reads or fragment pairs");
			System.err.println("#  -L  <int>     min length for an assembled sequence to be reported");
			System.err.println("#  -F  <int>     average fragment length");
			System.err.println("#  -C  <string>  prefix for component/reads file");
			System.err.println("#  ");
			System.err.println("# Optional:");
			System.err.println("#  -V <int>                      verbosity level ");
			System.err.println("#                                   (default: 10 - progress of method + some stats)");
			System.err.println("#                                   (15 - like (10) + final paths to be added)");
			System.err.println("#                                   (20 - maximum verbosity)");
			System.err.println("#");
			System.err.println("########################################################################################");
			System.err.println("");
			System.exit(1);

		}
		String outfileName = file+"_trulyAll";
		ERR_STREAM = new PrintStream(new FileOutputStream(outfileName + ".err"));

		debugMes("Started",10);

		//		MAX_PAIR_FOR_LONG_MIDDLE = MAX_PAIR_DISTANCE-100;
		Vector<Integer> rootIDs = new Vector<Integer>();

		HashMap<Integer,Integer> outFlow = new HashMap<Integer, Integer>();
		HashMap<Integer,Integer> inFlow = new HashMap<Integer, Integer>();
		HashMap<Integer,String> firstLetter = new HashMap<Integer, String>();

		preProcessGraphFile(file+".out",outFlow,inFlow,firstLetter);

		DirectedSparseGraph<SeqVertex, SimpleEdge> graph = buildNewGraphFirstLetter(file+".out",rootIDs,outFlow,inFlow,firstLetter); 


		//		DirectedSparseGraph<SeqVertex, SimpleEdge> graph = buildNewGraphFirstLetter(file+".out",rootIDs); 
		LAST_REAL_ID = LAST_ID;
		debugMes("Graph is built",10);

		String[] tmpFile = file.split("/");
		String graphName = tmpFile[tmpFile.length-1];

		PrintStream p;
		boolean createMiddleDotFiles = false;

		if (createMiddleDotFiles) 
		{
			p= new PrintStream(new FileOutputStream(outfileName + ".dot"));
			writeDotFile(graph,p,graphName);
			p.close();
		}

		fixExtremelyHighSingleEdges(graph,outFlow,inFlow);

		removeLightEdges(graph);

		if (createMiddleDotFiles)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_removedLightEdges.dot"));
			writeDotFile(graph,p,graphName);
			p.close();
		}

		compactLinearPaths(graph);
		if (createMiddleDotFiles)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_compactonly.dot"));
			writeDotFile(graph,p,graphName);
			p.close();
		}
		
		boolean generateFullSeqGraph = false;
		if (generateFullSeqGraph)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_compactonly_fullSeq.dot"));
			writeDotFile(graph,p,graphName,generateFullSeqGraph);
			p.close();
		}


		compactPrefixesBottomUp(graph);

		if (createMiddleDotFiles)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_compact1.dot"));
			writeDotFile(graph,p,graphName);
			p.close();
		}

		int i=1;
		boolean runAgain = removeLightEdges(graph);
		boolean changend_suffices = false,changend_rmEdges=false;
		while (runAgain)
		{
			changend_suffices = compactPrefixesBottomUp(graph);
			i++;
			changend_rmEdges = removeLightEdges(graph);

			runAgain = changend_suffices || changend_rmEdges;

			if (createMiddleDotFiles)
			{
				p = new PrintStream(new FileOutputStream(outfileName + "_compact"+i+".dot"));
				writeDotFile(graph,p,graphName);
				p.close();
			}

		}

		if (createMiddleDotFiles)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_compcatDone.dot"));
			writeDotFile(graph,p,graphName);
			p.close();
		}

		if (COLLAPSE_SNPs)
			if (USE_DEGENERATE_CODE)
				removeSingleNtBubblesWithDegenerateCode(graph);
			else
				removeSingleNtBubbles(graph);

		compactLinearPaths(graph);

		if (createMiddleDotFiles)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_compcatDone_withoutBubbles.dot"));
			writeDotFile(graph,p,graphName);
			p.close();
		}

		//remove small components
		calcSubComponentsStats(graph);

		if (createMiddleDotFiles)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_goodComp.dot"));
			writeDotFile(graph,p,graphName);
			p.close();
		}

		//		PrintBestSeqs(graph,file);

		HashMap<Integer, LocInGraph> originalVerIDsMapping = getOriginalVerIDsMappingHash(graph);

		p = new PrintStream(new FileOutputStream(outfileName + "_finalComps.dot"));
		writeDotFile(graph,p,graphName);
		p.close();

		if (generateFullSeqGraph)
		{
			p = new PrintStream(new FileOutputStream(outfileName + "_finalComps.fullSeq.dot"));
			writeDotFile(graph,p,graphName,generateFullSeqGraph);
			p.close();
		}


		//		p = new PrintStream(new FileOutputStream(file + "_origID.txt"));
		//		p.println(originalVerIDsMapping);
		//		p.close();

		int numXstructs = countNumOfXstructures(graph);
		if (numXstructs>0)
			debugMes("number X structures = "+numXstructs,10);


		DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis = new DijkstraDistance<SeqVertex, SimpleEdge>(graph, true);

		HashMap<String, List<Read>> readNameHash = getReadStarts(graph,file+".reads",originalVerIDsMapping,rootIDs);
		HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash = getSuffStats_wPairs(graph,readNameHash,dijkstraDis);

//		debugMes(combinedReadHash+"",15);

		//start working on one sub component at a time:
		// look for loops, try to solve them
		// if loops remain, move on to the next subComp.

		Set<Set<SeqVertex>> comps = divideIntoComponents(graph);
		debugMes("total number of components = "+comps.size(),10);
		//		pDot.println("digraph "+name+"{");
		int compID = -1;

		PrintStream pout_trulyAll = new PrintStream(new FileOutputStream(file+"_trulyAllProbPaths.fasta"));

		
		String[] pathName = file.split("/");

		int totalNumPaths = 0;
		int totalNumSuccComps = 0;
		for (Set<SeqVertex> comp : comps)
		{
			compID++;
			debugMes("working on subcomponent "+compID,10);

			dealWithSimpleLoops(graph,comp,combinedReadHash,false);
			HashMap<Integer,Integer> nLoops = dealWithSimpleLoops(graph,comp,combinedReadHash,true);

			int numLoops = 0;
			for (Integer len : nLoops.keySet())
			{
				if (nLoops.get(len)>0)
				{
					debugMes("number of loops of length "+len+" after simplification: "+nLoops.get(len),10);
					numLoops = numLoops + nLoops.get(len).intValue();
				}
			}
			//			if (createMiddleDotFiles)
			//			{
			//				p = new PrintStream(new FileOutputStream(file + "_comp"+compID+"_finalCompsWOloops.dot"));
			//				writeDotFile(graph,comp,p,graphName);
			//				p.close();
			//			}
			if (numLoops>0)
			{
				debugMes("stopped because of the loops, comp"+compID,10);
				continue;
			}

			addSandT(graph,comp,combinedReadHash);

			//			if (createMiddleDotFiles)
			//			{
			//				p = new PrintStream(new FileOutputStream(file + "_comp"+compID+ "_withSandT.dot"));
			//				writeDotFile(graph,comp,p,graphName);
			//				p.close();
			//			}



			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer = new DijkstraDistanceWoVer<SeqVertex, SimpleEdge>(graph);
			

//			HashMap<List<Integer>,Pair<Integer>> FinalPaths_diff = getAllProbablePaths(graph,comp,
//					combinedReadHash,dijkstraDis,dijkstraDisWoVer,FinalPaths_all);

			Pair<HashMap<List<Integer>,Pair<Integer>>> FinalPathsPair = getAllProbablePaths(graph,comp,
					combinedReadHash,dijkstraDis,dijkstraDisWoVer);

			HashMap<List<Integer>,Pair<Integer>> FinalPaths_trulyAll = FinalPathsPair.getFirst();

			String name = pathName[pathName.length-1]+"_c"+compID;

			if (FinalPaths_trulyAll==null)
				continue;

			totalNumPaths+=FinalPaths_trulyAll.keySet().size();
			if (FinalPaths_trulyAll.keySet().size()>0)
				totalNumSuccComps++;

			for (List<Integer> path : FinalPaths_trulyAll.keySet())
			{
				debugMes("FinalPath: "+path+" with support "+FinalPaths_trulyAll.get(path),15);
			}


			printFinalPaths(FinalPaths_trulyAll,graph,compID,pout_trulyAll,name,totalNumReads);

			numXstructs = countNumOfXstructuresResolved(graph,comp,FinalPaths_trulyAll);
			if (numXstructs>0)
				debugMes("number X structures resolved = "+numXstructs,10);

			//
			//		String fastaFile = outfileName+"_nodesSeq.fasta";
			//		printLongNodes(graph,fastaFile);
			removeAllEdgesOfSandT(graph);

		}
		pout_trulyAll.close();
		
		debugMes("total number of paths reported = "+totalNumPaths+" from "+totalNumSuccComps +" components",10);

		p = new PrintStream(new FileOutputStream(outfileName + "_finalCompsWOloops.dot"));
		writeDotFile(graph,p,graphName);
		p.close();

		debugMes("Done",10);
		ERR_STREAM.close();

	}




	/**
	 * given the graph, find all single nt bubbles, and choose the majority vote. 
	 * add the weights to the majority path, and add the prevID  
	 * v -> v1 -> vend
	 * v -> v2 -> vend
	 * @param graph
	 */
	private static void removeSingleNtBubbles(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {

		SeqVertex v1=null ,v2 = null, vend = null;
		SeqVertex vToKeep=null ,vToRemove = null;
		SimpleEdge e1ToKeep = null, e1ToRemove = null;
		SimpleEdge e2ToKeep = null, e2ToRemove = null;

		Vector<SeqVertex> removeV = new Vector<SeqVertex>();
		Collection<SeqVertex> allV = new HashSet<SeqVertex>();
		allV.addAll(graph.getVertices());

		for (SeqVertex v : allV)
		{
			if (removeV.contains(v))
				continue;
			if (graph.getSuccessorCount(v)==2)
			{
				Collection<SeqVertex> children = graph.getSuccessors(v);
				Iterator<SeqVertex> iter = children.iterator();
				v1 = iter.next();
				v2 = iter.next();


				int len1 = v1.getName().length();
				int len2 = v2.getName().length();

				if (len1==1 && len2==1 && 
						graph.getSuccessorCount(v1)==1 && 
						graph.getSuccessorCount(v2)==1 &&
						getSingleSuccessor(graph,v2).equals(getSingleSuccessor(graph,v1)))
				{
					vend = getSingleSuccessor(graph,v1);
					if (graph.findEdge(v, v1).getWeight() > graph.findEdge(v, v2).getWeight())
					{ //keep v1, loose v2
						vToKeep = v1;
						vToRemove = v2;
					}else
					{ //keep v2, loose v1
						vToKeep = v2;
						vToRemove = v1;
					}
					e1ToKeep = graph.findEdge(v, vToKeep);
					e2ToKeep = graph.findEdge(vToKeep, vend);
					e1ToRemove = graph.findEdge(v, vToRemove);
					e2ToRemove = graph.findEdge(vToRemove, vend);
					debugMes("merging the node "+vToRemove.getID()+" to the node "+vToKeep.getID(),20);

					SeqVertex newV = new SeqVertex(getNextID(), vToKeep.getName());
					newV.copyTheRest(vToKeep);
					newV.addToPrevIDs(vToKeep,vToRemove,LAST_REAL_ID);

					graph.addVertex(newV);
					graph.addEdge(new SimpleEdge(e1ToKeep.getWeight() + e1ToRemove.getWeight()), v, newV);
					graph.addEdge(new SimpleEdge(e2ToKeep.getWeight() + e2ToRemove.getWeight()), newV,vend);

					removeV.add(vToRemove);
					removeV.add(vToKeep);


				}
			}
		}

		for (SeqVertex rv : removeV)
		{
			debugMes("removing the single nt variation vertex "+rv.getID(),20);
			graph.removeVertex(rv);
		}

	}




	/**
	 * given the graph, find all single nt bubbles, and choose the majority vote. 
	 * add the weights to the majority path, and add the prevID  
	 * v -> v1 -> vend
	 * v -> v2 -> vend
	 * @param graph
	 * @throws Exception 
	 */
	private static void removeSingleNtBubblesWithDegenerateCode(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) throws Exception {

		SeqVertex v1=null ,v2 = null, vend = null;
		SimpleEdge eTop1 = null, eTop2 = null;
		SimpleEdge eBottom1 = null, eBottom2 = null;

		Vector<SeqVertex> removeV = new Vector<SeqVertex>();
		Collection<SeqVertex> allV = new HashSet<SeqVertex>();
		allV.addAll(graph.getVertices());

		for (SeqVertex v : allV)
		{
			if (removeV.contains(v))
				continue;
			if (graph.getSuccessorCount(v)==2)
			{
				Collection<SeqVertex> children = graph.getSuccessors(v);
				Iterator<SeqVertex> iter = children.iterator();
				v1 = iter.next();
				v2 = iter.next();


				int len1 = v1.getName().length();
				int len2 = v2.getName().length();

				if (len1==1 && len2==1 && 
						graph.getSuccessorCount(v1)==1 && 
						graph.getSuccessorCount(v2)==1 &&
						getSingleSuccessor(graph,v2).equals(getSingleSuccessor(graph,v1)))
				{
					vend = getSingleSuccessor(graph,v1);
					String key;
					if (String.CASE_INSENSITIVE_ORDER.compare(v1.getName(),v2.getName())<0)
						key =  v1.getName()+v2.getName();
					else
						key =  v2.getName()+v1.getName();
					String name = getDegenerateRepresentation(key);
					SeqVertex newV = new SeqVertex(getNextID(), name);

					if (graph.findEdge(v, v1).getWeight() > graph.findEdge(v, v2).getWeight())
						newV.copyTheRest(v1);
					else
						newV.copyTheRest(v2);

					eTop1 = graph.findEdge(v, v1);
					eBottom1 = graph.findEdge(v1, vend);
					eTop2 = graph.findEdge(v, v2);
					eBottom2 = graph.findEdge(v2, vend);
					debugMes("merging the nodes "+v1.getID()+" and the node "+v2.getID()+" to the node "+newV,20);

					newV.addToPrevIDs(v1,v2,LAST_REAL_ID);

					graph.addVertex(newV);
					graph.addEdge(new SimpleEdge(eTop1.getWeight() + eTop2.getWeight()), v, newV);
					graph.addEdge(new SimpleEdge(eBottom1.getWeight() + eBottom2.getWeight()), newV,vend);

					removeV.add(v1);
					removeV.add(v2);


				}
			}
		}

		for (SeqVertex rv : removeV)
		{
			debugMes("removing the single nt variation vertex "+rv.getID(),20);
			graph.removeVertex(rv);
		}

	}








	/**
	 * return the single successor of this node in this graph
	 * @param graph
	 * @param v2
	 * @return
	 */
	private static SeqVertex getSingleSuccessor(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v) {
		Collection<SeqVertex> children = graph.getSuccessors(v);
		if (children.size()!=1)
			return null;

		SeqVertex vout = children.iterator().next();
		return vout;

	}




	/**
	 * find edges that are extremely high compared to both side (a single very abundant kmer, and fix their support
	 * @param graph
	 * @param inFlow 
	 * @param outFlow 
	 */
	private static void fixExtremelyHighSingleEdges(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, HashMap<Integer,Integer> outFlow, HashMap<Integer,Integer> inFlow) {

		for (SimpleEdge e : graph.getEdges())
		{
			double supp =e.getWeight(); 
			Integer sourceID = graph.getSource(e).getID();
			Integer targetID = graph.getDest(e).getID();
			Integer inFlowToSource = inFlow.get(sourceID);
			Integer outFlowOfTarget = outFlow.get(targetID);

			if (inFlowToSource!= null && outFlowOfTarget!= null && 
					supp > inFlowToSource*EXTREME_EDGE_FLOW_FACTOR && supp > outFlowOfTarget*EXTREME_EDGE_FLOW_FACTOR)
			{
				double newSupp = Math.max(inFlowToSource, outFlowOfTarget);
				debugMes("the support of edge "+sourceID+"->"+targetID+" has changed from "+supp+" to "+newSupp,20);
				e.setWeight(newSupp);
			}

		}
	}






	/**
	 * given the graph and the final paths, find x structures that belong to only two paths, which resolve this structure.
	 * @param graph
	 * @param comp 
	 * @param finalPaths
	 * @return
	 */
	private static int countNumOfXstructuresResolved(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Set<SeqVertex> comp, HashMap<List<Integer>,Pair<Integer>> finalPaths) {

		int res = 0;

		//		for (SeqVertex v : graph.getVertices())
		for (SeqVertex v : comp)

		{
			if (graph.inDegree(v)>1 && graph.outDegree(v)>1)
			{
				//this is an x-structure
				int maxPaths = Math.max(graph.inDegree(v), graph.outDegree(v));
				Integer bef,after;
				int vid = v.getID();
				HashMap<Pair<Integer>,Integer> befAndAfterNodes = new HashMap<Pair<Integer>, Integer>();
				Pair<Integer> key;
				for (List<Integer> path : finalPaths.keySet())
				{
					int index = path.indexOf(vid);

					if (index!=-1 && index!=0 && index!=path.size()-1) // vid is not the first or the last
					{

						bef = path.get(index-1);
						after = path.get(index+1);
						key = new Pair<Integer>(bef,after);
						if (!befAndAfterNodes.containsKey(key))
							befAndAfterNodes.put(key,1);
						else
							befAndAfterNodes.put(key,befAndAfterNodes.get(key)+1);


						//						numPaths++;
						//						if (bef==-1)
						//						{
						//						} else {
						//							if (!path.get(index-1).equals(bef) && !path.get(index+1).equals(after))
						//							{
						//								foundDisjointPaths = true;
						//							}
						//						}
					}
				}

				if (befAndAfterNodes.keySet().size()==maxPaths)
				{
					debugMes("vertex "+v.getID()+" is resolved in an X-structure",10);
					res++;
				}
			}
		}

		return res;
	}

	/**
	 * Print all final paths
	 * @param finalPaths
	 * @param graph
	 * @param compID 
	 * @param p 
	 * @param name 
	 * @param totalNumReads 
	 * @throws FileNotFoundException 
	 */
	private static void printFinalPaths(
			HashMap<List<Integer>,Pair<Integer>> finalPaths,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, int compID, PrintStream p, String name, int totalNumReads) throws FileNotFoundException {


		int i=1;
		double fpkm_all = 0;
		double fpkm_rel = 0;

		DecimalFormat df = new DecimalFormat("#.###");
		for (List<Integer> path : finalPaths.keySet())
		{
			//print this path
			String seqName = name+"_seq"+i;

			String seq = getPathSeq(graph,path);


			if (seq.length()>=MIN_OUTPUT_SEQ) 
			{
				//				rpkm formula = #reads*1e9/(length*totalMapped);
				//				fpkm = (finalPaths.get(path)/(double)seq.length()) *(1e9/totalNumReads);
				fpkm_all = (finalPaths.get(path).getFirst()/(double)seq.length()) *(1e9/totalNumReads);
				fpkm_rel = (finalPaths.get(path).getSecond()/(double)seq.length()) *(1e9/totalNumReads);

				seqName = seqName + "_FPKM_all:" +df.format(fpkm_all)+ "_FPKM_rel:" +df.format(fpkm_rel)+ "_len:"+seq.length()+"_path:"+ path.subList(1, path.size()-1);
				seqName = seqName.replaceAll(" ", "");
				p.print(getSeqFasta(seq, seqName));
				i++;
			}

		}
	}


	/**
	 * given a path in the graph, return its sequence
	 * @param graph
	 * @param path
	 * @return
	 */
	private static String getPathSeq(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, List<Integer> path) {
		String seq = "";
		for (Integer nodeID : path)
			if (nodeID>=0)
				seq = seq + getSeqVertex(graph, nodeID).getName();	
		return seq;
	}




	/**
	 * For each path of a read pair, ask how many reads support it.
	 * @param graph
	 * @param readNameHash
	 * @param dijkstraDis 
	 * @return
	 */
	private static HashMap<Integer, HashMap<PairPath, Integer>> getSuffStats_wPairs(
			DirectedSparseGraph<SeqVertex,SimpleEdge> graph, HashMap<String, List<Read>> readNameHash, DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis) {
		HashMap<Integer,HashMap<PairPath,Integer>>   combinedReadHash = new HashMap<Integer,HashMap<PairPath,Integer>>  ();

		Set<String> usedReads = new HashSet<String>();
		List<Read> curList = null;
		int numReadsUsed = 0;
		//		debugMes(readNameHash.keySet()+"",10);
		for (String name : readNameHash.keySet())
		{
			//			if (name.equals(">SRR039231.7789180_FC42DB6AAXX:2:43:712:1680"))
			//				debugMes("stop here",10);
			if (usedReads.contains(name))
				continue;

			curList = readNameHash.get(name);
			if (curList.size()==1)
			{//single read
				Read r = curList.get(0);
				PairPath path = new PairPath(r.getPathIDs());
				Integer firstV = path.getFirstID();

				if (!combinedReadHash.containsKey(firstV))
					combinedReadHash.put(firstV, new HashMap<PairPath,Integer>());

				if (!combinedReadHash.get(firstV).containsKey(path))
					combinedReadHash.get(firstV).put(path, 0);

				Integer counts = combinedReadHash.get(firstV).get(path);
				combinedReadHash.get(firstV).put(path,++counts);
				numReadsUsed++;
				debugMes("we have "+combinedReadHash.get(firstV).get(path)+" reads supporting the path: "+path,20);

			}else {// paired read
				Read r1 = curList.get(0);
				List<Integer> path1 = r1.getPathIDs();

				Read r2 = curList.get(1);
				List<Integer> path2 = r2.getPathIDs();


				PairPath  combinedPath = combinePaths(graph,path1,path2,dijkstraDis);
				if (combinedPath.isEmpty())
				{
					debugMes("the paths "+path1+" and "+path2+" couldn't be combined",15);
					continue;
				}

				Integer firstV = combinedPath.getFirstID();

				if (!combinedReadHash.containsKey(firstV))
					combinedReadHash.put(firstV, new HashMap<PairPath,Integer>());

				if (!combinedReadHash.get(firstV).containsKey(combinedPath))
					combinedReadHash.get(firstV).put(combinedPath, 0);

				Integer counts = combinedReadHash.get(firstV).get(combinedPath);
				combinedReadHash.get(firstV).put(combinedPath,++counts);
				debugMes("we have "+combinedReadHash.get(firstV).get(combinedPath)+" reads supporting the path: "+combinedPath,20);

				numReadsUsed++;


			}
			usedReads.add(name);
		}
		debugMes("number of reads used = "+numReadsUsed,15);
		return combinedReadHash;
	}

	/**
	 * Given the graph, and two paths of the two reads, combine them into a single path
	 * @param graph
	 * @param path1
	 * @param path2
	 * @param dijkstraDis
	 * @return
	 */
	private static PairPath combinePaths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			List<Integer> path1, List<Integer> path2, DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis) {
		Integer firstV1 = path1.get(0);
		Integer lastV1 = path1.get(path1.size()-1);
		Integer firstV2 = path2.get(0);
		Integer lastV2 = path2.get(path2.size()-1);
		PairPath path  = new PairPath();

		if (path1.containsAll(path2))
		{
			path.setPath1(path1);
			//			return path;
		}
		else if (path2.containsAll(path1))
		{
			path.setPath2(path2);
			//			return path;
		}

		//path1 --> path2
		else if (isAncestral(getSeqVertex(graph, lastV1),getSeqVertex(graph, firstV2),dijkstraDis)>0)
		{
			path.setPath1(path1);
			path.setPath2(path2);
			//			return path;
		}

		//path2 --> path1
		else if (isAncestral(getSeqVertex(graph, lastV2),getSeqVertex(graph, firstV1),dijkstraDis)>0)
		{
			path.setPath1(path2);
			path.setPath2(path1);
			//			return path;
		}

		else if (isAncestral(getSeqVertex(graph, firstV2),getSeqVertex(graph, firstV1),dijkstraDis)==0 && 
				isAncestral(getSeqVertex(graph, lastV2),getSeqVertex(graph, lastV1),dijkstraDis)==0)
		{
			//there is no consistent path between read1 and read2
			//			return path;
		}

		//path1(partial) -> path2
		else if (isAncestral(getSeqVertex(graph, firstV1),getSeqVertex(graph, firstV2),dijkstraDis)>0 && 
				path1.indexOf(firstV2)>=0)
		{
			int i = path1.indexOf(firstV2);
			path.setPath1(path1.subList(0, i));
			path.addToPath1(path2);
		}

		//path2(partial) -> path1
		else if (isAncestral(getSeqVertex(graph, firstV2),getSeqVertex(graph, firstV1),dijkstraDis)>0 &&
				path2.indexOf(firstV1)>=0)
		{
			int i = path2.indexOf(firstV1);
			path.setPath1(path2.subList(0, i));
			path.addToPath1(path1);
		}

		if (path.getPath1().isEmpty() && !path.getPath2().isEmpty())
			path.movePath2To1();

		return path;

	}

	/** 
	 * using the dfs discovery and finishing times, figure out which is ancestral to which
	 * @param v1
	 * @param v2
	 * @param dijkstraDis 
	 * @return 1 if v1 is ancestral to v2
	 * @return -1 if v2 is ancestral to v1
	 * @return 0 if there is no path (v1,v2) and no path (v2, v1)
	 */
	private static int isAncestral(SeqVertex v1, SeqVertex v2, DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis) {

		// v1 ---> v2
		if (dijkstraDis.getDistance(v1, v2)!=null)
			return 1;
		// v2 ---> v1
		if (dijkstraDis.getDistance(v2, v1)!=null)
			return -1;

		return 0;
	}

	/**
	 * Count how many vertices we have with in degree >1 & out degree >1
	 * @param graph
	 * @return
	 */
	private static int countNumOfXstructures(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		int res = 0;

		for (SeqVertex v : graph.getVertices())
		{
			if (graph.inDegree(v)>1 && graph.outDegree(v)>1)
				res++;
		}

		return res;
	}

	/**
	 * Given the graph and the hash with all reads, find all probable paths from S to T.
	 * @param graph
	 * @param comp 
	 * @param combinedReadHash
	 * @param dijkstraDis 
	 * @param dijkstraDisWoVer 
	 */
	private static Pair<HashMap<List<Integer>, Pair<Integer>>> getAllProbablePaths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Set<SeqVertex> comp, HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash, 
			DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis, 
			DijkstraDistanceWoVer<SeqVertex,SimpleEdge> dijkstraDisWoVer) {

		HashMap<SeqVertex,List<List<Integer>>> Paths = new HashMap<SeqVertex,List<List<Integer>>>();
		HashMap<List<Integer>,HashMap<PairPath,Integer>> PathReads = new HashMap<List<Integer>,HashMap<PairPath,Integer>>();
		HashMap<List<Integer>,Boolean> Extensions = new HashMap<List<Integer>,Boolean>();

		HashMap<List<Integer>,Pair<Integer>> FinalPaths_diff = new HashMap<List<Integer>,Pair<Integer>>();
		HashMap<List<Integer>,Pair<Integer>> FinalPaths_all = new HashMap<List<Integer>,Pair<Integer>>();

		//initiation
		getTopologicalOrder(graph);
		ROOT.setDFS_FinishingTime(Integer.MAX_VALUE);
		T_VERTEX.setDFS_FinishingTime(-1);

		SeqVertexFinishTimeComparator finishingTimeComparator = new SeqVertexFinishTimeComparator();
		PriorityQueue<SeqVertex> C = new PriorityQueue<SeqVertex>(comp.size(),finishingTimeComparator  );

		C.add(ROOT);
		List<Integer> tmpL = new ArrayList<Integer>();
		tmpL.add(ROOT.getID());
		ArrayList<List<Integer>> tmpPathList = new ArrayList<List<Integer>>();
		tmpPathList.add(tmpL);
		Paths.put(ROOT, tmpPathList);
		SeqVertex v;

//		debugMes(""+combinedReadHash,15);
		while (!C.isEmpty())
		{
			debugMes("getPathsSize(Paths)="+getPathsSize(Paths)+" PathReads.size()="+ PathReads.size()+
					" Extensions.size()="+Extensions.size()+" FinalPaths.size()="+FinalPaths_diff.size(),20);
			if (getPathsSize(Paths)> MAX_PATH_SIZE)
			{
				debugMes("paths size has increased too much",10);
				return null;
			}	
			v = C.poll(); 
			debugMes("C = "+C,20);
			debugMes("the next node in the queue is "+v.getID(),20);
			HashMap<PairPath,Integer> readsStartingAtV = combinedReadHash.get(v.getID());

			//try to combine paths that are very similar and reach v
//			combineSimilarPathsThatEndAtV(graph,v,Paths,PathReads,Extensions);

			// go over all paths of P[v], add all reads that start at v
			for (List<Integer> path : Paths.get(v))
			{
				if (!PathReads.containsKey(path))
					PathReads.put(path, new HashMap<PairPath,Integer>());

				if (readsStartingAtV!=null && !readsStartingAtV.isEmpty())
				{
					debugMes("adding the reads " +readsStartingAtV +" to the path "+ path, 20);
					PathReads.get(path).putAll(readsStartingAtV);
				}

				//keep track of all extensions
				Extensions.put(path, false);

			}
			// go over all descendants of v
			for (SeqVertex u : graph.getSuccessors(v))
			{
				boolean vExtendedToU = false;
				for (List<Integer> path : Paths.get(v))
				{

					HashMap<PairPath,Integer> readsOfPathUntilV = PathReads.get(path);

					//					if (lastTripletHasEnoughReadSupport(readsOfPathUntilV,w_Ver,v,u,graph,dijkstraDisWoVer))
					if (pathHasEnoughReadSupport(readsOfPathUntilV,path,u,graph,dijkstraDisWoVer))

					{
						// add [path,u] to paths of u
						if (!Paths.containsKey(u))
							Paths.put(u, new ArrayList<List<Integer>>());

						List<Integer> pathWu = new ArrayList<Integer>();
						pathWu.addAll(path);
						pathWu.add(u.getID());
						debugMes("adding the path " +pathWu +" to the paths of "+ u.getID()+": "+Paths.get(u), 20);
						Paths.get(u).add(pathWu);

						//update reads of [path,u]
						updateReadsOfPath(PathReads,pathWu,readsOfPathUntilV,u.getID(),graph,dijkstraDis);

						//update extension
						Extensions.put(path, true);
						vExtendedToU = true;
					}

				}
				if (!C.contains(u))
				{
					debugMes(u.getID()+" was added to the queue",20);
					C.add(u);

				}
				//if v didn't extend to u, and we have an edge there, add (v,u) as a new path
				if (!vExtendedToU)
				{
					debugMes("the edge (v-u) was not used in any extension: "+v.getID()+"->"+u.getID(),15);
					if (!Paths.containsKey(u))
						Paths.put(u, new ArrayList<List<Integer>>());
					List<Integer> vuPath = new ArrayList<Integer>();
					vuPath.add(v.getID());
					vuPath.add(u.getID());

					Paths.get(u).add(vuPath);

					//add the reads
					if (!PathReads.containsKey(vuPath))
						PathReads.put(vuPath, new HashMap<PairPath,Integer>());

					if (readsStartingAtV!=null && !readsStartingAtV.isEmpty())
					{
						debugMes("adding the reads " +readsStartingAtV +" to the path "+ vuPath, 20);
						PathReads.get(vuPath).putAll(readsStartingAtV);
						updateReadsOfPath(PathReads,vuPath,readsStartingAtV,u.getID(),graph,dijkstraDis);

					}


				}
				//try to combine paths that are very similar and reach u
//				combineSimilarPathsThatEndAtV(graph,v,Paths,PathReads,Extensions);

			}

			//report the paths that were not extended AND remove them from Paths
			List<List<Integer>> removePaths = new ArrayList<List<Integer>>();
			for (List<Integer> path : Paths.get(v))
			{
				SeqVertex lastV = getSeqVertex(graph, path.get(path.size()-1));
				if (!lastV.equals(T_VERTEX) && !Extensions.get(path))
				{
					if (getSeqPathLength(graph,path)>MIN_OUTPUT_SEQ)
					{
						FinalPaths_diff.put(path,new Pair<Integer>(getSuppCalculation(PathReads.get(path)),0));
						debugMes("the unextended path: "+path+" was added to the final paths, with "+getSuppCalculation(PathReads.get(path)) +" support",1);
					} 
					removePaths.add(path);
				}
			}

			for (List<Integer> path : removePaths)
			{
				Paths.get(v).remove(path);
				Extensions.remove(path);
				//				PathReads.remove(path); // can't remove the pathReads, since we're using it later for expression profiling.
			}
		}

		for (List<Integer> path : Paths.get(T_VERTEX))
		{
			if (getSeqPathLength(graph,path)>MIN_OUTPUT_SEQ)
			{
				FinalPaths_diff.put(path,new Pair<Integer>(getSuppCalculation(PathReads.get(path)),0));
				if (path.get(0).intValue() == ROOT.getID())
					debugMes("the finished path: "+ path+" was added to the final paths, with "+getSuppCalculation(PathReads.get(path))+" support",15);
				else
					debugMes("the finished (from middle unextended) path: "+ path+" was added to the final paths, with "+getSuppCalculation(PathReads.get(path)) +" support",15);
			}
		}

		//remove similar final paths
		
//		combineSimilarPaths(graph,FinalPaths_diff,PathReads);

		// calc expression better
		calcExpressionOfFinalPaths(FinalPaths_diff,PathReads);


		return new Pair<HashMap<List<Integer>, Pair<Integer>>>(FinalPaths_diff,FinalPaths_all);

	}


	/**
	 * given these paths, and reads, re-calc the FPKM of each path
	 * @param FinalPaths
	 * @param PathReads
	 */
	private static void calcExpressionOfFinalPaths(
			HashMap<List<Integer>, Pair<Integer>> FinalPaths,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {
		HashMap<PairPath,Pair<Integer>> ReadPerPathCounts = new HashMap<PairPath, Pair<Integer>>();
		for (List<Integer> path : FinalPaths.keySet())
			for (PairPath read : PathReads.get(path).keySet())
			{
				Integer count = PathReads.get(path).get(read);
				if (!ReadPerPathCounts.containsKey(read))
					ReadPerPathCounts.put(read,new Pair<Integer>(1,count));
				else
					ReadPerPathCounts.put(read,new Pair<Integer>(ReadPerPathCounts.get(read).getFirst()+1,count));
			}

		for (List<Integer> path : FinalPaths.keySet())
		{
			Integer supp = 0;
			Integer totalCounts = 0;
			for (PairPath read : PathReads.get(path).keySet())
			{
				Integer numPaths = ReadPerPathCounts.get(read).getFirst();
				Integer count = ReadPerPathCounts.get(read).getSecond();
				supp += count/numPaths;
				totalCounts += count;
			}
			FinalPaths.put(path, new Pair<Integer>(totalCounts,supp));
		}
		
	}



//
//	/**
//	 * Go over all final paths, and combine those that are too similar.
//	 * @param graph
//	 * @param FinalPaths
//	 * @param PathReads
//	 */
//	private static void combineSimilarPaths(
//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
//			HashMap<List<Integer>, Pair<Integer>> FinalPaths,
//			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {
//		
//		List<List<Integer>> removeSimilarPaths = new ArrayList<List<Integer>>();
//
//		
//		Iterator<List<Integer>> i1,i2;
//		String path1S="", path2S="";
//		int index1;
//		for (i1=FinalPaths.keySet().iterator() ; i1.hasNext() ; )
//		{
//			List<Integer> path1 = i1.next();
//			path1S = getPathSeq(graph, path1);
//			
//			boolean gotToi1 = false;
//			for (i2=FinalPaths.keySet().iterator() ; i2.hasNext() ; )
//			{
//				List<Integer> path2 = i2.next();
//				
//				while (!gotToi1 && i2.hasNext())
//				{
//					if (path2.equals(path1))
//						gotToi1 = true;
//					path2 = i2.next();
//
//				}
//				
//				if (path2.equals(path1))
//					break;
//				
//				index1=path1S.length();
//
//				path2S = getPathSeq(graph, path2);
//
//				boolean noOverlappingVers = true;
//				int v1 = -1,v2 = -1, index2=-1;
//				for (int j1=path1.size()-1; j1>0 && noOverlappingVers ; j1--)
//				{
//					v1 = path1.get(j1);
//					index1 -= getSeqVertex(graph, v1).getName().length();
//					if (v1!=T_VERTEX.getID())
//					{
//						index2=path2S.length();
//						for (int j2=path2.size()-1; j2>0 && noOverlappingVers ; j2--)
//						{
//							v2 = path2.get(j2);
//							index2 -= getSeqVertex(graph, v2).getName().length();
//
//							if (v1==v2)
//								//update noOverlappingVers, so we'll get out of the loop:
//								noOverlappingVers = false;
//						}
//					}
//				}
//
//				if (!noOverlappingVers) //check only paths that share vertices
//				{
//					index1 += getSeqVertex(graph, v1).getName().length();
//					index2 += getSeqVertex(graph, v2).getName().length();
//					debugMes("checking paths: "+path1+ 
//							"(len="+path1S.length()+") and "+path2+"(len="+path2S.length()+")",15);
//
//					if (path1.lastIndexOf(T_VERTEX.getID())==-1)
//						index1--;
//					
//					if (path2.lastIndexOf(T_VERTEX.getID())==-1)
//						index2--;
//					if (twoPathsAreTooSimilar(path1S, path2S, index1, index2))
//					{
//						debugMes("they are too similar!",15);	
//						//remove the shorter path
//						removeTheShorterPath(path1S,path2S,path1,path2,removeSimilarPaths,PathReads);
//					}
//				}
//			}
//		}
//		
//		for (List<Integer> path2Remove : removeSimilarPaths)
//		{
//			debugMes("The final path "+path2Remove+" was removed becuase it was too close to another path",15);
//			FinalPaths.remove(path2Remove);
//		}
//	}
//
//
//
//	/**
//	 * check for similar paths that end at V, and start at different nodes
//	 * remove the shortest of the two
//	 * @param graph
//	 * @param v
//	 * @param Paths
//	 * @param PathReads
//	 * @param Extensions
//	 */
//	private static void combineSimilarPathsThatEndAtV(
//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,  SeqVertex v,
//			HashMap<SeqVertex, List<List<Integer>>> Paths,
//			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads,
//			HashMap<List<Integer>, Boolean> Extensions) {
//
//		List<List<Integer>> removeSimilarPaths = new ArrayList<List<Integer>>();
//		String path1S="", path2S="";
//		Iterator<List<Integer>> i1, i2;
//		int index1, index2;
//		for (i1=Paths.get(v).iterator() ; i1.hasNext() ; )
//		{
//			List<Integer> path1 = i1.next();
//			path1S = getPathSeq(graph, path1);
//			index1 = path1S.length()-1;
//
//			boolean gotToi1 = false;
//			for (i2=Paths.get(v).iterator() ; i2.hasNext() ; )
//			{
//				List<Integer> path2 = i2.next();
//				
//				while (!gotToi1 && i2.hasNext())
//				{
//					if (path2.equals(path1))
//						gotToi1 = true;
//					path2 = i2.next();
//
//				}
//				
//				if (path2.equals(path1))
//					break;
//
//				path2S = getPathSeq(graph, path2);
//				index2 = path2S.length()-1;
//				debugMes("checking for similarity the two paths: "+path1+ 
//						"(len="+path1S.length()+");"+path2+"(len="+path2S.length()+")",15);
//
//				if (twoPathsAreTooSimilar(path1S,path2S,index1,index2))
//				{
//					debugMes("they are too similar!",15);	
//					//remove the shorter path
//					removeTheShorterPath(path1S,path2S,path1,path2,removeSimilarPaths,PathReads);
//				}
//			}
//		}
//
//		for (List<Integer> path2Remove : removeSimilarPaths)
//		{
//			debugMes("The path "+path2Remove+" was removed becuase it was too close to another path",15);
//
//			Paths.get(v).remove(path2Remove);
//			Extensions.remove(path2Remove);
//		}
//
//
//
//	}

//
//	/** 
//	 * compare the sequences of the two paths, and return true if they are more than MAX_PATH_IDENTITY. 
//	 * @param path1s
//	 * @param path2s
//	 * @return
//	 */
//	private static boolean twoPathsAreTooSimilar(String path1s,
//			String path2s,int index1, int index2) {
//		int numMatch = 0, numMM = 0,longestMMstretch = 0;
//		boolean prevMM = true;
//		int curMMstretch = 0;
//		// zipper alignment from the given indices forward
//		for (int i1=index1+1,i2=index2+1 ; i1<path1s.length() && i2<path2s.length(); i1++, i2++)
//		{
//			if (path1s.charAt(i1) == path2s.charAt(i2)) //FIXME - should we use the degenrate code here?
//			{
//				numMatch++;
//				prevMM = false;
//			}
//			else
//			{
//				numMM++;
//				if (prevMM)
//					curMMstretch++;
//				else
//				{
//					if (curMMstretch>longestMMstretch)
//						longestMMstretch = curMMstretch;
//					curMMstretch=1;
//					prevMM = true;
//				}
//
//
//			}
//		}
//		if (curMMstretch>longestMMstretch)
//			longestMMstretch = curMMstretch;
//
//		// zipper alignment from the given indices backwards
//		curMMstretch = 0;
//		prevMM = true;
//		for (int i1=index1,i2=index2 ; i1>=0 && i2>=0; i1--, i2--)
//		{
//			if (path1s.charAt(i1) == path2s.charAt(i2)) //FIXME - should we use the degenrate code here?
//			{
//				numMatch++;
//				prevMM = false;
//			}
//			else
//			{
//				numMM++;
//				if (prevMM)
//					curMMstretch++;
//				else
//				{
//					if (curMMstretch>longestMMstretch)
//						longestMMstretch = curMMstretch;
//					curMMstretch=1;
//					prevMM = true;
//				}
//
//
//			}
//		}
//		if (curMMstretch>longestMMstretch)
//			longestMMstretch = curMMstretch;
//
//		int shortestLen = numMM+numMatch;
//		return isThisTooSimilar(numMM,longestMMstretch,shortestLen);
//	}
//
//	/**
//	 * given all the params, decide if the two seqs are too similar
//	 * FIXME - find a better criteria.
//	 * @param numMM - number of mismatches
//	 * @param longestMMstretch
//	 * @param shortestLen
//	 * @return
//	 */
//	private static boolean isThisTooSimilar(int numMM, 
//			int longestMMstretch, int shortestLen) {
//		double identity = ((shortestLen-numMM)/(double)(shortestLen))*100;
//		DecimalFormat df = new DecimalFormat("#.##");
//
//		debugMes("numMM = "+numMM+" identity = "+df.format(identity)+"% and logest MM stretch is "+longestMMstretch,15);
//
//		double minialMM = MINIMAL_MM_BETWEEN_PATHS;
////		if (shortestLen<10)
////			minialMM = 1;
//		
//		debugMes("the two paths have these stats: numMM="+numMM+
//				" numMatches="+(shortestLen-numMM)+" longestMMstretch="+longestMMstretch+" identity="+df.format(identity)+"%",15);
//		
//		return 
////		longestMMstretch < LONGEST_MM_STRETCH || //longest patch of mismatches is short enough -> similar paths 
//		numMM<minialMM ||                       //number of mismatches is low -> similar seqs
//		identity >= MAX_PATH_IDENTITY;        //identity is too high -> similar seqs
//	}
//
//
//
//
//	/**
//	 * given two paths (and their seqs) remove the shorter path, and add its reads to the other one.
//	 * if the are equal in length, remove the lighter one.
//	 * @param path1S
//	 * @param path2S
//	 * @param path1
//	 * @param path2
//	 * @param removeSimilarPaths
//	 * @param PathReads
//	 */
//	private static void removeTheShorterPath(String path1S, String path2S,
//			List<Integer> path1, List<Integer> path2, List<List<Integer>> removeSimilarPaths,
//			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {
//		List<Integer> path2remove,path2keep;
//		if (path1S.length() > path2S.length())
//		{
//			path2remove = path2;
//			path2keep = path1;
//		}
//		else if (path1S.length() < path2S.length())
//		{
//			path2remove = path1;
//			path2keep = path2;
//		}
//		else // they are equal
//		{
//			int sum1=0,sum2=0;
//			for (Integer s : PathReads.get(path1).values())
//				sum1+=s;
//			for (Integer s : PathReads.get(path2).values())
//				sum2+=s;
//
//			if (sum1<=sum2)
//			{
//				path2remove = path1;
//				path2keep = path2;
//			} 
//			else
//			{
//				path2remove = path2;
//				path2keep = path1;
//			}
//
//		}
//			
//		debugMes("removing path "+path2remove+" and keeping path "+path2keep,15);
//
//		if (!removeSimilarPaths.contains(path2remove))
//			removeSimilarPaths.add(path2remove);
//		if (PathReads.get(path2remove)!=null)
//		{
//			if (PathReads.get(path2keep)==null)
//				PathReads.put(path2keep, new HashMap<PairPath,Integer>());
//
//			PathReads.get(path2keep).putAll(PathReads.get(path2remove));
//			PathReads.remove(path2remove);
//		}								
////		}else
////		{
////			if (!removeSimilarPaths.contains(path1))
////				removeSimilarPaths.add(path1);
////			if (PathReads.get(path1)!=null)
////			{
////				if (PathReads.get(path2)==null)
////					PathReads.put(path2, new HashMap<PairPath,Integer>());
////
////				PathReads.get(path2).putAll(PathReads.get(path1));
////				PathReads.remove(path1);
////			}
////
////		}		
//	}




	/** 
	 * Given this path, ask whether it has enough support, either by last triplet, or by length
	 * @param readsOfPathUntilV - reads of this path, so far
	 * @param path - the path so far
	 * @param u - the extension to the path
	 * @param graph
	 * @param dijkstraDisWoVer
	 * @return
	 */
	private static boolean pathHasEnoughReadSupport(
			HashMap<PairPath, Integer> readsOfPathUntilV, List<Integer> path,
			SeqVertex u, DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer) {


		List<Integer> subPath = new ArrayList<Integer>();
		subPath.add(0, u.getID());
		if (USE_TRIPLETS)
		{
			SeqVertex v = getSeqVertex(graph,path.get(path.size()-1));
			if (v.getName().length() >= MAX_PAIR_DISTANCE)
				return true;
			subPath.add(0, v.getID());
			if (path.size()>1)
				subPath.add(0,path.get(path.size()-2));
			return subPathHasEnoughReadSupport(readsOfPathUntilV,subPath,graph,dijkstraDisWoVer);
		}else{
			int lookBack = MAX_PAIR_DISTANCE -50; 
			int lenSoFar = u.getName().length();
			for (int j = path.size()-1 ; j>=0 && lenSoFar < lookBack; j--){
				SeqVertex vLast = getSeqVertex(graph, path.get(j));
				subPath.add(0, vLast.getID());
				lenSoFar += vLast.getName().length();
			}
			return subPathHasEnoughReadSupport(readsOfPathUntilV,subPath,graph,dijkstraDisWoVer);
		}
	}



	/**
	 * Check that the given sub-path has N supporting reads or more.
	 * A supporting read is a read that enforces this triplet
	 * @param readsOfPathUntilV
	 * @param subPath
	 * @param graph
	 * @param dijkstraDisWoVer
	 * @return
	 */
	private static boolean subPathHasEnoughReadSupport(
			HashMap<PairPath, Integer> readsOfPathUntilV,
			List<Integer> subPath,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer) {

		int numberReadsSupporting = 0;
		for (PairPath pPath : readsOfPathUntilV.keySet())
		{
			boolean thisReadOK = true;
			for (Integer vTempID : subPath)
				if (thisReadOK)
					thisReadOK = thisReadOK && 	
					readEnforcesVertex(graph, dijkstraDisWoVer, pPath, getSeqVertex(graph, vTempID));

			if (thisReadOK)
			{
				numberReadsSupporting+=readsOfPathUntilV.get(pPath);
				debugMes("the read "+pPath+"("+readsOfPathUntilV.get(pPath)+") enforces the sub-path ("+subPath+")",20);
			} else
				debugMes("the read "+pPath+"("+readsOfPathUntilV.get(pPath)+") does not enforce the sub-path ("+subPath+")",20);
		}

		boolean res = (numberReadsSupporting>=MIN_TRIPLET_SUPPORT_THR);
		if (res)
			debugMes("the sub-path ("+subPath+") has PASSED",20);
		else
			debugMes("the sub-path ("+subPath+") has NOT PASSED",15);

		return res;	
	}






	/**
	 * given the graph and a list of nodes, calc the length of the seq of this path
	 * @param graph
	 * @param path
	 * @return
	 */
	private static int getSeqPathLength(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, List<Integer> path) {
		int len = 0;
		for (Integer vid : path)
			len +=getSeqVertex(graph, vid).getName().length();

		return len;
	}






	private static int getPathsSize(
			HashMap<SeqVertex, List<List<Integer>>> paths) {
		int res = 0;
		for (SeqVertex key : paths.keySet())
		{
			res+=paths.get(key).size();
		}
		return res;
	}


	private static boolean areTwoNucleotidesEqual(String n1, String n2)
	{
		if (n1.equals(n2))
			return true;

		if (USE_DEGENERATE_CODE && 
				((DEGENERATE_CODE_REV.containsKey(n1) && DEGENERATE_CODE_REV.get(n1).contains(n2)) ||
						(DEGENERATE_CODE_REV.containsKey(n2) && DEGENERATE_CODE_REV.get(n2).contains(n1))))
			return true;

		return false;
	}
	private static String getDegenerateRepresentation(String key) throws Exception {
		if (DEGENERATE_CODE.containsKey(key))
			return DEGENERATE_CODE.get(key);
		else
			throw new Exception("the letters "+key+" do not have a degenerate representation\n");
	}



	/**
	 * sum the counts of all the reads in this hash
	 * @param readHash
	 * @return
	 */
	private static Integer getSuppCalculation(HashMap<PairPath, Integer> readHash) {
		Integer res = 0;
		for (PairPath key : readHash.keySet())
			res = res + readHash.get(key);
		return res;
	}



	/**
	 * Given the new path (with u), and the set of reads that supported the path until v
	 * update the set of reads that support the new path
	 * @param PathReads 
	 * @param pathWu
	 * @param readsOfPathUntilV
	 * @param i 
	 * @param dijkstraDis 
	 * @param graph 
	 */
	private static void updateReadsOfPath(HashMap<List<Integer>,HashMap<PairPath,Integer>> PathReads, List<Integer> pathWu,
			HashMap<PairPath, Integer> readsOfPathUntilV, Integer i, DirectedSparseGraph<SeqVertex, SimpleEdge> graph, DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis) {

		if (!PathReads.containsKey(pathWu))
			PathReads.put(pathWu, new HashMap<PairPath,Integer>());

		for (PairPath pPath : readsOfPathUntilV.keySet())
		{
			if (!PathReads.get(pathWu).containsKey(pPath))  // only if this read doesn't exist in the PathReads for this pathWu
				// if this read is consistent with pathWu, then add it
				if (readIsConsistentWithNode(pPath,i,graph,dijkstraDis))
				{
					debugMes("read "+pPath+" is consistent with "+i, 20);
					PathReads.get(pathWu).put(pPath,readsOfPathUntilV.get(pPath));
				}else{
					debugMes("read "+pPath+" is not consistent with "+i, 20);
				}
		}

	}

	/**
	 * Check whether the pairPath is consistent with the node i
	 * @param pPath
	 * @param i
	 * @param graph
	 * @param dijkstraDis
	 * @return
	 */
	private static boolean readIsConsistentWithNode(PairPath pPath, Integer i,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis) {

		if (pPath.containsID(i) || i<0)
			return true;
		SeqVertex vI = getSeqVertex(graph, i);
		SeqVertex firstV = getSeqVertex(graph, pPath.getFirstID());
		// i --> firstV
		if (isAncestral(vI, firstV, dijkstraDis)>0)
			return true;

		SeqVertex lastV = getSeqVertex(graph, pPath.getLastID());
		// lastV --> i
		if (isAncestral(lastV,vI,dijkstraDis)>0)
			return true;

		if (pPath.hasSecondPath())
		{
			//last of first path
			lastV = getSeqVertex(graph, pPath.getLastID_path1());
			//first of second path
			firstV = getSeqVertex(graph, pPath.getFirstID_path2());

			// lastV --> i --> firstV
			if (isAncestral(lastV,vI,dijkstraDis)>0 && isAncestral(vI, firstV, dijkstraDis)>0)
				return true;
		}
		return false;
	}

	/**
	 * Check whether there are at least N reads enforcing 
	 * @param graph
	 * @param dijkstraDis
	 * @param pPath
	 * @param v
	 * @return
	 */
	private static boolean readEnforcesVertex(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer,
			PairPath pPath, SeqVertex v) {
		if (v==null || pPath.containsID(v.getID()) || v.equals(ROOT) || v.equals(T_VERTEX))
			return true;

		SeqVertex firstV = getSeqVertex(graph, pPath.getFirstID());
		if (dijkstraDisWoVer.getDistanceWoVer(ROOT, firstV,v)==null)
			return true;

		SeqVertex lastV = getSeqVertex(graph, pPath.getLastID());
		if (dijkstraDisWoVer.getDistanceWoVer(lastV, T_VERTEX,v)==null)
			return true;

		if (pPath.hasSecondPath())
		{
			//last of first path
			lastV = getSeqVertex(graph, pPath.getLastID_path1());
			//first of second path
			firstV = getSeqVertex(graph, pPath.getFirstID_path2());
			if (dijkstraDisWoVer.getDistanceWoVer(lastV,firstV,v)==null)
				return true;
		}
		return false;
	}

	/**
	 * Return the reads, hashed by their starting vertex
	 * @param graph
	 * @param filename
	 * @param originalVerIDsMapping
	 * @param rootIDs
	 * @return
	 * @throws IOException
	 */
	private static HashMap<String, List<Read>> getReadStarts(DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			String filename,
			HashMap<Integer, LocInGraph> originalVerIDsMapping, Vector<Integer> rootIDs) throws IOException {
		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 
		HashMap<String, List<Read>> readNameHash = new HashMap<String, List<Read>>();
		String l = fileB.readLine(); // read header of component
		int numReadsNotMapped = 0;
		int numReadsMapped = 0;

		while (fileB.ready())

		{
			l = fileB.readLine();
			if (l.isEmpty())
				continue;
			//			Component 0
			//			>@42MRYAAXX100104:7:100:1000:103#0      11      101393  36      101418          GAAAGACTGTCACCCTTGAGGTGGAGTCCTCTGACACTATTGACAATGTCAAGAGCAAAATCCAAGACAAGGAAGG
			String[] fields = l.split("\t");

			List<Integer> pathIDS = null;
			Read r = new Read();
			pathIDS = readSingleRead(fields,originalVerIDsMapping,graph,r,false);

			if (pathIDS==null || (pathIDS!=null && pathIDS.isEmpty()))
			{
				numReadsNotMapped++;
			}else
			{

				//add to readNameHash
				if (!readNameHash.containsKey(r.getName()))
					readNameHash.put(r.getName(), new ArrayList<Read>());

				readNameHash.get(r.getName()).add(r);
				numReadsMapped++;

			}
		}	
		//		debugMes("number of reads not found in graph = "+numReadsNotMapped +" of a total of "+(numReadsNotMapped+numReadsMapped),10);
		debugMes("number of reads found = "+numReadsMapped+" (from total of "+(numReadsNotMapped+numReadsMapped)+") which came from "+ readNameHash.keySet().size() + " pairs",10);

		if (numReadsNotMapped > .5*(numReadsNotMapped+numReadsMapped))
			debugMes("PROBLEM: less than half of the reads were mapped to this graph ("+numReadsMapped+"/"+(numReadsNotMapped+numReadsMapped)+")",10);
		return readNameHash;
	}


	/**
	 * given this read, try and map it to the graph. if rev= true, do it in reverse.
	 * @param fields
	 * @param originalVerIDsMapping
	 * @param graph
	 * @param r 
	 * @param rev
	 * @return
	 */
	private static List<Integer> readSingleRead(String[] fields,
			HashMap<Integer, LocInGraph> originalVerIDsMapping,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, Read r, boolean rev) {
		List<Integer> pathIDS = new ArrayList<Integer>();
		LocInGraph fromV;
		Integer startInRead,endInRead,fromOrigV;

		String name;
		String seq;

		name = fields[0];
		if (name.endsWith("/1") || name.endsWith("/2") || name.endsWith("\1") || name.endsWith("\2"))
			name = name.substring(0, name.length()-2);
		//		if (name.equals(">SRR039231.19456143_FC42DB6AAXX:2:108:138:1484"))
		//			debugMes("stop here",10);

		startInRead = Integer.parseInt(fields[1]);
		endInRead = Integer.parseInt(fields[3])+K-1;

		fromOrigV = Integer.parseInt(fields[2]);
		fromV = originalVerIDsMapping.get(fromOrigV);
		seq = fields[6]; //there is an empty field before the seq.
		seq = seq.substring(startInRead, endInRead);

		if (fromV!=null)// && toV!=null)
		{
			pathIDS = findPathInGraph(graph,seq,fromV,name);
			if (!pathIDS.isEmpty())
				r.init(name,seq, fromV, startInRead, endInRead,pathIDS);
			//			else
			//				debugMes("read "+name+" was not mapped to graph. couldn't match sequence.",15);

		}else
			debugMes("read "+name+" was not mapped to graph. original node doesn't exist anymore ("+fromOrigV+")",15);

		return pathIDS;
	}



	//	/**
	//	 * reverse complement the given string
	//	 * @param string
	//	 * @return
	//	 */
	//	private static String RC_seq(String seq) {
	//		String res = "";
	//		for (int i=seq.length()-1; i>=0 ; i--)
	//			res = res.concat(REV.get(""+seq.charAt(i)));
	//		return res;
	//	}





	/**
	 * Given the graph, and the read, find the path of the read in the graph
	 * @param graph
	 * @param seq
	 * @param fromV
	 * @return
	 */
	private static List<Integer> findPathInGraph(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, String seq,
			LocInGraph fromV, String readName) {


		List<Integer> path = new ArrayList<Integer>();

		SeqVertex fromVer = getSeqVertex(graph, fromV.getNodeID());
		List<SeqVertex> continueVers = new ArrayList<SeqVertex>();
		continueVers.add(fromVer);
		debugMes("trying to start the mapping to node "+fromVer.getID(),20);
		Integer totalNumMM = 0;
		updatePathRecursively(path,graph,continueVers,seq,fromV.getIndexInNode(),totalNumMM,readName);
		return path;
	}


	private static void updatePathRecursively(List<Integer> path,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			List<SeqVertex> fromVers, String seq,int locInNode, Integer totalNumMM,String readName) {


		for (SeqVertex fromV : fromVers)
		{
			debugMes("trying to continue the mapping to node "+fromV.getID(),20);
			String verSeq = fromV.getName();
			int startI = locInNode;
			int j=0, i = startI;
			//			int marginForSkipping = 3;
			for (; i>=0 && i<verSeq.length() && j<seq.length() ; i++,j++)
			{
				String readLetter = ""+seq.charAt(j);
				String verLetter = ""+verSeq.charAt(i);
				if (!areTwoNucleotidesEqual(readLetter,verLetter)) 
				{
					//we have a mismatch
					totalNumMM++;
					if (totalNumMM>=MAX_MM_ALLOWED)
						break;

					if (seq.length()>j+1 && verSeq.length()>i+1)
					{
						j++;
						i++;
					}
				}
			}

			if (totalNumMM>=MAX_MM_ALLOWED)
			{
				debugMes("read "+readName+" has too many mismatches ("+totalNumMM+")",15);
				path.clear();
				break;
			} else if (j==seq.length())
			{
				path.add(fromV.getID());
				break;
			}else if (i==verSeq.length()) // move to the next ver
			{
				path.add(fromV.getID());

				List<SeqVertex> continueVers = new ArrayList<SeqVertex>();

				Collection<SimpleEdge> outE = graph.getOutEdges(fromV);
				List<SimpleEdge> outE_list = new ArrayList<SimpleEdge>();
				outE_list.addAll(outE);
				SimpleEdgeComparator edgeComp = new SimpleEdgeComparator();
				Collections.sort(outE_list, edgeComp);

				for (SimpleEdge e : outE_list)
					//				for (SeqVertex v2 : graph.getSuccessors(fromV))
				{
					SeqVertex v2 = graph.getDest(e);
					if (v2.getName().charAt(0)==seq.charAt(j))
						continueVers.add(v2);
				}
				updatePathRecursively(path,graph,continueVers,seq.substring(j),0,totalNumMM,readName);
				break;
			}else // the seq hasn't ended, and the vertex hasn't ended either, wrong mapping 
			{


			}
		}
	}

	/**
	 * create a hash that hold all the original vertices ids and the new ones 	
	 * @param graph
	 * @param rootIDs 
	 * @return the hash
	 */
	private static HashMap<Integer, LocInGraph> getOriginalVerIDsMappingHash(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {

		// clear double entries in the prevID list - not sure why they happen?
		for (SeqVertex v : graph.getVertices())
			v.clearDoubleEntriesToPrevIDs();


		HashMap<Integer, LocInGraph> hash = new HashMap<Integer,LocInGraph>();
		for (SeqVertex v : graph.getVertices())
		{

			Integer loc = 0;
			Integer vid = v.getID();
			// if the node id is new, than the real start is in the vector 
			if (vid>LAST_REAL_ID)
				loc = loc-1;
			else
			{
				debugMes("adding to "+vid+": Location of original node "+vid+" in index "+loc,20);
				hash.put(vid, new LocInGraph(vid,loc));

			}

			for (Vector<Integer> vec :  v.getPrevVerIDs())
			{
				loc++;
				for (Integer id : vec)
				{
					debugMes("adding to "+id+": Location of original node "+v.getID()+" in index "+loc,20);
					hash.put(id, new LocInGraph(v.getID(),loc));
				}
			}
		}
		return hash;
	}

	//
	//	/**
	//	 * given the filename, make a graph out of the connected components
	//	 * @param filename
	//	 * @param rootIDs 
	//	 * @return
	//	 * @throws IOException
	//	 */
	//	private static DirectedSparseGraph<SeqVertex, SimpleEdge> buildNewGraph(String filename, Vector<Integer> rootIDs) 
	//	throws IOException
	//	{
	//
	//		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 
	//		DirectedSparseGraph<SeqVertex, SimpleEdge> graph = 
	//			new DirectedSparseGraph<SeqVertex,SimpleEdge>();
	//		String l = fileB.readLine(); // read header of component
	//		Integer from, to;
	//		double supp;
	//		while (fileB.ready())
	//
	//		{
	//			l = fileB.readLine();
	//			//	0       -1      3       ATTGAAAGCAAGTTTTCTTCGAAT        0
	//			//	1       0       3       TTGAAAGCAAGTTTTCTTCGAATT        0
	//			//	to		from	supp	kmer							stam       
	//			String[] fields = l.split("\t");
	//			from = Integer.parseInt(fields[1]);
	//			to = Integer.parseInt(fields[0]);
	//			supp = Double.parseDouble((fields[2]));
	//
	//			if (supp < INITIAL_EDGE_ABS_THR)
	//				continue;
	//
	//			if (from>LAST_ID)
	//				LAST_ID = from;
	//
	//			if (to>LAST_ID)
	//				LAST_ID = to;
	//
	//			String kmer = fields[3];
	//			K = kmer.length();
	//			SeqVertex toV = getSeqVertex(graph, to);
	//			SeqVertex fromV = getSeqVertex(graph, from);
	//
	//			boolean isRoot = (from<0 || fromV==null);
	//
	//			if (isRoot)
	//			{
	//				if (toV==null)
	//				{
	//					toV = new SeqVertex(to, kmer,supp);
	//					graph.addVertex(toV);
	//					//					debugMes("added vertex " + toV + " with sequence: " + kmer);
	//					rootIDs.add(to);
	//				}
	//			}
	//			else
	//			{
	//				kmer = kmer.substring(kmer.length()-1, kmer.length());
	//				if (toV==null)
	//				{
	//					toV = new SeqVertex(to, kmer);
	//					graph.addVertex(toV);
	//					//					debugMes("added vertex " + toV + " with sequence: " + kmer);
	//				}
	//				SimpleEdge e = new SimpleEdge(supp); 
	//				graph.addEdge(e, fromV, toV);
	//				//				debugMes("added edge: " + fromV + "->" + toV + "(" + supp + ")"grep  );
	//
	//			}
	//
	//		}
	//		return graph;
	//	}



	/**
	 * go over the graph file, and count the in flow and out flow of each node
	 * @param firstLetter 
	 * @throws IOException 
	 */
	private static void preProcessGraphFile(String filename,
			HashMap<Integer, Integer> outFlow, HashMap<Integer, Integer> inFlow, HashMap<Integer,String> firstLetter) throws IOException {

		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 
		String l = fileB.readLine(); // read header of component
		Integer from, to, supp;

		while (fileB.ready())

		{
			l = fileB.readLine();
			//	0       -1      3       ATTGAAAGCAAGTTTTCTTCGAAT        0
			//	1       0       3       TTGAAAGCAAGTTTTCTTCGAATT        0
			//	to		from	supp	kmer							stam       
			String[] fields = l.split("\t");
			from = Integer.parseInt(fields[1]);
			to = Integer.parseInt(fields[0]);
			supp = Integer.parseInt((fields[2]));
			String kmer = fields[3];

			if (!outFlow.containsKey(from))
				outFlow.put(from, supp);
			else
				outFlow.put(from, outFlow.get(from)+supp);

			if (!inFlow.containsKey(to))
				inFlow.put(to, supp);
			else
				inFlow.put(to, inFlow.get(to)+supp);

			firstLetter.put(to,kmer.substring(0,1));

		}


	}


	/**
	 * given the filename, make a graph out of the connected components
	 * This time, keep the first letter of each kmer:
	 * keep the whole kmer, and then if there is an edge out, leave only first letter 
	 * @param filename
	 * @param rootIDs 
	 * @param inFlow in flow for all vertices
	 * @param outFlow out flow for all vertices
	 * @param firstLetter 
	 * @return
	 * @throws IOException
	 */
	private static DirectedSparseGraph<SeqVertex, SimpleEdge> buildNewGraphFirstLetter(String filename, 
			Vector<Integer> rootIDs, HashMap<Integer,Integer> outFlow, HashMap<Integer,Integer> inFlow, HashMap<Integer, String> firstLetter) 
			throws IOException
			{

		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 
		DirectedSparseGraph<SeqVertex, SimpleEdge> graph = 
			new DirectedSparseGraph<SeqVertex,SimpleEdge>();
		String l = fileB.readLine(); // read header of component
		Integer from, to;
		double supp;
		while (fileB.ready())

		{
			l = fileB.readLine();
			//	0       -1      3       ATTGAAAGCAAGTTTTCTTCGAAT        0
			//	1       0       3       TTGAAAGCAAGTTTTCTTCGAATT        0
			//	to		from	supp	kmer							stam       
			String[] fields = l.split("\t");
			from = Integer.parseInt(fields[1]);
			to = Integer.parseInt(fields[0]);
			supp = Double.parseDouble((fields[2]));
			if (supp < INITIAL_EDGE_ABS_THR )
				continue;

			// not using this now, since I want to fix the extremely high edges first
			//				if  (supp<inFlow.get(from)*EDGE_THR)
			//				{
			//					debugMes("edge ("+from+")->("+to+"):"+supp+" was not parsed since "+supp+"<"+inFlow.get(from)+"*"+EDGE_THR,20);
			//					continue;
			//				}	
			//			if (supp<outFlow.get(to)*EDGE_THR)
			//			{
			//				debugMes("edge ("+from+")->("+to+"):"+supp+" was not parsed since "+supp+"<"+outFlow.get(to)+"*"+EDGE_THR,20);
			//				continue;
			//			}

			if (from>LAST_ID)
				LAST_ID = from;

			if (to>LAST_ID)
				LAST_ID = to;

			String kmer = fields[3];
			K = kmer.length();
			SeqVertex toV = getSeqVertex(graph, to);
			SeqVertex fromV = getSeqVertex(graph, from);
			if (fromV==null && from>=0)
			{
				fromV = new SeqVertex(from,firstLetter.get(from)+""+kmer.substring(0,K-1));
				graph.addVertex(fromV);
			}

			boolean isRoot = (from<0 || fromV==null);

			if (isRoot)
			{
				if (toV==null)
				{
					toV = new SeqVertex(to, kmer,supp);
					graph.addVertex(toV);
					rootIDs.add(to);
				}
			}
			else
			{
				if (toV==null)
				{
					toV = new SeqVertex(to, kmer);
					graph.addVertex(toV);
				}
				SimpleEdge e = new SimpleEdge(supp); 
				graph.addEdge(e, fromV, toV);
			}

		}

		//Go over the whole graph, and if there edges coming out, leave only first letter
		for (SeqVertex v : graph.getVertices())
		{
			if (graph.outDegree(v)>0)
			{
				v.removeAllButFirstLetter();
			}
		}

		return graph;
			}



	/**
//	 * given the filename, make a graph out of the connected components
//	 * This time, keep the first letter of each kmer:
//	 * keep the whole kmer, and then if there is an edge out, leave only first letter 
//	 * @param filename
//	 * @param rootIDs 
//	 * @return
//	 * @throws IOException
//	 */
	//	private static DirectedSparseGraph<SeqVertex, SimpleEdge> buildNewGraphFirstLetter(String filename, 
	//			Vector<Integer> rootIDs) 
	//	throws IOException
	//	{
	//	
	//		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 
	//		DirectedSparseGraph<SeqVertex, SimpleEdge> graph = 
	//			new DirectedSparseGraph<SeqVertex,SimpleEdge>();
	//		String l = fileB.readLine(); // read header of component
	//		Integer from, to;
	//		double supp;
	//		while (fileB.ready())
	//	
	//		{
	//			l = fileB.readLine();
	//			//	0       -1      3       ATTGAAAGCAAGTTTTCTTCGAAT        0
	//			//	1       0       3       TTGAAAGCAAGTTTTCTTCGAATT        0
	//			//	to		from	supp	kmer							stam       
	//			String[] fields = l.split("\t");
	//			from = Integer.parseInt(fields[1]);
	//			to = Integer.parseInt(fields[0]);
	//			supp = Double.parseDouble((fields[2]));
	//	
	//			if (supp < INITIAL_EDGE_ABS_THR)
	//				continue;
	//	
	//			if (from>LAST_ID)
	//				LAST_ID = from;
	//	
	//			if (to>LAST_ID)
	//				LAST_ID = to;
	//	
	//			String kmer = fields[3];
	//			K = kmer.length();
	//			SeqVertex toV = getSeqVertex(graph, to);
	//			SeqVertex fromV = getSeqVertex(graph, from);
	//	
	//			boolean isRoot = (from<0 || fromV==null);
	//	
	//			if (isRoot)
	//			{
	//				if (toV==null)
	//				{
	//					toV = new SeqVertex(to, kmer,supp);
	//					graph.addVertex(toV);
	//					//					debugMes("added vertex " + toV + " with sequence: " + kmer);
	//					rootIDs.add(to);
	//				}
	//			}
	//			else
	//			{
	////				kmer = kmer.substring(kmer.length()-1, kmer.length());
	//				if (toV==null)
	//				{
	//					toV = new SeqVertex(to, kmer);
	//					graph.addVertex(toV);
	//					//					debugMes("added vertex " + toV + " with sequence: " + kmer);
	//				}
	//				SimpleEdge e = new SimpleEdge(supp); 
	//				graph.addEdge(e, fromV, toV);
	//				//				debugMes("added edge: " + fromV + "->" + toV + "(" + supp + ")"grep  );
	//	
	//			}
	//	
	//		}
	//		
	//		//Go over the whole graph, and if there edges coming out, leave only first letter
	//		for (SeqVertex v : graph.getVertices())
	//		{
	//			if (graph.outDegree(v)>0)
	//			{
	//				v.removeAllButFirstLetter();
	//			}
	//		}
	//		
	//		return graph;
	//	}
	private static void writeDotFile(DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			PrintStream p, String name)
	{
		writeDotFile(graph,p,name,false);
	}

	/**
	 * Write to dot file, where the list of paths are colored red -> blue
	 * @param graph
	 * @param p
	 * @param name
	 * @param vertices which vertices to print
	 */
	private static void writeDotFile(DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			PrintStream p, String name,boolean printFullSeq)
	{

		p.println("digraph "+name+"{");
		SeqVertex toVertex;
		int weight;

		//for each edge decide it's color
		for (SeqVertex vertex : graph.getVertices())
		{ //go over all vertices

			String verDesc = ""+vertex.getID()+" [label=\"";
			if (printFullSeq)
				verDesc = verDesc.concat(""+vertex.getLongtSeqWID() + "["+vertex.getName().length()+"]\"");
			//				verDesc = verDesc.concat(""+vertex.getShortSeqWID() + "["+vertex.getName().length()+"]{"+vertex.getPrevVerIDs()+"}\"");
			//				p.println(vertex.getID()+" [label=\"" + vertex.getShortSeqWID() + "["+vertex.getName().length()+"]\"]"); // prints also id
			else
				//				verDesc = verDesc.concat(""+vertex.getShortSeqWID() + "["+vertex.getName().length()+"]\"");
//				verDesc = verDesc.concat(""+vertex.getShortSeqWID() + "["+vertex.getName().length()+"]{"+vertex.getPrevVerIDs()+"}\"");
				verDesc = verDesc.concat(""+vertex.getShortSeqWID() + "["+vertex.getName().length()+"]\"");


			//				p.println(vertex.getID()+" [label=\"" + vertex.getShortSeqWID() + "["+vertex.getName().length()+"]\"]");

			if (vertex.getWeightAvg()>25)
				verDesc = verDesc.concat(" ,style=bold,color=\"#AF0000\"");

			verDesc = verDesc.concat("]");
			p.println(verDesc);


			for (SimpleEdge edge : graph.getOutEdges(vertex)) //get all edges of vertex->?
			{
				toVertex = graph.getDest(edge);

				weight = (int) Math.round(edge.getWeight());
				String edgeStyle = "[label="+ weight +"]";

				if (toVertex.equals(T_VERTEX) || vertex.equals(ROOT))
					edgeStyle = "[style=dotted,label="+ weight+"]";
				if (weight>20)
					edgeStyle = "[style=bold,label="+ weight +",color=\"#AF0000\"]";

				p.println(vertex.getID() + "->" + toVertex.getID() +edgeStyle);

			}
		}


		p.println("}");

	}




	//	/** 
	//	 * Print the given graph into the given stream in a dot format
	//	 * highlight given path by bold.
	//	 */
	//	private static void writeDotFile(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
	//			PrintStream p, String name, boolean fullGraph)
	//	{
	//		ArrayList<List<SimpleEdge>> pathList = new ArrayList<List<SimpleEdge>>(1);
	//		pathList.add(path);
	//		writeDotFile(graph, p, name,pathList,COLORS,fullGraph);
	//	}

	//	/**
	//	 * write dot file only of the given comp
	//	 * @param graph
	//	 * @param comp
	//	 * @param p
	//	 * @param graphName
	//	 */
	//	private static void writeDotFile(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
	//			Set<SeqVertex> comp, PrintStream p, String graphName) {
	//		writeDotFile(graph, p, graphName,comp);
	//
	//	}

	/**
	 * Compact the given graph:
	 * for each vertex, if degree out = degree in = 1, and nextVertexIn ==1, remove this vertex, and connect edges
	 * @param graph
	 */
	private static boolean compactLinearPaths(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		debugMes("=================\nCOMPACTING THE GRAPH\n=================",10);
		//compact vertices
		Vector<SeqVertex> removeVertices = new Vector<SeqVertex>();
		Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();
		boolean changed = false;
		for (SeqVertex v1 : graph.getVertices())
		{
			//			debugMes("looking at vertex: "+v1);
			while (!v1.equals(ROOT) && graph.outDegree(v1)==1 )
			{
				SimpleEdge e = null;
				for (SimpleEdge ei : graph.getOutEdges(v1))
					e = ei;
				SeqVertex v2 = graph.getDest(e);

				if (graph.inDegree(v2)!=1 || removeVertices.contains(v2) || v2.equals(T_VERTEX))
					break;
				debugMes("Found potential edge: "+e +" between "+v1 +" and "+v2,20);
				v1.concatVertex(v2, e.getWeight(),LAST_REAL_ID);
				debugMes("removing vertex "+v2+" was concatenated into "+v1,20);
				removeVertices.add(v2);
				changed = true;
				removeEdges.clear();
				for (SimpleEdge e2 : graph.getOutEdges(v2))
				{
					SeqVertex v3 = graph.getDest(e2);
					debugMes("Want to move edge " + e2 + "("+v2 +"->"+v3+") to ("+v1+"->"+v3,20);

					SimpleEdge newEdge = new SimpleEdge(e2);
					graph.addEdge(newEdge, v1, v3);

					removeEdges.add(e2);
				}

				for (SimpleEdge re : removeEdges)
				{
					debugMes("removing edge " + re + "("+graph.getSource(re) +"->"+graph.getDest(re)+")",20);
					graph.removeEdge(re);
				}
				debugMes("removing edge " + e + "("+v1 +"->"+v2+")",20);
				graph.removeEdge(e);

			}
		}
		//remove all vertices that we don't want
		for (SeqVertex v : removeVertices)
		{
			graph.removeVertex(v);
		}
		return changed;
	}



	/**
	 * remove light edges from the graph. return true if something has changed
	 * @param graph
	 * @return
	 */
	private static boolean removeLightEdges(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		boolean comp = false ; //removeLightCompEdges(graph);
		boolean in = removeLightInEdges(graph);
		boolean out = removeLightOutEdges(graph);
		boolean flow = removeLightFlowEdges(graph);
		return comp || in || out || flow;
	}


	/**
//	 * given the graph, go over each component, and calc the average coverage of that component
//	 * then go over all edges, and if they are too light, remove them.
//	 * @param graph
//	 * @return
//	 */
	//	private static boolean removeLightCompEdges(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
	//		debugMes("=================\nREMOVING LIGHT COMP FLOW EDGES\n=================",10);
	//
	//		boolean changed = false;
	//		Set<Set<SeqVertex>> comps = divideIntoComponents(graph);
	//		
	//		for (Set<SeqVertex> comp : comps)
	//		{
	//			Set<SimpleEdge> compEdges = getCompEdges(graph,comp);
	//			double curSum = 0;
	//			double curLen = 0;
	//			for (SeqVertex v : comp)
	//			{
	//				if (v.getName().length()>1)
	//				{
	//					curLen += v.getWeights().size();
	//					curSum += v.getWeightSum();
	//				}
	//				
	//				for (SimpleEdge e : compEdges)
	//				{
	//					curLen ++;
	//					curSum += e.getWeight();
	//				}
	//			}
	//			
	//			double avgComp = curSum/curLen;
	//			Set<SimpleEdge> removeE = new HashSet<SimpleEdge>();
	//			for (SimpleEdge e : compEdges)
	//				if (e.getWeight()<EDGE_THR*avgComp)
	//					removeE.add(e);
	//			
	//			for (SimpleEdge e : removeE)
	//			{
	//				changed = true;
	//				debugMes("removing low comp flow edge "+e+" from "+ graph.getSource(e)+" to "+graph.getDest(e),20);
	//				graph.removeEdge(e);
	//			}
	//		}
	//		return changed;
	//	}




	//
	//	/**
	//	 * Return all the edges of this component 
	//	 * @param graph
	//	 * @param comp
	//	 * @return
	//	 */
	//	private static Set<SimpleEdge> getCompEdges(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
	//			Set<SeqVertex> comp) {
	//		
	//		Set<SimpleEdge> edges = new HashSet<SimpleEdge>();
	//		for (SeqVertex v : comp)
	//		{
	//			for (SimpleEdge e : graph.getOutEdges(v))
	//				edges.add(e);
	//			
	//			for (SimpleEdge e : graph.getInEdges(v))
	//			{
	//				if (!comp.contains(graph.getSource(e)))
	//					edges.add(e);
	//			}
	//		}
	//		
	//		return edges;
	//	}






	/**
	 * Given a graph, go over all vertices and remove incoming or outgoing edges that do not match the flow (<2% coverage) see FLOW_THR
	 * When considering flow, this considers both the incoming and outgoing edges, but also the average node coverage.
	 * @param graph
	 * @return true if graph was changed.
	 */
	private static boolean removeLightFlowEdges(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		debugMes("=================\nREMOVING LIGHT FLOW EDGES\n=================",10);

		boolean changed = false;
		for (SeqVertex v : graph.getVertices())
		{
			if (graph.inDegree(v)==0 && graph.outDegree(v)==0)
				continue;

			//			if (v.getID()==38725)
			//			{
			//				System.err.println("FINDME - "+graph.getOutEdges(v));
			//			}

			int totalIn = 0, totalOut = 0;
			for (SimpleEdge e : graph.getInEdges(v))
				totalIn+=e.getWeight();

			for (SimpleEdge e : graph.getOutEdges(v))
				totalOut+=e.getWeight();

			debugMes("FLOW: total in for vertex "+v+" is "+totalIn + " total out  is "+totalOut+" averageCov="+v.getWeightAvg(),20);

			Collection<SimpleEdge> removeEdges = new HashSet<SimpleEdge>();
			for (SimpleEdge e : graph.getOutEdges(v))
			{
				if (e.getWeight()<totalIn* FLOW_THR || e.getWeight()<v.getWeightAvg()*FLOW_THR)
					removeEdges.add(e);
			}

			for (SimpleEdge e : graph.getInEdges(v))
			{
				if (e.getWeight()<totalOut*FLOW_THR || e.getWeight()<v.getWeightAvg()*FLOW_THR)
					removeEdges.add(e);
			}

			for (SimpleEdge re : removeEdges)
			{
				debugMes("removing low flow edge "+re+" from "+ graph.getSource(re)+" to "+graph.getDest(re),20);
				graph.removeEdge(re);
				changed = true;
			}
		}
		return changed;
	}
	/**
	 * go over the graph, and remove edges that are less than EDGE_THR (5%) from the rest of the entry flow
	 * @param graph
	 */
	private static boolean removeLightInEdges(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		debugMes("=================\nREMOVING LIGHT In EDGES\n=================",10);
		boolean somethingChanged = false;
		Queue<SeqVertex> allCurVers = new LinkedList<SeqVertex>(graph.getVertices());
		SeqVertex v = null;
		while ((v = allCurVers.poll())!=null)
		{
			if (graph.inDegree(v)<=1)
				continue;
			Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();
			int totalIn = 0;
			for (SimpleEdge inE : graph.getInEdges(v))
			{
				totalIn+=inE.getWeight();
			}

			for (SimpleEdge inE : graph.getInEdges(v))
			{
				if (inE.getWeight() <= totalIn*EDGE_THR)
				{
					debugMes("removing the edge: "+graph.getSource(inE)+"->"+graph.getDest(inE)+" ("+inE.getWeight()+" <= "+totalIn*EDGE_THR+")",20);
					removeEdges.add(inE);
					somethingChanged = true;
				}
			}
			for (SimpleEdge e : removeEdges)
				graph.removeEdge(e);
		}
		return somethingChanged;
	}


	/**
	 * go over the graph, and remove edges that are less than EDGE_THR (10%) from the rest of the exit flow
	 * @param graph
	 */
	private static boolean removeLightOutEdges(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		debugMes("=================\nREMOVING LIGHT OUT EDGES\n=================",10);
		boolean somethingChanged = false;

		Queue<SeqVertex> allCurVers = new LinkedList<SeqVertex>(graph.getVertices());
		SeqVertex v = null;
		while ((v = allCurVers.poll())!=null)
		{
			if (graph.outDegree(v)<=1)
				continue;
			Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();
			int totalOut = 0;
			for (SimpleEdge outE : graph.getOutEdges(v))
			{
				totalOut+=outE.getWeight();
			}

			for (SimpleEdge outE : graph.getOutEdges(v))
			{
				if (outE.getWeight() <= totalOut*EDGE_THR)
				{
					debugMes("removing the edge: "+graph.getSource(outE)+"->"+graph.getDest(outE),20);
					removeEdges.add(outE);
					somethingChanged = true;
				}
			}
			for (SimpleEdge e : removeEdges)
				graph.removeEdge(e);

		}
		return somethingChanged;
	}
	/*
	 * Return the SeqVertex with the given id within the given graph.
	 */
	private static SeqVertex getSeqVertex(DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			int id)
	{
		for (SeqVertex v : graph.getVertices())
		{
			if (v.getID() == id)
				return v;
		}
		return null;
	}


	/**
	 * Given the string seq, return it in fasta format
	 * @param seq - seq
	 * @param name - seq name
	 * @return
	 */
	private static String getSeqFasta(String seq,String name){
		String res = "";
		res = res.concat(">"+name+"\n");

		int i=0;
		for (; i<seq.length()-LINE_LEN ; i+=LINE_LEN)
		{
			res = res.concat(seq.substring(i, i+LINE_LEN)+"\n");
		}
		res = res.concat(seq.substring(i)+"\n");
		return res;
	}

	//	/**
	//	 * Given the graph, go over all components, and print out the best sequence in each component
	//	 * Removes components with short sequences
	//	 * @param graph
	//	 * @param file
	//	 * @throws FileNotFoundException 
	//	 */
	//	private static void PrintBestSeqs(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, String file) throws FileNotFoundException {
	//		PrintStream p = new PrintStream(new FileOutputStream(file+"_bestPaths.fasta"));
	//		PrintStream pLong = new PrintStream(new FileOutputStream(file+"_bestLongPaths.fasta"));
	//
	//		//		PrintStream pDot = new PrintStream(new FileOutputStream(file+"_finalComps.dot"));
	//
	//		String[] pathName = file.split("/");
	//		String name =pathName[pathName.length-1];
	//		Set<Set<SeqVertex>> comps = divideIntoComponents(graph);		
	//		//		pDot.println("digraph "+name+"{");
	//
	//		for (Set<SeqVertex> comp : comps)
	//		{
	//
	//			Pair<Pair<String>> res = findBestSeq(graph,comp);
	//			if (res==null)
	//				continue;
	//			String seq = res.getFirst().getFirst();
	//			String seqName = res.getFirst().getSecond(); 
	//
	//			String seqLong = res.getSecond().getFirst();
	//			String seqLongName = res.getSecond().getSecond(); 
	//			if (seq.length()>=MIN_OUTPUT_SEQ)
	//			{
	//				String curName = name+":"+seqName+"["+seq.length()+"]";
	//				p.print(getSeqFasta(seq, curName));
	//
	//			}
	//			
	//			
	//			if (seqLong.length()>=MIN_OUTPUT_SEQ)
	//			{
	//				String curName = name+":"+seqLongName+"["+seqLong.length()+"]";
	//				pLong.print(getSeqFasta(seqLong, curName));
	//			}
	////			else
	////			{
	////				//remove this component
	////				for (SeqVertex vr : comp)
	////				{
	////					graph.removeVertex(vr);
	////				}
	////			}
	//		}
	//		//		pDot.println("}");
	//
	//		p.close();
	//		pLong.close();
	//		//		pDot.close();
	//	}

	//
	//	/**
	//	 * combine suffices:
	//	 * calc for each v it's "height" in terms of length of strings (from them on)
	//	 * draw all v's with the same height
	//	 * sort on their set of children
	//	 * draw all v's with same height and same set of children
	//	 * find subsets of those with same suffix
	//	 * create new node with suffix, connect accordingly.
	//	 * add the rest (those that removed the suffix) back into queue, with new heights
	//	 */
	//	private static boolean compactSufficesBottomUp(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	//	{
	//		//		try {
	//					PrintStream p = new PrintStream(new FileOutputStream("justBefore.dot"));
	//					writeDotFile(graph,p,"tmp");
	//					p.close();
	//		
	//		//		} catch (FileNotFoundException e) {
	//		//			e.printStackTrace();
	//		//		}
	//
	//		setVerticesHeights(graph);
	//		Comparator<Object> hieghtComparator = new SeqVertexHeightComparator();
	//		PriorityQueue<SeqVertex> hQueue = new PriorityQueue<SeqVertex>(graph.getVertexCount(),hieghtComparator );
	//		for (SeqVertex v : graph.getVertices())
	//			hQueue.add(v);
	//
	//		int curH;
	//		SeqVertex v;
	//		ListComparator listComp = new ListComparator();
	//		TreeMap<List<SeqVertex>,Collection<SeqVertex>> curParents = new TreeMap<List<SeqVertex>,Collection<SeqVertex>>(listComp );
	//		boolean changed = false;
	//		for (curH=0 ; curH<=MAX_HEIGHT ; curH++) 
	//		{
	//			curParents.clear();
	//
	//			while (!hQueue.isEmpty() && hQueue.peek().getDepth()==curH)
	//			{
	//				v = hQueue.poll();
	//				//				if (v.getID()==38940 || v.getID()==38922)
	//				//					System.err.println("stop here");
	//				//hQueue.contains(getSeqVertex(graph, 96))
	//				List<SeqVertex> chil = getSortedChildrenList(graph,v); 
	//				if (!chil.isEmpty())
	//				{
	//					debugMes("curChilder: "+curParents,20);
	//					if (!curParents.containsKey(chil))
	//					{
	//						debugMes(chil +" doesn't appear in curChildren",20);
	//						curParents.put(chil,new HashSet<SeqVertex>());
	//					}
	//					debugMes("adding "+ v +" to "+curParents.get(chil),20);
	//					curParents.get(chil).add(v);
	//				}	
	//
	//			}
	//			//			curChildren.get(getSeqVertex(graph, 38840))
	//			//look for subsets with identical children
	//			for (Collection<SeqVertex> children : curParents.keySet())
	//			{
	//				// this collection has vertices with the same children
	//				Collection<SeqVertex> candidateNodes = curParents.get(children);
	//				//				if (candidateNodes==null)
	//				//					System.err.println("stop");
	//
	//				if (candidateNodes.size()==1)
	//					continue;
	//
	//				// look for shared suffix
	//				changed = compactSuffixRecursive(graph,candidateNodes);
	//
	//				boolean updateQueue = false;
	//				for (SeqVertex ver : candidateNodes)
	//				{
	//					if (ver.getName().isEmpty())
	//					{
	//						debugMes("Need to update the queue. candidateNodes = "+candidateNodes,20);
	//						updateQueue = true;
	//					}
	//				}
	//				if (updateQueue)
	//				{
	//					for (SeqVertex ver : candidateNodes)
	//					{
	//						if (!ver.getName().isEmpty())
	//						{
	//							hQueue.add(ver);
	//							debugMes("adding "+ver+" to the queue, with height "+ver.getDepth(),20);
	//						}
	//					}
	//				}
	//			}
	//
	//		}
	//		// run compactGraph after the suffices are done. 
	//		if (compactLinearPaths(graph))
	//			changed = true;
	//		return changed;
	//	}

	//	/**
	//	 * Given the graph, and the vertex v, return a sorted list of its children
	//	 * @param graph
	//	 * @param v
	//	 * @return
	//	 */
	//	private static List<SeqVertex> getSortedChildrenList(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v) {
	//
	//		List<SeqVertex> res = new ArrayList<SeqVertex>(graph.getSuccessors(v));
	//		SeqComparator verComp = new SeqComparator();
	//		Collections.sort(res, verComp);
	//
	//		return res;
	//	}
	//	/**
	//	 * Given the graph, and the candidate nodes, look for shared suffices of a single letter, 
	//	 * and move on.
	//	 * @param graph
	//	 * @param candidateNodes
	//	 */
	//	private static boolean compactSuffixRecursive(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
	//			Collection<SeqVertex> 	candidateNodes) {
	//		boolean changed = false;
	//		for (String l : LETTERS)
	//		{
	//			Collection<SeqVertex> vWithL = getVerticesWithLastLetter(candidateNodes,l);
	//			if (vWithL.size()<=1)
	//				continue;
	//			// if there is a letter that has more than one vertex, create a new vertex with this letter
	//			changed = true;
	//			SeqVertex newV = new SeqVertex(getNextID(), l,vWithL);
	//			graph.addVertex(newV);
	//			debugMes("pulled the last letter from all vertices in "+vWithL+" to the new vertex "+newV,20);
	//			Vector<SeqVertex> removeVertices = new Vector<SeqVertex>();
	//			for (SeqVertex v1 : vWithL)
	//			{
	//				Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();
	//				v1.increaseHeightByOne();
	//				// handle outgoing edges (newE2)
	//				for (SimpleEdge edgeToRemove : graph.getOutEdges(v1))
	//				{
	//
	//
	//					double w2 = edgeToRemove.getWeight();
	//					SimpleEdge newE2 = null;
	//					SeqVertex v3 = graph.getDest(edgeToRemove);
	//
	//					if (graph.findEdge(newV, v3)==null)
	//					{
	//						newE2 = new SimpleEdge(w2);
	//						graph.addEdge(newE2, newV, v3); 
	//						debugMes("adding edge "+newE2+" between "+newV+" and "+v3,20);
	//					}else
	//					{
	//						newE2 = graph.findEdge(newV, v3);
	//						if (w2>newE2.getWeight())
	//						{
	//							debugMes("setting edge "+newE2+"'s weight from "+newE2.getWeight()+" to "+w2,20); 
	//							newE2.setWeight(w2);
	//						}
	//					}
	//
	//					removeEdges.add(edgeToRemove);
	//					debugMes("removed edge "+edgeToRemove+" between "+graph.getSource(edgeToRemove)+" and "+graph.getDest(edgeToRemove),20);
	//
	//				}
	//
	//				// handle incoming edges (newE1)
	//				if (v1.getName().length()==1)
	//				{
	//					v1.removeLastLetter();
	//					//go over all edges going into v1, and move then to newV
	//					for (SeqVertex v0 : graph.getPredecessors(v1))
	//					{
	//						double w = graph.findEdge(v0, v1).getWeight();
	//						graph.addEdge(new SimpleEdge(w), v0, newV);
	//						debugMes("adding edge "+w+" between "+v0+" and "+newV,20);
	//					}
	//					removeVertices.add(v1);
	//				}else
	//				{
	//					double w = v1.removeLastLetter();
	//					SimpleEdge newE1 = new SimpleEdge(w);
	//					graph.addEdge(newE1, v1, newV);
	//					debugMes("adding edge "+newE1+" between "+v1+" and "+newV,20);
	//
	//				}
	//
	//
	//				try {
	//					PrintStream p = new PrintStream(new FileOutputStream("tmp.dot"));
	//					writeDotFile(graph,p,"tmp");
	//					p.close();
	//
	//				} catch (FileNotFoundException e) {
	//					e.printStackTrace();
	//				}
	//
	//				for (SimpleEdge re : removeEdges)
	//				{
	//					graph.removeEdge(re);
	//				}
	//			}
	//			for (SeqVertex rv : removeVertices)
	//			{
	//				graph.removeVertex(rv);
	//				debugMes("removed vertex "+rv,20);
	//			}
	//
	//			compactSuffixRecursive(graph, vWithL);
	//
	//		}
	//		return changed;
	//	}

	/**
	 * return the next available vertex id.
	 * @return
	 */
	private static int getNextID() {
		LAST_ID++;
		return LAST_ID;
	}
	//	/**
	//	 * Given the set of nodes, return a set of nodes that has the given letter l as a final letter
	//	 * @param candidateNodes
	//	 * @param l
	//	 * @return
	//	 */
	//	private static Collection<SeqVertex> getVerticesWithLastLetter(
	//			Collection<SeqVertex> candidateNodes, String l) {
	//		Collection<SeqVertex> res = new HashSet<SeqVertex>();
	//		for (SeqVertex v : candidateNodes)
	//		{
	//			if (v.getName().endsWith(l))
	//				res.add(v);
	//		}
	//		return res;
	//	}


	//	/**
	//	 * Given the graph, go over all vertices, and calculate their height, as in distance from the leaves 
	//	 * (maximal or minimal??) = doesn't matter as long as it's consistent. I chose maximal. 
	//	 * @param graph
	//	 */
	//	private static void setVerticesHeights(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
	//		List<SeqVertex> bottomUp = getTopologicalOrder(graph);
	//		Collections.reverse(bottomUp);
	//		for (SeqVertex v : bottomUp)
	//		{
	//			if (graph.outDegree(v)==0)
	//			{
	//				v.setHeight(0);
	//				//				debugMes("The height of "+v.getID()+" is set to 0");
	//			}
	//			else
	//			{
	//				int h = -1;
	//				for (SeqVertex tv : graph.getSuccessors(v))
	//				{
	//					if (tv.getDepth() + tv.getName().length() >h)
	//						h=tv.getDepth() + tv.getName().length();
	//				}
	//				v.setHeight(h);
	//				//				debugMes("The height of "+v.getID()+" is set to "+h);
	//				if (h>MAX_HEIGHT)
	//					MAX_HEIGHT = h;
	//			}
	//		}
	//
	//	}

	//	/**
	//	 * Given the graph and a single component, perform dynamic programming to 
	//	 * find the best path (which has highest average coverage 
	//	 * @param graph
	//	 * @param comp
	//	 * @return sequence of best path
	//	 */
	//	private static Pair<Pair<String>> findBestSeq(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
	//			Set<SeqVertex> comp) {
	//
	//		int n = comp.size();
	//		if (n==1)  //single node, just get the name
	//		{
	//			for (SeqVertex v : comp) 
	//			{
	//				String seq = v.getName();
	//				String name = v.getWeightAvg()+"("+v.getID()+")";
	//				if (seq==null)
	//					debugMes("Error: seq == null at node "+v.getID(),20);
	//				
	//				debugMes("best seq: "+seq,10);
	//				return new Pair<Pair<String>>(new Pair<String>(seq,name), new Pair<String>(seq,name));
	//			}
	//		}
	//
	//		// more than one vertex
	//		// create hash tables
	//		HashMap<Number,String> names = new HashMap<Number,String>();
	//		HashMap<Number,Number> lenSoFar = new HashMap<Number,Number>();
	//		HashMap<Number,Number> traceback = new HashMap<Number,Number>();
	//		HashMap<Number,Number> bestAvgSoFar = new HashMap<Number,Number>();
	//
	//		Set<SeqVertex> roots = new HashSet<SeqVertex>();
	//
	//		for (SeqVertex v : comp)
	//		{
	//			int vid = v.getID();
	//			names.put(vid, v.getName());
	//			lenSoFar.put(vid, -1);
	//			traceback.put(vid,-1);
	//			bestAvgSoFar.put(vid,-1);
	//
	//			if (graph.inDegree(v)==0)
	//			{
	//				roots.add(v);
	//				bestAvgSoFar.put(vid,v.getWeightAvg());
	//				lenSoFar.put(vid, Math.sqrt(v.getWeights().size()));
	//			} 
	//		
	//		}
	//		List<SeqVertex> topoList = getTopologicalOrder(graph,roots);
	//		int bestFinalID = -1;
	//		int bestFinalID_leaves = -1;
	//		double bestFinalAvg=-1;
	//		double bestFinalAvg_leaves=-1;
	//
	//		for (SeqVertex v : topoList)
	//		{
	//			if (graph.inDegree(v)==0)
	//				continue;
	//			int vid = v.getID();
	//
	//			int bestID = -1;
	//			double bestAvg = -1;
	//			for (SeqVertex p : graph.getPredecessors(v))
	//			{
	//				int pid = p.getID();
	//				int pnum = lenSoFar.get(pid).intValue();
	//				int vnum = v.getWeights().size();
	//				//				int thisBest =(bestAvgSoFar.get(pid).intValue()*pnum + 
	//				//				(int)graph.findEdge(p, v).getWeight() +
	//				//				v.getWeightAvg()*vnum)/(pnum+vnum+1);
	//
	//				double thisBest =(bestAvgSoFar.get(pid).intValue()*Math.sqrt(pnum) + 
	//						(int)graph.findEdge(p, v).getWeight() +
	//						v.getWeightAvg()*vnum)/Math.sqrt(pnum+vnum+1);
	//
	//
	//				if (thisBest>bestAvg)
	//				{
	//					bestAvg = thisBest;
	//					bestID = pid;
	//					lenSoFar.put(vid, pnum+vnum+1);
	//					bestAvgSoFar.put(vid, bestAvg);
	//				}
	//				
	//			}
	//			traceback.put(vid, bestID);
	//			//			if (graph.outDegree(v)==0 && bestAvg>bestFinalAvg) //doesn't have to be out degree 0...
	//			if (bestAvg>bestFinalAvg) 
	//			{
	//				debugMes("vertex "+vid+" has replaced vertex "+bestFinalID+" ( "+bestAvg+" > "+bestFinalAvg+")",10);
	//				bestFinalID = vid;
	//				bestFinalAvg = bestAvg;
	//			}
	//			
	//			if (graph.outDegree(v)==0 && bestAvg > bestFinalAvg_leaves)
	//			{
	//				debugMes("leaf vertex "+vid+" has replaced leaf vertex "+bestFinalID_leaves+" ( "+bestAvg+" > "+bestFinalAvg_leaves+")",10);
	//				bestFinalID_leaves = vid;
	//				bestFinalAvg_leaves = bestAvg;
	//
	//			}
	//
	//			
	//		}
	//		debugMes("found best path, ending at: "+bestFinalID,10);
	//		debugMes("found best path, ending at leaf: "+bestFinalID_leaves,10);
	//
	//		if (bestFinalID==-1)
	//			return null;
	//		
	//		//follow traceback:
	//		String seq = getSeqFromTraceBack(bestFinalID,names,traceback);
	//		String name = bestAvgSoFar.get(bestFinalID).intValue()+"("+bestFinalID+")";
	//
	//		Pair<String> res = new Pair<String>(seq,name);
	//		Pair<Pair<String>> bothSeqs=null;
	//		
	//		String seqLeaf=null,nameLeaf=null;
	//		if (bestFinalID_leaves!=-1)
	//		{
	//			seqLeaf = getSeqFromTraceBack(bestFinalID_leaves,names,traceback);
	//			nameLeaf = bestAvgSoFar.get(bestFinalID_leaves).intValue()+"("+bestFinalID_leaves+")";
	//			debugMes("best seqLeaf: "+seqLeaf,10);
	//			Pair<String> resLeaf = new Pair<String>(seqLeaf,nameLeaf);
	//			bothSeqs = new Pair<Pair<String>>(res, resLeaf);
	//
	//		} else
	//		{
	//			debugMes("Problem: The graph has no leaves, so seqLeaf == null at comp: "+comp,10);
	//			bothSeqs = new Pair<Pair<String>>(res, res);
	//		}
	//		
	//		if (seq==null)
	//			debugMes("Error: seq == null at comp: "+comp,10);
	//		debugMes("best seq: "+seq,10);
	//		
	//		//		String avgCov = bestAvgSoFar.get(bestFinalID).intValue()+"";
	//		
	//		return bothSeqs;
	//	}


	//	private static String getSeqFromTraceBack(int bestFinalID,
	//			HashMap<Number, String> names, HashMap<Number, Number> traceback) {
	//		String seq = "";
	//		int prevID = bestFinalID;
	//		HashMap<Number,Boolean> seen = new HashMap<Number,Boolean>(); 
	//
	//		while (prevID!=-1)  {  
	//			if (seen.containsKey(prevID)) {  
	//				debugMes("Breaking out of cycle in graph",10);  
	//				break;  
	//			}  
	//
	//			seen.put(prevID, true);  
	//			seq = names.get(prevID)+seq;
	//			prevID = traceback.get(prevID).intValue();
	//		}
	//		return seq;
	//	}


	/**
	 * return a topological order on the graph's vertices.
	 * @param graph
	 * @return list of nodes.
	 */
	private static List<SeqVertex> getTopologicalOrder(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		Set<SeqVertex> roots = new HashSet<SeqVertex>();
		for (SeqVertex v : graph.getVertices())
		{
			if (graph.inDegree(v)==0)
				roots.add(v);
		}
		return getTopologicalOrder(graph, roots);
	}

	/**
	 * return a topological order on the graph's vertices, given these roots.
	 * @param graph
	 * @param rs
	 * @return list of nodes.
	 */
	private static List<SeqVertex> getTopologicalOrder(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,Set<SeqVertex> roots){
		DFS dfs = new DFS(graph,roots);
		Visitor vst = new Visitor();
		dfs.visit(vst);
		Map<SeqVertex,Number> finished = vst.getFinishing();
		SeqVertexFinishTimeComparator finishingTimeComparator = new SeqVertexFinishTimeComparator();
		PriorityQueue<SeqVertex> fQueue = new PriorityQueue<SeqVertex>(graph.getVertexCount(),finishingTimeComparator  );

		for (SeqVertex v : finished.keySet())
		{
			fQueue.add(v);
		}
		List<SeqVertex> order = new ArrayList<SeqVertex>();
		while (!fQueue.isEmpty())
		{
			order.add(fQueue.poll());
		}
		return order;
	}



	/**
	 * Go over each sub component of the given graph, and calc the following:
	 * total coverage (sum of weights)
	 * average coverage
	 * number of paths
	 * @param graph
	 */
	private static void calcSubComponentsStats(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		//		BFSDistanceLabeler<SeqVertex, SimpleEdge> bfs = new BFSDistanceLabeler<SeqVertex, SimpleEdge>();

		Set<Set<SeqVertex>> comps = divideIntoComponents(graph);
		int numComp = comps.size();
		for (Set<SeqVertex> comp : comps)
		{
			//now we have one comp in hand
			Vector<Double> allW = new Vector<Double>();
			int compID = -1;
			for (SeqVertex v : comp)
			{
				if (compID==-1)
					compID = v.getID();

				allW.addAll(0, v.getWeights());
				for (SimpleEdge outE : graph.getOutEdges(v))
				{
					allW.add(0, outE.getWeight());
				}
			}
			SeqVertex v1 = getSeqVertex(graph, compID);
			if (allW.size()==0 || (comp.size()==1 && v1.getName().length()<MIN_OUTPUT_SEQ))
			{
				//this is a single node with a single letter
				debugMes("removing component with node "+compID+" which has only one node with short seq "+v1.getName(),20);
				graph.removeVertex(v1);
				numComp = numComp-1;
				continue;

			}

			int t=0;
			for (Double w: allW)
				t+=w;
			//			int pathNum = calcPathComplexity(graph,r);
			//			System.err.println(r.getID());
			//			if (r.getID()==3113)
			//				System.err.println(r.getID());
			int avgCov = t/allW.size();
			debugMes("SubComp: "+compID+" has "+ comp.size() +" nodes; total coverage: "+t+" average: "+avgCov,20);

			if (avgCov<COMP_AVG_COV_THR)
			{
				debugMes("removing component with node "+compID+" which has only average coverage of "+
						avgCov+ " < "+COMP_AVG_COV_THR,20);
				for (SeqVertex v : comp)
					graph.removeVertex(v);
				numComp = numComp-1;
			}
		}
		debugMes("number of good components: "+numComp,10);
	}

	private static Set<Set<SeqVertex>> divideIntoComponents(DirectedSparseGraph<SeqVertex, SimpleEdge> graph) 
	{

		WeakComponentClusterer<SeqVertex, SimpleEdge> compClus = new WeakComponentClusterer<SeqVertex, SimpleEdge>();
		Set<Set<SeqVertex>> comps = compClus.transform(graph);
		return comps;

	}

	/**
	 * connect the source node to each node with indegree=0,
	 * connect each node with outdegree=0 to the target node 
	 * Also add reads from the root to each of the nodes, and from the ends too.
	 * @param graph
	 * @param comp the current component
	 * @param combinedReadHash 
	 */
	private static void addSandT(DirectedSparseGraph<SeqVertex, SimpleEdge> graph, Set<SeqVertex> comp, HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash)
	{
		//		debugMes("=================\nADDING S AND T\n=================",10);

		graph.addVertex(ROOT);
		graph.addVertex(T_VERTEX);
		SimpleEdge e=null;
		//		for (SeqVertex v : graph.getVertices())
		for (SeqVertex v : comp)
		{
			if (graph.inDegree(v)==0 && !v.equals(ROOT) && !v.equals(T_VERTEX)) // connect S to this vertex
			{
				double w = v.getFirstWeight();
				if (w==-1) // single letter node?
				{
					debugMes("got a single letter node here.. "+v,20);
					w = 1;
				}
				e = new SimpleEdge(w);
				graph.addEdge(e, ROOT, v);

				debugMes("Adding edge from S to "+v,20);

				for (SeqVertex v2 : graph.getSuccessors(v))
				{
					PairPath pathD = new PairPath();
					pathD.addToPath1(ROOT.getID());
					pathD.addToPath1(v.getID());
					pathD.addToPath1(v2.getID());
					if (!combinedReadHash.containsKey(ROOT.getID()))
						combinedReadHash.put(ROOT.getID(), new HashMap<PairPath,Integer>());
					combinedReadHash.get(ROOT.getID()).put(pathD, MIN_TRIPLET_SUPPORT_THR);

				}				
			}

			if (graph.outDegree(v)==0 && !v.equals(T_VERTEX) && !v.equals(ROOT)) // connect this vertex to T
			{
				double w = v.getLastWeight();
				if (w==-1)
					w=1;
				e = new SimpleEdge(w);
				graph.addEdge(e, v, T_VERTEX);
				debugMes("Adding edge from "+v+" to T",20);


				for (SeqVertex v2 : graph.getPredecessors(v))
				{
					PairPath pathD = new PairPath();
					pathD.addToPath1(v2.getID());
					pathD.addToPath1(v.getID());
					pathD.addToPath1(T_VERTEX.getID());
					if (!combinedReadHash.containsKey(v2.getID()))
						combinedReadHash.put(v2.getID(), new HashMap<PairPath,Integer>());
					combinedReadHash.get(v2.getID()).put(pathD, MIN_TRIPLET_SUPPORT_THR);


				}
			}
		}
	}

	/**
	 * given the graph, remove all edges of S and T
	 * @param graph
	 */
	private static void removeAllEdgesOfSandT(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		Set<SimpleEdge> removeEdges = new HashSet<SimpleEdge>();
		if (graph.containsVertex(ROOT))
			for (SimpleEdge e : graph.getOutEdges(ROOT))
				removeEdges.add(e);
		if (graph.containsVertex(T_VERTEX))
			for (SimpleEdge e : graph.getInEdges(T_VERTEX))
				removeEdges.add(e);

		for (SimpleEdge re : removeEdges)
			graph.removeEdge(re);

	}


	/**
	 * Given the graph, and the reads, solve simple loops (self and of length 2)
	 * @param graph
	 * @param comp current component
	 * @param combinedReadHash all mapped reads
	 * @return
	 */
	private static HashMap<Integer,Integer> dealWithSimpleLoops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			Set<SeqVertex> comp, HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash, boolean justCount) {
		HashMap<Integer,Integer> res = new HashMap<Integer,Integer>(); 
		DijkstraShortestPath<SeqVertex, SimpleEdge> dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(graph);
		Set<SeqVertex> dontCheckVers = new HashSet<SeqVertex>(); 

		Set<SeqVertex> selfLoops = new HashSet<SeqVertex>();
		HashMap<SeqVertex,SeqVertex> doubleLoops = new HashMap<SeqVertex,SeqVertex>();

		Set<SeqVertex> newVers = new HashSet<SeqVertex>();

		//		for (SeqVertex v : graph.getVertices())
		for (SeqVertex v : comp)

		{

			for (SeqVertex v2 : graph.getSuccessors(v))
			{
				if (!dontCheckVers.contains(v2) && dp.getDistance(v2, v)!=null)
				{

					if (v.equals(v2)) // self loop
					{
						selfLoops.add(v);
						dontCheckVers.add(v);

					}else if (dp.getDistance(v2, v).intValue()==1) // length 2 loop
					{
						doubleLoops.put(v, v2);
						dontCheckVers.add(v);
						dontCheckVers.add(v2);

					} else // longer than 2
					{
						Integer len = dp.getDistance(v2, v).intValue()+1;
						if (!res.containsKey(len))
							res.put(len, 0);
						res.put(len,res.get(len)+1);

						List<SimpleEdge> path = dp.getPath(v2, v);
						dontCheckVers.add(v);
						for (SimpleEdge e : path)
							dontCheckVers.add(graph.getDest(e));
					}
				}
			}
		}
		if (justCount) // don't deal with the loops, just count how many
		{
			res.put(1,selfLoops.size());
			res.put(2,doubleLoops.size());
		}else 
		{
			for (SeqVertex v : selfLoops)
			{
				dealWithSelfLoops(graph,v,combinedReadHash,newVers);
			}

			for (SeqVertex v : doubleLoops.keySet())
			{
				dealWithDoubleLoops(graph,v,doubleLoops.get(v),combinedReadHash,newVers);
			}
		}


		for (SeqVertex nv : newVers)
			comp.add(nv);
		return res;
	}

	/**
	 * given the graph and the node with the self loop,
	 * find the reads that support this loop, and multiply this vertex as many times as needed, and then remap these reads.
	 * @param graph
	 * @param v
	 * @param combinedReadHash
	 * @param newVers 
	 */
	private static void dealWithSelfLoops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash, Set<SeqVertex> newVers) {


		int vid = v.getID();
		int maxNumOfOccurrences = 0;

		HashMap<PairPath, Integer> relaventReads = new HashMap<PairPath, Integer>();
		for (Integer startV : combinedReadHash.keySet())
		{
			for (PairPath path: combinedReadHash.get(startV).keySet())
			{
				int numOcc = path.numOccurrences(vid);
				if (numOcc>0)
				{					
					Integer count = combinedReadHash.get(startV).get(path);
					if (count == null)
						debugMes("stop here",10);
					relaventReads.put(path,count);
				}

				if ( numOcc> maxNumOfOccurrences) //this read includes this vertex
				{
					debugMes("the read "+path+" includes the vertex "+vid+" "+numOcc+" times",20);
					maxNumOfOccurrences = numOcc;
				}
			}
		}

		// remove the self loop
		SimpleEdge removeE = graph.findEdge(v, v);
		double oldW = removeE.getWeight();
		List<Integer> newVerIDs = new ArrayList<Integer>();
		newVerIDs.add(vid);

		graph.removeEdge(removeE);
		debugMes("removing the edge between "+ v +" and itself",20);
		// multiply this node maxNumOfOccurrences times
		int upID = vid;
		int downID = -1;

		for (int i=2; i<=maxNumOfOccurrences; i++)
		{
			if (downID!=-1)
				upID = downID;

			downID = getNextID();
			newVerIDs.add(downID);

			SeqVertex newV = new SeqVertex(downID, v);
			graph.addVertex(newV);
			SeqVertex upV = getSeqVertex(graph, upID);

			newVers.add(newV);

			for (SeqVertex vOut : graph.getSuccessors(upV))
			{
				debugMes("adding an edge between "+newV.getID()+" and "+vOut.getID(),20);
				graph.addEdge(new SimpleEdge(graph.findEdge(v, vOut)), newV, vOut);
			}
			debugMes("adding an edge between "+upID+" and "+newV.getID(),20);
			graph.addEdge(new SimpleEdge(oldW), upV, newV);

		}
		List<Integer> loopVIDs = new ArrayList<Integer>();
		loopVIDs.add(vid);
		List<List<Integer>> newVerIDsList = new ArrayList<List<Integer>>();
		newVerIDsList.add(newVerIDs);
		updateReadsAfterLoopOpening(combinedReadHash,relaventReads,loopVIDs,newVerIDsList,maxNumOfOccurrences);

	}

	/**
	 * Given the combinedReadHash, and the relevant reads, update their paths.
	 * @param combinedReadHash
	 * @param relevantReads
	 * @param loopVIDs
	 * @param newVerIDs
	 * @param maxNumOfOccurrences
	 */
	private static void updateReadsAfterLoopOpening(
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash,
			HashMap<PairPath, Integer> relevantReads, List<Integer> loopVIDs,
			List<List<Integer>> newVerIDs, int maxNumOfOccurrences) {

		for (PairPath path: relevantReads.keySet())
		{
			Integer origFirstV = path.getFirstID();
			Integer origCount = combinedReadHash.get(origFirstV).get(path);
			List<Integer> newPath1 = new ArrayList<Integer>(path.getPath1());
			List<Integer> newPath2 = new ArrayList<Integer>(path.getPath2());

			if (loopVIDs.size()==1)
			{
				updatePathOfSelfLoop(newPath1,loopVIDs,newVerIDs.get(0),maxNumOfOccurrences);
				updatePathOfSelfLoop(newPath2,loopVIDs,newVerIDs.get(0),maxNumOfOccurrences);
			} else
			{
				updatePathOfDoubleLoop(newPath1,loopVIDs,newVerIDs.get(0),newVerIDs.get(1),maxNumOfOccurrences);
				updatePathOfDoubleLoop(newPath2,loopVIDs,newVerIDs.get(0),newVerIDs.get(1),maxNumOfOccurrences);
			}
			// path hasn't changed
			if (path.getPath1().equals(newPath1) && path.getPath2().equals(newPath2))
				continue;

			// both are empty now
			if (newPath1.isEmpty() && newPath2.isEmpty())
				combinedReadHash.get(origFirstV).remove(path);

			// at least one has changed
			PairPath newKey;
			if (newPath1.isEmpty())	
				newKey = new PairPath(newPath2,new ArrayList<Integer>());
			else if (newPath2.isEmpty())
				newKey = new PairPath(newPath1,new ArrayList<Integer>());
			else
				newKey = new PairPath(newPath1,newPath2);


			Integer firstV = newKey.getFirstID();
			if (!combinedReadHash.containsKey(firstV))
				combinedReadHash.put(firstV, new HashMap<PairPath, Integer>());

			if (combinedReadHash.get(firstV).containsKey(newKey))
			{
				Integer oldCount = combinedReadHash.get(firstV).get(newKey);
				combinedReadHash.get(firstV).put(newKey,oldCount+origCount);
				combinedReadHash.get(firstV).remove(path);
			}else
			{
				combinedReadHash.get(firstV).put(newKey,origCount);
			}
		}

	}



	/**
	 * given a path, the vid of the self loop, and the new vertices' id, update the path
	 * if the path starts of ends inside the loop, trim this part of the path, and leave only the outside info.
	 * @param path
	 * @param vid
	 * @param newVerIDs
	 * @return
	 */
	private static void updatePathOfSelfLoop(List<Integer> path, List<Integer> loopVIDs,
			List<Integer> newVerIDs,int maxNumOcc) {
		int vid = loopVIDs.get(0).intValue();
		String origPath = ""+path;
		Set<Integer> loopVs = new HashSet<Integer>();
		loopVs.add(vid);
		boolean changed = false;
		if (path.contains(vid))
		{
			if (path.get(0).intValue()==vid)
			{ //starts inside the loop
				changed = true;
				if (path.get(path.size()-1).intValue()==vid)
					//starts and ends inside the loop
					if (path.size()==maxNumOcc)
					{
						for (int i=1 ; i<=path.size()-1 ; i++)
							path.set(i,newVerIDs.get(i));

						changed = true;
					}else
						path.clear();
				else
					updatePathToRemoveLoopNodes(path,loopVs);
				//			} else if (path.get(path.size()-1).intValue()==vid)
				//			{ //ends inside the loop
				//				updatePathToRemoveLoopNodes(path,loopVs);
				//				changed = true;
			}else
			{ // starts and ends outside the loop
				for (int i=1 ; i<=path.size()-1 ; i++)
				{
					if (path.get(i).intValue()==vid) // i>0
					{
						int j = newVerIDs.indexOf(path.get(i-1));
						if (j>=0)
						{
							path.set(i, newVerIDs.get(j+1));
							changed = true;
						}
					}

				}
			}
		}
		if (changed)
			debugMes("path changed from "+origPath+" to "+path,20);
	}

	/**
	 * remove the integers that are inside the loop
	 * @param path
	 * @param loopVs
	 */
	private static void updatePathToRemoveLoopNodes(List<Integer> path,
			Set<Integer> loopVs) {
		List<Integer> indicesToRemove = new ArrayList<Integer>();
		for (int i=0 ; i<=path.size()-1 ; i++)
			if (loopVs.contains(path.get(i)))
				indicesToRemove.add(i);
		Collections.sort(indicesToRemove);
		Collections.reverse(indicesToRemove);
		for (Integer i : indicesToRemove)
			path.remove(i.intValue());
	}



	/**
	 * given the graph and the node with the self loop,
	 * find the reads that support this loop, and multiply this vertex as many times as needed
	 * @param graph
	 * @param t_v1
	 * @param t_v2
	 * @param combinedReadHash
	 * @param newVers 
	 */
	private static void dealWithDoubleLoops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex t_v1,
			SeqVertex t_v2,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash, Set<SeqVertex> newVers) 
	{
		int v1_id=-1; // the one inside the regular flow
		int v2_id=-1; // the addition
		if (graph.getSuccessorCount(t_v1)>1)
		{
			v1_id = t_v1.getID();
			v2_id = t_v2.getID();
		} else if (graph.getSuccessorCount(t_v2)>1)
		{
			v1_id = t_v2.getID();
			v2_id = t_v1.getID();
		}

		if (v1_id==-1)
			return;

		HashMap<PairPath, Integer> relaventReads = new HashMap<PairPath, Integer>();

		int maxNumOfOccurrences = 0;
		for (Integer startV : combinedReadHash.keySet())
		{
			for (PairPath path: combinedReadHash.get(startV).keySet())
			{
				int numOcc2 = path.numOccurrences(v2_id);
				if (numOcc2>0)
				{
					Integer count = combinedReadHash.get(startV).get(path);
					if (count == null)
					{
						//						debugMes("stop here",10);
						count = combinedReadHash.get(startV).get(path);

						for (PairPath path2: combinedReadHash.get(startV).keySet())
						{
							debugMes("path: "+path2+" with hashCode "+path2.hashCode(),15);
							debugMes("path: "+path2+" with value "+combinedReadHash.get(startV).get(path2),15);
						}
					}
					relaventReads.put(path,count);

					if ( numOcc2> maxNumOfOccurrences) //this read includes this vertex
					{
						debugMes("MAX:the read "+path+" includes the vertex "+v1_id+" "+numOcc2+" times",15);
						maxNumOfOccurrences = numOcc2;
					}

				}
			}
		}
		// the loop is v1 (v2,v1)*
		//if we count how many v2 appears, then the number appearances of v1 is one more.
		SeqVertex v1 = getSeqVertex(graph, v1_id);
		SeqVertex v2 = getSeqVertex(graph, v2_id);
		List<Integer> newVerIDs_v1 = new ArrayList<Integer>();
		List<Integer> newVerIDs_v2 = new ArrayList<Integer>();
		newVerIDs_v1.add(v1_id);
		newVerIDs_v2.add(v2_id);

		// remove the self loop
		SimpleEdge removeE = graph.findEdge(v2, v1);
		double oldW = removeE.getWeight();
		double oldW2 = graph.findEdge(v1, v2).getWeight();

		graph.removeEdge(removeE);
		debugMes("removing the edge between "+ v2_id +" and "+v1_id,20);
		// multiply this node maxNumOfOccurrences times

		int up_v1 = v1_id;
		if (maxNumOfOccurrences>=1) //multiply only v1
		{
			SeqVertex newV = new SeqVertex(getNextID(), v1);
			graph.addVertex(newV);
			newVerIDs_v1.add(newV.getID());

			newVers.add(newV);

			for (SeqVertex vOut : graph.getSuccessors(v1))
			{
				if (!vOut.equals(v2))
				{
					debugMes("adding an edge between "+newV.getID()+" and "+vOut.getID(),20);
					graph.addEdge(new SimpleEdge(graph.findEdge(v1, vOut)), newV, vOut);
				}
			}
			debugMes("adding an edge between "+v2_id+" and "+newV.getID(),20);
			graph.addEdge(new SimpleEdge(oldW), v2, newV);

			up_v1 = newV.getID();
		}

		int up_v2 = v2_id;
		int down_v1 = -1;
		int down_v2 = -1;

		for (int i=2; i<=maxNumOfOccurrences; i++) // multiple v2-v1
		{
			if (down_v1!=-1)
			{
				up_v1 = down_v1;
			}

			down_v1 = getNextID();
			down_v2 = getNextID();
			newVerIDs_v1.add(down_v1);
			newVerIDs_v2.add(down_v2);

			SeqVertex newV1 = new SeqVertex(down_v1, v1);
			SeqVertex newV2 = new SeqVertex(down_v2, v2);

			graph.addVertex(newV1);
			graph.addVertex(newV2);

			newVers.add(newV1);
			newVers.add(newV2);

			SeqVertex upV = getSeqVertex(graph, up_v1);
			for (SeqVertex vOut : graph.getSuccessors(upV))
			{
				if (!newVerIDs_v2.contains(vOut.getID()))
				{
					debugMes("adding an edge between "+newV1.getID()+" and "+vOut.getID(),20);
					graph.addEdge(new SimpleEdge(graph.findEdge(upV, vOut)), newV1, vOut);
				}
			}

			for (SeqVertex vIn : graph.getPredecessors(getSeqVertex(graph, up_v2)))
			{
				if (!newVerIDs_v1.contains(vIn.getID()))
				{
					debugMes("adding an edge between "+vIn.getID()+" and "+down_v2,20);
					graph.addEdge(new SimpleEdge(graph.findEdge(vIn, getSeqVertex(graph, up_v2))), vIn, newV2);
				}
			}


			debugMes("adding an edge between "+up_v1+" and "+newV2.getID(),20);
			graph.addEdge(new SimpleEdge(oldW), upV, newV2);
			//			debugMes("adding an edge between "+newV1.getID()+" and "+newV2.getID(),20);
			//			graph.addEdge(new SimpleEdge(oldW2), newV1, newV2);
			debugMes("adding an edge between "+newV2.getID()+" and "+newV1.getID(),20);
			graph.addEdge(new SimpleEdge(oldW2), newV2, newV1);


		}

		List<Integer> loopVIDs = new ArrayList<Integer>();
		loopVIDs.add(v1_id);
		loopVIDs.add(v2_id);
		List<List<Integer>> newVerIDs = new ArrayList<List<Integer>>();
		newVerIDs.add(newVerIDs_v1);
		newVerIDs.add(newVerIDs_v2);

		updateReadsAfterLoopOpening(combinedReadHash,relaventReads,loopVIDs,newVerIDs,maxNumOfOccurrences);


	}




	private static void updatePathOfDoubleLoop(List<Integer> path, List<Integer> loopVIDs,
			List<Integer> newVerIDsV1, List<Integer> newVerIDsV2, int maxNumOfOccurrences) {
		int v1_id = loopVIDs.get(0).intValue();
		int v2_id = loopVIDs.get(1).intValue();

		if (path.isEmpty())
			return;
		boolean changed = false;

		String origPath = ""+path;

		Set<Integer> loopVs = new HashSet<Integer>();
		loopVs.add(v1_id);
		loopVs.add(v2_id);

		int firstV = path.get(0).intValue();
		int lastV = path.get(path.size()-1).intValue();
		if (path.contains(v2_id))
		{
			if (firstV==v1_id || firstV==v2_id)
			{
				changed = true;
				if (firstV==v1_id || lastV==v2_id)

					// the whole path is inside the loop
					if ((firstV==v1_id && lastV==v1_id && path.size()==maxNumOfOccurrences*2+1) ||
							(firstV==v2_id && lastV==v2_id && path.size()==maxNumOfOccurrences*2-1) || 
							(firstV==v1_id && lastV==v2_id && path.size()==maxNumOfOccurrences*2) || 
							(firstV==v2_id && lastV==v1_id && path.size()==maxNumOfOccurrences*2) ) // all path is in the loop, but there is only one new path that matches
					{
						changed = updateSinglePathWithDoubleLoopNodes(path,v1_id,v2_id,newVerIDsV1,newVerIDsV2);
					}else
						path.clear();
				else
				{// only the start is inside the loop
					updatePathToRemoveLoopNodes(path, loopVs); 
					changed = true;
				}
				//			}else if (lastV==v1_id || lastV==v2_id) 
				//			{ // only the end is inside the loop
				//				updatePathToRemoveLoopNodes(path, loopVs); 
				//				changed = true;
			}else
			{ // start and ends outside the loop
				changed = updateSinglePathWithDoubleLoopNodes(path,v1_id,v2_id,newVerIDsV1,newVerIDsV2);
			}
		}
		if (changed)
			debugMes("path changed from "+origPath+" to "+path,20);
	}





	/**
	 * given this path, and the loop info, update the path to its single option.
	 * @param path
	 * @param v1_id
	 * @param v2_id
	 * @param newVerIDsV1
	 * @param newVerIDsV2
	 * @return
	 */
	private static boolean updateSinglePathWithDoubleLoopNodes(
			List<Integer> path, int v1_id, int v2_id, List<Integer> newVerIDsV1,
			List<Integer> newVerIDsV2) {
		boolean changed = false;
		for (int i=1 ; i<=path.size()-1 ; i++)
		{
			if (path.get(i).intValue()==v1_id) 
			{
				int j = newVerIDsV2.indexOf(path.get(i-1));
				if (j>=0)
				{
					path.set(i, newVerIDsV1.get(j+1));
					changed = true;
				}
			} else if (path.get(i).intValue()==v2_id) 
			{
				int j = newVerIDsV1.indexOf(path.get(i-1));
				if (j>=1)
				{
					path.set(i, newVerIDsV2.get(j));
					changed = true;
				}
			}
		}
		return changed;
	}



	//	/**
	//	 * 
	//	 * @param graph
	//	 * @param file
	//	 * @throws FileNotFoundException 
	//	 */
	//	private static void printLongNodes(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, String fileName) throws FileNotFoundException {
	//		PrintStream p = new PrintStream(new FileOutputStream(fileName));
	//		String[] pathName = fileName.split("/");
	//		String name =pathName[pathName.length-1];
	//
	//		for (SeqVertex v: graph.getVertices())
	//		{
	//			if (v.getName().length()>MIN_SEQ_LENTH)
	//				p.print(getSeqFasta(v.getName(),name+"_"+v.getID()));
	//		}
	//		p.close();
	//
	//	}
	/**
	 * print out the given error message, only if DEBUG=true
	 * @param mes Message
	 */
	private static void debugMes(String mes, int verbosityLevel)
	{
		if (DEBUG && verbosityLevel<=VERBOSE_LEVEL)
			ERR_STREAM.println(mes);
	}
	/**
	 * remove the given vertex and its outgoing edges
	 * @param graph
	 * @param curVer
	 */

	/**
	 * combine suffices:
	 * calc for each v it's "depth" in terms of length of strings (from them on)
	 * draw all v's with the same height
	 * sort on their set of children
	 * draw all v's with same height and same set of children
	 * find subsets of those with same suffix
	 * create new node with suffix, connect accordingly.
	 * add the rest (those that removed the suffix) back into queue, with new heights
	 */

	private static boolean compactPrefixesBottomUp(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{

		setVerticesDepths(graph);
		Comparator<Object> depthComparator = new SeqVertexDepthComparator();
		PriorityQueue<SeqVertex> dQueue = new PriorityQueue<SeqVertex>(graph.getVertexCount(),depthComparator );
		for (SeqVertex v : graph.getVertices())
			dQueue.add(v);

		int curD;
		SeqVertex v;
		ListComparator listComp = new ListComparator();
		TreeMap<List<SeqVertex>,Collection<SeqVertex>> curParents = new TreeMap<List<SeqVertex>,Collection<SeqVertex>>(listComp );
		boolean changed = false;
		for (curD=0 ; curD<=MAX_DEPTH ; curD++) 
		{
			curParents.clear();

			while (!dQueue.isEmpty() && dQueue.peek().getDepth()==curD)
			{
				v = dQueue.poll();
				if (!graph.containsVertex(v))
					continue;

				List<SeqVertex> parents = getSortedParentList(graph,v); 
				if (!parents.isEmpty())
				{
					debugMes("curParents: "+curParents,20);
					if (!curParents.containsKey(parents))
					{
						debugMes(parents +" doesn't appear in curParents",20);
						curParents.put(parents,new HashSet<SeqVertex>());
					}
					debugMes("adding "+ v +" to "+curParents.get(parents),20);
					curParents.get(parents).add(v);
				}	

			}
			//			curChildren.get(getSeqVertex(graph, 38840))
			//look for subsets with identical children
			for (Collection<SeqVertex> parents : curParents.keySet())
			{
				// this collection has vertices with the same children
				Collection<SeqVertex> candidateNodes = curParents.get(parents);

				if (candidateNodes.size()==1)
					continue;

				// look for shared suffix
				boolean updateQueue = false;
				Collection<SeqVertex> updatedNodes = new HashSet<SeqVertex>();
				changed = compactPrefixRecursive(graph,candidateNodes,updatedNodes);

				for (SeqVertex ver : updatedNodes)
				{
					if (ver.getName().isEmpty())
					{
						debugMes("Need to update the queue. candidateNodes = "+updatedNodes,20);
						updateQueue = true;
					}
				}
				if (updateQueue)
				{
					for (SeqVertex ver : updatedNodes)
					{
						if (!ver.getName().isEmpty())
						{
							dQueue.add(ver);
							debugMes("adding "+ver+" to the queue, with depth "+ver.getDepth(),20);
						}
					}
				}
			}

		}
		// run compactGraph after the suffices are done. 
		if (compactLinearPaths(graph))
			changed = true;
		return changed;
	}

	/**
	 * Given the graph, go over all vertices, and calculate their depth, as in distance from the roots 
	 * (maximal or minimal??) = doesn't matter as long as it's consistent. I chose maximal. 
	 * @param graph
	 */
	private static void setVerticesDepths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		List<SeqVertex> topBottom = getTopologicalOrder(graph);

		for (SeqVertex v : topBottom)
		{
			if (graph.inDegree(v)==0)
			{
				v.setDepth(0);
			}
			else
			{
				int d = -1;
				for (SeqVertex tv : graph.getPredecessors(v))
				{
					if (tv.getDepth() + tv.getName().length() >d)
						d=tv.getDepth() + tv.getName().length();
				}
				v.setDepth(d);
				if (d>MAX_DEPTH)
					MAX_DEPTH = d;
			}
		}

	}

	/**
	 * Given the graph, and the vertex v, return a sorted list of its parents
	 * @param graph
	 * @param v
	 * @return
	 */
	private static List<SeqVertex> getSortedParentList(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v) {
		List<SeqVertex> res = new ArrayList<SeqVertex>(graph.getPredecessors(v));
		SeqComparator verComp = new SeqComparator();
		Collections.sort(res, verComp);

		return res;
	}

	/**
	 * Given the graph, and the candidate nodes, look for shared prefixes of a single letter, 
	 * and move on.
	 * @param graph
	 * @param candidateNodes
	 * @param updateQueue 
	 */
	private static boolean compactPrefixRecursive(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Collection<SeqVertex> 	candidateNodes, Collection<SeqVertex> updatedNodes) {
		boolean changed = false;
		//		try {
		//			PrintStream p = new PrintStream(new FileOutputStream("justBefore.dot"));
		//			writeDotFile(graph,p,"tmp");
		//			p.close();
		//
		//		} catch (FileNotFoundException e) {
		//			e.printStackTrace();
		//		}

		for (String l : LETTERS)
		{
			Collection<SeqVertex> vWithL = getVerticesWithFirstLetter(candidateNodes,l);
			if (vWithL.size()<=1)
				continue;


			// if there is a letter that has more than one vertex, create a new vertex with this letter
			changed = true;
			SeqVertex newV = new SeqVertex(getNextID(), l);
			newV.addIDsAsFirstPrevIDs(vWithL,LAST_REAL_ID);
			Collection<SeqVertex> new_vWithL = new HashSet<SeqVertex>(); 
			Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();

			for (SeqVertex v_withL : vWithL)
			{
				if (!graph.containsVertex(v_withL))
					continue;

				// create a new vertex with the first prevID as id
				SeqVertex newReplaceV_withL;
				if (!v_withL.getPrevVerIDs().isEmpty() && v_withL.getPrevVerIDs().firstElement().size()>1)
				{
					newReplaceV_withL = new SeqVertex(getNextID(), v_withL.getName());
					newReplaceV_withL.copyTheRest(v_withL);
				}else
					newReplaceV_withL = v_withL.generateNewVerWithFirstIDasID(); 
				// move all edges from and to the orig, to the new
				if (!newReplaceV_withL.equals(v_withL)) // they will be equal if the v_withL has no prevIDs, and only his original id
				{
					for (SimpleEdge e : graph.getOutEdges(v_withL))
					{
						removeEdges.add(e);
						graph.addEdge(new SimpleEdge(e.getWeight()), newReplaceV_withL, graph.getDest(e));
					}
					for (SimpleEdge e : graph.getInEdges(v_withL))
					{
						removeEdges.add(e);
						graph.addEdge(new SimpleEdge(e.getWeight()), graph.getSource(e), newReplaceV_withL);
					}
				}
				//replace it's location within vWithL
				new_vWithL.add(newReplaceV_withL);
			}

			for (SimpleEdge re : removeEdges)
			{
				debugMes("removing edge "+re+" between "+graph.getSource(re)+" and "+graph.getDest(re),20);
				graph.removeEdge(re);
			}

			for (SeqVertex rv : vWithL)
			{
				if (!new_vWithL.contains(rv))
				{
					debugMes("removing vertex (because new_vWithL doesn't contain it) "+rv,20);
					graph.removeVertex(rv);
				}
			}

			vWithL = new_vWithL;
			graph.addVertex(newV);
			debugMes("pulled the first letter from all vertices in "+vWithL+" to the new vertex "+newV,20);
			Vector<SeqVertex> removeVertices = new Vector<SeqVertex>();
			for (SeqVertex v1 : vWithL)
			{

				if (removeVertices.contains(v1) || !graph.containsVertex(v1))
					continue;
				removeEdges.clear();
				v1.increaseDepthByOne();
				// handle outgoing edges (newE2)
//				if (!graph.containsVertex(v1))
//					debugMes("v1 = "+v1,10);

				for (SimpleEdge edgeToRemove : graph.getInEdges(v1))
				{
					double w2 = edgeToRemove.getWeight();
					SimpleEdge newE2 = null;
					SeqVertex v3 = graph.getSource(edgeToRemove);

					if (graph.findEdge(v3,newV)==null)
					{
						newE2 = new SimpleEdge(w2);
						graph.addEdge(newE2, v3,newV); 
						debugMes("adding edge "+newE2+" between "+v3+" and "+newV,20);
					}else
					{
						newE2 = graph.findEdge(v3,newV);
						if (w2>newE2.getWeight())
						{
							//FIXME ?? do we want to add up the weights?
							debugMes("setting edge "+newE2+"'s weight from "+newE2.getWeight()+" to "+w2,20); 
							newE2.setWeight(w2);
						}
					}

					removeEdges.add(edgeToRemove);
					debugMes("removed edge "+edgeToRemove+" between "+graph.getSource(edgeToRemove)+" and "+graph.getDest(edgeToRemove),20);

				}
				// handle outgoing edges (newE1)
				SeqVertex newV1; // needed only if this node is less than K in length
				if  (v1.getName().length()==1)
				{
					v1.removeFirstLetter();

					//go over all edges going out of v1, and move them to exit newV
					for (SeqVertex v0 : graph.getSuccessors(v1))
					{
						double w = graph.findEdge(v1,v0).getWeight();
						graph.addEdge(new SimpleEdge(w), newV,v0);
						debugMes("adding edge "+w+" between "+newV+" and "+v0,20);
					}
					debugMes("vertex "+v1+" is going to be removed",20);
					removeVertices.add(v1);
				}else if (v1.getName().length()<=K && graph.outDegree(v1)==0) 
				{
					v1.removeFirstLetter();
					Collection<SeqVertex> upV = graph.getPredecessors(v1);

					if (v1.getID()<=LAST_REAL_ID)
					{
						newV1 = new SeqVertex(getNextID(),v1.getName());
						graph.addVertex(newV1);
						removeVertices.add(v1);
					} else
						newV1 = v1;

					//go over all edges going into v1, and move them to exit newV
					if (upV.size()==1)
					{
						for (SeqVertex upV1 : upV)
						{
							SimpleEdge oldE = graph.findEdge(upV1, v1);
							double w = oldE.getWeight();
							graph.addEdge(new SimpleEdge(w), newV,newV1);
							removeEdges.add(oldE);
							debugMes("adding edge "+w+" between "+newV+" and "+newV1,20);
							debugMes("removing edge "+w+" between "+upV1+" and "+v1,20);
							graph.addEdge(new SimpleEdge(1), v1, newV1);
						}
					}
				}else
				{
					double w = v1.removeFirstLetter();
					SimpleEdge newE1 = new SimpleEdge(w);
					graph.addEdge(newE1, newV,v1);
					debugMes("adding edge "+newE1+" between "+newV+" and "+v1,20);

				}

				for (SimpleEdge re : removeEdges)
				{
					graph.removeEdge(re);
				}
			}
			//try this out
			updatedNodes.clear();
			Set<SeqVertex> toAddTo_vWithL = new HashSet<SeqVertex>();
			int curDepth = -1;
			// use this curDepth to decide if to add the children or not.
			if (!removeVertices.isEmpty())
				for (SeqVertex ver : vWithL)
					if (!removeVertices.contains(ver))
						curDepth = ver.getDepth();
			for (SeqVertex rv : removeVertices)
			{
				for (SeqVertex vChild : graph.getSuccessors(newV))
					if (!vWithL.contains(vChild) && vChild.getDepth()==curDepth)
						toAddTo_vWithL.add(vChild);

				graph.removeVertex(rv);
				debugMes("removed vertex "+rv,20);
				if (vWithL.contains(rv))
					vWithL.remove(rv);
				if (candidateNodes.contains(rv))
					candidateNodes.remove(rv);
			}

			for (SeqVertex vToAdd : toAddTo_vWithL)
				vWithL.add(vToAdd);

			for (SeqVertex vToAdd : vWithL)
			{
				updatedNodes.add(vToAdd);
				//				if (vToAdd.getName().isEmpty() || vToAdd.getFirstWeight()==-1)
				//					for (SeqVertex vChild : graph.getSuccessors(newV))
				//						if (!vWithL.contains(vChild))
				//							toAddTo_vWithL.add(vChild);
			}


			if (vWithL.size()>1)
				compactPrefixRecursive(graph, vWithL,updatedNodes);

		}

		return changed;
	}

	/**
	 * Given the set of nodes, return a set of nodes that has the given letter l as a final letter
	 * @param candidateNodes
	 * @param l
	 * @return
	 */
	private static Collection<SeqVertex> getVerticesWithFirstLetter(
			Collection<SeqVertex> candidateNodes, String l) {
		Collection<SeqVertex> res = new HashSet<SeqVertex>();
		for (SeqVertex v : candidateNodes)
		{
			if (v.getName().startsWith(l))
				res.add(v);
		}
		return res;
	}

	

	//	/**
	//		 * Check that the last triplet of the path (w,v,u) has N supporting reads or more
	//		 * A supporting read is a read that enforces this triplet
	//		 * if w==null, check only (v,u)
	//		 * @param readsOfPathUntilV
	//		 * @param w
	//		 * @param v
	//		 * @param u
	//		 * @param dijkstraDis 
	//		 * @param graph 
	//		 * @return
	//		 */
	//		private static boolean lastTripletHasEnoughReadSupport(
	//				HashMap<PairPath, Integer> readsOfPathUntilV, SeqVertex w,
	//				SeqVertex v, SeqVertex u, DirectedSparseGraph<SeqVertex,SimpleEdge> graph, DijkstraDistanceWoVer<SeqVertex,SimpleEdge> dijkstraDisWoVer) {
	//	
	//			// if the middle node is too long, we don't expect to see connections.
	//			int vlen = v.getName().length();
	//			if (vlen > MAX_PAIR_DISTANCE)
	//				return true;
	//
	////			int ulen = u.getName().length();
	////			int wlen = 0; 
	////			if (w!=null)
	////				wlen = w.getName().length();
	////			if (vlen>MIN_SINGLE_NODE && (wlen+vlen+ulen)<MAX_PAIR_FOR_LONG_MIDDLE)
	////			{
	////				if (w==null)
	////					debugMes("the triplet (-,"+v.getID()+","+u.getID()+") is enforced due to long middle node",10);
	////				else
	////					debugMes("the triplet ("+w.getID()+","+v.getID()+","+u.getID()+") is enforced due to long middle node",10);
	////				return true;
	////			}
	//			int numberReadsSupporting = 0;
	//			for (PairPath pPath : readsOfPathUntilV.keySet())
	//			{
	//				if (readEnforcesVertex(graph, dijkstraDisWoVer, pPath, w) &&
	//						readEnforcesVertex(graph, dijkstraDisWoVer, pPath, v) &&
	//						readEnforcesVertex(graph, dijkstraDisWoVer, pPath, u))
	//				{
	//					numberReadsSupporting+=readsOfPathUntilV.get(pPath);
	//					if (w==null)
	//						debugMes("the read "+pPath+"("+readsOfPathUntilV.get(pPath)+") enforces the triplet (-,"+v.getID()+","+u.getID()+")",20);
	//					else
	//						debugMes("the read "+pPath+"("+readsOfPathUntilV.get(pPath)+") enforces the triplet ("+w.getID()+","+v.getID()+","+u.getID()+")",20);
	//				} else
	//					if (w==null)
	//						debugMes("the read "+pPath+" DOES NOT enforce the triplet (-,"+v.getID()+","+u.getID()+")",20);
	//					else
	//						debugMes("the read "+pPath+" DOES NOT enforce the triplet ("+w.getID()+","+v.getID()+","+u.getID()+")",20);
	//	
	//			}
	//	
	//			boolean res = (numberReadsSupporting>=MIN_TRIPLET_SUPPORT_THR);
	//			if (res)
	//				if (w==null)
	//					debugMes("the triplet (-,"+v.getID()+","+u.getID()+") has PASSED",20);
	//				else
	//					debugMes("the triplet ("+w.getID()+","+v.getID()+","+u.getID()+") has PASSED",20);
	//			else
	//				if (w==null)
	//					debugMes("the triplet (-,"+v.getID()+","+u.getID()+") has NOT PASSED",15);
	//				else
	//					debugMes("the triplet ("+w.getID()+","+v.getID()+","+u.getID()+") has NOT PASSED",15);
	//	
	//			return res;
	//		}


	
	//	/**
	//	 * Given the graph and the hash with all reads, find all probable paths from S to T.
	//	 * @param graph
	//	 * @param combinedReadHash
	//	 */
	//	private static void printAllProbablePathsWithReads(
	//			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
	//			HashMap<Integer,HashMap<List<List<Integer>>,Integer>> combinedReadHash) {
	//	
	//		Queue<SeqVertex> C = new LinkedList<SeqVertex>();
	//		HashMap<SeqVertex,List<List<Integer>>> Paths = new HashMap<SeqVertex,List<List<Integer>>>();
	//		HashMap<List<Integer>,List<Read>> PathReads = new HashMap<List<Integer>,List<Read>>();
	//		
	//		//initiation
	//		C.add(ROOT);
	//		List<Integer> tmpL = new ArrayList<Integer>();
	//		tmpL.add(ROOT.getID());
	//		ArrayList<List<Integer>> tmpPathList = new ArrayList<List<Integer>>();
	//		tmpPathList.add(tmpL);
	//		Paths.put(ROOT, tmpPathList);
	//		
	//		SeqVertex v;
	//		while (!C.isEmpty())
	//		{
	//			v = C.poll(); 
	//			debugMes("the next node in the queue is "+v.getID(),20);
	//			
	//			// go over all paths of P[v]
	//			for (List<Integer> path : Paths.get(v))
	//			{
	//				if (!PathReads.containsKey(path))
	//					PathReads.put(path, new ArrayList<Read>());
	//				
	//				Set<Read> reads = combinedReadHash.get(v.getID());
	//				if (reads!=null && !reads.isEmpty())
	//				{
	//					debugMes("adding the reads " +reads +" to the path "+ path, 20);
	//					PathReads.get(path).addAll(reads);
	//				}
	//			}
	//			
	//			for (SeqVertex u : graph.getSuccessors(v))
	//			{
	//				for (List<Integer> path : Paths.get(v))
	//				{
	//					// add [path,u] to paths of u
	//					if (!Paths.containsKey(u))
	//						Paths.put(u, new ArrayList<List<Integer>>());
	//	
	//					List<Integer> pathWu = new ArrayList<Integer>();
	//					pathWu.addAll(path);
	//					pathWu.add(u.getID());
	//					debugMes("adding the path " +pathWu +" to the paths of "+ u.getID()+": "+Paths.get(u), 20);
	//					Paths.get(u).add(pathWu);
	//					
	//					//update reads of [path,u]
	//					if (!PathReads.containsKey(pathWu))
	//						PathReads.put(pathWu, new ArrayList<Read>());
	//					
	//					List<Read> readsOfPath = PathReads.get(path);
	//					for (Read r : readsOfPath)
	//					{
	//						// if this read is consistent with pathWu, then add it
	//						if (r.getPathIDs().contains(u.getID()))
	//						{
	//							debugMes("read "+r+" with the path "+r.getPathIDs()+" is consistent with "+u, 20);
	//							PathReads.get(pathWu).add(r);
	//						}
	//					}
	//					
	//					if (PathReads.get(pathWu).isEmpty())
	//					{
	//						debugMes("removing path "+pathWu + " from the paths of "+u.getID(),20);
	//						Paths.get(u).remove(pathWu);
	//					}
	//	
	//				}
	//				if (!C.contains(u))
	//					C.add(u);
	//			}
	//		}
	//		
	//		debugMes("readHash = \n"+combinedReadHash,10);
	//		for (List<Integer> path : Paths.get(getSeqVertex(graph, 202)))
	//			debugMes("Final Paths: "+ path+" with "+PathReads.get(path).size() +" support: "+PathReads.get(path),20);
	//		
	//		
	//	}

}

// End TransAssembly.java
//
//		



