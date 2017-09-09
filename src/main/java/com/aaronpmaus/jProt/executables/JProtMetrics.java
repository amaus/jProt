package com.aaronpmaus.jProt.executables;

import com.aaronpmaus.jProt.protein.*;
import com.aaronpmaus.jProt.metrics.*;
import com.aaronpmaus.jProt.io.*;
import com.aaronpmaus.jMath.graph.*;

import java.io.FileNotFoundException;
import java.io.InputStream;

import java.util.Scanner;
import java.util.ArrayList;
import java.util.Date;

/**
* <p>This program calculates the angular distance, regions of local similarity, and regions under
* the global distance test. For the latter two tasks, it prints out the pymol scripts to color
* those  regions of the structure. </p>
*
* Pass -h as a command line argument for usage information.
*
* @version 0.6.0
* @since 0.5.0
*/
public class JProtMetrics{
  private static boolean runAngularDistance = false;
  private static boolean runLocalSimilarity = false;
  private static boolean runGDT = false;
  private static boolean usePDBs = false;
  private static boolean useCSVs = false;
  private static double localSimilarityThreshold = 1.0;
  private static double[] gdtThresholds;
  private static String mol1FileName;
  private static String mol2FileName;
  private static String differencesFileName;
  private static boolean diffFileProvided = false;
  private static boolean mol1FileProvided = false;
  private static boolean mol2FileProvided = false;
  private static Metrics theTool;
  /**
  * Runs the metrics of this program.
  * run <p>AngularDistanceJProt -h</p> for help on how to run it.
  * @param arguments the command line arguments. use -h for usage info.
  */
  public static void main(String[] arguments){
    gdtThresholds = new double[4];
    gdtThresholds[0] = 1.0;
    gdtThresholds[1] = 2.0;
    gdtThresholds[2] = 4.0;
    gdtThresholds[3] = 8.0;
    localSimilarityThreshold = 1.0;
    CommandLineParser args = new CommandLineParser(arguments);
    if(arguments.length == 0 || args.contains("-h")){
      InputStream stream = JProtMetrics.class.getResourceAsStream("JProtMetricsUsage.txt");
      Scanner in = new Scanner(stream);
      while(in.hasNextLine()){
        System.out.println(in.nextLine());
      }
      System.exit(1);
    } else {
      if(args.contains("-a") || args.contains("--angular-distance")){ //set set runAngularDistance flag to true
        runAngularDistance = true;
      }
      if(args.contains("--gdt")){
        runGDT = true;
      }
      if(args.contains("--gdt-ha")){
        runGDT = true;
        gdtThresholds[0] = 0.5;
        gdtThresholds[1] = 1.0;
        gdtThresholds[2] = 2.0;
        gdtThresholds[3] = 4.0;
      }
      if(args.contains("--ls") || args.contains("--local-similarity")){
        runLocalSimilarity = true;
      }
      if(args.contains("--ls-t")){
        runLocalSimilarity = true;
        localSimilarityThreshold = Double.parseDouble(args.getValue("--ls-t"));
      }
      if(args.contains("--mol1-f")){
        mol1FileName = args.getValue("--mol1-f");
        mol1FileProvided = true;
      }
      if(args.contains("--mol2-f")){
        mol2FileName = args.getValue("--mol2-f");
        mol2FileProvided = true;
      }
    }

    try {
      theTool = null;
      if(mol1FileProvided && mol2FileProvided){
        String[] mol1FileParts = mol1FileName.split("\\.");
        String[] mol2FileParts = mol2FileName.split("\\.");
        String mol1Ext = mol1FileParts[mol1FileParts.length - 1].trim().toLowerCase();
        String mol2Ext = mol2FileParts[mol2FileParts.length - 1].trim().toLowerCase();
        
        if(mol1Ext.equals("pdb") && mol2Ext.equals("pdb")){
          usePDBs = true;
        } else if(mol1Ext.equals("csv") && mol2Ext.equals("csv")){
          useCSVs = true;
        } else {
          System.out.println("Both files must either be PDBs with extension .pdb"
            + " or CSVs with extension .csv");
          System.out.println("You provided files: \n" + mol1FileName + "\n" + mol2FileName);
          System.exit(1);
        }

        if(usePDBs){
          PDBFileIO pdb = new PDBFileIO();
          Protein prot1 = pdb.readInPDBFile(mol1FileName);
          Protein prot2 = pdb.readInPDBFile(mol2FileName);
          theTool = new Metrics(prot1, prot2);
        } else if(useCSVs){
          theTool = new Metrics(mol1FileName, mol2FileName);
        }
      } else {
        System.out.println("You must provide the two protein data files.");
        System.exit(1);
      }
      if(runAngularDistance) angularDistance();
      if(runLocalSimilarity) localSimilarity(localSimilarityThreshold);
      if(runGDT) globalDistanceTest(gdtThresholds);
    } catch (FileNotFoundException e){
      System.out.println("Could not open required files. Check for existence.");
      System.exit(1);
    }

  }

  /**
  * Calculates and prints the Angular Distance between the two structures.
  * Range: (0,100)
  * 0 is identical.
  * @since 0.5.0
  */
  public static void angularDistance(){
    System.out.println("############################### AngularDistance ###############################");
    System.out.println();
    System.out.printf("AngularDistance: %.4f\n",theTool.angularDistance());
    System.out.println();
  }

  /**
  * Finds the local similarity covering of the graph. It finds a set of regions that are all
  * internally consistent and cover the graph. This can be used to find the domains of a structure
  * that don't changes internally but do change in orientation of position relative to each other.
  * @param threshold the threshold to use for finding the regions of local similarity
  * @since 0.5.0
  */
  public static void localSimilarity(double threshold){
    System.out.println("########################### Local Similarity Covering #########################");
    System.out.printf("\nUnder a threshold of: %.2f\n",threshold);
    long start = new Date().getTime();
    ArrayList<UndirectedGraph<Integer>> localSimilarityRegions;
    localSimilarityRegions = theTool.getLocalSimilarityRegions(threshold);

    ArrayList<String> pymolScript = theTool.getPymolColoringScript(localSimilarityRegions);

    System.out.println("\n#Pymol Script:");
    for(String cmd: pymolScript){
      System.out.println(cmd);
    }
    System.out.println("#End of Pymol Script\n");

    double[][] globalDistanceTest = theTool.getGlobalDistanceTestScore(localSimilarityRegions);
    String cliquesStr = "Cliques:";
    String resNumStr = "Num Res:";
    String percentsStr = "Percents:";
    double percentInTopFourCliques = 0;
    for(int i = 0; i < localSimilarityRegions.size(); i++){
      cliquesStr += String.format("\tclique%d",i+1);
      resNumStr += String.format("\t%.0f",globalDistanceTest[i][0]);
      percentsStr += String.format("\t%.2f%%",globalDistanceTest[i][1]*100);
      if(i < 4){
        percentInTopFourCliques += globalDistanceTest[i][1];
      }
    }
    System.out.println(cliquesStr);
    System.out.println(resNumStr);
    System.out.println(percentsStr);
    System.out.println("The Local Similarity Score is the percent of residues in the top four regions.");
    System.out.printf("LS Score: %.2f\n", percentInTopFourCliques*100);
    long end = new Date().getTime();
    System.out.println("Total Time for Local Similarity Covering: " + (end - start) + " milleseconds.");
  }

  /**
  * Performs the Global Distance Test given a set of thresholds.
  * It calculates the largest region that is internally consistent for each threshold.
  * It prints out the number and percent of residues under each threshold and it prints
  * out the average of the percent of residues under all thresholds.
  * @param thresholds the thresholds to use in the global distance test.
  * @since 0.5.0
  */
  public static void globalDistanceTest(double[] thresholds){
    System.out.println("############################# Global Distance Test ############################");
    long start = new Date().getTime();
    ArrayList<UndirectedGraph<Integer>> globalDistanceRegions;
    int numThresholds = thresholds.length;
    globalDistanceRegions = theTool.getGlobalDistanceRegions(thresholds);

    ArrayList<String> pymolScript = theTool.getPymolColoringScript(globalDistanceRegions);

    System.out.println("\n#Pymol Script:");
    for(String cmd: pymolScript){
      System.out.println(cmd);
    }
    System.out.println("#End of Pymol Script\n");

    double[][] globalDistanceTest = theTool.getGlobalDistanceTestScore(globalDistanceRegions);
    System.out.println("Global Distance Test:");
    System.out.println("Total num residues: " + theTool.getNumResidues());
    String thresholdsStr = "Thresholds:";
    String resNumStr = "Num Res:";
    String percentsStr = "Percents:";
    for(int i = 0; i < thresholds.length; i++){
      thresholdsStr += String.format("\t%.2f",thresholds[i]);
      resNumStr += String.format("\t%.0f",globalDistanceTest[i][0]);
      percentsStr += String.format("\t%.2f%%",globalDistanceTest[i][1]*100);
    }
    System.out.println(thresholdsStr);
    System.out.println(resNumStr);
    System.out.println(percentsStr);
    //System.out.printf("Thresholds:\t%.2f\t%.2f\t%.2f\t%.2f\n",
    //                    thresholds[0], thresholds[1], thresholds[2], thresholds[3]);
    //System.out.printf("Num Res:\t%.0f\t%.0f\t%.0f\t%.0f\n",
    //                    globalDistanceTest[0][0], globalDistanceTest[1][0],
    //                    globalDistanceTest[2][0], globalDistanceTest[3][0]);
    //System.out.printf("Percents:\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",
    //                    globalDistanceTest[0][1]*100, globalDistanceTest[1][1]*100,
    //                    globalDistanceTest[2][1]*100, globalDistanceTest[3][1]*100);

    // The last row in the array holds the averages. If there are 4 thresholds,
    // rows 0-3 hold the number and percents of residues for each threshold. Row
    // 4 holds the averages.
    System.out.printf("Score: %.4f\n",globalDistanceTest[numThresholds][1]*100);
    long end = new Date().getTime();
    System.out.println("Total Time for Global Distance Test: " + (end - start) + " milleseconds.");
  }
}
