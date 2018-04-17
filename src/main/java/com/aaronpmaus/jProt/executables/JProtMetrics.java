package com.aaronpmaus.jProt.executables;

import com.aaronpmaus.jProt.protein.*;
import com.aaronpmaus.jProt.metrics.*;
import com.aaronpmaus.jProt.io.*;
import com.aaronpmaus.jMath.graph.*;

import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.File;

import java.util.Scanner;
import java.util.ArrayList;
import java.util.Date;

/**
* <p>This program can calculate the angular distance, regions of local similarity, and regions under
* the global distance test for two protein structures. For the latter two tasks, it prints out the
* pymol scripts to color those  regions of the structure. </p>
*
* <pre>
* <code>
* {@literal Usage: JProtMetrics [<options>] <mol1-f fname> <mol2-f fname>}
*  options :
*      -h
*          Display the usage file.
*      -a, --angular-distance
*          Calculate and print the angular distance.
*      --mol1-f fname
*          This file must be a PDB file. It must conform to the PDB File
*          Format Version 3.30.
*      --mol2-f fname
*          The file for molecule two. The format is the same as mol1-f.
*      --ros
*          Find the regions of similarity between the structures using  a
*          default threshold of 1.0 angstroms. These are the sets of
*          residues that are internally consistent. That is, the intra
*          structure distances between the CA's of all pairs of residues
*          this this set are the same (under a default tolerance of 1.0
*          angstroms) in both structures.
*      --ros-t d
*          Find the regions of similarity between the structures using
*          the tolerance d.
*      --diff
*          Find the regions of dissimilarity between the structures.
*          These are the sets of residues where none of their intra
*          structure distances are preserved between the structures. By
*          default, the similarity threshold is 1.0 angstroms.
*      --diff-t d
*          Find the regions of dissimilarity between the structures using
*          the threshold d.
*      --gdt
*          Global Distance Test: Iteratively find the largest regions of
*          similarity between the two structures under increasing
*          thresholds. The thresholds are {1.0, 2.0, 4.0, 8.0}.
*      --gdt-ha
*          Global Distance Test - High Accuracy: the same as gdt except
*          with thresholds: {0.5, 1.0, 2.0, 4.0}.
*      --gdt-plot
*          Global Distance Test: Operates as --gdt but with thresholds
*          {0.5, 1.0, 1.5, ..., 9.5, 10.0}.
*      --chimera
*          Print out a chimera script to color the regions of similarity
*          found by ros, gdt, gdt-ha, or diff.
*      --pymol
*          Print out a pymol script to color the regions of similarity
*          found by ros, gdt, gdt-ha, or diff.
* </code>
* </pre>
*
* @version 0.6.0
* @since 0.5.0
*/
public class JProtMetrics{
  private static boolean runAngularDistance = false;
  private static boolean runLocalSimilarity = false;
  private static boolean runGDT = false;
  private static boolean printGDTPlotData = false;
  private static boolean runGDTHA = false;
  private static boolean usePDBs = false;
  private static boolean printChimera = false;
  private static boolean printPymol = false;
  private static double localSimilarityThreshold = 1.0;
  private static double dissimilarityThreshold = 1.0;
  private static double[] gdtThresholds;
  private static String mol1FilePath;
  private static String mol2FilePath;
  private static String differencesFileName;
  private static boolean diffFileProvided = false;
  private static boolean mol1FileProvided = false;
  private static boolean mol2FileProvided = false;
  private static boolean runningGDTHA = false;
  private static boolean runningGDTPlot = false;
  private static boolean runDissimilarity = false;
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
    dissimilarityThreshold = 1.0;
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
        runGDTHA = true;
      }
      if(args.contains("--gdt-plot")){
        printGDTPlotData = true;
      }
      if(args.contains("--ros")){
        runLocalSimilarity = true;
      }
      if(args.contains("--ros-t")){
        runLocalSimilarity = true;
        localSimilarityThreshold = Double.parseDouble(args.getValue("--ros-t"));
      }
      if(args.contains("--diff")){
        runDissimilarity = true;
      }
      if(args.contains("--diff-t")){
        runDissimilarity = true;
        dissimilarityThreshold = Double.parseDouble(args.getValue("--diff-t"));
      }
      if(args.contains("--mol1-f")){
        mol1FilePath = args.getValue("--mol1-f");
        mol1FileProvided = true;
      }
      if(args.contains("--mol2-f")){
        mol2FilePath = args.getValue("--mol2-f");
        mol2FileProvided = true;
      }
      if(args.contains("--chimera")){
        printChimera = true;
      }
      if(args.contains("--pymol")){
        printPymol = true;
      }
    }

    try {
      theTool = null;
      if(mol1FileProvided && mol2FileProvided){
        File mol1File = new File(mol1FilePath);
        File mol2File = new File(mol2FilePath);
        String mol1FileName = mol1File.getName();
        String mol2FileName = mol2File.getName();

        Protein prot1 = new PDBFileIO().readInPDBFile(new FileInputStream(mol1FileName),mol1FileName);
        Protein prot2 = new PDBFileIO().readInPDBFile(new FileInputStream(mol2FileName),mol2FileName);
        theTool = new Metrics(prot1, prot2);
      } else {
        System.out.println("You must provide the two protein data files.");
        System.exit(1);
      }
      if(runAngularDistance) angularDistance();
      if(runLocalSimilarity) localSimilarity(localSimilarityThreshold);
      if(runDissimilarity)   dissimilarity(dissimilarityThreshold);
      if(runGDT) globalDistanceTest(gdtThresholds);
      if(runGDTHA) {
        gdtThresholds[0] = 0.5;
        gdtThresholds[1] = 1.0;
        gdtThresholds[2] = 2.0;
        gdtThresholds[3] = 4.0;
        runningGDTHA = true;
        globalDistanceTest(gdtThresholds);
        runningGDTHA = false;
      }
      if(printGDTPlotData){
        gdtThresholds = new double[20];
        for(int i = 0; i < 20; i++){
          gdtThresholds[i] = (i / 2.0) + 0.5;
        }
        runningGDTPlot = true;
        globalDistanceTest(gdtThresholds);
        runningGDTPlot = false;
      }
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
  private static void angularDistance(){
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
  private static void localSimilarity(double threshold){
    System.out.println("############################## Regions of Similarity ##########################");
    System.out.printf("\nUnder a threshold of: %.2f\n",threshold);
    long start = new Date().getTime();
    ArrayList<UndirectedGraph<Integer>> localSimilarityRegions;
    localSimilarityRegions = theTool.getLocalSimilarityRegions(threshold);

    // print out the pymol and/or chimera scripts if the user has specified to do so
    printScripts(localSimilarityRegions);

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
    System.out.printf("RoS-LS Score: %.2f\n", percentInTopFourCliques*100);
    long end = new Date().getTime();
    System.out.println("Total Time for RoS-LS: " + (end - start) + " milleseconds.");
  }

  private static void dissimilarity(double threshold){
    System.out.println("############################ Regions of Dissimilarity #########################");
    System.out.printf("\nUnder a threshold of: %.2f\n",threshold);
    long start = new Date().getTime();
    ArrayList<UndirectedGraph<Integer>> localSimilarityRegions;
    localSimilarityRegions = theTool.getRegionsOfDissimilarity(threshold);

    // print out the pymol and/or chimera scripts if the user has specified to do so
    printDiffScripts(localSimilarityRegions);

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
    System.out.println("The Local Dissimilarity Score is the percent of residues in the top four regions.");
    System.out.printf("RoS-DIFF Score: %.2f\n", percentInTopFourCliques*100);
    long end = new Date().getTime();
    System.out.println("Total Time for RoS-DIFF: " + (end - start) + " milleseconds.");
  }

  /**
  * Performs the Global Distance Test given a set of thresholds.
  * It calculates the largest region that is internally consistent for each threshold.
  * It prints out the number and percent of residues under each threshold and it prints
  * out the average of the percent of residues under all thresholds.
  * @param thresholds the thresholds to use in the global distance test.
  * @since 0.5.0
  */
  private static void globalDistanceTest(double[] thresholds){
    if(runningGDTHA){
      System.out.println("##################### Global Distance Test - High Accuracy ####################");
    } else if (runningGDTPlot){
      System.out.println("####################### Global Distance Test - Plot Data ######################");
    } else {
      System.out.println("############################# Global Distance Test ############################");
    }
    long start = new Date().getTime();
    ArrayList<UndirectedGraph<Integer>> globalDistanceRegions;
    int numThresholds = thresholds.length;
    globalDistanceRegions = theTool.getGlobalDistanceRegions(thresholds);

    // print out the pymol and/or chimera scripts if the user has specified to do so
    printScripts(globalDistanceRegions);

    double[][] globalDistanceTest = theTool.getGlobalDistanceTestScore(globalDistanceRegions);
    if(runningGDTHA){
      System.out.println("RoS - Global Distance Test - High Accuracy:");
    } else if(runningGDTPlot){
      System.out.println("RoS - Global Distance Test - Print Plot Data:");
    } else {
      System.out.println("RoS - Global Distance Test:");
    }
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

    if(runningGDTPlot){
      System.out.println();
      System.out.println("RoS-GDT Plot Data:");
      System.out.println("Percent Res, Threshold");
      for(int i = 0; i < thresholds.length; i++){
        System.out.printf("%.2f, %.2f\n", globalDistanceTest[i][1]*100, thresholds[i]);
      }
    } else {
      // The last row in the array holds the averages. If there are 4 thresholds,
      // rows 0-3 hold the number and percents of residues for each threshold. Row
      // 4 holds the averages.
      if(runningGDTHA){
        System.out.printf("RoS-GDT-HA Score: %.4f\n",globalDistanceTest[numThresholds][1]*100);
      } else {
        System.out.printf("RoS-GDT Score: %.4f\n",globalDistanceTest[numThresholds][1]*100);
      }
    }
    long end = new Date().getTime();
    if(runningGDTHA){
      System.out.println("Total Time for RoS-GDT-HA: " + (end - start) + " milleseconds.");
    } else if(runningGDTPlot){
      System.out.println("Total Time for RoS-GDT-Plot: " + (end - start) + " milleseconds.");
    } else {
      System.out.println("Total Time for RoS-GDT: " + (end - start) + " milleseconds.");
    }
  }

  /**
  * Private helper method to print out the pymol and/or chimera scripts if those options are
  * specified.
  */
  private static void printScripts(ArrayList<UndirectedGraph<Integer>> regionsOfSimilarity){
    ArrayList<String> script;

    if(printChimera){
      script = theTool.getChimeraColoringScript(regionsOfSimilarity);
      System.out.println("\n#Chimera Script:");
      for(String cmd: script){
        System.out.println(cmd);
      }
      System.out.println("#End of Chimera Script\n");
    }

    if(printPymol){
      script = theTool.getPymolColoringScript(regionsOfSimilarity);
      System.out.println("\n#Pymol Script:");
      for(String cmd: script){
        System.out.println(cmd);
      }
      System.out.println("#End of Pymol Script\n");
    }
  }

  /**
  * Private helper method to print out the pymol and/or chimera scripts if those options are
  * specified.
  */
  private static void printDiffScripts(ArrayList<UndirectedGraph<Integer>> regionsOfSimilarity){
    ArrayList<String> script;

    if(printChimera){
      script = theTool.getDiffChimeraColoringScript(regionsOfSimilarity);
      System.out.println("\n#Chimera Script:");
      for(String cmd: script){
        System.out.println(cmd);
      }
      System.out.println("#End of Chimera Script\n");
    }

    if(printPymol){
      script = theTool.getPymolColoringScript(regionsOfSimilarity);
      System.out.println("\n#Pymol Script:");
      for(String cmd: script){
        System.out.println(cmd);
      }
      System.out.println("#End of Pymol Script\n");
    }
  }

}
