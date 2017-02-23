package com.aaronpmaus.jProt.executables;

import com.aaronpmaus.jProt.metrics.*;
import com.aaronpmaus.jProt.io.*;
import com.aaronpmaus.jMath.graph.*;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Date;

/**
 * <p>This program calculates the angular distance, regions of local similarity, and regions under
 * the global distance test. For the latter two tasks, it prints out the pymol scripts to color
 * those  regions of the structure. </p>
 * Pass -h as a command line argument for usage information.
 * @since 0.1.2
*/
public class AngularDistanceJProt{
    private static boolean runAngularDistance = false;
    private static boolean runLocalSimilarity = false;
    private static boolean runGDT = false;
    private static double localSimilarityThreshold = 2.0;
    private static double[] gdtThresholds;
    private static String mol1FileName;
    private static String mol2FileName;
    private static String differencesFileName;
    private static boolean diffFileProvided = false;
    private static boolean mol1FileProvided = false;
    private static boolean mol2FileProvided = false;
    private static MausMetrics theTool;
    /**
     * Runs the metrics of this program.
     * run <p>AngularDistanceJProt -h</p> for help on how to run it.
     * @param arguments the command line arguments. use -h for usage info.
    */
    public static void main(String[] arguments){
        gdtThresholds = new double[3];
        gdtThresholds[0] = 2.0;
        gdtThresholds[1] = 4.0;
        gdtThresholds[2] = 8.0;
        CommandLineParser args = new CommandLineParser(arguments);
        if(arguments.length == 0 || args.contains("-h")){
            System.out.println(" Usage: AngularDistanceJProt [<options>] <mol1-f fname> <mol2-f fname>");
            System.out.println("    options :");
            System.out.println("              -h");
            System.out.println("              \tdisplay this help file.");
            System.out.println("              -a, --angular-distance");
            System.out.println("              \tcalculate the angular distance.");
            System.out.println("              --mol1-f fname");
            System.out.println("              \tthe file containing the distance matrix for molecule one. This file is csv.");
            System.out.println("              \tThe first row contains the residue IDs. The rest of the file contains the");
            System.out.println("              \tcarbon alpha distance matrix. Each row contains that residues distances to");
            System.out.println("              \tevery other residue. Only the upper right hand side will be used.");
            System.out.println("              --mol2-f fname");
            System.out.println("              \tthe file containing the distance matrix for molecule two. The format is the");
            System.out.println("              \tsame as mol1-f.");
            System.out.println("              --diff-f fname");
            System.out.println("              \tthe file containing the differences matrix for the molecules. The format is the");
            System.out.println("              \tsame as mol1-f.");
            System.out.println("              --ls, --local-similarity");
            System.out.println("              \tfind the local similarity regions of the structures. These are the sets of");
            System.out.println("              \tresidues that are internally consistent, that is, all residues within the set");
            System.out.println("              \tare the same distance apart (within a threshold t, default 2.0 angstroms) in");
            System.out.println("              \tboth structures. This in effect find the regions of the structures that are");
            System.out.println("              \tthe same in both structures but that may oriented differently between them.");
            System.out.println("              \tThis prints to the screen a pymol script to select and color these regions of");
            System.out.println("              \tsimilarity.");
            System.out.println("              --ls_t d");
            System.out.println("              \tset the threshold used to determine the regions of similarity to d angstroms");
            System.out.println("              --gdt");
            System.out.println("              \tGlobal Distance Test: finds the largest region of similarity between the structures");
            System.out.println("              \tunder increading thresholds. It prints out a pymol script to select and color these");
            System.out.println("              \tregions and prints out the number and percent of residues within each region along");
            System.out.println("              \twith the average of the percents. By default the thresholds are {2.0, 4.0, 8.0}");
            System.out.println("              --gdt-ha");
            System.out.println("              \tGlobal Distance Test - High Accuracy: the same as gdt except with thresholds:");
            System.out.println("              \t{1.0,2.0,4.0}. Warning, this may take a LONG time to run depending on the");
            System.out.println("              \tstructures.");
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
                gdtThresholds[0] = 1.0;
                gdtThresholds[1] = 2.0;
                gdtThresholds[2] = 4.0;
            }
            if(args.contains("--ls") || args.contains("--local-similarity")){
                runLocalSimilarity = true;
            }
            if(args.contains("--ls-t")){
                localSimilarityThreshold = Double.parseDouble(args.getValue("--ls-t"));
            }
            if(args.contains("--mol1-f")){
                mol1FileName = args.getValue("--mol1-f");
                mol1FileProvided = true;
            }/* else {
                System.out.println("Both mol1-f and mol2-f are required arguments. Rerun the command passing them both in.");
                System.exit(1);
            }*/
            if(args.contains("--mol2-f")){
                mol2FileName = args.getValue("--mol2-f");
                mol2FileProvided = true;
            } /*else {
                System.out.println("Both mol1-f and mol2-f are required arguments. Rerun the command passing them both in.");
                System.exit(1);
            }*/
            if(args.contains("--diff-f")){
                differencesFileName = args.getValue("--diff-f");
                diffFileProvided = true;
            }
        }

        try {
            theTool = null;
            if(mol1FileProvided && mol2FileProvided){
                theTool = new MausMetrics(mol1FileName, mol2FileName);
            } else {
                System.out.println("You must provide the two molecule CA distance matrix files");
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
        for(int i = 0; i < localSimilarityRegions.size(); i++){
            cliquesStr += String.format("\tclique%d",i+1);
            resNumStr += String.format("\t%.0f",globalDistanceTest[i][0]);
            percentsStr += String.format("\t%.2f%%",globalDistanceTest[i][1]*100);
        }
        System.out.println(cliquesStr);
        System.out.println(resNumStr);
        System.out.println(percentsStr);
        long end = new Date().getTime();
        System.out.println("Total Time for Local Similarity Covering: " + (end - start) + " milleseconds.");
    }

    /**
     * Performs the Global Distance Test given a set of thresholds.
     * It calculates the largest region that is internally consistent for each threshold.
     * It prints out the number and percent of residues under each threshold and it prints
     * out the average of the percent of residues under all thresholds.
     * @param thresholds the thresholds to use in the global distance test.
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
