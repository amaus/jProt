package com.aaronpmaus.jProtein.executables;

import com.aaronpmaus.jProtein.metrics.*;
import com.aaronpmaus.jMath.graph.*;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Date;

/**
 * <p>Usage: AngularDistanceJProt alphaDistancesFile betaDistancesFile differencesFile</p>
 * 
 * <p>This program calculates the angular distance, regions of local similarity, and regions under
 * the global distance test. For the latter two tasks, it prints out the pymol scripts to color
 * those  regions of the structure. </p>
 * @since 0.1.1
*/
public class AngularDistanceJProt{
    private static MausMetrics theTool;
    /**
     * <p>Usage: AngularDistanceJProt alphaDistancesFile betaDistancesFile differencesFile</p>
     *
     * Runs the metrics of this program.
     * @param args the command line arguments. See usage above.
    */
    public static void main(String[] args){
        if(args.length != 3){
            System.out.println("Usage: AngularDistanceJProt alphaDistancesFile betaDistancesFile differencesFile");
            System.exit(1);
        }

        String alphaFileName = args[0];
        String betaFileName = args[1];
        String differencesFileName = args[2];

        try {
            theTool = new MausMetrics(alphaFileName, betaFileName, differencesFileName);
            angularDistance();
            localSimilarity(2.0);
            double[] thresholds = {1.0, 2.0, 4.0, 8.0};
            globalDistanceTest(thresholds);
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
        long start = new Date().getTime();
        ArrayList<UndirectedGraph<Integer>> localSimilarityRegions;
        localSimilarityRegions = theTool.getLocalSimilarityRegions(threshold);

        ArrayList<String> pymolScript = theTool.getPymolColoringScript(localSimilarityRegions);

        System.out.println("\n#Pymol Script:");
        for(String cmd: pymolScript){
            System.out.println(cmd);
        }
        System.out.println("#End of Pymol Script\n");
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
        System.out.printf("Thresholds:\t%.2f\t%.2f\t%.2f\t%.2f\n",
                            thresholds[0], thresholds[1], thresholds[2], thresholds[3]);
        System.out.printf("Num Res:\t%.0f\t%.0f\t%.0f\t%.0f\n",
                            globalDistanceTest[0][0]*100, globalDistanceTest[1][0]*100, 
                            globalDistanceTest[2][0]*100, globalDistanceTest[3][0]*100);
        System.out.printf("Percents:\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",
                            globalDistanceTest[0][1]*100, globalDistanceTest[1][1]*100, 
                            globalDistanceTest[2][1]*100, globalDistanceTest[3][1]*100);
        System.out.printf("Score: %.4f\n",globalDistanceTest[4][1]*100);
        long end = new Date().getTime();
        System.out.println("Total Time for Global Distance Test: " + (end - start) + " milleseconds.");
    }
}
