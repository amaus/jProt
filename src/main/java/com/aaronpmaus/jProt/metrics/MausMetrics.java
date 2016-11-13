package com.aaronpmaus.jProt.metrics;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.linearAlgebra.*;
import java.util.Date;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;

/**
 * <p>MausMetrics is a collection of protein similarity metrics. It consists of Angular Distance,
 * Local Similarity, and Global Similarity.</p>
 * 
 * <p>To calculate Angular Distance, we vectorize the Carbon Alpha Distance matrices of two proteins
 * and calculate the angle between those two vectors.</p>
 *
 * <p>Local similarity takes the difference of the two distance matrices, builds a graph out of the
 * differences, finds the max clique (that is, max region where all residues are the same distance
 * apart in both structures) of that graph, removes all those nodes from the graph, finds
 * the next max clique and continues to find a clique covering of the graph. That is, a set of 
 * regions that are all internally consistent and cover the protein structure. To build the graph,
 * if two residues are the same distance apart in the two structures, we draw an edge between vertices
 * representing those residues. In practice, for the two residues, if the difference of their distance
 * from the two distance matrices is below a threshold, we say they are the same and draw
 * an edge.</p>
 *
 * <p>Global similarity relies on the same ideas for building a graph from the differences matrix and
 * finding regions that are internally consistent, but it operates as follows. It uses a series of
 * thresholds. eg. 1.0, 2.0, 4.0, 8.0. It builds a graph using each threshold and finds the max
 * clique for each graph. We get a set of regions that are subsequently larger (since we are relaxing
 * out idea of "the same"). Taking the average of the percent of residues under each threshold gives
 * us a score (credit to Zemla for the recipe for this score).</p>
 *
 * <p>A note on the file formats. All the distance matrix files use the same format.
 * format: csv. The first line is a comma separated list all all residue
 * IDs. Each line after that is a row containing that residues distances
 * to every other residue. This file should only have values in the upper right
 * hand triangle of the matrix. every other value should be a single space.
 * There should be no value for a residues distance to itself. The differences file
 * is formatted the same, except the values contains the difference of the two
 * elements from the distances matrices.</p>
 * @since 0.1.1
*/
public class MausMetrics{
    private Double[][] differencesMatrix;
    private Double[][] alphaDistancesMatrix;
    private Double[][] betaDistancesMatrix;
    private String[] residueIDs;

    /**
     * A constructor that takes the differences file. If this constructor is used,
     * the AngularDistance method can not be called because it relies on the
     * distance matrices for each structure.
     * @param differencesMatrixFileName the file containing the differences matrix for two
     *          structures.
     * @throws FileNotFoundException if the file is not found
    */
    public MausMetrics(String differencesMatrixFileName) throws FileNotFoundException{
        Scanner reader = new Scanner(new File(differencesMatrixFileName));
        this.residueIDs = reader.nextLine().split(",");
        this.differencesMatrix = readInDistanceFile(differencesMatrixFileName);
    }

    /**
     * A constructor that both the distance matrix files and the differences file.
     * @param alphaDistancesFileName the file containing the distances matrix for the first structure
     * @param betaDistancesFileName the file containing the distances matrix for the second structure
     * @param differencesMatrixFileName the file containing the differences matrix for two
     *          structures.
     * @throws FileNotFoundException if any of the files are not found
    */
    public MausMetrics(String alphaDistancesFileName, String betaDistancesFileName, String differencesMatrixFileName) throws FileNotFoundException{
        this(differencesMatrixFileName);
        this.alphaDistancesMatrix = readInDistanceFile(alphaDistancesFileName);
        this.betaDistancesMatrix = readInDistanceFile(betaDistancesFileName);
    }

    /**
     * A constructor that both the distance matrix files. It will calculate the values for the differences.
     * @param alphaDistancesFileName the file containing the distances matrix for the first structure
     * @param betaDistancesFileName the file containing the distances matrix for the second structure
     * @throws FileNotFoundException if any of the files are not found
    */
    public MausMetrics(String alphaDistancesFileName, String betaDistancesFileName) throws FileNotFoundException{
        Scanner reader = new Scanner(new File(alphaDistancesFileName));
        this.residueIDs = reader.nextLine().split(",");
        this.alphaDistancesMatrix = readInDistanceFile(alphaDistancesFileName);
        this.betaDistancesMatrix = readInDistanceFile(betaDistancesFileName);
        this.differencesMatrix = new Double[getNumResidues()][getNumResidues()];
        for(int i = 0; i < this.alphaDistancesMatrix.length; i++){
            for(int j = 0; j < this.alphaDistancesMatrix[i].length; j++){
                if(this.alphaDistancesMatrix[i][j].isNaN()){       
                    this.differencesMatrix[i][j] = Double.NaN;
                } else {
                    this.differencesMatrix[i][j] = this.alphaDistancesMatrix[i][j] - this.betaDistancesMatrix[i][j];
                }
            }
        }
    }

    /*
     * a private helper method that takes in a distances file and returns a 2D array with the
     * values. 
    */
    private Double[][] readInDistanceFile(String fileName) throws FileNotFoundException{
        Scanner fileReader = new Scanner(new File(fileName));
        fileReader.nextLine(); // throw away the first line. it's the residue IDs.

        ArrayList<String> dataLines = new ArrayList<String>();
        while(fileReader.hasNext()){
            dataLines.add(fileReader.nextLine());
        }
        int numResidues = dataLines.size();
        Double[][] data = new Double[numResidues][numResidues];
        for(int i = 0; i < dataLines.size(); i++){
            String[] tokens = dataLines.get(i).split(",");
            for(int j = 0; j < tokens.length; j++){
                if(tokens[j].equals(" ")){
                    data[i][j] = Double.NaN;
                } else {
                    data[i][j] = Double.parseDouble(tokens[j]);
                }
            }
        }
        return data;
    }

    /**
     * A query to get the differences matrix
     * @return a 2D array of Doubles containing the differences values
    */
    public Double[][] getDifferencesMatrix(){
        return this.differencesMatrix;
    }

    /**
     * A query to get the alpha distance matrix
     * @return a 2D array of Doubles containing the values
    */
    public Double[][] getAlphaDistancesMatrix(){
        return this.alphaDistancesMatrix;
    }

    /**
     * A query to get the beta distance matrix
     * @return a 2D array of Doubles containing the values
    */
    public Double[][] getBetaDistancesMatrix(){
        return this.betaDistancesMatrix;
    }

    /**
     * A query to get the residue IDs
     * @return a array of Strings containing the IDs
    */
    public String[] getResidueIDs(){
        return this.residueIDs;
    }

    /**
     * A query to get the number of residues
     * @return the number of residues
    */
    public int getNumResidues(){
        return getResidueIDs().length;
    }

    /**
     * A method to calculate the Angular Distance. Given two structure distance matrices, it
     * flattens them into a vector and calculates the angle between them. It scales the angle
     * to be in the range 0-100 where 0 is identical and returns that value.
     * @return the angular distance between the two structures. range 0-100. 0 is identical.
    */
    public double angularDistance( ) throws NullPointerException{
        if(getAlphaDistancesMatrix() == null || getBetaDistancesMatrix() == null){
            throw new NullPointerException("Can not calculate angular distance. No alpha distances or beta distances available.");
        }
        Vector alphaVec = buildVector(getAlphaDistancesMatrix());
        Vector betaVec = buildVector(getBetaDistancesMatrix());
        return alphaVec.angle(betaVec)*100/90; // adjust the range to be 0-100
    }

    /*
     * A private helper method to flatten a matrix into a vector.
    */
    private Vector buildVector(Double[][] matrix){
        int numRes = getNumResidues();
        int numValues = ((numRes - 1)*numRes)/2;
        double[] vals = new double[numValues];
        int counter = 0;
        for(int i = 0; i < matrix.length; i++){
            for(int j = i+1; j < matrix[i].length; j++){
                vals[counter] = matrix[i][j];
                counter++;
            }
        }
        return new Vector(vals);
    }

    /**
     * Returns an ArrayList of UndirectedGraphs representing the similarity regions
     * under the default threshold of 1.0 angstroms. That is, if the difference of
     * the distances between two residues is less than 1.0 angstroms, we say they are
     * "the same" distance apart in both structures. We draw an edge between them in
     * the similarity graph of the structures and use it to find regions of internal
     * consistency.
     * @return an ArrayList of the graphs of the regions of local similarity. Each graph
     *          consists of nodes representing the residues in the region of similarity.
    */
    public ArrayList<UndirectedGraph<Integer>> getLocalSimilarityRegions(){
        return getLocalSimilarityRegions(1.0);
    }

    /**
     * Returns an ArrayList of UndirectedGraphs representing the similarity regions
     * under the threshold passed in as an argument. That is, if the difference of
     * the distances between two residues is less than the threshold, we say they are
     * "the same" distance apart in both structures. We draw an edge between them in
     * the similarity graph of the structures and use it to find regions of internal
     * consistency.
     * @param threshold the threshold to use when building the similarity graph.
     * @return an ArrayList of the graphs of the regions of local similarity. Each graph
     *          consists of nodes representing the residues in the region of similarity.
    */
    public ArrayList<UndirectedGraph<Integer>> getLocalSimilarityRegions(double threshold){
        UndirectedGraph<Integer> graph = buildSimilarityGraph(threshold);
        ArrayList<UndirectedGraph<Integer>> cliques = new ArrayList<UndirectedGraph<Integer>>();
        do {
            UndirectedGraph<Integer> clique = graph.findMaxClique(graph);
            cliques.add(clique);
            for(Node<Integer> node : clique.getNodes()){
                // need to pass in a node from the original graph, not one from the clique
                // that has the same object
                graph.removeNodeFromGraph(graph.getNode(node.get()));
            }
            //System.out.println(getNodesString(clique,"+"));
        } while(graph.size() > 0);
        return cliques;
    }

    /**
     * Returns an ArrayList of the graphs of the regions of similarity under the default
     * global distance thresholds of 1.0, 2.0, 4.0, and 8.0 angstroms.
     * @return an ArrayList of the graphs of the regions of similarity under the default
     *          global distance thresholds of 1.0, 2.0, 4.0, and 8.0 angstroms.
    */
    public ArrayList<UndirectedGraph<Integer>> getGlobalDistanceRegions(){
        double[] thresholds = {1.0, 2.0, 4.0, 8.0};
        return getGlobalDistanceRegions(thresholds);
    }

    /**
     * Returns an ArrayList of the graphs of the regions of similarity under the
     * global distance thresholds passed in.
     * @param thresholds the thresholds to use to find the regions of similarity.
     * @return an ArrayList of the graphs of the regions of similarity under the
     *          global distance thresholds passed in.
    */
    public ArrayList<UndirectedGraph<Integer>> getGlobalDistanceRegions(double[] thresholds){
        ArrayList<UndirectedGraph<Integer>> regions = new ArrayList<UndirectedGraph<Integer>>();
        UndirectedGraph<Integer> lastClique = null;
        for(int i = 0; i < thresholds.length; i++){
            double threshold = thresholds[i];
            UndirectedGraph<Integer> graph = buildSimilarityGraph(threshold);
            UndirectedGraph<Integer> clique;
            if(i == 0){
                clique = graph.findMaxClique(graph);
            } else {
                graph = graph.getNeighborhood(lastClique.getNodes());
                clique = graph.findMaxClique(graph);
            }
            regions.add(clique);
            lastClique = clique;
        }
        return regions;
    }
    
    /**
     * Calculates and returns the percent of residues in each region of similarity along with
     * the average of those percents.
     * @param regions an ArrayList of the graphs of the regions of similarity. Each graph
     *                  consists of nodes residues in that region of similarity.
     * @return a 2D array. The number of rows is regions.size()+1. Each row has 2 columns.
     *          The i-th row contains at index 0 the number of residues in that region and
     *          at index 1 the percent of residues in that regions.
     *          The last row contains at index 0 the average of the number of residues in
     *          each region and at index 1 the average of the percent of residues in each region.
    */
    public double[][] getGlobalDistanceTestScore(ArrayList<UndirectedGraph<Integer>> regions){
        // returns an array of 5 numbers, the first 4 are the percent residues under the
        // thresholds. the last is the score.
        double numRes = getNumResidues();
        int numRegions = regions.size();
        double[][] ret = new double[numRegions+1][2];
        double numResAve = 0;
        double percentAve = 0;
        for(int i = 0; i < numRegions; i++){
            ret[i][0] = regions.get(i).size();
            ret[i][1] = ret[i][0] / numRes;
            numResAve += ret[i][0];
            percentAve += ret[i][1];
        }
        ret[numRegions][0] = numResAve/(double)numRegions;
        ret[numRegions][1] = percentAve/(double)numRegions;
        //ret[0] = regions.get(0).size() / numRes;
        //ret[1] = regions.get(1).size() / numRes;
        //ret[2] = regions.get(2).size() / numRes;
        //ret[3] = regions.get(3).size() / numRes;
        //ret[4] = (ret[0] + ret[1] + ret[2] + ret[3])/4;
        return ret;
    }

    /**
     * Builds and returns an ArrayList containing the pymol commands to color a protein given
     * a set of regions of similarity.
     * @param regions an ArrayList of the graphs representing the regions of similarity.
     * @return the set of pymol commands to properly color the structure.
    */
    public ArrayList<String> getPymolColoringScript(ArrayList<UndirectedGraph<Integer>> regions){
        ArrayList<String> pymolScript = new ArrayList<String>();
        int numRes = getNumResidues();
        int numRegions = regions.size();
        pymolScript.add("hide everything");
        pymolScript.add("show cartoon");
        int counter = 1;
        for(UndirectedGraph<Integer> region : regions){
            pymolScript.add("select clique"+counter+", resi " + getNodesString(region, "+"));
            counter++;
        }
        if(regions.size() > 3){
            pymolScript.add("color red");
            pymolScript.add("color orange, clique4");
            pymolScript.add("color yellow, clique3");
            pymolScript.add("color cyan, clique2");
            pymolScript.add("color green, clique1");
        } else {
            pymolScript.add("color red");
            String[] colors = {"green", "cyan", "yellow", "orange"};
            for(int i = regions.size(); i > 0; i--){
                pymolScript.add("color " + colors[i-1] + ", clique"+i);
            }
        }
        return pymolScript;
    }

    /*
     * a private helper method to build a similarity graph from the differences matrix given
     * a threshold.
     * @param threshold the threshold to use when determining whether to add an edge or not
     * @return the similarity graph built from that threshold. It consists of a node for each
     *          residue and an edge between two residues if the difference of their distance 
     *          within each of the structures is less than the threshold.
    */
    private UndirectedGraph<Integer> buildSimilarityGraph( double threshold ){
        UndirectedGraph<Integer> graph = new UndirectedGraph<Integer>();
        Double[][] data = getDifferencesMatrix();
        //String dataString = ""; // building the String for debugging purposes is VERY slow (minutes!).
        for(int row = 0; row < data.length; row++){
            for(int col = 0; col < data[row].length; col++){
                if(data[row][col].isNaN()){
                    //dataString += " N";
                } else {
                    double val = data[row][col];
                    if(val < threshold){
                        //dataString += " -";
                        //dataString += " "+val;
                        graph.addEdge(new Node<Integer>(row), new Node<Integer>(col));
                    } else if(val >= threshold){
                        //dataString += " O";
                        //dataString += " "+val;
                    }
                }
            }
            //dataString += "\n";
        } // finished reading in the data
        //System.out.print(dataString);
        return graph;
    }

    /*
     * a private helper method to build a string containing all the residue IDs of the nodes
     * in an similarity graph.
     * @param graph the similarity graph
     * @param delimiter the delimiter to use between residue IDs.
     * @return the string that contains all the residue IDs of the residues in this graph separated
     *          by the delimiter
    */
    private String getNodesString(UndirectedGraph<Integer> graph, String delimiter){
        String[] resIDs = getResidueIDs();
        Collection<Node<Integer>> nodes = graph.getNodes();
        String ret = "";
        int i = 0;
        for(Node<Integer> node : nodes){
            if(i == 0){
                i++;
                if(node.get() instanceof Integer){
                    ret += resIDs[node.get()].trim();
                }
            } else {
                if(node.get() instanceof Integer){
                    ret += delimiter+resIDs[node.get()].trim();
                }
            }
        }
        return ret;
    }
}
