package com.aaronpmaus.jProt.io;

import com.aaronpmaus.jProt.energy.*;

import java.util.Collections;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import java.io.FileInputStream;
import java.io.File;
import java.io.InputStream;
import java.io.PrintWriter;

import java.math.BigDecimal;

/**
* PMFFileIO provides an easy to use interface to PMF Files. This class allows for simple and
* efficient access into any of the energy functions listed in a pmf.
* <p>
* To read in a PMF:
* <p>
* {@code InputStream in = new FileInputStream(new File(pathToFile));}
* <br>
* {@code PotentialOfMeanForce pmf = PMFFileIO.readInPMF(in);}
* <p>
* The values in this class are immutable. It is intended to be used as the tool to read in
* or write out a PMF.
* <p>
* To read in a PMF,
* @version 0.6.0
* @since 0.6.0
*/
public class PMFFileIO{

  /**
  * Read in a Potential of Mean Force from an InputStream.
  * @param inputStream an input stream that points to an open PMF file
  * @return a PotentialOfMeanForce containing the inform held in that file
  */
  public static PotentialOfMeanForce readInPMF(InputStream inputStream){
    HashMap<String, EnergyRecord> energyRecords = new HashMap<String, EnergyRecord>();
    HashMap<String, Integer> atomTypeIndices = new HashMap<String, Integer>();
    String encadFileName = null;
    int fullResidueSkip = -1;
    double minDistance = -1;
    double maxDistance = -1;
    BinRecord binRecord = null;
    Scanner in = new Scanner(inputStream);
    while(in.hasNextLine()){
      String line = in.nextLine();
      String[] tokens = line.trim().split("\\s+");
      switch(tokens[0].trim()){
        case "ENCADFILENAME":
          encadFileName = in.nextLine();
          break;
        case "FULLRESIDUESKIP":
          fullResidueSkip = Integer.parseInt(tokens[1]);
          break;
        case "SHORTESTDISTANCE":
          minDistance = Double.parseDouble(tokens[1]);
          break;
        case "LONGESTDISTANCE":
          maxDistance = Double.parseDouble(tokens[1]);
          break;
        case "PMFSTART":
          line = in.nextLine();
          binRecord = parseBinsLine(line);
          break;
        case "%RS":
          EnergyRecord record = parseEnergyRecord(line, binRecord.getNumBins());
          energyRecords.put(record.getAtomOne()+"-"+record.getAtomTwo(), record);
          break;
      }
    }

    for(EnergyRecord record: energyRecords.values()){
      atomTypeIndices.put(record.getAtomOne(), record.getAtomOneIndex());
    }

    return new PotentialOfMeanForce(encadFileName, fullResidueSkip, minDistance, maxDistance,
                                    binRecord.getBinWidth(), binRecord.getNumBins(),
                                    getEnergies(energyRecords), atomTypeIndices,
                                    binRecord.getBinCenters());
  }

  /**
  * Write a PMF out to the given PrintWriter. The format is the standard pmf format as used by
  * encad.
  * @param out the PrintWriter to write the the pmf out to
  * @param pmf the PotentialOfMeanForce to write to file
  */
  public static void writePMF(PrintWriter out, PotentialOfMeanForce pmf) {
    List<String[]> pairs = pmf.getAtomPairs();
    out.printf("ENCADFILENAME\n");
    out.printf("%s\n", pmf.getEncadFileName());
    out.printf("FULLRESIDUESKIP\t%d\n", pmf.getFullResidueSkip());
    out.printf("SHORTESTDISTANCE\t%.2f\n", pmf.getMinDistance());
    out.printf("LONGESTDISTANCE\t%.2f\n", pmf.getMaxDistance());
    out.printf("DISTANCE\n");
    double low =  pmf.getMinDistance();
    double high = low +  pmf.getBinWidth();
    for(int i = 0; i <  pmf.getNumBins(); i++){
      out.printf("%10s\t%s\n",getNumInSigFigs(low, 6), getNumInSigFigs(high, 6));
      low +=  pmf.getBinWidth();
      high +=  pmf.getBinWidth();
    }
    out.printf("PMFSTART\n");
    // print out the bins line
    out.printf("%s\n",new BinRecord(pmf.getBinWidth(), pmf.getNumBins(), pmf.getBinCenters()));
    for(String[] pair : pairs){
      String atomOne = pair[0];
      String atomTwo = pair[1];
      int atomOneIndex = pmf.getAtomTypeIndex(atomOne);
      int atomTwoIndex = pmf.getAtomTypeIndex(atomTwo);
      out.printf("%s\n", new EnergyRecord(atomOneIndex, atomTwoIndex, atomOne, atomTwo,
                                          pmf.getEnergyFunction(atomOne, atomTwo)));
    }
  }

  private static String getNumInSigFigs(double num, int significantFigures){
    BigDecimal bd = BigDecimal.valueOf(num);
    return String.format("%."+significantFigures+"G",bd);
  }

  /**
  * Build a map of energy functions from a map of EnergyRecords. The key is an atom type pair, the
  * atom types concatenated together with a `-` between them. For example: "AN-VO" for "CSG-CSG".
  * Each atom-type pair is only contained once. That is, if there is a key for "AN-VO", then there
  * is no key for "VO-AN".
  * @return a HashMap where the key are the atom type pairs and the values are ArrayLists of the
  * values for energy function for that atom type pair.
  */
  private static HashMap<String, List<Double>> getEnergies(HashMap<String, EnergyRecord> records){
    HashMap<String, List<Double>> energies = new HashMap<String, List<Double>>();
    for(String key : records.keySet()){
      energies.put(key, records.get(key).getEnergies());
    }
    return energies;
  }

  private static BinRecord parseBinsLine(String line){
    //example beginning of line:
    //%RS del_r=0.1 n_r=201R=0.000   0.100   0.200   ...
    String[] tokens = line.trim().split("\\s+");
    //get the del_r=0.1 or bin size value
    double binWidth = Double.parseDouble((tokens[1].split("="))[1].trim());
    System.out.println("bin size: " + binWidth);
    //work with the n_r=201R=0.000 part
    String[] mess = tokens[2].split("R");
    //mess[0] is n_r=201
    //mess[1] is =0.000
    //get n_r=201 or numBins value:
    int numBins = Integer.parseInt(mess[0].split("=")[1]);
    System.out.println("num bins: " + numBins);
    ArrayList<Double> binCenters = new ArrayList<Double>(numBins);
    double firstBin = Double.parseDouble(mess[1].split("=")[1]);
    System.out.println("first bin: " + firstBin);
    binCenters.add(firstBin);

    for(int i = 3; i < tokens.length; i++){
      binCenters.add(Double.parseDouble(tokens[i]));
    }
    return new BinRecord(binWidth, numBins, binCenters);
  }

  private static EnergyRecord parseEnergyRecord(String line, int numBins){
    String[] tokens = line.split("\\s+");
    int atomOneIndex = Integer.parseInt(tokens[1]);
    int atomTwoIndex = Integer.parseInt(tokens[2]);
    String atomOne = tokens[3];
    String atomTwo = tokens[4];

    ArrayList<Double> energies = new ArrayList<Double>(numBins);
    for(int i = 5; i < tokens.length; i++){
      energies.add(Double.parseDouble(tokens[i]));
    }

    return new EnergyRecord(atomOneIndex, atomTwoIndex, atomOne, atomTwo, energies);
  }


  private static class EnergyRecord{
    private int atomOneIndex;
    private int atomTwoIndex;
    private String atomOne;
    private String atomTwo;
    private List<Double> energies;

    public EnergyRecord(int atomOneIndex, int atomTwoIndex,
                        String atomOne, String atomTwo,
                        List<Double> energies){

      this.atomOneIndex = atomOneIndex;
      this.atomTwoIndex = atomTwoIndex;
      this.atomOne = atomOne;
      this.atomTwo = atomTwo;
      this.energies = energies;
    }

    public int getAtomOneIndex(){ return this.atomOneIndex; }
    public int getAtomTwoIndex(){ return this.atomTwoIndex; }
    public String getAtomOne(){ return this.atomOne; }
    public String getAtomTwo(){ return this.atomTwo; }
    public List<Double> getEnergies(){ return new ArrayList<Double>(this.energies); }

    @Override
    public String toString(){
      String str = String.format("%%RS %3d %3d %-4s %-4s  ",
          this.getAtomOneIndex(), this.getAtomTwoIndex(), this.getAtomOne(), this.getAtomTwo());
      for(Double energy : energies){
        str += String.format("%6.3f  ",energy);
      }
      return str;
    }
  }

  private static class BinRecord{
    private double binWidth;
    private int numBins;
    private List<Double> binCenters;

    public BinRecord(double binWidth, int numBins, List<Double> binCenters){
      this.binWidth = binWidth;
      this.numBins = numBins;
      this.binCenters = binCenters;
    }
    public double getBinWidth(){ return this.binWidth; }
    public int getNumBins(){ return this.numBins; }
    public List<Double> getBinCenters(){ return this.binCenters; }

    @Override
    public String toString(){
      String str = String.format("%%RS del_r=%2.1f n_r=%dR=",getBinWidth(), getNumBins());
      for(Double center : getBinCenters()){
        str += String.format("%.3f   ",center);
      }
      return str;
    }
  }
}
