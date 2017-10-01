package com.aaronpmaus.jProt.energy;

import com.aaronpmaus.jProt.io.*;

import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Collection;

import java.io.FileInputStream;
import java.io.File;
import java.io.InputStream;
import java.io.PrintWriter;

/**
* A PotentialOfMeanForce is a set of energy functions, one for each possible atom type pair. Each
* energy function represents the non bonded energy for that pair.
* <p>
* To read in a PMF from file, use PMFFileIo:
* <p>
* {@code InputStream in = new FileInputStream(new File(pathToFile));}
* <br>
* {@code PotentialOfMeanForce pmf = PMFFileIO.readInPMF(in);}
* <p>
*
* @version 0.6.0
* @since 0.6.0
*/
public class PotentialOfMeanForce {
  private String encadFileName;
  private int fullResidueSkip;
  private double minDistance;
  private double maxDistance;
  private double binWidth;
  private int numBins;
  private HashMap<String, List<Double>> energies;
  private HashMap<String, Integer> atomTypeIndices;
  private List<Double> binCenters;
  private HashMap<Integer, String> atomTypeNames;

  /**
  * Construct a Potential of Mean Force with all the information it requires.
  *
  * @param encadFileName the name of the encad force field file that corresponds to this PMF
  * @param fullResidueSkip the number of residues skipped when counting interactions for this PMF
  * @param minDistance the closest distance any contact can be counted at
  * @param maxDistance the farthest distance any contact can be counted at
  * @param binWidth the width of the bins, a physical distances.
  * @param numBins the number of bins in this PMF.
  * @param energies a map containing an energy function for each pair of atoms. The keys are atom
  * type pairs formatted as: "AN-AN", "AN-"VO", etc. The values are the set of energy values
  * for that atom type pair ordered by the bin that energy is for. The zeroth value is the energy
  * of the closest bin.
  * @param atomTypeIndices a map where the key is an atom type name and the value is its index in
  * the PMF.
  * @param binCenters a list of the bin centers.
  */
  public PotentialOfMeanForce(String encadFileName, int fullResidueSkip,
                              double minDistance, double maxDistance,
                              double binWidth, int numBins,
                              HashMap<String, List<Double>> energies,
                              HashMap<String, Integer> atomTypeIndices,
                              List<Double> binCenters){

    this.encadFileName = encadFileName;
    this.fullResidueSkip = fullResidueSkip;
    this.minDistance = minDistance;
    this.maxDistance = maxDistance;
    this.binWidth = binWidth;
    this.numBins = numBins;
    this.energies = energies;
    this.atomTypeIndices = atomTypeIndices;
    this.binCenters = binCenters;
    this.atomTypeNames = buildAtomTypeNames();
  }

  /**
  * @return the encad force field file for this PFM.
  */
  public String getEncadFileName(){
    return this.encadFileName;
  }

  /**
  * The full residue skip is the number of residues skipped when counting up the interactions. For
  * example, a value of 1 indicates that only nonbonded interactions outside of an Atom's
  * residue and immediately neighboring residues are counted in this PMF
  * @return the number of neighboring residues skipped when counting nonbonded interactions.
  */
  public int getFullResidueSkip(){
    return this.fullResidueSkip;
  }

  /**
  * @return the closest distance any pair interaction can be counted at
  */
  public double getMinDistance(){
    return this.minDistance;
  }

  /**
  * @return the farthest distance any pair interaction can be counted at
  */
  public double getMaxDistance(){
    return this.maxDistance;
  }

  /**
  * Return the width of the bins in this PMF. A single bin is a range of distances at which to
  * count interactions for a given atom pair. The bin centers and bin width completely define
  * these ranges.
  * @return the bin width for all bins in this PMF
  */
  public double getBinWidth(){
    return this.binWidth;
  }

  /**
  * @return the number of bins for each energy function in this PMF
  */
  public int getNumBins(){
    return this.numBins;
  }

  /**
  * @return a List of the bin centers in order from smallest to largest
  */
  public List<Double> getBinCenters(){
    return new ArrayList<Double>(this.binCenters);
  }

  /**
  * Return a list of the energy values for the given pair.
  * <p>
  * This list is the same size as
  * {@code getBinCenters()} and each value represents the energy at the corresponding bin
  * as defined in the Bin Centers List. For example, if the 21 element in the bin centers list
  * is 2.1, and the bin width is 0.1, then the 21st element in this energy function is the
  * energy for an interaction of these two atoms in the distance range 2.05 - 2.15.
  * @param atomOne one of these atoms in this pair
  * @param atomTwo the other atom in this pair
  * @return the energy list of values of the energy function for this pair of atoms
  */
  public List<Double> getEnergyFunction(String atomOne, String atomTwo){
    if(!containsEnergyFunctionFor(atomOne, atomTwo)){
      throw new IllegalStateException("PMFFileIO: no energy function for atom pair: " +
          atomOne + " - " + atomTwo);
    } else {
      String key = buildPairString(atomOne, atomTwo);
      return this.energies.get(key);
    }
  }

  /**
  * Set an energy function for two atoms to hold the energies passed in. If there is already
  * an energy function for this atom pair, it is replaced by the one passed in.
  * @param atomOne one of the atoms of the pair
  * @param atomTwo the other atom of the pair
  * @param energies the energy function to set for this atom pair
  */
  public void setEnergyFunction(String atomOne, String atomTwo, List<Double> energies){
    String key = buildPairString(atomOne, atomTwo);
    this.energies.put(key, energies);
  }

  public void write(PrintWriter fileWriter){
    validateEnergies();
    PMFFileIO.writePMF(fileWriter, this);
  }

  private HashMap<Integer, String> buildAtomTypeNames(){
    HashMap<Integer, String> atomTypeNames = new HashMap<Integer, String>();
    Collection<String> names = atomTypeIndices.keySet();
    for(String name : names){
      //System.out.println("looking for index of " + name);
      atomTypeNames.put(getAtomTypeIndex(name), name);
    }
    return atomTypeNames;
  }

  /**
  * Return the PMF atom type index for a given atom type.
  * @param atomName one of the atom types of the PMF
  * @return the index of that atom type as listed in the PMF
  */
  public int getAtomTypeIndex(String atomName){
    if(this.atomTypeIndices.containsKey(atomName)){
      return this.atomTypeIndices.get(atomName);
    }
    throw new IllegalArgumentException(atomName + " is not a valid atom type for this PMF");
  }

  /**
  * Return the PMF atom type name for a given atom type index.
  * @param atomIndex a valid index into from the PMF
  * @return the name of that atom type as listed in the PMF
  */
  private String getAtomTypeName(int atomIndex){
    if(this.atomTypeNames.containsKey(atomIndex)){
      return this.atomTypeNames.get(atomIndex);
    }
    throw new IllegalArgumentException(atomIndex + " is not a valid atom type index for this PMF");
  }

  /**
  * Ensure that there are energies for every atom type pair, and ensure that every
  * set of energies has enough values.
  */
  private void validateEnergies(){
    List<String[]> pairs = getAtomPairs();

    for(String[] pair : pairs){
      String atomOne = pair[0];
      String atomTwo = pair[1];
      List<Double> energyValues = getEnergyFunction(atomOne, atomTwo);
      if(energyValues.size() != getNumBins()){
        throw new IllegalStateException("PMFFileIO: " + pair + " has " + energyValues.size()
            + "energy values, requires " + getNumBins() + " energy values.");
      }
    }
  }

  /**
  * Check if this PMF contains an energy function for the two atom types
  * @param atomOne one of the atoms in the pair
  * @param atomTwo the other atom in this pair
  * @return true if this PMF contains an energy function for this atom pair
  */
  private boolean containsEnergyFunctionFor(String atomOne, String atomTwo){
    int atomOneIndex = getAtomTypeIndex(atomOne);
    int atomTwoIndex = getAtomTypeIndex(atomTwo);
    String pair = buildPairString(atomOne, atomTwo);
    return energies.containsKey(pair);
  }

  /**
  * @return a list of the atom pairs in this PMF. Each pair is a String[] pair of size 2 containing
  * the two atom types in this pair.
  */
  public List<String[]> getAtomPairs( ){
    List<Integer> indices = new ArrayList<Integer>(this.atomTypeNames.keySet());
    Collections.sort(indices);
    // get all atom names in ascending order by their index (as they would be listed in a PMF)
    ArrayList<String> names = new ArrayList<String>(indices.size());
    for(Integer index : indices){
      names.add(getAtomTypeName(index));
    }

    // Build a list of the pair strings "AN-AN", "AN-ACA",... "GO-GO"
    ArrayList<String[]> pairs = new ArrayList<String[]>();
    for(int i = 0; i < names.size(); i++){
      for(int j = i; j < names.size(); j++){
        String[] pair = new String[2];
        pair[0] = names.get(i);
        pair[1] = names.get(j);
        pairs.add(pair);
      }
    }
    return pairs;
  }

  private String buildPairString(String atomOne, String atomTwo){
    String pair = "";
    int atomOneIndex = getAtomTypeIndex(atomOne);
    int atomTwoIndex = getAtomTypeIndex(atomTwo);
    if(atomOneIndex < atomTwoIndex){
      pair = atomOne+"-"+atomTwo;
    } else {
      pair = atomTwo+"-"+atomOne;
    }
    return pair;
  }
}
