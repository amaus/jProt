package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.io.*;
import com.aaronpmaus.jMath.graph.*;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.LinkedList;
import java.util.Collection;

import java.io.InputStream;

/**
* A Protein has one or more PolypeptideChains
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class Protein {
  private ArrayList<PolypeptideChain> chains;
  private String pdbName;
  private ArrayList<Bond> disulfideBonds;
  private UndirectedGraph<Atom> atoms;
  private PDBFileIO pdbIO;

  /**
  * Construct a protein from an input stream pointing to a pdb file.
  * @param pdbFileInputStream an InputStream pointing to a pdb file
  * @param pdbName the PDB filename base.
  */
  public Protein(InputStream pdbFileInputStream, String pdbName){
    this.chains = new ArrayList<PolypeptideChain>();
    this.atoms = new UndirectedGraph<Atom>();
    this.disulfideBonds = new ArrayList<Bond>();
    this.pdbName = pdbName;
    this.pdbIO = new PDBFileIO(pdbFileInputStream);
    buildProteinFromPDB();
  }

  /**
  * Return the name of the PDB file for this Protein
  * @return the name of the PDB file without the .pdb extension.
  */
  public String getPDBName(){
    return this.pdbName;
  }

  /**
  * Add a chain to this protein
  * @param chain the PolypeptideChain to add to this protein
  */
  public void addChain(PolypeptideChain chain){
    this.chains.add(chain);
    for(Bond bond : chain.getBonds()){
      atoms.addEdge(bond.getAtomOne(), bond.getAtomTwo());
    }
  }

  /**
  * Create a disulfide bond in this protein.
  * <p>
  * Atoms must both be sulfur atoms in cysteine and must not be the same atom.
  * @param chainID1 the chain that the first sulfur is in
  * @param resID1 the residue id of the residue that the first sulfur is in
  * @param chainID2 the chain that the second sulfur is in
  * @param resID2 the residue id of the residue that the second sulfur is in
  * @throws IllegalArgumentException if the cysteine sulfurs are not found or are the same atom
  */
  public void addDisulfideBond(String chainID1, int resID1, String chainID2, int resID2){
    //System.out.println(getChain(chainID1).getResidue(resID1));
    //System.out.println(getChain(chainID2).getResidue(resID2));
    Atom atomOne = getChain(chainID1).getResidue(resID1).getAtom("SG");
    Atom atomTwo = getChain(chainID2).getResidue(resID2).getAtom("SG");
    if(atomOne.equals(atomTwo)){
      throw new IllegalArgumentException("Atoms in disulfide bond must be different atoms.");
    }
    Bond disulfideBond = new Bond(atomOne, atomTwo);
    disulfideBonds.add(disulfideBond);
    atoms.addEdge(atomOne, atomTwo);
    // if the disulfide bond is within the same chain, add it to that chain.
    if(chainID1.equals(chainID2)){
      getChain(chainID1).addBond(disulfideBond);
    }
  }

  /**
  * Returns a reference to a chain in the protein
  * @param chainID the ID of the chain to get
  * @return the PolypeptideChain with the ID provided as an argument
  */
  public PolypeptideChain getChain(String chainID){
    for(PolypeptideChain chain : chains){
      if(chain.getChainID().equals(chainID)){
        return chain;
      }
    }
    throw new IllegalArgumentException("There is no chain "
    +chainID+ " in the protein.");
  }

  /**
  * Returns the number of Chains in this Protein
  * @return the number of chains in this protein
  */
  public int getNumChains(){
    return this.chains.size();
  }

  /**
  * Return the number of bonds between the two atoms via the shortest path from a
  * starting atom to an ending atom.
  *
  * @param atomStart one of the end points in the path
  * @param atomEnd the other end point in the path
  * @return the number of bonds between them
  */
  public int getBondSeparation(Atom atomStart, Atom atomEnd){
    int pathSize = this.getShortestPath(atomStart, atomEnd).size();
    if(pathSize == 0){
      return 0;
    }
    return pathSize-1;
  }

  /**
  * Return the shortest path from one atom to another in this protein.
  *
  * @param atomStart the atom to calculate the path from
  * @param atomEnd the atom to calculate the path to.
  * @return the collection of atoms starting with atomStart and ending with atomEnd or an empty
  * collection if no path exists.
  */
  public List<Atom> getShortestPath(Atom atomStart, Atom atomEnd){
    System.out.printf("%s\n%s",atomStart, atomEnd);
    System.out.println(atoms.contains(atomStart));
    System.out.println(atoms.contains(atomEnd));
    LinkedList<Atom> list = new LinkedList<Atom>();
    for(Node<Atom> n : this.atoms.shortestPath(atomStart, atomEnd)){
      list.add(n.get());
    }
    return list;
  }

  /**
  * Returns the number of Residues in this Protein
  * @return the number of Residues in this protein
  */
  public int getNumResidues(){
    int numResidues = 0;
    for(PolypeptideChain chain : this.chains){
      numResidues += chain.getNumResidues();
    }
    return numResidues;
  }

  /**
  * Returns the number of Atoms in this Protein
  * @return the number of Atoms in this Protein
  */
  public int getNumAtoms(){
    int numAtoms = 0;
    for(PolypeptideChain chain : this.chains){
      numAtoms += chain.getNumAtoms();
    }
    return numAtoms;
  }

  /**
  * Return the sequence of this Protein in the form of a String of the single letter
  * residue names.
  *
  * @return a String containing the single letter residue names with no spaces
  */
  public String getSequence(){
    String seq = "";
    for(PolypeptideChain chain : this.chains){
      seq += chain.getSequence();
    }
    return seq;
  }

  /**
  * Calculate and return the CA Distance Matrix of this protein
  *
  * @return a 2D array of Double containing the CA distances
  */
  public Double[][] calculateCarbonAlphaDistanceMatrix(){
    boolean[] mask = new boolean[getNumResidues()];
    Arrays.fill(mask, true);
    return calculateCarbonAlphaDistanceMatrix(mask);
  }

  /**
  * Calculate and return the CA Distance Matrix of this protein using only the residues
  * specified by the mask.
  *
  * The mask must contain as many elements as there are residues in this protein. Each
  * residue in order will have a value true or false. If true, include that residue
  * in the distance matrix, otherwise, ignore it.
  *
  * @param mask an array of boolean containing as many values as there are residues, each
  *  indicating whether to include that residue in the distance matrix.
  * @return a 2D array of Double containing the CA distances
  * @since 0.6.0
  */
  public Double[][] calculateCarbonAlphaDistanceMatrix(boolean[] mask){
    ArrayList<Residue> residues = new ArrayList<Residue>();
    int maskIndex = 0;
    for(PolypeptideChain chain : this.chains){
      for(Residue res : chain){
        if(mask[maskIndex]){
          residues.add(res);
        }
        maskIndex++;
      }
    }
    int numResidues = residues.size();
    Double[][] distanceMatrix = new Double[numResidues][numResidues];
    int i = 0;
    int j = 0;
    for(Residue residueOne : residues){
      for(Residue residueTwo : residues){
        Atom carbonAlphaOne = residueOne.getAtom("CA");
        Atom carbonAlphaTwo = residueTwo.getAtom("CA");
        distanceMatrix[i][j] = carbonAlphaOne.distance(carbonAlphaTwo);
        j++;
      }
      j = 0;
      i++;
    }
    return distanceMatrix;
  }

  /**
  * Calculates and returns the CA Distance Matrix of one of the chains in
  * this protein
  *
  * @param chainID the ID of the chain to get the CA Distance Matrix of
  * @return a 2D array of Double containing the CA distances
  */
  public Double[][] calculateCarbonAlphaDistanceMatrix(String chainID){
    PolypeptideChain chain = getChain(chainID);
    return chain.calculateCarbonAlphaDistanceMatrix();
  }

  /**
  * Return the residue IDs of all the residues in this Protein
  *
  * @return an array of Integers holding the residue IDs
  */
  public Integer[] getResidueIDs(){
    boolean[] mask = new boolean[getNumResidues()];
    Arrays.fill(mask, true);
    return getResidueIDs(mask);
  }

  /**
  * Return the residue IDs of all the residues in this Protein specified by the mask.
  *
  * The mask must contain as many elements as there are residues in this protein. Each
  * residue in order will have a value true or false. If true, include that residue
  * ID in the array, otherwise, ignore it.
  *
  * @param mask an array of boolean containing as many values as there are residues, each
  *  indicating whether to include that residue in the IDs list.
  * @return an array of Integers holding the residue IDs
  */
  public Integer[] getResidueIDs(boolean[] mask){
    int numResidues = 0;
    for(boolean val : mask){
      if(val){
        numResidues++;
      }
    }
    Integer[] residueIDs = new Integer[numResidues];
    int residueArrayIndex = 0;
    int maskIndex = 0;
    for(PolypeptideChain chain : this.chains){
      for(Integer id : getResidueIDs(chain.getChainID())){
        if(mask[maskIndex]){
          residueIDs[residueArrayIndex] = id;
          residueArrayIndex++;
        }
        maskIndex++;
      }
    }
    return residueIDs;
  }

  /**
  * Returns the residue IDs of all the residues in the chain with the
  * chainID provided as an argument
  * @param chainID the letter ID of the chain to get the residue IDs from
  * @return an array of Integers holding the residue IDs
  */
  public Integer[] getResidueIDs(String chainID){
    return getChain(chainID).getResidueIDs();
  }

  public boolean contains(Atom atom){
    return this.atoms.contains(atom);
  }

  //*************** File IO Helper Methods *******************//
  private void buildProteinFromPDB(){
    Collection<String> chainIDs = this.pdbIO.getListOfChainIDs();
    for(String chainID : chainIDs){

      PolypeptideChain chain = new PolypeptideChain(chainID);
      // get a list of all the atomRecords that have chainID
      Collection<AtomRecord> chainAtomRecords = this.pdbIO.getAtomRecordsInChain(chainID);

      HashMap<Integer, ArrayList<AtomRecord>> residueRecordsLists =
          this.pdbIO.getResidueRecordsLists(chainAtomRecords);

      // for every residue's list of AtomRecords, build a list of the atoms in that residue from the
      // AtomRecords. Use that list to construct and add that new residue with those atoms to the
      // chain.
      for(ArrayList<AtomRecord> residueAtomRecords : residueRecordsLists.values()){
        String resName = residueAtomRecords.get(0).getResName();
        int resSeq = residueAtomRecords.get(0).getResSeq();

        Collection<Atom> residueAtoms = constructAtomsInResidue(residueAtomRecords);

        chain.addResidue(new Residue(resName, resSeq, residueAtoms));
      }
      this.addChain(chain);
    }

    addDisulfideBonds();
  }

  private void addDisulfideBonds(){
    for(SSBondRecord ssBondRec: this.pdbIO.getSSBondRecords()){
      String chainID1 = ssBondRec.getChainID1();
      String chainID2 = ssBondRec.getChainID2();
      int resID1 = ssBondRec.getResID1();
      int resID2 = ssBondRec.getResID2();
      this.addDisulfideBond(chainID1, resID1, chainID2, resID2);
    }
  }


  /**
  * From a Collection of Atom Records, build and return a collection of those Atoms.
  */
  private Collection<Atom> constructAtomsInResidue(Collection<AtomRecord> residueAtomRecords){
    ArrayList<Atom> residueAtoms = new ArrayList<Atom>();
    for(AtomRecord rec : residueAtomRecords){
      residueAtoms.add(constructAtom(rec));
    }
    return residueAtoms;
  }

  private Atom constructAtom(AtomRecord rec){
    return new Atom(rec.getAtomName(), rec.getSerial(), rec.getOccupancy(),
        rec.getTempFactor(), rec.getCharge(),
        rec.getX(), rec.getY(), rec.getZ());
  }
}
