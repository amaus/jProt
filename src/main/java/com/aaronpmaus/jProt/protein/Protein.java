package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.io.*;
import com.aaronpmaus.jMath.graph.*;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.LinkedList;
import java.util.Collection;
import java.util.Iterator;

import java.io.InputStream;

/**
* A Protein is composed of PolypeptideChains which are composed of Residues. Proteins can either be
* constructed from PDB files or they can be constructed residue by residue ab initio style. <p>
*
* To read in a Protein from a PDB file:
* <p>
* {@code InputStream in = new FileInputStream(new File(pathToPDBFile));}
* <br>
* {@code Protein prot = PDBFileIO.readInPDBFile(in, proteinName)}
* <p>
*
* A protein by default has hydrogen atoms disabled. They can be enabled via the {@code
* enableHydrogens()} method. When {@code enableHydrogens()} is called, all hydrogens read in from
* the pdb structure will be added to the protein. TODO: Automatically build all hydrogens when when
* a protein is constructed so that when they are enabled ({@code enableHydrogens()}), all hydrogens
* in the structure are added to the protein. Hydrogens can also be disabled with the method {@code
* disableHydrogens()}.
* <p>
* @version 0.6.0
* @since 0.6.0
*/
public class Protein implements Iterable<PolypeptideChain>{
  private ArrayList<PolypeptideChain> chains;
  private String proteinName;
  private ArrayList<Bond> disulfideBonds;
  private UndirectedGraph<Atom> atoms;
  private PDBFileIO pdbIO;
  private boolean hydrogensEnabled = false;

  /**
  * Create an instance of a Protein. All chains will then have to be added to this instance.
  * @param proteinName the name of the protein, recommended: the PDB filename base (the part before
  * the extension).
  * @param pdbIO a PDBFileIO instance that this protein can use to write out to file
  */
  public Protein(String proteinName, PDBFileIO pdbIO){
    this.chains = new ArrayList<PolypeptideChain>();
    this.atoms = new UndirectedGraph<Atom>();
    this.disulfideBonds = new ArrayList<Bond>();
    this.proteinName = proteinName;
    this.pdbIO = pdbIO;
  }

  /**
  * Create an instance of a Protein. All chains will then have to be added to this instance.
  * @param proteinName the name of the protein, recommended: the PDB filename base (the part before
  * the extension).
  */
  public Protein(String proteinName){
    this(proteinName, new PDBFileIO());
  }

  /**
  * Return the name of the PDB file for this Protein
  * @return the name of the PDB file without the .pdb extension.
  */
  public String getProteinName(){
    return this.proteinName;
  }

  /**
  * Add a chain to this protein
  * @param chain the PolypeptideChain to add to this protein
  */
  public void addChain(PolypeptideChain chain){
    this.chains.add(chain);
    addToProteinGraph(chain);
  }

  private void addToProteinGraph(PolypeptideChain chain){
    Collection<Bond> bonds = chain.getBonds();
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
  private List<Atom> getShortestPath(Atom atomStart, Atom atomEnd){
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

  public int getNumBonds(){
    return this.atoms.numEdges();
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
      for(Integer id : chain.getResidueIDs()){
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
  * @param atom the atom to check if it is in this protein
  * @return true if this protein contains the atom
  */
  private boolean contains(Atom atom){
    return this.atoms.contains(atom);
  }

  public void enableHydrogens(){
    // set the hydrogensEnabled flag in Protein (instance variable)
    this.hydrogensEnabled = true;
    // tell every chain to enable hydrogens, then rebuild the atoms graph
    // getting all the bonds from each chain
    this.atoms = new UndirectedGraph<Atom>();
    for(PolypeptideChain chain : this){
      chain.enableHydrogens();
      addToProteinGraph(chain);
    }
    // Add all edges for inter-chain disulfide bonds.
    for(Bond bond : this.disulfideBonds){
      String atomOneChain = getAtomChainID(bond.getAtomOne());
      String atomTwoChain = getAtomChainID(bond.getAtomTwo());
      if(!atomOneChain.equals(atomTwoChain)){
        atoms.addEdge(bond.getAtomOne(), bond.getAtomTwo());
      }
    }
  }

  public void disableHydrogens(){
    // set the hydrogensEnabled flag in Protein (instance variable)
    this.hydrogensEnabled = false;
    // tell every chain to enable hydrogens, then rebuild the atoms graph
    // getting all the bonds from each chain
    this.atoms = new UndirectedGraph<Atom>();
    for(PolypeptideChain chain : this){
      chain.disableHydrogens();
      addToProteinGraph(chain);
    }
  }

  private String getAtomChainID(Atom atom){
    String chainID;
    for(PolypeptideChain chain : this){
      if(chain.containsAtom(atom)){
        return chain.getChainID();
      }
    }
    throw new IllegalStateException(String.format("Atom \n%sNot in Protein",atom));
  }

  /**
  * @return true if hydrogens are enabled for this protein
  */
  public boolean hydrogensEnabled(){
    return this.hydrogensEnabled;
  }

  @Override
  public Iterator<PolypeptideChain> iterator(){
    return this.chains.iterator();
  }

}
