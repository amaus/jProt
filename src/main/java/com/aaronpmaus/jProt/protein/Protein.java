package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.io.*;
import com.aaronpmaus.jProt.sequence.*;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.transformations.Transformable;
import com.aaronpmaus.jMath.transformations.Transformation;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.LinkedList;
import java.util.Collection;
import java.util.Iterator;

import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.FileWriter;
import java.io.IOException;


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
* @version 0.6.0
* @since 0.6.0
*/
public class Protein implements Iterable<PolypeptideChain>, Transformable{
  private ArrayList<PolypeptideChain> chains;
  private String proteinName;
  private ArrayList<Bond> disulfideBondsBetweenChains;
  private ArrayList<Bond> disulfideBondsWithinChains;
  private UndirectedGraph<Atom> atoms;
  private PDBFileIO pdbIO;

  /**
  * Create an instance of a Protein. All chains will then have to be added to this instance.
  * @param proteinName the name of the protein, recommended: the PDB filename base (the part before
  * the extension).
  * @param pdbIO a PDBFileIO instance that this protein can use to write out to file
  */
  public Protein(String proteinName, PDBFileIO pdbIO){
    this.chains = new ArrayList<PolypeptideChain>();
    this.disulfideBondsBetweenChains = new ArrayList<Bond>();
    this.disulfideBondsWithinChains = new ArrayList<Bond>();
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
    // if the disulfide bond is within the same chain, add it to that chain.
    if(chainID1.equals(chainID2)){
      getChain(chainID1).addBond(disulfideBond);
      this.disulfideBondsWithinChains.add(disulfideBond);
    } else {
      this.disulfideBondsBetweenChains.add(disulfideBond);
    }
  }

  /**
  * Return all the disulfide bonds in this Protein.
  * @return a {@code List<Bond>} of all the disulfide bonds.
  * @since 0.7.0
  */
  public List<Bond> getDisulfideBonds(){
    ArrayList<Bond> allDisulfideBonds = new ArrayList<Bond>(this.disulfideBondsWithinChains);
    allDisulfideBonds.addAll(this.disulfideBondsBetweenChains);
    return allDisulfideBonds;
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
    throw new IllegalArgumentException("There is no chain " + chainID + " in the protein.");
  }

  /**
  * Get the Chain ID of the chain that contains a Residue.
  * @param residue a Residue in this Protein
  * @return the chain ID of the chain that contains residue
  * @throws IllegalStateException if no chain contains the residue
  * @since 0.7.0
  */
  public String getChainID(Residue residue){
    String chainID = null;
    for(PolypeptideChain chain : this){
      if(chain.contains(residue)){
        chainID = chain.getChainID();
      }
    }
    if(chainID == null){
      throw new IllegalStateException(
          String.format("Residue: %s Not contained in any chain.",residue));
    }
    return chainID;
  }

  /**
  * Get the Chain ID of the chain that contains an Atom.
  * @param atom an Atom in this Protein
  * @return the chain ID of the chain that contains atom
  * @throws IllegalStateException if no chain contains the residue
  * @since 0.7.0
  */
  public String getChainID(Atom atom){
    String chainID = null;
    for(PolypeptideChain chain : this){
      if(chain.contains(atom)){
        chainID = chain.getChainID();
      }
    }
    if(chainID == null){
      throw new IllegalStateException(
          String.format("Atom: %s Not contained in any chain.",atom));
    }
    return chainID;
  }

  /**
  * Returns the number of Chains in this Protein
  * @return the number of chains in this protein
  */
  public int getNumChains(){
    return this.chains.size();
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
  * @return the number of covalent bonds in this Protein
  */
  public int getNumBonds(){
    int numBonds = 0;
    for(PolypeptideChain chain : this){
      numBonds += chain.getNumBonds();
    }
    numBonds += this.disulfideBondsBetweenChains.size();
    return numBonds;
  }

  /**
  * @return a Collection of all the covalent bonds in this Protein
  */
  public List<Bond> getBonds(){
    ArrayList<Bond> bonds = new ArrayList<Bond>();
    for(PolypeptideChain chain : this){
      bonds.addAll(chain.getBonds());
    }
    bonds.addAll(this.disulfideBondsBetweenChains);
    return bonds;
  }

  /**
  * Return the sequence of this Protein in the form of a String of the single letter
  * residue names.
  *
  * @return a String containing the single letter residue names with no spaces
  */
  public ProteinSequence getSequence(){
    String seq = "";
    for(PolypeptideChain chain : this.chains){
      seq += chain.getSequence().getSequenceString();
    }
    return new ProteinSequence(seq);
  }

  /**
  * Return the residue IDs of all the residues in this Protein
  *
  * @return an array of Integers holding the residue IDs in order
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
    for(PolypeptideChain chain : this){
      if(chain.contains(atom)){
        return true;
      }
    }
    return false;
  }

  /**
  * @param chain the PolypeptideChain to check for
  * @return true if this protein contains the chain
  */
  public boolean contains(PolypeptideChain chain){
    return this.chains.contains(chain);
  }

  public boolean contains(String chainID){
    for(PolypeptideChain chain: this){
      if(chain.getChainID().equals(chainID)){
        return true;
      }
    }
    return false;
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
  * Write this protein out to a PDB file.
  * @param fileName the name of the file to write out to, will create new file if doesn't exist,
  * will override if does exist.
  * @throws IOException if the output file can not be opened
  */
  public void writeToFile(String fileName) throws IOException{
    OutputStreamWriter outputStream = new FileWriter(fileName);
    this.pdbIO.writeToPDB(outputStream, this);
    outputStream.close();
  }

  @Override
  public void applyTransformation(Transformation t){
    for(PolypeptideChain chain : this ){
      chain.applyTransformation(t);
    }
  }

  @Override
  public Iterator<PolypeptideChain> iterator(){
    return this.chains.iterator();
  }

}
