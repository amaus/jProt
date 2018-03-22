package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.sequence.*;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.linearAlgebra.Vector3D;
import com.aaronpmaus.jMath.transformations.Transformable;
import com.aaronpmaus.jMath.transformations.Transformation;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;

/**
* A PolypeptideChain is a single chain of amino acids.
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class PolypeptideChain extends Molecule implements Iterable<Residue>, Transformable{
  // ordered by residueIDs
  private ArrayList<Residue> residues;
  // key: residueID, value: residue index in residues list
  private HashMap<Integer, Integer> residueIndices;
  // key: Bond, value: residueID
  private HashMap<Bond, Integer> bondResidueIDs;
  // An efficient way of looking up which residue an atom is in
  private HashMap<Atom, Residue> residuesByAtom;
  private ArrayList<Bond> rotatableBonds;
  private String chainID;

  /**
  * A constructor for a PolypeptideChain.
  * @param chainID the letter ID of this chain
  */
  public PolypeptideChain(String chainID){
    super();
    this.residues = new ArrayList<Residue>();
    this.residueIndices = new HashMap<Integer,Integer>();
    this.bondResidueIDs = new HashMap<Bond, Integer>();
    this.rotatableBonds = new ArrayList<Bond>();
    this.residuesByAtom = new HashMap<Atom, Residue>();
    this.chainID = chainID;
  }

  /**
  * Returns the Chain ID of this chain
  * @return the String that is the Chain ID: A, B, C, etc...
  */
  public String getChainID(){
    return this.chainID;
  }

  /**
  * Add a residue to this chain
  * @param residue the residue to add
  */
  public void addResidue(Residue residue){
    residues.add(residue);
    Collections.sort(residues, new ResidueComparator());
    buildResidueIndices();
    addBondsFromResidue(residue);
    for(Atom atom : residue){
      this.residuesByAtom.put(atom, residue);
    }
  }

  private void buildResidueIndices(){
    residueIndices = new HashMap<Integer,Integer>();
    int index = 0;
    for(Residue res : this){
      residueIndices.put(res.getResidueID(), index);
      index++;
    }
  }

  /**
  * Get and add all Bonds from the Residue to this Chain.
  * @param residue the residue whose bonds need to be added to this chain.
  */
  private void addBondsFromResidue(Residue residue){
    Collection<Bond> bonds = residue.getBonds();
    for(Bond b : bonds){
      this.addBond(b); // inherited from Molecule
      if(b.containsAtom("N") && b.containsAtom("CA")){
        this.bondResidueIDs.put(b, residue.getResidueID());
        this.rotatableBonds.add(b);
      }
      if(b.containsAtom("CA") && b.containsAtom("C")){
        this.bondResidueIDs.put(b, residue.getResidueID());
        this.rotatableBonds.add(b);
      }
    }
    // if the previous residue - the residue with the residueID immediately
    // preceeding this residue's ID - exists, then add a bond between it's
    // C and this residue's N.
    if(this.contains(residue.getResidueID()-1)){
      Residue prevResidue = getResidue(residue.getResidueID()-1);
      Bond b = new Bond(prevResidue.getAtom("C"), residue.getAtom("N"), 1);
      this.addBond(b);
      this.bondResidueIDs.put(b, residue.getResidueID());
      this.rotatableBonds.add(b);
    }
    // if the next residue - the residue with the residueID immediately
    // after this residue's ID - exists, then add a bond between this
    // residue's C and its N.
    if(this.contains(residue.getResidueID()+1)){
      Residue nextResidue = getResidue(residue.getResidueID()+1);
      Bond b = new Bond(residue.getAtom("C"), nextResidue.getAtom("N"), 1);
      this.addBond(b);
      this.bondResidueIDs.put(b, nextResidue.getResidueID());
      this.rotatableBonds.add(b);
    }
  }

  /**
  * Get the residue with the given residueID.
  * @param residueID the residueID of the residue to get
  * @return the Residue with residueID
  * @throws IllegalStateException if there is no residue with residueID
  */
  public Residue getResidue(int residueID){
    if(!contains(residueID)){
      throw new IllegalArgumentException(
          String.format("There is no Residue %d in chain %s", residueID, getChainID()));
    }
    int residueIndex = residueIndices.get(residueID);
    return this.residues.get(residueIndex);
  }

  /**
  * Get the residue that contains atom
  * @param atom the atom we whose containing Residue we want to know
  * @return the Residue that contains atom
  * @throws IllegalStateException if this Atom is not in this chain
  */
  public Residue getResidue(Atom atom){
    if(!contains(atom)){
      throw new IllegalArgumentException(
          String.format("There is no Atom %s in chain %s", atom, getChainID()));
    }
    return this.residuesByAtom.get(atom);
  }

  /**
  * @return the residue with the max residue ID.
  */
  public Residue getLastResidue(){
    return this.residues.get(this.residues.size()-1);
  }

  /**
  * Return a Collection of Residues from the Residue at residueID to the end of the chain.
  * @param residueID The ID of the first residue in the Collection that spans from it to the end of
  * the Chain.
  * @return a Collection of Residues from the Residue at residueID to the end of the chain
  * @throws IllegalArgumentException if the residueID is greater than the maxResidueID in this
  * chain.
  */
  public Collection<Residue> getResiduesToEnd(int residueID){
    // if there is no Residue with this ID in the chain, find the ID of the closest (in sequence)
    // following Residue.
    if(!contains(residueID)){
      int maxResidueID = this.residues.get(residues.size()-1).getResidueID();
      if(residueID > maxResidueID){
        return new ArrayList<Residue>();
        //throw new IllegalArgumentException(String.format("There is no Residue with ID %d in chain."
        //    + " The max residue ID is %d.", residueID, maxResidueID));
      }
      while(residueID >= maxResidueID){
        if(contains(residueID)){
          break;
        } else {
          residueID++;
        }
      }
    }
    int residueIndex = this.residueIndices.get(residueID);
    //System.out.printf("Index of Residue %d: %d\n",residueID, residueIndex);
    return this.residues.subList(residueIndex, this.residues.size());
  }

  /**
  * Checks if there is a residue with the given ID in this PolypeptideChain
  * @param residueID the numeric residue ID of the residue to check for.
  * @return true if there is a residue with the given ID in this
  *         PolypeptideChain, false otherwise
  */
  public boolean contains(int residueID){
    return residueIndices.containsKey(residueID);
  }

  /**
  * @param atom an Atom, maybe in this chain, maybe not
  * @return true if this chain contains atom, false otherwise
  */
  public boolean contains(Atom atom){
    return this.residuesByAtom.containsKey(atom);
  }

  /**
  * @param residue a Residue, maybe in this chain, maybe not
  * @return true if this chain contains residue, false otherwise
  */
  public boolean contains(Residue residue){
    if(contains(residue.getResidueID())){
      Residue resInChain = getResidue(residue.getResidueID());
      return (residue == resInChain);
    }
    return false;
  }

  /**
  * Returns the number of residues in this PolypeptideChain
  * @return the  number of residues in this PolypeptideChain
  */
  public int getNumResidues(){
    return residues.size();
  }

  /**
  * Returns the number of Atoms in this PolypeptideChain
  * @return the number of Atoms in this PolypeptideChain
  */
  public int getNumAtoms(){
    int numAtoms = 0;
    for(Residue res : this.residues){
      numAtoms += res.getNumAtoms();
    }
    return numAtoms;
  }

  /**
  * Returns an array of all the residue IDs of the residues in this chain.
  * @return an array of the residue IDs of the residues in this chain
  */
  public Integer[] getResidueIDs(){
    Integer[] resIDs = new Integer[getNumResidues()];
    int i = 0;
    for(Residue res : this.getResidues()){
      resIDs[i] = res.getResidueID();
      i++;
    }
    return resIDs;
  }

  /**
  * Return the sequence of this PolypeptideChain in the form of a String of the single letter
  * residue names.
  *
  * @return a String containing the single letter residue names with no spaces
  */
  public ProteinSequence getSequence(){
    String seq = "";
    for(Residue res : this){
      seq += res.getOneLetterName();
    }
    return new ProteinSequence(seq);
  }

  private Collection<Residue> getResidues(){
    return this.residues;
  }

  /**
  * Return the Omega Angle for a residue.
  * <p>
  * This residue and the previous residue must both be present.
  * @param residueID the residue ID of the residue to return the omega angle of
  * @return the omega angle - the angle about the (C-1) -- (N) bond, in degrees, or 1000 if
  * undefined
  * @since 0.7.0
  */
  public double getOmegaAngle(int residueID){
    if(!contains(residueID) || !contains(residueID-1)){
      return 1000;
    }
    // get the CA and C from residueID-1 and the N and CA from residueID.
    Residue prev = getResidue(residueID-1);
    Residue current = getResidue(residueID);
    return calculateDihedralAngle(prev.getAtom("CA").getCoordinates(),
                                  prev.getAtom("C").getCoordinates(),
                                  current.getAtom("N").getCoordinates(),
                                  current.getAtom("CA").getCoordinates());
  }

  /**
  * Return the Phi Angle for a residue.
  * <p>
  * This residue and the previous residue must both be present.
  * @param residueID the residue ID of the residue to return the omega angle of
  * @return the phi angle - the angle about the (N) -- (CA) bond, in degrees, or 1000 if undefined
  * @since 0.7.0
  */
  public double getPhiAngle(int residueID){
    if(!contains(residueID) || !contains(residueID-1)){
      return 1000;
      //throw new IllegalStateException(String.format("Both residues %d and %d must be present\n"
      //    , residueID, residueID-1));
    }
    // get the C from residueID-1 and the N, CA, and C from residueID.
    Residue prev = getResidue(residueID-1);
    Residue current = getResidue(residueID);
    return calculateDihedralAngle(prev.getAtom("C").getCoordinates(),
                                  current.getAtom("N").getCoordinates(),
                                  current.getAtom("CA").getCoordinates(),
                                  current.getAtom("C").getCoordinates());
  }

  /**
  * Return the Psi Angle for a residue.
  * <p>
  * This residue and the next residue must both be present.
  * @param residueID the residue ID of the residue to return the omega angle of
  * @return the psi angle - the angle about the (CA) -- (C) bond, in degrees, or 1000 if undefined
  * @since 0.7.0
  */
  public double getPsiAngle(int residueID){
    if(!contains(residueID) || !contains(residueID+1)){
      return 1000;
    }
    // get the N, CA and C from residueID and the N from residueID+1.
    Residue current = getResidue(residueID);
    Residue next = getResidue(residueID+1);
    return calculateDihedralAngle(current.getAtom("N").getCoordinates(),
                                  current.getAtom("CA").getCoordinates(),
                                  current.getAtom("C").getCoordinates(),
                                  next.getAtom("N").getCoordinates());
  }

  /**
  * @return an {@code List<Double>} of all the dihedral angles in this PolypeptideChain. The
  * angles are listed in order by residue from N to C terminus and for each residue in the order of
  * omega, phi, and psi.
  */
  public List<Double> getBackboneDihedralAngles(){
    ArrayList<Double> dihedralAngles = new ArrayList<Double>(this.getNumResidues()*3);
    for(Residue res : this){
      dihedralAngles.add(getOmegaAngle(res.getResidueID()));
      dihedralAngles.add(getPhiAngle(res.getResidueID()));
      dihedralAngles.add(getPsiAngle(res.getResidueID()));
    }
    return dihedralAngles;
  }

  /*
  * Calculate and return the angle between the ab vector and the cd vector about the bc vector.
  */
  private static double calculateDihedralAngle(Vector3D a, Vector3D b, Vector3D c, Vector3D d){
    return -1.0 * Vector3D.calculateDihedralAngle(a,b,c,d);
  }


  @Override
  public void applyTransformation(Transformation t){
    for(Residue residue : this ){
      residue.applyTransformation(t);
    }
  }

  @Override
  public Iterator<Residue> iterator(){
    return this.getResidues().iterator();
  }

  private class ResidueComparator implements Comparator<Residue>{
    public int compare(Residue one, Residue two){
      return one.getResidueID() - two.getResidueID();
    }
  }
}
