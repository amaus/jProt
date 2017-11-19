package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.sequence.*;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.linearAlgebra.Vector3D;
import com.aaronpmaus.jMath.transformations.Transformable;
import com.aaronpmaus.jMath.transformations.Transformation;

import java.util.HashMap;
import java.util.ArrayList;
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
  private ArrayList<Residue> residues;
  // key: residueID, value: residue index in residues list
  private HashMap<Integer, Integer> residueIndices;
  private String chainID;

  /**
  * A constructor for a PolypeptideChain.
  * @param chainID the letter ID of this chain
  */
  public PolypeptideChain(String chainID){
    super();
    residues = new ArrayList<Residue>();
    residueIndices = new HashMap<Integer,Integer>();
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
    }
    // if the previous residue - the residue with the residueID immediately
    // preceeding this residue's ID - exists, then add a bond between it's
    // C and this residue's N.
    if(this.containsResidue(residue.getResidueID()-1)){
      Residue prevResidue = getResidue(residue.getResidueID()-1);
      this.addBond(new Bond(prevResidue.getAtom("C"), residue.getAtom("N"), 1));
    }
    // if the next residue - the residue with the residueID immediately
    // after this residue's ID - exists, then add a bond between this
    // residue's C and its N.
    if(this.containsResidue(residue.getResidueID()+1)){
      Residue nextResidue = getResidue(residue.getResidueID()+1);
      this.addBond(new Bond(residue.getAtom("C"), nextResidue.getAtom("N"), 1));
    }
  }

  /**
  * Gets the residue with the given residueID.
  * @param residueID the residueID of the residue to get
  * @return the Residue if it is in this PolypeptideChain, null otherwise
  */
  public Residue getResidue(int residueID){
    if(this.containsResidue(residueID)){
      int residueIndex = residueIndices.get(residueID);
      return this.residues.get(residueIndex);
    }
    return null;
  }

  /**
  * Checks if there is a residue with the given ID in this PolypeptideChain
  * @param residueID the numeric residue ID of the residue to check for.
  * @return true if there is a residue with the given ID in this
  *         PolypeptideChain, false otherwise
  */
  public boolean containsResidue(int residueID){
    return residueIndices.containsKey(residueID);
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

  /**
  * Add hydrogens to this PolypeptideChain.
  */
  public void enableHydrogens(){
    // for each residue in this chain, get add all its bonds (now containing hydrogens) to
    // this chain.
    for(Residue residue : this){
      residue.enableHydrogens();
      addBondsFromResidue(residue);
    }
  }

  /**
  * Remove all hydrogens from this chain.
  */
  public void disableHydrogens(){
    for(Residue residue : this){
      for(Bond bond : residue.getBondsToHydrogens()){
        this.removeBond(bond); // from superclass Molecule
      }
      for(Atom hydrogen : residue.getHydrogens()){
        this.removeAtom(hydrogen); // from superclass Molecule
      }
      residue.disableHydrogens();
    }
  }

  private Collection<Residue> getResidues(){
    return this.residues;
  }

  /**
  * Return the Phi Angle for a residue.
  * <p>
  * This residue and the previous residue must both be present.
  * @param residueID the residue ID of the residue to return the omega angle of
  * @return the phi angle - the angle about the (N) -- (CA) bond, in degrees
  * @throws IllegalStateException if either residue at residueID or residueID-1 are not present
  */
  public double getPhiAngle(int residueID){
    if(!containsResidue(residueID) || !containsResidue(residueID-1)){
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
  * @return the psi angle - the angle about the (CA) -- (C) bond, in degrees
  * @throws IllegalStateException if either residue at residueID or residueID+1 are not present
  */
  public double getPsiAngle(int residueID){
    if(!containsResidue(residueID) || !containsResidue(residueID+1)){
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
  * Return the Omega Angle for a residue.
  * <p>
  * This residue and the previous residue must both be present.
  * @param residueID the residue ID of the residue to return the omega angle of
  * @return the omega angle - the angle about the (C-1) -- (O) bond, in degrees
  * @throws IllegalStateException if either residue at residueID or residueID-1 are not present
  */
  public double getOmegaAngle(int residueID){
    if(!containsResidue(residueID) || !containsResidue(residueID-1)){
      return 1000;
    }
    // get the CA and C from residueID-1 and the N and CA from residueID.
    Residue prev = getResidue(residueID-1);
    Residue current = getResidue(residueID);
    /*return calculateDihedralAngle(prev.getAtom("CA").getCoordinates(),
                                  prev.getAtom("C").getCoordinates(),
                                  current.getAtom("N").getCoordinates(),
                                  current.getAtom("CA").getCoordinates());*/
    return calculateDihedralAngle(current.getAtom("CA").getCoordinates(),
                                  current.getAtom("N").getCoordinates(),
                                  prev.getAtom("C").getCoordinates(),
                                  prev.getAtom("CA").getCoordinates());
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
