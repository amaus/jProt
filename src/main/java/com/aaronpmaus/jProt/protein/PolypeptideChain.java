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
  // key: Bond, value: residueID
  private HashMap<Bond, Integer> bondResidueIDs;
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
    if(this.containsResidue(residue.getResidueID()-1)){
      Residue prevResidue = getResidue(residue.getResidueID()-1);
      Bond b = new Bond(prevResidue.getAtom("C"), residue.getAtom("N"), 1);
      this.addBond(b);
      this.bondResidueIDs.put(b, residue.getResidueID());
      this.rotatableBonds.add(b);
    }
    // if the next residue - the residue with the residueID immediately
    // after this residue's ID - exists, then add a bond between this
    // residue's C and its N.
    if(this.containsResidue(residue.getResidueID()+1)){
      Residue nextResidue = getResidue(residue.getResidueID()+1);
      Bond b = new Bond(residue.getAtom("C"), nextResidue.getAtom("N"), 1);
      this.addBond(b);
      this.bondResidueIDs.put(b, nextResidue.getResidueID());
      this.rotatableBonds.add(b);
    }
  }

  /**
  * Gets the residue with the given residueID.
  * @param residueID the residueID of the residue to get
  * @return the Residue if it is in this PolypeptideChain, null otherwise
  */
  public Residue getResidue(int residueID){
    if(containsResidue(residueID)){
      int residueIndex = residueIndices.get(residueID);
      return this.residues.get(residueIndex);
    }
    return null;
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
    if(!containsResidue(residueID)){
      int maxResidueID = this.residues.get(residues.size()-1).getResidueID();
      if(residueID > maxResidueID){
        return new ArrayList<Residue>();
        //throw new IllegalArgumentException(String.format("There is no Residue with ID %d in chain."
        //    + " The max residue ID is %d.", residueID, maxResidueID));
      }
      while(residueID >= maxResidueID){
        if(containsResidue(residueID)){
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
  * Set the Omega Angle of a Residue. Valid for every residue except the N-terminus.
  * @param residueID the residue to set the angle of, must not be the first residue in the chain.
  * @param degrees the value to set it to
  * @since 0.7.0
  */
  public void setOmegaAngle(int residueID, double degrees){
    double currentAngle = getOmegaAngle(residueID);
    double delta = degrees - currentAngle;
    rotateOmegaAngle(residueID, delta);
  }

  /**
  * Set the Phi Angle of a Residue. Valid for every residue except the N-terminus.
  * @param residueID the residue to set the angle of, must not be the first residue in the chain.
  * @param degrees the value to set it to
  * @since 0.7.0
  */
  public void setPhiAngle(int residueID, double degrees){
    double currentAngle = getPhiAngle(residueID);
    double delta = degrees - currentAngle;
    rotatePhiAngle(residueID, delta);
  }

  /**
  * Set the Psi Angle of a Residue. Valid for every residue except the C-terminus.
  * @param residueID the residue to set the angle of, must not be the last residue in the chain.
  * @param degrees the value to set it to
  * @since 0.7.0
  */
  public void setPsiAngle(int residueID, double degrees){
    double currentAngle = getPsiAngle(residueID);
    double delta = degrees - currentAngle;
    rotatePsiAngle(residueID, delta);
  }

  /**
  * Rotate the Omega angle of the specified Residue by the given degrees.
  * <p>
  * Both the Residue and the Residue immediately before it must be present in the Chain.
  * @param residueID the Residue whose Omega Angle to rotate
  * @param degrees the number of degrees to rotate that angle
  * @since 0.7.0
  */
  public void rotateOmegaAngle(int residueID, double degrees){
    if(!containsResidue(residueID) || !containsResidue(residueID-1)){
      throw new IllegalArgumentException(
          String.format("Can not rotate about the Omega Angle of residue %d."
          + "Both residues %d and %d must be in chain", residueID, residueID-1));
    }
    Residue prev = getResidue(residueID-1);
    Residue res = getResidue(residueID);
    //System.out.printf("Rotating about the Omega Bond of Residue %d, Prev: %d\n", res.getResidueID(), prev.getResidueID());
    //System.out.printf("C: %d, N: %d\n", prev.getAtom("C").getSerialNumber(), res.getAtom("N").getSerialNumber());
    rotateAboutBond(prev.getAtom("C"), res.getAtom("N"), residueID, "Omega", degrees);
  }

  /**
  * Rotate the Phi angle of the specified Residue by the given degrees.
  * <p>
  * The Residue must be present in the Chain.
  * @param residueID the Residue whose Phi Angle to rotate
  * @param degrees the number of degrees to rotate that angle
  * @since 0.7.0
  */
  public void rotatePhiAngle(int residueID, double degrees){
    if(!containsResidue(residueID)){
      throw new IllegalArgumentException(
          String.format("Can not rotate Phi Angle of residue %d, not in chain", residueID));
    }
    Residue res = getResidue(residueID);
    rotateAboutBond(res.getAtom("N"), res.getAtom("CA"), residueID, "Phi", degrees);
  }

  /**
  * Rotate the Psi angle of the specified Residue by the given degrees.
  * <p>
  * The Residue must be present in the Chain.
  * @param residueID the Residue whose Psi Angle to rotate
  * @param degrees the number of degrees to rotate that angle
  * @since 0.7.0
  */
  public void rotatePsiAngle(int residueID, double degrees){
    if(!containsResidue(residueID)){
      throw new IllegalArgumentException(
          String.format("Can not rotate Psi Angle of residue %d, not in chain", residueID));
    }
    Residue res = getResidue(residueID);
    rotateAboutBond(res.getAtom("CA"), res.getAtom("C"), residueID, "Psi", degrees);
  }

  /**
  * Rotate a PolypeptideChain about one of its rotatable bonds, one of its defined dihedral angles.
  * <p>
  * All atoms after that bond are transformed.
  * @param bond the bond to rotate about
  * @param degrees the number of degrees to rotate the residue
  * @throws IllegalArgumentException if the bond is not a rotatable bond, an omega, phi, or psi bond
  * @since 0.7.0
  */
  public void rotateAboutBond(Bond bond, double degrees){
    // break this into two methods, one that takes a bond, and one that takes two atoms
    Integer residueID = this.bondResidueIDs.get(bond);
    if(residueID == null){
      throw new IllegalArgumentException("%s not a rotatable bond. Must be C-N, N-CA, or CA-C.");
    }
    Atom atomOne = bond.getAtomOne();
    Atom atomTwo = bond.getAtomTwo();
    // The axis of the bond is the vector from atomOne to atomTwo.
    // Ensure that the atom assignments are correct so the axis points in the correct direction.
    String bondType = null;
    if(bond.containsAtom("C") && bond.containsAtom("N")){
      bondType = "Omega";
      if(atomOne.getAtomName().equals("N")){
        Atom temp = atomOne;
        atomOne = atomTwo;
        atomTwo = temp;
      }
    } else if(bond.containsAtom("N") && bond.containsAtom("CA")){
      bondType = "Phi";
      if(atomOne.getAtomName().equals("CA")){
        Atom temp = atomOne;
        atomOne = atomTwo;
        atomTwo = temp;
      }
    } else if(bond.containsAtom("CA") && bond.containsAtom("C")){
      bondType = "Psi";
      if(atomOne.getAtomName().equals("C")){
        Atom temp = atomOne;
        atomOne = atomTwo;
        atomTwo = temp;
      }
    }

    rotateAboutBond(atomOne, atomTwo, residueID, bondType, degrees);
  }

  private void rotateAboutBond(Atom atomOne, Atom atomTwo, int residueID,
      String bondType, double degrees){

    Vector3D coor1 = atomOne.getCoordinates();
    Vector3D coor2 = atomTwo.getCoordinates();
    Vector3D axis = coor2.subtract(coor1);
    Transformation t = new Transformation();
    t.addRotationAboutAxis(coor1, axis, degrees);
    // Different bonds require different atoms to be rotated.
    switch(bondType){
      case "Omega":
        applyTransformationToOmegaBond(t, residueID);
        break;
      case "Phi":
        applyTransformationToPhiBond(t, residueID);
        break;
      case "Psi":
        applyTransformationToPsiBond(t, residueID);
        break;
    }
  }

  private void applyTransformationToOmegaBond(Transformation t, int residueID){
    // In the Omega case, apply this transformation to this residue and every one after
    // it.
    for(Residue res : getResiduesToEnd(residueID)){
      //System.out.printf("Rotating Residue %d\n", res.getResidueID());
      res.applyTransformation(t);
    }
  }

  private void applyTransformationToPhiBond(Transformation t, int residueID){
    // In the Phi case, apply the transformation to all atoms in this residue except for
    // N, H, and CA. Then apply the transformation to every residue after.
    ArrayList<Atom> atomsInThisRes = new ArrayList<Atom>(getResidue(residueID).getHeavyAtoms());
    atomsInThisRes.addAll(getResidue(residueID).getHydrogens());
    Iterator<Atom> it = atomsInThisRes.iterator();
    while(it.hasNext()){
      Atom atom = it.next();
      if(atom.getAtomName().equals("H")
          || atom.getAtomName().equals("N")
          || atom.getAtomName().equals("CA")){
        it.remove();
      }
    }
    for(Atom atom : atomsInThisRes){
      atom.applyTransformation(t);
    }
    // If the next residue is in a gap, getResiduesToEnd will return all Residues after the gap.
    for(Residue res : getResiduesToEnd(residueID+1)){
      res.applyTransformation(t);
    }
  }

  private void applyTransformationToPsiBond(Transformation t, int residueID){
    // in the Psi case, apply the transformation to the O and OXT (if present) atoms,
    // and then to every residue after this one.
    Residue residue = getResidue(residueID);
    if(residue.contains("O")){
      residue.getAtom("O").applyTransformation(t);
    }
    if(residue.isCarboxylTerminus()){
      residue.getAtom("OXT").applyTransformation(t);
    }
    // If the next residue is in a gap, getResiduesToEnd will return all Residues after the gap.
    for(Residue res : getResiduesToEnd(residueID+1)){
      res.applyTransformation(t);
    }
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
    if(!containsResidue(residueID) || !containsResidue(residueID-1)){
      return 1000;
    }
    // get the CA and C from residueID-1 and the N and CA from residueID.
    Residue prev = getResidue(residueID-1);
    Residue current = getResidue(residueID);
    return calculateDihedralAngle(prev.getAtom("CA").getCoordinates(),
                                  prev.getAtom("C").getCoordinates(),
                                  current.getAtom("N").getCoordinates(),
                                  current.getAtom("CA").getCoordinates());
    /*return calculateDihedralAngle(current.getAtom("CA").getCoordinates(),
                                  current.getAtom("N").getCoordinates(),
                                  prev.getAtom("C").getCoordinates(),
                                  prev.getAtom("CA").getCoordinates());*/
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
  * @return the psi angle - the angle about the (CA) -- (C) bond, in degrees, or 1000 if undefined
  * @since 0.7.0
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
