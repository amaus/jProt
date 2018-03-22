package com.aaronpmaus.jProt.manipulators;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.linearAlgebra.Vector3D;
import com.aaronpmaus.jMath.transformations.Transformation;

import com.aaronpmaus.jProt.protein.*;

import java.util.HashSet;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Comparator;
import java.util.Collections;
import java.io.IOException;

/**
* A ConformationManipulator operates on a Protein and allows that Protein's conformation to be
* altered. The three types of manipulations that are allowed are: <br>
* - Bond Length Manipulations <br>
* - Angle Manipulations <br>
* - Dihedral Angle Manipulations <p>
* bond length manipulations allow the legnth of a Bond to be altered, angle manipulations allow the
* angle formed by 3 consecutive atoms to be modified, and dihedral angle manipulations allow
* dihedral angles to be altered.
* <p>
* In each of the manipulations, the ordering of the atoms that define the length, angle, or dihedral
* angle is important because it will affect the Transformation needed to perform the manipulation
* and it will determine which atoms to apply that Transformation to. The ordered imposed on all
* atoms is the conventional ordering within a protein chain, from the N-Terminus to the C-Terminus,
* and if the atoms are within a side chain, from CA on back.
* <p>
* For Angle and Dihedral Angle manipulations, the Atom(s) further from the N-Terminus will be the
* ones transformed. For example, if a simple angle is defined by N-CA-C, C will be the atom rotated
* to modify that angle.
* <p>
* IMPORTANT NOTE: if a ConformationManipulator has been constructed and the Protein it operates on
* has been modified, then the update() method must be called or the manipulator's behavior becomes
* undefined. Modifications that necessitate a call to update() are the insertion or removal of
* Polypeptide Chains, Residues, or Atoms.
* <p>
* Usage:<br>
* {@code // Construct a custom protein from the sequence IAMSTARSTFF}<br>
* {@code Protein prot = VirtualRibosome.synthesizeProtein(new ProteinSequence("IAMSTARSTFF"), "strstf");}<br>
* {@code // and build a manipulator to operate on it.}<br>
* {@code ConformationManipulator manip = new CascadeConformationManipulator(prot);}<br>
* {@code // Set the phi angle of residue 2 in chain A to 141.42 degrees}<br>
* {@code manip.setPhiAngle("A", 2, 141.42);}<br>
* {@code // Get the penultimate phenylalanine and rotate its benzyl group -90 degrees}<br>
* {@code Residue phe = chain.getResidue(10);}<br>
* {@code // rotate all atoms after the CB-CG bond about the axis formed by the bond by -90 degrees}<br>
* {@code manip.rotateAboutBond(phe.getAtom("CB"), phe.getAtom("CG"), -90);}<br>
* @version 0.7.0
* @since 0.7.0
*/
public abstract class ConformationManipulator {
  // The Protein that this ConformationManipulator manipulator manipulates.
  private Protein protein;
  // these indices are of the atoms from N-terminus to C-terminus and from alpha level back within
  // residues. They allow ConformationManipulator to ensure consistent ordering when calculating
  // angles. This will include all atoms from all chains: all atoms in the first chain, then all
  // atoms in the second chain and so on. The ordering of the chains doesn't matter. These indices
  // will be used to compare consecutive covalently bonded atoms.
  private HashMap<Atom,Integer> atomIndices;
  // All valid angle triplets
  private HashSet<AngleTriplet> validTriplets;
  // All valid dihedral angle quartets
  private HashSet<DihedralQuartet> validQuartets;
  // a mapping between each rotatable bond and the quartet of atoms that define its dihedral angle.
  private HashMap<OrderedBond, DihedralQuartet> dihedrals;
  // All valid bonds.
  private HashSet<OrderedBond> validBonds;

  /**
  * Construct a ConformationManipulator to manipulate a Protein.
  * @param protein the Protein to manipulate
  */
  public ConformationManipulator(Protein protein){
    this.protein = protein;
    this.atomIndices = constructAtomIndices();
    this.validTriplets = buildAngleTriplets();
    this.validQuartets = new HashSet<DihedralQuartet>();
    this.dihedrals = buildDihedralQuartets();
    this.validBonds = buildBonds();
  }

  /**
  * This method must be called if the Protein this manipulator operates on has been modified
  * (adding Polypeptide Chains, Residues, or Atoms), otherwise the behavior of this manipulator
  * becomes undefined.
  */
  public void update(){
    this.atomIndices = constructAtomIndices();
    this.validTriplets = buildAngleTriplets();
    this.validQuartets = new HashSet<DihedralQuartet>();
    this.dihedrals = buildDihedralQuartets();
    this.validBonds = buildBonds();
  }

  //###########################################################################################\\
  //#############################  THREE ATOM ANGLE MANIPULATIONS  ############################\\
  //###########################################################################################\\

  /**
  * Set the angle of the three atoms to that specified by degrees.
  * <p>
  * The three Atoms must a be valid Angle Triplet. A valid angle triplet is a set of three
  * consecutive atoms that obey the ordering of atoms in the protein and where at most two of them
  * can be in a ring. The ordering of atoms is from N to C terminus and from CA back to the furthest
  * atom in a side chain. Put another way, the Atoms can be ordered by their distance from the
  * N-Terminus.
  * <p>
  * Examples: CA-C-O is valid angle triplet. O-C-N is not a valid angle triplet because
  * C is closer to the N-terminus than O is. In phenylalanine, CB-CG-CD1 is a valid triplet because
  * only CG and CD1 are part of the ring. CG-CD1-CE1 is not valid because all three are part of the
  * ring.
  * @param atomOne the first atom of the three
  * @param atomTwo the second atom of the three
  * @param atomThree the third atom of the three
  * @param degrees the degrees to set the angle to.
  * @throws IllegalArgumentException if the three atoms are not a valid angle triplet
  */
  public void setAngle(Atom atomOne, Atom atomTwo, Atom atomThree, double degrees){
    double currentAngle = getAngle(atomOne, atomTwo, atomThree);
    double delta = degrees - currentAngle;
    modifyAngle(atomOne, atomTwo, atomThree, delta);
  }

  /**
  * Manipulate the three atoms that define an angle to adjust that angle by degrees.
  * <p>
  * The three Atoms must a be valid Angle Triplet. A valid angle triplet is a set of three
  * consecutive atoms that obey the ordering of atoms in the protein and where at most two of them
  * can be in a ring. The ordering of atoms is from N to C terminus and from CA back to the furthest
  * atom in a side chain. Put another way, the Atoms can be ordered by their distance from the
  * N-Terminus.
  * <p>
  * Examples: CA-C-O is valid angle triplet. O-C-N is not a valid angle triplet because
  * C is closer to the N-terminus than O is. In phenylalanine, CB-CG-CD1 is a valid triplet because
  * only CG and CD1 are part of the ring. CG-CD1-CE1 is not valid because all three are part of the
  * ring.
  * @param atomOne the first atom of the three
  * @param atomTwo the second atom of the three
  * @param atomThree the third atom of the three
  * @param degrees the degrees to adjust the angle by
  * @throws IllegalArgumentException if the three atoms are not a valid angle triplet
  */
  public abstract void modifyAngle(Atom atomOne, Atom atomTwo, Atom atomThree, double degrees);

  /**
  * Package Private (intended only accessible to this class and subclasses of this class, but java
  * provides no such access modifier) method to construct a Transformation to rotate the third atom
  * in a three atom triplet in order to manipulate the angle formed by those three atoms. In more
  * technical terms, the third atom is rotated about the normal formed by the three by the number of
  * degrees specified.
  * <p>
  * The directionality of the atoms is important. They should be specified as they would appear in
  * the Protein from N-Terminus to C-Terminus and if within a side chain from CA on back. If not,
  * they will be reorderd and the angle calculated from the required ordering.
  * <p>
  * For example, if the atoms specified are the C, CA, and N or Residue 42, the atom rotated will be
  * C. It will be rotated about the normal formed by N, CA, and C by the number of degrees
  * specified.
  * @param atomOne the first atom of the triplet
  * @param atomTwo the second atom of the triplet
  * @param atomThree the third atom of the triplet
  * @param degrees the degrees to rotate the third atom by
  * @return a Transformation that can be applied to the third atom in the angle triplet to rotate it
  * about the normal formed by the three atoms by the specified degrees
  * @throws IllegalArgumentException if the three atoms are not a valid angle triplet
  */
  Transformation buildAngleTransformation(Atom atomOne, Atom atomTwo,
                                          Atom atomThree, double degrees){
    if(!isValidAngleTriplet(atomOne, atomTwo, atomThree)){
      throw new IllegalArgumentException(
          String.format("%s-%d, %s-%d, %s-%d is not a valid defined angle triplet.\n",
                        atomOne.getName(), atomOne.getSerialNumber(),
                        atomTwo.getName(), atomTwo.getSerialNumber(),
                        atomThree.getName(), atomThree.getSerialNumber()));
    }
    AngleTriplet triplet = new AngleTriplet(atomOne, atomTwo, atomThree);
    Vector3D atomOneCoor = triplet.getAtomOne().getCoordinates();
    Vector3D atomTwoCoor = triplet.getAtomTwo().getCoordinates();
    Vector3D atomThreeCoor = triplet.getAtomThree().getCoordinates();

    Vector3D normal = atomOneCoor.subtract(atomTwoCoor)
        .crossProduct(atomThreeCoor.subtract(atomTwoCoor));
    Transformation t = new Transformation();
    t.addRotationAboutAxis(atomTwoCoor, normal, degrees);

    return t;
  }


  //###########################################################################################\\
  //##############################  DIHEDRAL ANGLE MANIPULATIONS  #############################\\
  //###########################################################################################\\

  /**
  * Given two atom names that specify a rotatable bond, set that bond's dihedral angle as specified.
  * <p>
  * To calculate the dihedral angle, the order of the atoms used will be by their distance from the
  * N-Terminus.
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @param degrees the degrees to set the dihedral angle to
  * @throws IllegalArgumentException if the bond is not a rotatable bond
  */
  public void setDihedralAngle(Atom atomOne, Atom atomTwo, double degrees){
    DihedralQuartet quartet = getDihedralQuartet(atomOne, atomTwo);
    setDihedralAngle(quartet.getAtomOne(), quartet.getAtomTwo(),
                     quartet.getAtomThree(), quartet.getAtomFour(), degrees);
  }

  /**
  * Given four atoms of a dihedral angle, set that dihedral angle as specified.
  * <p>
  * To calculate the dihedral angle, the order of the atoms used will be by their distance from the
  * N-Terminus.
  * @param atomOne the first atom of the angle
  * @param atomTwo the second atom of the angle
  * @param atomThree the third atom of the angle
  * @param atomFour the fourth atom of the angle
  * @param degrees the degrees to set the dihedral angle to
  * @throws IllegalArgumentException if the atoms do not define a valid dihedral angle
  */
  public void setDihedralAngle(Atom atomOne, Atom atomTwo,
                               Atom atomThree, Atom atomFour, double degrees){
    double currentAngle = getDihedralAngle(atomOne, atomTwo, atomThree, atomFour);
    double delta = degrees - currentAngle;
    modifyDihedralAngle(atomOne, atomTwo, atomThree, atomFour, delta);
  }

  /**
  * Set the backbone Omega Dihedral Angle of a Residue.
  * @param chainID the ID of the chain the Residue belongs to
  * @param residueID the ID of the Residue
  * @param degrees the degrees to set the angle to
  * @throws IllegalStateException if the chain, residue, or previous residue are not in the Protein.
  */
  public void setOmegaAngle(String chainID, int residueID, double degrees){
    if(!this.protein.contains(chainID)){
      throw new IllegalStateException("Chain " + chainID + " not in Protein.");
    }
    PolypeptideChain chain = this.protein.getChain(chainID);
    if(!chain.contains(residueID - 1) || !chain.contains(residueID)){
      throw new IllegalStateException(String.format("Both residues %d and %d must be in Chain %s."
                                                  + " Can't rotate omega angle of Residue %d.\n"
                                                  ,residueID-1, residueID, chainID));
    }
    Residue previousResidue = chain.getResidue(residueID - 1);
    Residue nextResidue = chain.getResidue(residueID);
    if(!previousResidue.contains("C") || !nextResidue.contains("N")){
      throw new IllegalStateException(
          String.format("Both Atoms C and N must be in Protein to rotate about that bond."));
    }
    Atom atomOne = previousResidue.getAtom("C");
    Atom atomTwo = nextResidue.getAtom("N");
    setDihedralAngle(atomOne, atomTwo, degrees);
  }

  /**
  * Set the backbone Phi Dihedral Angle of a Residue.
  * @param chainID the ID of the chain the Residue belongs to
  * @param residueID the ID of the Residue
  * @param degrees the degrees to set the angle to
  * @throws IllegalStateException if the chain, residue, or previous residue are not in the Protein.
  */
  public void setPhiAngle(String chainID, int residueID, double degrees){
    if(!this.protein.contains(chainID)){
      throw new IllegalStateException("Chain " + chainID + " not in Protein.");
    }
    PolypeptideChain chain = this.protein.getChain(chainID);
    if(!chain.contains(residueID - 1) || !chain.contains(residueID)){
      throw new IllegalStateException(String.format("Both residues %d and %d must be in Chain %s."
                                                  + " Can't rotate omega angle of Residue %d.\n"
                                                  ,residueID-1, residueID, chainID));
    }
    Residue residue = chain.getResidue(residueID);
    if(!residue.contains("N") || !residue.contains("CA")){
      throw new IllegalStateException(
          String.format("Both Atoms N and CA must be in Protein to rotate about that bond."));
    }
    Atom atomOne = residue.getAtom("N");
    Atom atomTwo = residue.getAtom("CA");
    setDihedralAngle(atomOne, atomTwo, degrees);
  }

  /**
  * Set the backbone Psi Dihedral Angle of a Residue.
  * @param chainID the ID of the chain the Residue belongs to
  * @param residueID the ID of the Residue
  * @param degrees the degrees to set the angle to
  * @throws IllegalStateException if the chain, residue, or next residue are not in the Protein.
  */
  public void setPsiAngle(String chainID, int residueID, double degrees){
    if(!this.protein.contains(chainID)){
      throw new IllegalStateException("Chain " + chainID + " not in Protein.");
    }
    PolypeptideChain chain = this.protein.getChain(chainID);
    if(!chain.contains(residueID - 1) || !chain.contains(residueID)){
      throw new IllegalStateException(String.format("Both residues %d and %d must be in Chain %s."
                                                  + " Can't rotate omega angle of Residue %d.\n"
                                                  ,residueID-1, residueID, chainID));
    }
    Residue residue = chain.getResidue(residueID);
    if(!residue.contains("CA") || !residue.contains("C")){
      throw new IllegalStateException(
          String.format("Both Atoms CA and C must be in Protein to rotate about that bond."));
    }
    Atom atomOne = residue.getAtom("CA");
    Atom atomTwo = residue.getAtom("C");
    setDihedralAngle(atomOne, atomTwo, degrees);
  }

  /**
  * Rotate a residue about one of its rotatable bonds, one of its defined dihedral angles.
  * <p>
  * To calculate the dihedral angle, the order of the atoms used will be by their distance from the
  * N-Terminus.
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @param degrees the number of degrees to rotate about the bond
  * @throws IllegalArgumentException if the bond is not a rotatable bond
  */
  public void rotateAboutBond(Atom atomOne, Atom atomTwo, double degrees){
    DihedralQuartet quartet = getDihedralQuartet(atomOne, atomTwo);
    modifyDihedralAngle(quartet.getAtomOne(), quartet.getAtomTwo(),
                        quartet.getAtomThree(), quartet.getAtomFour(), degrees);
  }

  /**
  * Manipulate one of the Protein's defined dihedral angles.
  * <p>
  * To calculate the dihedral angle, the order of the atoms used will be by their distance from the
  * N-Terminus.
  * @param atomOne the first atom in the quartet
  * @param atomTwo the second atom in the quartet
  * @param atomThree the third atom in the quartet
  * @param atomFour the fourth atom in the quartet
  * @param degrees the number of degrees to rotate about the bond
  * @throws IllegalArgumentException if the bond is not a rotatable bond
  */
  public abstract void modifyDihedralAngle(Atom atomOne, Atom atomTwo,
                                       Atom atomThree, Atom atomFour,
                                       double degrees);

  /**
  * Package Private method to construct a Transformation to rotate the bond of a dihedral angle
  * by the number of degrees specified.
  * <p>
  * The directionality of the atoms is important. They should be specified as they would appear in
  * the Protein from N-Terminus to C-Terminus and if within a side chain from CA on back. In other
  * words, they will be ordered by their distance from the N-Terminus. If not, they will be
  * re-orderd and the angle calculated from the required ordering.
  * <p>
  * For example, if the bond is specified by the atoms CA and N of Residue 42, the N-CA bond will
  * form the axis of rotation.
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @param degrees the number of degrees to rotate about the bond
  * @return a Transformation that can be applied to the atoms after this bond to rotate them about
  * the axis of the bond by the specified degrees
  * @throws IllegalArgumentException if the bond is not a rotatable bond
  */
  Transformation buildDihedralRotationTransformation(Atom atomOne, Atom atomTwo,
                                                     Atom atomThree, Atom atomFour,
                                                     double degrees){
    DihedralQuartet quartet = new DihedralQuartet(atomOne, atomTwo, atomThree, atomFour);
    if(!isValidDihedralQuartet(quartet)){
      throw new IllegalArgumentException(
          String.format("%s-%s-%s-%s is not a valid dihedral quartet.\n Either the Atoms are not"
                      + "consecutive or one of them is missing from the Protein.",
                      atomOne.getName(), atomTwo.getName(),
                      atomThree.getName(), atomFour.getName()));
    }
    Transformation t = new Transformation();
    Vector3D atomTwoCoor = quartet.getAtomTwo().getCoordinates();
    Vector3D atomThreeCoor = quartet.getAtomThree().getCoordinates();
    t.addRotationAboutAxis(atomTwoCoor, atomThreeCoor.subtract(atomTwoCoor), degrees);
    return t;
  }

  //###########################################################################################\\
  //###################################  BOND MANIPULATIONS  ##################################\\
  //###########################################################################################\\

  /**
  * Set the length of one of the Bonds of the Protein.
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @param length the length in Angstroms to set the bond to
  * @throws IllegalArgumentException if the atoms are not a valid bond in the Protein
  */
  public void setBondLength(Atom atomOne, Atom atomTwo, double length){
    if(!isValidBond(atomOne, atomTwo)){
      throw new IllegalArgumentException(
          String.format("%s-%d, %s-%d is not a valid bond.",
              atomOne.getName(), atomOne.getSerialNumber(),
              atomTwo.getName(), atomTwo.getSerialNumber()));
    }
    OrderedBond bond = new OrderedBond(atomOne, atomTwo);
    double currentDistance = atomOne.distance(atomTwo);
    double delta = length - currentDistance;
    modifyBondLength(bond.getAtomOne(), bond.getAtomTwo(), delta);
  }

  /**
  * Adjust the length of one of the Bonds of the Protein.
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @param delta the delta in Angstroms to adjust the bond length by
  */
  public abstract void modifyBondLength(Atom atomOne, Atom atomTwo, double delta);

  /**
  * Package Private method to construct a Transformation to modify the length of a bond.
  * <p>
  * This method translates atomTwo by delta towards or away from atomOne depending on the sign of
  * delta.
  * <p>
  * In the other build transformation methods in this class, an ordering on the atoms is imposed.
  * That is not the case in this method. As long as the two atoms are bonded, the second one
  * specified will be the one transformed.
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @param delta the delta in Angstroms to adjust the bond length by
  * @return a Transformation that can be applied to atomTwo to translation it by delta distance
  * towards or away from atomOne
  * @throws IllegalArgumentException if the atoms are not a valid bond in the Protein
  */
  Transformation buildBondManipulationTransformation(Atom atomOne, Atom atomTwo, double delta){
    if(!isValidBond(atomOne, atomTwo)){
      throw new IllegalArgumentException(
          String.format("%s-%d, %s-%d is not a valid bond.",
              atomOne.getName(), atomOne.getSerialNumber(),
              atomTwo.getName(), atomTwo.getSerialNumber()));
    }
    Transformation t = new Transformation();
    Vector3D atomOneCoor = atomOne.getCoordinates();
    Vector3D atomTwoCoor = atomTwo.getCoordinates();
    Vector3D translation = atomOneCoor.subtract(atomTwoCoor)
                                      .toUnitVector()
                                      .multiply(-1.0*delta);
    t.addTranslation(translation);
    return t;
  }

  //###########################################################################################\\
  //##############################  VALIDATION & HELPER METHODS  ##############################\\
  //###########################################################################################\\

  /**
  * @param atomOne one of the atoms in the triplet
  * @param atomTwo another of the atoms in the triplet
  * @param atomThree a third of the atoms in the triplet
  * @return true if the three Atoms passed in are a valid angle triplet in this Protein, false
  * otherwise.
  */
  private boolean isValidAngleTriplet(Atom atomOne, Atom atomTwo, Atom atomThree){
    return this.validTriplets.contains(new AngleTriplet(atomOne, atomTwo, atomThree));
  }

  /**
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @return true if the two Atoms passed in are a valid rotatable bond in this Protein, false
  * otherwise.
  */
  public boolean isRotatableBond(Atom atomOne, Atom atomTwo){
    return isRotatableBond(new OrderedBond(atomOne, atomTwo));
  }

  /**
  * @param bond an OrderedBond to check if rotatable
  * @return true if the OrderedBond passed in is a valid rotatable bond in this Protein, false
  * otherwise.
  */
  private boolean isRotatableBond(OrderedBond bond){
    return this.dihedrals.containsKey(bond);
  }

  /**
  * @param atomOne the first atom in the quartet
  * @param atomTwo the second atom in the quartet
  * @param atomThree the third atom in the quartet
  * @param atomFour the fourth atom in the quartet
  * @return true if the four Atoms passed in define a valid dihedral angle, false otherwise.
  */
  private boolean isValidDihedralQuartet(Atom atomOne, Atom atomTwo, Atom atomThree, Atom atomFour){
    return isValidDihedralQuartet(new DihedralQuartet(atomOne, atomTwo, atomThree, atomFour));
  }

  /**
  * @param quartet a DihedralQuartet to check if valid
  * @return true if the DihedralQuartet passed in defines a valid Dihedral angle in this Protein,
  * false otherwise.
  */
  private boolean isValidDihedralQuartet(DihedralQuartet quartet){
    return this.validQuartets.contains(quartet);
  }

  /**
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @return true if the two Atoms passed in are a valid bond in this Protein, false otherwise.
  */
  private boolean isValidBond(Atom atomOne, Atom atomTwo){
    return this.validBonds.contains(new OrderedBond(atomOne, atomTwo));
  }

  /**
  * Get the Dihedral Quartet that the rotatable bond formed by atomOne and atomTwo are a part of.
  * @param atomOne an Atom in the bond
  * @param atomTwo the other Atom in the bond
  * @return a DihedralQuartet, the four atoms that define the dihedral angle of the bond
  * @throws IllegalArgumentException if atomOne and atomTwo are not a rotatable bond
  */
  private DihedralQuartet getDihedralQuartet(Atom atomOne, Atom atomTwo){
    if(!isRotatableBond(atomOne, atomTwo)){
      throw new IllegalArgumentException(
          String.format("%s-%d, %s-%d is not a valid rotatable bond.\n"
              + "Either the atoms are not consecutive or either of their neighbors are missing.",
              atomOne.getName(), atomOne.getSerialNumber(),
              atomTwo.getName(), atomTwo.getSerialNumber()));
    }
    OrderedBond bond = new OrderedBond(atomOne, atomTwo);
    return this.dihedrals.get(bond);
  }

  /**
  * Get a List of the pairs of Atoms that comprise all the Rotatable Bonds in the Protein.
  * @return a List containing all the pairs of Atoms that make up all the Rotatable Bonds in the
  * Protein
  */
  public List<Atom[]> getRotatableBonds(){
    ArrayList<Atom[]> rotatableBonds = new ArrayList<Atom[]>();
    for(OrderedBond bond : this.dihedrals.keySet()){
      rotatableBonds.add(new Atom[]{bond.getAtomOne(), bond.getAtomTwo()});
    }
    return rotatableBonds;
  }

  /**
  * Get a List of the pairs of Atoms that comprise all the valid Dihedral Quartets in the Protein.
  * @return a List containing all the quartets of Atoms that define all the Dihedral Angles in this
  * Protein
  */
  public List<Atom[]> getDihedralQuartets(){
    ArrayList<Atom[]> quartets = new ArrayList<Atom[]>();
    for(DihedralQuartet quartet : this.validQuartets){
      quartets.add(new Atom[]{quartet.getAtomOne(), quartet.getAtomTwo(),
                              quartet.getAtomThree(), quartet.getAtomFour()});
    }
    return quartets;
  }

  /**
  * Get a List of the pairs of Atoms that comprise the backbone Rotatable Bonds in the Protein.
  * @return a List containing all the pairs of Atoms that make up the backbone  Rotatable Bonds in
  * the Protein
  */
  public List<Atom[]> getBackboneRotatableBonds(){
    ArrayList<Atom[]> rotatableBonds = new ArrayList<Atom[]>();
    for(OrderedBond bond : this.dihedrals.keySet()){
      if(bond.isOnBackbone()){
        rotatableBonds.add(new Atom[]{bond.getAtomOne(), bond.getAtomTwo()});
      }
    }
    return rotatableBonds;
  }

  /**
  * Get a List of the pairs of Atoms that comprise all the Bonds in the Protein.
  * @return a List containing all the pairs of Atoms that make up all the Bonds in the Protein
  */
  public List<Atom[]> getBonds(){
    ArrayList<Atom[]> allBonds = new ArrayList<Atom[]>();
    for(OrderedBond bond : this.validBonds){
      allBonds.add(new Atom[]{bond.getAtomOne(), bond.getAtomTwo()});
    }
    return allBonds;
  }

  /**
  * Get a List of the pairs of Atoms that comprise all the backbone Bonds in the Protein.
  * @return a List containing all the pairs of Atoms that make up all the backbone Bonds in the
  * Protein
  */
  public List<Atom[]> getBackboneBonds(){
    ArrayList<Atom[]> allBonds = new ArrayList<Atom[]>();
    for(OrderedBond bond : this.validBonds){
      if(bond.isOnBackbone()){
        allBonds.add(new Atom[]{bond.getAtomOne(), bond.getAtomTwo()});
      }
    }
    return allBonds;
  }

  /**
  * Get a List of the triplets of Atoms that define all angles in the Protein.
  * @return a List containing all the triplets of Atoms that define all the angles in the Protein
  */
  public List<Atom[]> getAngleTriplets(){
    ArrayList<Atom[]> angleTriplets = new ArrayList<Atom[]>();
    for(AngleTriplet triplet : this.validTriplets){
      angleTriplets.add(new Atom[]{triplet.getAtomOne(),
                                   triplet.getAtomTwo(),
                                   triplet.getAtomThree()});
    }
    return angleTriplets;
  }

  /**
  * Get a List of the triplets of Atoms that define all angles in the Protein.
  * @return a List containing all the triplets of Atoms that define all the angles in the Protein
  */
  public List<Atom[]> getBackboneAngleTriplets(){
    ArrayList<Atom[]> angleTriplets = new ArrayList<Atom[]>();
    for(AngleTriplet triplet : this.validTriplets){
      if(triplet.isOnBackbone()){
        angleTriplets.add(new Atom[]{triplet.getAtomOne(),
                                     triplet.getAtomTwo(),
                                     triplet.getAtomThree()});
      }
    }
    return angleTriplets;
  }

  public void writePDB(String fileName){
    try{
      this.protein.writeToFile(fileName);
    }catch(IOException e){
      System.out.println(String.format("Could not write to file %s. File not accessible.",fileName));
    }
  }

  //###########################################################################################\\
  //###################################  ANGLE CALCULATIONS  ##################################\\
  //###########################################################################################\\

  /**
  * Return the angle formed by 3 atoms in this residue.
  * @param atomOne one of the atoms in the angle
  * @param atomTwo the pivot atom
  * @param atomThree the other atom in the angle
  * @return the angle defined by atomOne-atomTwo-atomThree, in degrees
  * @throws IllegalArgumentException if the three atoms are not a valid defined angle triplet
  */
  public double getAngle(Atom atomOne, Atom atomTwo, Atom atomThree){
    if(!isValidAngleTriplet(atomOne, atomTwo, atomThree)){
      throw new IllegalArgumentException(
          String.format("%s-%d, %s-%d, %s-%d is not a valid defined angle triplet.\n",
                        atomOne.getName(), atomOne.getSerialNumber(),
                        atomTwo.getName(), atomTwo.getSerialNumber(),
                        atomThree.getName(), atomThree.getSerialNumber()));
    }
    AngleTriplet triplet = new AngleTriplet(atomOne, atomTwo, atomThree);
    Vector3D one = triplet.getAtomOne().getCoordinates();
    Vector3D two = triplet.getAtomTwo().getCoordinates();
    Vector3D three = triplet.getAtomThree().getCoordinates();
    return one.subtract(two).angle(three.subtract(two));
  }

  /**
  * Return the dihedral angle of a rotatable bond.
  * @param atomOne the first atom of the rotatable bond
  * @param atomTwo the second atom of the rotatable bond
  * @return the degree value of the dihedral angle about this bond
  * @throws IllegalArgumentException if the two atoms do not define a rotatable bond
  */
  public double getDihedralAngle(Atom atomOne, Atom atomTwo){
    if(!isRotatableBond(atomOne, atomTwo)){
      throw new IllegalArgumentException(
          String.format("%s-%d, %s-%d is not a valid rotatable bond.\n"
              + "Either the atoms are not consecutive or either of their neighbors are missing.",
              atomOne.getName(), atomOne.getSerialNumber(),
              atomTwo.getName(), atomTwo.getSerialNumber()));
    }
    DihedralQuartet quartet = getDihedralQuartet(atomOne, atomTwo);
    return getDihedralAngle(quartet.getAtomOne(), quartet.getAtomTwo(),
                            quartet.getAtomThree(), quartet.getAtomFour());
  }

  /**
  * Return the dihedral angle formed by the four atoms.
  * @param atomOne the first atom of the dihedral
  * @param atomTwo the second atom of the dihedral
  * @param atomThree the third atom of the dihedral
  * @param atomFour the fourth atom of the dihedral
  * @return the degree value of the dihedral angle about this bond, or 1000 if any of the atoms
  * required to calculate the angle are missing.
  * @throws IllegalArgumentException if the four atoms do not define a valid dihedral angle
  */
  public double getDihedralAngle(Atom atomOne, Atom atomTwo, Atom atomThree, Atom atomFour){
    DihedralQuartet quartet = new DihedralQuartet(atomOne, atomTwo, atomThree, atomFour);
    if(!isValidDihedralQuartet(quartet)){
      throw new IllegalArgumentException(
        String.format("%s-%s-%s-%s is not a valid dihedral quartet.\n Either the Atoms are not"
                    + "consecutive or one of them is missing from the Protein.",
                    atomOne.getName(), atomTwo.getName(),
                    atomThree.getName(), atomFour.getName()));
    }

    Vector3D one = quartet.getAtomOne().getCoordinates();
    Vector3D two = quartet.getAtomTwo().getCoordinates();
    Vector3D three = quartet.getAtomThree().getCoordinates();
    Vector3D four = quartet.getAtomFour().getCoordinates();
    return -1 * Vector3D.calculateDihedralAngle(one, two, three, four);
  }

  //###########################################################################################\\
  //##############################  BOOKKEEPING & INITIALIZATION  #############################\\
  //###########################################################################################\\

  private HashMap<Atom,Integer> constructAtomIndices(){
    HashMap<Atom,Integer> atomIndices = new HashMap<Atom,Integer>();
    int index = 0;
    for(PolypeptideChain chain : this.protein){
      for(Residue res : chain){
        List<Atom> atomsInOrder = res.getAtomsInOrder();
        for(Atom atom : atomsInOrder){
          atomIndices.put(atom, index);
          index++;
        }
      }
    }
    return atomIndices;
  }

  /**
  * Build all of the valid Angle Triplets: sets of three consecutive atoms that define an angle.
  * This set can be used to efficiently validate arguments passed in when manipulating angles.
  * Importantly, since order matters when calculating angles, it also imposes an ordering the atoms
  * in the triplet - from N-Terminus to C-Terminus and within the side chain, from CA on back.
  * @return a {@code HashSet<AngleTriplet>} containing the valid angle triplets for the protein.
  */
  private HashSet<AngleTriplet> buildAngleTriplets(){
    HashSet<AngleTriplet> triplets = new HashSet<AngleTriplet>();
    for(PolypeptideChain chain : this.protein){
      for(Residue res : chain){
        ArrayList<Atom[]> resTriplets = res.getAngleTriplets();
        for(Atom[] resTriplet : resTriplets){
          triplets.add(new AngleTriplet(resTriplet[0], resTriplet[1], resTriplet[2]));
        }
        if(chain.contains(res.getResidueID()-1)){
          Residue prev = chain.getResidue(res.getResidueID()-1);
          if(prev.contains("C") && res.contains("N") && res.contains("CA")){
            triplets.add(new AngleTriplet(prev.getAtom("C"), res.getAtom("N"), res.getAtom("CA")));
          }
        }
        if(chain.contains(res.getResidueID()+1)){
          Residue next = chain.getResidue(res.getResidueID()+1);
          if(res.contains("CA") && res.contains("C") && next.contains("N")){
            triplets.add(new AngleTriplet(res.getAtom("CA"), res.getAtom("C"), next.getAtom("N")));
          }
          if(res.contains("C") && next.contains("N") && next.contains("H")){
            triplets.add(new AngleTriplet(res.getAtom("C"), next.getAtom("N"), next.getAtom("H")));
          }
        }
      }
    }
    return triplets;
  }

  /**
  * Build all of the valid Dihedral Quartets: sets of four consecutive atoms that define a dihedral
  * angle.
  * This set can be used to efficiently validate arguments passed in when manipulating angles.
  * Importantly, since order matters when calculating angles, it also imposes an ordering the atoms
  * in the triplet - from N-Terminus to C-Terminus and within the side chain, from CA on back.
  * @return a {@code HashMap<OrderedBond, DihedralQuartet>}, a mapping from the rotatable bonds
  * of the Protein to the DihedralQuartets containing the atoms definind the dihedral angles of
  * those dihedral angles.
  */
  private HashMap<OrderedBond, DihedralQuartet> buildDihedralQuartets(){
    HashMap<OrderedBond, DihedralQuartet> dihedrals = new HashMap<OrderedBond, DihedralQuartet>();
    for(PolypeptideChain chain : this.protein){
      for(Residue res : chain){
        ArrayList<Atom[]> resQuartets = res.getDihedralQuartets();
        for(Atom[] resQuartet : resQuartets){
          DihedralQuartet quartet = new DihedralQuartet(resQuartet[0], resQuartet[1],
                                                        resQuartet[2], resQuartet[3]);
          OrderedBond rotatableBond = new OrderedBond(quartet.getAtomTwo(), quartet.getAtomThree());
          dihedrals.put(rotatableBond, quartet);
          this.validQuartets.add(quartet);
        }
        if(chain.contains(res.getResidueID()-1)){
          Residue prev = chain.getResidue(res.getResidueID()-1);
          if(prev.contains("CA") && prev.contains("C") && res.contains("N") && res.contains("CA")){
            DihedralQuartet quartet = new DihedralQuartet(prev.getAtom("CA"), prev.getAtom("C"),
                                                          res.getAtom("N"), res.getAtom("CA"));
            OrderedBond rotatableBond = new OrderedBond(prev.getAtom("C"), res.getAtom("N"));
            dihedrals.put(rotatableBond, quartet);
            this.validQuartets.add(quartet);
          }
          if(prev.contains("C") && res.contains("N") && res.contains("CA") && res.contains("C")){
            DihedralQuartet quartet = new DihedralQuartet(prev.getAtom("C"), res.getAtom("N"),
                                                          res.getAtom("CA"), res.getAtom("C"));
            OrderedBond rotatableBond = new OrderedBond(res.getAtom("N"), res.getAtom("CA"));
            dihedrals.put(rotatableBond, quartet);
            this.validQuartets.add(quartet);
          }
          // Side Chain Rotation
          if(prev.contains("C") && res.contains("N") && res.contains("CA") && res.contains("CB")){
            DihedralQuartet quartet = new DihedralQuartet(prev.getAtom("C"), res.getAtom("N"),
                                                          res.getAtom("CA"), res.getAtom("CB"));
            this.validQuartets.add(quartet);
          }
        }
        if(chain.contains(res.getResidueID()+1)){
          Residue next = chain.getResidue(res.getResidueID()+1);
          if(res.contains("N") && res.contains("CA") && res.contains("C") && next.contains("N")){
            DihedralQuartet quartet = new DihedralQuartet(res.getAtom("N"), res.getAtom("CA"),
                                                          res.getAtom("C"), next.getAtom("N"));
            OrderedBond rotatableBond = new OrderedBond(res.getAtom("CA"), res.getAtom("C"));
            dihedrals.put(rotatableBond, quartet);
            this.validQuartets.add(quartet);
          }
        }
        if(res.contains("N") && res.contains("CA") && res.contains("C") && res.contains("O")){
          DihedralQuartet quartet = new DihedralQuartet(res.getAtom("N"), res.getAtom("CA"),
                                                        res.getAtom("C"), res.getAtom("O"));
          this.validQuartets.add(quartet);
        }
      }
    }
    return dihedrals;
  }

  /**
  * Build all of the valid Bonds: sets of two consecutive atoms that define a Bond.
  * This set can be used to efficiently validate arguments passed in when manipulating bonds.
  * Importantly, since order matters when calculating angles, it also imposes an ordering the atoms
  * in the triplet - from N-Terminus to C-Terminus and within the side chain, from CA on back.
  * @return a {@code HashSet<OrderedBond>} containing the valid Bond Atom pairs for the protein.
  */
  private HashSet<OrderedBond> buildBonds(){
    HashSet<OrderedBond> validBonds = new HashSet<OrderedBond>();
    for(PolypeptideChain chain : this.protein){
      for(Bond bond : chain.getBonds()){
        validBonds.add(new OrderedBond(bond.getAtomOne(), bond.getAtomTwo()));
      }
    }
    for(Bond disulfideBond : this.protein.getDisulfideBonds()){
      validBonds.add(new OrderedBond(disulfideBond.getAtomOne(), disulfideBond.getAtomTwo()));
    }
    return validBonds;
  }

  //###########################################################################################\\
  //######################################  Innerclasses  #####################################\\
  //###########################################################################################\\

  /**
  * Inner Class OrderedBond is a pair of atoms: consecutive, covalently bonded atoms
  * that form a single bond. This class imposes an ordering on the atoms in the pair, from
  * N-terminus to C-terminus and within residues from CA on back.
  */
  class OrderedBond extends OrderedAtoms{
    private boolean onBackbone;

    public OrderedBond(Atom one, Atom two){
      super(one, two);
      this.onBackbone = false;
      if(getAtomTwo().getName().equals("N")){
        onBackbone = true;
      }
      if(getAtomTwo().getName().equals("CA")){
        onBackbone = true;
      }
      if(getAtomTwo().getName().equals("C")){
        onBackbone = true;
      }
    }

    public Atom getAtomOne(){
      return super.getAtom(0);
    }

    public Atom getAtomTwo(){
      return super.getAtom(1);
    }

    public boolean isOnBackbone(){
      return this.onBackbone;
    }
  }

  /**
  * Inner Class AngleTriplet is a triplet of atoms: consecutive, covalently bonded atoms
  * that form a simple angle. This class imposes an ordering on the atoms in the triple, from
  * N-terminus to C-terminus and within residues from CA on back.
  */
  class AngleTriplet extends OrderedAtoms{
    private boolean onBackbone;

    public AngleTriplet(Atom one, Atom two, Atom three){
      super(one, two, three);
      this.onBackbone = false;
      if(getAtomTwo().getName().equals("C") && getAtomThree().getName().equals("N")){
        onBackbone = true;
      }
      if(getAtomTwo().getName().equals("N") && getAtomThree().getName().equals("CA")){
        onBackbone = true;
      }
      if(getAtomTwo().getName().equals("CA") && getAtomThree().getName().equals("C")){
        onBackbone = true;
      }
    }

    public Atom getAtomOne(){
      return super.getAtom(0);
    }

    public Atom getAtomTwo(){
      return super.getAtom(1);
    }

    public Atom getAtomThree(){
      return super.getAtom(2);
    }

    public boolean isOnBackbone(){
      return this.onBackbone;
    }
  }

  /**
  * Inner Class DihedralQuartet is a quartet of atoms: consecutive, covalently bonded atoms
  * that form a dihedral angle. This class imposes an ordering on the atoms in the quartet, from
  * N-terminus to C-terminus and within residues from CA on back.
  */
  class DihedralQuartet extends OrderedAtoms{
    private boolean onBackbone;

    public DihedralQuartet(Atom one, Atom two, Atom three, Atom four){
      super(one, two, three, four);
      this.onBackbone = false;
      if(getAtomTwo().getName().equals("C") && getAtomThree().getName().equals("N")){
        onBackbone = true;
      }
      if(getAtomTwo().getName().equals("N") && getAtomThree().getName().equals("CA")
          && !getAtomFour().getName().equals("CB")){
        onBackbone = true;
      }
      if(getAtomTwo().getName().equals("CA") && getAtomThree().getName().equals("C")){
        onBackbone = true;
      }
    }

    public final Atom getAtomOne(){
      return super.getAtom(0);
    }

    public final Atom getAtomTwo(){
      return super.getAtom(1);
    }

    public final Atom getAtomThree(){
      return super.getAtom(2);
    }

    public final Atom getAtomFour(){
      return super.getAtom(3);
    }

    public boolean isOnBackbone(){
      return this.onBackbone;
    }
  }

  /**
  * Private Inner Class OrderedAtoms is a set of atoms: consecutive, covalently bonded atoms. This
  * class imposes an ordering on the atoms in the set, from N-terminus to C-terminus and within
  * residues from CA on back.
  */
  private abstract class OrderedAtoms{
    private ArrayList<Atom> orderedAtoms;

    public OrderedAtoms(Atom... atoms){
      orderedAtoms = new ArrayList<Atom>();
      HashSet<Atom> set = new HashSet<Atom>();
      for(Atom atom : atoms){
        set.add(atom);
        this.orderedAtoms.add(atom);
      }
      if(set.size() != atoms.length){
        String atomsString = "";
        for(Atom atom : atoms){
          atomsString += atom.getName()+ " ";
        }
        throw new IllegalArgumentException(
            String.format("Recieved atoms not all unique. Atoms must be different.\n"
                        + "Atoms recieved: %s", atomsString));
      }
      Collections.sort(orderedAtoms, new AtomComparator());
    }

    public Atom getAtom(int index){
      return this.orderedAtoms.get(index);
    }

    public int getNumAtoms(){
      return this.orderedAtoms.size();
    }

    private boolean contains(Atom atom){
      return orderedAtoms.contains(atom);
    }

    /**
    * Returns true if a transformation applied to this set of atoms would modify the backbone of
    * the protein.
    */
    public abstract boolean isOnBackbone();

    @Override
    public int hashCode(){
      long hash = 0;
      for(Atom atom : this.orderedAtoms){
        hash += (long) atom.hashCode();
      }
      hash %= Integer.MAX_VALUE;
      return (int) hash;
    }

    @Override
    public boolean equals(Object obj){
      if(obj instanceof OrderedAtoms){
        OrderedAtoms other = (OrderedAtoms)obj;
        if(this.getNumAtoms() != other.getNumAtoms()){
          return false;
        }
        boolean containsAllSameAtoms = true;
        for(int i = 0; i < this.getNumAtoms(); i++){
          containsAllSameAtoms = containsAllSameAtoms && this.getAtom(i).equals(other.getAtom(i));
        }
        return containsAllSameAtoms;
      }
      return false;
    }

  }

  private class AtomComparator implements Comparator<Atom>{
    public int compare(Atom a, Atom b){
      int indexA = atomIndices.get(a);
      int indexB = atomIndices.get(b);
      return indexA - indexB;
    }
  }
}
