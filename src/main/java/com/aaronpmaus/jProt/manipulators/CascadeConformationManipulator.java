package com.aaronpmaus.jProt.manipulators;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.transformations.Transformation;
import com.aaronpmaus.jProt.protein.*;

import java.util.Collection;

/**
* A CascadeConformationManipulator is a ConformationManipulator where the manipulation on a bond,
* angle, or dihedral angle affects all atoms after that bond, angle, or dihedral angle.  The
* direction that the manipulation is applied in is the N to C terminus direction. For example,  if a
* bond is being manipulated, a cascade manipulation would translate the second atom in that bond and
* all atoms after it.
* <p>
* For example, if the second atom is the C of Residue 2, all atoms from it to
* the C-Terminus of that chain will be modified. If the second atom is the CB of Residue 2, all
* atoms from it to the end of the side chain of Residue 2 will be modified, but no atoms on the
* backbone or in following residues will be modified.
* @version 0.7.0
* @since 0.7.0
*/
public class CascadeConformationManipulator extends ConformationManipulator {
  // a directed graph of all the atoms in this protein. Not necessarily connected if the chains in
  // this protein are not covalently connected. The edges are directed from N-terminus to
  // C-terminus. This gives CascadeConformationManipulator an easy way of determining all the atoms to
  // apply a cascate transformation to: a depth first search from the initial atom on.
  private Graph<Atom> atomsGraph;
  // A reference to the Protein that this ConformationManipulator manipulator manipulates.
  private Protein protein;


  public CascadeConformationManipulator(Protein protein){
    super(protein);
    this.protein = protein;
  }

  @Override
  public void modifyAngle(Atom atomOne, Atom atomTwo, Atom atomThree, double degrees){

    Transformation t = buildAngleTransformation(atomOne, atomTwo, atomThree, degrees);

    ConformationManipulator.AngleTriplet triplet =
        new ConformationManipulator.AngleTriplet(atomOne, atomTwo, atomThree);

    String chainID = protein.getChainID(triplet.getAtomTwo());
    PolypeptideChain chain = protein.getChain(chainID);
    Residue res = chain.getResidue(triplet.getAtomThree());

    applyTransformation(t, triplet.isOnBackbone(), chain, res, triplet.getAtomThree());
  }

  @Override
  public void modifyDihedralAngle(Atom atomOne, Atom atomTwo,
                              Atom atomThree, Atom atomFour, double degrees){

    Transformation t = buildDihedralRotationTransformation(atomOne, atomTwo,
                                                           atomThree, atomFour, degrees);

    ConformationManipulator.DihedralQuartet quartet =
        new ConformationManipulator.DihedralQuartet(atomOne, atomTwo, atomThree, atomFour);

    String chainID = protein.getChainID(quartet.getAtomThree());
    PolypeptideChain chain = protein.getChain(chainID);
    Residue res = chain.getResidue(quartet.getAtomThree());

    // Apply the transformation to all atoms after bond.getAtomTwo() within its residue
    // If the bond is on the backbone, apply the transformation to all following residues.
    if(quartet.getAtomOne().getName().equals("C")
        && quartet.getAtomTwo().getName().equals("N")
        && quartet.getAtomThree().getName().equals("CA")
        && quartet.getAtomFour().getName().equals("CB")){
      // if this is a side chain rotation, only apply the Transformation to the atoms in the side
      // chain.
      res = chain.getResidue(quartet.getAtomFour());
      applyTransformation(t, quartet.isOnBackbone(), chain, res, quartet.getAtomFour());
    } else if(quartet.getAtomOne().getName().equals("C")
              && quartet.getAtomTwo().getName().equals("N")
              && quartet.getAtomThree().getName().equals("CA")
              && quartet.getAtomFour().getName().equals("C")){
      // If this is a Phi Angle rotation, apply the Transformation to all atoms after N. This will
      // include the H bonded to the N.
      res = chain.getResidue(quartet.getAtomTwo());
      applyTransformation(t, quartet.isOnBackbone(), chain, res, quartet.getAtomTwo());
    } else {
      applyTransformation(t, quartet.isOnBackbone(), chain, res, quartet.getAtomThree());
    }
  }

  @Override
  public void modifyBondLength(Atom atomOne, Atom atomTwo, double delta){

    ConformationManipulator.OrderedBond bond =
        new ConformationManipulator.OrderedBond(atomOne, atomTwo);

    Transformation t = buildBondManipulationTransformation(bond.getAtomOne(),
                                                           bond.getAtomTwo(), delta);
    // The rotatable bonds are those that form the axis of rotation of a defined dihedral angle.
    // This means that there are atoms after this bond that will be affected by any changes to
    // it. Non rotatable bonds don't have atoms after them (such as the C-O bond, no atoms after
    // O), or they are in parts of structures where dihedral rotations can't occur such as within
    // rings. In those parts of the structures, follow on translations from Bond Length
    // manipulations aren't easily definable so don't cascade the manipulation.
    // By definition, rotatable bonds are missing important bonds on the backbone where the cascade
    // should occur: all N-CA and CA-C bonds where there is not a residue immediately before or
    // following that residue. This includes the N and C termini and all residues before or after
    // a gap. The bonds that missing are all backbone bonds, so apply the cascade to all rotatable
    // bonds and all bonds on the backbone.
    if(isRotatableBond(atomOne, atomTwo) || bond.isOnBackbone()){
      String chainID = protein.getChainID(bond.getAtomTwo());
      PolypeptideChain chain = protein.getChain(chainID);
      Residue res = chain.getResidue(atomTwo);
      //  Apply the transformation to all atoms after bond.getAtomTwo() within its residue
      //  If the bond is on the backbone, apply the transformation to all following residues.
      applyTransformation(t, bond.isOnBackbone(), chain, res, bond.getAtomTwo());

    } else {
      bond.getAtomTwo().applyTransformation(t);
    }
  }

  /**
  * Apply the Transformation to all atoms after and including atom within its residue. If the atom
  * is on the backbone, apply the transformation to all residues after res in the chain.
  * @param t the Transformation to apply
  * @param isOnBackbone true if the Transformation will affect backbone atoms, false otherwise
  * @param chain the PolypeptideChain that this Tranformation is manipulating
  * @param res the Residue that contains atom, the first atom that needs to be Transformed
  * @param atom the first Atom that needs to be Transformed: the pivot in a three atom angle, the
  * second atom in the rotatable bond of a dihedral quartet, or the second atom in a bond.
  */
  private void applyTransformation(Transformation t, boolean isOnBackbone,
                                   PolypeptideChain chain, Residue res, Atom atom){
    int resID = res.getResidueID();
    // get all atoms after atom in the Residue and apply the transformation to them
    Collection<Atom> atomsToTransform = res.getAtomsAfterAndIncluding(atom);
    //System.out.printf("There are %d atoms to transform.\n", atomsToTransform.size());
    for(Atom atomToTransform : res.getAtomsAfterAndIncluding(atom)){
      //System.out.printf("   Applying Transformation to %s.\n", atomToTransform.getName());
      atomToTransform.applyTransformation(t);
    }
    if(isOnBackbone){
      for(Residue residueToTransform : chain.getResiduesToEnd(resID + 1)){
        residueToTransform.applyTransformation(t);
      }
    }
  }
}
