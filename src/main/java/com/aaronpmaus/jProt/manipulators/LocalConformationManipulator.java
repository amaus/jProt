package com.aaronpmaus.jProt.manipulators;

import com.aaronpmaus.jProt.protein.*;
import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.transformations.Transformation;

/**
* A LocalConformationManipulator is a ConformationManipulator where the manipulation on a bond,
* angle, or dihedral angle only affects the atoms of that bond, angle, or dihedral angle.   For
* example, if a bond is being manipulated, a local manipulation would alter only the atoms in that
* bond.
* @version 0.7.0
* @since 0.7.0
*/
public class LocalConformationManipulator extends ConformationManipulator {
  // A reference to the Protein that this ConformationManipulator manipulator manipulates.
  private Protein protein;

  public LocalConformationManipulator(Protein protein){
    super(protein);
    this.protein = protein;
  }

  @Override
  public void modifyAngle(Atom atomOne, Atom atomTwo, Atom atomThree, double degrees){
    Transformation t = buildAngleTransformation(atomOne, atomTwo, atomThree, degrees);
    ConformationManipulator.AngleTriplet triplet =
        new ConformationManipulator.AngleTriplet(atomOne, atomTwo, atomThree);
    triplet.getAtomThree().applyTransformation(t);
  }

  @Override
  public void modifyDihedralAngle(Atom atomOne, Atom atomTwo,
                                  Atom atomThree, Atom atomFour, double degrees){
    Transformation t = buildDihedralRotationTransformation(atomOne, atomTwo,
                                                           atomThree, atomFour, degrees);
    ConformationManipulator.DihedralQuartet quartet =
        new ConformationManipulator.DihedralQuartet(atomOne, atomTwo, atomThree, atomFour);
    quartet.getAtomFour().applyTransformation(t);
  }

  @Override
  public void modifyBondLength(Atom atomOne, Atom atomTwo, double delta){
    Transformation t = buildBondManipulationTransformation(atomOne, atomTwo, delta);
    atomTwo.applyTransformation(t);
  }
}
