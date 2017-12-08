package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.sequence.ProteinSequence;
import com.aaronpmaus.jMath.transformations.Transformation;
import com.aaronpmaus.jMath.linearAlgebra.*;

import java.math.BigDecimal;
import java.math.MathContext;


public class VirtualRibosome {
  /**
  * Construct a Protein with a single chain. This chain is unfolded, linear, all phi, psi, and
  * omega angles are either 180 or -180.
  * @param sequence the sequence of the protein to construct
  * @return a Protein with a single chain built out of the residues in that sequence
  */
  public static Protein synthesizeProtein(ProteinSequence sequence){
    Protein prot = new Protein("Synthaprot");
    prot.addChain(synthesizeChain(sequence, "A"));
    return prot;
  }

  private static PolypeptideChain synthesizeChain(ProteinSequence sequence, String chainID){
    PolypeptideChain chain = new PolypeptideChain(chainID);
    Residue prevResidue = null;
    Residue newResidue = null;
    int resID = 1;
    // For every residue, rotate it and translate it so that its peptide bond with the previous
    // residue is physically reasonable
    for(Character resString : sequence){
      if(resID == 1){
        prevResidue = new Residue(resString, resID);
        //geometryExperiments(prevResidue);
        chain.addResidue(prevResidue);
      } else {
        newResidue = new Residue(resString, resID);
        moveResidueIntoPlace(prevResidue, newResidue);
        chain.addResidue(newResidue);
        chain.setOmegaAngle(resID, 180);
        prevResidue = newResidue;
      }
      resID++;
    }
    return chain;
  }

  /*
  * Move the next Residue into place next to the previous Residue.
  *
  * After the Transformations are applied, the O-C-N angle is 125, N is 1.32 Angstroms away from
  * C, and the C-N-CA angle is 123 degrees.
  */
  private static void moveResidueIntoPlace(Residue prevResidue, Residue nextResidue){
    Atom c = prevResidue.getAtom("C");
    Atom o = prevResidue.getAtom("O");
    Atom ca = prevResidue.getAtom("CA");
    Atom n = nextResidue.getAtom("N");
    Atom ca2 = nextResidue.getAtom("CA");

    Vector3D c_o = o.getCoordinates().subtract(c.getCoordinates());
    Vector3D c_ca = ca.getCoordinates().subtract(c.getCoordinates());
    Vector3D normal = c_ca.crossProduct(c_o);

    // Construct and apply transformations such that the N of the next Residue is in a 125 degree
    // angle O-C-N and is 1.32 Angstroms away from C.
    Transformation rotateCO = new Transformation();
    rotateCO.addRotationAboutAxis(c.getCoordinates(), normal, 125);
    Vector3D nitrogenOrientation = new Vector3D(o.getCoordinates());
    // Rotate the coordinates of oxygen about the axis passing through C in the
    // direction of the normal 125 degrees.
    nitrogenOrientation.applyTransformation(rotateCO);
    // move nitrogen's location out from C so that it is 1.32 Angstroms away
    Vector3D c_n = nitrogenOrientation.subtract(c.getCoordinates());
    c_n = c_n.toUnitVector().multiply(1.32);
    Vector3D nitrogenLocation = c.getCoordinates().add(c_n);

    // translate the next residue so that the nitrogen is in the correct place
    Vector3D translation = nitrogenLocation.subtract(n.getCoordinates());
    Transformation translateResIntoPlace = new Transformation();
    translateResIntoPlace.addTranslation(translation);
    nextResidue.applyTransformation(translateResIntoPlace);

    // Rotate the next Residue so that the CA-N-C angle is 123 degrees
    Vector3D n_c = c.getCoordinates().subtract(n.getCoordinates());
    Vector3D n_ca2 = ca2.getCoordinates().subtract(n.getCoordinates());
    normal = n_c.crossProduct(n_ca2);
    double angle = n_c.angle(n_ca2);
    // calculate the angle to rotate
    double angleToRotate = 123 - angle;
    // construct a rotation about the CA-N-C normal for the needed amount of degrees
    // and rotate the residue that much
    Transformation rotateAboutNitrogen = new Transformation();
    rotateAboutNitrogen.addRotationAboutAxis(n.getCoordinates(), normal, angleToRotate);
    nextResidue.applyTransformation(rotateAboutNitrogen);
  }
}
