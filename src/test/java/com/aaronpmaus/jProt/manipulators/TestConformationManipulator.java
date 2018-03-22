package com.aaronpmaus.jProt;

import com.aaronpmaus.jProt.protein.*;
import com.aaronpmaus.jProt.manipulators.*;
import com.aaronpmaus.jProt.sequence.*;
import com.aaronpmaus.jMath.linearAlgebra.*;

import static org.junit.Assert.*;
import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Rule;
import org.junit.rules.ExpectedException;

import java.io.IOException;
import java.util.List;

/*
 * @Test flags a method as a test method.
 * @Before indicates that a method will be run before every
 *test method is run.
 * @BeforeClass indicates that a method will be run once before
 *  any of the other methods in the test suite are run.
 * @After indicates that a method will be run after every
 *  test method is run.
 * @AfterClass indicates that a method will be run once after
 *  all the other methods in the test suite finish..
*/
public class TestConformationManipulator{
  Protein prot;
  PolypeptideChain chain;
  CascadeConformationManipulator cascadeManip;

  @Rule
  public final ExpectedException exception = ExpectedException.none();

  @Before
  public void setup() throws IOException{
    prot = VirtualRibosome.synthesizeProtein(new ProteinSequence("IAMSTARSTFF"),
                                                     "strstf");
    chain = prot.getChain("A");
    cascadeManip = new CascadeConformationManipulator(prot);


    Residue phe = chain.getResidue(10);
    cascadeManip.rotateAboutBond(phe.getAtom("CA"), phe.getAtom("CB"), 30);
    cascadeManip.rotateAboutBond(phe.getAtom("CB"), phe.getAtom("CG"), -90);
    Residue thr = chain.getResidue(5);
    cascadeManip.rotateAboutBond(thr.getAtom("CA"), thr.getAtom("CB"), -90);
    thr = chain.getResidue(9);
    cascadeManip.rotateAboutBond(thr.getAtom("CA"), thr.getAtom("CB"), -160);

  }

  @Test
  public void testCascadeSetBackboneDihedralAngles(){

    assertTrue(Math.abs(Math.abs(chain.getOmegaAngle(2)) - 180.0) < 0.00000001);
    cascadeManip.setOmegaAngle("A", 2, 30.0);
    assertTrue(Math.abs(chain.getOmegaAngle(2) - 30.0) < 0.00000001);

    assertTrue(Math.abs(Math.abs(chain.getPhiAngle(2)) - 180.0) < 0.00000001);
    cascadeManip.setPhiAngle("A", 2, 60.0);
    assertTrue(Math.abs(chain.getPhiAngle(2) - 60.0) < 0.00000001);

    assertTrue(Math.abs(Math.abs(chain.getPsiAngle(2)) - 180.0) < 0.00000001);
    cascadeManip.setPsiAngle("A", 2, 90.0);
    assertTrue(Math.abs(chain.getPsiAngle(2) - 90.0) < 0.00000001);

    cascadeManip.setPsiAngle("A", 2, 180.0);
    assertTrue(Math.abs(Math.abs(chain.getPsiAngle(2)) - 180.0) < 0.00000001);
    cascadeManip.setPsiAngle("A", 2, -90.0);
    assertTrue(Math.abs(chain.getPsiAngle(2) - -90.0) < 0.00000001);

    //cascadeManip.setOxygenDihedralAngle(2, 180.0);
    //assertTrue(Math.abs(Math.abs(chain.getOxygenDihedralAngle(2)) - 180.0) < 0.00000001);
    //cascadeManip.setOxygenDihedralAngle(2, -90.0);
    //assertTrue(Math.abs(chain.getOxygenDihedralAngle(2) - -90.0) < 0.00000001);

  }

  @Test
  public void testCascadeSetBondLength(){

    // ensure that the cascade manipulation works at the beginning of the chain
    Atom atomOne = chain.getResidue(1).getAtom("N");
    Atom atomTwo = chain.getResidue(1).getAtom("CA");
    Atom atomThree = chain.getResidue(1).getAtom("C");

    double nextBondLength = atomTwo.distance(atomThree);

    cascadeManip.setBondLength(atomOne, atomTwo, 1.00);

    assertTrue(Math.abs(atomOne.distance(atomTwo) - 1.00) < 0.000000001);
    assertTrue(Math.abs(nextBondLength - atomTwo.distance(atomThree)) < 0.000000001);

    // ensure that the cascade manipulation works at the end of the chain
    atomOne = chain.getResidue(11).getAtom("CA");
    atomTwo = chain.getResidue(11).getAtom("C");
    atomThree = chain.getResidue(11).getAtom("O");

    nextBondLength = atomTwo.distance(atomThree);

    cascadeManip.setBondLength(atomOne, atomTwo, 1.00);

    assertTrue(Math.abs(atomOne.distance(atomTwo) - 1.00) < 0.000000001);
    assertTrue(Math.abs(nextBondLength - atomTwo.distance(atomThree)) < 0.000000001);
  }

  @Test
  public void testCascadeSetAngle(){
    // Cascade angle manipulations should affect every atom after the 3rd atom that defines
    // the angle.
    Atom atomOne = chain.getResidue(1).getAtom("N");
    Atom atomTwo = chain.getResidue(1).getAtom("CA");
    Atom atomThree = chain.getResidue(1).getAtom("C");
    Atom atomFour = chain.getResidue(1).getAtom("CB");

    // save to ensure that the side chain angle isn't changed by altering the backbone angle.
    double beforeAngle = cascadeManip.getAngle(atomOne, atomTwo, atomFour);

    cascadeManip.setAngle(atomOne, atomTwo, atomThree, 30.0);
    assertTrue(Math.abs(cascadeManip.getAngle(atomOne, atomTwo, atomThree) - 30.0) < 0.000000001);

    // ensure that the side chain angle isn't modified
    double afterAngle = cascadeManip.getAngle(atomOne, atomTwo, atomFour);
    assertTrue(Math.abs(beforeAngle - afterAngle) < 0.000000001);

    // check that for the last residue, if the N-CA-C angle is modified, that both the OXT and O
    // atoms are transformed, that the CA-C-OXT and the CA-C-O angles are invariant to the
    // manipulation.
    atomOne = chain.getResidue(11).getAtom("N");
    atomTwo = chain.getResidue(11).getAtom("CA");
    atomThree = chain.getResidue(11).getAtom("C");
    atomFour = chain.getResidue(11).getAtom("O");
    Atom atomFive = chain.getResidue(11).getAtom("OXT");

    // save to ensure that the following angles aren't changed by altering the backbone angle.
    double beforeAngle1 = cascadeManip.getAngle(atomTwo, atomThree, atomFour);
    double beforeAngle2 = cascadeManip.getAngle(atomTwo, atomThree, atomFive);

    cascadeManip.setAngle(atomOne, atomTwo, atomThree, 45.0);
    assertTrue(Math.abs(cascadeManip.getAngle(atomOne, atomTwo, atomThree) - 45.0) < 0.000000001);

    double afterAngle1 = cascadeManip.getAngle(atomTwo, atomThree, atomFour);
    double afterAngle2 = cascadeManip.getAngle(atomTwo, atomThree, atomFive);
    assertTrue(Math.abs(beforeAngle1 - afterAngle1) < 0.000000001);
    assertTrue(Math.abs(beforeAngle2 - afterAngle2) < 0.000000001);
  }

  @Test
  public void testCascadeSetAngleDoesNotAlterDihedralAngles(){
    List<Double> beforeAngles = chain.getBackboneDihedralAngles();
    Atom atomOne = chain.getResidue(1).getAtom("N");
    Atom atomTwo = chain.getResidue(1).getAtom("CA");
    Atom atomThree = chain.getResidue(1).getAtom("C");

    cascadeManip.setAngle(atomOne, atomTwo, atomThree, 30.0);

    List<Double> afterAngles = chain.getBackboneDihedralAngles();

    for(int i = 0; i < beforeAngles.size(); i++){
      assertTrue(Math.abs(beforeAngles.get(i) - afterAngles.get(i)) < 0.000000001);
    }
  }

  @Test
  public void testCascadeSideChainRotation(){
    Atom atomOne = chain.getResidue(1).getAtom("C");
    Atom atomTwo = chain.getResidue(2).getAtom("N");
    Atom atomThree = chain.getResidue(2).getAtom("CA");
    Atom atomFour = chain.getResidue(2).getAtom("CB");

    double beforePsi1 = chain.getPsiAngle(1);
    double beforeOmega2 = chain.getOmegaAngle(2);
    double beforePhi2 = chain.getPhiAngle(2);
    double beforePsi2 = chain.getPsiAngle(2);

    cascadeManip.setDihedralAngle(atomOne, atomTwo, atomThree, atomFour, 15.0);
    assertTrue(Math.abs(
        cascadeManip.getDihedralAngle(atomOne, atomTwo, atomThree, atomFour) - 15.0) < 0.000000001);

    double afterPsi1 = chain.getPsiAngle(1);
    double afterOmega2 = chain.getOmegaAngle(2);
    double afterPhi2 = chain.getPhiAngle(2);
    double afterPsi2 = chain.getPsiAngle(2);
    assertTrue(Math.abs(beforePsi1 - afterPsi1) < 0.000000001);
    assertTrue(Math.abs(beforeOmega2 - afterOmega2) < 0.000000001);
    assertTrue(Math.abs(beforePhi2 - afterPhi2) < 0.000000001);
    assertTrue(Math.abs(beforePsi2 - afterPsi2) < 0.000000001);
  }


  private static double getAngle(Atom atom1, Atom atom2, Atom atom3){
    Vector3D one = atom1.getCoordinates();
    Vector3D two = atom2.getCoordinates();
    Vector3D three = atom3.getCoordinates();

    return one.subtract(two).angle(three.subtract(two));
  }
}
