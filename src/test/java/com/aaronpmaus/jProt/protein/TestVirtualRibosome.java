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
public class TestVirtualRibosome{
  Protein prot;
  PolypeptideChain chain;

  @Before
  public void setup() throws IOException{
    ProteinSequence seq = new ProteinSequence("IAMSTARSTFF");
    //ProteinSequence seq = new ProteinSequence("GGGGGG");
    prot = VirtualRibosome.synthesizeProtein(seq, "strstf");
    chain = prot.getChain("A");

    // Pose pretty for picture
    CascadeConformationManipulator manip = new CascadeConformationManipulator(prot);
    Residue phe = chain.getResidue(10);
    manip.rotateAboutBond(phe.getAtom("CA"), phe.getAtom("CB"), 30);
    manip.rotateAboutBond(phe.getAtom("CB"), phe.getAtom("CG"), -90);
    Residue thr = chain.getResidue(5);
    manip.rotateAboutBond(thr.getAtom("CA"), thr.getAtom("CB"), -90);
    thr = chain.getResidue(9);
    manip.rotateAboutBond(thr.getAtom("CA"), thr.getAtom("CB"), -160);
    /*
    int numResidues = chain.getNumResidues();
    System.out.printf("Res ID   | Phi      Psi      Omega    |   CA1-C-O  CA1-C-N  O-C-N    C-N-CA2  C-N-H    H-N-CA2\n");
    System.out.printf("---------|----------------------------|-------------------------------------------------------\n");
    for(int i = 1; i <= numResidues; i++){
      double phi = Math.abs(chain.getPhiAngle(i));
      double psi = Math.abs(chain.getPsiAngle(i));
      double omega = Math.abs(chain.getOmegaAngle(i));
      System.out.printf("Res %-5d| %7.2f  %7.2f  %7.2f  |%s\n", i, phi, psi, omega, getAnglesString(getPeptideBondAngles(chain, i)));
    }
    */
    //prot.writeToFile("out.pdb"); // take picture
  }

  @Test
  public void TestDihedralAngles(){
    int numResidues = chain.getNumResidues();
    //System.out.printf("Res ID   | Phi      Psi      Omega    |   CA1-C-N  CA1-C-O  O-C-N    C-N-CA2  C-N-H    H-N-CA2\n");
    //System.out.printf("---------|----------------------------|-------------------------------------------------------\n");
    for(int i = 2; i < numResidues; i++){
      double phi = Math.abs(chain.getPhiAngle(i));
      double psi = Math.abs(chain.getPsiAngle(i));
      double omega = Math.abs(chain.getOmegaAngle(i));
      //System.out.printf("Res %-5d| %7.2f  %7.2f  %7.2f  |%s\n", i, phi, psi, omega, getAnglesString(getPeptideBondAngles(chain, i)));
      assertTrue(Math.abs(omega - 180) < 0.000000001);
      assertTrue(Math.abs(phi - 180) < 0.000000001);
      assertTrue(Math.abs(psi - 180) < 0.000000001);
    }

  }

  @Test
  public void TestPeptideBondAngles(){
    for(int i = 1; i < chain.getNumResidues(); i++){
      double[] angles = getPeptideBondAngles(chain, i);
      assertTrue(Math.abs(angles[0] - 121.00) < 0.01);
      assertTrue(Math.abs(angles[1] - 114.00) < 0.01);
      assertTrue(Math.abs(angles[2] - 125.00) < 0.01);
      assertTrue(Math.abs(angles[3] - 123.00) < 0.01);
      assertTrue(Math.abs(angles[4] - 118.67) < 0.01);
      assertTrue(Math.abs(angles[5] - 118.33) < 0.01);
    }
  }

  private static double[] getPeptideBondAngles(PolypeptideChain chain, int resID){
    Atom ca1 = chain.getResidue(resID).getAtom("CA");
    Atom c = chain.getResidue(resID).getAtom("C");
    Atom o = chain.getResidue(resID).getAtom("O");
    double[] angles = new double[6];
    angles[0] = getAngle(ca1, c, o);
    if(chain.contains(resID+1)){
      Atom n = chain.getResidue(resID+1).getAtom("N");
      Atom h = chain.getResidue(resID+1).getAtom("H");
      Atom ca2 = chain.getResidue(resID+1).getAtom("CA");
      angles[1] = getAngle(ca1, c, n);
      angles[2] = getAngle(o, c, n);
      angles[3] = getAngle(c, n, ca2);
      angles[4] = getAngle(c, n, h);
      angles[5] = getAngle(h, n, ca2);
    } else {
      angles[1] = 1000;
      angles[2] = 1000;
      angles[3] = 1000;
      angles[4] = 1000;
      angles[5] = 1000;
    }
    return angles;
  }

  private String getAnglesString(double[] angles){
    String str = "";
    for(double angle : angles){
      str += String.format("  %7.2f", angle);
    }
    return str;
  }

  private static double getAngle(Atom atom1, Atom atom2, Atom atom3){
    Vector3D one = atom1.getCoordinates();
    Vector3D two = atom2.getCoordinates();
    Vector3D three = atom3.getCoordinates();

    return one.subtract(two).angle(three.subtract(two));
  }
}
