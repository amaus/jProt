package com.aaronpmaus.jProt;

import com.aaronpmaus.jProt.protein.*;
import com.aaronpmaus.jProt.io.*;
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
    prot = VirtualRibosome.synthesizeProtein(seq);
    prot.enableHydrogens();
    chain = prot.getChain("A");

    Residue phe = chain.getResidue(10);
    phe.rotateAboutBond("CA", "CB", 30);
    phe.rotateAboutBond("CB", "CG", -90);
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
    for(int i = 2; i < chain.getNumResidues(); i++){
      double[] angles = getPeptideBondAngles(chain, i);
      assertTrue(Math.abs(angles[0] - 114.50) < 0.01);
      assertTrue(Math.abs(angles[1] - 120.50) < 0.01);
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
    Atom n = chain.getResidue(resID+1).getAtom("N");
    Atom h = chain.getResidue(resID+1).getAtom("H");
    Atom ca2 = chain.getResidue(resID+1).getAtom("CA");
    double[] angles = new double[6];
    angles[0] = getAngle(ca1, c, n);
    angles[1] = getAngle(ca1, c, o);
    angles[2] = getAngle(o, c, n);
    angles[3] = getAngle(c, n, ca2);
    angles[4] = getAngle(c, n, h);
    angles[5] = getAngle(h, n, ca2);
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
