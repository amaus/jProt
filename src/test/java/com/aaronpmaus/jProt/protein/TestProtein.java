package com.aaronpmaus.jProt;

import com.aaronpmaus.jProt.protein.*;
import com.aaronpmaus.jProt.io.*;

import static org.junit.Assert.*;
import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Rule;
import org.junit.rules.ExpectedException;

import java.io.InputStream;
import java.util.Scanner;
import java.util.ArrayList;

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

public class TestProtein{
  private Protein rop;
  private Protein m2j;
  private ArrayList<Double[]> dihedrals;

  @Before
  public void setup(){
    InputStream stream = TestProtein.class.getResourceAsStream("1rop.pdb");
    rop = new PDBFileIO().readInPDBFile(stream, "1rop");

    stream = TestProtein.class.getResourceAsStream("5m2j.pdb");
    m2j = new PDBFileIO().readInPDBFile(stream, "5m2j");
    dihedrals = readInDihedrals();
  }

  public ArrayList<Double[]> readInDihedrals(){
    ArrayList<Double[]> dihedrals = new ArrayList<Double[]>();
    InputStream stream = TestProtein.class.getResourceAsStream("1rop.dihedrals");
    Scanner in = new Scanner(stream);
    while(in.hasNext()){
      String[] tokens = in.nextLine().split(",");
      Double[] angles = {Double.parseDouble(tokens[0]),
                         Double.parseDouble(tokens[1]),
                         Double.parseDouble(tokens[2])};
      dihedrals.add(angles);
    }
    return dihedrals;
  }

  @Test
  public void testSequence(){
    PolypeptideChain chain = rop.getChain("A");
    //System.out.println(chain.getResidue(1));
    String sequence = "MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARF";
    //System.out.println(rop.getSequence());
    assertTrue(rop.getSequence().toString().equals(sequence));
    assertTrue(rop.getNumAtoms() == 447);

    chain = m2j.getChain("A");
    sequence = "SDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTK"
             + "VNLLSAIKSPCQPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL";
    assertTrue(chain.getSequence().toString().equals(sequence));

    chain = m2j.getChain("D");
    sequence = "QVQLVESGGGLVQPGGSLRLSCAASGFTFSNYWMYWVRQAPGKGLEWVSEINTNGLITKYPDSVKGRFTISRDNAKNTLYLQ"
             + "MNSLKPEDTALYYCARSPSGFNRGQGTQVTVSS";
    assertTrue(chain.getSequence().toString().equals(sequence));
  }

  @Test
  public void testAddingDisulfideBond(){
    // Ensure that disulfide bonds are handled correctly whether they are
    // First Case: between residues in separate chains
    Protein prot = new Protein("CC");
    PolypeptideChain chainA = new PolypeptideChain("A");
    PolypeptideChain chainB = new PolypeptideChain("B");
    chainA.addResidue(new Residue('C',1));
    chainB.addResidue(new Residue('C',1));

    prot.addChain(chainA);
    prot.addChain(chainB);
    prot.addDisulfideBond("A", 1, "B", 1);
    assertEquals(prot.getNumBonds(), 19);

    // Second Case: between residues in the same chain
    prot = new Protein("CC");
    chainA = new PolypeptideChain("A");
    chainA.addResidue(new Residue('C',1));
    chainA.addResidue(new Residue('C',2));
    prot.addChain(chainA);
    prot.addDisulfideBond("A", 1, "A", 2);
    assertEquals(prot.getNumBonds(), 20);
  }

  @Test
  public void testAtomCharges(){
    PolypeptideChain chain = rop.getChain("A");
    Atom atomOne = chain.getResidue(1).getAtom("O"); // charge   -1
    Atom atomTwo = chain.getResidue(2).getAtom("O"); // charge   1-
    Atom atomThree = chain.getResidue(2).getAtom("OG1"); // charge +2
    Atom atomFour = chain.getResidue(3).getAtom("O"); // charge   2+
    Atom atomFive = chain.getResidue(4).getAtom("O"); // charge   3
    Atom atomSix = chain.getResidue(4).getAtom("OE1"); // charge  3
    assertTrue(Math.abs(atomOne.getCharge() - -1.0) < 0.000000001);
    assertTrue(Math.abs(atomTwo.getCharge() - -1.0) < 0.000000001);
    assertTrue(Math.abs(atomThree.getCharge() - 2.0) < 0.000000001);
    assertTrue(Math.abs(atomFour.getCharge() - 2.0) < 0.000000001);
    assertTrue(Math.abs(atomFive.getCharge() - 3.0) < 0.000000001);
    assertTrue(Math.abs(atomSix.getCharge() - 3.0) < 0.000000001);
  }

  @Test
  public void testDihedralAngles(){
    PolypeptideChain chain = rop.getChain("A");
    //System.out.printf("Res ID   | Omega   Phi     Psi\n");
    for(int resID : chain.getResidueIDs()){
      double omega = chain.getOmegaAngle(resID);
      double phi = chain.getPhiAngle(resID);
      double psi = chain.getPsiAngle(resID);
      double expectedOmega = dihedrals.get(resID-1)[0];
      double expectedPhi = dihedrals.get(resID-1)[1];
      double expectedPsi = dihedrals.get(resID-1)[2];
      //System.out.printf("Res %d    | %.2f  %.2f  %.2f\n", resID, omega, phi, psi);
      assertTrue(Math.abs(omega - expectedOmega) < 0.001);
      assertTrue(Math.abs(phi - expectedPhi) < 0.001);
      assertTrue(Math.abs(psi - expectedPsi) < 0.001);
    }
  }
}
