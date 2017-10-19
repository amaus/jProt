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

  @Before
  public void setup(){
    InputStream stream = TestProtein.class.getResourceAsStream("1rop.pdb");
    //PDBFileIO pdb = new PDBFileIO(stream);
    rop = new PDBFileIO().readInPDBFile(stream, "1rop");

    stream = TestProtein.class.getResourceAsStream("5m2j.pdb");
    //pdb = new PDBFileIO(stream);
    m2j = new PDBFileIO().readInPDBFile(stream, "5m2j");
  }

  @Test
  public void testSequence(){
    PolypeptideChain chain = rop.getChain("A");
    //System.out.println(chain.getResidue(1));
    String sequence = "MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARF";
    //System.out.println(rop.getSequence());
    assertTrue((rop.getSequence()).equals(sequence));
    assertTrue(rop.getNumAtoms() == 447);

    chain = m2j.getChain("A");
    sequence = "SDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL";
    assertTrue(chain.getSequence().equals(sequence));

    chain = m2j.getChain("D");
    sequence = "QVQLVESGGGLVQPGGSLRLSCAASGFTFSNYWMYWVRQAPGKGLEWVSEINTNGLITKYPDSVKGRFTISRDNAKNTLYLQMNSLKPEDTALYYCARSPSGFNRGQGTQVTVSS";
    assertTrue(chain.getSequence().equals(sequence));
  }

  @Test
  public void testBondSeparation(){
    PolypeptideChain chainA = m2j.getChain("A");
    Residue ser = chainA.getResidue(9);
    Residue cys69 = chainA.getResidue(69);
    Residue cys101 = chainA.getResidue(101);

    Atom atomOne = ser.getAtom("N");
    Atom atomTwo = ser.getAtom("C");
    assertEquals(m2j.getBondSeparation(atomOne, atomTwo), 2);

    atomOne = cys69.getAtom("SG");
    atomTwo = cys101.getAtom("SG");
    assertEquals(m2j.getBondSeparation(atomOne, atomTwo), 1);

    atomOne = cys69.getAtom("CB");
    atomTwo = cys101.getAtom("CB");
    assertEquals(m2j.getBondSeparation(atomOne, atomTwo), 3);

    Residue one = chainA.getResidue(10);
    Residue two = chainA.getResidue(11);
    atomOne = one.getAtom("O");
    atomTwo = two.getAtom("N");
    assertEquals(m2j.getBondSeparation(atomOne, atomTwo), 2);

    one = chainA.getResidue(10);
    two = chainA.getResidue(12);
    atomOne = one.getAtom("C");
    atomTwo = two.getAtom("C");
    assertEquals(m2j.getBondSeparation(atomOne, atomTwo), 6);
  }

  @Test
  public void testEnableDisableHydrogens(){
    // TODO Test that whether a pdb includes hydrogens or not,
    // when a protein is constructed, the hydrogens are constructed so that they can
    // be enabled. Ensure that after reading in 1rop (which does not include hydrogens),
    // when hydrogens are enabled, it has the correct number of atoms and bonds.
    PolypeptideChain chain = new PolypeptideChain("A");
                                          //{HeavyAtoms, Hydrogens} {main-bonds, h-bonds}
    chain.addResidue(new Residue('I',1)); //  {8,11}  {7,11}

    chain.addResidue(new Residue('A',6)); //  {5,5}   {4,5}
    chain.addResidue(new Residue('M',3)); //  {8,9}   {7,9}

    chain.addResidue(new Residue('S',8)); //  {6,5}   {5,5}
    chain.addResidue(new Residue('T',9)); //  {7,7}   {6,7}
    chain.addResidue(new Residue('A',2)); //  {5,5}   {4,5}
    chain.addResidue(new Residue('R',7)); //  {11,13} {10,13}

    chain.addResidue(new Residue('S',4)); //  {6,5}   {5,5}
    chain.addResidue(new Residue('T',5)); //  {7,7}   {6,7}
    Residue phe = new Residue('F',11); //     {12,9}  {12,9}
    phe.setAsCarboxylTerminus();
    chain.addResidue(phe);
    chain.addResidue(new Residue('F',10)); // {11,9}  {11,9}
                              // TOTALS       {86,85} {77,85}
                              // TOTAL of      171     172 (includes 10 peptide bonds)

    Protein starStuff = new Protein("starStuff");
    starStuff.addChain(chain);

    // Hydrogens start out disabled.
    assertEquals(starStuff.getNumAtoms(), 86);
    assertEquals(starStuff.getNumBonds(), 87);
    // Protein keeps an atoms graph where the edges are bonds
    // ensure that the Proteins number of Bonds matches the
    // number of bonds in the chain
    assertEquals(starStuff.getNumBonds(), chain.getNumBonds());
    starStuff.enableHydrogens();
    assertEquals(starStuff.getNumAtoms(), 171);
    assertEquals(starStuff.getNumBonds(), 172);
    assertEquals(starStuff.getNumBonds(), chain.getNumBonds());

    starStuff.disableHydrogens();
    assertEquals(starStuff.getNumAtoms(), 86);
    assertEquals(starStuff.getNumBonds(), 87);
    assertEquals(starStuff.getNumBonds(), chain.getNumBonds());
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
    assertEquals(prot.getNumBonds(), 11);
    prot.enableHydrogens();
    assertEquals(prot.getNumBonds(), 19);

    // Second Case: between residues in the same chain
    prot = new Protein("CC");
    chainA = new PolypeptideChain("A");
    chainA.addResidue(new Residue('C',1));
    chainA.addResidue(new Residue('C',2));
    prot.addChain(chainA);
    prot.addDisulfideBond("A", 1, "A", 2);
    assertEquals(prot.getNumBonds(), 12);
    prot.enableHydrogens();
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
}
