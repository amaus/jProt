package com.aaronpmaus.jProt.protein;

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
}
