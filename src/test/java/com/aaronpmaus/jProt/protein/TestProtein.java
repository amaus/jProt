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

  @Before
  public void setup(){
    PDBFileIO pdb = new PDBFileIO();
    InputStream stream = TestProtein.class.getResourceAsStream("1rop.pdb");
    rop = pdb.readInPDBFile(stream, "1rop.pdb");
  }

  @Test
  public void testSequence(){
    PolypeptideChain chain = rop.getChain("A");
    System.out.println(chain.getResidue(1));
    String sequence = "MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARF";
    System.out.println(rop.getSequence());
    assertTrue((rop.getSequence()).equals(sequence));
    assertTrue(rop.getNumAtoms() == 447);
  }
}
