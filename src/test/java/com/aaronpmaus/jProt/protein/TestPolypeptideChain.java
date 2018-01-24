package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.protein.*;
import com.aaronpmaus.jProt.io.*;
import com.aaronpmaus.jProt.sequence.*;

import static org.junit.Assert.*;
import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Rule;
import org.junit.rules.ExpectedException;

import java.util.Collection;
import java.util.ArrayList;
import java.util.Scanner;
import java.math.BigDecimal;

import java.io.InputStream;
import java.io.File;
import java.io.FileNotFoundException;
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

public class TestPolypeptideChain{
  PolypeptideChain chain;

  @Before
  public void setup(){
    chain = new PolypeptideChain("A");
                                          // {main-bonds, h-bonds}
    chain.addResidue(new Residue('I',1)); //  {7,11}

    chain.addResidue(new Residue('A',6)); //  {4,5}
    chain.addResidue(new Residue('M',3)); //  {7,9}

    chain.addResidue(new Residue('S',8)); //  {5,5}
    chain.addResidue(new Residue('T',9)); //  {6,7}
    chain.addResidue(new Residue('A',2)); //  {4,5}
    chain.addResidue(new Residue('R',7)); //  {10,13}

    chain.addResidue(new Residue('S',4)); //  {5,5}
    chain.addResidue(new Residue('T',5)); //  {6,7}
    Residue phe = new Residue('F',11); // {12,9}
    phe.setAsCarboxylTerminus();
    chain.addResidue(phe);
    chain.addResidue(new Residue('F',10)); // {11,9}
                              // TOTALS {77,85}
                              // TOTAL of 172 (incudes 10 peptide bonds)
  }

  @Test
  public void testGetSequence(){
    assertEquals(chain.getSequence().toString(), "IAMSTARSTFF");
  }

  @Test
  public void testGetResidueIDs(){
    // even though the residues were added out of order
    // this still returns them in the proper order
    Integer[] ids = chain.getResidueIDs();
    int counter = 1;
    for(Integer id : ids){
      assertEquals(id, new Integer(counter));
      counter++;
    }
  }

  @Test
  public void testGetResidue(){
    // ensure that getting a residue by ID returns the residue with the
    // proper ID
    assertEquals(chain.getResidue(1).getResidueID(), 1);
    assertEquals(chain.getResidue(2).getResidueID(), 2);
    assertEquals(chain.getResidue(3).getResidueID(), 3);
    assertEquals(chain.getResidue(4).getResidueID(), 4);
    assertEquals(chain.getResidue(5).getResidueID(), 5);
    assertEquals(chain.getResidue(6).getResidueID(), 6);
    assertEquals(chain.getResidue(7).getResidueID(), 7);
    assertEquals(chain.getResidue(8).getResidueID(), 8);
    assertEquals(chain.getResidue(9).getResidueID(), 9);
    assertEquals(chain.getResidue(10).getResidueID(), 10);
    assertEquals(chain.getResidue(11).getResidueID(), 11);
  }

  @Test
  public void testResidueIterator(){
    // ensure that the iterator iterates over the residues in order
    int id = 1;
    for(Residue res : chain){
      assertEquals(res.getResidueID(),id);
      id++;
    }
  }

  @Test
  public void testHydrogensToggle(){
    chain.enableHydrogens();
    Collection<Bond> bonds = chain.getBonds();
    assertEquals(bonds.size(), 172);

    chain.disableHydrogens();
    bonds = chain.getBonds();
    assertEquals(bonds.size(), 87);
    for(Bond bond : bonds){
      assertFalse(bond.containsHydrogen());
    }
  }

  @Test
  public void testDihedralAngles(){
    Protein prot = VirtualRibosome.synthesizeProtein(new ProteinSequence("IAMSTARSTFF"),
                                                     "strstf");
    PolypeptideChain chain = prot.getChain("A");

    assertTrue(Math.abs(Math.abs(chain.getOmegaAngle(2)) - 180.0) < 0.00000001);
    chain.setOmegaAngle(2, 30.0);
    assertTrue(Math.abs(chain.getOmegaAngle(2) - 30.0) < 0.00000001);

    assertTrue(Math.abs(Math.abs(chain.getPhiAngle(2)) - 180.0) < 0.00000001);
    chain.setPhiAngle(2, 60.0);
    assertTrue(Math.abs(chain.getPhiAngle(2) - 60.0) < 0.00000001);

    assertTrue(Math.abs(Math.abs(chain.getPsiAngle(2)) - 180.0) < 0.00000001);
    chain.setPsiAngle(2, 90.0);
    assertTrue(Math.abs(chain.getPsiAngle(2) - 90.0) < 0.00000001);

    chain.setPsiAngle(2, 180.0);
    assertTrue(Math.abs(chain.getPsiAngle(2) - 180.0) < 0.00000001);

    chain.setPsiAngle(2, -90.0);
    assertTrue(Math.abs(chain.getPsiAngle(2) - -90.0) < 0.00000001);
  }
}
