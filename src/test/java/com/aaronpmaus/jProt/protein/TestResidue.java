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

import java.util.ArrayList;
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

public class TestResidue{

  /**
  * Test that they residue has the proper number of atoms and bonds before and after both
  * setting it as the Carboxyl Terminus and enabling hydrogens
  *
  * @param res the Residue to test
  * @param numAllAtoms the number of atoms including the carboxyl terminus oxygen and all hydrogens
  * @param numAllBonds the number of bonds including the carboxyl terminus oxygen and all hydrogens
  */
  private void testResidueNumAtomsAndBonds(Residue res, int numAllAtoms, int numAllBonds){

    // res starts out not as Carboxyl Terminus
    assertFalse(res.contains("OXT"));
    assertEquals(res.getNumAtoms(), numAllAtoms-1);
    assertEquals(res.getBonds().size(),numAllBonds-1);

    // set as Carboxyl Terminus and test num atoms, bonds, and that it contains OXT
    res.setAsCarboxylTerminus();
    assertEquals(res.getNumAtoms(), numAllAtoms);
    assertEquals(res.getBonds().size(),numAllBonds);
    assertTrue(res.contains("OXT"));


  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // A-ALA {5,5} {4,5}
  @Test
  public void testAlanine(){
    Residue res = new Residue("ALA",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("A"));
    assertTrue(res.getThreeLetterName().equals("ALA"));
    assertTrue(res.getName().equals("Alanine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 11, 10);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("3HB"));

  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // R-ARG {11,13} {10,13}
  @Test
  public void testArginine(){
    Residue res = new Residue("ARG",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("R"));
    assertTrue(res.getThreeLetterName().equals("ARG"));
    assertTrue(res.getName().equals("Arginine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 25, 24);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD"));
    assertTrue(res.contains("NE"));
    assertTrue(res.contains("CZ"));
    assertTrue(res.contains("NH1"));
    assertTrue(res.contains("NH2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("2HG"));
    assertTrue(res.contains("1HD"));
    assertTrue(res.contains("2HD"));
    assertTrue(res.contains("HE"));
    assertTrue(res.contains("1HH1"));
    assertTrue(res.contains("2HH1"));
    assertTrue(res.contains("1HH2"));
    assertTrue(res.contains("2HH2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // N-ASN {8,6} {7,6}
  @Test
  public void testAsparagine(){
    Residue res = new Residue("ASN",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("N"));
    assertTrue(res.getThreeLetterName().equals("ASN"));
    assertTrue(res.getName().equals("Asparagine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 15, 14);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("OD1"));
    assertTrue(res.contains("ND2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HD2"));
    assertTrue(res.contains("2HD2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // D-ASP {8,4} {7,4}
  @Test
  public void testAsparticAcid(){
    Residue res = new Residue("ASP",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("D"));
    assertTrue(res.getThreeLetterName().equals("ASP"));
    assertTrue(res.getName().equals("Aspartic Acid"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 13, 12);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("OD1"));
    assertTrue(res.contains("OD2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // D-ASH {8,5} {7,5}
  @Test
  public void testProtonatedAsparticAcid(){
    Residue res = new Residue("ASH",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("D"));
    assertTrue(res.getThreeLetterName().equals("ASH"));
    assertTrue(res.getName().equals("Aspartic Acid"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 14, 13);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("OD1"));
    assertTrue(res.contains("OD2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HD2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // C-CYS {6,4} {5,4}
  @Test
  public void testCysteine(){
    Residue res = new Residue("CYS",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("C"));
    assertTrue(res.getThreeLetterName().equals("CYS"));
    assertTrue(res.getName().equals("Cysteine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 11, 10);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("SG"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // C-CYH {6,5} {5,5}
  @Test
  public void testProtonatedCysteine(){
    Residue res = new Residue("CYH",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("C"));
    assertTrue(res.getThreeLetterName().equals("CYH"));
    assertTrue(res.getName().equals("Cysteine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 12, 11);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HG"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // Q-GLN {9,8} {8,8}
  @Test
  public void testGlutamine(){
    Residue res = new Residue("GLN",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("Q"));
    assertTrue(res.getThreeLetterName().equals("GLN"));
    assertTrue(res.getName().equals("Glutamine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 18, 17);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD"));
    assertTrue(res.contains("OE1"));
    assertTrue(res.contains("NE2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("2HG"));
    assertTrue(res.contains("1HE2"));
    assertTrue(res.contains("2HE2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // E-GLU {9,6} {8,6}
  @Test
  public void testGlutamicAcid(){
    Residue res = new Residue("GLU",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("E"));
    assertTrue(res.getThreeLetterName().equals("GLU"));
    assertTrue(res.getName().equals("Glutamic Acid"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 16, 15);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD"));
    assertTrue(res.contains("OE1"));
    assertTrue(res.contains("OE2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // E-GLH {9,7} {8,7}
  @Test
  public void testProtonatedGlutamicAcid(){
    Residue res = new Residue("GLH",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("E"));
    assertTrue(res.getThreeLetterName().equals("GLH"));
    assertTrue(res.getName().equals("Glutamic Acid"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 17, 16);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD"));
    assertTrue(res.contains("OE1"));
    assertTrue(res.contains("OE2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("HE2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // G-GLY {4,3} {3,3}
  @Test
  public void testGlycine(){
    Residue res = new Residue("GLY",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("G"));
    assertTrue(res.getThreeLetterName().equals("GLY"));
    assertTrue(res.getName().equals("Glycine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 8, 7);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("1HA"));
    assertTrue(res.contains("2HA"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // H-HIS {10,7} {10,7}
  @Test
  public void testHistidine(){
    Residue res = new Residue("HIS",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("H"));
    assertTrue(res.getThreeLetterName().equals("HIS"));
    assertTrue(res.getName().equals("Histidine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 18, 18);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("ND1"));
    assertTrue(res.contains("CD2"));
    assertTrue(res.contains("CE1"));
    assertTrue(res.contains("NE2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HE1"));
    assertTrue(res.contains("HE2"));
    assertTrue(res.contains("HD2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // H-HID {10,7} {10,7}
  @Test
  public void testProtonatedHistidinD(){
    Residue res = new Residue("HID",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("H"));
    assertTrue(res.getThreeLetterName().equals("HID"));
    assertTrue(res.getName().equals("Histidine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 18, 18);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("ND1"));
    assertTrue(res.contains("CD2"));
    assertTrue(res.contains("CE1"));
    assertTrue(res.contains("NE2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HD1"));
    assertTrue(res.contains("HD2"));
    assertTrue(res.contains("HE1"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // H-HIE {10,7} {10,7}
  @Test
  public void testProtonatedHistidinE(){
    Residue res = new Residue("HIE",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("H"));
    assertTrue(res.getThreeLetterName().equals("HIE"));
    assertTrue(res.getName().equals("Histidine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 18, 18);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("ND1"));
    assertTrue(res.contains("CD2"));
    assertTrue(res.contains("CE1"));
    assertTrue(res.contains("NE2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HD2"));
    assertTrue(res.contains("HE1"));
    assertTrue(res.contains("HE2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // H-HIP {10,8} {10,8}
  @Test
  public void testProtonatedHistidinP(){
    Residue res = new Residue("HIP",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("H"));
    assertTrue(res.getThreeLetterName().equals("HIP"));
    assertTrue(res.getName().equals("Histidine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 19, 19);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("ND1"));
    assertTrue(res.contains("CD2"));
    assertTrue(res.contains("CE1"));
    assertTrue(res.contains("NE2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HD1"));
    assertTrue(res.contains("HD2"));
    assertTrue(res.contains("HE1"));
    assertTrue(res.contains("HE2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // I-ILE {8,11} {7,11}
  @Test
  public void testIsoleucine(){
    Residue res = new Residue("ILE",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("I"));
    assertTrue(res.getThreeLetterName().equals("ILE"));
    assertTrue(res.getName().equals("Isoleucine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 20, 19);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG1"));
    assertTrue(res.contains("CG2"));
    assertTrue(res.contains("CD1"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("HB"));
    assertTrue(res.contains("1HG1"));
    assertTrue(res.contains("2HG1"));
    assertTrue(res.contains("1HG2"));
    assertTrue(res.contains("2HG2"));
    assertTrue(res.contains("3HG2"));
    assertTrue(res.contains("1HD1"));
    assertTrue(res.contains("2HD1"));
    assertTrue(res.contains("3HD1"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // L-LEU {8,11} {7,11}
  @Test
  public void testLeucine(){
    Residue res = new Residue("LEU",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("L"));
    assertTrue(res.getThreeLetterName().equals("LEU"));
    assertTrue(res.getName().equals("Leucine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 20, 19);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD1"));
    assertTrue(res.contains("CD2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HG"));
    assertTrue(res.contains("1HD1"));
    assertTrue(res.contains("2HD1"));
    assertTrue(res.contains("3HD1"));
    assertTrue(res.contains("1HD2"));
    assertTrue(res.contains("2HD2"));
    assertTrue(res.contains("3HD2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // K-LYN {9,12} {8,12}
  @Test
  public void testProtonatedLysine(){
    Residue res = new Residue("LYN",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("K"));
    assertTrue(res.getThreeLetterName().equals("LYN"));
    assertTrue(res.getName().equals("Lysine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 22, 21);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD"));
    assertTrue(res.contains("CE"));
    assertTrue(res.contains("NZ"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("2HG"));
    assertTrue(res.contains("1HD"));
    assertTrue(res.contains("2HD"));
    assertTrue(res.contains("1HE"));
    assertTrue(res.contains("2HE"));
    assertTrue(res.contains("1HZ"));
    assertTrue(res.contains("2HZ"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // K-LYS {9,13} {8,13}
  @Test
  public void testLysine(){
    Residue res = new Residue("LYS",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("K"));
    assertTrue(res.getThreeLetterName().equals("LYS"));
    assertTrue(res.getName().equals("Lysine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 23, 22);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD"));
    assertTrue(res.contains("CE"));
    assertTrue(res.contains("NZ"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("2HG"));
    assertTrue(res.contains("1HD"));
    assertTrue(res.contains("2HD"));
    assertTrue(res.contains("1HE"));
    assertTrue(res.contains("2HE"));
    assertTrue(res.contains("1HZ"));
    assertTrue(res.contains("2HZ"));
    assertTrue(res.contains("3HZ"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // M-MET {8,9} {7,9}
  @Test
  public void testMethionine(){
    Residue res = new Residue("MET",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("M"));
    assertTrue(res.getThreeLetterName().equals("MET"));
    assertTrue(res.getName().equals("Methionine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 18, 17);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("SD"));
    assertTrue(res.contains("CE"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("2HG"));
    assertTrue(res.contains("1HE"));
    assertTrue(res.contains("2HE"));
    assertTrue(res.contains("3HE"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // F-PHE {11,9} {11,9}
  @Test
  public void testPhenylalanine(){
    Residue res = new Residue("PHE",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("F"));
    assertTrue(res.getThreeLetterName().equals("PHE"));
    assertTrue(res.getName().equals("Phenylalanine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 21, 21);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD1"));
    assertTrue(res.contains("CD2"));
    assertTrue(res.contains("CE1"));
    assertTrue(res.contains("CE2"));
    assertTrue(res.contains("CZ"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HD1"));
    assertTrue(res.contains("HD2"));
    assertTrue(res.contains("HE1"));
    assertTrue(res.contains("HE2"));
    assertTrue(res.contains("HZ"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // P-PRO {7,7} {7,7}
  @Test
  public void testProline(){
    Residue res = new Residue("PRO",1);
    assertTrue(!res.isMissingAtoms());
    assertTrue(res.getOneLetterName().equals("P"));
    assertTrue(res.getThreeLetterName().equals("PRO"));
    assertTrue(res.getName().equals("Proline"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 15, 15);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD"));

    assertFalse(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("2HG"));
    assertTrue(res.contains("1HD"));
    assertTrue(res.contains("2HD"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // S-SER {6,5} {5,5}
  @Test
  public void testSerine(){
    Residue res = new Residue("SER",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("S"));
    assertTrue(res.getThreeLetterName().equals("SER"));
    assertTrue(res.getName().equals("Serine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 12, 11);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("OG"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HG"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // T-THR {7,7} {6,7}
  @Test
  public void testThreonine(){
    Residue res = new Residue("THR",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("T"));
    assertTrue(res.getThreeLetterName().equals("THR"));
    assertTrue(res.getName().equals("Threonine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 15, 14);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("OG1"));
    assertTrue(res.contains("CG2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("HB"));
    assertTrue(res.contains("HG1"));
    assertTrue(res.contains("1HG2"));
    assertTrue(res.contains("2HG2"));
    assertTrue(res.contains("3HG2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // W-TRP {14,10} {15,10}
  @Test
  public void testTryptophan(){
    Residue res = new Residue("TRP",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("W"));
    assertTrue(res.getThreeLetterName().equals("TRP"));
    assertTrue(res.getName().equals("Tryptophan"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 25, 26);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD1"));
    assertTrue(res.contains("CD2"));
    assertTrue(res.contains("NE1"));
    assertTrue(res.contains("CE2"));
    assertTrue(res.contains("CE3"));
    assertTrue(res.contains("CZ2"));
    assertTrue(res.contains("CZ3"));
    assertTrue(res.contains("CH2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HD1"));
    assertTrue(res.contains("HE1"));
    assertTrue(res.contains("HE3"));
    assertTrue(res.contains("HZ2"));
    assertTrue(res.contains("HZ3"));
    assertTrue(res.contains("HH2"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // Y-TYR {12,9} {12,9}
  @Test
  public void testTyrosine(){
    Residue res = new Residue("TYR",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("Y"));
    assertTrue(res.getThreeLetterName().equals("TYR"));
    assertTrue(res.getName().equals("Tyrosine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 22, 22);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG"));
    assertTrue(res.contains("CD1"));
    assertTrue(res.contains("CD2"));
    assertTrue(res.contains("CE1"));
    assertTrue(res.contains("CE2"));
    assertTrue(res.contains("CZ"));
    assertTrue(res.contains("OH"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("HD1"));
    assertTrue(res.contains("HD2"));
    assertTrue(res.contains("HE1"));
    assertTrue(res.contains("HE2"));
    assertTrue(res.contains("HH"));
  }

  // RES {HeavyAtoms, Hydrogens} {main-bonds, h-bonds)}
  // Doesn't include OXT or bond to OXT
  // V-VAL {7,9} {6,9}
  @Test
  public void testValine(){
    Residue res = new Residue("VAL",1);
    assertTrue(!res.isMissingAtoms());

    assertTrue(res.getOneLetterName().equals("V"));
    assertTrue(res.getThreeLetterName().equals("VAL"));
    assertTrue(res.getName().equals("Valine"));

    // (Residue, numAllAtoms, numAllBonds)
    testResidueNumAtomsAndBonds(res, 17, 16);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG1"));
    assertTrue(res.contains("CG2"));

    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("HB"));
    assertTrue(res.contains("1HG1"));
    assertTrue(res.contains("2HG1"));
    assertTrue(res.contains("3HG1"));
    assertTrue(res.contains("1HG2"));
    assertTrue(res.contains("2HG2"));
    assertTrue(res.contains("3HG2"));
  }
  
}
