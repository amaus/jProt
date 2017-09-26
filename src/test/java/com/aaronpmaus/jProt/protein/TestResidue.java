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

public class TestResidue{

  /*
  @Before
  public void setup(){
    Residue res = new Residue("PRO",3);
    Residue res = new Residue("PHE",4);
    Residue res = new Residue("TYR",5);
    Residue res = new Residue("GLU",6);
    Residue res = new Residue("HIS",7);
    Residue res = new Residue("SER",8);
    Residue res = new Residue("CYS",9);
    Residue res = new Residue("ASN",10);
    Residue res = new Residue("GLY",11);
    Residue res = new Residue("LEU",12);
    Residue res = new Residue("VAL",13);
    Residue res = new Residue("TRP",14);
    Residue res = new Residue("ASP",15);
    Residue res = new Residue("ARG",16);
    Residue res = new Residue("LYS",17);
    Residue res = new Residue("THR",18);
    Residue res = new Residue("MET",19);
    Residue res = new Residue("GLN",20);
  }
  */
  @Test
  public void testAlanine(){
    Residue res = new Residue("ALA",1);
    assertTrue(!res.isMissingAtoms());
    assertTrue(res.getOneLetterName().equals("A"));
    assertTrue(res.getThreeLetterName().equals("ALA"));
    assertTrue(res.getName().equals("Alanine"));
    assertEquals(res.getNumAtoms(), 5);
    assertFalse(res.contains("OXT"));
    res.setAsCarboxylTerminus();
    assertEquals(res.getNumAtoms(), 6);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("OXT"));

    assertFalse(res.contains("H"));
    assertFalse(res.contains("HA"));
    assertFalse(res.contains("1HB"));
    assertFalse(res.contains("2HB"));
    assertFalse(res.contains("3HB"));

    Residue.enableHydrogens();
    assertEquals(res.getNumAtoms(), 11);
    assertTrue(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("3HB"));
    Residue.disableHydrogens();

    for(Bond bond : res.getBonds()){
      assertFalse(bond.containsHydrogen());
    }
  }

  @Test
  public void testArginine(){
    Residue res = new Residue("ARG",1);
    assertTrue(!res.isMissingAtoms());
    assertTrue(res.getOneLetterName().equals("R"));
    assertTrue(res.getThreeLetterName().equals("ARG"));
    assertTrue(res.getName().equals("Arginine"));
    assertEquals(res.getNumAtoms(), 11);
    assertFalse(res.contains("OXT"));
    res.setAsCarboxylTerminus();
    assertEquals(res.getNumAtoms(), 12);

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
    assertTrue(res.contains("OXT"));

    assertFalse(res.contains("H"));
    assertFalse(res.contains("HA"));
    assertFalse(res.contains("1HB"));
    assertFalse(res.contains("2HB"));
    assertFalse(res.contains("1HG"));

    Residue.enableHydrogens();
    assertEquals(res.getNumAtoms(), 25);
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
    Residue.disableHydrogens();

    for(Bond bond : res.getBonds()){
      assertFalse(bond.containsHydrogen());
    }
  }

  @Test
  public void testIsoleucine(){
    Residue res = new Residue("ILE",1);
    assertTrue(!res.isMissingAtoms());
    assertTrue(res.getOneLetterName().equals("I"));
    assertTrue(res.getThreeLetterName().equals("ILE"));
    assertTrue(res.getName().equals("Isoleucine"));
    assertEquals(res.getNumAtoms(), 8);
    assertFalse(res.contains("OXT"));
    res.setAsCarboxylTerminus();
    assertEquals(res.getNumAtoms(), 9);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("CG1"));
    assertTrue(res.contains("CG2"));
    assertTrue(res.contains("CD1"));
    assertTrue(res.contains("OXT"));

    assertFalse(res.contains("H"));
    assertFalse(res.contains("HA"));
    assertFalse(res.contains("1HB"));
    assertFalse(res.contains("2HB"));
    assertFalse(res.contains("3HB"));

    Residue.enableHydrogens();
    assertEquals(res.getNumAtoms(), 20);
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
    Residue.disableHydrogens();

    for(Bond bond : res.getBonds()){
      assertFalse(bond.containsHydrogen());
    }
  }

  @Test
  public void testProline(){
    Residue res = new Residue("PRO",1);
    assertTrue(!res.isMissingAtoms());
    assertTrue(res.getOneLetterName().equals("P"));
    assertTrue(res.getThreeLetterName().equals("PRO"));
    assertTrue(res.getName().equals("Proline"));
    assertEquals(res.getNumAtoms(), 7);
    assertFalse(res.contains("OXT"));
    res.setAsCarboxylTerminus();
    assertEquals(res.getNumAtoms(), 8);

    assertTrue(res.contains("C"));
    assertTrue(res.contains("CA"));
    assertTrue(res.contains("O"));
    assertTrue(res.contains("N"));
    assertTrue(res.contains("CB"));
    assertTrue(res.contains("OXT"));

    assertFalse(res.contains("H"));
    assertFalse(res.contains("HA"));
    assertFalse(res.contains("1HB"));
    assertFalse(res.contains("2HB"));
    assertFalse(res.contains("3HB"));

    Residue.enableHydrogens();
    assertEquals(res.getNumAtoms(), 15);
    assertFalse(res.contains("H"));
    assertTrue(res.contains("HA"));
    assertTrue(res.contains("1HB"));
    assertTrue(res.contains("2HB"));
    assertTrue(res.contains("1HG"));
    assertTrue(res.contains("2HG"));
    assertTrue(res.contains("1HD"));
    assertTrue(res.contains("2HD"));
    Residue.disableHydrogens();

    for(Bond bond : res.getBonds()){
      assertFalse(bond.containsHydrogen());
    }
  }

}
