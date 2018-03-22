package com.aaronpmaus.jProt.sequence;

import com.aaronpmaus.jProt.sequence.*;

import static org.junit.Assert.*;
import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Rule;
import org.junit.rules.ExpectedException;

import java.util.Arrays;

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

public class TestSequenceAligner{
  Alignment alignment;

  @Rule
  public final ExpectedException exception = ExpectedException.none();

  @Test
  public void testDifferentCharacterMatchAlignment(){
    Sequence seq1 = new ProteinSequence("SHAKE");
    Sequence seq2 = new ProteinSequence("SPEARE");
    alignment = seq1.align(seq2);
    // The H should align with the E because histidine is more similar to Glutamic Acid than
    // Proline. Lysine and Arginine are also similar enough for the K and R to align.
    assertEquals(alignment.getAlignment(seq1), "S-HAKE");
    assertEquals(alignment.getAlignment(seq2), "SPEARE");
  }

  @Test
  public void testSameSequenceAlignment(){
    ProteinSequence seq1 = new ProteinSequence("IAMSTARSTT");
    ProteinSequence seq2 = new ProteinSequence("IAMSTARSTT");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "IAMSTARSTT");
    assertEquals(alignment.getAlignment(seq2), "IAMSTARSTT");
  }

  @Test
  public void testInvalidProteinSequenceExpectException(){
    exception.expect(IllegalArgumentException.class);
    ProteinSequence seq1 = new ProteinSequence("IAMSTARSUTT");
    ProteinSequence seq2 = new ProteinSequence("IAMSTARSTT");
    alignment = seq1.align(seq2);
  }

  @Test
  public void testInvalidDNASequenceExpectException(){
    exception.expect(IllegalArgumentException.class);
    DNASequence seq1 = new DNASequence("GAUUACA");
    DNASequence seq2 = new DNASequence("GAUUACA");
    alignment = seq1.align(seq2);
  }

  @Test
  public void testInvalidRNASequenceExpectException(){
    exception.expect(IllegalArgumentException.class);
    RNASequence seq1 = new RNASequence("GATTACA");
    RNASequence seq2 = new RNASequence("GATTACA");
    alignment = seq1.align(seq2);
  }

  @Test
  public void testEmptySequenceExpectException(){
    ProteinSequence seq1 = new ProteinSequence("");
    ProteinSequence seq2 = new ProteinSequence("MT");
    exception.expect(IllegalArgumentException.class);
    alignment = seq1.align(seq2);
  }

  @Test
  public void testSingleGapInMiddle(){
    ProteinSequence seq1 = new ProteinSequence("MTKQ");
    ProteinSequence seq2 = new ProteinSequence("MTAKQ");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "MT-KQ");
    assertEquals(alignment.getAlignment(seq2), "MTAKQ");

    seq1 = new ProteinSequence("MTAKQ");
    seq2 = new ProteinSequence("MTKQ");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "MTAKQ");
    assertEquals(alignment.getAlignment(seq2), "MT-KQ");
  }

  @Test
  public void testGapAtStartAndEnd(){
    ProteinSequence seq1 = new ProteinSequence("M");
    ProteinSequence seq2 = new ProteinSequence("MT");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "M-");
    assertEquals(alignment.getAlignment(seq2), "MT");

    seq1 = new ProteinSequence("MT");
    seq2 = new ProteinSequence("M");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "MT");
    assertEquals(alignment.getAlignment(seq2), "M-");

    seq1 = new ProteinSequence("T");
    seq2 = new ProteinSequence("MT");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "-T");
    assertEquals(alignment.getAlignment(seq2), "MT");

    seq1 = new ProteinSequence("MT");
    seq2 = new ProteinSequence("T");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "MT");
    assertEquals(alignment.getAlignment(seq2), "-T");
  }

  @Test
  public void testAlignmentSingleGapExpected(){
    DNASequence seq1 = new DNASequence("AAATCA");
    DNASequence seq2 = new DNASequence("AAAGAATTCA");
    alignment = seq1.align(seq2);
    assertEquals(alignment.getAlignment(seq1), "AAA----TCA");
    assertEquals(alignment.getAlignment(seq2), "AAAGAATTCA");
  }

  @Test
  public void testGetSequenceMatchMasks(){
    boolean[][] masks = SequenceAligner.getSequenceMatchMasks("AA-EYE","AAP-E-");
    boolean[] mask0 = {true, true, false, true, false};
    boolean[] mask1 = {true, true, false, true};
    assertTrue(Arrays.equals(masks[0],mask0));
    assertTrue(Arrays.equals(masks[1],mask1));
    String seq1 = "AAA----TCA";
    String seq2 = "AAAGAATTCA";
    masks = SequenceAligner.getSequenceMatchMasks(seq1,seq2);
    boolean[] mask2 = {true, true, true, true, true, true};
    boolean[] mask3 = {true, true, true, false, false, false, false, true, true, true};
    assertTrue(Arrays.equals(masks[0],mask2));
    assertTrue(Arrays.equals(masks[1],mask3));
  }
}
