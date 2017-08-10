package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jProt.tools.*;

import static org.junit.Assert.*;
import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Rule;
import org.junit.rules.ExpectedException;

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
  String[] alignment;
  String seq1;
  String seq2;

  @Rule
  public final ExpectedException exception = ExpectedException.none();

  @Test
  public void testDifferentCharacterMatchAlignment(){
    seq1 = "SHAKE";
    seq2 = "SPEARE";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    // The H should align with the E because histidine is more similar to Glutamic Acid than
    // Proline. Lysine and Arginine are also similar enough for the K and R to align.
    assertEquals(alignment[0], "S-HAKE");
    assertEquals(alignment[1], "SPEARE");
  }

  @Test
  public void testSameSequenceAlignment(){
    seq1 = "IAMSTARSTT";
    seq2 = "IAMSTARSTT";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    assertEquals(alignment[0], "IAMSTARSTT");
    assertEquals(alignment[1], "IAMSTARSTT");
  }

  @Test
  public void testInvalidProteinSequenceExpectException(){
    seq1 = "IAMSTARSTUFF";
    seq2 = "IAMSTARSTFF";
    exception.expect(IllegalArgumentException.class);
    SequenceAligner.alignProteinSequences(seq1,seq2);
  }

  @Test
  public void testInvalidDNASequenceExpectException(){
    seq1 = "GAUUACA";
    seq2 = "GAUUACA";
    exception.expect(IllegalArgumentException.class);
    SequenceAligner.alignDNASequences(seq1,seq2);
  }

  @Test
  public void testInvalidRNASequenceExpectException(){
    seq1 = "GATTACA";
    seq2 = "GATTACA";
    exception.expect(IllegalArgumentException.class);
    SequenceAligner.alignRNASequences(seq1,seq2);
  }

  @Test
  public void testEmptySequenceExpectException(){
    seq1 = "";
    seq2 = "MT";
    exception.expect(IllegalArgumentException.class);
    SequenceAligner.alignProteinSequences(seq1,seq2);
  }

  @Test
  public void testSingleGapInMiddle(){
    seq1 = "MTKQ";
    seq2 = "MTAKQ";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    assertEquals(alignment[0], "MT-KQ");
    assertEquals(alignment[1], "MTAKQ");

    seq1 = "MTAKQ";
    seq2 = "MTKQ";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    assertEquals(alignment[0], "MTAKQ");
    assertEquals(alignment[1], "MT-KQ");
  }

  @Test
  public void testGapAtStartAndEnd(){
    seq1 = "M";
    seq2 = "MT";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    assertEquals(alignment[0], "M-");
    assertEquals(alignment[1], "MT");

    seq1 = "MT";
    seq2 = "M";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    assertEquals(alignment[0], "MT");
    assertEquals(alignment[1], "M-");

    seq1 = "T";
    seq2 = "MT";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    assertEquals(alignment[0], "-T");
    assertEquals(alignment[1], "MT");

    seq1 = "MT";
    seq2 = "T";
    alignment = SequenceAligner.alignProteinSequences(seq1,seq2);
    assertEquals(alignment[0], "MT");
    assertEquals(alignment[1], "-T");
  }

  @Test
  public void testAlignmentSingleGapExpected(){
    seq1 = "AAATCA";
    seq2 = "AAAGAATTCA";
    alignment = SequenceAligner.alignDNASequences(seq1,seq2);
    assertEquals(alignment[0], "AAA----TCA");
    assertEquals(alignment[1], "AAAGAATTCA");
  }
}
