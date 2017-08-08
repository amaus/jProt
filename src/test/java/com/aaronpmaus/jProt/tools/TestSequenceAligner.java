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

  @Test
  public void testAlign(){
    String sequence = "MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARF";
    String seq1 = "MTKQ";
    String seq2 = "MTAKQ";
    SequenceAligner aligner = new SequenceAligner();
    String[] alignment = aligner.alignProteinSequences(seq1,seq2);
    System.out.println(alignment[0]);
    System.out.println(alignment[1]);
    alignment = aligner.alignProteinSequences(seq1,seq1);
    System.out.println(alignment[0]);
    System.out.println(alignment[1]);

    seq1 = "SHAKE";
    seq2 = "SPEARE";
    alignment = aligner.alignProteinSequences(seq1,seq2);
    System.out.println(alignment[0]);
    System.out.println(alignment[1]);

    seq1 = "ATGGCGT";
    seq2 = "ATGAGT";
    alignment = aligner.alignNucleotideSequences(seq1,seq2);
    System.out.println(alignment[0]);
    System.out.println(alignment[1]);

    seq1 = "AAAGAATTCA";
    seq2 = "AAATCA";
    alignment = aligner.alignNucleotideSequences(seq1,seq2);
    System.out.println(alignment[0]);
    System.out.println(alignment[1]);


  }
}
