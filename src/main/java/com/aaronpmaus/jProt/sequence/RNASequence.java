package com.aaronpmaus.jProt.sequence;

import java.util.ArrayList;

/**
* A RNASequence is a String of single letter identifier specifying the sequence of the bases of a
* strand of RNA
* <p>
* Valid characters are only G A U C.
* @since 0.7.0
*/
public class RNASequence extends Sequence {

  /**
  * Construct a RNA Sequence. The sequences must contain only the characters G A U C.
  * @param seq a String containing only valid RNA base characters
  */
  public RNASequence(String seq){
    super(seq);
    validateSequence(seq);
  }

  @Override
  public Alignment align(Sequence other){
    return SequenceAligner.align(this, other, "RNA");
  }

  /**
  * Verify that the seq passed in is a valid RNA sequence.
  *
  * A valid RNA sequence can only contain G A U C
  *
  * @param seq the sequence to check
  * @throws IllegalArgumentException if the sequence is invalid
  */
  private void validateSequence(String seq){
    ArrayList<Character> validCharacters = new ArrayList<Character>();
    validCharacters.add('G');
    validCharacters.add('A');
    validCharacters.add('U');
    validCharacters.add('C');

    for(Character nucleotide : seq.toCharArray()){
      if(!validCharacters.contains(nucleotide)){
        throw new IllegalArgumentException("Seq invalid. " + nucleotide + " not a valid RNA base."
            + "\n only G A U C allowed.");
      }
    }
  }
}
