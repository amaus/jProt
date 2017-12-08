package com.aaronpmaus.jProt.sequence;

import java.util.ArrayList;

/**
* A DNASequence is a String of single letter identifier specifying the sequence of the bases of a
* strand of DNA
* <p>
* Valid characters are only G A T C.
* @since 0.7.0
*/
public class DNASequence extends Sequence {

  /**
  * Construct a DNA Sequence. The sequences must contain only the characters G A T C.
  * @param seq a String containing only valid DNA base characters
  */
  public DNASequence(String seq){
    super(seq);
    validateSequence(seq);
  }

  /**
  * Calculate and return an alignment of the two DNA sequences.
  * @param other the other sequence to align to this one
  * @return an object of type Alignment which can be queried to get the results of this alignment
  */
  public Alignment align(DNASequence other){
    return SequenceAligner.align(this, other, "DNA.txt");
  }

  /**
  * Verify that the seq passed in is a valid DNA sequence.
  *
  * A valid DNA sequence can only contain G A T C
  *
  * @param seq the sequence to check
  * @throws IllegalArgumentException if the sequence is invalid
  */
  private void validateSequence(String seq){
    ArrayList<Character> validCharacters = new ArrayList<Character>();
    validCharacters.add('G');
    validCharacters.add('A');
    validCharacters.add('T');
    validCharacters.add('C');

    for(Character nucleotide : seq.toCharArray()){
      if(!validCharacters.contains(nucleotide)){
        throw new IllegalArgumentException("Seq invalid. " + nucleotide + " not a valid DNA base."
        + "\n only G A T C allowed.");
      }
    }
  }
}
