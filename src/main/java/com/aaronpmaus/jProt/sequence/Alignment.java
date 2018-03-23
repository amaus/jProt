package com.aaronpmaus.jProt.sequence;

/**
* An Alignment holds the results of aligning two sequences, and provides methods for accessing those
* results.
* <p>
* It can not be constructed by clients outside of this package. Rather it is the object returned as
* the result of Sequence::align().
* @see com.aaronpmaus.jProt.sequence.Sequence
* @see com.aaronpmaus.jProt.sequence.SequenceAligner
*/
public class Alignment {
  private final Sequence seqOne;
  private final Sequence seqTwo;
  private final String seqOneAlignment;
  private final String seqTwoAlignment;
  private final double score;
  private final boolean[] seqOneAlignmentMask;
  boolean[] seqTwoAlignmentMask;

  Alignment(Sequence seqOne, String seqOneAlignment, boolean[] seqOneAlignmentMask,
      Sequence seqTwo, String seqTwoAlignment, boolean[] seqTwoAlignmentMask,
      double score){
    this.seqOne = seqOne;
    this.seqTwo = seqTwo;
    this.seqOneAlignment = seqOneAlignment;
    this.seqTwoAlignment = seqTwoAlignment;
    this.seqOneAlignmentMask = seqOneAlignmentMask;
    this.seqTwoAlignmentMask = seqTwoAlignmentMask;
    this.score = score;
  }

  /**
  * Return the alignment string corresponding to one of the input sequences.
  * @param seq a Sequence, one of the sequences input for the alignment
  * @return a String representing that sequences alignment to the other
  * @throws IllegalArgumentException if seq is not one of the input Sequences for the alignment
  */
  public String getAlignment(Sequence seq){
    if(seqOne.equals(seq)){
      return this.seqOneAlignment;
    } else if(seqTwo.equals(seq)){
      return this.seqTwoAlignment;
    } else {
      throw new IllegalArgumentException("Seq not one of the sequences Aligment was built from");
    }
  }

  /**
  * Return the alignment mask for one of the input sequences. This is a mask that indicates which
  * elements in the input sequence had a match in the other sequence.
  * @param seq a Sequence, one of the sequences input for the alignment
  * @return a boolean array where mask[i] is true if the element at index i in both sequences was a
  * match
  * @throws IllegalArgumentException if seq is not one of the input Sequences for the alignment
  */
  public boolean[] getAlignmentMask(Sequence seq){
    if(seqOne.equals(seq)){
      return this.seqOneAlignmentMask;
    } else if(seqTwo.equals(seq)){
      return this.seqTwoAlignmentMask;
    } else {
      throw new IllegalArgumentException("Seq not one of the sequences Aligment was built from");
    }
  }

  /**
  * @return the alignment score as calculated by the Needleman-Wuncsh alignment algorithm
  */
  public double getScore(){
    return this.score;
  }
}
