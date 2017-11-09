package com.aaronpmaus.jProt.sequence;

public class Alignment {
  private final Sequence seqOne;
  private final Sequence seqTwo;
  private final String seqOneAlignment;
  private final String seqTwoAlignment;
  private final double score;
  private final boolean[] seqOneAlignmentMask;
  boolean[] seqTwoAlignmentMask;

  public Alignment(Sequence seqOne, String seqOneAlignment, boolean[] seqOneAlignmentMask,
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

  public String getAlignment(Sequence seq){
    if(seqOne.equals(seq)){
      return this.seqOneAlignment;
    } else if(seqTwo.equals(seq)){
      return this.seqTwoAlignment;
    } else {
      throw new IllegalArgumentException("Seq not one of the sequences Aligment was built from");
    }
  }

  public boolean[] getAlignmentMask(Sequence seq){
    if(seqOne.equals(seq)){
      return this.seqOneAlignmentMask;
    } else if(seqTwo.equals(seq)){
      return this.seqTwoAlignmentMask;
    } else {
      throw new IllegalArgumentException("Seq not one of the sequences Aligment was built from");
    }
  }

  public double getScore(){
    return this.score;
  }
}
