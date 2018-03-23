package com.aaronpmaus.jProt.sequence;

import java.util.Iterator;
import java.util.LinkedList;

/**
* A biological sequence. A Sequence can specify either the bases that make up DNA or RNA or
* the Amino Acids that make up Proteins.
* <p>
* Sequences can be aligned:
* <p>
* {@code InputStream prot1Stream = new FileInputStream(new File(pathToProt1));}<br>
* {@code InputStream prot2Stream = new FileInputStream(new File(pathToProt2));}<br>
* {@code Protein prot1 = PDBFileIO.readInPDBFile(prot1Stream, "prot1");}<br>
* {@code Protein prot2 = PDBFileIO.readInPDBFile(prot2Stream, "prot2");}<br>
* <br>
* {@code ProteinSequence prot1Seq = prot1.getSequence();}<br>
* {@code ProteinSequence prot2Seq = prot2.getSequence();}<br>
* {@code // Align the two sequences}<br>
* {@code Alignment alignment = prot1Seq.align(prot2Seq);}<br>
* {@code String prot1Alignment = alignment.getAlignment(prot1Seq;}<br>
* {@code String prot1Alignment = alignment.getAlignment(prot2Seq;}<br>
* @see com.aaronpmaus.jProt.sequence.SequenceAligner
* @see com.aaronpmaus.jProt.sequence.Alignment
* @version 0.7.0
* @since 0.6.0
*/
public abstract class Sequence implements Iterable<Character>{
  String seq;

  public Sequence(String seq){
    this.seq = seq;
  }

  /**
  * @return the String representation of this sequence
  */
  public String getSequenceString(){
    return this.seq;
  }

  /**
  * @return the number of characters in the sequence
  */
  public int getLength(){
    return seq.length();
  }

  /**
  * Calculate and return an alignment of the two sequences.
  * @param other the other sequence to align to this one
  * @return an object of type Alignment which can be queried to get the results of this alignment
  */
  public abstract Alignment align(Sequence other);

  @Override
  public Iterator<Character> iterator(){
    LinkedList<Character> list = new LinkedList<Character>();
    for(Character c :  seq.toCharArray()){
      list.add(c);
    }
    return list.iterator();
  }

  @Override
  public String toString(){
    return this.seq;
  }

  @Override
  public int hashCode(){
    return toString().hashCode();
  }

  @Override
  public boolean equals(Object obj){
    if(obj instanceof Sequence){
      Sequence other = (Sequence)obj;
      return this.toString().equals(other.toString());
    }
    return false;
  }
}
