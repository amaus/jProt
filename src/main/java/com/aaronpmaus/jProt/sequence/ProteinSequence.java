package com.aaronpmaus.jProt.sequence;

import java.util.ArrayList;

/**
* A ProteinSequence is a String of single letter residue identifier specifying the sequence of a
* protein.
* <p>
* Valid residue IDs include the IDs for the standard 20 amino acids plus:
* - B: Aspartic Acid or Asparagine
* - Z: Glutamic Acid or Glutamine
* - X: Undetermined
* - *: Any
* <p>
* Protein Sequences can be aligned. To perform a sequence alignment:
* <p>
* {@code InputStream prot1Stream = new FileInputStream(new File(pathToProt1));}<br>
* {@code InputStream prot2Stream = new FileInputStream(new File(pathToProt2));}<br>
* {@code Protein prot1 = PDBFileIO.readInPDBFile(prot1Stream, "prot1");}<br>
* {@code Protein prot2 = PDBFileIO.readInPDBFile(prot2Stream, "prot2");}<br>
* <br>
* {@code ProteinSequence prot1Seq = prot1.getSequence();}<br>
* {@code ProteinSequence prot2Seq = prot2.getSequence();}<br>
* {@code Alignment alignment = prot1Seq.align(prot2Seq);}<br>
* {@code String prot1Alignment = alignment.getAlignment(prot1Seq;}<br>
* {@code String prot1Alignment = alignment.getAlignment(prot2Seq;}<br>
* @see com.aaronpmaus.jProt.sequence.Alignment
* @since 0.7.0
*/
public class ProteinSequence extends Sequence {

  /**
  * The sequence must not be empty and must contain only valid capital single letter IDs.
  * These are the standard IDs for the standard 20 amino acids plus:
  * - B: Aspartic Acid or Asparagine
  * - Z: Glutamic Acid or Glutamine
  * - X: Undetermined
  * - *: Any
  * Undetermined indicates undetermined via the crystallographic process. This itself offers
  * clues about which amino acid it may be.
  * @param seq a String containing only valid residue identifiers
  */
  public ProteinSequence(String seq){
    super(seq);
    validateSequence(seq);
  }

  @Override
  public Alignment align(Sequence other){
    return SequenceAligner.align(this, other, "BLOSUM62");
  }

  /**
  * Verify that the seq passed in is a valid protein sequence.
  *
  * The sequences must contain only valid capital single letter IDs.
  * These are the standard IDs for the standard 20 amino acids plus:
  * - B: Aspartic Acid or Asparagine
  * - Z: Glutamic Acid or Glutamine
  * - X: Undetermined
  * - *: Any
  * Undetermined indicates undetermined via the crystallographic process. This itself offers
  * clues about which amino acid it may be.
  *
  * @param seq the sequence to check
  * @throws IllegalArgumentException if the sequence is invalid
  */
  private void validateSequence(String seq){
    ArrayList<Character> validCharacters = new ArrayList<Character>();
    validCharacters.add('A');
    validCharacters.add('R');
    validCharacters.add('N');
    validCharacters.add('D');
    validCharacters.add('C');
    validCharacters.add('Q');
    validCharacters.add('E');
    validCharacters.add('G');
    validCharacters.add('H');
    validCharacters.add('I');
    validCharacters.add('L');
    validCharacters.add('K');
    validCharacters.add('M');
    validCharacters.add('F');
    validCharacters.add('P');
    validCharacters.add('S');
    validCharacters.add('T');
    validCharacters.add('W');
    validCharacters.add('Y');
    validCharacters.add('V');
    validCharacters.add('B');
    validCharacters.add('Z');
    validCharacters.add('X');
    validCharacters.add('*');

    for(Character res : seq.toCharArray()){
      if(!validCharacters.contains(res)){
        throw new IllegalArgumentException("Seq invalid. " + res + " not a valid residue.");
      }
    }
  }

}
