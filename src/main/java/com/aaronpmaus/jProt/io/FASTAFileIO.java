package com.aaronpmaus.jProt.io;

import com.aaronpmaus.jProt.sequence.ProteinSequence;

import java.io.InputStream;
import java.io.IOException;;
import java.util.Scanner;
import java.util.ArrayList;

/**
* Provides the ability to read in FASTA Files.
* <p>
* Example Usage:
* <p>
* To read in a FASTA File:
* <p>
* {@code FASTAFileIO fastaReader = new FASTAFileIO();} <br>
* {@code InputStream inputStream = new FileInputStream(new File(pathToFASTAFile));} <br>
* {@code ArrayList<ProteinSequence> sequences = fastaReader.readInFASTAFile(inputStream);} <br>
* {@code for(ProteinSequence seq : sequences) System.out.println(seq);} <br>
* @since 0.7.0
*/
public class FASTAFileIO{
  private ArrayList<ProteinSequence> sequences;
  private ArrayList<String> comments;
  private boolean fastaReadIn = false;

  /**
  * Construct an instance of FASTAFileIO. This instance can be used to read in a single FASTA File.
  */
  public FASTAFileIO(){
    sequences = new ArrayList<ProteinSequence>();
    comments = new ArrayList<String>();
  }

  /**
  * Read in a FASTA File and return an ArrayList of the ProteinSequences in it.
  * <p>
  * An instance of FASTAFileIO can only call this method once.
  * @param inputStream the inputStream to read from
  */
  public void readInFASTAFile(InputStream inputStream){
    if(this.fastaReadIn){
      throw new IllegalStateException("A FASTAFileIO Object can only read in a single FASTA File. "
      + " If you wish to read in a second FASTA file, you must instantiate another FASTAFileIO.");
    }
    this.fastaReadIn = true;
    Scanner in = new Scanner(inputStream);
    boolean lastLineComment = true;
    String seqString = "";
    String commentString = "";
    boolean readingFirstLine = true;
    while(in.hasNextLine()){
      String line = in.nextLine().trim();
      char firstChar = line.toCharArray()[0];
      if(firstChar == '>' || firstChar == ';'){
        // line is a comment
        if(!lastLineComment){
          // if the last line was not a comments, then build a new ProteinSequence out of the
          // seqString and add it to the sequences.
          this.sequences.add(new ProteinSequence(seqString));
          seqString = "";
          this.comments.add(line);
        } else if(lastLineComment){
          if(readingFirstLine){
            this.comments.add(line);
          } else {
            // if the last line was a comment, concatenate onto it in the comments ArrayList this
            // comment as well.
            int lastIndex = comments.size() - 1; //((comments.size() - 1) < 0  ?  0  :  comments.size() - 1);
            this.comments.set(lastIndex, comments.get(lastIndex) + "\n" + line);
          }
        }
        lastLineComment = true;
      } else {
        // if the line is not a comment, add it onto seqString
        seqString += line;
        lastLineComment = false;
      }
    }
    this.sequences.add(new ProteinSequence(seqString));
    try{
      in.close();
      inputStream.close();
    } catch(IOException e){

    }
    //return this.sequences;
  }

  /**
  * @return the number of sequences read in from the FASTA File
  */
  public int getNumSequences(){
    return this.sequences.size();
  }

  /**
  * @param index a number in the range [0, getNumSequences()-1]
  * @return a ProteinSequence
  * @throws IllegalArgumentException if index is not in range of [0,getNumSequences()-1]
  */
  public ProteinSequence getSequence(int index) throws IllegalArgumentException {
    if(index < 0 || index >= getNumSequences()){
      throw new IllegalArgumentException("FASTAFileIO::getSequence() - Index must be in range "
          + "[0, getNumSequences()-1], provided: " + index);
    }
    return this.sequences.get(index);
  }

  /**
  * @param index a number in the range [0, getNumSequences()-1]
  * @return a comment. Comments correspond to ProteinSequences by index. The comment at index
  * 0 corresponds to the ProteinSequence at index 0.
  * @throws IllegalArgumentException if index is not in range of [0,getNumSequences()-1]
  */
  public String getComment(int index){
    if(index < 0 || index >= getNumSequences()){
      throw new IllegalArgumentException("FASTAFileIO::getComment() - Index must be in range "
          + "[0, getNumSequences()-1], provided: " + index);
    }
    return this.comments.get(index);
  }
}
