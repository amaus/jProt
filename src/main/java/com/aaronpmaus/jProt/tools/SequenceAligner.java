package com.aaronpmaus.jProt.tools;

import com.aaronpmaus.jProt.protein.*;
import java.util.ArrayList;
import java.util.Scanner;

import java.io.InputStream;

/**
 * @version 0.10.0
 * @since 0.10.0
*/
public class SequenceAligner{
  private boolean debug = false;
  public String[] alignProteinSequences(String seq1, String seq2){
    validateProteinSequence(seq1);
    validateProteinSequence(seq2);

    if(debug) System.out.println("Sequences Validated. Reading in BLOSUM Matrix...");
    ScoringMatrix scoringMatrix = new ScoringMatrix("BLOSUM62.txt");
    if(debug) System.out.println("BLOSUM62 read in. aligning...");
    return align(seq1, seq2, scoringMatrix);
  }

  public String[] alignNucleotideSequences(String seq1, String seq2){
    validateNucleotideSequence(seq1);
    validateNucleotideSequence(seq2);

    if(debug) System.out.println("Sequences Validated. Reading in DNA Matrix...");
    ScoringMatrix scoringMatrix = new ScoringMatrix("DNA.txt");
    if(debug) System.out.println("DNA read in. aligning...");
    return align(seq1, seq2, scoringMatrix);
  }

  private String[] align(String upSeq, String leftSeq, ScoringMatrix scoringMatrix){
    // seq2 goes across the top and so determines the number of columns.
    // seq1 goes down the side and likewise determines the number of rows.
    int numRows = leftSeq.length() + 1;
    int numCols = upSeq.length() + 1;
    Cell[][] scores = new Cell[numRows][numCols];
    // initialize all Cells in first row and column;
    if(debug) System.out.println("Initializing Scores Matrix...");
    initializeScoresMatrix(scores, leftSeq, upSeq);
    if(debug) System.out.println("Scores Initialized. Calculating Matrix Values...");
    // fill in values for all remaining scores
    calculateMatrixValues(scores, scoringMatrix, leftSeq, upSeq);
    if(debug) System.out.println("Matrix values calculated. Tracing back optimal alignment...");

    return traceBackOptimalAlignments(scores[numRows-1][numCols-1]);
  }

  public void initializeScoresMatrix(Cell[][] scores, String leftSeq, String upSeq){
    scores[0][0] = new Cell(0,null,null); // top left corner is 0 with no traceback
    // initialize all values on the first row
    int numCols = scores[0].length;
    for(int j = 1; j < numCols; j++){
      int value = j * gapPenalty();
      Cell cell  = new Cell(value, leftSeq.charAt(0), upSeq.charAt(j-1));
      cell.addLeft(scores[0][j-1]);
      scores[0][j] = cell;
    }
    // initialize all values on the first column
    int numRows = scores.length;
    for(int i = 1; i < numRows; i++){
      int value = i * gapPenalty();
      Cell cell  = new Cell(value, leftSeq.charAt(i-1), upSeq.charAt(0));
      cell.addUp(scores[i-1][0]);
      scores[i][0] = cell;
    }
  }

  public void calculateMatrixValues(Cell[][] scores, ScoringMatrix scoringMatrix,
                                        String leftSeq, String upSeq){
    int numRows = scores.length;
    int numCols = scores[0].length;
    for(int i = 1; i < numRows; i++){
      for(int j = 1; j < numCols; j++){
        // calculate score from diagonal, up, and left
        // match between x and y
        int qDiag = scores[i-1][j-1].getValue()
            + scoringMatrix.getSimilarityScore(leftSeq.charAt(i-1), upSeq.charAt(j-1));
        // gap in y
        int qUp = scores[i-1][j].getValue() + gapPenalty();
        // gap in x
        int qLeft = scores[i][j-1].getValue() + gapPenalty();
        // save the max and construct the cell
        int max = max(qDiag, qUp, qLeft);
        Cell cell = new Cell(max, leftSeq.charAt(i-1), upSeq.charAt(j-1));
        // for every direction that is equal to the max, add a pointer to it
        if(qDiag == max){
          cell.addDiag(scores[i-1][j-1]);
        }
        if(qUp == max){
          cell.addUp(scores[i-1][j]);
        }
        if(qLeft == max){
          cell.addLeft(scores[i][j-1]);
        }
        // save the cell in the scores matrix
        scores[i][j] = cell;
      }
    }
  }

  public String[] traceBackOptimalAlignments(Cell start){
    String seq1 = "";
    String seq2 = "";
    Cell current = start;
    while(!current.isRoot()){
      if(current.hasDiag()){
        seq1 += current.getUpId();
        seq2 += current.getLeftId();
        current = current.getDiag();
      } else if(current.hasLeft()){
        seq1 += current.getUpId();
        seq2 += "-";
        current = current.getLeft();
      } else if(current.hasUp()){
        seq1 += "-";
        seq2 += current.getLeftId();
        current = current.getUp();
      } else {
        if(debug) System.out.println("Should NOT BE HERE");
        throw new NullPointerException("Cell is not Root but has no Pointers.");
      }
    }
    String[] alignment = new String[2];
    alignment[0] = new StringBuilder(seq1).reverse().toString();
    alignment[1] = new StringBuilder(seq2).reverse().toString();
    return alignment;
  }

  public int gapPenalty(){
    return -10;
  }

  private int max(int a, int b, int c){
    int max = a;
    if(b > max){
      max = b;
    }
    if(c > max){
      max = c;
    }
    return max;
  }

  private void validateProteinSequence(String seq){
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

  private void validateNucleotideSequence(String seq){
    ArrayList<Character> validCharacters = new ArrayList<Character>();
    validCharacters.add('G');
    validCharacters.add('A');
    validCharacters.add('T');
    validCharacters.add('C');

    for(Character nucleotide : seq.toCharArray()){
      if(!validCharacters.contains(nucleotide)){
        throw new IllegalArgumentException("Seq invalid. " + nucleotide + " not a valid nucleotide."
            + "\n only G A T C allowed.");
      }
    }
  }

  private class Cell {
    private Cell diag = null;
    private Cell up = null;
    private Cell left = null;
    private String leftId;
    private String upId;
    private int value;

    public Cell(int value, Character leftId, Character upId){
      this.value = value;
      this.leftId = ""+leftId;
      this.upId = ""+upId;
    }

    public int getValue(){
      return this.value;
    }

    public void addDiag(Cell cell){
      this.diag = cell;
    }

    public void addUp(Cell cell){
      this.up = cell;
    }

    public void addLeft(Cell cell){
      this.left = cell;
    }

    public String getLeftId(){
      return this.leftId;
    }

    public String getUpId(){
      return this.upId;
    }

    public boolean isRoot(){
      if(value == 0 && this.numPointers() == 0){
        return true;
      }
      return false;
    }

    public Cell getDiag(){
      return this.diag;
    }

    public Cell getLeft(){
      return this.left;
    }

    public Cell getUp(){
      return this.up;
    }

    public boolean hasDiag(){
      if(this.getDiag() == null){
        return false;
      }
      return true;
    }

    public boolean hasLeft(){
      if(this.getLeft() == null){
        return false;
      }
      return true;
    }

    public boolean hasUp(){
      if(this.getUp() == null){
        return false;
      }
      return true;
    }

    public int numPointers(){
      int num = 0;
      if(this.hasDiag()){
        num++;
      }
      if(this.hasUp()){
        num++;
      }
      if(this.hasLeft()){
        num++;
      }
      return num;
    }
  }

  private class ScoringMatrix {
    private ArrayList<Character> ids;
    private int[][] scores;

    public ScoringMatrix(String filename){
      if(debug) System.out.println("Building InputStream to read in BLOSUM file.");
      InputStream stream = SequenceAligner.class.getResourceAsStream(filename);
      if(debug) System.out.println("InputStream built. Constructing the Scanner.");
      Scanner in = new Scanner(stream);
      if(debug) System.out.println("Scanner pointing to BLOSUM file obtained. Reading in the matrix...");
      ids = new ArrayList<Character>();
      scores = readInMatrix(in);
      if(debug) {
        System.out.println("Matrix read in. Values: ");
        for(int i = 0; i < scores.length; i++){
          for(int j = 0; j < scores[i].length; j++){
            System.out.print(scores[i][j] + " ");
          }
          System.out.println();
        }
      }
    }

    private int[][] readInMatrix(Scanner in){
      int[][] scores = null;
      while(in.hasNext()){
        String line = in.nextLine();
        if(debug) System.out.println(line);
        String[] tokens = (line.trim()).split("\\s+");

        if(line.charAt(0) == '#'){
          // this is a comment, do nothing
        } else if(line.charAt(0) == ' '){
          // this is the first row of the matrix, it contains all the element IDs
          // add them to the arrayList of ids.
          for(String id : tokens){
            if(debug) System.out.println("|"+id+"|");
            ids.add(id.charAt(0));
          }
          scores = new int[tokens.length][tokens.length];
        } else {
          // this is a row containing similarity values
          int colIndex = 0;
          Character rowId = tokens[0].charAt(0);
          int rowIndex = ids.indexOf(rowId);
          for(int i = 1; i < tokens.length; i++){
            scores[rowIndex][colIndex] = Integer.parseInt(tokens[i]);
            colIndex++;
          }
        }
      }
      return scores;
    }

    public int getSimilarityScore(Character id1, Character id2){
      int index1 = ids.indexOf(id1);
      int index2 = ids.indexOf(id2);
      return scores[index1][index2];
    }
  }

}
