package com.aaronpmaus.jProt.sequence;

import com.aaronpmaus.jProt.protein.*;
import java.util.ArrayList;
import java.util.Scanner;

import java.io.InputStream;

/**
* A SequenceAligner can be used to align Protein, DNA, or RNA sequences.
* <p>
* It is designed as a class with only static methods. An individual alignment is completely
* self-contained so there need be no state to this class.
* <p>
* It uses the Needleman-Wuncsh algorithm with an affine gap penalty. The default values for the gap
* penalty are -10 to start a gap and -2 to extend it. The similarity scores matchMatrix used in the
* BLOSUM62 matrix.
* <p>
* It is intended that sequences are aligned using the align method from the subclasses
* of Sequence. For example, to perform a sequence alignment of two protein sequences:
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
* @see com.aaronpmaus.jProt.sequence.ProteinSequence
* @see com.aaronpmaus.jProt.sequence.Alignment
* @version 0.6.0
* @since 0.6.0
*/
public class SequenceAligner{
  private static int gapExtendPenalty = -2; //-2
  private static int gapStartPenalty = -10; //-10

  /**
  * Calculate and return an alignment of the two sequences.
  *
  * @param seq1 one of the sequences to align
  * @param seq2 the other sequence to align
  * @param matrixName the name of the Scoring Matrix to use, must be either "BLOSUM62.txt",
  * "DNA.txt", or "RNA.txt"
  * @return an object of type Alignment which can be queried to get the results of this alignment
  */
  public static Alignment align(Sequence seq1, Sequence seq2, String matrixName){
    ScoringMatrix scoringMatrix = new ScoringMatrix(matrixName);
    return align(seq1, seq2, scoringMatrix);
  }

  /**
  * A private helper method to perform the alignment given the two sequences and a scoringMatrix
  * for matches.
  *
  * This method uses the Needleman-Wuncsh algorithm. It's implementation was derived from
  * information on Carl Kingsford's Gap Penalties slides
  * {@link http://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf}.
  *
  * @param upSeq the sequence that goes across the top of the dynamically generated matrices
  * @param leftSeq the sequence that goes down the left side of the dynamically generated matrices
  * @param scoringMatrix a matrix containing the values of the values of the scores of alignment
  *   matches.
  * @return an array of String where the 0th String is the alignment of upSeq and the 1st String
  *     is the alignment of leftSeq
  */
  private static Alignment align(Sequence seq1, Sequence seq2, ScoringMatrix scoringMatrix){
    String upSeq = seq1.getSequenceString();
    String leftSeq = seq2.getSequenceString();
    if(upSeq.length() == 0 || leftSeq.length() == 0){
      throw new IllegalArgumentException("Sequences must not be empty.\n"
          + "seq1: |" + upSeq + "|\n"
          + "seq2: |" + leftSeq + "|\n");
    }
    // seq2 goes across the top and so determines the number of columns.
    // seq1 goes down the side and likewise determines the number of rows.
    int numRows = leftSeq.length() + 1;
    int numCols = upSeq.length() + 1;

    // These matrices will be used by several helper methods below and could be made instance
    // variables to save from having to pass them around, but since an instance of this class is
    // designed to be able to be used to align many sequence pairs, the working data of an
    // individual alignment should not be part of the state of this class
    Cell[][] matchMatrix = new Cell[numRows][numCols];
    Cell[][] xGapMatrix = new Cell[numRows][numCols];
    Cell[][] yGapMatrix = new Cell[numRows][numCols];

    // initialize all Cells in first row and column;
    initializeScoresMatrices(leftSeq, upSeq, matchMatrix, xGapMatrix, yGapMatrix);
    // fill in values for all remaining scores
    calculateMatrixValues(scoringMatrix, leftSeq, upSeq, matchMatrix, xGapMatrix, yGapMatrix);

    // the traceback needs to start at the bottom right cell with the largest value. Compare
    // the bottom right cell of all three matrices to find the largest.
    Cell max = matchMatrix[numRows-1][numCols-1];
    Cell xBottomRightCell = xGapMatrix[numRows-1][numCols-1];
    Cell yBottomRightCell = yGapMatrix[numRows-1][numCols-1];
    if(xBottomRightCell.getValue() > max.getValue()){
      max = xBottomRightCell;
    }
    if(yBottomRightCell.getValue() > max.getValue()){
      max = yBottomRightCell;
    }
    String[] alignments = traceBackOptimalAlignments(max);
    boolean[][] alignmentMasks = getSequenceMatchMasks(alignments[0], alignments[1]);
    return new Alignment(seq1, alignments[0], alignmentMasks[0],
                         seq2, alignments[1], alignmentMasks[1],
                         max.getValue());
  }

  /**
  * A private helper method to initialize the dynamically generated matrices.
  *
  * There are 3: the Matches Matrix, the X Gap Matrix, and the Y Gap Matrix. The Matches Matrix
  * holds the sequence scores and pointers for sequence matches. The X Gap Matrix holds scores and
  * pointers for gaps in the left sequence. The Y Gap Matrix holds scores and pointers for gaps in
  * the up Matrix. Every cell in the matrices has one or more pointers to other cells that were
  * used to calculate it's value.
  *
  * This helper method fills in the first row and column of each of these matrices.
  *
  * For an explanation of this matrices, see
  * {@link https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf}.
  *
  * The top left cell of each matrix should contain the value 0 and map to no sequence members.
  * The traceback will arrive at one of these, and at that point is complete.
  *
  * For the Matches Matrix, each of the remaining cells in the first row and column get a value
  * approximating -INFINITY because they represent the score of the best alignment of 0 characters
  * of one of the sequences and i characters of the other sequence that end in a match. -INFINITY
  * because this is impossibe. To avoid int overflow, Integer.MIN_VALUE/2 is used.
  *
  * For the X Gap Matrix, the values in the first row are calculated as gap_start + j*gap_extend.
  * where j is the column index. This is because the first row indicates the best alignment of
  * 0 elements of the left sequence to j elements of the top sequence, that is, a gap of length
  * j in the left sequence. The values of the first column are -INFINITY because they represent
  * a situation which is impossible, an alignment of 0 characters of upSeq with i characters of
  * leftSeq that ends in a gap.
  *
  * Likewise, the first row and column values of the Y Gap Matrix are like the X, but reversed
  * for similar logic.
  *
  * @param leftSeq the sequence that goes down the left side of the scores matrices
  * @param upSeq the sequence that goes across the top of the scores matrices
  */
  private static void initializeScoresMatrices(String leftSeq, String upSeq,
                                        Cell[][] matchMatrix,
                                        Cell[][] xGapMatrix,
                                        Cell[][] yGapMatrix){
    matchMatrix[0][0] = new Cell(0,null,null,1); // top left corner is 0 with no traceback
    xGapMatrix[0][0] = new Cell(0,null,null,2); // the 4th parameter indicates which matrix this
    yGapMatrix[0][0] = new Cell(0,null,null,3); // cell is in.
    // initialize all values on the first row
    int numCols = matchMatrix[0].length;
    for(int j = 1; j < numCols; j++){
      // Fill in first row of matchMatrix
      int value = Integer.MIN_VALUE/2; //sufficiently low number to approximate -INFINITY
      Cell cell  = new Cell(value, leftSeq.charAt(0), upSeq.charAt(j-1), 1);
      cell.addLeft(matchMatrix[0][j-1]);
      matchMatrix[0][j] = cell;
      // Fill in first row of xGapMatrix
      value = getGapStartPenalty() + (j) * getGapExtendPenalty();
      cell  = new Cell(value, leftSeq.charAt(0), upSeq.charAt(j-1), 2);
      cell.addLeft(xGapMatrix[0][j-1]);
      xGapMatrix[0][j] = cell;
      // Fill in first row of yGapMatrix
      value = Integer.MIN_VALUE/2;
      cell  = new Cell(value, leftSeq.charAt(0), upSeq.charAt(j-1), 3);
      cell.addLeft(yGapMatrix[0][j-1]);
      yGapMatrix[0][j] = cell;

    }
    // initialize all values on the first column
    int numRows = matchMatrix.length;
    for(int i = 1; i < numRows; i++){
      // Fill in first column of matchMatrix
      int value = Integer.MIN_VALUE/2;
      Cell cell  = new Cell(value, leftSeq.charAt(i-1), upSeq.charAt(0), 1);
      cell.addUp(matchMatrix[i-1][0]);
      matchMatrix[i][0] = cell;
      // Fill in first column of xGapMatrix
      value = Integer.MIN_VALUE/2;
      cell = new Cell(value, leftSeq.charAt(i-1), upSeq.charAt(0), 2);
      cell.addUp(xGapMatrix[i-1][0]);
      xGapMatrix[i][0] = cell;
      // Fill in first column of yGapMatrix
      value = getGapStartPenalty() + (i) * getGapExtendPenalty();
      cell  = new Cell(value, leftSeq.charAt(i-1), upSeq.charAt(0), 3);
      cell.addUp(yGapMatrix[i-1][0]);
      yGapMatrix[i][0] = cell;
    }
  }

  /**
  * Private helper method to loop over the matrices, dynamically calculating the values.
  *
  * @param ScoringMatrix an Object that allows match scores to be queried.
  * @param leftSeq the sequence that goes down the left side of the scores matrices
  * @param upSeq the sequence that goes across the top of the scores matrices
  */
  private static void calculateMatrixValues(ScoringMatrix scoringMatrix,
                                        String leftSeq, String upSeq,
                                        Cell[][] matchMatrix,
                                        Cell[][] xGapMatrix,
                                        Cell[][] yGapMatrix){
    int numRows = matchMatrix.length;
    int numCols = matchMatrix[0].length;
    for(int i = 1; i < numRows; i++){
      for(int j = 1; j < numCols; j++){
        // Calculate xGapMatrix value
        calculateXGapValue(i, j, leftSeq.charAt(i-1), upSeq.charAt(j-1), matchMatrix, xGapMatrix, yGapMatrix);
        // Calculate yGapMatrix value
        calculateYGapValue(i, j, leftSeq.charAt(i-1), upSeq.charAt(j-1), matchMatrix, xGapMatrix, yGapMatrix);
        // Calculate matchMatrix value
        calculateMatchValue(i, j, leftSeq.charAt(i-1), upSeq.charAt(j-1), scoringMatrix, matchMatrix, xGapMatrix, yGapMatrix);
      }
    }
  }

  /**
  * Calculate an individual value for the Matches Matrix.
  *
  * M[i][j] = scoringMatrix[i][j] + MAX(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1])
  * This calculates the score of that cell. A pointer are then added from this cell to a cell
  * that was used to calculate the Max value. For example, if both M[i-1][j-1] resulted in the max
  * then this cell will have a pointer back to it.
  * Likewise, if X[i-1][j-1] resulted in the max score, but not the corresponding values in M or Y
  * then this cell will only have one pointer back to it.
  *
  * @param i the row index of the cell that is being calculated
  * @param j the col index of the cell that is being calculated
  * @param leftSeqChar the character in the left sequence at this row
  * @param upSeqChar the character in the up sequence at this col
  * @param scoringMatrix the matrix to use to get the score of a match between leftSeqChar and
  *   upSeqChar
  */
  private static void calculateMatchValue(int i, int j,
                                   Character leftSeqChar, Character upSeqChar,
                                   ScoringMatrix scoringMatrix,
                                   Cell[][] matchMatrix,
                                   Cell[][] xGapMatrix,
                                   Cell[][] yGapMatrix){
    // calculate score from diagonal, up, and left
    // match between x and y
    int matchScore = scoringMatrix.getSimilarityScore(leftSeqChar, upSeqChar);
    int matchVal = matchMatrix[i-1][j-1].getValue() + matchScore;
    // gap in y
    int yVal = yGapMatrix[i-1][j-1].getValue() + matchScore;
    // gap in x
    int xVal = xGapMatrix[i-1][j-1].getValue() + matchScore;
    // save the max and construct the cell
    int max = max(matchVal, yVal, xVal);
    Cell cell = new Cell(max, leftSeqChar, upSeqChar, 1);
    // for every direction that is equal to the max, add a pointer to it
    if(matchVal == max){
      cell.addDiag(matchMatrix[i-1][j-1]);
    } else if(yVal == max){
      cell.addDiag(yGapMatrix[i-1][j-1]);
    } else if(xVal == max){
      cell.addDiag(xGapMatrix[i-1][j-1]);
    }
    // save the cell in the matchMatrix matrix
    matchMatrix[i][j] = cell;

  }

  /**
  * Calculate an individual value for the X Gap Matrix.
  *
  *                          Gap in X - Same row previous column
  *              { gap_start + gap_extend + M[i][j-1]
  * X[i][j] = MAX{ gap_extend + X[i][j-1]
  *              { gap_start + gap_extend + Y[i][j-1]
  *
  * The pointers are likewise calculated as the match matrix above
  *
  * @param i the row index of the cell that is being calculated
  * @param j the col index of the cell that is being calculated
  * @param leftSeqChar the character in the left sequence at this row
  * @param upSeqChar the character in the up sequence at this col
  */
  private static void calculateXGapValue(int i, int j, Character leftSeqChar, Character upSeqChar,
                                  Cell[][] matchMatrix, Cell[][] xGapMatrix, Cell[][] yGapMatrix){
    int matchVal = getGapStartPenalty() + getGapExtendPenalty() + matchMatrix[i][j-1].getValue();
    int xVal = getGapExtendPenalty() + xGapMatrix[i][j-1].getValue();
    int yVal = getGapStartPenalty() + getGapExtendPenalty() + yGapMatrix[i][j-1].getValue();
    int max = max(matchVal, xVal, yVal);
    Cell cell = new Cell(max, leftSeqChar, upSeqChar, 2);
    if(matchVal == max){
      cell.addLeft(matchMatrix[i][j-1]);
    } else if(xVal == max){
      cell.addLeft(xGapMatrix[i][j-1]);
    } else if(yVal == max){
      cell.addLeft(yGapMatrix[i][j-1]);
    }
    xGapMatrix[i][j] = cell;
  }

  /**
  * Calculate an individual value for the Y Gap Matrix.
  *
  *                          Gap in Y - Same column previous row
  *              { gap_start + gap_extend + M[i-1][j]
  * Y[i][j] = MAX{ gap_start + gap_extend + X[i-1][j]
  *              { gap_extend + Y[i-1][j]
  *
  * The pointers are likewise calculated as the match matrix above
  *
  * @param i the row index of the cell that is being calculated
  * @param j the col index of the cell that is being calculated
  * @param leftSeqChar the character in the left sequence at this row
  * @param upSeqChar the character in the up sequence at this col
  */
  private static void calculateYGapValue(int i, int j, Character leftSeqChar, Character upSeqChar,
                                  Cell[][] matchMatrix, Cell[][] xGapMatrix, Cell[][] yGapMatrix){
    int matchVal = getGapStartPenalty() + getGapExtendPenalty() + matchMatrix[i-1][j].getValue();
    int xVal = getGapStartPenalty() + getGapExtendPenalty() + xGapMatrix[i-1][j].getValue();
    int yVal = getGapExtendPenalty() + yGapMatrix[i-1][j].getValue();
    int max = max(matchVal, xVal, yVal);
    Cell cell = new Cell(max, leftSeqChar, upSeqChar, 3);
    if(matchVal == max){
      cell.addUp(matchMatrix[i-1][j]);
    } else if(xVal == max){
      cell.addUp(xGapMatrix[i-1][j]);
    } else if(yVal == max){
      cell.addUp(yGapMatrix[i-1][j]);
    }
    yGapMatrix[i][j] = cell;
  }

  /**
  * Trace back the optimum alignment starting at the start cell.
  *
  * For each cell, until it reaches a root, append the proper characters for the pointer
  * it has. Then follow that pointer back to the next previous cell.
  *
  * @param start the cell to start tracing the alignment from.
  * @return an array of String where the 0th String is the alignment of the left seq String
  *         is the alignment of right seq
  */
  private static String[] traceBackOptimalAlignments(Cell start){
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
        throw new NullPointerException("Cell is not Root but has no Pointers.");
      }
    }
    String[] alignment = new String[2];
    alignment[0] = new StringBuilder(seq1).reverse().toString();
    alignment[1] = new StringBuilder(seq2).reverse().toString();
    return alignment;
  }

  /**
  * Get a mask for each alignment string indicating the matches in its alignment to the other.
  *
  * The alignment strings are those produced by the alignment methods above and must be the
  * same length.
  *
  * The size of the mask is the number of elements in that alignment string exluding gaps. The
  * idea is to return a mask for the original sequence indicating which parts of it were able to
  * find a match in the other sequence.
  *
  * @param alignmentSequenceOne the top sequence in the alignment. Can contain element ids and gaps
  * @param alignmentSequenceTwo the bottom sequence in the alignment. Can contain element ids and
  *   gaps
  * @return a array containing 2 arrays of booleans. The first is the mask for sequenceOne, and
  *   the second the mask for sequenceTwo
  */
  public static boolean[][] getSequenceMatchMasks(String alignmentSequenceOne,
    String alignmentSequenceTwo){
    boolean[][] masks = new boolean[2][];
    masks[0] = getMask(alignmentSequenceOne, alignmentSequenceTwo);
    masks[1] = getMask(alignmentSequenceTwo, alignmentSequenceOne);
    return masks;
  }

  /**
  * From the alignment Strings seq1 and seq2, get a mask for all the elements in seq1 that
  * have a match in seq2. The elements are the original elements in the sequence that produced
  * seq1 alignment - that is, ignore all gaps in seq1
  *
  * given an alignment, for example:
  * AA-EYE
  * AAP-E-
  * we want to perform a comparison where there are matches.
  * The masks should be as long as their corresponding protein sequences, and include true
  * everywhere there is a match and false otherwise.
  *
  * In the case above, we are expecting
  * TTFTF
  * because AA have matches, E does not, Y does, and E again doesn't
  *
  * If we wanted to find the mask for the bottom sequence, that is for the case:
  * AAP-E-
  * AA-EYE
  * It would be:
  * TTFT
  *
  */
  private static boolean[] getMask(String seq1, String seq2){
    int seq1Length = 0; // the length of seq1 with no gaps
    for(char c : seq1.toCharArray()){
      if(c != '-'){
        seq1Length++;
      }
    }
    boolean[] mask = new boolean[seq1Length];
    int maskIndex = 0;
    for(int alignmentIndex = 0; alignmentIndex < seq1.length(); alignmentIndex++){
      char seq1Char = seq1.charAt(alignmentIndex);
      char seq2Char = seq2.charAt(alignmentIndex);
      if((seq1.charAt(alignmentIndex) != '-')
          && (seq2.charAt(alignmentIndex) != '-')){
        mask[maskIndex] = true;
        maskIndex++;
      } else if((seq1.charAt(alignmentIndex) != '-')
          && (seq2.charAt(alignmentIndex) == '-')){
        mask[maskIndex] = false;
        maskIndex++;
      } else if(seq1.charAt(alignmentIndex) == '-'){
        // ignore gaps in prot1;
      }
    }
    return mask;
  }

  /**
  * Return the max of 3 ints.
  *
  * @param a first int
  * @param b second int
  * @param c third int
  * @return the mac of these three ints
  */
  private static int max(int a, int b, int c){
    int max = a;
    if(b > max){
      max = b;
    }
    if(c > max){
      max = c;
    }
    return max;
  }

  /**
  * The gap start penalty is the given to starting a new gap when using affine gap penalties
  *
  * @return the gap start penalty
  */
  public static int getGapStartPenalty(){
    return SequenceAligner.gapStartPenalty;
  }

  /**
  * The gap start penalty is the given to extending a gap by one space when using affine gap
  * penalties
  *
  * @return the gap extend penalty
  */
  public static int getGapExtendPenalty(){
    return SequenceAligner.gapExtendPenalty;
  }

  /**
  * The default value of the gap extend penalty is -2. I don't have informed recommendations
  * on reasonable values.
  *
  * @param penalty the value to set the gap extend penalty to
  */
  public static void setGapExtendPenalty(int penalty){
    if(penalty > 0){
      throw new IllegalArgumentException("Gap Extend Penalty must be less than or equal 0.");
    }
    SequenceAligner.gapExtendPenalty = penalty;
  }

  /**
  * The default value of the gap start penalty is -10. I don't have informed recommendations
  * on reasonable values.
  *
  * @param penalty the value to set the gap start penalty to
  */
  public static void setGapStartPenalty(int penalty){
    if(penalty > 0){
      throw new IllegalArgumentException("Gap Start Penalty must be less than or equal 0.");
    }
    SequenceAligner.gapStartPenalty = penalty;
  }

  /**
  * Print out the values in a matrix of Cells. This is used for debugging purposes.
  *
  * The default -INFINITY values make this table look ugly. If debugging, it may help to decrease
  * the -INFINITY value to value still lower than any possible value, but small enough in
  * magnitude for reasonable printing.
  *
  * @param mat the matrix of Cells to print
  * @param leftSeq the sequence that goes down the left side of the scores matrices
  * @param upSeq the sequence that goes across the top of the scores matrices
  */
  private static void printMatrix(Cell[][] mat, String leftSeq, String upSeq){
    System.out.printf("%12s","");
    for(Character c : upSeq.toCharArray()){
      System.out.printf("%6s", c);
    }
    System.out.println();
    for(int i = 0; i < mat.length; i++){
      if(i == 0){
        System.out.printf("%6s","");
      } else {
        System.out.printf("%6s", leftSeq.charAt(i-1));
      }
      for(int j = 0; j < mat[i].length; j++){
        System.out.printf("%6d",mat[i][j].getValue());
      }
      System.out.println();
    }
  }

  /**
  * This private inner class Cell allows for elements of the matrices to have their values
  * but also pointers to other cells so that the traceback can be performed.
  */
  private static class Cell {
    private Cell diag = null;
    private Cell up = null;
    private Cell left = null;
    private String leftId;
    private String upId;
    private int value;
    private int matrixID;

    /**
    * Construct a cell using its value, the sequence character of its row and column, and
    * an int that represents which matrix it is in (1 - matches, 2 - X gap, 3 - Y gap)
    *
    * @param value the value of this cell
    * @param leftId the character in the sequence at this row of the matrix
    * @param upId the character in the sequence at this col of the matrix
    * @param matrixId an int representing which matrix this cell is in: 1 - matches, 2 - X gap,
    *     3 - Y gap
    */
    public Cell(int value, Character leftId, Character upId, int matrixId){
      this.value = value;
      this.leftId = ""+leftId;
      this.upId = ""+upId;
      this.matrixID = matrixId;
    }

    public String getMatrixName(){
      if(this.matrixID == 1){
        return "Match";
      } else if(this.matrixID == 2){
        return "X";
      } else if(this.matrixID == 3){
        return "Y";
      }
      return null;
    }

    public int getValue(){
      return this.value;
    }

    /**
    * If this cell was assigned its value because of a value in a cell diagonal to it, then
    * this method allows this cell to remember that it came from that diagonal cell. This means
    * that in the optimal path, the sequence goes from that cell to this one.
    *
    * @param cell the cell to the diagonal of this one to remember.
    */
    public void addDiag(Cell cell){
      this.diag = cell;
    }

    /**
    * If this cell was assigned its value because of a value in a cell above it, then
    * this method allows this cell to remember that it came from the above cell. This means
    * that in the optimal path, the sequence goes from that cell to this one.
    *
    * @param cell the cell above this one to remember.
    */
    public void addUp(Cell cell){
      this.up = cell;
    }

    /**
    * If this cell was assigned its value because of a value in a cell below it, then
    * this method allows this cell to remember that it came from the below cell. This means
    * that in the optimal path, the sequence goes from that cell to this one.
    *
    * @param cell the cell below this one to remember.
    */
    public void addLeft(Cell cell){
      this.left = cell;
    }

    public String getLeftId(){
      return this.leftId;
    }

    public String getUpId(){
      return this.upId;
    }

    /**
    * A Cell is root if it has no pointers to other cells.
    *
    * @return true if this cell is root, false otherwise
    */
    public boolean isRoot(){
      if(this.numPointers() == 0){
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

  /**
  * This class provides the ability to read in a Matches Scores Matrix, such as any BLOSUM matrix,
  * and be queried given sequence characters to find out their match value.
  */
  private static class ScoringMatrix {
    private ArrayList<Character> ids;
    private int[][] scores;

    /**
    * Constructing a ScoringMatrix requires the name of the file holding the Scores Matrix values.
    * This file will be read in from the resources folder for this package in this project.
    *
    * @param filename the name of the file holding the scores values. This filename does not include
    * the path to the file because it must be within this project.
    */
    public ScoringMatrix(String filename){
      InputStream stream = SequenceAligner.class.getResourceAsStream(filename);
      Scanner in = new Scanner(stream);
      ids = new ArrayList<Character>();
      scores = readInMatrix(in);
    }

    private int[][] readInMatrix(Scanner in){
      int[][] scores = null;
      while(in.hasNext()){
        String line = in.nextLine();
        String[] tokens = (line.trim()).split("\\s+");

        if(line.charAt(0) == '#'){
          // this is a comment, do nothing
        } else if(line.charAt(0) == ' '){
          // this is the first row of the matrix, it contains all the element IDs
          // add them to the arrayList of ids.
          for(String id : tokens){
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

    /**
    * Given two characters, return the score of their match in the sequence alignment.
    *
    * Examples include (from BLOSUM62), P vs. F has a score of -4 because in an aligned sequence,
    * a phenylalanine is not expected to be where a proline is. P vs. P has a score of 7 because
    * that is expected.
    *
    * @param id1 a Character in a match from one of the sequences
    * @param id2 the other Character in the match from the other sequence
    * @return the score of a match of these two characters.
    */
    public int getSimilarityScore(Character id1, Character id2){
      int index1 = ids.indexOf(id1);
      int index2 = ids.indexOf(id2);
      return scores[index1][index2];
    }
  }

}
