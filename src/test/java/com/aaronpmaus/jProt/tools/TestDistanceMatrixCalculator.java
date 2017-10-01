package com.aaronpmaus.jProt;

import com.aaronpmaus.jProt.protein.Protein;
import com.aaronpmaus.jProt.protein.PolypeptideChain;
import com.aaronpmaus.jProt.io.PDBFileIO;
import com.aaronpmaus.jProt.tools.DistanceMatrixCalculator;

import static org.junit.Assert.*;
import org.junit.Test;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Rule;
import org.junit.rules.ExpectedException;

import java.util.Collection;
import java.util.ArrayList;
import java.util.Scanner;
import java.math.BigDecimal;

import java.io.InputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

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

public class TestDistanceMatrixCalculator{
  Double[][] ropMatrix;
  Double[][] ropOddResMatrix;
  Protein rop;

  @Before
  public void setup() throws FileNotFoundException, IOException{
    InputStream stream = TestDistanceMatrixCalculator.class.getResourceAsStream("1ropDistanceMatrix.csv");
    ropMatrix = readInDistanceFile(stream);
    stream.close();

    stream = TestDistanceMatrixCalculator.class.getResourceAsStream("1rop_OddResiduesDistanceMatrix.csv");
    ropOddResMatrix = readInDistanceFile(stream);
    stream.close();

    stream = TestDistanceMatrixCalculator.class.getResourceAsStream("1rop.pdb");
    //PDBFileIO pdb = new PDBFileIO(stream);
    rop = new PDBFileIO().readInPDBFile(stream, "1rop");
    stream.close();
  }

  @Test
  public void testChainDistanceMatrix(){
    PolypeptideChain chain = rop.getChain("A");

    Double[][] calculatedMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(chain);
    verifyMatrix(calculatedMatrix, ropMatrix);
  }

  @Test
  public void testChainDistanceMatrixWithMask(){
    PolypeptideChain chain = rop.getChain("A");

    boolean[] mask = new boolean[chain.getNumResidues()];
    for(int i = 0; i < mask.length; i++){
      if(i%2 == 0){
        mask[i] = true;
      } else {
        mask[i] = false;
      }
    }

    Double[][] calculatedMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(chain, mask);
    verifyMatrix(calculatedMatrix, ropOddResMatrix);
  }

  @Test
  public void testProteinDistanceMatrix(){
    Double[][] calculatedMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(rop);
    verifyMatrix(calculatedMatrix, ropMatrix);
  }

  @Test
  public void testProteinDistanceMatrixWithMask(){
    boolean[] mask = new boolean[rop.getNumResidues()];
    for(int i = 0; i < mask.length; i++){
      if(i%2 == 0){
        mask[i] = true;
      } else {
        mask[i] = false;
      }
    }

    Double[][] calculatedMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(rop, mask);
    verifyMatrix(calculatedMatrix, ropOddResMatrix);
  }

  private void verifyMatrix(Double[][] calculatedMatrix, Double[][] referenceMatrix){
    boolean different = false;

    for(int i = 1; i < calculatedMatrix.length; i++){
      for(int j = i+1; j < calculatedMatrix[i].length; j++){
        if(Math.abs(referenceMatrix[i][j] - getNumInSigFigs(calculatedMatrix[i][j],4)) > 0.001){
          different = true;
        }
      }
    }
    assertFalse(different);
  }

  private Double[][] readInDistanceFile(InputStream stream) throws FileNotFoundException{
    Scanner fileReader = new Scanner(stream);
    fileReader.nextLine(); // throw away the first line. it's the residue IDs.

    ArrayList<String> dataLines = new ArrayList<String>();
    while(fileReader.hasNext()){
      dataLines.add(fileReader.nextLine());
    }
    int numResidues = dataLines.size();
    Double[][] data = new Double[numResidues][numResidues];
    for(int i = 0; i < numResidues; i++){
      for(int j = 0; j < numResidues; j++){
        data[i][j] = Double.NaN;
      }
    }
    // only use the values from the upper right hand side of the matrix.
    for(int i = 0; i < dataLines.size(); i++){
      String[] tokens = dataLines.get(i).split(",");
      for(int j = i+1; j < tokens.length; j++){
        data[i][j] = Double.parseDouble(tokens[j]);
      }
    }
    return data;
  }

  private static double getNumInSigFigs(double num, int significantFigures){
    BigDecimal bd = BigDecimal.valueOf(num);
    return new Double(String.format("%."+significantFigures+"G",bd));
  }

}
