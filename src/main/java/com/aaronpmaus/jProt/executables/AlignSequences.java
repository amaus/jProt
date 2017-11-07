package com.aaronpmaus.jProt.executables;

import com.aaronpmaus.jProt.io.CommandLineParser;
import com.aaronpmaus.jProt.io.PDBFileIO;
import com.aaronpmaus.jProt.io.FASTAFileIO;

import com.aaronpmaus.jProt.protein.Protein;
import com.aaronpmaus.jProt.sequence.ProteinSequence;
import com.aaronpmaus.jProt.sequence.Alignment;

import java.io.InputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;

import java.util.Scanner;

/**
*
* @since 0.7.0
*/
public class AlignSequences {
  private static boolean mol1FileProvided = false;
  private static boolean mol2FileProvided = false;
  private static String mol1FilePath;
  private static String mol2FilePath;
  private static boolean usePDBs = false;
  private static boolean useFASTAs = false;


  public static void main(String[] arguments){
    CommandLineParser args = new CommandLineParser(arguments);
    if(arguments.length == 0 || args.contains("-h")){
      InputStream stream = JProtMetrics.class.getResourceAsStream("AlignSequencesUsage.txt");
      Scanner in = new Scanner(stream);
      while(in.hasNextLine()){
        System.out.println(in.nextLine());
      }
      System.exit(1);
    }
    if(args.contains("--mol1-f")){
      mol1FilePath = args.getValue("--mol1-f");
      mol1FileProvided = true;
    }
    if(args.contains("--mol2-f")){
      mol2FilePath = args.getValue("--mol2-f");
      mol2FileProvided = true;
    }

    try {
      if(mol1FileProvided && mol2FileProvided){
        File mol1File = new File(mol1FilePath);
        File mol2File = new File(mol2FilePath);
        String mol1FileName = mol1File.getName();
        String mol2FileName = mol2File.getName();
        String[] fileOneTokens = mol1FileName.split("\\.");
        String[] fileTwoTokens = mol2FileName.split("\\.");
        String mol1Ext = fileOneTokens[fileOneTokens.length-1];
        String mol2Ext = fileTwoTokens[fileTwoTokens.length-1];
        int extLength = mol1Ext.length();
        String mol1Base = mol1FileName.substring(0,mol1FileName.length()-extLength+1);
        String mol2Base = mol2FileName.substring(0,mol2FileName.length()-extLength+1);

        if(mol1Ext.equals("pdb") && mol2Ext.equals("pdb")){
          usePDBs = true;
        } else if(mol1Ext.equals("fasta") && mol2Ext.equals("fasta")){
          useFASTAs = true;
        } else {
          System.out.println("Both files must either be PDBs with extension .pdb"
            + " or FASTAs with extension .fasta");
          System.out.println("You provided files: \n" + mol1FileName + "\n" + mol2FileName);
          System.exit(1);
        }

        ProteinSequence seq1 = null;
        ProteinSequence seq2 = null;
        if(usePDBs){
          Protein prot1 = new PDBFileIO().readInPDBFile(new FileInputStream(mol1FileName),mol1Base);
          Protein prot2 = new PDBFileIO().readInPDBFile(new FileInputStream(mol2FileName),mol2Base);
          seq1 = prot1.getSequence();
          seq2 = prot2.getSequence();
        } else if(useFASTAs){
          seq1 = new FASTAFileIO().readInFASTAFile(new FileInputStream(mol1FileName)).get(0);
          seq2 = new FASTAFileIO().readInFASTAFile(new FileInputStream(mol2FileName)).get(0);
        }
        Alignment alignment = seq1.align(seq2);
        System.out.println(alignment.getAlignment(seq1));
        System.out.println(alignment.getAlignment(seq2));
      } else {
        System.out.println("You must provide the two files, either PDB or FASTA.");
        System.exit(1);
      }

    } catch (FileNotFoundException e){
      System.out.println("Could not open required files. Check for existence.");
      System.exit(1);
    }
  }
}
