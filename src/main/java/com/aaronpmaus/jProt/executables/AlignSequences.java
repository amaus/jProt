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
* This program aligns two sequences and prints out the alignment.
*
* <pre>
* <code>
* {@literal Usage: AlignSequences [<options>] <mol1-f fname> <mol2-f fname>}
* options :
*     -h
*         Display the usage file.
*     --mol1-f fname
*         This file must either be a PDB file or a fasta file. If it is
*         a PDB file it must end in the .pdb extension and conform to
*         the PDB File Format Version 3.30. Otherwise, the file must be
*         a fasta file ending in the .fasta extension.
*     --mol2-f fname
*         The file for molecule two. The format is the same as mol1-f.
*     --mol1-chain chainID
*         The chain from mol1 to use in the alignment, capital letter
*         indicating the chain. This argument is only applicable if PDBs
*         are provided as input, otherwise it will be ignored. If
*         omitted, the entire sequence will be used.
*     --mol2-chain chainID
*         The chain from mol2 to use in the alignment, capital letter
*         indicating the chain. This argument is only applicable if PDBs
*         are provided as input, otherwise it will be ignored. If
*         omitted, the entire sequence will be used.
*     --mol1-entry entryNumber
*         A number indicating which entry to use in the alignment.
*         Possible values: [1,N], where N is the number of entries in
*         the FASTA File. If omitted the first entry will be used.
*     --mol2-entry entryNumber
*         A number indicating which entry to use in the alignment.
*         Possible values: [1,N], where N is the number of entries in
*         the FASTA File. If omitted the first entry will be used.
* </code>
* </pre>
*
* @since 0.7.0
*/
public class AlignSequences {
  private static boolean mol1FileProvided = false;
  private static boolean mol2FileProvided = false;
  private static boolean mol1ChainProvided = false;
  private static boolean mol2ChainProvided = false;
  private static boolean mol1EntryProvided = false;
  private static boolean mol2EntryProvided = false;

  private static String mol1FilePath;
  private static String mol2FilePath;
  private static String mol1Chain;
  private static String mol2Chain;
  private static int mol1Entry;
  private static int mol2Entry;
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
    if(args.contains("--mol1-chain")){
      mol1ChainProvided = true;
      mol1Chain = args.getValue("--mol1-chain");
    }
    if(args.contains("--mol2-chain")){
      mol2ChainProvided = true;
      mol2Chain = args.getValue("--mol2-chain");
    }
    if(args.contains("--mol1-entry")){
      mol1EntryProvided = true;
      mol1Entry = Integer.parseInt(args.getValue("--mol1-entry"));
    }
    if(args.contains("--mol2-entry")){
      mol2EntryProvided = true;
      mol2Entry = Integer.parseInt(args.getValue("--mol2-entry"));
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
          if(mol1ChainProvided){
            seq1 = prot1.getChain(mol1Chain).getSequence();
          } else {
            seq1 = prot1.getSequence();
          }
          if(mol2ChainProvided){
            seq2 = prot2.getChain(mol2Chain).getSequence();
          } else {
            seq2 = prot2.getSequence();
          }
        } else if(useFASTAs){

          if(mol1EntryProvided){
            try {
              FASTAFileIO io = new FASTAFileIO();
              io.readInFASTAFile(new FileInputStream(mol1FileName));
              seq1 = io.getSequence(mol1Entry-1);
            } catch (IllegalArgumentException e){
              System.out.println("mol1-entry must be between 1 and the number "
                  + "of entries in the fasta file. Provided: " + mol1Entry);
              System.exit(1);
            }
          } else {
            FASTAFileIO io = new FASTAFileIO();
            io.readInFASTAFile(new FileInputStream(mol1FileName));
            seq1 = io.getSequence(0);
          }

          if(mol2EntryProvided){ // use the entry specified
            try {
              FASTAFileIO io = new FASTAFileIO();
              io.readInFASTAFile(new FileInputStream(mol2FileName));
              seq2 = io.getSequence(mol2Entry-1);
            } catch (IllegalArgumentException e){
              System.out.println("mol2-entry must be between 1 and the number "
                  + "of entries in the fasta file. Provided: " + mol2Entry);
              System.exit(1);
            }
          } else { // Use the first entry
            FASTAFileIO io = new FASTAFileIO();
            io.readInFASTAFile(new FileInputStream(mol2FileName));
            seq2 = io.getSequence(0);
          }
        }
        Alignment alignment = seq1.align(seq2);
        System.out.println(alignment.getAlignment(seq1));
        System.out.println(alignment.getAlignment(seq2));
        System.out.println("Alignment Score: " + alignment.getScore());
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
