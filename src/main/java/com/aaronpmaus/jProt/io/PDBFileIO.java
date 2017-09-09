package com.aaronpmaus.jProt.io;

import com.aaronpmaus.jProt.protein.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Collection;
import java.util.HashSet;

import java.io.FileInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;

/**
* Provides the ability to read in and write out PDB Files. Every protein has its own
* instance of PDBFileIO, even proteins that aren't read in from a PDB, but are built from
* sequence data alone. Those proteins will need the ability to write themselves out to file.
* <p>
* Example Usage:
* <p>
* To build a protein from a pdb file:
* {@code Protein prot = new Protein(pdbInputStream, proteinName);}
* To build a protein from sequence data:
* {@code Protein prot = new Protein(new FastaFileIO(fastaInputStream), proteinName);}
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class PDBFileIO{
  private ArrayList<String> fileLines;
  private ArrayList<SSBondRecord> ssBonds;
  private ArrayList<AtomRecord> atomRecords;
  private HashMap<String,ArrayList<AtomRecord>> defaultResidues;
  private InputStream pdbInputStream;

  public PDBFileIO(InputStream pdbInputStream){
    fileLines = new ArrayList<String>();
    ssBonds = new ArrayList<SSBondRecord>();
    atomRecords = new ArrayList<AtomRecord>();
    defaultResidues = readInDefaultResidues();
    this.pdbInputStream = pdbInputStream;
    readInPDBFile(pdbInputStream);
  }

  /**
  * Read in a PDB file and return a Protein
  *
  * @param fileName the name of the PDB file
  * @throws FileNotFoundException if the file can not be read from
  */
  private void readInPDBFile(String fileName) throws FileNotFoundException{
    readInPDBFile(new FileInputStream(new File(fileName)));
  }

  /**
  * Read in a PDB file and return a Protein
  *
  * @param inputStream the inputStream to read from
  * @param fileName the name of the PDB file
  */
  private void readInPDBFile(InputStream inputStream){
    Scanner in = new Scanner(inputStream);

    while(in.hasNextLine()){
      String line = in.nextLine();
      if(line.length() < 80){
        line = padRight(line);
      }
      String recordName = line.substring(0,6).trim();
      switch(recordName){
        case "SSBOND":
          ssBonds.add(parseSSBondLine(line));
          break;
        case "ATOM":
          atomRecords.add(parseAtomLine(line));
          break;
        case "TER":
          break;
        case "END":
          break;
      }
    }
    //return buildProtein(fileBase);
  }

  /*
  * Build HashMap of the lists of AtomRecords for each residue in a chain.
  * the key is the residueID, the value is an ArrayList holding all the atomRecords for
  * that residue.
  */
  public HashMap<Integer, ArrayList<AtomRecord>> getResidueRecordsLists(
      Collection<AtomRecord> atomRecords){

    HashMap<Integer, ArrayList<AtomRecord>> residueRecordsLists =
        new HashMap<Integer, ArrayList<AtomRecord>>();
    // add all atom records to their respective residue atomRecords lists.
    // At this point, the residueRecordsLists is empty. After this loop, it will
    // contain an arrayList of atomRecords for every residueID.
    for(AtomRecord rec : atomRecords){
      int resID = rec.getResSeq();
      // if the arraylist for the residue already exists, add the record to it.
      if(residueRecordsLists.containsKey(resID)){
        residueRecordsLists.get(resID).add(rec);
      } else { //otherwise, build the arraylist and add the record to it.
        residueRecordsLists.put(resID, new ArrayList<AtomRecord>());
        residueRecordsLists.get(resID).add(rec);
      }
    }
    return residueRecordsLists;
  }

  /**
  * Build HashMap of the lists of AtomRecords for each of the default residues.
  *
  * The key is the residue name, the value is an ArrayList holding all the atomRecords for
  * that residue.
  */
  private HashMap<String, ArrayList<AtomRecord>> buildDefaultResiduesLists(
      ArrayList<AtomRecord> list){

    HashMap<String, ArrayList<AtomRecord>> defaultResiduesLists =
        new HashMap<String, ArrayList<AtomRecord>>();

    for(AtomRecord rec : list){
      String resThreeLetterID = rec.getResName();

      if(defaultResiduesLists.containsKey(resThreeLetterID)){
        defaultResiduesLists.get(resThreeLetterID).add(rec);
      } else { //otherwise, build the arraylist and add the record to it.
        defaultResiduesLists.put(resThreeLetterID, new ArrayList<AtomRecord>());
        defaultResiduesLists.get(resThreeLetterID).add(rec);
      }
    }
    return defaultResiduesLists;
  }

  /**
  * Return all the atom records for a defaultResidue with all atomSerialNumber resID fields
  * set to -1.
  * @param residueThreeLetterID a three letter ID of a residue in all caps.
  * @return an ArrayList containing all the atomRecords for that residue. includes hydrogens.
  * @throws IllegalArgumentException if the three letter ID is not a valid ID.
  */
  public ArrayList<AtomRecord> getDefaultResidueRecords(String residueThreeLetterID){
    if(defaultResidues.containsKey(residueThreeLetterID.toUpperCase())){
      return defaultResidues.get(residueThreeLetterID.toUpperCase());
    }
    throw new IllegalArgumentException(residueThreeLetterID + " not a valid residue name.");
  }


  /**
  * From the list of all AtomRecords in the pdb, return a list of those with the given chainID.
  * @param chainID a chainID that must be present in the PDB
  * @return a collection of all the AtomRecords that contain that chainID
  */
  public Collection<AtomRecord> getAtomRecordsInChain(String chainID){
    ArrayList<AtomRecord> atoms = new ArrayList<AtomRecord>();
    for(AtomRecord rec : this.atomRecords){
      if(rec.getChainID().equals(chainID)){
        atoms.add(rec);
      }
    }
    return atoms;
  }

  /**
  * Get a list all the chainIDs in this PDB.
  * @return a list of all the chainIDs in this PDB.
  */
  public Collection<String> getListOfChainIDs(){
    HashSet<String> chainIDs = new HashSet<String>();
    for(AtomRecord rec : this.atomRecords){
      chainIDs.add(rec.getChainID());
    }
    return chainIDs;
  }

  private HashMap<String, ArrayList<AtomRecord>> readInDefaultResidues(){
    InputStream stream = PDBFileIO.class.getResourceAsStream("residues.pdb");
    Scanner in = new Scanner(stream);
    ArrayList<AtomRecord> allResAtomRecords = new ArrayList<AtomRecord>();
    while(in.hasNextLine()){
      String line = in.nextLine();
      if(line.length() < 80){
        line = padRight(line);
      }
      String recordName = line.substring(0,6).trim();
      allResAtomRecords.add(parseAtomLine(line));
    }

    return buildDefaultResiduesLists(allResAtomRecords);
  }

  private void addDisulfideBonds(Protein protein){
    for(SSBondRecord ssBondRec: ssBonds){
      String chainID1 = ssBondRec.getChainID1();
      String chainID2 = ssBondRec.getChainID2();
      int resID1 = ssBondRec.getResID1();
      int resID2 = ssBondRec.getResID2();
      protein.addDisulfideBond(chainID1, resID1, chainID2, resID2);
    }
  }

  private String padRight(String line){
    return String.format("%1$-" + (80-line.length()) + "s", line);
  }

  private AtomRecord parseAtomLine(String line){
    int serial = Integer.parseInt(line.substring(6,11).trim());
    String atomName = line.substring(12,16).trim();
    String altLoc = line.substring(16,17).trim();
    String resName = line.substring(17,20).trim();
    String chainID = line.substring(21,22).trim();
    int resSeq = Integer.parseInt(line.substring(22,26).trim());
    String iCode = line.substring(26,27).trim(); // code for insertion of residues
    String x = line.substring(30,38).trim();
    String y = line.substring(38,46).trim();
    String z = line.substring(46,54).trim();

    // make sure there is a value for occupancy. if so, assign.
    // else, default value of -1.0.
    String occupancyString = line.substring(54,60).trim();
    double occupancy = -1.0;
    if(!occupancyString.equals("")){
      occupancy = Double.parseDouble(occupancyString);
    }

    // make sure there is a value for temp factor. if so, assign.
    // else, default value of -1.0.
    String tempFactorString = line.substring(60,66).trim();
    double tempFactor = -1.0;
    if(!tempFactorString.equals("")){
      tempFactor = Double.parseDouble(tempFactorString);
    }
    String element = line.substring(76,78).trim();

    // make sure there is a value for charge. if so, parse and assign.
    // else default charge of 0.
    // format for charge is 1+ or 2+ or 1-.
    String chargeString = line.substring(78,80).trim();
    double charge = 0;
    if(!chargeString.equals("")){
      charge = Double.parseDouble(chargeString.substring(0,1));
      if(chargeString.substring(1,2).equals("-")){
        charge *= -1;
      }
    }
    return new AtomRecord(serial, atomName, altLoc, resName, chainID, resSeq, iCode,
                          x, y, z, occupancy, tempFactor, element, charge);
  }

  private SSBondRecord parseSSBondLine(String line){
    String chainID1 = line.substring(15,16).trim();
    int resID1 = Integer.parseInt(line.substring(17,21).trim());
    String chainID2 = line.substring(29,30).trim();
    int resID2 = Integer.parseInt(line.substring(31,35).trim());
    return new SSBondRecord(chainID1, resID1, chainID2, resID2);
  }

}
