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
* Provides the ability to read in and write out PDB Files
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class PDBFileIO{
  private ArrayList<String> fileLines;
  private ArrayList<SSBondRecord> ssBonds;
  private ArrayList<AtomRecord> atomRecords;

  public PDBFileIO(){
    fileLines = new ArrayList<String>();
    ssBonds = new ArrayList<SSBondRecord>();
    atomRecords = new ArrayList<AtomRecord>();
  }

  /**
  * Read in a PDB file and return a Protein
  *
  * @param fileName the name of the PDB file
  * @throws FileNotFoundException if the file can not be read from
  * @return a Protein built from that file
  */
  public Protein readInPDBFile(String fileName) throws FileNotFoundException{
    return readInPDBFile(new FileInputStream(new File(fileName)), fileName);
  }

  /**
  * Read in a PDB file and return a Protein
  *
  * @param inputStream the inputStream to read from
  * @param fileName the name of the PDB file
  * @return a Protein built from that file
  */
  public Protein readInPDBFile(InputStream inputStream, String fileName){
    Scanner in = new Scanner(inputStream);
    String[] fileNameParts = fileName.split("\\.");
    String fileBase = fileNameParts[0];

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
    return buildProtein(fileBase);
  }

  private Atom constructAtom(AtomRecord rec){
    return new Atom(rec.getAtomName(), rec.getSerial(), rec.getOccupancy(),
                    rec.getTempFactor(), rec.getCharge(),
                    rec.getX(), rec.getY(), rec.getZ());
  }

  private Protein buildProtein(String fileBase){
    Protein protein = new Protein(fileBase);
    Collection<String> chains = getListOfChains();
    for(String chainID : chains){
      PolypeptideChain chain = new PolypeptideChain(chainID);
      // get a list of all atoms in that chain
      Collection<AtomRecord> atoms = getAtomsInChain(chainID);
      // build HashMap of residues
      HashMap<Integer, ArrayList<AtomRecord>> residues =
          new HashMap<Integer, ArrayList<AtomRecord>>();
      // add all atom records to their respective atom lists in the map.
      for(AtomRecord rec : atoms){
        // if the arraylist for the residue is already built, add the record to it.
        int resID = rec.getResSeq();
        if(residues.containsKey(resID)){
          residues.get(resID).add(rec);
        } else { //otherwise, build the arraylist and add the record to it.
          residues.put(resID, new ArrayList<AtomRecord>());
          residues.get(resID).add(rec);
        }
      }
      // for every residue atoms list in the map
      for(ArrayList<AtomRecord> residueAtomRecords : residues.values()){
        String resName = residueAtomRecords.get(0).getResName();
        int resSeq = residueAtomRecords.get(0).getResSeq();
        ArrayList<Atom> residueAtoms = new ArrayList<Atom>();
        for(AtomRecord rec : residueAtomRecords){
          residueAtoms.add(constructAtom(rec));
        }
        chain.addResidue(new Residue(resName, resSeq, residueAtoms));
      }
      protein.addChain(chain);
    }
    return protein;
  }

  private Collection<AtomRecord> getAtomsInChain(String chain){
    ArrayList<AtomRecord> atoms = new ArrayList<AtomRecord>();
    for(AtomRecord rec : this.atomRecords){
      if(rec.getChainID().equals(chain)){
        atoms.add(rec);
      }
    }
    return atoms;
  }

  private Collection<String> getListOfChains(){
    HashSet<String> chains = new HashSet<String>();
    for(AtomRecord rec : this.atomRecords){
      chains.add(rec.getChainID());
    }
    return chains;
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

  private String padRight(String line){
    return String.format("%1$-" + (80-line.length()) + "s", line);
  }

  /**
  * Private Inner Class SSBondRecord
  */
  private class SSBondRecord{
    private String chainID1;
    private int resID1;
    private String chainID2;
    private int resID2;

    public SSBondRecord(String chainID1, int resID1, String chainID2, int resID2){
      this.chainID1 = chainID1;
      this.resID1 = resID1;
      this.chainID2 = chainID2;
      this.resID2 = resID2;
    }
    public String getChainID1(){ return this.chainID1; }
    public String getChainID2(){ return this.chainID2; }
    public int getResID1(){ return this.resID1; }
    public int getResID2(){ return this.resID2; }

  }

  /**
  * Private Inner Class AtomRecord
  */
  private class AtomRecord{
    private int serial;
    private String atomName;
    private String altLoc;
    private String resName;
    private String chainID;
    private int resSeq;
    private String iCode;
    private String x;
    private String y;
    private String z;
    private double occupancy;
    private double tempFactor;
    private String element;
    private double charge;

    public AtomRecord(int serial, String atomName, String altLoc,
                      String resName, String chainID, int resSeq,
                      String iCode, String x, String y, String z,
                      double occupancy, double tempFactor,
                      String element, double charge){
      this.serial = serial;
      this.atomName = atomName;
      this.altLoc = altLoc;
      this.resName = resName;
      this.chainID = chainID;
      this.resSeq = resSeq;
      this.iCode = iCode;
      this.x = x;
      this.y = y;
      this.z = z;
      this.occupancy = occupancy;
      this.tempFactor = tempFactor;
      this.element = element;
      this.charge = charge;
    }
    public int getSerial(){ return this.serial; }
    public String getAtomName(){ return this.atomName; }
    public String getAltLoc(){ return this.altLoc; }
    public String getResName(){ return this.resName; }
    public String getChainID(){ return this.chainID; }
    public int getResSeq(){ return this.resSeq; }
    public String getICode(){ return this.iCode; }
    public String getX(){ return this.x; }
    public String getY(){ return this.y; }
    public String getZ(){ return this.z; }
    public double getOccupancy(){ return this.occupancy; }
    public double getTempFactor(){ return this.tempFactor; }
    public String getElement(){ return this.element; }
    public double getCharge(){ return this.charge; }

  }

}
