package com.aaronpmaus.jProt.io;

import com.aaronpmaus.jProt.protein.*;

import com.aaronpmaus.jMath.linearAlgebra.Vector3D;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Collection;
import java.util.HashSet;

import java.io.FileInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.IOException;

/**
* Provides the ability to read in and write out PDB Files. Every protein has its own instance of
* PDBFileIO. Even proteins that aren't read in from a PDB file, but are built from sequence data
* alone, have an instance of PDBFileIO because they will need the ability to write themselves out to
* file.
* <p>
* Example Usage:
* <p>
* To build a protein from a pdb file:
* <p>
* {@code InputStream in = new FileInputStream(new File(pathToPDBFile));}
* <br>
* {@code Protein prot = PDBFileIO.readInPDBFile(in, proteinName)}
* <p>
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class PDBFileIO{
  private ArrayList<String> fileLines;
  private ArrayList<SSBondRecord> ssBonds;
  private ArrayList<AtomRecord> atomRecords;
  private static HashMap<String,ArrayList<AtomRecord>> defaultResiduesRecords =
      readInDefaultResidues();
  private boolean proteinReadIn = false;

  public PDBFileIO(){
    fileLines = new ArrayList<String>();
    ssBonds = new ArrayList<SSBondRecord>();
    atomRecords = new ArrayList<AtomRecord>();
  }

  /**
  * Read in a PDB file and return a Protein. A PDBFileIO Object can only call this method once.
  * If you wish to read in multiple pdb file, you must instantiate a PDBFileIO object for each.
  * <p>
  * This method does not close the InputStream. It is the clients job to do so.
  * <p>
  * @param inputStream the inputStream to read from
  * @param pdbFileNameBase the base name of the PDB file (the part before the extension).
  * @return the Protein built from the information in this PDB
  * @throws IllegalStateException if this method is attempted to be called more than once.
  */
  public Protein readInPDBFile(InputStream inputStream, String pdbFileNameBase){
    if(this.proteinReadIn){
      throw new IllegalStateException("A PDBFileIO Object can only read in a single PDB. " +
          " If you wish to read in a second PDB file, you must instantiate another PDBFileIO.");
    }
    this.proteinReadIn = true;
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
    return constructProtein(pdbFileNameBase);
  }

  /**
  * Write prot out to file.
  * @param outputStream the OutputStreamWriter to write out to
  * @param prot the Protein to write out to file.
  * @throws IOException if an I/O error occurs
  */
  public void writeToPDB(OutputStreamWriter outputStream, Protein prot) throws IOException{
    String chainID = null;
    int resID = -1;
    String resName = null;
    int serialNum = -1;
    for(PolypeptideChain chain : prot){
      chainID = chain.getChainID();
      for(Residue res : chain){
        resID = res.getResidueID();
        resName = res.getThreeLetterName();
        for(Atom atom : res){
          serialNum = atom.getSerialNumber();
          String altLoc = " ";
          String iCode = " ";
          Vector3D loc = atom.getCoordinates();
          AtomRecord atomRecord = new AtomRecord(serialNum, atom.getName(), altLoc,
              resName, chainID, resID, iCode, loc.getX().doubleValue(), loc.getY().doubleValue(),
              loc.getZ().doubleValue(), atom.getOccupancy(), atom.getTempFactor(),
              atom.getElement(), atom.getCharge());
          outputStream.write(atomRecord.toString());
          outputStream.flush();
        }
      }
      String ter = String.format("%-6s%5d      %s %s%4s\n",
          "TER", serialNum+1, resName, chainID, resID);
      outputStream.write(ter);
      outputStream.flush();
    }
    String end = "END\n";
    outputStream.write(end);
    outputStream.flush();
  }

  private static HashMap<String, ArrayList<AtomRecord>> readInDefaultResidues(){
    InputStream stream = PDBFileIO.class.getResourceAsStream("residues.pdb");
    Scanner in = new Scanner(stream);
    ArrayList<AtomRecord> allResAtomRecords = new ArrayList<AtomRecord>();
    while(in.hasNextLine()){
      String line = in.nextLine();
      if(line.length() < 80){
        line = padRight(line);
      }
      allResAtomRecords.add(parseAtomLine(line));
    }
    in.close();
    return buildDefaultResiduesLists(allResAtomRecords);
  }

  /**
  * Build HashMap of the lists of AtomRecords for each of the default residues.
  *
  * The key is the residue name, the value is an ArrayList holding all the atomRecords for
  * that residue.
  */
  private static HashMap<String, ArrayList<AtomRecord>> buildDefaultResiduesLists(
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
  * Return all the atom records for a defaultResidue with all atomSerialNumber resID fields
  * set to -1.
  * @param residueThreeLetterID a three letter ID of a residue in all caps.
  * @return an ArrayList containing all the atomRecords for that residue. includes hydrogens.
  * @throws IllegalArgumentException if the three letter ID is not a valid ID.
  */
  private static ArrayList<AtomRecord> getDefaultResidueRecords(String residueThreeLetterID){
    if(defaultResiduesRecords.containsKey(residueThreeLetterID.toUpperCase())){
      return defaultResiduesRecords.get(residueThreeLetterID.toUpperCase());
    }
    throw new IllegalArgumentException(residueThreeLetterID + " not a valid residue name.");
  }

  /**
  * @param resThreeLetterName the three letter name in all caps of the residue to get the atoms of.
  * @return the default atoms for the residue corresponding to the three letter name.
  * This residue is at arbitrary coordinates. It will have to be translates and rotated
  * to a desired location
  */
  public static Collection<Atom> getDefaultResidueAtoms(String resThreeLetterName){
        return constructAtoms(getDefaultResidueRecords(resThreeLetterName));
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

  public Collection<SSBondRecord> getSSBondRecords(){
    return new ArrayList<SSBondRecord>(this.ssBonds);
  }

  private static String padRight(String line){
    return String.format("%-" + 80 + "s", line);
  }

  private static AtomRecord parseAtomLine(String line){
    int serial = Integer.parseInt(line.substring(6,11).trim());
    String atomName = line.substring(12,16).trim();
    String altLoc = line.substring(16,17).trim();
    String resName = line.substring(17,20).trim();
    String chainID = line.substring(21,22).trim();
    int resSeq = Integer.parseInt(line.substring(22,26).trim());
    String iCode = line.substring(26,27).trim(); // code for insertion of residues
    double x = Double.parseDouble(line.substring(30,38).trim());
    double y = Double.parseDouble(line.substring(38,46).trim());
    double z = Double.parseDouble(line.substring(46,54).trim());

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
    // PDB format specifies charge in form of 1+ or 2+ or 1-, but it could also
    // be +1, +2, or -1. It could also omit the sign.
    String chargeString = line.substring(78,80).trim();
    double charge = parseChargeString(chargeString);

    return new AtomRecord(serial, atomName, altLoc, resName, chainID, resSeq, iCode,
                          x, y, z, occupancy, tempFactor, element, charge);
  }

  private static double parseChargeString(String chargeString){
    // By PDB Specification, the chargeString should be 2+, 1-, etc,
    // but they could also be +2, -1, etc, or even simply 1, 2.
    double charge = 0;
    if(!chargeString.equals("")){
      // The charge string can be 1 or 2 characters long.
      // If the chargeString is only 1 character long and that character
      // is a digit, save it as the charge.
      if(chargeString.length() == 1){
        if(chargeString.matches("\\d")){
          charge = Double.parseDouble(chargeString);
        }
      } else if(chargeString.length() == 2){
        // The charge string must be 2 characters long
        String firstChar = chargeString.substring(0,1);
        String secondChar = chargeString.substring(1,2);
        // If the first character is a digit, save it as the charge, and if the secondChar
        // is "-", make the charge negative.
        // Otherwise, if the second char is a digit, use it for the charge and the firstChar
        // to determine the sign.
        if(firstChar.matches("\\d")){
          charge = Double.parseDouble(firstChar);
          if(secondChar.equals("-")){
            charge *= -1;
          }
        } else if(secondChar.matches("\\d")){
          charge = Double.parseDouble(secondChar);
          if(firstChar.equals("-")){
            charge *= -1;
          }
        }
      }
    }
    return charge;
  }

  private static SSBondRecord parseSSBondLine(String line){
    String chainID1 = line.substring(15,16).trim();
    int resID1 = Integer.parseInt(line.substring(17,21).trim());
    String chainID2 = line.substring(29,30).trim();
    int resID2 = Integer.parseInt(line.substring(31,35).trim());
    return new SSBondRecord(chainID1, resID1, chainID2, resID2);
  }

  //*************** Protein Construction Methods *******************//
  private Protein constructProtein(String pdbName){
    Protein protein = new Protein(pdbName, this);
    Collection<String> chainIDs = getListOfChainIDs();
    for(String chainID : chainIDs){
      PolypeptideChain chain = new PolypeptideChain(chainID);
      // get a list of all the atomRecords that have chainID
      Collection<AtomRecord> chainAtomRecords = getAtomRecordsInChain(chainID);

      HashMap<Integer, ArrayList<AtomRecord>> residueRecordsLists =
          getResidueRecordsLists(chainAtomRecords);

      // for every residue's list of AtomRecords, build a list of the atoms in that residue from the
      // AtomRecords. Use that list to construct and add that new residue with those atoms to the
      // chain.
      for(ArrayList<AtomRecord> residueAtomRecords : residueRecordsLists.values()){
        String resName = residueAtomRecords.get(0).getResName();
        int resSeq = residueAtomRecords.get(0).getResSeq();

        Collection<Atom> residueAtoms = constructAtoms(residueAtomRecords);
        // only build and add the residue if it contains the backbone atoms N CA C
        if(resContainsBackboneAtoms(residueAtoms)){
          try{
            Residue res = new Residue(resName, resSeq, residueAtoms);
            if(containsCarboxylOxygen(residueAtoms)){
              res.setAsCarboxylTerminus();
            }
            chain.addResidue(res);
          } catch (IllegalArgumentException e){
            throw new IllegalStateException(
                String.format("Problem when trying to construct Residue %s-%d for Protein %s",
                              resName, resSeq, pdbName), e);
          }
        }
      }
      protein.addChain(chain);
    }

    addDisulfideBonds(protein);
    return protein;
  }

  /**
  * @param atoms a Collection of atoms that make up a residue
  * @return true if atoms contains atoms with the names N, CA, and C. All three must be present
  */
  private boolean resContainsBackboneAtoms(Collection<Atom> atoms){
    boolean n = false;
    boolean ca = false;
    boolean c = false;
    for(Atom atom : atoms){
      if(atom.getName().equals("N")){
        n = true;
      }
      if(atom.getName().equals("CA")){
        ca = true;
      }
      if(atom.getName().equals("C")){
        c = true;
      }
    }
    return (n && ca && c);
  }

  private boolean containsCarboxylOxygen(Collection<Atom> atoms){
    for(Atom a: atoms){
      if(a.getName().equals("OXT")){
        return true;
      }
    }
    return false;
  }

  private void addDisulfideBonds(Protein protein){
    for(SSBondRecord ssBondRec: getSSBondRecords()){
      String chainID1 = ssBondRec.getChainID1();
      String chainID2 = ssBondRec.getChainID2();
      int resID1 = ssBondRec.getResID1();
      int resID2 = ssBondRec.getResID2();
      if(protein.contains(chainID1) && protein.getChain(chainID1).contains(resID2)){
        if(protein.contains(chainID2) && protein.getChain(chainID2).contains(resID2)){
          if(!chainID1.equals(chainID2) || resID1 != resID2){
            protein.addDisulfideBond(chainID1, resID1, chainID2, resID2);
          }
        }
      }
    }
  }


  /**
  * From a Collection of Atom Records, build and return a collection of those Atoms.
  */
  private static Collection<Atom> constructAtoms(Collection<AtomRecord> atomRecords){
    ArrayList<Atom> atoms = new ArrayList<Atom>();
    for(AtomRecord rec : atomRecords){
      atoms.add(constructAtom(rec));
    }
    return atoms;
  }

  private static Atom constructAtom(AtomRecord rec){
    return new Atom(rec.getName(), rec.getSerial(), rec.getOccupancy(),
        rec.getTempFactor(), rec.getCharge(),
        rec.getX(), rec.getY(), rec.getZ());
  }

  private static class AtomRecord{
    private int serial;
    private String atomName;
    private String altLoc;
    private String resName;
    private String chainID;
    private int resSeq;
    private String iCode;
    private double x;
    private double y;
    private double z;
    private double occupancy;
    private double tempFactor;
    private String element;
    private double charge;

    public AtomRecord(int serial, String atomName, String altLoc,
        String resName, String chainID, int resSeq,
        String iCode, double x, double y, double z,
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
    public String getName(){ return this.atomName; }
    public String getAltLoc(){ return this.altLoc; }
    public String getResName(){ return this.resName; }
    public String getChainID(){ return this.chainID; }
    public int getResSeq(){ return this.resSeq; }
    public String getICode(){ return this.iCode; }
    public double getX(){ return this.x; }
    public double getY(){ return this.y; }
    public double getZ(){ return this.z; }
    public double getOccupancy(){ return this.occupancy; }
    public double getTempFactor(){ return this.tempFactor; }
    public String getElement(){ return this.element; }
    public double getCharge(){ return this.charge; }
    public String getChargeString(){
      int charge = (int)getCharge();
      return ""+charge;
    }

    @Override
    public String toString(){
      String record = String.format("%-6s%5d %-4s%s%3s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
          "ATOM", getSerial(), getName(), getAltLoc(), getResName(), getChainID(),
          getResSeq(), getICode(), getX(), getY(), getZ(), getOccupancy(), getTempFactor(),
          getElement(), getChargeString());
      return record;
    }
  }

  private static class SSBondRecord{
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

}
