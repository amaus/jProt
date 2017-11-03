package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.transformations.Transformable;
import com.aaronpmaus.jMath.transformations.Transformation;
import com.aaronpmaus.jProt.io.PDBFileIO;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Collection;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.util.Iterator;

import java.io.File;
import java.io.InputStream;

/**
* A Residue is one of the 20 standard Amino acids.
* <p>
* By default, residues are not terminal. If a residue is to be the c-terminus, the method
* setAsCarboxylTerminus() must be called BEFORE the Residue is added to the chain, otherwise the
* carboxyl Oxygen will be missing from the Chain.
*
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class Residue implements Iterable<Atom>, Transformable{
  private final String name;
  private final String threeLetterName;
  private final String oneLetterName;
  private final int residueID;
  private boolean residueComplete = true;
  private static int maxResidueID = 0;
  private boolean hydrogensEnabled = false;

  // valid keys are the atomNames: CA, CB, CD, CD1, CD2, CE, C, O, N, etc...
  private HashMap<String, Atom> heavyAtoms;
  private HashMap<String, Atom> hydrogens;
  private Atom carboxylOxygen;
  private HashSet<Bond> bonds;
  // bonds_hydrogens: unusual var name format to avoid confusion with hydrogen bonds
  private HashSet<Bond> bondsToHydrogens;

  private static String[][] resNames = {
    {"ALA","A","Alanine"},        {"GLY","G","Glycine"},
    {"ILE","I","Isoleucine"},     {"LEU","L","Leucine"},
    {"PRO","P","Proline"},        {"VAL","V","Valine"},
    {"PHE","F","Phenylalanine"},  {"TRP","W","Tryptophan"},
    {"TYR","Y","Tyrosine"},       {"ASP","D","Aspartic Acid"},
    {"GLU","E","Glutamic Acid"},  {"ARG","R","Arginine"},
    {"HIS","H","Histidine"},      {"LYS","K","Lysine"},
    {"SER","S","Serine"},         {"THR","T","Threonine"},
    {"CYS","C","Cysteine"},       {"MET","M","Methionine"},
    {"ASN","N","Asparagine"},     {"GLN","Q","Glutamine"},
    {"HID","H","Histidine"},      {"HIE","H","Histidine"},
    {"HIP","H","Histidine"},      {"CYH","C","Cysteine"},
    {"ASH","D","Aspartic Acid"},  {"LYN","K","Lysine"},
    {"GLH","E","Glutamic Acid"},
  };

  protected void addAtom(Atom atom){
    this.heavyAtoms.put(atom.getAtomName(), atom);
  }
  /**
  * A constructor for a residue that builds a default residue of this type
  * @param threeLetterID the three letter name of the residue to builds
  * @param residueID the numeric ID of the residue to build
  */
  public Residue(String threeLetterID, int residueID){
    this(threeLetterID, residueID, PDBFileIO.getDefaultResidueAtoms(threeLetterID));
  }

  /**
  * Construct a default residue
  * @param oneLetterName the one letter name of the residue
  * @param residueID the numeric ID of the residue to build
  */
  public Residue(char oneLetterName, int residueID){
    this(lookUpThreeLetterName(String.format("%c",oneLetterName)),residueID);
  }

  /**
  * A constructor for an amino acid.
  * @param threeLetterName the one letter ID of the amino acid to build.
  * @param residueID the numeric residue ID of the residue being built.
  * @param atoms the atoms in this residue
  */
  public Residue(String threeLetterName, int residueID, Collection<Atom> atoms){
    this.hydrogensEnabled = true;
    this.threeLetterName = threeLetterName;
    this.oneLetterName = Residue.lookUpOneLetterName(this.threeLetterName);
    this.name = Residue.lookUpFullName(this.threeLetterName);
    if(residueID == -1){
      this.residueID = Residue.maxResidueID + 1;
    } else {
      this.residueID = residueID;
    }
    Residue.maxResidueID = Math.max(this.residueID, Residue.maxResidueID);
    this.heavyAtoms = new HashMap<String, Atom>();
    this.bonds = new HashSet<Bond>();
    this.hydrogens = new HashMap<String, Atom>();
    this.bondsToHydrogens = new HashSet<Bond>();
    initializeAminoAcid(this.threeLetterName+".dat", atoms);
    this.hydrogensEnabled = false;
  }

  // For Initializing all Residues:
  // use the atoms passed in. add all to hashMap and create all bonds.
  // What to do about missing atoms? Flag this residue as incomplete. Then have the protein
  // constructor resolve all incomplete residues after they've all been added.
  private void initializeAminoAcid(String dataFileName, Collection<Atom> atoms){
    this.residueComplete = true;
    // if atoms is null, build amino acid from default values
    if(atoms == null){

    } else {
      if(atoms.isEmpty()){
        throw new IllegalArgumentException("Collection of Atoms must not be empty");
      }
      // use the atoms passed in
      for(Atom a: atoms){
        //System.out.printf("Adding %s to residue atoms\n",a.getAtomName());
        if(a.getElement().equals("H")){
          this.hydrogens.put(a.getAtomName(), a);
        } else if(a.getAtomName().equals("OXT")){
          this.carboxylOxygen = a;
        } else {
          this.heavyAtoms.put(a.getAtomName(), a);
        }
      }
      // read in residue data file. use it to add all main bonds.
      // Then add all hydrogens
      InputStream stream = Residue.class.getResourceAsStream(dataFileName);
      Scanner in = new Scanner(stream);
      boolean readingInMainBonds = false;
      boolean readingInHydrogens = false;
      while(in.hasNext()){
        String line = in.nextLine().trim();
        //System.out.println(line);
        String firstChar = line.substring(0,1);
        if(line.equals("!main bonds")){
          readingInMainBonds = true;
        }
        if(line.equals("!hydrogens")){
          readingInMainBonds = false;
          readingInHydrogens = true;
        }
        if(line.equals("!defined dihedral angles")){
          readingInMainBonds = false;
          readingInHydrogens = false;
        }
        if(line.equals("!atom types")){
          readingInMainBonds = false;
          readingInHydrogens = false;
        }
        if((!firstChar.equals("!")) && readingInMainBonds){
          String[] tokens = line.trim().split("\\s+");
          String atomOne = tokens[0].trim();
          String atomTwo = tokens[1].trim();
          //System.out.printf("[%s-%s]\n",atomOne, atomTwo);
          // TODO specify single or double bond depending on
          // atoms and residue
          if(atomTwo.equals("OXT")){
            // do nothing. The bond will be added if this residue is setAsCarboxylTerminus()
          } else {
            residueComplete = addBond(atomOne, atomTwo) & residueComplete;
          }
        }
        if((!firstChar.equals("!")) && readingInHydrogens){
          String[] tokens = line.trim().split("\\s+");
          String atomOne = tokens[0];
          String atomTwo = tokens[1];
          // TODO specify single or double bond depending on
          // atoms and residue

          addBondToHydrogen(atomOne, atomTwo);
        }
      }
    }
  }

  /**
  * gets an atom by its name. C, CA, CB, CD, etc...
  * @param atomName the name of the atom, C, CA, CB, CD, etc...
  * @return the Atom in this residue with that name
  */
  public Atom getAtom(String atomName){
    Atom atom = this.heavyAtoms.get(atomName);
    if(this.hydrogensEnabled() && atom == null){
      atom = this.hydrogens.get(atomName);
    }
    if(atom == null){
      throw new NoSuchElementException("No atom of name " + atomName
          + " in residue " + getResidueID() +": " +getThreeLetterName());
    }
    return atom;
  }

  /**
  * @param atomName the name of an Atom to check for. These names are formatted as the name
  * field in the Atom record of the PDB File Format.
  * @return true if this Residue contains an Atom with this name
  */
  public final boolean contains(String atomName){
    if(this.heavyAtoms.containsKey(atomName)){
      return true;
    }
    if(this.hydrogensEnabled() && this.hydrogens.containsKey(atomName)){
      return true;
    }
    return false;
  }

  /**
  * Return the numeric residue ID of this residue
  *
  * @return the residue ID of this residue
  */
  public int getResidueID(){
    return this.residueID;
  }

  /**
  * Return the single letter residue name of this residue
  *
  * @return the single letter residue name of this residue
  */
  public String getOneLetterName(){
    return this.oneLetterName;
  }

  /**
  * Return the three letter residue name of this residue
  *
  * @return the three letter residue name of this residue
  */
  public String getThreeLetterName(){
    return this.threeLetterName;
  }

  /**
  * Return the full name of this residue
  *
  * @return the full name of this residue
  */
  public String getName(){
    return this.name;
  }

  /**
  * Get a Collection of the covalent bonds in this residue
  *
  * @return the bonds in this residue
  */
  public Collection<Bond> getBonds(){
    ArrayList<Bond> bonds = new ArrayList<Bond>(this.bonds);
    if(this.hydrogensEnabled()){
      bonds.addAll(getBondsToHydrogens());
    }
    return bonds;
  }

  /**
  * Return a Collection of the Bonds to hydrogens in this residue. It will return these bonds
  * whether hydrogens are enabled or not.
  *
  * @return a {@code Collection<Bond>} containing all bonds with a hydrogen
  */
  public Collection<Bond> getBondsToHydrogens(){
    return new ArrayList<Bond>(this.bondsToHydrogens);
  }

  /**
  * @return a {@code Collection<Atom>} containing all the atoms in this residue. If hydrogens are
  * enabled, it will include the hydrogens
  */
  public Collection<Atom> getAtoms(){
    ArrayList<Atom> atoms = new ArrayList<Atom>(this.heavyAtoms.values());
    if(this.hydrogensEnabled()){
      atoms.addAll(this.getHydrogens());
    }
    return atoms;
  }

  /**
  * Return all the heavy atoms in this residue
  * @return a {@code Collection<Atom>} containing the heavy atoms in this residue
  */
  protected Collection<Atom> getHeavyAtoms(){
    return this.heavyAtoms.values();
  }

  /**
  * Return all the hydrogens in this residue whether they are enabled or not.
  * @return a {@code Collection<Atom>} containing the hydrogens in this residue
  */
  protected Collection<Atom> getHydrogens(){
    return this.hydrogens.values();
  }

  /**
  * Return the number of Atoms in this Residue. This includes all heavy atoms and hydrogens.
  * @return the number of Atoms in this Residue
  */
  public int getNumAtoms(){
    int numAtoms = getNumHeavyAtoms();
    if(this.hydrogensEnabled()){
      numAtoms += getNumHydrogens();
    }
    return numAtoms;
  }

  /**
  * Return the number of Heavy Atoms in this Residue.
  * @return the number of Heavy Atoms in this Residue
  */
  public int getNumHeavyAtoms(){
    return this.heavyAtoms.size();
  }

  /**
  * Return the number of Hydrogens in this Residue.
  * @return the number of Hydrogens in this Residue
  */
  public int getNumHydrogens(){
    return this.hydrogens.size();
  }

  /**
  * Check if this residue is missing any of it's heavy atoms - non hydrogen atoms.
  * This is useful after a protein has been read in from PDB to know which residues have
  * missing atoms.
  * @return true if this residue is missing any of its heavy atoms.
  */
  public boolean isMissingAtoms(){
    return !(this.residueComplete);
  }

  /**
  * Set this Residue as the Carboxyl Terminus. This adds the OXT atom to this residue.
  */
  public void setAsCarboxylTerminus(){
    if(carboxylOxygen != null){
      this.heavyAtoms.put(carboxylOxygen.getAtomName(), carboxylOxygen);
      if(contains("C")){
        Atom carbon = getAtom("C");
        this.bonds.add(new Bond(carbon, carboxylOxygen, 1));
      }
    }
  }

  private boolean addBond(String atomNameOne, String atomNameTwo){
    return addBond(atomNameOne, atomNameTwo,1);
  }

  private boolean addBond(String atomNameOne, String atomNameTwo, int bondStrength){
    if(contains(atomNameOne) && contains(atomNameTwo)){
      this.bonds.add(new Bond(getAtom(atomNameOne), getAtom(atomNameTwo), bondStrength));
      return true;
    }
    return false;
  }

  private void addBondToHydrogen(String heavyAtom, String hydrogen){
    if(contains(heavyAtom) && contains(hydrogen)){
      this.bondsToHydrogens.add(new Bond(getAtom(heavyAtom), getAtom(hydrogen), 1));
    }
  }

  /**
  * Static lookup method to get a residue one letter name from
  * a three letter name.
  * @param threeLetterName the three letter name of the residue
  * @return the one letter name of the residue
  */
  private static String lookUpOneLetterName(String threeLetterName){
    threeLetterName = threeLetterName.toUpperCase();
    for(String[] triplet: resNames){
      if(triplet[0].equals(threeLetterName)){
        return triplet[1];
      }
    }
    throw new IllegalArgumentException(threeLetterName + ": not a valid residue name.");
  }

  /**
  * Static lookup method to get a residue three letter name from
  * a one letter name.
  * @param oneLetterName the one letter name of the residue
  * @return the three letter name of the residue
  */
  private static String lookUpThreeLetterName(String oneLetterName){
    oneLetterName = oneLetterName.toUpperCase();
    for(String[] triplet: resNames){
      if(triplet[1].equals(oneLetterName)){
        return triplet[0];
      }
    }
    throw new IllegalArgumentException(oneLetterName + ": not a valid residue name.");
  }

  /**
  * Static lookup method to get a residue full name from
  * a three letter name.
  * @param threeLetterName the three letter name of the residue
  * @return the full name of the residue
  */
  private static String lookUpFullName(String threeLetterName){
    threeLetterName = threeLetterName.toUpperCase();
    for(String[] triplet: resNames){
      if(triplet[0].equals(threeLetterName)){
        return triplet[2];
      }
    }
    throw new IllegalArgumentException(threeLetterName + ": not a valid residue name.");
  }

  /**
  * Enable Hydrogens for this Residue. The methods {@code getAtoms()} and  {@code getBonds()} will
  * now return Collections which include Hydrogens and Bonds containing hydrogens respectively.
  */
  public void enableHydrogens(){
    this.hydrogensEnabled = true;
  }

  public void disableHydrogens(){
    this.hydrogensEnabled = false;
  }

  public boolean hydrogensEnabled(){
    return this.hydrogensEnabled;
  }

  @Override
  public void applyTransformation(Transformation t){
    for(Atom atom : getHeavyAtoms()){
      atom.applyTransformation(t);
    }
    for(Atom atom : getHydrogens()){
      atom.applyTransformation(t);
    }
  }

  @Override
  public String toString(){
    String str = "";
    str += getName() + "\n";
    str += getOneLetterName() + "\n";
    str += getResidueID() + "\n";
    str += "Atoms:\n";
    for(Atom atom : this){
      str += atom.getAtomName()+"\n";
    }
    return str;
  }

  @Override
  public Iterator<Atom> iterator(){
    return getAtoms().iterator();
  }
}
