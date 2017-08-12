package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Collection;
import java.util.Scanner;
import java.util.ArrayList;

import java.io.FileNotFoundException;
import java.io.File;
import java.io.InputStream;
import java.net.URL;

/**
* This class represents in general an AminoAcid.
*
* An AminoAcid consists of a collection of Atoms bonded together
* is a particular way.
*
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class Residue{
  private final String name;
  private final String threeLetterName;
  private final String oneLetterName;
  private final int residueID;
  private boolean residueComplete = true;

  // valid keys are the atomNames: CA, CB, CD, CD1, CD2, CE, C, O, N, etc...
  private HashMap<String, Atom> atoms;
  private HashSet<Bond> bonds;

  private static String[][] resNameLookupTable = {
    {"ALA","A","Alanine"},        {"GLY","G","Glycine"},
    {"ILE","I","Isoleucine"},     {"LEU","L","Leucine"},
    {"PRO","P","Proline"},        {"VAL","V","Valine"},
    {"PHE","F","Phenylalanine"},  {"TRP","W","Tryptophan"},
    {"TYR","Y","Tyrosine"},       {"ASP","D","Aspartic Acid"},
    {"GLU","E","Glutamic Acid"},  {"ARG","R","Arginine"},
    {"HIS","H","Histidine"},      {"LYS","K","Lysine"},
    {"SER","S","Serine"},         {"THR","T","Threonine"},
    {"CYS","C","Cysteine"},       {"MET","M","Methionine"},
    {"ASN","N","Asparagine"},     {"GLN","Q","Glutamine"}
  };

  /**
  * A constructor for a residue that builds a default residue of this type
  * @param threeLetterID the three letter name of the residue to builds
  * @param residueID the numeric ID of the residue to build
  */
  public Residue(String threeLetterID, int residueID){
    this(threeLetterID,residueID,null);
  }

  /**
  * A constructor for an amino acid.
  * @param threeLetterName the one letter ID of the amino acid to build.
  * @param residueID the numeric residue ID of the residue being built.
  * @param atoms the atoms in this residue
  */
  public Residue(String threeLetterName, int residueID, Collection<Atom> atoms){
    this.threeLetterName = threeLetterName;
    this.oneLetterName = Residue.lookUpOneLetterName(this.threeLetterName);
    this.name = Residue.lookUpFullName(this.threeLetterName);
    this.residueID = residueID;
    this.atoms = new HashMap<String, Atom>();
    this.bonds = new HashSet<Bond>();
    initializeAminoAcid(this.threeLetterName+".dat", atoms);
  }

  /**
  * gets an atom by its name. C, CA, CB, CD, etc...
  * @param atomName the name of the atom, C, CA, CB, CD, etc...
  * @return the Atom in this residue with that name
  */
  public Atom getAtom(String atomName){
    return atoms.get(atomName);
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
  * @return a {@code HashSet<Bond>} of the bonds in this residue
  */
  public Collection<Bond> getBonds(){
    return new ArrayList<Bond>(this.bonds);
  }

  // For Initializing all Residues:
  // use the atoms passed in. add all to hashMap and create all bonds.
  // What to do if any atoms are missing? Easy solution: treat it
  // as a missing residue. Throw an exception if missing any atoms.
  // then the client will have to ignore this residue.
  // Hard Solution: fill in the missing atoms. remember, no clashes.
  // Umm, flag this residue as incomplete. Then have the protein
  // constructor resolve all incomplete residues after they've all
  // been added?
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
        this.atoms.put(a.getAtomName(), a);
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
        if((!firstChar.equals("!")) && readingInMainBonds){
          String[] tokens = line.trim().split(" ");
          String atomOne = tokens[0];
          String atomTwo = tokens[1];
          // TODO specify single or double bond depending on
          // atoms and residue
          if(atomTwo.equals("OXT")){
            addBond(atomOne, atomTwo);
          } else {
            residueComplete = addBond(atomOne, atomTwo)
            & residueComplete;
          }
        }
        if((!firstChar.equals("!")) && readingInHydrogens){
          String[] tokens = line.trim().split(" ");
          String atomOne = tokens[0];
          String atomTwo = tokens[1];
          // TODO specify single or double bond depending on
          // atoms and residue
          addBond(atomOne, atomTwo);
        }
      }
    }
  }

  private boolean addBond(String atomNameOne, String atomNameTwo){
    return addBond(atomNameOne, atomNameTwo,1);
  }

  private boolean addBond(String atomNameOne, String atomNameTwo, int bondStrength){
    if(this.atoms.containsKey(atomNameOne) && this.atoms.containsKey(atomNameTwo)){
      this.bonds.add(new Bond(getAtom(atomNameOne), getAtom(atomNameTwo), bondStrength));
      return true;
    }
    return false;
  }

  /**
  * Static lookup method to get a residue one letter name from
  * a three letter name.
  * @param threeLetterName the three letter name of the residue
  * @return the one letter name of the residue
  */
  public static String lookUpOneLetterName(String threeLetterName){
    threeLetterName = threeLetterName.toUpperCase();
    for(String[] pair: resNameLookupTable){
      if(pair[0].equals(threeLetterName)){
        return pair[1];
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
  public static String lookUpThreeLetterName(String oneLetterName){
    oneLetterName = oneLetterName.toUpperCase();
    for(String[] pair: resNameLookupTable){
      if(pair[1].equals(oneLetterName)){
        return pair[0];
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
  public static String lookUpFullName(String threeLetterName){
    threeLetterName = threeLetterName.toUpperCase();
    for(String[] pair: resNameLookupTable){
      if(pair[0].equals(threeLetterName)){
        return pair[2];
      }
    }
    throw new IllegalArgumentException(threeLetterName + ": not a valid residue name.");
  }

  /**
  * Returns the number of Atoms in this Residue
  * @return the number of Atoms in this Residue
  */
  public int getNumAtoms(){
    return this.atoms.size();
  }

  public String toString(){
    String str = "";
    str += getName() + "\n";
    str += getOneLetterName() + "\n";
    str += getResidueID() + "\n";
    return str;
  }
}
