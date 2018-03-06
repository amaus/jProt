package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.linearAlgebra.Vector3D;
import com.aaronpmaus.jMath.transformations.Transformable;
import com.aaronpmaus.jMath.transformations.Transformation;
import com.aaronpmaus.jProt.io.PDBFileIO;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Iterator;
import java.util.Comparator;

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
* @version 0.7.0
* @since 0.6.0
*/
public class Residue implements Iterable<Atom>, Transformable{
  private final String name;
  private final String threeLetterName;
  private final String oneLetterName;
  private final int residueID;
  private boolean residueComplete = true;
  private static int maxResidueID = 0;
  private boolean isCarboxylTerminus = false;

  // valid keys are the atomNames: CA, CB, CD, CD1, CD2, CE, C, O, N, etc...
  private HashMap<String, Atom> heavyAtoms;
  private HashMap<String, Atom> hydrogens;
  private HashMap<String, ArrayList<String>> heavyAtomHydrogens;
  private Atom carboxylOxygen;
  private ArrayList<Bond> bonds;
  private ArrayList<Bond> rotatableBonds;
  private ArrayList<AngleTriple> angleTriplets;
  // bonds_hydrogens: unusual var name format to avoid confusion with hydrogen bonds
  private ArrayList<Bond> bondsToHydrogens;
  // A list of atom types, every consecutive 4 define a dihedral angle
  // Eg.
  // N CA CB CG CD NE CZ
  // There are 4 dihedral angles:
  // N-CA-CB-CG, CA-CB-CG-CD, CB-CG-CD-NE, and CG-CD-NE-CZ
  private ArrayList<String> definedDihedralAngles;
  // A list of all the atom types in the residue, ordered by the level they occupy.
  // For example, all atoms at the delta level are before the atoms at the epsilon level.
  // excludes backbone atoms.
  private ArrayList<String> atomsInOrder;
  // a directed graph of all the atoms in this residue. contains edges as specified by the .dat
  // files in the resources directory for this package.
  private Graph<Atom> atomsGraph;

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
    {"GLH","E","Glutamic Acid"},  {"UNK","X","Unknown"}
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
  public Residue(String threeLetterName, int residueID, Collection<Atom> atoms) {
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
    this.bonds = new ArrayList<Bond>();
    this.hydrogens = new HashMap<String, Atom>();
    this.bondsToHydrogens = new ArrayList<Bond>();
    this.definedDihedralAngles = new ArrayList<String>();
    this.heavyAtomHydrogens = new HashMap<String, ArrayList<String>>();
    this.rotatableBonds = new ArrayList<Bond>();
    this.angleTriplets = new ArrayList<AngleTriple>();
    this.atomsGraph = new Graph<Atom>();
    this.atomsInOrder = new ArrayList<String>();
    this.atomsInOrder.add("N");
    this.atomsInOrder.add("CA");
    initializeAminoAcid(this.threeLetterName+".dat", atoms);
    this.atomsInOrder.add("C");
    this.atomsInOrder.add("O");
    this.atomsInOrder.add("OXT");
    Collections.sort(this.bonds, new ResidueBondComparator());
    Collections.sort(this.bondsToHydrogens, new ResidueBondComparator());
  }

  // For Initializing all Residues:
  // use the atoms passed in. add all to hashMap and create all bonds.
  // What to do about missing atoms? Flag this residue as incomplete. Then have the protein
  // constructor resolve all incomplete residues after they've all been added.
  private void initializeAminoAcid(String dataFileName, Collection<Atom> atoms){
    this.residueComplete = true;
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

    while(in.hasNext()){
      String line = in.nextLine().trim();
      //System.out.println(line);
      String firstChar = line.substring(0,1);
      if(line.equals("!main bonds")){
        while(in.hasNext()) {
          line = in.nextLine().trim();
          if(!line.substring(0,1).equals("!")){
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
            if(!atomOne.equals("N") && !atomOne.equals("CA") && !atomOne.equals("C")
                && !atomOne.equals("O") && !atomOne.equals("OXT")){

              if(!this.atomsInOrder.contains(atomOne)){
                this.atomsInOrder.add(atomOne);
              }
            }
            if(!atomTwo.equals("N") && !atomTwo.equals("CA") && !atomTwo.equals("C")
                && !atomTwo.equals("O") && !atomTwo.equals("OXT")){

              if(!this.atomsInOrder.contains(atomTwo)){
                this.atomsInOrder.add(atomTwo);
              }
            }
          } else { // read next section header
            break;
          }
        }
      }
      if(line.equals("!hydrogens")){
        while(in.hasNext()) {
          line = in.nextLine().trim();
          if(!line.substring(0,1).equals("!")){
            String[] tokens = line.trim().split("\\s+");
            String heavyAtom = tokens[0];
            String hydrogen = tokens[1];
            addBondToHydrogen(heavyAtom, hydrogen);
            if(!this.atomsInOrder.contains(hydrogen)){
              int index = this.atomsInOrder.indexOf(heavyAtom);
              this.atomsInOrder.add(index+1, hydrogen);
            }
            if(heavyAtomHydrogens.containsKey(heavyAtom)){
              heavyAtomHydrogens.get(heavyAtom).add(hydrogen);
            } else {
              heavyAtomHydrogens.put(heavyAtom, new ArrayList<String>());
            }
          } else { // read next section header
            break;
          }
        }
      }
      if(line.equals("!defined dihedral angles")){
        line = in.nextLine().trim();
        if(!line.substring(0,1).equals("!")){
          String[] atomTypes = line.split("\\s+");
          for(String atomType : atomTypes){
            definedDihedralAngles.add(atomType);
          }
        }
      }
      if(line.equals("!defined angle triplets")){
        while(in.hasNext()) {
          line = in.nextLine().trim();
          if(!line.substring(0,1).equals("!")){
            if(!line.substring(0,1).equals("#")){
              String[] tokens = line.trim().split("\\s+");
              if(this.contains(tokens[0]) && this.contains(tokens[1]) && this.contains(tokens[2])){
                this.angleTriplets.add(new AngleTriple(getAtom(tokens[0]),
                                                      getAtom(tokens[1]),
                                                      getAtom(tokens[2])));
              }
            }
          } else { // read next section header
            break;
          }
        }
      }
      if(line.equals("!atom types")){
        // place holder for future use of this information
      }
    }
    for(Bond bond : this.getBonds()){
      if(isRotatableBond(bond.getAtomOne().getAtomName(), bond.getAtomTwo().getAtomName())){
        this.rotatableBonds.add(bond);
      }
    }
  }

  /**
  * gets an atom by its name. C, CA, CB, CD, etc...
  * @param atomName the name of the atom, C, CA, CB, CD, etc...
  * @return the Atom in this residue with that name
  * @throws NoSuchElementException if the atom is not present in this residue
  */
  public Atom getAtom(String atomName){
    Atom atom = this.heavyAtoms.get(atomName);
    if(atom == null){
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
    if(this.hydrogens.containsKey(atomName)){
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
    bonds.addAll(getBondsToHydrogens());
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
    atoms.addAll(this.getHydrogens());
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
    numAtoms += getNumHydrogens();
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
  * @return a List of the rotatable Bonds in this residue
  */
  public List<Bond> getRotatableBonds(){
    return this.rotatableBonds;
  }

  /**
  * Set this Residue as the Carboxyl Terminus. This adds the OXT atom to this residue.
  */
  public void setAsCarboxylTerminus(){
    if(carboxylOxygen != null){
      this.isCarboxylTerminus = true;
      this.heavyAtoms.put(carboxylOxygen.getAtomName(), carboxylOxygen);
      if(contains("C")){
        Atom carbon = getAtom("C");
        this.bonds.add(new Bond(carbon, carboxylOxygen, 1));
        this.atomsGraph.addEdge(carbon, carboxylOxygen);
      }
      if(this.contains("CA") && this.contains("C") && this.contains("OXT")){
        this.angleTriplets.add(new AngleTriple(getAtom("CA"), getAtom("C"), getAtom("OXT")));
      }
    }
  }

  /**
  * @return true if this residue is the carboxyl terminus
  */
  public boolean isCarboxylTerminus(){
    return this.isCarboxylTerminus;
  }

  private boolean addBond(String atomNameOne, String atomNameTwo){
    return addBond(atomNameOne, atomNameTwo,1);
  }

  private boolean addBond(String atomNameOne, String atomNameTwo, int bondStrength){
    if(contains(atomNameOne) && contains(atomNameTwo)){
      Atom atomOne = getAtom(atomNameOne);
      Atom atomTwo = getAtom(atomNameTwo);
      Bond bond = new Bond(atomOne, atomTwo, bondStrength);
      this.atomsGraph.addEdge(atomOne,atomTwo);
      this.bonds.add(bond);
      return true;
    }
    return false;
  }

  private void addBondToHydrogen(String heavyAtomName, String hydrogenName){
    if(contains(heavyAtomName) && contains(hydrogenName)){
      Atom heavyAtom = getAtom(heavyAtomName);
      Atom hydrogen = getAtom(hydrogenName);
      this.bondsToHydrogens.add(new Bond(heavyAtom, hydrogen, 1));
      this.atomsGraph.addEdge(heavyAtom, hydrogen);
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
  * Given one of the rotatable bonds, set its dihedral angle to degrees.
  * @param bond a rotatable bond in this residue
  * @param degrees the degrees to set the angle to
  * @throws IllegalArgumentException if the bond is not one of the rotatable bonds
  * @since 0.7.0
  */
  public void setDihedralAngle(Bond bond, double degrees){
    setDihedralAngle(bond.getAtomOne().getAtomName(), bond.getAtomTwo().getAtomName(), degrees);
  }

  /**
  * Given two atom names that specify a rotatable bond, set that bond's dihedral angle as specified.
  * <p>
  * The order the atoms are passed in does not matter.
  * @param atomOne one of the atoms in the bond
  * @param atomTwo the other atom in the bond
  * @param degrees the degrees to set the angle to
  * @throws IllegalArgumentException if the bond is not one of the rotatable bonds
  * @since 0.7.0
  */
  public void setDihedralAngle(String atomOne, String atomTwo, double degrees){
    double currentAngle = getDihedralAngle(atomOne, atomTwo);
    double delta = degrees - currentAngle;
    rotateAboutBond(atomOne, atomTwo, delta);
  }

  /**
  * Set the angle of these three atoms in this residue. The atoms must be consecutive.
  * @param atomOneName the first atom of the three
  * @param atomTwoName the second atom (the pivot) of the three
  * @param atomThreeName the third atom of the three
  * @param degrees the degrees to set the angle to.
  */
  public void setAngle(String atomOneName, String atomTwoName,
                        String atomThreeName, double degrees){
    // these must be 3 consecutive atoms. one pair must be a rotatable bond.
    // all atoms from and including the middle will be transformed.
    // build list of angle triplets at residue initialization. then check agains that list.
    // using graph and depthFirstSearch, N-CA-CB would result in the Transformation being applied
    // to C and O as well. Be careful.
    if(!contains(atomOneName) || !contains(atomTwoName) || !contains(atomThreeName)){
      throw new NoSuchElementException(
          String.format("%s, %s, and %s must all be present in the residue.\n%s",
                        atomOneName, atomTwoName, atomThreeName, toString()));
    }
    Atom atomOne = getAtom(atomOneName);
    Atom atomTwo = getAtom(atomTwoName);
    Atom atomThree = getAtom(atomThreeName);
    if(!angleTriplets.contains(new AngleTriple(atomOne, atomTwo, atomThree))){
      throw new IllegalArgumentException(String.format("%s-%s-%s not a valid triple.",
          atomOneName, atomTwoName, atomThreeName));
    }
    Transformation t = buildAngleTransformation(atomOne, atomTwo, atomThree, degrees);
    // special cases: N-CA-C, N-CA-CB,
    Collection<Atom> atomsToTransform = getAtomsAfterAtom(atomTwoName);
    if(atomOneName.equals("N") && atomTwoName.equals("CA") && atomThreeName.equals("CB")){
      // remove C, O, and OXT from atomsToTransform
      Iterator<Atom> it = atomsToTransform.iterator();
      while(it.hasNext()){
        Atom next = it.next();
        if(next.getAtomName().equals("C")
            || next.getAtomName().equals("O")
            || next.getAtomName().equals("OXT")){
          it.remove();
        }
      }
    }
    for(Atom atom : atomsToTransform){
      if(!atom.getAtomName().equals(atomOneName)) {
        atom.applyTransformation(t);
      }
    }
  }

  /**
  * @return a list of all valid angle triplets. The are the sets of three consecutive atoms whose
  * angle can be modified. They are returned as an ArrayList of Arrays of size 3.
  */
  public ArrayList<String[]> getAngleTriplets(){
    ArrayList<String[]> triplets = new ArrayList<String[]>();
    for(AngleTriple trip : angleTriplets){
      String[] names = new String[3];
      names[0] = trip.getAtomOne().getAtomName();
      names[1] = trip.getAtomTwo().getAtomName();
      names[2] = trip.getAtomThree().getAtomName();
      if(names[0].equals("N") && names[1].equals("CA") && names[2].equals("C")){
        continue;
      }
      if(names[0].equals("CA") && names[1].equals("C") && names[2].equals("O")){
        continue;
      }
      if(names[0].equals("CA") && names[1].equals("C") && names[2].equals("OXT")){
        continue;
      }
      triplets.add(names);
    }
    return triplets;
  }

  private Transformation buildAngleTransformation(Atom atomOne, Atom atomTwo,
                                                  Atom atomThree, double degrees){
    // AngleTriple sorts the atoms by their relative position from N
    AngleTriple trip = new AngleTriple(atomOne, atomTwo, atomThree);
    Vector3D atomOneCoor = trip.getAtomOne().getCoordinates();
    Vector3D atomTwoCoor = trip.getAtomTwo().getCoordinates();
    Vector3D atomThreeCoor = trip.getAtomThree().getCoordinates();
    double angle = getAngle(atomOne, atomTwo, atomThree);
    //System.out.printf("Current Angle: %.2f\n", angle);
    //System.out.printf("Goal: %.2f\n", degrees);
    double delta = degrees - angle;
    //System.out.printf("Delta: %.2f\n", delta);
    Vector3D normal = atomOneCoor.subtract(atomTwoCoor)
        .crossProduct(atomThreeCoor.subtract(atomTwoCoor));
    Transformation t = new Transformation();
    t.addRotationAboutAxis(atomTwoCoor, normal, delta);

    return t;
  }


  /**
  * Given one of the rotatable bonds, set its dihedral angle to degrees.
  * @param bond a rotatable bond in this residue
  * @param length the length (in Angstroms) to set the bond distance to
  * @throws IllegalArgumentException if the bond is not one of the rotatable bonds
  * @since 0.7.0
  */
  public void setRotatableBondLength(Bond bond, double length){
    setRotatableBondLength(bond.getAtomOne().getAtomName(), bond.getAtomTwo().getAtomName(), length);
  }

  /**
  * Given two atom names that specify a rotatable bond, set that bond's dihedral angle as specified.
  * <p>
  * The order the atoms are passed in does not matter.
  * @param atomOneName one of the atoms in the bond
  * @param atomTwoName the other atom in the bond
  * @param length the length (in Angstroms) to set the distance between these atoms to
  * @throws IllegalStateException if either atomOne or atomTwo are missing from this residue
  * @throws IllegalArgumentException if the bond is not one of the rotatable bonds
  * @since 0.7.0
  */
  public void setRotatableBondLength(String atomOneName, String atomTwoName, double length){
    if(!contains(atomOneName)){
      throw new IllegalStateException(String.format("Residue %d: Can't rotate about %s-%s bond,"
          + "missing atom %s", getResidueID(), atomOneName, atomTwoName, atomOneName));
    }
    if(!contains(atomTwoName)){
      throw new IllegalStateException(String.format("Residue %d: Can't rotate about %s-%s bond,"
          + "missing atom %s", getResidueID(), atomOneName, atomTwoName, atomTwoName));
    }
    if(!isRotatableBond(atomOneName, atomTwoName)){
      throw new IllegalArgumentException(String.format("%s-%s are not a rotable bond in %s",
          atomOneName, atomTwoName, this.getThreeLetterName()));
    }
    // ensure that the atoms are in the correct order, the second atom further from the backbone
    int indexOfAtomOne = this.definedDihedralAngles.indexOf(atomOneName);
    int indexOfAtomTwo = this.definedDihedralAngles.indexOf(atomTwoName);
    if(indexOfAtomOne > indexOfAtomTwo){
      String temp = atomOneName;
      atomOneName = atomTwoName;
      atomTwoName = temp;
    }
    // construct the Translation Transformation
    Transformation t = new Transformation();
    Vector3D a1 = getAtom(atomOneName).getCoordinates();
    Vector3D a2 = getAtom(atomTwoName).getCoordinates();
    Vector3D newAtomTwoPosition = a2.subtract(a1)
                                    .toUnitVector()
                                    .multiply(length)
                                    .add(a1);
    t.addTranslation(newAtomTwoPosition.subtract(a2));
    // get the atoms that need to be rotated, all those after this bond away from the backbone
    Collection<Atom> atomsToTransform = getAtomsAfterAtom(atomTwoName);
    for(Atom atom : atomsToTransform){
      atom.applyTransformation(t);
    }

  }

  /**
  * Rotate a residue about one of its rotatable bonds, one of its defined dihedral angles.
  * <p>
  * All atoms after that bond are transformed.
  * @param bond the bond to rotate about
  * @param degrees the number of degrees to rotate the residue
  * @throws IllegalStateException if either atomOne or atomTwo are missing from this residue
  * @throws IllegalArgumentException if the atomOne and atomTwo are not a rotatable bond in this
  * residue
  * @since 0.7.0
  */
  public void rotateAboutBond(Bond bond, double degrees){
    rotateAboutBond(bond.getAtomOne().getAtomName(), bond.getAtomTwo().getAtomName(), degrees);
  }

  /**
  * Rotate a residue about one of its rotatable bonds, one of its defined dihedral angles.
  * <p>
  * All atoms after that bond are transformed.
  * @param atomOneName the name of the first atom in the bond, the atom closer to the backbone
  * @param atomTwoName the name of the second atom in the bond, the atom further from the backbone
  * @param degrees the number of degrees to rotate the residue
  * @throws IllegalStateException if either atomOne or atomTwo are missing from this residue
  * @throws IllegalArgumentException if the atomOne and atomTwo are not a rotatable bond in this
  * residue
  * @since 0.7.0
  */
  public void rotateAboutBond(String atomOneName, String atomTwoName, double degrees){
    if(!contains(atomOneName)){
      throw new IllegalStateException(String.format("Residue %d: Can't rotate about %s-%s bond,"
          + "missing atom %s", getResidueID(), atomOneName, atomTwoName, atomOneName));
    }
    if(!contains(atomTwoName)){
      throw new IllegalStateException(String.format("Residue %d: Can't rotate about %s-%s bond,"
          + "missing atom %s", getResidueID(), atomOneName, atomTwoName, atomTwoName));
    }
    if(!isRotatableBond(atomOneName, atomTwoName)){
      throw new IllegalArgumentException(String.format("%s-%s are not a rotable bond in %s",
          atomOneName, atomTwoName, this.getThreeLetterName()));
    }
    // ensure that the atoms are in the correct order, the second atom further from the backbone
    int indexOfAtomOne = this.definedDihedralAngles.indexOf(atomOneName);
    int indexOfAtomTwo = this.definedDihedralAngles.indexOf(atomTwoName);
    if(indexOfAtomOne > indexOfAtomTwo){
      String temp = atomOneName;
      atomOneName = atomTwoName;
      atomTwoName = temp;
    }
    // construct the Rotation Transformation
    Transformation t = new Transformation();
    Vector3D a1 = getAtom(atomOneName).getCoordinates();
    Vector3D a2 = getAtom(atomTwoName).getCoordinates();
    t.addRotationAboutAxis(a1, a2.subtract(a1), degrees);
    // get the atoms that need to be rotated, all those after this bond away from the backbone
    Collection<Atom> atomsToTransform = getAtomsAfterAtom(atomTwoName);
    for(Atom atom : atomsToTransform){
      atom.applyTransformation(t);
    }
  }

  private boolean isRotatableBond(String atomOneName, String atomTwoName){
    // if the atoms of the bond are not in definedDihedralAngles, then they are not a rotatable bond
    if(!this.definedDihedralAngles.contains(atomOneName)
        || !this.definedDihedralAngles.contains(atomTwoName)){
      return false;
    }
    int indexOfAtomOne = this.definedDihedralAngles.indexOf(atomOneName);
    int indexOfAtomTwo = this.definedDihedralAngles.indexOf(atomTwoName);
    // the first and last elements in definedDihedralAngles are not rotatable bond atoms
    String first = this.definedDihedralAngles.get(0);
    String last = this.definedDihedralAngles.get(this.definedDihedralAngles.size()-1);
    // if the atoms are either the first or last atoms in definedDihedralAngles, then they are not
    // a rotable bond
    if(atomOneName.equals(first) || atomTwoName.equals(first)
        || atomOneName.equals(last) || atomTwoName.equals(last)){
      return false;
    }
    // if they are not consecutive in definedDihedralAngles, then they are not a rotatable bond
    if(Math.abs(indexOfAtomOne - indexOfAtomTwo) != 1){
      return false;
    }
    // if the atom before and the atom after this bond is not in this residue, then this is
    // not a rotatable bond
    if(indexOfAtomOne < indexOfAtomTwo){
      String before = this.definedDihedralAngles.get(indexOfAtomOne-1);
      String after = this.definedDihedralAngles.get(indexOfAtomTwo+1);
      if(!this.contains(before) || !this.contains(after)){
        return false;
      }
    } else {
      String before = this.definedDihedralAngles.get(indexOfAtomTwo-1);
      String after = this.definedDihedralAngles.get(indexOfAtomOne+1);
      if(!this.contains(before) || !this.contains(after)){
        return false;
      }
    }
    return true;
  }

  /**
  * Return a Collection of atoms after and including the atom provided as an argument
  * @param atomName the name of the atom to get all the atoms after it in the residue. These are the
  * atoms bonded to it extending away from CA.
  * @return a collection of those atoms.
  * @throws NoSuchElementException if there is no Atom with atomName in this residue
  */
  private Collection<Atom> getAtomsAfterAtom(String atomName){
    if(!this.contains(atomName)){
      throw new NoSuchElementException(String.format("Atom %s not in residue.", atomName));
    }
    return this.atomsGraph.depthFirstSearch(this.getAtom(atomName));
  }

  /**
  * Get the distance between two atoms in this residue
  * @param atomOneName one atom in this residue
  * @param atomTwoName the other atom in this residue
  * @return the distance between these two atoms
  * @since 0.7.0
  */
  public double getDistance(String atomOneName, String atomTwoName){
    return this.getAtom(atomOneName).distance(this.getAtom(atomTwoName));
  }

  /**
  * Get the value for one of the dihedral angles in this Residue.
  * @param atomOneName the name of the first atom in the bond, the atom closer to the backbone
  * @param atomTwoName the name of the second atom in the bond, the atom further from the backbone
  * @return the degree value of the dihedral angle about this bond, or 1000 if any of the atoms
  * required to calculate the angle are missing.
  * @since 0.7.0
  */
  public double getDihedralAngle(String atomOneName, String atomTwoName){
    try{
      int indexOfAtomOne = this.definedDihedralAngles.indexOf(atomOneName);
      int indexOfAtomTwo = this.definedDihedralAngles.indexOf(atomTwoName);
      int indexOfB = indexOfAtomOne;
      if(indexOfAtomOne > indexOfAtomTwo){
        indexOfB = indexOfAtomTwo;
        String temp = atomOneName;
        atomOneName = atomTwoName;
        atomTwoName = temp;
      }
      Vector3D a = this.getAtom(this.definedDihedralAngles.get(indexOfB-1)).getCoordinates();
      Vector3D b = this.getAtom(atomOneName).getCoordinates();
      Vector3D c = this.getAtom(atomTwoName).getCoordinates();
      Vector3D d = this.getAtom(this.definedDihedralAngles.get(indexOfB+2)).getCoordinates();
      return -1 * Vector3D.calculateDihedralAngle(a,b,c,d);
    } catch(NoSuchElementException e){
      // one of the atoms is missing from the residue, can't calculate a dihedral angle
      // return 1000
      return 1000;
    }
  }

  /**
  * Return the angle formed by 3 atoms in this residue. The angle is that between atomOne and
  * atomThree using atomTwo as the pivot.
  * @param atomOneName one of the atoms in the angle
  * @param atomTwoName the pivot atom
  * @param atomThree the other atom in the angle
  * @return the angle defined by atomOne-atomTwo-atomThree, in degrees
  * @since 0.7.0
  */
  public double getAngle(String atomOneName, String atomTwoName, String atomThree){
    Vector3D one = getAtom(atomOneName).getCoordinates();
    Vector3D two = getAtom(atomTwoName).getCoordinates();
    Vector3D three = getAtom(atomThree).getCoordinates();
    return one.subtract(two).angle(three.subtract(two));
  }

  private double getAngle(Atom atomOne, Atom atomTwo, Atom atomThree){
    return getAngle(atomOne.getAtomName(), atomTwo.getAtomName(), atomThree.getAtomName());
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

  private class ResidueBondComparator implements Comparator<Bond>{
    public int compare(Bond a, Bond b){
      int indexA1 = atomsInOrder.indexOf(a.getAtomOne());
      int indexA2 = atomsInOrder.indexOf(a.getAtomTwo());
      int indexB1 = atomsInOrder.indexOf(b.getAtomOne());
      int indexB2 = atomsInOrder.indexOf(b.getAtomTwo());
      if(indexA1 - indexB1 == 0){
        return indexA2 - indexB2;
      } else {
        return indexA1 - indexB1;
      }
    }
  }

  private class AngleTriple{
    private ArrayList<Atom> atoms;
    public AngleTriple(Atom one, Atom two, Atom three){
      if(one.equals(two) || one.equals(three) || two.equals(three)){
        throw new IllegalArgumentException(
            String.format("Recieved %s, %s, %s\nAtoms must be different",
                          one.getAtomName(), two.getAtomName(), three.getAtomName()));
      }
      atoms = new ArrayList<Atom>();
      atoms.add(one);
      atoms.add(two);
      atoms.add(three);
      Collections.sort(atoms, new Comparator<Atom>(){
        public int compare(Atom a, Atom b){
          int indexA = atomsInOrder.indexOf(a);
          int indexB = atomsInOrder.indexOf(b);
          return indexA - indexB;
        }
      });
    }

    public Atom getAtomOne(){
      return atoms.get(0);
    }

    public Atom getAtomTwo(){
      return atoms.get(1);
    }

    public Atom getAtomThree(){
      return atoms.get(2);
    }

    private boolean contains(Atom atom){
      return atoms.contains(atom);
    }

    @Override
    public int hashCode(){
      return String.format("%s%s%s",getAtomOne(), getAtomTwo(), getAtomThree()).hashCode();
    }
    @Override
    public boolean equals(Object obj){
      if(obj instanceof AngleTriple){
        AngleTriple other = (AngleTriple)obj;
        if(other.contains(this.getAtomOne())
            && other.contains(this.getAtomTwo())
            && other.contains(this.getAtomThree())){
          return true;
        }
      }
      return false;
    }
  }
}
