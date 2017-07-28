package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import java.util.HashMap;

/**
 * This class represents in general an AminoAcid.
 * <p>
 * An AminoAcid consists of a collection of Atoms bonded together
 * is a particular way.
 * <p>
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
*/

public class Residue{
    private final String name;
    private final String threeLetterID;
    private final String oneLetterID;
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
                    {"ASN","N","Asparagine"},     {"GLN","Q","Glutamine"}};

    public Residue(String threeLetterID, int residueID){
        this(threeLetterID,residueID,null);
    }
    /**
     * A constructor for an amino acid.
     * @param threeLetterID the one letter ID of the amino acid to build.
    */
    public Residue(String threeLetterID, int residueID, Collection<Atom> atoms){
        this.threeLetterID = threeLetterID;
        this.oneLetterID = Residue.lookUpOneLetterName(this.threeLetterID);
        this.name = Residue.lookUpFullName(this.threeLetterID);
        this.residueID = residueID;
        initializeAminoAcid(this.threeLetterID+".dat", atoms);
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
     * get the residue ID of this residue
     * @return the residue ID of this residue
    */
    public int getResidueID(){
        return this.residueID;
    }

    /**
     * get the set of covalent bonds in this residue
     * @return a {@code HashSet<Bond>} of the bonds in this residue
    */
    public HashSet<Bond> getBonds(){
        return this.bonds;
    }

    private void initializeAminoAcid(String dataFileName, Collection<Atom> atoms){
        this.residueComplete = true;
        // if atoms is null, build amino acid from default values
        if(atoms == null){

        } else {
            // use the atoms passed in
            for(Atom a: atoms){
                this.atoms.put(a.getName(), a);
            }
            residueComplete = addBond("N","CA",1) & residueComplete;
            residueComplete = addBond("CA","CB",1) & residueComplete;
            residueComplete = addBond("CA","C",1) & residueComplete;
            residueComplete = addBond("C","O",2) & residueComplete;
        }
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
    private void initializeAlanine(Collection<Atom> atoms){
        this.residueComplete = true;
        // if atoms is null, build amino acid from default values
        if(atoms == null){

        } else {
            // use the atoms passed in
            for(Atom a: atoms){
                this.atoms.put(a.getName(), a);
            }
            residueComplete = addBond("N","CA",1) & residueComplete;
            residueComplete = addBond("CA","CB",1) & residueComplete;
            residueComplete = addBond("CA","C",1) & residueComplete;
            residueComplete = addBond("C","O",2) & residueComplete;
        }

    }

    private void initializeGlycine(Collection<Atom> atoms){
        this.residueComplete = true;
        // if atoms is null, build amino acid from default values
        if(atoms == null){

        } else {
            // use the atoms passed in
            for(Atom a: atoms){
                this.atoms.put(a.getName(), a);
            }
            residueComplete = addBond("N","CA",1) & residueComplete;
            residueComplete = addBond("CA","C",1) & residueComplete;
            residueComplete = addBond("C","O",2) & residueComplete;
        }

    }

    private void initializeIsoluecine(Collection<Atom> atoms){
        this.residueComplete = true;
        // if atoms is null, build amino acid from default values
        if(atoms == null){

        } else {
            // use the atoms passed in
            for(Atom a: atoms){
                this.atoms.put(a.getName(), a);
            }
            residueComplete = addBond("N","CA",1) & residueComplete;
            residueComplete = addBond("CA","CB",1) & residueComplete;
            residueComplete = addBond("CB","CG1",1) & residueComplete;
            residueComplete = addBond("CB","CG2",1) & residueComplete;
            residueComplete = addBond("CG1","CD1",1) & residueComplete;
            residueComplete = addBond("CA","C",1) & residueComplete;
            residueComplete = addBond("C","O",2) & residueComplete;
        }

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
}
