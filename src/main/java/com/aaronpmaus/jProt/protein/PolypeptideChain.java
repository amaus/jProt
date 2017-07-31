package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Collection;

/**
 * A PolypeptideChain is a single chain of amino acids
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
*/

public class PolypeptideChain extends Molecule{
    private ArrayList<Residue> residues;
    // key: residueID, value: residue index in residues list
    private HashMap<Integer, Integer> residueIndexLookupTable;
    private String chainID;

    /**
     * A constructor for a PolypeptideChain.
     * @param chainID the letter ID of this chain
    */
    public PolypeptideChain(String chainID){
        super();
        residues = new ArrayList<Residue>();
        residueIndexLookupTable = new HashMap<Integer,Integer>();
        this.chainID = chainID;
    }

    /**
     * Returns the Chain ID of this chain
     * @return the String that is the Chain ID: A, B, C, etc...
    */
    public String getChainID(){
        return this.chainID;
    }

    /**
     * Add a residue to this chain
     * @param residue the residue to add
    */
    public void addResidue(Residue residue){
        int residueIndex = residues.size();
        residues.add(residue);
        residueIndexLookupTable.put(residue.getResidueID(), residueIndex);
        Collection<Bond> bonds = residue.getBonds();
        for(Bond b : bonds){
            this.addBond(b); // inherited from Molecule
        }
        // if the previous residue - the residue with the residueID immediately
        // preceeding this residue's ID - exists, then add a bond between it's
        // C and this residue's N.
        if(this.containsResidue(residue.getResidueID()-1)){
            Residue prevResidue = getResidue(residue.getResidueID()-1);
            this.addBond(new Bond(prevResidue.getAtom("C"), residue.getAtom("N"), 1));
        }
        // if the next residue - the residue with the residueID immediately
        // after this residue's ID - exists, then add a bond between this
        // residue's C and its N.
        if(this.containsResidue(residue.getResidueID()+1)){
            Residue nextResidue = getResidue(residue.getResidueID()+1);
            this.addBond(new Bond(residue.getAtom("C"), nextResidue.getAtom("N"), 1));
        }
    }

    /**
     * Gets the residue with the given residueID.
     * @param residueID the residueID of the residue to get
     * @return the Residue if it is in this PolypeptideChain, null otherwise
    */
    public Residue getResidue(int residueID){
        if(this.containsResidue(residueID)){
            int residueIndex = residueIndexLookupTable.get(residueID);
            return this.residues.get(residueIndex);
        }
        return null;
    }

    /**
     * Checks if there is a residue with the given ID in this PolypeptideChain
     * @param residueID the numeric residue ID of the residue to check for.
     * @return true if there is a residue with the given ID in this
     *         PolypeptideChain, false otherwise
    */
    public boolean containsResidue(int residueID){
        return residueIndexLookupTable.containsKey(residueID);
    }

    /**
     * Returns the number of residues in this PolypeptideChain
     * @return the  number of residues in this PolypeptideChain
    */
    public int getNumResidues(){
        return residues.size();
    }

    /**
     * Returns the number of Atoms in this PolypeptideChain
     * @return the number of Atoms in this PolypeptideChain
    */
    public int getNumAtoms(){
        int numAtoms = 0;
        for(Residue res : this.residues){
            numAtoms += res.getNumAtoms();
        }
        return numAtoms;
    }
}
