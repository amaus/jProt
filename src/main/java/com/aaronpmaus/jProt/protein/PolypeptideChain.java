package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import java.util.HashMap;
import java.util.ArrayList;
/**
 * A PolypeptideChain is a single chain of amino acids
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
*/

public class PolypeptideChain extends Molecule{
    private ArrayList<AminoAcid> residues;
    // key - resID, value - index into residues list
    private HashMap<Integer,Integer> residueIDLookUpTable;
    private String chainID;

    /**
     * A constructor for a PolypeptideChain.
    */
    public PolypeptideChain(){
        super();
        residues = new ArrayList<AminoAcid>();
        residueIDLookUpTable = new HashMap<Integer, Integer>();
        chainID = "A";
    }

    /**
     * Returns the Chain ID of this chain
     * @return the String that is the Chain ID: A, B, C, etc...
    */
    public String getChainID(){
        return this.chainID;
    }
}
