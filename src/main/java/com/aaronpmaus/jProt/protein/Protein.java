package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import java.util.HashMap;
/**
 * A Protein has one or more PolypeptideChains
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
*/

public class Protein {
    private ArrayList<PolypeptideChain> chains;
    // key - chainID, value - index into residues list
    private HashMap<String,Integer> chainIDLookUpTable;

    /**
     * A constructor for a PolypeptideChain.
    */
    public Protein(){
        super();
        chains = new ArrayList<PolypeptideChain>();
        chainIDLookUpTable = new HashMap<String, Integer>();
    }

    /**
     * A constructor for a protein that takes in the chains of this protein.
     * @param chains one or more chains that make up this protein
    */
    public Protein(PolypeptideChain... chains){
        int index = 0;
        for(PolypeptideChain chain: chains){
            chainIDLookUpTable.put(chain.getChainID(), index);
            chains.add(chain);
            index++;
        }
    }
}
