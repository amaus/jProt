package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import java.util.HashMap;
import java.util.ArrayList;
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
            this.chains.add(chain);
            index++;
        }
    }

    /**
     * Add a chain to this protein
     * @param chain the PolypeptideChain to add to this protein
    */
    public void addChain(PolypeptideChain chain){
        this.chains.add(chain);
    }

    /**
     * Returns a reference to a chain in the protein
     * @param chainID the ID of the chain to get
     * @return the PolypeptideChain with the ID provided as an argument
    */
    public PolypeptideChain getChain(String chainID){
        for(PolypeptideChain chain : chains){
            if(chain.getChainID().equals(chainID)){
                return chain;
            }
        }
        throw new IllegalArgumentException("There is no chain "
                    +chainID+ " in the protein.");
    }

    /**
     * Returns the number of Chains in this Protein
     * @return the number of chains in this protein
    */
    public int getNumChains(){
        return this.chains.size();
    }

    /**
     * Returns the number of Residues in this Protein
     * @return the number of Residues in this protein
    */
    public int getNumResidues(){
        int numResidues = 0;
        for(PolypeptideChain chain : chains){
            numResidues += chain.getNumResidues();
        }
        return numResidues;
    }

    /**
     * Returns the number of Atoms in this Protein
     * @return the number of Atoms in this Protein
    */
    public int getNumAtoms(){
        int numAtoms = 0;
        for(PolypeptideChain chain : chains){
            numAtoms += chain.getNumAtoms();
        }
        return numAtoms;
    }
}
