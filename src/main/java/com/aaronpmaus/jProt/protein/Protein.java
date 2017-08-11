package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.graph.*;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * A Protein has one or more PolypeptideChains
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
 */
 public class Protein {
   private ArrayList<PolypeptideChain> chains;
   private String pdbName;

   /**
   * A constructor for a PolypeptideChain.
   * @param pdbName the name of the PDB file for this protein
   */
   public Protein(String pdbName){
     chains = new ArrayList<PolypeptideChain>();
     this.pdbName = pdbName;
   }

   /**
   * A constructor for a protein that takes in the chains of this protein.
   * @param pdbName the name of the PDB file for this protein
   * @param chains one or more chains that make up this protein
   */
   public Protein(String pdbName, PolypeptideChain... chains){
     this(pdbName);
     int index = 0;
     for(PolypeptideChain chain: chains){
       this.chains.add(chain);
     }
   }

   /**
   * Return the name of the PDB file for this Protein
   * @return the name of the PDB file without the .pdb extension.
   */
   public String getPDBName(){
     return this.pdbName;
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
     for(PolypeptideChain chain : this.chains){
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
     for(PolypeptideChain chain : this.chains){
       numAtoms += chain.getNumAtoms();
     }
     return numAtoms;
   }

   /**
   * Return the sequence of this Protein in the form of a String of the single letter
   * residue names.
   *
   * @return a String containing the single letter residue names with no spaces
   */
   public String getSequence(){
     String seq = "";
     for(PolypeptideChain chain : this.chains){
       seq += chain.getSequence();
     }
     return seq;
   }

   /**
   * Calculate and return the CA Distance Matrix of this protein
   *
   * @return a 2D array of Double containing the CA distances
   */
   public Double[][] calculateCarbonAlphaDistanceMatrix(){
     boolean[] mask = new boolean[getNumResidues()];
     Arrays.fill(mask, true);
     return calculateCarbonAlphaDistanceMatrix(mask);
   }

   /**
   * Calculate and return the CA Distance Matrix of this protein using only the residues
   * specified by the mask.
   *
   * The mask must contain as many elements as there are residues in this protein. Each
   * residue in order will have a value true or false. If true, include that residue
   * in the distance matrix, otherwise, ignore it.
   *
   * @param mask an array of boolean containing as many values as there are residues, each
   *  indicating whether to include that residue in the distance matrix.
   * @return a 2D array of Double containing the CA distances
   * @since 0.10.0
   */
   public Double[][] calculateCarbonAlphaDistanceMatrix(boolean[] mask){
     ArrayList<Residue> residues = new ArrayList<Residue>();
     int maskIndex = 0;
     for(PolypeptideChain chain : this.chains){
       for(Residue res : chain){
         if(mask[maskIndex]){
           residues.add(res);
         }
         maskIndex++;
       }
     }
     int numResidues = residues.size();
     Double[][] distanceMatrix = new Double[numResidues][numResidues];
     int i = 0;
     int j = 0;
     for(Residue residueOne : residues){
       for(Residue residueTwo : residues){
         Atom carbonAlphaOne = residueOne.getAtom("CA");
         Atom carbonAlphaTwo = residueTwo.getAtom("CA");
         distanceMatrix[i][j] = carbonAlphaOne.distance(carbonAlphaTwo);
         j++;
       }
       j = 0;
       i++;
     }
     return distanceMatrix;
   }

   /**
   * Calculates and returns the CA Distance Matrix of one of the chains in
   * this protein
   *
   * @param chainID the ID of the chain to get the CA Distance Matrix of
   * @return a 2D array of Double containing the CA distances
   */
   public Double[][] calculateCarbonAlphaDistanceMatrix(String chainID){
     PolypeptideChain chain = getChain(chainID);
     return chain.calculateCarbonAlphaDistanceMatrix();
   }

   /**
   * Return the residue IDs of all the residues in this Protein
   *
   * @return an array of Integers holding the residue IDs
   */
   public Integer[] getResidueIDs(){
     boolean[] mask = new boolean[getNumResidues()];
     Arrays.fill(mask, true);
     return getResidueIDs(mask);
   }

   /**
   * Return the residue IDs of all the residues in this Protein specified by the mask.
   *
   * The mask must contain as many elements as there are residues in this protein. Each
   * residue in order will have a value true or false. If true, include that residue
   * ID in the array, otherwise, ignore it.
   *
   * @param mask an array of boolean containing as many values as there are residues, each
   *  indicating whether to include that residue in the IDs list.
   * @return an array of Integers holding the residue IDs
   */
   public Integer[] getResidueIDs(boolean[] mask){
     int numResidues = 0;
     for(boolean val : mask){
       if(val){
         numResidues++;
       }
     }
     Integer[] residueIDs = new Integer[numResidues];
     int residueArrayIndex = 0;
     int maskIndex = 0;
     for(PolypeptideChain chain : this.chains){
       for(Integer id : getResidueIDs(chain.getChainID())){
         if(mask[maskIndex]){
           residueIDs[residueArrayIndex] = id;
           residueArrayIndex++;
         }
         maskIndex++;
       }
     }
     return residueIDs;
   }

   /**
   * Returns the residue IDs of all the residues in the chain with the
   * chainID provided as an argument
   * @param chainID the letter ID of the chain to get the residue IDs from
   * @return an array of Integers holding the residue IDs
   */
   public Integer[] getResidueIDs(String chainID){
     return getChain(chainID).getResidueIDs();
   }
 }
