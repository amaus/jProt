package com.aaronpmaus.jProt.tools;

import com.aaronpmaus.jProt.protein.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collection;

/**
* This tool can be used to calculate distance matrices for Proteins and Polypeptide Chains.
* <p>
* To calculate a distance matrix for a protein or a chain:
* <pre>
* <code>
* Protein prot1 = PDBFileIO.readInPDBFile(streamToProtFile, "1ABC");
* // get a CA-distance matrix for prot1
* Double[][] protMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot1);
* // get a CA-distance matrix for chain "A" of prot1
* Double[][] chainMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot1.getChain("A"));
* </code>
* </pre>
* <p>
* The masked versions of the methods allow for residues to be selected to be used in the
* distance calculations. This is intended to be used in concert with the SequenceAligner
* to be able to calculate distance matrices for separate proteins using only those
* residues in each that have been aligned. These matrices can then be used to directly compare
* the two structures.
* <p>
* For example:
* <pre>
* <code>
* Protein prot1 = PDBFileIO.readInPDBFile(streamToProt1, "prot1");
* Protein prot2 = PDBFileIO.readInPDBFile(streamToProt2, "prot2");
*
* // first, get an alignment of the sequences of these proteins
* String[] alignment = SequenceAligner.alignProteinSequences(prot1.getSequence(),
*                                                           prot2.getSequence());
* // get masks which indicate which residues in each protein have a match in the other protein.
* boolean[][] masks = SequenceAligner.getSequenceMatchMasks(alignment[0], alignment[1]);
* boolean[] prot1Mask = masks[0];
* boolean[] prot2Mask = masks[1];
*
* Double[][] prot1DistanceMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot1, prot1Mask);
* Double[][] prot2DistanceMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot2, prot2Mask);
* // Now prot1DistanceMatrix and prot2DistanceMatrix can be used directly to compare the two
* // structures
* </code>
* </pre>
*
* @version 0.6.0
* @since 0.6.0
*/
public class DistanceMatrixCalculator{

  private static Double[][] calculateDistanceMatrix(Collection<Atom> atoms){
    int numAtoms = atoms.size();
    Double[][] distanceMatrix = new Double[numAtoms][numAtoms];
    int i = 0;
    int j = 0;
    for(Atom atomOne : atoms){
      for(Atom atomTwo : atoms){
        distanceMatrix[i][j] = atomOne.distance(atomTwo);
        j++;
      }
      j = 0;
      i++;
    }
    return distanceMatrix;
  }

  /**
  * Calculate and return the CA Distance Matrix of the protein using all Residues
  *
  * @param prot the protein to calculate the CA Distance Matrix of
  * @return a 2D array of Double containing the CA distances
  * @since 0.6.0
  */
  public static Double[][] calculateDistanceMatrix(Protein prot){
    boolean[] mask = new boolean[prot.getNumResidues()];
    Arrays.fill(mask, true);
    return calculateDistanceMatrix(prot, mask);
  }

  /**
  * Calculate and return the CA Distance Matrix of the protein using only the residues
  * specified by the mask.
  *
  * The mask must contain as many elements as there are residues in the protein. Each
  * residue in order will have a value true or false. If true, include that residue
  * in the distance matrix, otherwise, ignore it.
  *
  * @param prot the protein to calculate the CA Distance Matrix of
  * @param mask an array of boolean containing as many values as there are residues, each
  *  indicating whether to include that residue in the distance matrix.
  * @return a 2D array of Double containing the CA distances
  * @since 0.6.0
  */
  public static Double[][] calculateDistanceMatrix(Protein prot, boolean[] mask){
    ArrayList<Atom> atoms = new ArrayList<Atom>(prot.getNumResidues());
    int maskIndex = 0;
    for(PolypeptideChain chain : prot){
      for(Residue res : chain){
        if(mask[maskIndex]){
          atoms.add(res.getAtom("CA"));
        }
        maskIndex++;
      }
    }
    return calculateDistanceMatrix(atoms);
  }

  /**
  * Calculate and return the CA Distance Matrix of the Chain using all Residues
  *
  * @param chain the Chain to calculate the CA Distance Matrix of
  * @return a 2D array of Double containing the CA distances
  * @since 0.6.0
  */
  public static Double[][] calculateDistanceMatrix(PolypeptideChain chain){
    boolean[] mask = new boolean[chain.getNumResidues()];
    Arrays.fill(mask, true);
    return calculateDistanceMatrix(chain, mask);
  }

  /**
  * Calculate and return the CA Distance Matrix of the Chain using only the residues
  * specified by the mask.
  *
  * The mask must contain as many elements as there are residues in the Chain. Each
  * residue in order will have a value true or false. If true, include that residue
  * in the distance matrix, otherwise, ignore it.
  *
  * @param chain the protein to calculate the CA Distance Matrix of
  * @param mask an array of boolean containing as many values as there are residues, each
  *  indicating whether to include that residue in the distance matrix.
  * @return a 2D array of Double containing the CA distances
  * @since 0.6.0
  */
  public static Double[][] calculateDistanceMatrix(PolypeptideChain chain, boolean[] mask){
    ArrayList<Atom> atoms = new ArrayList<Atom>(chain.getNumResidues());
    int maskIndex = 0;
    for(Residue res : chain){
      if(mask[maskIndex]){
        atoms.add(res.getAtom("CA"));
      }
      maskIndex++;
    }
    return calculateDistanceMatrix(atoms);
  }

}
