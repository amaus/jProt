package com.aaronpmaus.jProt.tools;

import com.aaronpmaus.jProt.protein.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collection;

/**
* This tool can be used to calculate distance matrices for Proteins and Polypeptide Chains.
* <p>
* To calculate a distance matrix for a protein or a chain:
* <p>
* {@code InputStream in = new FileInputStream(new File(pathToPDBFile));}<br>
* {@code Protein prot1 = PDBFileIO.readInPDBFile(in, "1ABC");}<br>
* {@code // get a CA-distance matrix for prot1}<br>
* {@code Double[][] protMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot1);}<br>
* {@code // get a CA-distance matrix for chain "A" of prot1}<br>
* {@code Double[][] chainMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot1.getChain("A"));}<br>
* <p>
* The masked versions of the methods allow for residues to be selected to be used in the
* distance calculations. This is intended to be used in concert with the SequenceAligner
* to be able to calculate distance matrices for separate proteins using only those
* residues in each that have been aligned. These matrices can then be used to directly compare
* the two structures.
* <p>
* For example:
* <p>
* {@code InputStream prot1Stream = new FileInputStream(new File(pathToProt1));}<br>
* {@code InputStream prot2Stream = new FileInputStream(new File(pathToProt2));}<br>
* {@code Protein prot1 = PDBFileIO.readInPDBFile(prot1Stream, "prot1");}<br>
* {@code Protein prot2 = PDBFileIO.readInPDBFile(prot2Stream, "prot2");}<br>
* <br>
* {@code // first, get an alignment of the sequences of these proteins}<br>
* {@code String[] alignment = SequenceAligner.alignProteinSequences(prot1.getSequence(),}<br>
* {@code     prot2.getSequence());}<br>
* {@code // get masks which indicate which residues in each protein have a match in the other protein.}<br>
* {@code boolean[][] masks = SequenceAligner.getSequenceMatchMasks(alignment[0], alignment[1]);}<br>
* {@code boolean[] prot1Mask = masks[0];}<br>
* {@code boolean[] prot2Mask = masks[1];}<br>
* <br>
* {@code Double[][] prot1DistanceMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot1, prot1Mask);}<br>
* {@code Double[][] prot2DistanceMatrix = DistanceMatrixCalculator.calculateDistanceMatrix(prot2, prot2Mask);}<br>
* {@code // Now prot1DistanceMatrix and prot2DistanceMatrix can be used directly to compare the two}<br>
* {@code // structures}<br>
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
