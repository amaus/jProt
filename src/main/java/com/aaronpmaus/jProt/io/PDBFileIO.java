package com.aaronpmaus.jProt.io;

import com.aaronpmaus.jProt.protein.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

/**
 * Provides the ability to read in and write out PDB Files
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
*/
public class PDBFileIO{
    ArrayList<String> fileLines;


    public PDBFileIO(){
        fileLines = new ArrayList<String>();
    }

    /**
     * Method to read in a PDB file and return a Protein
     * @param fileName the name of the PDB file
     * @return a Protein built from that file
     * @throws FileNotFoundException if the file can not be read from
    */
    public Protein readInPDBFile(String fileName) throws FileNotFoundException{
        Scanner in = new Scanner(new File(fileName));
        PolypeptideChain currentChain = null;
        boolean firstAtomOfChain = true;
        int currentResidueID = -1;
        int resSeq = -1;
        String resName = null;
        ArrayList<Atom> residueAtoms = null;
        String[] fileNameParts = fileName.split("\\.");
        String fileBase = fileNameParts[0];
        System.out.println("reading in PDB: " + fileBase);
        Protein protein = new Protein(fileBase);

        // bookkeeping to make sure that the number of chains that has
        // been instantiated is the same number that has been added
        // by the time the whole PDB is read in. This is because some malformed
        // PDBs leave off the final TER record. GRRR
        int numChainsInstantiated = 0;
        while(in.hasNextLine()){
            String line = in.nextLine();
            if(line.length() < 80){
                line = padRight(line);
            }
            String recordName = line.substring(0,6).trim();
            switch(recordName){
                case "ATOM":
                    // if currentResidueID is -1, then this is the first atom
                    // of the first residue in this chain. Get the chainID
                    // and build the chain
                    if(currentResidueID == -1){
                        String chainID = line.substring(21,22).trim();
                        currentChain = new PolypeptideChain(chainID);
                        numChainsInstantiated++;
                    }
                    // what residue does this Atom belong to?
                    resSeq = Integer.parseInt(line.substring(22,26).trim());
                    resName = line.substring(17,20).trim();
                    // if this is a new residue number,
                    if(resSeq != currentResidueID){
                        // and if this isn't the first atom of the first residue
                        if(currentResidueID != -1){
                            // build a new residue from the atoms read in and
                            // add it to the chain
                            Residue res = new Residue(resName, resSeq, residueAtoms);
                            currentChain.addResidue(res);
                        }
                        // save the new residue number
                        currentResidueID = resSeq;
                        // create an empty ArrayList to hold this residue's
                        // atoms
                        residueAtoms = new ArrayList<Atom>();
                        // build this atom and add it to the list
                        Atom atom = parseAtomLine(line);
                        residueAtoms.add(atom);

                    } else { // this atom belongs to the same residue as prev
                        // build the atom and add it to the list
                        Atom atom = parseAtomLine(line);
                        residueAtoms.add(atom);
                    }
                    //chain.addAtom(atom);
                    break;
                case "TER":
                    // end of a chain
                    // add last residue to this chain
                    Residue res = new Residue(resName, resSeq, residueAtoms);
                    currentChain.addResidue(res);
                    // add this chain to the protein
                    protein.addChain(currentChain);
                    // reset the currentResidueID back to -1 to indicate
                    // that no residues have been read in for the current chain
                    currentResidueID = -1;
                    break;
                case "END":
                    // if the number of chains in the protein don't match the number of
                    // chains instantiated, then we need to add the last chain to the
                    // protein
                    if(protein.getNumChains() != numChainsInstantiated){
                        // add the last residue to the last chain
                        res = new Residue(resName, resSeq, residueAtoms);
                        currentChain.addResidue(res);
                        // and add the last chain to the protein
                        protein.addChain(currentChain);
                    }
            }
        }

        return protein;
    }

    private Atom parseAtomLine(String line){
        int serial = Integer.parseInt(line.substring(6,11).trim());
        String atomName = line.substring(12,16).trim();
        String altLoc = line.substring(16,17).trim();
        String resName = line.substring(17,20).trim();
        String chainID = line.substring(21,22).trim();
        int resSeq = Integer.parseInt(line.substring(22,26).trim());
        String iCode = line.substring(26,27).trim(); // code for insertion of residues
        String x = line.substring(30,38).trim();
        String y = line.substring(38,46).trim();
        String z = line.substring(46,54).trim();

        // make sure there is a value for occupancy. if so, assign.
        // else, default value of -1.0.
        String occupancyString = line.substring(54,60).trim();
        double occupancy = -1.0;
        if(!occupancyString.equals("")){
            occupancy = Double.parseDouble(occupancyString);
        }

        // make sure there is a value for temp factor. if so, assign.
        // else, default value of -1.0.
        String tempFactorString = line.substring(60,66).trim();
        double tempFactor = -1.0;
        if(!tempFactorString.equals("")){
            tempFactor = Double.parseDouble(tempFactorString);
        }
        String element = line.substring(76,78).trim();

        // make sure there is a value for charge. if so, parse and assign.
        // else default charge of 0.
        // format for charge is 1+ or 2+ or 1-.
        String chargeString = line.substring(78,80).trim();
        double charge = 0;
        if(!chargeString.equals("")){
            charge = Double.parseDouble(chargeString.substring(0,1));
            if(chargeString.substring(1,2).equals("-")){
                charge *= -1;
            }
        }
        return new Atom(atomName, serial, occupancy, tempFactor, charge, x, y, z);
    }

    private String padRight(String line){
        return String.format("%1$-" + (80-line.length()) + "s", line);
    }
}
