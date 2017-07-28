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
    */
    public Protein readInPDBFile(String fileName) throws FileNotFoundException{
        Scanner in = new Scanner(new File(fileName));
        PolypeptideChain chain = new PolypeptideChain();
        boolean firstAtomOfChain = true;
        while(in.hasNextLine()){
            String line = in.nextLine();
            if(line.length() < 80){
                line = padRight(line);
            }
            String recordName = line.substring(0,6).trim();
            switch(recordName){
                case "ATOM":
                    Atom atom = parseAtomLine(line);
                    //chain.addAtom(atom);
                    break;
            }
        }
        return null;
    }

    private Atom parseAtomLine(String line){
        int serial = Integer.parseInt(line.substring(6,11).trim());
        String atomName = line.substring(12,16).trim();
        String altLoc = line.substring(16,17).trim();
        String resName = line.substring(17,20).trim();
        String chainID = line.substring(21,22).trim();
        int resSeq = Integer.parseInt(line.substring(22,26).trim());
        String iCode = line.substring(26,27).trim(); // code for insertion of residues
        double x = Double.parseDouble(line.substring(30,38).trim());
        double y = Double.parseDouble(line.substring(38,46).trim());
        double z = Double.parseDouble(line.substring(46,54).trim());

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
