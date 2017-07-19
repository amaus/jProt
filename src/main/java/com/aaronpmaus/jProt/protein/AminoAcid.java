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

public class AminoAcid extends Molecule{
    private final String name;
    private final String threeLetterID;
    private final String oneLetterID;

    private boolean hydrophobic;
    private HashMap<String, Atom> atoms;

    /**
     * A constructor for an amino acid.
     * @param oneLetterID the one letter ID of the amino acid to build.
    */
    public AminoAcid(String oneLetterID){
        super();
        this.name = null;
        this.threeLetterID = null;
        this.oneLetterID = oneLetterID;
    }
}
