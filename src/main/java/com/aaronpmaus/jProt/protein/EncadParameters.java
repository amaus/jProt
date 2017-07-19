package com.aaronpmaus.jProt.protein;
import java.util.HashMap;
/**
 * A list of commonly used parameters for computational chemistry.
 * Includes bond lengths.
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
*/
public class EncadParameters{
    private static HashMap<String,Double> bondLengths = initializeBondLengths();

    /*
     * A helper method to initialize the bond lengths.
     * @return a HashMap where the keys are atom bond strings, and the values
     * are their bond lengths.
    */
    private static HashMap<String,Double> initializeBondLengths(){
        HashMap<String,Double> bondLengths = new HashMap<String,Double>();
        bondLengths.put("O-H", 1.0000);
        bondLengths.put("N-H", 1.0000);

        bondLengths.put("C-H", 1.09);
        bondLengths.put("C-O", 1.437);
        bondLengths.put("C-N", 1.4670);
        bondLengths.put("C-C", 1.5250);
        bondLengths.put("C-S", 1.8080);
        return bondLengths;

    }

    /**
     * Returns the bond length of the two Atoms
     * @param a1 one of the atoms in the bond
     * @param a2 the other atom in the bond
     * @return the length of their covalent bond
    */
    public static double getBondLength(Atom a1, Atom a2){
        return bondLengths.get( a1.getElement()+"-"+a2.getElement());
    }

    /**
     * Returns the bond force constant of a bond containing these two atoms.
     * @param a1 one of the atoms in the bond
     * @param a2 the other atom in the bond
     * @return the bond force constant of their covalent bond in kcal/mol-A**2
    */
    public static double getBondForceConstant(Atom a1, Atom a2){
        return 250.0;
    }
}
