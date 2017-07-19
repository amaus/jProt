package com.aaronpmaus.jProt.protein;
import java.lang.IllegalArgumentException;

// make this an Observer of the atoms in it? That way if one of them
// is moved, the bond will atomatically update its distance and energy?
/**
 * Class Bond to represent an atomic bond.
 * @author Aarom Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.6.0
*/
public class Bond {
    private Atom a1;
    private Atom a2;
    private double bondLength;
    private double energy;
    private int bondStrength; // single, double, triple bond, etc.

    /**
     * Create a bond containing the two atoms provided as arguments.
     * @param a1 one of the atoms in the bond
     * @param a2 the other atom in the bond.
    */
    public Bond(Atom a1, Atom a2){
        this(a1, a2, 1);
    }

    /**
     * Create a bond containing the two atoms and with a bond strength
     * as provided as arguments.
     * @param a1 one of the atoms in the bond
     * @param a2 the other atom in the bond.
     * @param bondStrength 1 for single bond, 2 for double bond, etc...
    */
    public Bond(Atom a1, Atom a2, int bondStrength){
        // Use an ordering when assigning atoms into the bond.
        // The smaller atom goes into atomOne and the larger into atomTwo.
        if(a1.compareTo(a2) < 0){
            this.a1 = a1;
            this.a2 = a2;
        } else if(a1.compareTo(a2) > 0){
            this.a2 = a1;
            this.a1 = a2;
        } else {
            throw new IllegalArgumentException("Atoms added to a bond must be different.");
        }
        this.bondStrength = bondStrength;
        this.bondLength = a1.distance(a2);
        this.energy = calculateEnergy();
    }

    /**
     * Returns the length of the bond
     * @return the length of the bond in Angstroms.
    */
    public double getBondLength(){
        return this.bondLength;
    }

    /**
     * Returns the strength of the bond.
     * @return the strength of the bond, 1 for single, 2 for double, etc...
    */
    public int getBondStrength(){
        return this.bondStrength;
    }

    /**
     * Returns the potential energy of the bond. The potential energy is modeled
     * as a spring. The bond wants to be at a particular length, corresponding
     * to a minimum energy. Any deviation from that length increases the energy.
     * @return the potential energy of this bond.
    */
    public double getEnergy(){
        return this.energy;
    }

    /**
     * Returns one of the atoms in this bond.
     * @return one of the Atoms in this bond.
    */
    public Atom getAtomOne(){
        return this.a1;
    }

    /**
     * Returns the other Atom in this bond.
     * @return  the other Atom in this bond.
    */
    public Atom getAtomTwo(){
        return this.a2;
    }

    /**
     * Checks whether a given atom is in this bond.
     * @param atom the atom to check if in bond.
     * @return true if atom is in this bond
    */
    public boolean containsAtom(Atom atom){
        if(getAtomOne().equals(atom) || getAtomTwo().equals(atom)){
            return true;
        }
        return false;
    }

    /**
     * A private helper method to calculate the energy of this bond.
     * @return the energy of this bond in kcal/mol
    */
    private double calculateEnergy(){
        double bondForceConstant = EncadParameters.getBondForceConstant(a1, a2);
        double idealBondLength = EncadParameters.getBondLength(a1, a2);
        return (bondForceConstant * Math.pow(getBondLength() - idealBondLength, 2));
    }

    /**
     * The hashCode of a bond is the concatenation of the hash codes of the
     * atoms in the bond.
     * @return an int, the concatenations of the hash codes of the atoms in
     *         the bond.
    */
    @Override
    public int hashCode(){
        return Integer.parseInt("" + getAtomOne().hashCode() + getAtomTwo().hashCode());
    }

    /**
     * Overriden equals method. Returns true if the bonds contain the same
     * atoms. That is,
    */
    @Override
    public boolean equals(Object o){
        if(o instanceof Bond){
            Bond other = (Bond)o;
            if(this.getAtomOne().equals(other.getAtomOne())
              && this.getAtomTwo().equals(other.getAtomTwo())){
                return true;
            }
            if(this.getAtomOne().equals(other.getAtomTwo())
              && this.getAtomTwo().equals(other.getAtomOne())){
                return true;
            }
        }
        return false;
    }

}
