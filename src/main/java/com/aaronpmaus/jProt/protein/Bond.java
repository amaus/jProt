package com.aaronpmaus.jProt.protein;
import java.lang.IllegalArgumentException;

/**
 * A bond is the covalent bond between two atoms.
 * @version 0.6.0
 * @since 0.6.0
*/
public class Bond {
    private Atom a1;
    private Atom a2;
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
    }

    /**
     * Returns the length of the bond
     * @return the length of the bond in Angstroms.
    */
    public double getBondLength(){
        return this.getAtomOne().distance(this.getAtomTwo());
    }

    /**
     * Returns the strength of the bond.
     * @return the strength of the bond, 1 for single, 2 for double, etc...
    */
    public int getBondStrength(){
        return this.bondStrength;
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

    public boolean containsHydrogen(){
      if(getAtomOne().getElement().equals("H") || getAtomTwo().getElement().equals("H")){
        return true;
      }
      return false;
    }

    @Override
    public int hashCode(){
      long hash = Long.parseLong("" + getAtomOne().hashCode() + getAtomTwo().hashCode());
      hash %= Integer.MAX_VALUE;
      return (int) hash;
    }

    /**
    * Overriden equals method.
    * @return true if the bonds contain the same atoms.
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

    @Override
    public String toString(){
      return getAtomOne().getAtomName() + " - " + getAtomTwo().getAtomName();
    }
  }
