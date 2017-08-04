package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.linearAlgebra.*;

/**
 * The classical building block of all matter, the Atom.
 * <p>
 * An atom is defined as having a mass, charge, radius, and coordinates.
 * It also knows its element, which is the string representations of its
 * element name. Eg. "H", "C", "N", "O", "S". In fact, these are the only
 * elements allowed at the moment.
 * <p>
 * TODO: Add input file with atom parameters to allow for an atom of
 * any element to be instantiated.
 * @author Aaron Maus aaron@aaronpmaus.com
 * @version 0.6.0
 * @since 0.1.0
*/
public class Atom implements Comparable<Atom> {
    private static int maxSerialNum = 0;
    private final double mass;
    private double charge;
    private final String element;
    private final String atomName;
    private final double radius;
    private Vector itsCoordinates;
    private int serialNumber;
    private double occupancy; // -1.0 means no value
    private double tempFactor; // -1.0 means no value

    /**
     * Builds an atom of type element located at the origin.
     * <p>
     * The charge is initialized to 0. The mass and radius are assigned
     * according to the element type. The radius is determined according
     * to the encad protein energy function parameters.
     * @param atomName   The String representation of the atom. One
     * of "CA", "CB", "CG", "O", "C", "NH1", etc...
     * @param serialNumber The Serial Number for this Atom.
     * @param occupancy the Occupancy for this Atom.
     * @param tempFactor the Temperature Factor for this Atom.
     * @param charge the charge of the atom.
     * @throws IllegalArgumentException thrown if the element is not one of
     *                                  the allowed elements.
     * @since 0.6.0
    */
    public Atom(String atomName, int serialNumber, double occupancy,
                double tempFactor, double charge){
        this(atomName, serialNumber, occupancy, tempFactor, charge, "0.0", "0.0", "0.0");
    }

    /**
     * Builds an atom of type element located at the coordinates specified.
     * <p>
     * The charge is initialized to 0. The mass and radius are assigned
     * according to the element type. The radius is determined according
     * to the encad protein energy function parameters.
     * @param atomName   The String representation of the atom. One
     * of "CA", "CB", "CG", "O", "C", "NH1", etc...
     * @param serialNumber The Serial Number for this Atom.
     * @param occupancy the Occupancy for this Atom.
     * @param tempFactor the Temperature Factor for this Atom.
     * @param charge the charge of the atom.
     * @param x         The x coordinate of this atom.
     * @param y         The y coordinate of this atom.
     * @param z         The z coordinate of this atom.
     * @throws IllegalArgumentException thrown if the element is not one of
     *                                  the allowed elements.
     * @since 0.6.0
    */
    public Atom(String atomName, int serialNumber, double occupancy,
                double tempFactor, double charge, String x, String y, String z){
        itsCoordinates = new Vector(x,y,z);
        this.atomName = atomName.toUpperCase();
        this.charge = 0;
        this.serialNumber = serialNumber;
        this.occupancy = occupancy;
        this.tempFactor = tempFactor;

        // set the mass and radius depending on the element
        if(this.atomName.charAt(0) == 'H'){
            this.mass = 1.008;
            this.radius = 1.10;
            this.element = "H";
        } else if(this.atomName.charAt(0) == 'C'){
            this.mass = 12.011;
            this.radius = 1.85;
            this.element = "C";
        } else if(this.atomName.charAt(0) == 'N'){
            this.mass = 14.007;
            this.radius = 1.65;
            this.element = "N";
        } else if(this.atomName.charAt(0) == 'O'){
            this.mass = 15.999;
            this.radius = 1.60;
            this.element = "O";
        } else if(this.atomName.charAt(0) == 'S'){
            this.mass = 32.064;
            this.radius = 1.85;
            this.element = "S";
        } else { // must be one of the above. otherwise, throw exception
            throw new IllegalArgumentException("Atom() invalid element type. Must be H, C, N, O, or S.");
        }
    }

    /**
     * Returns the coordinates of this atom.
     * @return  a reference to this atom's Vector.
     * @since 0.1.0
    */
    public Vector getCoordinates(){
        return itsCoordinates;
    }

    /**
     * Returns the serial number for this atom. Every atom is required
     * to have a unique serial number.
     * @return the serial number for this atom.
     * @since 0.6.0
    */
    public int getSerialNumber(){
        return this.serialNumber;
    }

    /**
     * Returns the occupancy for this atom. The Occupancy is a percent, and
     * it represents the number of models from crystallography that had
     * the atom at these coordinates. If the atom has been placed at these
     * coordinates by simulation, the value will be -1.
     * @return the occupancy of this atom
     * @since 0.6.0
    */
    public double getOccupancy(){
        return this.occupancy;
    }

    /**
     * Returns the temperature factor for this atom. The Temperature Factor
     * is the degree to which the electron density of this atom was spread
     * out in the x-ray scattering data. It can represent the true static
     * or dynamic mobility of the atom. It also can indicate where there are
     * errors in the model building and can be used to predict protein structure
     * disorder.
     * @return the Temperature Factor of this Atom.
     * @since 0.6.0
    */
    public double getTempFactor(){
        return this.tempFactor;
    }

    /**
     * Returns the charge of the atom. This treats atoms as point charges which
     * is a simplifying assumption. A charge on an atom is not necessarily
     * uniform.
     * @return the charge of the atom.
     * @since 0.6.0
    */
    public double getCharge(){
        return this.charge;
    }

    /**
     * Returns the mass of the atom.
     * @return the mass of the atom
     * @since 0.6.0
    */
    public double getMass(){
        return this.mass;
    }

    /**
     * Returns the chemical abbriviation for this element, C for carbon,
     * H for hydrogren, Fe for iron, and so on.
     * @return returns the chemical abbriviation for this element.
     * @since 0.6.0
    */
    public String getElement(){
        return this.element;
    }

    /**
     * Returns the name of this atom. eg. CA, CB, C, O, N, CD, CD1, etc...
     * @return the name of this atom.
     * @since 0.6.0
    */
    public String getAtomName(){
        return this.atomName;
    }

    /**
     * Calculates the distance between two atoms.
     * @param otherAtom     The other atom to calculate the distance to.
     * @return The euclidean distance between these two atoms.
     * @since 0.1.0
    */
    public double distance(Atom otherAtom) {
        return getCoordinates().distance(otherAtom.getCoordinates());
    }

    @Override
    /**
     * Compares the two atoms based on their serialNumber
     * @return a negative number if this.getSerialNumber() < other.getSerialNumber()
     *      0 if this.getSerialNumber() == other.getSerialNumber()
     *      a positive number if this.getSerialNumber() > other.getSerialNumber()
     * @since 0.6.0
    */
    public int compareTo(Atom other){
        return this.getSerialNumber() - other.getSerialNumber();
    }

    @Override
    /**
     * Returns a String representation of this atom.
     * @return a string containing the atom's element, mass, charge,
     * and coordinates.
     * @since 0.1.0
    */
    public String toString(){
        String str = String.format("Atom Name: %s\n",getAtomName()) +
                    String.format("Serial Number: %d\n", getSerialNumber()) +
                    String.format("Mass: %.2f\n",getMass()) +
                    String.format("Charge: %.2f\n", getCharge())+
                    String.format("Coords: %s\n", getCoordinates());
                    //String.format("Occupancy: %.2f\n", getOccupancy()) +
                    //String.format("Temp Factor: %.2f\n", getTempFactor());
        //str += "Coords: " + getCoordinates().toString();
        return str;
    }

    @Override
    /**
     * Returns a hash code value for this Atom.
     * The hash code is the atom's serial number
     * @return the hash code
     * @since 0.6.0
    */
    public int hashCode(){
        return getSerialNumber();
    }

    @Override
    /**
     * Idicates whether two atoms are equal. Equality is based on the serial
     * number. Every atom must have a different serial number
     * @return returns true if the two atoms have the same serial number
     * @since 0.6.0
    */
    public boolean equals(Object obj){
        if(obj instanceof Atom){
            Atom other = (Atom)obj;
            if(other.getSerialNumber() == this.getSerialNumber()){
                return true;
            }
        }
        return false;
    }

} // end class Atom
