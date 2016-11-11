package com.aaronpmaus.jProtein.protein;

import java.lang.IllegalArgumentException;
import com.aaronpmaus.jMath.*;

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
 * @version 1.0.0
 * @since 1.0.0
*/
public class Atom{
    private final double mass;
    private double charge;
    private final String element;
    private final double radius;
    private CartesianCoordinates itsCoordinates;

    /**
     * Builds an atom of type element located at the origin.
     * <p>
     * The charge is initialized to 0. The mass and radius are assigned
     * according to the element type. The radius is determined according
     * to the encad protein energy function parameters.
     * @param element   The String representation of the atom element. Must
     * be one of "H", "C", "N", "O", or "S".
     * @throws IllegalArgumentException thrown if the element is not one of
     *                                  the allowed elements.
    */
    public Atom(String element){
        this(element, 0, 0, 0);
    }

    /**
     * Builds an atom of type element located at the coordinates specified.
     * <p>
     * The charge is initialized to 0. The mass and radius are assigned
     * according to the element type. The radius is determined according
     * to the encad protein energy function parameters.
     * @param element   The String representation of the atom element. Must
     * be one of "H", "C", "N", "O", or "S".
     * @param x         The x coordinate of this atom.
     * @param y         The y coordinate of this atom.
     * @param z         The z coordinate of this atom.
     * @throws IllegalArgumentException thrown if the element is not one of
     *                                  the allowed elements.
    */
    public Atom(String element, double x, double y, double z){
        itsCoordinates = new CartesianCoordinates(x,y,z);
        this.element = element;
        this.charge = 0;

        //UndirectedGraph<String, DefaultEdge> stringGraph = new SimpleGraph<String, DefaultEdge>(DefaultEdge.class);

        // set the mass and radius depending on the element
        if(this.element == "H"){
            this.mass = 1.008;
            this.radius = 1.10;
        } else if(this.element == "C"){
            this.mass = 12.011;
            this.radius = 1.85;
        } else if(this.element == "N"){
            this.mass = 14.007;
            this.radius = 1.65;
        } else if(this.element == "O"){
            this.mass = 15.999;
            this.radius = 1.60;
        } else if(this.element == "S"){
            this.mass = 32.064;
            this.radius = 1.85;
        } else { // must be one of the above. otherwise, throw exception
            throw new IllegalArgumentException("Atom() invalid element type. Must be H, C, N, O, or S.");
        }
    }

    /**
     * Returns the coordinates of this atom.
     * @return  a reference to this atom's CartesianCoordinates.
    */
    public CartesianCoordinates getCoordinates(){
        return itsCoordinates;
    }

    /**
     * Calculates the distance between two atoms.
     * @param otherAtom     The other atom to calculate the distance to.
     * @return The euclidean distance between these two atoms.
    */
    public double distance(Atom otherAtom) {
        return getCoordinates().distance(otherAtom.getCoordinates());
    }

    /**
     * Changes the location of this atoms coordinates to be those specified.
     * @param x    The x coordinate to move this atom to.
     * @param y    The y coordinate to move this atom to.
     * @param z    The z coordinate to move this atom to.
    */
    public void moveTo(double x, double y, double z) {
        getCoordinates().moveTo(x, y, z);
    }

    @Override
    /**
     * Returns a String representation of this atom.
     * <p>
     * Contains the atom's element, mass, charge, and coordinates.
    */
    public String toString(){
        String str = this.element +
                    "\nMass: " + this.mass + 
                    "\nCharge: " + this.charge;
        str += "\nCoords: " + getCoordinates().toString();
        return str;
    }

} // end class Atom
