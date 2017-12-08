package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.linearAlgebra.Vector3D;
import com.aaronpmaus.jMath.transformations.Transformable;
import com.aaronpmaus.jMath.transformations.Transformation;

/**
* The classical building block of all matter, the Atom.
*
* An atom is defined as having a mass, charge, radius, and coordinates.
* It also knows its element, which is the string representations of its
* element name. Eg. "H", "C", "N", "O", "S".
*
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.1.0
*/
public class Atom implements Comparable<Atom>, Transformable{
  // maxSerialNum is used for constructing Atoms when they haven't been read in from
  // a PDB file. Every atom is given (maxSerialNum+1)
  private static int maxSerialNum = 0;
  private final double mass;
  private double charge;
  private final String element;
  private final String atomName;
  private final double radius;
  private Vector3D itsCoordinates;
  private int serialNumber;
  private double occupancy; // -1.0 means no value
  private double tempFactor; // -1.0 means no value

  /**
  * Build an atom of type element at the origin.
  *
  * The charge is initialized to 0. The mass and radius are assigned
  * according to the element type. The radius is determined according
  * to the encad protein energy function parameters.
  *
  * @param atomName   The String representation of the atom. One
  * of "CA", "CB", "CG", "O", "C", "NH1", etc...
  * @param serialNumber The Serial Number for this Atom.
  * @param occupancy the Occupancy for this Atom.
  * @param tempFactor the Temperature Factor for this Atom.
  * @param charge the charge of the atom.
  * @since 0.6.1
  */
  public Atom(String atomName, int serialNumber, double occupancy,
  double tempFactor, double charge){
    this(atomName, serialNumber, occupancy, tempFactor, charge, 0.0, 0.0, 0.0);
  }

  /**
  * Builds an atom of type element located at the coordinates specified.
  *
  * The charge is initialized to 0. The mass and radius are assigned
  * according to the element type. The radius is determined according
  * to the encad protein energy function parameters.
  *
  * @param atomName   The String representation of the atom. One
  * of "CA", "CB", "CG", "O", "C", "NH1", etc...
  * @param serialNumber The Serial Number for this Atom.
  * @param occupancy the Occupancy for this Atom.
  * @param tempFactor the Temperature Factor for this Atom.
  * @param charge the charge of the atom.
  * @param x         The x coordinate of this atom.
  * @param y         The y coordinate of this atom.
  * @param z         The z coordinate of this atom.
  * @since 0.6.1
  */
  public Atom(String atomName, int serialNumber, double occupancy,
  double tempFactor, double charge, double x, double y, double z){
    this.itsCoordinates = new Vector3D(x,y,z);
    this.atomName = atomName.toUpperCase().trim();
    this.charge = charge;
    // If the serialNumber is -1, then this atom was built from one of the default amino acids.
    // Assign the serianNumber to be be the max serial number + 1
    if(serialNumber == -1){
      this.serialNumber = Atom.maxSerialNum+1;
    } else {
      this.serialNumber = serialNumber;
    }
    Atom.maxSerialNum = Math.max(Atom.maxSerialNum, this.serialNumber);
    this.occupancy = occupancy;
    this.tempFactor = tempFactor;

    // set the mass and radius depending on the element
    switch(this.atomName.charAt(0)){
      case 'C':
        this.mass = 12.011;
        this.radius = 1.85;
        this.element = "C";
        break;
      case 'N':
        this.mass = 14.007;
        this.radius = 1.65;
        this.element = "N";
        break;
      case 'O':
        this.mass = 15.999;
        this.radius = 1.60;
        this.element = "O";
        break;
      case 'S':
        this.mass = 32.064;
        this.radius = 1.85;
        this.element = "S";
        break;
      default:
        this.mass = 1.008;
        this.radius = 1.10;
        this.element = "H";
    }
  }

  /**
  * Return the coordinates of this atom.
  *
  * @return a Vector containing the coordinates of this atom
  * @since 0.1.0
  */
  public Vector3D getCoordinates(){
    return itsCoordinates;
  }

  /**
  * Return the serial number for this atom. Every atom is required
  * to have a unique serial number.
  *
  * @return the serial number for this atom.
  * @since 0.6.0
  */
  public int getSerialNumber(){
    return this.serialNumber;
  }

  /**
  * Return the occupancy for this atom. The Occupancy is a percent, and
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
  * Return the temperature factor for this atom. The Temperature Factor
  * is the degree to which the electron density of this atom was spread
  * out in the x-ray scattering data. It can represent the true static
  * or dynamic mobility of the atom. It also can indicate where there are
  * errors in the model building and can be used to predict protein structure
  * disorder.
  *
  * @return the Temperature Factor of this Atom.
  * @since 0.6.0
  */
  public double getTempFactor(){
    return this.tempFactor;
  }

  /**
  * Return the charge of the atom. This treats atoms as point charges which
  * is a simplifying assumption. A charge on an atom is not necessarily
  * uniform.
  *
  * @return the charge of the atom.
  * @since 0.6.0
  */
  public double getCharge(){
    return this.charge;
  }

  /**
  * Return the mass of the atom.
  *
  * @return the mass of the atom
  * @since 0.6.0
  */
  public double getMass(){
    return this.mass;
  }

  /**
  * Return the chemical abbriviation for this element, C for carbon,
  * H for hydrogren, Fe for iron, and so on.
  *
  * @return returns the chemical abbriviation for this element.
  * @since 0.6.0
  */
  public String getElement(){
    return this.element;
  }

  /**
  * Return the name of this atom. eg. CA, CB, C, O, N, CD, CD1, etc...
  *
  * This name is the name of the Atom as listed in PDB files.
  *
  * @return the name of this atom.
  * @since 0.6.0
  */
  public String getAtomName(){
    return this.atomName;
  }

  /**
  * Calculate and return the distance between two atoms.
  *
  * @param otherAtom The other atom to calculate the distance to.
  * @return The euclidean distance between these two atoms.
  * @since 0.1.0
  */
  public double distance(Atom otherAtom) {
    return getCoordinates().distance(otherAtom.getCoordinates());
  }

  @Override
  public void applyTransformation(Transformation t){
    this.itsCoordinates.applyTransformation(t);
  }

  /**
  * Compare the two atoms based on their serialNumber
  *
  * @return a negative number if {@code this.getSerialNumber() < other.getSerialNumber()}
  *      0 if this.getSerialNumber() == other.getSerialNumber()
  *      a positive number if {@code this.getSerialNumber() > other.getSerialNumber()}
  * @since 0.6.0
  */
  @Override
  public int compareTo(Atom other){
    return this.getSerialNumber() - other.getSerialNumber();
  }

  /**
  * Return a String representation of this atom.
  *
  * @return a string containing the atom's element, mass, charge,
  * and coordinates.
  * @since 0.1.0
  */
  @Override
  public String toString(){
    String str = String.format("Atom Name: %s\n",getAtomName()) +
    String.format("Element: %s\n", getElement()) +
    String.format("Serial Number: %d\n", getSerialNumber()) +
    String.format("Mass: %.2f\n",getMass()) +
    String.format("Charge: %.2f\n", getCharge())+
    String.format("Coords: %s\n", getCoordinates());
    //String.format("Occupancy: %.2f\n", getOccupancy()) +
    //String.format("Temp Factor: %.2f\n", getTempFactor());
    //str += "Coords: " + getCoordinates().toString();
    return str;
  }

  /**
  * Return a hash code value for this Atom.
  *
  * The hash code is the atom's serial number
  *
  * @return the hash code
  * @since 0.6.0
  */
  @Override
  public int hashCode(){
    return getSerialNumber();
  }

  /**
  *
  * Equality is based on the serial number. Every atom must have a different serial number
  *
  * @return returns true if the two atoms have the same serial number
  * @since 0.6.0
  */
  @Override
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
