package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.*;
import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.linearAlgebra.*;
import java.util.Collection;
import java.util.HashSet;

/**
* This class represents a Molecule. It can be used for Amino Acids,
* PolypeptideChains, H20, anything that is a Molecule.
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class Molecule {
  private double mass;
  private Vector centerOfMass;
  private HashSet<Bond> covalentBonds;
  private HashSet<Atom> atoms;

  public Molecule(){
    this.atoms = new HashSet<Atom>();
    this.covalentBonds = new HashSet<Bond>();
    this.mass = 0.0;
    this.centerOfMass = new Vector(0.0, 0.0, 0.0);
  }
  /**
  * A constructor for a molecule. Takes a Collection of Bonds and builds
  * the Molecule from that.
  * @param covalentBonds the bonds in this Molecule
  */
  public Molecule(Collection<Bond> covalentBonds){
    this();
    for(Bond b: covalentBonds){
      this.covalentBonds.add(b);
    }
    for(Bond b : covalentBonds){
      atoms.add(b.getAtomOne());
      atoms.add(b.getAtomTwo());
    }
    calcCenterOfMass();
  }

  private void calcCenterOfMass(){
    double comX = 0.0;
    double comY = 0.0;
    double comZ = 0.0;
    this.mass = 0.0;
    for(Atom a: atoms){
      this.mass += a.getMass();
      comX += a.getMass() * a.getCoordinates().getValue(0).doubleValue();
      comY += a.getMass() * a.getCoordinates().getValue(1).doubleValue();
      comZ += a.getMass() * a.getCoordinates().getValue(2).doubleValue();
    }
    this.centerOfMass = new Vector(comX/mass, comY/mass, comZ/mass);
  }

  /**
  * A query to get the mass.
  * @return the mass of this molecule
  */
  public double getMass(){
    return this.mass;
  }

  /**
  * A query to get the center of mass of this molecule
  * @return a Vector holding the coordinates of the center of mass
  */
  public Vector getCenterOfMass(){
    return this.centerOfMass;
  }

  /**
  * Adds a bond to this molecule. It will add both atoms if they are not
  * already added. It will also add the bond if it hasn't been added.
  * @param bond the bond to add
  */
  protected void addBond(Bond bond){
    atoms.add(bond.getAtomOne());
    atoms.add(bond.getAtomTwo());
    this.covalentBonds.add(bond);
    calcCenterOfMass();
  }

  /**
  * @return a Collection of the covalent bonds in this molecule
  */
  public Collection<Bond> getBonds(){
    return new HashSet<Bond>(covalentBonds);
  }

  public boolean containsAtom(Atom atom){
    return this.atoms.contains(atom);
  }

  /**
  * Returns a String Representation of this Molecule.
  * @return a String representation of this molecule
  */
  @Override
  public String toString(){
    String str = "";
    for(Atom a : this.atoms){
      str += String.format(a.getAtomName() + " ");
    }
    str += "Center Of Mass: " + getCenterOfMass() + "\n";
    str += "Mass: " + getMass() + "\n";
    return str;
  }
}
