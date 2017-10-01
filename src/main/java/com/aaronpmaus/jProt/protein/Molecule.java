package com.aaronpmaus.jProt.protein;

import com.aaronpmaus.jMath.*;
import com.aaronpmaus.jMath.graph.*;
import com.aaronpmaus.jMath.linearAlgebra.*;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;

/**
* A Molecule is a set of atoms all connected by covalent bonds.
* @author Aaron Maus aaron@aaronpmaus.com
* @version 0.6.0
* @since 0.6.0
*/
public class Molecule {
  private HashSet<Bond> covalentBonds;
  private HashSet<Atom> atoms;

  public Molecule(){
    this.atoms = new HashSet<Atom>();
    this.covalentBonds = new HashSet<Bond>();
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
  }

  /**
  * A query to get the mass.
  * @return the mass of this molecule
  */
  public double getMass(){
    double mass = 0.0;
    for(Atom a: atoms){
      mass += a.getMass();
    }
    return mass;
  }

  /**
  * Calculate and return the center of mass of this molecule
  * @return a Vector holding the coordinates of the center of mass
  */
  public Vector getCenterOfMass(){
    double comX = 0.0;
    double comY = 0.0;
    double comZ = 0.0;
    double mass = 0.0;
    for(Atom a: atoms){
      mass += a.getMass();
      comX += a.getMass() * a.getCoordinates().getValue(0).doubleValue();
      comY += a.getMass() * a.getCoordinates().getValue(1).doubleValue();
      comZ += a.getMass() * a.getCoordinates().getValue(2).doubleValue();
    }
    return new Vector(comX/mass, comY/mass, comZ/mass);
  }

  /**
  * Adds a bond to this molecule. It will add both atoms if they are not
  * already added. It will also add the bond if it hasn't been added.
  * @param bond the bond to add
  */
  public void addBond(Bond bond){
    atoms.add(bond.getAtomOne());
    atoms.add(bond.getAtomTwo());
    this.covalentBonds.add(bond);
  }

  /**
  * Remove the bond from this molecule if it exists.
  * @param bond a bond to remove from this molecule
  */
  public void removeBond(Bond bond){
    this.covalentBonds.remove(bond);
  }

  /**
  * Remove the atom from this molecule if it exists.
  * @param atom an atom to remove from this molecule
  */
  public void removeAtom(Atom atom){
    this.atoms.remove(atom);
  }

  /**
  * @return a Collection of the covalent bonds in this molecule
  */
  public Collection<Bond> getBonds(){
    return new HashSet<Bond>(covalentBonds);
  }

  /**
  * @return the number of covalent bonds in this molecule
  */
  public int getNumBonds(){
    return this.covalentBonds.size();
  }

  public boolean containsAtom(Atom atom){
    return this.atoms.contains(atom);
  }

  /**
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
