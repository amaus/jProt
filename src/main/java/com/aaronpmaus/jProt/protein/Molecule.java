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
    private UndirectedGraph<Atom> atoms;
    private double mass;
    private Vector centerOfMass;
    private HashSet<Bond> covalentBonds;

    public Molecule(){
        this.atoms = new UndirectedGraph<Atom>();
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
        for(Bond b: covalentBonds){
            this.covalentBonds.add(b);
        }
        atoms = new UndirectedGraph<Atom>();
        for(Bond b : covalentBonds){
            atoms.addEdge(b.getAtomOne(), b.getAtomTwo());
        }
        calcCenterOfMass();
    }

    private void calcCenterOfMass(){
        Collection<Node<Atom>> atomsList = this.atoms.getNodes();
        double comX = 0.0;
        double comY = 0.0;
        double comZ = 0.0;
        this.mass = 0.0;
        for(Node<Atom> n: atomsList){
            Atom a = n.get();
            this.mass += a.getMass();
            comX += a.getMass()*a.getCoordinates().getValue(0);
            comY += a.getMass()*a.getCoordinates().getValue(1);
            comZ += a.getMass()*a.getCoordinates().getValue(2);
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
     * Adds an atom to this molecule. This will add this and all the atoms
     * bonded to it to this molecule if they are not already there. Will also
     * add the bonds if they haven't been added as well.
     * @param atomToAdd the atom to add to this molecule
     * @param atomsBondedToIt the atoms bonded to it.
    */
    public void addAtom(Atom atomToAdd, Atom... atomsBondedToIt){
        for(Atom atomBondedToIt : atomsBondedToIt){
            atoms.addEdge(atomToAdd, atomBondedToIt);
            Bond bond = new Bond(atomToAdd, atomBondedToIt);
            this.covalentBonds.add(bond);
        }
        calcCenterOfMass();
    }

    @Override
    /**
     * Returns a String Representation of this Molecule.
     * @return a String representation of this molecule
    */
    public String toString(){
        String str = "";
        for(Node<Atom> n : this.atoms.getNodes()){
            Atom a = n.get();
            str += String.format(a.getAtomName() + " ");
            str += String.format(n.toString() + "\n");
        }
        str += "Center Of Mass: " + getCenterOfMass() + "\n";
        str += "Mass: " + getMass() + "\n";
        return str;
    }
}
