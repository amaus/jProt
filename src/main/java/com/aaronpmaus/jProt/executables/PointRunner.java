package com.aaronpmaus.jProt.executables;
import com.aaronpmaus.jProt.protein.Atom;
import com.aaronpmaus.jMath.linearAlgebra.Matrix;

public class PointRunner{
    public static void main(String[] args){
        Atom atom1 = new Atom("H");
        atom1.moveTo(3,4,5);
        System.out.println(atom1.toString());

        Atom atom2 = new Atom("S");
        atom2.moveTo(3,4,6);
        System.out.println(atom2.toString());

        System.out.println("distance: " + atom1.distance(atom2));
        double[][] mat = { {1,2,3}, {4,5,6}, {0,9,8} };
        Matrix A = new Matrix(mat);
        double[][] mat2 = { {1}, {1}, {1} };
        Matrix B = new Matrix(mat2);
        System.out.println(A);
        System.out.println();
        System.out.println(B);
        System.out.println();
        System.out.println(A.multiply(B));
        System.out.println();

        double[][] mat3 = { {2,3,5}, {7,11,13}, {17,19,23} };
        Matrix C = new Matrix(mat3);
        double[][] mat4 = { {2,3,5}, {7,11,13} };
        Matrix D = new Matrix(mat4);
        System.out.println(C);
        System.out.println();
        System.out.println(D.transpose());
        System.out.println();
        System.out.println(C.multiply(D.transpose()));
        
    }
}
