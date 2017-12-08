package com.aaronpmaus.jProt.sequence;

import java.util.Iterator;
import java.util.LinkedList;

public class Sequence implements Iterable<Character>{
  String seq;

  public Sequence(String seq){
    this.seq = seq;
  }

  /**
  * @return the String representation of this sequence
  */
  public String getSequenceString(){
    return this.seq;
  }

  @Override
  public Iterator<Character> iterator(){
    LinkedList<Character> list = new LinkedList<Character>();
    for(Character c :  seq.toCharArray()){
      list.add(c);
    }
    return list.iterator();
  }

  @Override
  public String toString(){
    return this.seq;
  }

  @Override
  public int hashCode(){
    return toString().hashCode();
  }

  @Override
  public boolean equals(Object obj){
    if(obj instanceof Sequence){
      Sequence other = (Sequence)obj;
      return this.toString().equals(other.toString());
    }
    return false;
  }
}
