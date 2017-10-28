package com.aaronpmaus.jProt.sequence;

import java.util.Iterator;
import java.util.LinkedList;

public class Sequence implements Iterable<String>{
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
  public Iterator<String> iterator(){
    LinkedList<String> list = new LinkedList<String>();
    for(Character character : seq.toCharArray()){
      list.add(String.format("%s",character));
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
