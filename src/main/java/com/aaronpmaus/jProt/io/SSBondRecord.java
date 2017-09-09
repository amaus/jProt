package com.aaronpmaus.jProt.io;

public class SSBondRecord{
  private String chainID1;
  private int resID1;
  private String chainID2;
  private int resID2;

  public SSBondRecord(String chainID1, int resID1, String chainID2, int resID2){
    this.chainID1 = chainID1;
    this.resID1 = resID1;
    this.chainID2 = chainID2;
    this.resID2 = resID2;
  }
  public String getChainID1(){ return this.chainID1; }
  public String getChainID2(){ return this.chainID2; }
  public int getResID1(){ return this.resID1; }
  public int getResID2(){ return this.resID2; }

}
