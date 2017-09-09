package com.aaronpmaus.jProt.io;

public class AtomRecord{
  private int serial;
  private String atomName;
  private String altLoc;
  private String resName;
  private String chainID;
  private int resSeq;
  private String iCode;
  private String x;
  private String y;
  private String z;
  private double occupancy;
  private double tempFactor;
  private String element;
  private double charge;

  public AtomRecord(int serial, String atomName, String altLoc,
  String resName, String chainID, int resSeq,
  String iCode, String x, String y, String z,
  double occupancy, double tempFactor,
  String element, double charge){
    this.serial = serial;
    this.atomName = atomName;
    this.altLoc = altLoc;
    this.resName = resName;
    this.chainID = chainID;
    this.resSeq = resSeq;
    this.iCode = iCode;
    this.x = x;
    this.y = y;
    this.z = z;
    this.occupancy = occupancy;
    this.tempFactor = tempFactor;
    this.element = element;
    this.charge = charge;
  }
  public int getSerial(){ return this.serial; }
  public String getAtomName(){ return this.atomName; }
  public String getAltLoc(){ return this.altLoc; }
  public String getResName(){ return this.resName; }
  public String getChainID(){ return this.chainID; }
  public int getResSeq(){ return this.resSeq; }
  public String getICode(){ return this.iCode; }
  public String getX(){ return this.x; }
  public String getY(){ return this.y; }
  public String getZ(){ return this.z; }
  public double getOccupancy(){ return this.occupancy; }
  public double getTempFactor(){ return this.tempFactor; }
  public String getElement(){ return this.element; }
  public double getCharge(){ return this.charge; }

}
