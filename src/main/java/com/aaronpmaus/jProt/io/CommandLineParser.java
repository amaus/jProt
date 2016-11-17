package com.aaronpmaus.jProt.io;

import java.util.ArrayList;
import java.security.InvalidParameterException;

/**
 * Provides a set of methods to parse command line arguments for
 * executables. Arguments to set values are in the form flag:val
 * where flag is something like -d or --distance, and val is the
 * value to be used.
*/
public class CommandLineParser{
    private ArrayList<String> args;

    /**
     * Constructor requires the array of command line arguments
     * @param arguments the array of command line args
    */
    public CommandLineParser(String[] arguments){
        this.args = new ArrayList<String>( );
        for(String arg : arguments){
            this.args.add(arg);
        }
    }

    /**
     * Checks to see if the command line arguments contains the given
     * flag. A command line argument can consist of a flag and a value
     * separated by a ':'. eq. -d:2.6 or --distance:2.6.
     * This method checks if the flag of any of the arguments is
     * equal to the given flag.
     * @param flag the flag to check for
     * @return true if found, false otherwise.
    */
    public boolean contains(String flag){
        for(String arg : args){ // for every cmd argument
            String[] pair = arg.split(":");
            if(pair[0].equals(flag)){ // if it contains the flag
                return true; // return true
            }
        }
        // if none of the args contained the flag, return false
        return false;
    }

    /*
     * private helper method to get the index of the flag in the parameter list
    */
    private int getIndex(String flag){
        boolean found = false;
        int index = -1;
        for(int i = 0; i < args.size(); i++){
            String[] pair = args.get(i).split(":");
            if(pair[0].equals(flag)){
                found = true;
                index = i;
            }
        }
        if(!found){
            throw new InvalidParameterException("Flag " + flag + 
                        "not found in command line arguments.");
        }
        return index;
    }

    /**
     * returns the value for the given flag.
     * @param flag the flag to get the value of
     * @return the value for that flag
     * @throws InvalidParameterException if the flag is not in the cmd line args
    */
    public String getValue(String flag) throws InvalidParameterException {
        return args.get( getIndex(flag)+1 );
    }
/*
    public String getValue(String flag) throws InvalidParameterException {
        boolean found = false;
        String val = "";
        for(String arg : args){ // for every cmd argument
            char lastCharacter = arg.charAt(arg.length()-1);
            String[] tokens = arg.split(":");
            if(tokens[0].equals(flag)){ // if it contains the flag
                // just in case any values have a ':' in them,
                // concatenate all tokens after the flag back together
                for(int i = 1; i < tokens.length; i++){
                    val += tokens[i];
                    if(i != tokens.length -1){
                        val+=":";
                    }
                }
                if(lastCharacter == ':'){
                    val+=":";
                }
                found = true;
            }
        }
        if(!found){
            throw new InvalidParameterException("Flag " + flag + 
                        "not found in command line arguments.");
        }
        return val;
    }
*/
}
