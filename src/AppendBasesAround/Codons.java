import java.io.*;
import java.util.*;

public class Codons{

    private Hashtable<String, Character> codon2amino;

    /*
     * codonfile is a tab-delimited file where 1st column is codon triplet in UPPER characters
     * and 2nd column is a amino acid symbol in UPPER character. 
     * symbol, * is used to denote a STOP codon
     * 
     */
    public Codons(String codonfile){
	this.codon2amino = new Hashtable<String, Character>();
	this.loadCodon2amino(codonfile);
    }

    private void loadCodon2amino(String codonfile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(codonfile));
	    String curLine = "";
	    while((curLine=br.readLine())!=null){
		String[] tokens = curLine.split("\\t");
		this.codon2amino.put(tokens[0], tokens[1].charAt(0));
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    
    public char toAmino(String triplet){
	
	Character ch =  this.codon2amino.get(triplet);
	if(ch == null){
	    System.err.println(triplet);
	    return '!';
	    //System.exit(1);
	}
	    
	return ch.charValue();
	    
	/*if(ch !=null)
	    return ch.charValue();
	else
	return '@';*/
    }


    

    /* RETURNS 1 if Syn, 0 if non-syn, -1, if NG due to bases other than atgc */
    /*public int isSyn(String t1, String t2){
	//System.out.println(t1 + "\t" + t2);
	char c1 = this.toAmino(t1);
	char c2 = this.toAmino(t2);
	if(c1 == '!' || c2 == '!')
	    return -1;
	if(c1 == c2)
	    return 1;
	else
	    return 0;
	    }*/

    public MutType isSyn(String t1, String t2){
	//System.out.println(t1 + "\t" + t2);
	char c1 = this.toAmino(t1);
	char c2 = this.toAmino(t2);
	if(c1 == '!' || c2 == '!')
	    return new MutType(-1,c1,c2,t1,t2);
	if(c1 == c2)
	    return new MutType(1,c1,c2,t1,t2);
	else
	    return new MutType(0,c1,c2,t1,t2);
    }
}
