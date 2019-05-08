import java.util.*;

public class Consensus{

    private String contig;
    private int pos;
    private char ref;
    private ArrayList<Character> chars;
    private char consensus;
    private int mutType; //default -5: unset
    

    /* 
     * tokens are from .consensus file where first three columns are headers(ContigName, position, reference base)
     * the very last column is the consensus of all lines at that specific position.
     */
    public Consensus(String[] tokens){
	this.contig = tokens[0];
	this.pos = Integer.parseInt(tokens[1]);
	this.ref = tokens[2].charAt(0);
	this.chars = new ArrayList<Character>();
	for(int i=3;i<(tokens.length-1);i++){
	    this.chars.add(new Character(tokens[i].charAt(0)));
	}
	this.consensus = tokens[tokens.length-1].charAt(0);
	this.mutType = -5;
    }
    
    public char getConsensus(){
	return this.consensus;
    }

    public String toString(){
	StringBuffer buffer =  new StringBuffer(this.contig + " " + this.pos + " " + this.ref);
	for(int i=0;i<chars.size();i++){
	    buffer.append(" " + chars.get(i).charValue());
	}
	buffer.append(" " + consensus);
	return buffer.toString();
    }

    public String getHeader(){
	return this.contig + "\t" + this.pos;
    }

    /*
     * -2 : AT > GC
     * -1 : GC > AT
     *  1 : AT > TA 
     *  2 : GC > TA
     *  3 : AT > CG
     *  4 : GC > CG
     *
     *  0 : no mutation or error
     *
     * EDIT: sets the mutType field with proper mutType code returned from determineMutationtype();
     *       then, it returns the index number of line(sample) in which the mutation has been observed.
     * return -1 if no mutation is found.
     */
    public int mutationType(){
	for(int i=0; i<chars.size();i++){
	    if(chars.get(i).charValue() != '-' && chars.get(i).charValue() != this.consensus){
		this.mutType = determineMutationType(chars.get(i).charValue());
		return i;
	    }
	}
	this.mutType = 0;
	return -1;
    }

    public char getMutBase(){
	for(int i=0; i<chars.size();i++){
	    if(chars.get(i).charValue() != this.consensus){
		return chars.get(i).charValue();
	    }
	}
	return this.consensus;
    }

    public int getMutType(){
	return this.mutType;
    }

    public int determineMutationType(char c){
	if(this.consensus == 'A'){
	    if(c == 'C')
		return 3;
	    else if(c == 'G')
		return -2;
	    else if(c == 'T')
		return 1;
	    else
		return 0;
	}else if(this.consensus == 'C'){
	    if(c == 'A')
		return 2;
	    else if(c == 'G')
		return 4;
	    else if(c == 'T')
		return -1;
	    else
		return 0;
	}else if(this.consensus == 'G'){
	    if(c == 'A')
		return -1;
	    else if(c == 'C')
		return 4;
	    else if(c == 'T')
		return 2;
	    else
		return 0;
	}else if(this.consensus == 'T'){
	    if(c == 'A')
		return 1;
	    else if(c == 'C')
		return -2;
	    else if(c == 'G')
		return 3;
	    else
		return 0;
	}else
	    return 0;
	    
    }
}
