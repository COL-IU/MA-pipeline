import java.util.*;

class Mutation{
    private String contig;
    private int pos;
    private char ref;
    private char mut;

    public Mutation(String ctg, int pos, char ref, char mut){
	this.contig = ctg;
	this.pos = pos;
	this.ref = ref;
	this.mut = mut;
    }

    public String getContig(){
	return this.contig;
    }

    public int getPos(){
	return this.pos;
    }
    
    public char getRef(){
	return this.ref;
    }
    
    public char getMut(){
	return this.mut;
    }
}

public class Consensus{

    private String contig;
    private int pos;
    private char ref;
    private ArrayList<Character> chars;
    private char consensus;
    private int mutType; //default -5: unset
    private int firstMutIndex;

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
	this.firstMutIndex = -1;
    }

    public Consensus(String ctg, int pos, char ref, char consensus, int lineIndex, char baseForThatLine, int numLines){
	this.contig = ctg;
	this.pos = pos;
	this.ref = ref;
	this.consensus = consensus;
	this.chars =new ArrayList<Character>();
	for(int i=0;i<numLines;i++){
	    if( (i+1) != lineIndex)
		this.chars.add(new Character('-'));
	    else
		this.chars.add(baseForThatLine);
	}
	this.mutType = -5;
	this.firstMutIndex = -1;
    }

    public boolean charsHasAtLeast2Consensus(){
	int counter = 0;
	for(int i=0; i<chars.size(); i++){
	    if(chars.get(i).charValue() == consensus)
		counter++;
	}
	if(counter > 1)
	    return true;
	return false;
    }

    public boolean isPolymorphism(){
	if(isACGT(consensus) && ref != consensus){
	    if(charsHasAtLeast2Consensus())
		return true;
	}
	return false;
    }
    
    public String getSignatureString(){
	return this.contig + "_" + this.pos;
    }

    public String simpleString(){
	StringBuffer buffer = new StringBuffer(this.contig + "\t" + this.pos + "\t" + this.ref + "\t" + this.consensus + "\t" + this.getMutBase());
	return this.firstMutIndex + "\t" + buffer.toString();
    }
    
    public String simpleString(Hashtable lineNum2Sample){
	StringBuffer buffer = new StringBuffer(this.contig + "\t" + this.pos + "\t" + this.ref + "\t" + this.consensus + "\t" + this.getMutBase());
	return this.firstMutIndex + "\t" + lineNum2Sample.get(new Integer(this.firstMutIndex)) + "\t" + buffer.toString();
    }

    public String toGFFString(Hashtable lineNum2Sample,String synStr, String mutTypeStr, String trackName){
	char tmp = '.';
	if(synStr.equals("Syn"))
	    tmp = '+';
	else if(synStr.equals("N-Syn"))
	    tmp = '-';
	else if(synStr.equals("NC"))
	    tmp = '.';

	StringBuffer sb = new StringBuffer(this.contig + "\t\t" + trackName + "\t" + this.pos + "\t" + this.pos + "\t.\t" + tmp + "\t.\t" + this.getGFFAttribute(lineNum2Sample, mutTypeStr) + "Reference_Base=" + this.ref + ";Consensus=" + consensus + ";Modified_Call=" + this.getMutBase() +";snp_type=" + synStr + "\n");
	return sb.toString();
    }

    public String getGFFAttribute(Hashtable lineNum2Sample, String mutTypeStr){
	//return "ID=" + this.contig + "_" + this.pos + ";Name=" + this.contig + "_" + this.pos + ";"; 
	return "ID=" + this.contig+"_"+this.pos + ";Name="+ lineNum2Sample.get(new Integer(this.firstMutIndex)) + "_" + mutTypeStr + ";"; 
    }

    public String refConsString(){
	return this.contig + "\t" + this.pos + "\t" + this.ref + "\t" + this.consensus;
    }
    
    public String fakePolyString(int n){
	StringBuffer tmp = new StringBuffer(this.contig + " " + this.pos + " " + this.ref + " " + this.consensus);
	for(int i=0;i<n;i++)
	    tmp.append(" " + this.ref);
	return tmp.toString();
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

    public int getPos(){
	return this.pos;
    }

    public char getRef(){
	return this.ref;
    }

    public String getContig(){
	return this.contig;
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
     *  EDIT: sets the mutType field with proper mutType code returned from determineMutationtype();                                                                                                    
     *       then, it returns the index number of line(sample) in which the mutation has been observed.                                                                                                
     * return -1 if no mutation is found.
     */
    public int mutationType(){
	for(int i=0; i<chars.size();i++){
	    char tmp = chars.get(i).charValue();
	    if(tmp != '-' && isACGT(tmp) && tmp != this.consensus){
		this.mutType = determineMutationType(chars.get(i).charValue());
		this.firstMutIndex = (i+1);
		return i;
	    }
	}
	//System.out.println(this.getHeader());
	//return 0;
	this.mutType = 0;
	return -1;
    }

    public int getMutType(){
	return this.mutType;
    }

    public char getMutBase(){
	//System.out.println(getHeader());
	for(int i=0; i<chars.size();i++){
	    char tmp = chars.get(i).charValue();
	    if(tmp != '-' && isACGT(tmp) && tmp != this.consensus){
		this.firstMutIndex = (i+1);
		return chars.get(i).charValue();
	    }
	}
	return this.consensus;
    }

    private boolean isACGT(char c){
	if(c == 'A' || c == 'C'|| c == 'G' || c == 'T')
	    return true;
	return false;	    
    }

    public int determineMutationType(int c){
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

    public char[] getIndexArrayForLinesWithMutations(){
	char[] all = new char[chars.size()];
	int count = 0;
	for(int i=0;i<chars.size();i++){
	    if(this.chars.get(i).charValue() != this.consensus){
		all[i] = this.chars.get(i).charValue();
		count++;
	    }else
		all[i] = '-';
	}
	return all;
    }
}
