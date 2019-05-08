import java.util.*;
import java.io.*;

public class PileupRecord{
    
    private String contig; /* FIELD 1 */
    private int pos;       /* FIELD 2 */
    private char refBase;  /* FIELD 3 - always stored as upper case in constructor*/ 
    private short depth;   /* FIELD 4 */
    private String bases;  /* FIELD 5 */
    private String quals;  /* FIELD 6 */
    private ExtComposition eComp; //base counts for AaCcGgTt Insetion Deletion, where caps represent fwd reads and lower cases represent rev reads.
    private ArrayList<InDel> indels; 
    private double refFreq;

    private static int nextId = 0;
    private int readerid; 

    private int start;
    private int end;
    
    /* IMPORTANT TO SET RIGHT VALUE for MIN_COUNT*/
    private static final int MIN_COUNT = 1; //since we are couting from 1 not from 0


    public String getContig(){
	return this.contig;
    }
    
    public int getPos(){
	return this.pos;
    }

    public ArrayList<InDel> getIndels(){
	return this.indels;
    }
    
    public char getRef(){
	return this.refBase;
    }

    public String getPileupString(){
	return contig + "\t" + pos + "\t" + refBase + "\t" + depth + "\t" + this.bases + this.quals;
    }

    public int getReaderID(){
	return this.readerid;
    }

    public PileupRecord(String line){
	this(line, nextId);
	nextId++;
    }

    /* constructor to load each field */
    public PileupRecord(String line, int index){
	this.start = 0;
	this.end = 0;
	//System.err.println("@inputRecord : " + line);
	this.readerid = index;
	String[] tokens = line.split("\\t");
	this.contig = tokens[0];
	this.pos = Integer.parseInt(tokens[1]);
	this.refBase = Character.toUpperCase(tokens[2].charAt(0)); /*refBase is always UPPER-case*/
	this.depth = Short.parseShort(tokens[3]);
	this.bases = tokens[4];
	this.quals = tokens[5];
	this.eComp = new ExtComposition(); /* initialize aAcCgGtT to ZERO */
	this.load();
	//System.err.println("@parsedAs : " + this.toString());
    }

        private void load(){
	/* hashtable to hold indel records, key is indelString concatenated with +/-*/
	Hashtable<String, InDel> tmptable = new Hashtable<String, InDel>();
	char lowerRef = Character.toLowerCase(this.refBase);
	int refcount = 0;
	for(int i=0;i<this.bases.length(); i++){
	    char base = this.bases.charAt(i);
	    if(base == '$'){
		this.end++;//ignore
	    }else if(base == '+' || base == '-'){//indels
		InDel tmp = new InDel(this.bases.substring(i));
		String tmpSig = tmp.getSignatureStr();
		if(tmptable.get(tmpSig) == null)//if there is no such indel in hashtable, we put a new entry
		    tmptable.put(tmpSig, tmp);
		else
		    tmptable.get(tmpSig).update(tmp);//otherwise, we just update the counter(fwd or rev)
		i = i+(tmpSig.length()-1);//update i
	    }else if(base == '^'){//start of a read, we ignore ^ and mapping score
		this.start++;
		i++;//there is a ascii character for mapping score after ^ symbol, so need to skip over.
	    }else{
		if(base == '.'){// base == ',')
		    this.eComp.update(this.refBase);
		    refcount++;
		}else if(base == ','){
		    this.eComp.update(lowerRef);
		    refcount++;
		}else
		    this.eComp.update(base);
	    }
	    
	}
	this.refFreq = (refcount*1.0d)/(this.depth*1.0d);
	this.indels = new ArrayList<InDel>(tmptable.values());
    }

    public boolean isThereIndels(){
	if(indels.size() > 0 )
	    return true;
	return false;
    }
    
    /* adjustedDepth is depth - |#($) + #(^)|*/
    public String isThereSignificantIndels(){
	int adjustedDepth = this.depth - (this.start + this.end);
	ArrayList<String> sigIndels = new ArrayList<String>();
	for(int i=0;i<indels.size();i++){
	    if(indels.get(i).isSignificantPortion(adjustedDepth))
		sigIndels.add(indels.get(i).getSignatureStr());
	}
	if(sigIndels.size() == 1)
	    return sigIndels.get(0);
	else if(sigIndels.size() > 1)
	    return "STRANGE";
	return null;
    }

    public String toString(){
	return this.contig + "\t" + this.pos + "\t" + this.refBase
	    + "\t" + this.depth + "\t" + this.bases + "\t" + this.quals;
    }

    
    /*public int run(){
	return this.processConsensus(this.getConsensus());
	}*/

    /* 
     * result[0] = classifier for mutation type: 0-same, 1-singleMutation, 2-polymorphism 
     * 
     */
    public int[] run2(StringBuffer simple){
	int[] result = new int[2];
	result[0] = this.processConsensus(this.getConsensus(), simple);
	result[1] = this.pos;
	return result;
	
    }


    public int getEffectiveDepth(int[] baseCounts){
	int count = 0;
	for(int i=0; i<baseCounts.length;i++){
	    if(baseCounts[i] >= MIN_COUNT)
		count += baseCounts[i];
	}
	return count;
    }

    /*
     * prints the result to standard out in the following format
     *  contigsome [TAB] 1bpPOS [TAB] RefBase [TAB] depth [TAB] pileupbases + Comma-delimited observed bases
     *
     * RETURNS: 0 if no mutation
     *          1 if point-mutation 
     *          2 if polymorphic
     *
     *         -1 otherwise
     */
    public int processConsensus(int[] baseCounts, StringBuffer simple){
	char[] consensus = getConsensus(baseCounts); // array of candidate bases
	StringBuffer tmpBf = new StringBuffer();
	if(consensus.length > 0){
	    tmpBf.append(this.contig + "\t" + this.pos + "\t" + this.refBase + "\t" + getEffectiveDepth(baseCounts) + "\t");
	    System.out.print(this.contig + "\t" + this.pos 
	    		     + "\t" + this.refBase + "\t" + this.depth 
	    		     + "\t" + this.bases + "\t");
	    int i=0;
	    //boolean multi = false;
	    for(; i<(consensus.length-1);i++){
		tmpBf.append(consensus[i] + ",");
		System.out.print(consensus[i]+",");
		//multi = true;
	    }
	    //if(!multi)
	    tmpBf.append(consensus[i]);
	    System.out.print(consensus[i]);	
		//if(multi)
		//	System.out.print(consensus[i] + "\t" );
	    
	    if(consensus.length > 1){
		tmpBf.append("\tP\t");
		tmpBf.append(baseCounts[0] + "\t" + baseCounts[1] + "\t" + baseCounts[2] + "\t" + baseCounts[3] + "\n");
		simple.append(tmpBf);
		System.out.println("\tP");
		return 2;//polymorphic
	    }else{
		if(consensus[0] == Character.toLowerCase(this.refBase)){
		    tmpBf.append("\tS\t");//tmpBf.append("\n");
		    tmpBf.append(baseCounts[0] + "\t" + baseCounts[1] + "\t" + baseCounts[2] + "\t" + baseCounts[3] + "\n");
		    simple.append(tmpBf);
		    System.out.println("\tS");
		    return 0;//same as ref
		}else{
		    tmpBf.append("\tM\t");
		    tmpBf.append(baseCounts[0] + "\t" + baseCounts[1] + "\t" + baseCounts[2] + "\t" + baseCounts[3] + "\n");
		    simple.append(tmpBf);
		    System.out.println("\tM");
		    return 1;//point mutation;
		}
	    }
	}else
	    return -1;
    }



    

    /* simply return int array of base counts in order of ACGT. base count has to be greater than or equal to MIN_COUNT to be counted
     * Otherwise, the specific base count to 0. ex) if the real count A C G T is 5 3 2 10 and MIN_COUNT is 3 then 5 3 0 10 will be returned
     */
    private int[] getConsensus(){
	return this.eComp.getSimplifiedACGTspectrum(MIN_COUNT);
    }
    

    /*
     * returns char arrary of bases that have passed a filtering. so this is a list of cadidates for consensus bases
     */
    private char[] getConsensus(int[] baseCounts){
	int len = 0;
	for(int i=0; i<baseCounts.length; i++){
	    if(baseCounts[i] > MIN_COUNT)
		len++;
	}
	char[] result = new char[len];
	int j =0;
	for(int i =0; i<baseCounts.length;i++){
	    if(baseCounts[i] > MIN_COUNT)
		result[j] = indexToDNABase(i);
	}
	return result;
    }

    private char indexToDNABase(int i){
	if(i == 0)
	    return 'A';
	else if(i == 1)
	    return 'C';
	else if(i == 2)
	    return 'G';
	else if(i == 3)
	    return 'T';
	else
	    return 'N';
    }
    
}
