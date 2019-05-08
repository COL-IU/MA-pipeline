public class InDel{

    private String indelStr; /* indelString written in caps*/
    private int fwd; /* # fwd reads having this specific indel type */
    private int rev; /* # rev reads having this specific indel type */
    private boolean insertion; /* insertion if true, deleteion otherwise */

    /* input string here is substring of pileupline starting with either +/- 
     * this will only process the first indel record it sees.
     */
    public InDel(String input){
	boolean legalFormat = true;
	if(input.charAt(0) == '+')
	    this.insertion = true;
	else
	    this.insertion = false;
	int iLen = input.length();
	int numDig = 0;
	StringBuffer buffer = new StringBuffer("");
	for(int i = 1; i<iLen; i++){
	    if(!Character.isDigit(input.charAt(i)))
		break;
	    buffer.append(input.charAt(i));
	    numDig++;
	}
	int len = 0;
	if(buffer.length() > 0){
	    len = Integer.parseInt(buffer.toString()); //length of indel
	    String tmp = input.substring(numDig+1,numDig+1+len); //extracts the indel sequence
	    if(tmp.length() == len){//lengths match!
		this.indelStr = tmp.toUpperCase(); // store it all capital
		if(Character.isUpperCase(tmp.charAt(0)))
		    this.fwd = 1;
		else
		    this.rev = 1;
	    }
	    else//length does not match
		legalFormat = false;
	}else{ /* not a indel format illegal format set indelStr to null and return*/
	    legalFormat = false;
	}
	if(legalFormat == false)
	    this.indelStr = "";
    }

    public void update(InDel other){
	this.fwd += other.fwd;
	this.rev += other.rev;
    }

    /* adjustedDepth is depth - |#($) + #(^)|*/
    public boolean isSignificantPortion(int adjustedDepth){
	double majorIndelFreq = ((this.fwd + this.rev) / (adjustedDepth*1.0d));
	if(this.fwd >= 3 && this.rev >= 3){
	    if( majorIndelFreq >= 0.5){
		return true;
	    }
	}else if(majorIndelFreq >= 0.8 && (this.fwd + this.rev) >= 5){
	    return true;
	}
	return false;

    }

    public int getFwd(){
	return this.fwd;
    }
    public int getRev(){
	return this.rev;
    }
    public int getDepth(){
	return this.fwd + this.rev;
    }
    public String getIndelStr(){
	return this.indelStr;
    }
    /* Signature str is  +/- then indel sequence*/
    public String getSignatureStr(){
	return (this.insertion ? "+" : "-") + this.indelStr.length() + this.indelStr;
    }

    public boolean isInsertion(){
	return this.insertion;
    }

    
    public void increment(boolean fwd){
	if(fwd)
	    this.fwd++;
	else
	    this.rev++;
    }
    
    public boolean isSameIndel(InDel other){
	if(this.indelStr.equals(other.getIndelStr()) && (this.insertion==other.isInsertion()))
	    return true;
	return false;
    }
    
}
