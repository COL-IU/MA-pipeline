public class SAM{
    String readName;
    short flag;
    String genome;//this variable stores contig name in case of using multiple-contig-genomes
    int start;
    byte mQual;
    //int CIGAR;
    short readLength;
    byte editDistance;
    boolean unique;
    byte numMismatch;
    byte numGapOpen;
    byte numGapExtentions;

    double percentMismatch;
    
    String samline;

    final byte qThresh;
    
    public String toString(){
	return readName + "\t" + flag + "\t" + start;
    }

    public SAM(String samLine, byte mappingQualCutOff, double mismatchPer){
	this.samline = samLine;
	this.qThresh = mappingQualCutOff;
	this.percentMismatch = mismatchPer;
	this.parseSAMline(samLine);
    }

    /*public SAM(String[] tokens, byte mappingQualCutOff, double mismatchPer){
	this.samline = samLine;
	this.qThresh = mappingQualCutOff;
	this.percentMismatch = mismatchPer;
	this.parseSAMline(tokens);
	}*/

    public void parseSAMline(String sl){
	this.parseSAMline(sl.split("\\t"));
    }

    
    public void parseSAMline(String[] tokens){
	//String[] tokens = sl.split("\\t");
	byte mappingQual = Byte.parseByte(tokens[4]);
	this.flag = 0x0004;
	this.readLength = (short)tokens[9].length();
	this.readName = tokens[0];
	if(mappingQual >= this.qThresh){
	    //this.readName = tokens[0];
	    this.flag = Short.parseShort(tokens[1]);
	    this.genome = tokens[2];
	    this.start = Integer.parseInt(tokens[3]);
	    this.mQual = mappingQual;
	    this.editDistance = this.extractBWAtag(tokens, "NM");
	    this.numMismatch = this.extractBWAtag(tokens, "XM");
	    this.numGapOpen = this.extractBWAtag(tokens, "XO");
	    this.numGapExtentions = this.extractBWAtag(tokens, "XG");
	    this.unique = this.isUniqueMapping(tokens, true);
	    
	}
    }

    public boolean qcCompliant(){
	//double percentMismatch = 0.03;
	if( ((editDistance*1.0d)/(this.readLength*1.0d))
	    <= percentMismatch)
	    return true;
	return false;
    }

    public byte extractBWAtag(String[] tokens, String tag){
	for(int i=11; i<tokens.length; i++){
	    if(tokens[i].startsWith(tag))
		return Byte.parseByte(tokens[i].substring(tag.length()+3));
	}
	return 0;
    }
    
    public boolean isUniqueMapping(String[] tokens){
	for(int i=11; i<tokens.length; i++){
            if(tokens[i].startsWith("XT:A:"))
		if(tokens[i].substring("XT:A:".length()).equals("U"))
		    return true;
        }
	return false;
    }

    public boolean isUniqueMapping(String[] tokens, boolean mem){
	if(mem){
	    int alignscore = Integer.parseInt(tokens[13].substring(tokens[13].lastIndexOf(":")+1));
	    int suboptimalscore = Integer.parseInt(tokens[14].substring(tokens[13].lastIndexOf(":")+1));
	    if(alignscore > suboptimalscore)
		return true;
	}else{
	    return this.isUniqueMapping(tokens);
	}
	return false;
    }
}
