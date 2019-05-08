public class SAMPair{
    private SAM read1;
    private SAM read2;
    final short unmapped = 0x0004;
    final short reverse  = 0x0010;
    
    private boolean qcuPCOCS;
    /*
     * orientation encoding!
     *
     * 1: correct ----1----> <-----2-----  ,  3: Opposite <------1----- -----2----->,  5 rightAll ----1---->  ----2---->,  7: leftAll <----1---- <-----2-----, 
     * 2: correct ----2----> <-----1-----  ,  4: Opposite <------2----- -----1----->,  6 rightAll ----2---->  ----1---->,  8: leftAll <----2---- <-----1-----,
     * -1: ----1-----> , -2: <-----1----, -3: -----2-----> , -4 <-----2-----
     * 0 : unmapped
     */
    private byte orientation; //
    private int insertSize;
                       
    public SAMPair(SAM r1, SAM r2){
	this.read1 = r1;
	this.read2 = r2;
	this.orientation = this.calcOrientation();
	this.insertSize = this.calcInsertSize();
	qcuPCOCS = false;
    }

    public SAM getRead1(){
	return this.read1;
    }
    
    public SAM getRead2(){
	return this.read2;
    }

    public boolean getQcuPCOCS(){
	return this.qcuPCOCS;
    }

    public byte getOrientation(){
	return this.orientation;
    }
    
    public int getInsertSize(){
	return this.insertSize;
    }

    public String printReadLengths(){
	return this.read1.readLength + "\t" + this.read2.readLength;
    }

    public void setQcuPCOCS(boolean q){
	this.qcuPCOCS = q;
    }
    
    public int sumOfReadLengths(){
	return this.read1.readLength + this.read2.readLength;
    }

    public String getSAMLines(){
	return this.read1.samline + "\n" + this.read2.samline +"\n";
    }


    private int calcInsertSize(){
	if(orientation == 1)
	    return (int)(this.read2.start + this.read2.readLength - this.read1.start);
	else if(orientation == 2)
	    return (int)(this.read1.start + this.read1.readLength - this.read2.start);
	else{
	    return 0;
	}
    }

    public boolean qcCompliant(){
	if(this.read1.qcCompliant() && this.read2.qcCompliant())
	    return true;
	return false;
    }

    public String toString(){
	StringBuffer sb = new StringBuffer(read1.readName + "\t" + read2.readName + "\t"); 
	sb.append(flag2pic() + "\t" + this.insertSize);
	
	return sb.toString();
    }
    
    public String flag2pic(){
	if(orientation == 1)
	    return "----1----> <----2----";
	else if(orientation == 2)
	    return "----2----> <----1----";
	else if(orientation == 3)
	    return "<----1---- ----2---->";
	else if(orientation == 4)
	    return "<----2---- ----1---->";
	else if(orientation == 5)
	    return "----1---->  ----2---->";
	else if(orientation == 6)
	    return "----2---->  ----1---->";
	else if(orientation == 7)
	    return "<----1---- <-----2-----";
	else if(orientation == 8)
	    return "<----2---- <-----1-----";
	else
	    return orientation+"";
    }

    /*
      0x0004 : unmapped
      0x0010 : strand is reverse
      0x0040 : first Read in a pair?
      0x0080 : second read in a pair?
      0x0200 : QC failure
      0x0800 : optical or PCR duplicate
    */
    private byte calcOrientation(){
	//System.out.println("READ1 :\t" + this.read1.toString() + "\n"
	//		   +"READ2 :\t" + this.read2.toString() );
	boolean r1fwd = true;
	boolean r1mapped = false;
	boolean r2fwd = true;
	boolean r2mapped = false;
	if( (this.read1.flag & unmapped) == 0){ // if mapped
	    r1mapped = true;
	    r1fwd = getDirection(this.read1.flag);
	}
	
	if( (this.read2.flag & unmapped) == 0){ // if mapped
	    r2mapped = true;
	    r2fwd = getDirection(this.read2.flag);
	}

	if(r1mapped && r2mapped){
	    if(r1fwd ^ r2fwd){// if they are opposite direction
		if(r1fwd && (read1.start <= read2.start) )
		    return 1;
		else if(r1fwd)
		    return 4;
		else if(read2.start <= read1.start)
		    return 2;
		else
		    return 3;
	    }else{ // same direction
		if(r1fwd && (read1.start <= read2.start) )
		    return 5;
		else if(r1fwd)
		    return 6;
		else if(read2.start <= read1.start)
		    return 8;
		else
		    return 7;
	    }
	}else if(r1mapped){
	    if(r1fwd)
		return -1;
	    else
		return -2;
	}
	else if(r2mapped){
	    if(r2fwd)
		return -3;
	    else 
		return -4;
	}
	else // unmapped
	    return 0;
    }

    private boolean getDirection(short flag){
	if((flag & reverse) == 0) //if it's forward
	    return true;
	return false;
    }
   
}
