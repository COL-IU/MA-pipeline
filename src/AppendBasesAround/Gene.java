public class Gene{

    private String contig;
    private int start;//smaller number
    private int end;//larger number
    private boolean fwd; //true if fwd, false otherwise.
    private String gName;

    public Gene(String ct, String line){
	String[] tokens = line.split("\\t");
    }

    public Gene(String ct, int st, int end, boolean fwd, String geneName){
	this.contig = ct;
	this.start = st;
	this.end = end;
	this.fwd = fwd;
	this.gName = geneName;
    }

    public Gene(String ct, int st, int end, boolean fwd, int frame){
	if(frame == 3)
	    frame = 0;
	if( (end-st+1)%3 != 0 ){
	    if(fwd){
		if( st%3 == frame)
		    end = end - ((end-st+1)%3);
		else if((end + 1) %3 == frame)
		    st = st + ((end-st+1)%3);
		else{
		    int oldst = st;
		    st = ( (st+1)%3 == frame ? st+1 : st+2);
		    end = ( (oldst+1)%3 == frame ? (end - ((end-oldst+1)%3) + 1) : (end- ((end-oldst+1)%3) + 2));
		}
	    }else{
		if( (end+1)%3 == frame)
		    st = st + ((end-st+1)%3);
		else if(st%3 == frame)
		    end = end - ((end-st+1)%3);
		else{
		    int oldst = st;
		    st = ( (st+1)%3 == frame ? st+1 : st+2);
		    end = ( (oldst+1)%3 == frame ? (end - ((end-oldst+1)%3) + 1) : (end- ((end-oldst+1)%3) + 2));
		}
	    }
	    //System.err.println(ct + "\t" + st + "\t" + end + "\t" + fwd + "\t" + frame); 
	}
	
	if( (end-st+1)%3 != 0)
	    System.err.println("Are you  Kidding me?");
	this.contig = ct;
	this.start = st;
	this.end = end;
	this.fwd = fwd;
	
	
    }
    
    public Gene(String ct, String[] fragGeneScanTokens){
	this(ct, 
	     Integer.parseInt(fragGeneScanTokens[0]), 
	     Integer.parseInt(fragGeneScanTokens[1]), 
	     ( fragGeneScanTokens[2].equals("+") ? true : false )
	     , Integer.parseInt(fragGeneScanTokens[3]));
    }

    public String getGName(){
	return this.gName;
    }
    public String getContig(){
	return this.contig;
    }
    
    public int getStart(){
	return this.start;
    }
    
    public int getEnd(){
	return this.end;
    }
    public boolean fwd(){
	return this.fwd;
    }
    public String directionStr(){
	if(this.fwd)
	    return "fwd";
	else
	    return "rev";
    }

    public String toString(){
	return this.contig + "\t" + this.start + "\t" + this.end + "\t" + this.fwd;
    }
    
    public boolean contains(int pos){
	if( (pos >= this.start) &&
	    (pos <= this.end) )
	    return true;
	return false;
    }

}
