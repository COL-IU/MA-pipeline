import java.util.*;

public class SimplePileupIndelsRecord{
    
    private String contig;
    private int pos;
    private char refBase;
    private ArrayList<String> bases;

    private static int nextId = 0;
    private int readerid;
    
    public String getContig(){
	return this.contig;
    }
    
    public int getPos(){
	return this.pos;
    }
    
    public char getRef(){
	return this.refBase;
    }
    
    public String getPileupString(){
	return contig + "\t" + pos + "\t" + refBase + "\t" + this.bases;
    }

    public ArrayList<String> getBases(){
	return this.bases;
    }
    
    public int getReaderID(){
	return this.readerid;
    }

    public String getHeader(){
	return this.contig + "\t" + this.pos;
    }

    public SimplePileupIndelsRecord appendToBases(String b){
	this.bases.add(b);
	return this;
    }

    public SimplePileupIndelsRecord merge(SimplePileupIndelsRecord other){
	if(this.isSamePileupLocation(other))
	    appendToBases(other.getBases());
	return this;
    }

    public SimplePileupIndelsRecord appendToBases(ArrayList<String> basesL){
	this.bases.addAll(basesL);
	return this;
    }

    public StringBuffer toStringBuffer(){
	StringBuffer buffer = new StringBuffer(this.contig + "\t" + pos + "\t" + refBase);
	for(int i=0;i<this.bases.size();i++){
	    buffer.append("\t" + this.bases.get(i).toString());
	}
	return buffer;
    }

    public boolean isSamePileupLocation(SimplePileupIndelsRecord other){
	if(this.contig.equals(other.getContig()) && this.pos == other.getPos())
	    return true;
	return false;
    }

    public SimplePileupIndelsRecord(String contig, int pos, char ref){
	this.contig = contig;
	this.pos = pos;
	this.refBase = ref;
	this.bases = new ArrayList<String>();
	this.readerid = nextId;
	nextId++;
    }

    public SimplePileupIndelsRecord(String[] tokens, int index){
	this(tokens);
	this.readerid = index;
	nextId--;
    }
    
    public SimplePileupIndelsRecord(String[] tokens){
	this(tokens[0],Integer.parseInt(tokens[1]),tokens[2].charAt(0));
	if(tokens.length > 3){
	    bases.add(tokens[4]);
	}
    }

    
}
