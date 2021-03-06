import java.io.*;
import java.util.*;

public class MappingStat{

    
    //StringBuffer allPairBuffer;
    //int pairCount;
    ArrayList<SAMPair> allPairs;
    
    int uPairCOCS; //uniquely mapped Pair Correct Orientation Correct   insert Size (within 3 std?)
    int uPairCOIS; //uniquely mapped Pair Correct Orientation Incorrect insert Size
    int uPairIO;   //uniquely mapped Pair Incorrect Orientation
    int uSingle;   //uniquely mapped Single
    int rPair;     //repeatedly mapped Pair
    int rSingle;   //repeatedly mapped Single                                                                                                                 
    int splitSingle;
    
    int qcuPairCOCS; //uniquely mapped Pair Correct Orientation Correct   insert Size (within 3 std?)
    int qcuPairCOIS; //uniquely mapped Pair Correct Orientation Incorrect insert Size
    int qcuPairIO;   //uniquely mapped Pair Incorrect Orientation
    int qcuSingle;   //uniquely mapped Single
    int qcrPair;     //repeatedly mapped Pair
    int qcrSingle;   //repeatedly mapped Single                                                                                                               
    long totalLength;
    long genomeLength;// = 5163189L;

    String key;
    String rLengths;

    BufferedWriter bw = null;
    String header;
    
    String fileName; //ADDED [heewlee on 01/13/14] to merge PEs for multiple fasta referencew

    MappingStat(String fileName, String header, String key){
	
	//System.err.println(fileName + "===" + header + "===" + key);

	this.fileName = fileName; //ADDED [heewlee on 01/13/14] to merge PEs for multiple fasta referencew
	
	this.key = key;

	this.allPairs = new ArrayList<SAMPair>();
	
	this.uPairCOCS = 0;
        this.uPairCOIS = 0;
        this.uPairIO = 0;
        this.uSingle = 0;
        this.rPair = 0;
        this.rSingle = 0;
        this.splitSingle = 0;
    	
        this.qcuPairCOCS = 0;
        this.qcuPairCOIS = 0;
        this.qcuPairIO = 0;
        this.qcuSingle = 0;
        this.qcrPair = 0;
        this.qcrSingle = 0;

        this.totalLength = 0L;
        this.genomeLength = 5000000L;
	this.header = header;
	try{
	    bw = new BufferedWriter(new FileWriter(fileName));
	    bw.write(this.header);
	    bw.flush();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    MappingStat(String fileName, String header, String key, long gLen){
	this(fileName, header, key);
	this.genomeLength = gLen;
    }

    private static final int MAX_SIZE = 20000;

    //processing pair of records that hit to same Genome
    void updateStat(SAM s1, SAM s2, String gName1, String gName2, double[] medMAD){
	this.totalLength += (s1.readLength + s2.readLength);
	SAMPair curPair = new SAMPair(s1,s2);
	if(s1.unique && s2.unique){
	    int ori = curPair.getOrientation();
	    int curISize = curPair.getInsertSize();
	    if(ori == 1 || ori == 2){//correct orientation
		if( (curISize >= (medMAD[0] - (3*medMAD[1])))
		    && (curISize <= (medMAD[0] + (3*medMAD[1]))) ){// within 3MAD ~= 2SD away
		    this.uPairCOCS++;
		    if(curPair.qcCompliant()){
			this.qcuPairCOCS++;
			curPair.setQcuPCOCS(true);
			this.allPairs.add(curPair);
			//this.pairCount++;
			if(this.allPairs.size()%MAX_SIZE == 0)
			    flushOutToFile();
		    }else
			this.uPairCOIS++;
		}else{
		    this.uPairCOIS++;
		    if(curPair.qcCompliant())
			this.qcuPairCOIS++;
		}
	    }else{//unexpected orientation
		this.uPairIO++;
		if(curPair.qcCompliant())
		    this.qcuPairIO++;
	    }
	}else{//at least one end read is not unique
	    this.rPair++;
	    if(curPair.qcCompliant())
		this.qcrPair++;
	}
    }

    void flushOutToFile(){
	StringBuffer buffer = new StringBuffer();

	try{
	    for(int i=0; i<allPairs.size();i++){
		if(i > 0 && (i%5000 == 0)){
		    bw.write(buffer.toString());
		    buffer = new StringBuffer();
		}
		buffer.append(allPairs.get(i).getSAMLines());
	    }
	    bw.write(buffer.toString());
	    buffer = null;
	    //bw.close();
	    this.allPairs = new ArrayList<SAMPair>();

	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    void close(){
	StringBuffer buffer = new StringBuffer();
	if(allPairs.size() < 1){
	    System.out.println("#NO PAIRS mapped");//System.err.println("allPairs dont have anything!!");
	}else{
	this.rLengths = this.allPairs.get(0).printReadLengths();
	try{
	    for(int i=0; i<allPairs.size();i++){
		if(i > 0 && (i%5000 == 0)){
		    bw.write(buffer.toString());
		    buffer = new StringBuffer();
		}
		buffer.append(allPairs.get(i).getSAMLines());
	    }
	    bw.write(buffer.toString());
	    buffer = null;
	    bw.close();
	    this.allPairs = new ArrayList<SAMPair>();

	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    }

    public void printStat(String projname, String medStats){
	System.out.println("============ for " + this.key + "  ============================\n" 
			   //+ "Mean      : " + this.mean + "\nSTD      : " + this.std + "\n"
			   //+ "Median    : " + medMAD[0] + "\nMAD      : " + medMAD[1] + "\n"
			   + "uPairCOCS : " + this.uPairCOCS + "\tqcuPairCOCS : " + this.qcuPairCOCS + "\n"
			   + "uPairCOIS : " + this.uPairCOIS + "\tqcuPairCOIS : " + this.qcuPairCOIS + "\n"
			   + "uPairIO   : " + this.uPairIO + "\tqcuPairIO   : " + this.qcuPairIO + "\n"
			   + "uSingle   : " + this.uSingle + "\tqcuSingle   : " + this.qcuSingle + "\n"
			   + "rPair     : " + this.rPair + "\tqcrPair     : " + this.qcrPair + "\n"
			   + "rSinlge   : " + this.rSingle + "\tqcrSinlge   : " + this.qcrSingle + "\n"
			   + "Mapped Tot: " + this.totalLength + " bps" + "\n"
			   + "Coverage  : " + ((this.totalLength*1.0d) / this.genomeLength) //+ "\n"
			   //+ "totalPair : " + this.totalPair
			   );
	this.outToFile(projname, this.statToString(medStats));
    }
	    
    private void outToFile(String projname, String content){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(projname + "_" + key +"_stat2.txt"));
	    bw.write(projname + "\t" + content + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    private String statToString(String medStats){//double[] medMAD){
	
	return medStats + "\n" + this.key + "\t" + this.rLengths + "\t" /*+ medMAD[0] + "\t" + medMAD[1] + "\t" */+ this.uPairCOCS + "\t" + this.qcuPairCOCS
	    + "\t" + this.uPairCOIS + "\t" + this.qcuPairCOIS + "\t" + this.uPairIO + "\t" + this.qcuPairIO
	    + "\t" + this.uSingle + "\t" + this.qcuSingle + "\t" + this.rPair + "\t" + this.qcrPair
	    + "\t" + this.rSingle + "\t" + this.qcrSingle
	    //+ "\t" + (this.totalPair - (this.uPairCOCS + this.uPairCOIS+ this.uPairIO + this.uSingle + this.rPair + this.rSingle))
	    + "\t" + this.totalLength + "\t" + this.genomeLength + "\t" + ((this.totalLength*1.0d) / this.genomeLength);
    }

}
