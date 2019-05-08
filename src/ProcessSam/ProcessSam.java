import java.io.*;
import java.util.*;


class MappingStat{

    int uPairCOCS; //uniquely mapped Pair Correct Orientation Correct   insert Size (within 3 std?)
    int uPairCOIS; //uniquely mapped Pair Correct Orientation Incorrect insert Size 
    int uPairIO;   //uniquely mapped Pair Incorrect Orientation
    int uSingle;   //uniquely mapped Single
    int rPair;     //repeatedly mapped Pair
    int rSingle;   //repeatedly mapped Single
    int splitSingle;
    //int totalPair; 

    int qcuPairCOCS; //uniquely mapped Pair Correct Orientation Correct   insert Size (within 3 std?)
    int qcuPairCOIS; //uniquely mapped Pair Correct Orientation Incorrect insert Size 
    int qcuPairIO;   //uniquely mapped Pair Incorrect Orientation
    int qcuSingle;   //uniquely mapped Single
    int qcrPair;     //repeatedly mapped Pair
    int qcrSingle;   //repeatedly mapped Single
    //private int qctotalPair; 

    long totalLength;
    long genomeLength;// = 5163189L;

    MappingStat(){
	
	this.uPairCOCS = 0;
  	this.uPairCOIS = 0;
	this.uPairIO = 0;
	this.uSingle = 0;
	this.rPair = 0;
	this.rSingle = 0;
	this.splitSingle = 0;
	//this.totalPair = 0;

	this.qcuPairCOCS = 0;
  	this.qcuPairCOIS = 0;
	this.qcuPairIO = 0;
	this.qcuSingle = 0;
	this.qcrPair = 0;
	this.qcrSingle = 0;

	this.totalLength = 0L;
	this.genomeLength = 5000000L;
	
    }

    MappingStat(long gLen){
	this();
	this.genomeLength = gLen;
    }
}

public class ProcessSam{

    //private int numOfPairs;
    //private int readLength;
    public final byte qThresh;
    private Hashtable<String, Long> hashOfGenomeSize;
    private Hashtable<String, ArrayList<SAMPair>> hashOfAllPairs;
    private Hashtable<String, MappingStat> hashOfStats;
    private ArrayList<SAMPair> allPairs;
    private ArrayList<SAM> r1single;
    private ArrayList<SAM> r2single;

    double mean;
    double std;
    double median;

    int totalPair;

    String header;
    
    String projectName;

    byte readLength;
    
    double mismatchAllowed;

    /*
     * samfiles : list of all samfiles being processed 
     * paired   : flag for paired or single ended
     * mapQualCutoff : minimum mapping quality cutoff (currently not used so it's set at 0)
     * prob     : this is the percentile cutoff value to calculate insert size.
     *
     */
    public ProcessSam(String[] samfiles, boolean paired, byte mapQualCutoff, double prob, String genomeListFile, double mismatchPercent){
	//this.projectName = samfiles[0].substring(0,samfiles[0].indexOf("_")) + "_" + new File(samfiles[0]).getParentFile().getName();
	this.projectName = samfiles[0];
 	
	this.hashOfGenomeSize = new Hashtable<String, Long>();
	this.loadHashOfGenomeSize(genomeListFile);
	this.header = "";
	this.totalPair = 0;
	//this.qctotalPair = 0;

	this.mean = 0;
        this.std = 0;
        this.median = 0;

	this.mismatchAllowed = mismatchPercent;

	if(paired){
	    //this.numOfPairs = 0;
	    //this.readLength = 0;
	    qThresh = mapQualCutoff;
	    this.hashOfAllPairs = new Hashtable<String, ArrayList<SAMPair>>();
	    //this.allPairs = new ArrayList<SAMPair>();
	    this.hashOfStats = new Hashtable<String, MappingStat>();


	    this.r1single = new ArrayList<SAM>();
	    this.r2single = new ArrayList<SAM>();
	    parsePairsPE(samfiles);

	    this.calcStat(prob);
	    
	    qcPairToFiles(samfiles[0]);
	    
	    printStats(samfiles[0]);
	}else{
	    qThresh = 0;
	}

    }

    public void loadHashOfGenomeSize(String genomeListFile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(genomeListFile));
	    String curLine  = "";
	    while((curLine=br.readLine())!=null){
		String[] tokens = curLine.split("\\t");
		this.hashOfGenomeSize.put(tokens[0], new Long(tokens[3]));
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    private void printStat(MappingStat ms, String abrv){
	StringBuffer buffer = new StringBuffer();
	System.out.println("============ for " + abrv + "  ============================\n" 
			   + "Mean      : " + this.mean + "\nSTD      : " + this.std + "\n"
			   + "Median    : " + this.median + "\n"
			   + "uPairCOCS : " + ms.uPairCOCS + "\tqcuPairCOCS : " + ms.qcuPairCOCS + "\n"
			   + "uPairCOIS : " + ms.uPairCOIS + "\tqcuPairCOIS : " + ms.qcuPairCOIS + "\n"
			   + "uPairIO   : " + ms.uPairIO + "\tqcuPairIO   : " + ms.qcuPairIO + "\n"
			   + "uSingle   : " + ms.uSingle + "\tqcuSingle   : " + ms.qcuSingle + "\n"
			   + "rPair     : " + ms.rPair + "\tqcrPair     : " + ms.qcrPair + "\n"
			   + "rSinlge   : " + ms.rSingle + "\tqcrSinlge   : " + ms.qcrSingle + "\n"
			   + "Mapped Tot: " + ms.totalLength + " bps" + "\n"
			   + "Coverage  : " + ((ms.totalLength*1.0d) / ms.genomeLength) //+ "\n"
			   //+ "totalPair : " + ms.totalPair
			   );
    }

    public void printStats(String projName){
	Enumeration<String> msList = hashOfStats.keys();
	while(msList.hasMoreElements()){
	    String curKey = msList.nextElement();
	    this.printStat(hashOfStats.get( curKey ), curKey);
	    this.outToFile(projName, curKey, statToString(curKey, hashOfStats.get(curKey)));
	}
    }



    public void qcPairToFiles(String baseName){
	Enumeration<String> gList = hashOfAllPairs.keys();
	while(gList.hasMoreElements()){
	    String curKey = gList.nextElement(); 
	    qcPairToFile(baseName + "." + curKey + ".qc.PE", hashOfAllPairs.get(curKey), curKey);
	}
    }


    private void qcPairToFile(String fName, ArrayList<SAMPair> pairs, String key){
	BufferedWriter bw = null;
	int count = 0;
	try{
	    bw = new BufferedWriter(new FileWriter(fName));
	    bw.write(this.header);
	    StringBuffer buffer = new StringBuffer();
	    if(pairs.size() == 0)
		bw.write("#NO PAIRS HAVE BEEN FOUND!!!!!\n");
	    for(int i=0; i<pairs.size();i++){
		if(pairs.get(i).getQcuPCOCS()){
		    count++;
		    if(count >= 5000){
			bw.write(buffer.toString());
			count = 0;
			buffer = new StringBuffer();
		    }
		    buffer.append(pairs.get(i).getRead1().samline + "\n" + pairs.get(i).getRead2().samline + "\n");
		}
	    }
	    bw.write(buffer.toString());
	    buffer = null;
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void outToFile(String projname, String abrv, String content){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(projname+"."+abrv+"_stat2.txt"));
	    bw.write(projname + "\t" + content + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    
    private String statToString(String abrv, MappingStat ms){
	
	/*int readLen = -1;
	if(r1single.size() > 0)
	    readLen = r1single.get(0).readLength;
	else if(r2single.size() > 0)
	readLen = r2single.get(0).readLength;*/
    
	return abrv + "\t" + this.readLength + "\t" + this.mean + "\t" +this.std + "\t" +  this.median + "\t" + ms.uPairCOCS + "\t" + ms.qcuPairCOCS
	    + "\t" + ms.uPairCOIS + "\t" + ms.qcuPairCOIS + "\t" + ms.uPairIO + "\t" + ms.qcuPairIO
	    + "\t" + ms.uSingle + "\t" + ms.qcuSingle + "\t" + ms.rPair + "\t" + ms.qcrPair
	    + "\t" + ms.rSingle + "\t" + ms.qcrSingle
	    //+ "\t" + (this.totalPair - (this.uPairCOCS + this.uPairCOIS+ this.uPairIO + this.uSingle + this.rPair + this.rSingle))
	    + "\t" + ms.totalLength + "\t" + ms.genomeLength + "\t" + ((ms.totalLength*1.0d) / ms.genomeLength);
    }

    private void calcStat(double prob){
	int sum = 0;
	//int properCount = 0;
	//int curSize = 0;

	Enumeration<ArrayList<SAMPair>> contigList = hashOfAllPairs.elements();
	this.allPairs = new ArrayList<SAMPair>();
	while(contigList.hasMoreElements()){
	    this.allPairs.addAll(contigList.nextElement()); 
	}
	//this.hashOfAllPairs = null;
	
	if(allPairs.size() > 0){
	    int[] iSizes = new int[allPairs.size()];
	    for(int i=0; i<allPairs.size(); i++){
		//if( (curSize = allPairs.get(i).getInsertSize()) > 0 )
		//	properCount++;
		iSizes[i] = allPairs.get(i).getInsertSize();//curSize;
		//sum += curSize;
	    }
	    
	    Arrays.sort(iSizes);
	    
	    int qtIndexUpper = (int)(iSizes.length * prob + 0.5);
	    int qtIndexLower = (int)(iSizes.length * (1.0-prob) + 0.5);
	    //System.err.println("ERrrrr");
	    //System.err.println(iSizes.length + "\t" + qtIndexLower + "\t" + qtIndexUpper);
	    int qUpper = iSizes[qtIndexUpper];
	    int qLower = iSizes[qtIndexLower];
	    for(int i=qtIndexLower; i<qtIndexUpper+1;i++){
		sum += iSizes[i];
	    }
	    
	    this.mean = (sum*1.0d) / ((qtIndexUpper+1-qtIndexLower)*1.0d);
	    
	    int numDatas = qtIndexUpper - qtIndexLower + 1;
	    if(numDatas%2 == 0)
		this.median = (iSizes[qtIndexLower+numDatas/2] + iSizes[qtIndexLower+numDatas/2+1])*1.0d / 2.0d;
	    else
		this.median = iSizes[qtIndexLower+numDatas/2];
	    
	    double squareOfErrorSum = 0;
	    for(int i=qtIndexLower; i<qtIndexUpper+1;i++){
		squareOfErrorSum += Math.pow(iSizes[i]*1.0d - this.mean, 2.0d);
	    }
	    this.std = Math.sqrt( squareOfErrorSum / (qtIndexUpper+1-qtIndexLower) );
	    
	    for(int i=0; i<allPairs.size();i++){
		SAMPair curPair = allPairs.get(i);
		MappingStat ms = hashOfStats.get(getGenomeAbrv(curPair.getRead1().genome));
		if(curPair.getRead1().unique && curPair.getRead2().unique){
		    if(curPair.getOrientation() == 1 || curPair.getOrientation() == 2){
			//if( (curPair.getInsertSize() >= (this.mean - this.std*2) )
			//    && (curPair.getInsertSize() <= (this.mean + this.std*2) )){
			if( (curPair.getInsertSize() >= (this.median*0.8) )
			    && (curPair.getInsertSize() <= (this.median*1.2) )){
			    ms.uPairCOCS++;
			    if(curPair.qcCompliant()){
				ms.qcuPairCOCS++;
				curPair.setQcuPCOCS(true);
			    }
			}else{
			    ms.uPairCOIS++;
			    if(curPair.qcCompliant())
				ms.qcuPairCOIS++;
			}
		    }else{
			ms.uPairIO++;
			if(curPair.qcCompliant())
			    ms.qcuPairIO++;
		    }
		}else{//both mapped, not unique
		    ms.rPair++;
		    if(curPair.qcCompliant())
			ms.qcrPair++;
		}
	    }
	}else{
	    System.out.println("NO Pairs have been found in :\t" + projectName);
	}
    }

    
    public void parsePairsPE(String[] samfiles){
	for(int i=0;i<samfiles.length;i++){
	    System.out.println("processing " + samfiles[i] + " ...");
	    this.parsePairPE(samfiles[i]);
	}
    }

    public void parsePairsSE(String[] samfiles){
	for(int i=0;i<samfiles.length;i=i+2){
	    this.parsePair(samfiles[i], samfiles[i+1]);
	}
    }

    /*
     */
    public String getGenomeAbrv(String genome){
	if(genome.contains("_"))
	   return genome.substring(0,genome.indexOf("_"));
	return genome;
    }

    public void parsePairPE(String sf1){
	BufferedReader br1 = null;
	String curLine1 = "";
	String curLine2 = "";
	StringBuffer headerBuffer = new StringBuffer();
	try{
	    br1 = new BufferedReader(new FileReader(sf1));
	    SAM s1 = null;
	    SAM s2 = null;
	    while((curLine1=br1.readLine())!=null){
		if(!curLine1.startsWith("@")){
		    curLine2 = br1.readLine();
		    s1 = new SAM(curLine1, qThresh, this.mismatchAllowed);
		    s2 = new SAM(curLine2, qThresh, this.mismatchAllowed);
		    String gName1 = getGenomeAbrv(s1.genome); //extract mapped genome names for read1 and read2
		    String gName2 = getGenomeAbrv(s2.genome);
		    if( (s1.flag & 0x0004) == 0){//read 1 mapped
			if(hashOfAllPairs.get(gName1) == null){ //if we dont have this genome in our hash
			    hashOfAllPairs.put(gName1, new ArrayList<SAMPair>()); //create a ArrayList<SAMPair>
			    hashOfStats.put(gName1, new MappingStat(hashOfGenomeSize.get(gName1))); // storage for stat as well
			}
			//hashOfStats.get(gName1).totalLength += s1.readLength;
			
			if((s2.flag & 0x0004) == 0){//read 2 mapped
			
			    if(hashOfAllPairs.get(gName2) == null){// if we dont have this genome in our hash
				hashOfAllPairs.put(gName2, new ArrayList<SAMPair>());
				hashOfStats.put(gName2, new MappingStat(hashOfGenomeSize.get(gName2)));
			    }
			    //hashOfStats.get(gName2).totalLength += s2.readLength;
			    
			    if( (s1.genome).equals(s2.genome) ){ // if read1 & 2 are both mapped to same genome
				hashOfStats.get(gName1).totalLength += (2*s1.readLength);
				//System.out.println(gName1 + " " + gName2 + "\t" + s1.samline + "\t" + s2.samline);
				hashOfAllPairs.get(gName1).add(new SAMPair(s1, s2)); // then as add them as a pair
			    }else{ // else that means we have split mapping so no need to store for now but keep track of the stat
				//r1single.add(s1);
				//r2single.add(s2);
				hashOfStats.get(gName1).splitSingle++;
				hashOfStats.get(gName2).splitSingle++;
			    }
			    
			}else{//ONLY read1 MAPPED
			    //r1single.add(s1);
			    if(s1.unique){
				hashOfStats.get(gName1).uSingle++;
				if(s1.qcCompliant())
				    hashOfStats.get(gName1).qcuSingle++;
			    }else{// not UNIQUE
				hashOfStats.get(gName1).rSingle++;
				if(s1.qcCompliant())
				    hashOfStats.get(gName1).qcrSingle++;
			    }
			}
		    }else if((s2.flag & 0x0004) == 0){//ONLY READ2 MAPPED
			if(hashOfAllPairs.get(gName2) == null){
			    hashOfAllPairs.put(gName2, new ArrayList<SAMPair>());
			    hashOfStats.put(gName2, new MappingStat(hashOfGenomeSize.get(gName2)));
			}
			//hashOfStats.get(gName2).totalLength += s2.readLength;
			//r2single.add(s2);
			if(s2.unique){
			    hashOfStats.get(gName2).uSingle++;
			    if(s2.qcCompliant())
				hashOfStats.get(gName2).qcuSingle++;
			}else{
			    hashOfStats.get(gName2).rSingle++;
			    if(s2.qcCompliant())
				hashOfStats.get(gName2).qcrSingle++;
			}
		    }
		    this.totalPair++;
		}else
		    headerBuffer.append(curLine1+"\n");
		//this.header = curLine1;
	    }
	    this.readLength = s1.readLength;
	    this.header = headerBuffer.toString();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void parsePair(String sf1, String sf2){
	BufferedReader br1 = null;
	String curLine1 = "";
	BufferedReader br2 = null;
	String curLine2 = "";
   
	try{
	    br1 = new BufferedReader(new FileReader(sf1));
	    br2 = new BufferedReader(new FileReader(sf2));
	    SAM s1 = null;
	    SAM s2 = null;
	    while( ((curLine1 = br1.readLine()) != null) &&
		   ((curLine2 = br2.readLine()) != null) ){
		
		s1 = new SAM(curLine1, qThresh, this.mismatchAllowed);
		s2 = new SAM(curLine2, qThresh, this.mismatchAllowed);
		
		if( (s1.flag & 0x004) == 0){
		    if((s2.flag & 0x004) == 0)
			allPairs.add(new SAMPair(s1, s2));
		    else
			r1single.add(s1);
		}else if((s2.flag & 0x004) == 0)
		    r2single.add(s2);
				
	    }

	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    /*
     * ls -d MH* | nohup xargs -I{} -P5 sh -c 'cd {};java -jar /scratch/heewlee/ebi_gut/mapping/ProcessSam.jar *_PE 0 0.999' > tanktank.log &
     *
     * args --> list of sam files, mapQualcutoff, prob(percentile cutoff val), percentMismatchAllowed(double) (ex: 0.01 for 1% mismatch, 0.05 for 5%)
     */
    public static void main(String[] args){
	String[] fs = new String[args.length-4];
	for(int i=0; i<args.length-4; i++){
	    fs[i] = args[i];
	}
	new ProcessSam(fs, true, Byte.parseByte(args[args.length-4]), Double.parseDouble(args[args.length-3]), args[args.length-2], Double.parseDouble(args[args.length-1]));
    }

}
