import java.io.*;
import java.util.*;

public class MergePileups{
    
    private Comparator<SimplePileupRecord> posSorter;
    private String emptyRecord;
    private ArrayList<BufferedReader> brs;             /* Contains readers for each pileup-formatted file */
    private SimplePileupRecord[] curList;
    private ArrayList<SimplePileupRecord> mins;
    private int nullBrCounter;
    private ArrayList<SimplePileupRecord> totalList;
    private BufferedWriter bw;
    private String outFileName;
    
    public static void main(String[] args){
	if(args.length >= 3){
	    String[] files = new String[args.length-2];
	    for(int i=1; i<(args.length-1); i++){
		files[i-1] = args[i];
	    }
	    System.out.println("Processing : " + args[0]);
	    new MergePileups(args[0], files, args[args.length-1]).build(30000);
	}else
	    System.err.println("USAGE: java MergerPileups <headerFileList, headers are in SAMFORMAT> <x1.pileup> <x2.pileup> <xn.pileup <outFileName>");
    } 


    /*NOT USED since we use headerList file instead of headerLoc*/
    public String[] grabHeaderFiles(String headerLoc){
	File[] files = new File(headerLoc).listFiles();
	String[] fileNames = new String[files.length];
	for(int i=0;i<files.length;i++){
	    fileNames[i] = files[i].getAbsolutePath();
	}
	return fileNames;
    }

    public String[] getHeaderFiles(String headerListFile){
	BufferedReader br = null;
	ArrayList<String> files = new ArrayList<String>();
	try{
	    br = new BufferedReader(new FileReader(headerListFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		File curFile = new File(curline);
		if(curFile.exists() && curFile.isFile())
		    files.add(curline);
		else{
		    br.close();
		    System.err.println("[MergePrintBases.getHeaderFiles] : Invalid file --> " + curline);
		    System.exit(1);
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	String[] fileArrays = new String[files.size()];
	return files.toArray(fileArrays);
    }

    public MergePileups(String headerListFile, String[] files, String outFileName){
	File tmp = new File(outFileName);
	if(tmp.exists()){
	    tmp.delete();
	}
	emptyRecord = "-";
	this.posSorter = new PosSorter(this.getHeaderFiles(headerListFile));
	brs = new ArrayList<BufferedReader>();
	
	/* Loads all pileup files into brs arrayList */
	for(int i=0; i<files.length; i++){
	    try{
		brs.add(new BufferedReader(new FileReader(files[i])));
	    }catch(IOException ioe){
		ioe.printStackTrace();
	    }
	}
	
	this.curList = new SimplePileupRecord[brs.size()];
	//System.err.println("@curListSize:\t" + curList.length);
	nullBrCounter = 0;
	this.totalList = new ArrayList<SimplePileupRecord>();
	initCurList();
	bw = null;
	this.outFileName = outFileName;
    }


    private void initCurList(){
	BufferedReader br = null;
	try{
	    for(int i=0; i<this.curList.length; i++){
		br = brs.get(i);
		String curline = br.readLine();
		if(curline != null){
		    //System.out.println(curline);
		    this.curList[i] = new SimplePileupRecord(curline.split("\\t"));
		    System.out.println(this.curList[i].toStringBuffer().toString());
		}else{
		    SimplePileupRecord.incrementReaderId();//since the file is empty we still need to increment the nextId counter to correclty assign readerIds.
		    //System.err.print("\t$CLOSING A FILE and INCREMENTING nullBrCounter:\t" + nullBrCounter + "\t" );
		    br.close();
		    nullBrCounter++;
		    //System.err.println(nullBrCounter);
		    brs.set(i,null);
		    br = null;
		    this.curList[i] = null;
		}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public void build(int mod){
	int count = 1;
	while(nullBrCounter < curList.length){
	    if(count%mod == 0){
		System.out.print(".");
		flushTotalList();
	    }
	    processEachLine();
	    count++;
	}
	if(totalList.size() > 0)
	    flushTotalList();
	try{
	    this.bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	System.out.println("DONE\n");
    }

    /* this outputs the content of totalList and flushes the totalList. DONE not to use too much memory */
    private void flushTotalList(){
	StringBuffer buffer = new StringBuffer();
	for(int i=0; i<totalList.size();i++){
	    buffer.append(totalList.get(i).toStringBuffer() + "\n");
	}
	try{
	    bw = new BufferedWriter(new FileWriter(this.outFileName, true));
	    bw.write(buffer.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	totalList = new ArrayList<SimplePileupRecord>();
    }
    
    private void processEachLine(){
	updateMins();
	if(this.mins.size() == 0)
	    ;
	else
	    generateFromMins();
    }

    /*
     * will set this.mins's size as 0 if no mins were found other than null 
     *
     * returns int[] filled with indices that needed to be updated in this.curList
     *
     * mins is the arrayList that keeps track of all the minimum ties. 
     */
    private void updateMins(){
	mins = new ArrayList<SimplePileupRecord>();
	mins.add(curList[0]);
	//System.err.println(">>>curMin\t" + (curList[0]==null?null:curList[0].getHeader()));
	for(int i=1; i<curList.length;i++){
	    int compareVal = posSorter.compare(mins.get(0),curList[i]);
	    if(compareVal == 0)
		mins.add(curList[i]); /* keep on adding the ties */
	    else if(compareVal > 0){
		mins = null;
		mins = new ArrayList<SimplePileupRecord>(); /* we initialize the minimal tie array since new minimum pos has been found */
		mins.add(curList[i]);
		//System.err.println(">>>NEWcurMin\t" + (curList[i]==null?null:curList[i].getHeader()));
	    }
	}
	
	if(mins.get(0) == null){
	    mins = new ArrayList<SimplePileupRecord>();
	    //System.err.println("###RETURNING SIZE O mins");
	}else
	    ;//System.err.println("####REURNING SIZE " + mins.size() + " mins with value:\t" + mins.get(0).getHeader());
    }

    private void generateFromMins(){
	//System.err.println("@@@@@@@@ generateFromMins");
	int offset = 0;
	int curPos = 0;
	SimplePileupRecord pr = new SimplePileupRecord(mins.get(0).getContig()
					   , mins.get(0).getPos()
					   , mins.get(0).getRef());
	
	for(int i=0; i<mins.size();i++){/* for each min entries in mins array*/
	    offset = mins.get(i).getReaderID(); //take the reader index which is identical to the column order
	    for(;curPos<offset;curPos++){ /* for any samples that don't have any data for this point */
		pr.appendToBases(this.emptyRecord); /* we simply append empty record which is just '*' */
	    }
	    pr.merge(mins.get(i));//now we can enter the current min entry to SimplePileupRecord obj.
	    curPos++;
	}
	/* this is to ensure the samples w/o data at this locus to be filled with empty after the last min entry*/
	for(; curPos < this.curList.length; curPos++){
	    pr.appendToBases(this.emptyRecord);
	}
	this.totalList.add(pr);
	feedIn();
    }
    
    /* this loads the next line */
    private void feedIn(){
	//System.err.println("FEED IN : mins size of " + mins.size() + ":\t" + mins.get(0).getHeader());
	for(int i=0; i<mins.size();i++){
	    int index = mins.get(i).getReaderID();
	    //System.err.println(index + " : br's len : " + brs.size() );
	    
	    BufferedReader br = brs.get(index);
	    try{
		if(br != null){
		    String curline = br.readLine();
		    //System.err.println("NEXT LINE FOR CURMIN\t)" + curline);
		    if(curline!=null)
			curList[index] = new SimplePileupRecord(curline.split("\\t"), index);
		    else{
			//System.err.print("\t$CLOSING A FILE and INCREMENTING nullBrCounter:\t" + nullBrCounter + "\t" );
			curList[index] = null;
			nullBrCounter++;
			//System.err.println(nullBrCounter);
			br.close();
			br = null;
			brs.set(index, null);
		    }
		}
	    }catch(IOException ioe){
		ioe.printStackTrace();
	    }
	}
    }

}
