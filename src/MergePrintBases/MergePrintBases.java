import java.io.*;
import java.util.*;

public class MergePrintBases{

    private Comparator<PBCount> posSorter;
    private ExtComposition emptySet; /* placeholder for line with no data at given locus */
    private ArrayList<BufferedReader> brs; /* List of all the readers for each .printBases file */
    private PBCount[] curList; /* */
    private ArrayList<PBCount> mins; //keeps track of which lines have the minimum tieing record in terms of position
    private int nullBrCounter;
    private ArrayList<PBCount> totalList;
    private BufferedWriter bw;
    private String outFileName;


    public static void main(String[] args){
	if(args.length > 3){
	    String[] files = new String[args.length-2];
	    for(int i=1; i<(args.length-1); i++){
		files[i-1] = args[i];
	    }
	    System.out.println("Processing : " + args[0]);
	    new MergePrintBases(args[0], files, args[args.length-1]).build(30000);
	}else if(args.length == 3){
	    new MergePrintBases(args[0], loadFileNames(args[1]), args[2]).build(30000);
	}else{
	    System.err.println("USAGE: java MergePrintBases <headerListFile - headers are in SAMFORMAT> <x1.printBases> <x2.printBases> ... <xn.printBases> <outFileName>");
	    System.err.println("OR ");
	    System.err.println("java MergePrintBases <headerListFile - headers are in SAMFORMAT> <listingFile> <outFileName>");
	}
    }
    

    private static String[] loadFileNames(String listing){
	BufferedReader br = null;
	int size = 0;
	try{
	    br = new BufferedReader(new FileReader(listing));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		size++;
	    }
	    br.close();
	    String[] files = new String[size];
	    
	    int count = 0;
	    br = new BufferedReader(new FileReader(listing));
	    //String curline = "";
	    while((curline=br.readLine())!=null){
		
		files[count] = curline.trim() + File.separator + curline.trim() + ".printBases";
		count++;
	    }
	    br.close();
	    for(int i=0; i<files.length; i++){
		System.out.println(files[i]);
	    }
	    return files;
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return null;
    }


    public String[] grabHeaderFiles(String headerLoc){
	File[] files = new File(headerLoc).listFiles();
	ArrayList<String> fileNames = new ArrayList<String>();
	//String[] fileNames = new String[files.length];
	for(int i=0;i<files.length;i++){
	    if(files[i].isFile() && files[i].getName().endsWith("_header.txt"))
		fileNames.add(files[i].getAbsolutePath());
	}
	String[] fNames = new String[fileNames.size()];
	return fileNames.toArray(fNames);
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


    public MergePrintBases(String headerListFile, String[] files, String outFileName){
	this.deleteOutFile(outFileName);
	emptySet = new ExtComposition();
	this.posSorter = new PosSorter( this.getHeaderFiles(headerListFile) );
	brs = new ArrayList<BufferedReader>();
	
	for(int i=0; i<files.length;i++){
	    try{
		brs.add(new BufferedReader(new FileReader(files[i])));
	    }catch(IOException ioe){
		ioe.printStackTrace();
	    }
	}
	
	this.curList = new PBCount[brs.size()];
	nullBrCounter = 0;
	this.totalList = new ArrayList<PBCount>();
	initCurList();
	
	bw = null;
	this.outFileName = outFileName;;
	
    }

    private void deleteOutFile(String file){
	File f = new File(file);
	if(f.exists())
	    f.delete();
    }
    
    private void initCurList(){
	BufferedReader br = null;
	try{
	    for(int i=0; i<this.curList.length; i++){
		br = brs.get(i);
		String curline = br.readLine();
		if(curline !=null){
		    //String[] toks = curline.split("(\\t)|( )");
		    //if(toks.length == 3)
		    //this.curList[i] = new PBCount(toks[0], Integer.parseInt(toks[1]), toks[2].charAt(0));
		    //else{
		    //			System.out.println(curline);
		    this.curList[i] = new PBCount(curline.split("(\\t)|( )"));
		    //	System.out.println("PB id assigned : " + this.curList[i].getReaderID());
		    //}
		}else{
		    br.close();
		    nullBrCounter++;
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
	    if(count%mod == 0){//every "mod" generation, dump out what we have in the memory out to outputFile
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
	for(int i=0; i<totalList.size(); i++){
	    buffer.append(totalList.get(i).toStringBuffer() + "\n");
	}
	try{
	    bw = new BufferedWriter(new FileWriter(this.outFileName, true));
	    bw.write(buffer.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	totalList = new ArrayList<PBCount>();
    }

    private void processEachLine(){
	updateMins();
	if(this.mins.size() == 0)
	    ;
	else
	    generatePBCountFromMins();
    }
    
    /*
     * will set this.mins's size as 0 if no mins were found other than null 
     *
     * returns int[] filled with indices that needed to be updated in this.curList
     *
     * mins is the arrayList that keeps track of all the minimum ties. 
     */
    private void updateMins(){
	mins = new ArrayList<PBCount>();
	mins.add(curList[0]);
	for(int i=1; i<curList.length;i++){
	    int compareVal = posSorter.compare(mins.get(0),curList[i]);
	    if(compareVal == 0)
		mins.add(curList[i]);//keep on adding the ties
	    else if(compareVal > 0){
		mins = null;
		mins = new ArrayList<PBCount>();//we initialize the minimal tie array since new minimum pos has been found
		mins.add(curList[i]);
	    }
	}
	
	if(mins.get(0) == null)
	    mins = new ArrayList<PBCount>();
    }

    private void generatePBCountFromMins(){
	int offset = 0;
	int curPos = 0;
	PBCount pbct = new PBCount(mins.get(0).getContig(), mins.get(0).getPos(), mins.get(0).getRef());
	/* Note that mins array is sorted array in terms of reader indices */
	for(int i=0; i<mins.size();i++){/* for each min entries in mins array*/
	    offset = mins.get(i).getReaderID(); //take the reader index which is identical to the column order
	    for(;curPos<offset;curPos++){ /* for any samples that don't have any data for this point */
		pbct.appendToACGT(this.emptySet); /* we simply append empty set which is just 0 0 0 0 0 0 0 0 0 0*/
	    }
	    pbct.merge(mins.get(i));//now we can enter the current min entry to PBCount obj.
	    curPos++;
	}
	/* this is to ensure the samples w/o data at this locus to be filled with empty after the last min entry*/
	for(; curPos<this.curList.length;curPos++){
	    pbct.appendToACGT(this.emptySet);
	}
	this.totalList.add(pbct); 
	feedIn();
    }

    /* this loads the next line */
    private void feedIn(){
	//System.out.println("FEED IN");
	for(int i=0; i<mins.size();i++){
	    int index = mins.get(i).getReaderID();
	    // System.out.println(index + " : br's len : " + brs.size() );
	    
	    BufferedReader br = brs.get(index);
	    try{
		if(br != null){
		    String curline = br.readLine();
		    if(curline!=null)
			curList[index] = new PBCount(curline.split("(\\t)|( )"), index);
		    else{
			curList[index] = null;
			nullBrCounter++;
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
