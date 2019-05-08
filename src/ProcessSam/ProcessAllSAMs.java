import java.io.*;
import java.util.*;

public class ProcessAllSAMs{

    /*
     * ls -d MH* | nohup xargs -I{} -P5 sh -c 'cd {};java -jar /scratch/heewlee/ebi_gut/mapping/ProcessSam.jar *_PE 0 0.2 ~/gutGenomes gList.txt ' > tanktank.log &
     *
     * args --> list of sam files, mapQualcutoff, percentMismatchAllowed(double) (ex: 0.01 for 1% mismatch, 0.05 for 5%), headerfileList (ex: util/M101.headers), genomeListFile
     */

    public static void main(String[] args){
	boolean scanOnly = false;
	int offset = 0;
	if(args.length < 1){
	    System.err.println("USAGE: java ProcessAllSAMs [-S (for SCANONLY to output median stats)] <samfile 1> ... <samfile n> <mapQualCutoff> <mismatchPercent> <headerfileList> <genomeListFile>");
	    System.exit(1);
	}
	if(args[0].startsWith("-")){
	    if(args[0].equals("-S")){ //scan only mode, only scan the PE files to get median stats
		offset = 1;
		scanOnly = true;
	    }else{
		System.err.println("USAGE: java ProcessAllSAMs [-S (for SCANONLY to output median stats)] <samfile 1> ... <samfile n> <mapQualCutoff> <mismatchPercent> <headerfileList> <genomeListFile>");
		System.exit(1);
	    }
	}
	String[] fs = new String[args.length-4];
	String[] rest = new String[4];
	int count = 0;
	for(int i=offset;i<args.length; i++){
	    if(i<args.length-4)
		fs[i-offset] = args[i];
	    else{
		rest[count] = args[i];
		count++;
	    }
	}
	/* samfile, paired or not, mapQualCutoff, mistmatchLimit, headerFileList, genomeListFile*/
	new ProcessAllSAMs(fs, true, Byte.parseByte(rest[0]), Double.parseDouble(rest[1]), rest[2], rest[3], scanOnly);
	
    }



    public final byte qThresh;
    private Hashtable<String, Long> hashOfGenomeSize;
    private Hashtable<String, MappingStat> hashOfStats;
    private Hashtable<String, String> hashOfHeaders;
    private Hashtable<String, MedianScan> hashOfMedianStats;
    
    String projectName;
    double mismatchAllowed;

    
    String representativeGenomeName; //ADDED [heewlee on 01/13/14] to merge PEs for multiple fasta referencew

    //ADDED [heewlee on 01/13/14] to merge PEs for multiple fasta referencew
    public void setGenomeName(String genomeListFile){
	if(genomeListFile.indexOf(File.separator) < 0)
	    this.representativeGenomeName = genomeListFile.substring(0,genomeListFile.indexOf("_gList"));
	else{
	    String tmp = genomeListFile.substring(genomeListFile.lastIndexOf(File.separator)+1);
	    this.representativeGenomeName = tmp.substring(0,tmp.indexOf("_gList"));
	}
    }
    

    public ProcessAllSAMs(String[] samfiles, boolean paired, byte mapQualCutoff, double mismatchPercent, String headerListFile, String genomeListFile, boolean scanOnly){
	
	this.setGenomeName(genomeListFile);//ADDED [heewlee on 01/13/14] to merge PEs for multiple fasta referencew

	this.projectName = samfiles[0].substring(0,samfiles[0].indexOf("_"));
	
	this.hashOfGenomeSize = new Hashtable<String, Long>(); /* key : genomeCode */
	this.hashOfHeaders = new Hashtable<String, String>();  /* key : genomeCod */
	this.hashOfMedianStats = new Hashtable<String, MedianScan>(); /* key : samfiles[i] (samfile name) */
	this.loadHashOfGenomeSize(genomeListFile);
	this.loadHashOfHeaders(headerListFile);
	//this.header = "";
		
	this.mismatchAllowed = mismatchPercent;

        if(paired){
	    qThresh = mapQualCutoff;
            this.hashOfStats = new Hashtable<String, MappingStat>();
	    if(scanOnly){
		this.parsePairsPEScan(samfiles);//load all median
		System.out.println("#fileName\tmedian\tMAD");
		for(int i=0;i<samfiles.length;i++){
		    this.hashOfMedianStats.get(samfiles[i]).toString(samfiles[i]);
		}
	    }else
		this.parsePairsPE(samfiles);
            //printStats(samfiles[0]);
        }else{
            qThresh = 0;
        }
    }

    public void parsePairsPE(String[] samfiles){
	//double[] medianMADpair = this.parsePairsPEScan(samfiles);
	this.parsePairsPEScan(samfiles); //load all median
	StringBuffer medMADStringBuffer = new StringBuffer();
	for(int i=0; i<samfiles.length;i++){
	    if(this.hashOfMedianStats.get(samfiles[i]).isEmptySet()){
		this.hashOfMedianStats.remove(samfiles[i]);
		System.err.println("\tSkipping : " + samfiles[i] + " due to (either no data or no proper pairs)");
	    }else{
		System.err.println("\tProcessing : " + samfiles[i] + " . . . (2nd iteration)" );
		double[] tmpMM = this.hashOfMedianStats.get(samfiles[i]).getMedMAD();
		String tmpMMstr = this.hashOfMedianStats.get(samfiles[i]).medMADToString();
		this.hashOfMedianStats.remove(samfiles[i]);
		this.parsePairPE2(samfiles[i], tmpMM);//this.parsePairPE(samfiles[i], medianMADpair);
		medMADStringBuffer.append(samfiles[i] + "\t" + tmpMMstr + "\n");
	    }
	}
	this.printStats(medMADStringBuffer.toString());
    }
    
    /* ms[0] --> median, ms[1] --> Median Absolute Deviation */
    private void parsePairPE(String sf1, double[] ms){
	BufferedReader br1 = null;
	String curLine1 = "";
	String curLine2 = "";
	
	try{
	    br1 = new BufferedReader(new FileReader(sf1));
	    SAM s1 = null;
	    SAM s2 = null;
	    while((curLine1=br1.readLine())!=null){
		if(!curLine1.startsWith("@")){
		    curLine2 = br1.readLine();
		    s1 = new SAM(curLine1, qThresh, this.mismatchAllowed);
		    s2 = new SAM(curLine2, qThresh, this.mismatchAllowed);
		    String gName1 = getGenomeAbrv(s1.genome);
		    String gName2 = getGenomeAbrv(s2.genome);
		    
		    int s1UnmappedVal = s1.flag & 0x0004;
		    int s2UnmappedVal = s2.flag & 0x0004;

		    //System.err.println("@gName1: ["+gName1+"]");
		    //System.err.println("@gName1: ["+gName2+"]");

		    if(s1UnmappedVal == 0){
			if(hashOfStats.get(gName1) == null)         /* if we dont have this genome in our hash, we add */
                            hashOfStats.put(gName1, new MappingStat(this.projectName+"_"+gName1+".qc.PE" , hashOfHeaders.get(gName1) , gName1, hashOfGenomeSize.get(gName1) ));
		    }

		    if(s2UnmappedVal == 0){
			if(hashOfStats.get(gName2) == null)         /* if we dont have this genome in our hash, we add */
			    hashOfStats.put(gName2, new MappingStat(this.projectName+"_"+gName2+".qc.PE" , hashOfHeaders.get(gName2) , gName2, hashOfGenomeSize.get(gName2) ));
		    }
		    
		    if( s1UnmappedVal == 0 ){                /* READ 1 MAPPED */
			if( s2UnmappedVal == 0){              /* READ 1 & 2 BOTH MAPPED */
			    if(gName1.equals(gName2))
				hashOfStats.get(gName1).updateStat(s1,s2, gName1, gName2, ms);
			    else{
				hashOfStats.get(gName1).splitSingle++;
				hashOfStats.get(gName2).splitSingle++;
			    }
			}else{                                    /* ONLY READ 1 MAPPED */
			    if(s1.unique){
				hashOfStats.get(gName1).uSingle++;
				if(s1.qcCompliant())
				    hashOfStats.get(gName1).qcuSingle++;
			    }else{
				hashOfStats.get(gName1).rSingle++;
				if(s1.qcCompliant())
				    hashOfStats.get(gName1).qcrSingle++;
			    }
			}
		    }else if( s2UnmappedVal == 0){           /* ONLY READ 2 MAPPED */
			if(s2.unique){
			    hashOfStats.get(gName2).uSingle++;
			    if(s2.qcCompliant())
				hashOfStats.get(gName2).qcuSingle++;
			}else{
			    hashOfStats.get(gName2).rSingle++;
			    if(s2.qcCompliant())
				hashOfStats.get(gName2).qcrSingle++;
			}
		    }else                                        /* NONE MAPPED */
			;
		}else
		    ;//headerBuffer.append(curLine1+"\n");
	    }
	    //this.header = headerBuffer.toString();
	    br1.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
    }

    private void parsePairPE2(String sf1, double[] ms){
	BufferedReader br1 =null;
	String curline = "";
	ArrayList<SAM> sams = new ArrayList<SAM>();
	try{
	    br1 = new BufferedReader(new FileReader(sf1));
	    while((curline=br1.readLine()).startsWith("@")){
		;
	    }
	    sams.add(new SAM(curline,qThresh, this.mismatchAllowed));
	    
	    while((curline=br1.readLine())!=null){
		String[] tmptoks = curline.split("\\t");
		if(tmptoks[0].equals(sams.get(sams.size()-1).readName)){//we only add reads that are paired into sams
                    sams.add(new SAM(curline, qThresh, this.mismatchAllowed));
                    //curCount++;                                                                                                                                     
                }else{//if start of new pair.
                    if(sams.size() == 2)
                        this.processDoublePE(sams,ms);
                    else if(sams.size() < 2){
                        ;//singleBuffer.append(sams.get(0).samline+"\n");//this.processSingleton(sams);//System.out.println("1:\t" + preName);                           
                    }else if(sams.size() >2){
                        ;/*for(int i=0;i<sams.size();i++){
                            splitBuffer.append(sams.get(i).samline+"\n");
			    }//processSplits(sams);                                                                                                                        */
                    }
                    sams = new ArrayList<SAM>();//prevline = curline;                                                                                                    
                    sams.add(new SAM(curline, qThresh, this.mismatchAllowed));
                    //curCount = 1;                                                                                                                                   
                }
            }
	    if(sams.size() == 2)
		this.processDoublePE(sams, ms);
	    sams = null;
	    br1.close();
	    
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    private void processDoublePE(ArrayList<SAM> sams, double[] ms){
	SAM s1 = sams.get(0);
	SAM s2 = sams.get(1);
	String gName1 = getGenomeAbrv(s1.genome);
	String gName2 = getGenomeAbrv(s2.genome);
	
	int s1UnmappedVal = s1.flag & 0x0004;
	int s2UnmappedVal = s2.flag & 0x0004;
	
	//System.err.println("@gName1: ["+gName1+"]");
	//System.err.println("@gName1: ["+gName2+"]");

	if(s1UnmappedVal == 0){
	    if(hashOfStats.get(gName1) == null)         /* if we dont have this genome in our hash, we add */
		hashOfStats.put(gName1, new MappingStat(this.projectName+"_"+gName1+".qc.PE" , hashOfHeaders.get(gName1) , gName1, hashOfGenomeSize.get(gName1) ));
	}
	
	if(s2UnmappedVal == 0){
	    if(hashOfStats.get(gName2) == null)         /* if we dont have this genome in our hash, we add */
		hashOfStats.put(gName2, new MappingStat(this.projectName+"_"+gName2+".qc.PE" , hashOfHeaders.get(gName2) , gName2, hashOfGenomeSize.get(gName2) ));
	}
	
	if( s1UnmappedVal == 0 ){                /* READ 1 MAPPED */
	    if( s2UnmappedVal == 0){              /* READ 1 & 2 BOTH MAPPED */
		if(gName1.equals(gName2))
		    hashOfStats.get(gName1).updateStat(s1,s2, gName1, gName2, ms);
		else{
		    hashOfStats.get(gName1).splitSingle++;
		    hashOfStats.get(gName2).splitSingle++;
		}
	    }else{                                    /* ONLY READ 1 MAPPED */
		if(s1.unique){
		    hashOfStats.get(gName1).uSingle++;
		    if(s1.qcCompliant())
			hashOfStats.get(gName1).qcuSingle++;
		}else{
		    hashOfStats.get(gName1).rSingle++;
		    if(s1.qcCompliant())
			hashOfStats.get(gName1).qcrSingle++;
		}
	    }
	}else if( s2UnmappedVal == 0){           /* ONLY READ 2 MAPPED */
	    if(s2.unique){
		hashOfStats.get(gName2).uSingle++;
		if(s2.qcCompliant())
		    hashOfStats.get(gName2).qcuSingle++;
	    }else{
		hashOfStats.get(gName2).rSingle++;
		if(s2.qcCompliant())
		    hashOfStats.get(gName2).qcrSingle++;
	    }
	}else                                        /* NONE MAPPED */
	    ;
    }
    
    
    /* process one samfile to load medianstat*/
    private void parsePairPEScan(String sf1, MedianScan ms){
	BufferedReader br1 = null;
	String curLine1 = "";
	String curLine2 = "";
	
	try{
	    br1 = new BufferedReader(new FileReader(sf1));
	    SAM s1 = null;
	    SAM s2 = null;
	    while((curLine1=br1.readLine())!=null){
		if(!curLine1.startsWith("@")){
		    curLine2 = br1.readLine();
		    s1 = new SAM(curLine1, qThresh, this.mismatchAllowed);
		    s2 = new SAM(curLine2, qThresh, this.mismatchAllowed);
		    String gName1 = getGenomeAbrv(s1.genome);
		    String gName2 = getGenomeAbrv(s2.genome);
		    if( ((s1.flag & 0x0004) == 0) 
			&& ((s2.flag & 0x0004) == 0) 
			&& (gName1.equals(gName2))
			&& s1.unique
			&& s2.unique ){//read 1 && 2 mapped && mapped to same genomes
			ms.updateStat(s1,s2);
		    }//everything else is junk in our case.
		}
		
	    }
	    br1.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
    }


    private void parsePairPEScan2(String sf1, MedianScan ms){
	BufferedReader br1 = null;
	BufferedWriter bfSingle = null;
	BufferedWriter bfSplit = null;
	StringBuffer singleBuffer = new StringBuffer();
	StringBuffer splitBuffer = new StringBuffer();
	//String curLine1 = "";
	//String curLine2 = "";
	//String[] prevline = {""};
	//String[] curline = "";
	ArrayList<SAM> sams = new ArrayList<SAM>();
	try{
	    br1 = new BufferedReader(new FileReader(sf1));
	    String line;
	    while((line=br1.readLine()).startsWith("@")){
	    }
	    //String[] tmptoks = tmp.split("\\t");
	    sams.add(new SAM(line, qThresh, this.mismatchAllowed));
	    //SAM s1 = null;
	    //SAM s2 = null;
	    //ArrayList<SAM> sams = new ArrayList<SAM>();
	    //int curCount = 1;
	    
	    
	    while((line=br1.readLine())!=null){
		//if(!line.startsWith("@")){
		String[] tmptoks = line.split("\\t");
		if(tmptoks[0].equals(sams.get(sams.size()-1).readName)){
		    sams.add(new SAM(line, qThresh, this.mismatchAllowed));
		    //curCount++;
		}else{
		    if(sams.size() == 2)
			this.processDouble(sams,ms);
		    else if(sams.size() < 2){
			singleBuffer.append(sams.get(0).samline+"\n");//this.processSingleton(sams);//System.out.println("1:\t" + preName);
		    }else if(sams.size() >2){
			for(int i=0;i<sams.size();i++){
			    splitBuffer.append(sams.get(i).samline+"\n");
			}//processSplits(sams);
		    }
		    sams = new ArrayList<>();//prevline = curline;
		    sams.add(new SAM(line, qThresh, this.mismatchAllowed));
		    //curCount = 1;
		}
		
	    }
	    if(sams.size() == 2)
		this.processDouble(sams, ms);
	    else if(sams.size() < 2){
		singleBuffer.append(sams.get(0).samline + "\n");
	    }else if(sams.size() > 2){
		for(int i=0;i<sams.size();i++){
		    splitBuffer.append(sams.get(i).samline+"\n");
		}//processSplits(sams);
	    }
	    br1.close();
	    
	    if(singleBuffer.length() > 0){
		bfSingle = new BufferedWriter(new FileWriter(sf1+".singleton"));
		bfSingle.write(singleBuffer.toString());
		bfSingle.close();
	    }
	    if(splitBuffer.length() > 0){
		bfSplit = new BufferedWriter(new FileWriter(sf1+".split"));
		bfSplit.write(splitBuffer.toString());
		bfSplit.close();
	    }
	    
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    private void processDouble(ArrayList<SAM> sams, MedianScan ms){
	String gName1 = getGenomeAbrv(sams.get(0).genome);
	String gName2 = getGenomeAbrv(sams.get(1).genome);
	if( ((sams.get(0).flag & 0x0004) == 0) 
	    && ((sams.get(1).flag & 0x0004) == 0) 
	    && (gName1.equals(gName2))
	    && sams.get(0).unique
	    && sams.get(0).unique ){//read 1 && 2 mapped && mapped to same genomes
	    ms.updateStat(sams.get(0), sams.get(1));
	}//everything else is junk in our case.
		    //}
    }

    public void parsePairsPEScan(String[] samfiles){
	//MedianScan ms = new MedianScan();
	/* for each sam file */
	for(int i=0;i<samfiles.length;i++){
	    System.err.println("Processing : " + samfiles[i] + " . . . (1st iteration)");
	    this.hashOfMedianStats.put(samfiles[i], new MedianScan()); /* added to have separate medianScan stat for each experiment */
	    this.parsePairPEScan2(samfiles[i], this.hashOfMedianStats.get(samfiles[i])); /* process each samfile*/
	}
	for(int i=0;i<samfiles.length;i++){
	    this.hashOfMedianStats.get(samfiles[i]).finalizeStat(samfiles[i]); // this returs double[] containing median and MAD
	}
    }


    /*                                                                                                                                                        
     */
    public String getGenomeAbrv(String genome){
        //if(genome.contains("_"))
	//    return genome.substring(0,genome.indexOf("_"));
        return genome;
    }
    
    public void loadHashOfGenomeSize(String genomeListFile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(genomeListFile));
	    String curLine  = "";
	    while((curLine=br.readLine())!=null){
		String[] tokens = curLine.split("\\t");
		System.err.println("GSHash Putting: [" +tokens[1] + "][" + tokens[3] +"]");
		this.hashOfGenomeSize.put(tokens[1], new Long(tokens[3]));
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public File[] getHeaderFiles(String headerListFile){
	BufferedReader br = null;
	ArrayList<File> files = new ArrayList<File>();
	try{
	    br = new BufferedReader(new FileReader(headerListFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		File curFile = new File(curline);
		if(curFile.exists() && curFile.isFile())
		    files.add(new File(curline));
		else{
		    br.close();
		    System.err.println("[ProcessAllSAMs.getHeaderFiles] : Invalid file --> " + curline);
		    System.exit(1);
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	File[] fileArrays = new File[files.size()];
	return files.toArray(fileArrays);
    }

    public void loadHashOfHeaders(String headerListFile){
	BufferedReader br = null;
	try{
	    
	    File[] files =  this.getHeaderFiles(headerListFile); /* get File Array of all files listed in headerListFile */
	    for(int i=0;i<files.length;i++){
		String fName = files[i].getName();
		//if(files[i].isFile() && fName.endsWith("_header.txt")){    /* we only want header files! */
		br = new BufferedReader(new FileReader(files[i]));
		String curLine = "";
		StringBuffer buffer = new StringBuffer();
		while((curLine=br.readLine())!=null){
		    buffer.append(curLine + "\n");
		}
		br.close();
		hashOfHeaders.put(fName.substring( 0 , fName.indexOf("_header.txt") ),buffer.toString()); /* ADD EACH header content to hash under genome key*/
		//}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    /* Removed and updated with the function above that takes the filename that has a list of header files to allow multiple fasta reference [heewlee 01/14/14]*/
    /*
    public void loadHashOfHeaders(String headerLoc){
	BufferedReader br = null;
	try{
	    File[] files =  new File(headerLoc).listFiles(); //grab all the files in headerLoc directory 
	    for(int i=0;i<files.length;i++){
		String fName = files[i].getName();
		if(files[i].isFile() && fName.endsWith("_header.txt")){    // we only want header files! 
		    br = new BufferedReader(new FileReader(files[i]));
		    String curLine = "";
		    StringBuffer buffer = new StringBuffer();
		    while((curLine=br.readLine())!=null){
			buffer.append(curLine + "\n");
		    }
		    br.close();
		    hashOfHeaders.put(fName.substring( 0 , fName.indexOf("_header.txt") ),buffer.toString()); // ADD EACH header content to hash under genome key
		}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	}*/

    public void printStats(String medMADStrings){
	Enumeration<String> msList = hashOfStats.keys();
	System.out.println(medMADStrings);
	int count = 0;

	//for each mappingStat stored that is for each fasta sequence
	while(msList.hasMoreElements()){
	    count++;
	    String curKey = msList.nextElement();
	    hashOfStats.get(curKey).close();
	    hashOfStats.get(curKey).printStat(this.projectName, medMADStrings);
	}
    

	if(count > 1){
	    System.err.println("\n###MultipleFiles Detected: MERGING ...");
	    this.mergePEs();
	}else{
	    this.mergePEs();
	}
    }

    public void mergePEs(){
	Enumeration<String> msList = hashOfStats.keys();
	ArrayList<BufferedReader> readers = new ArrayList<BufferedReader>();
	StringBuffer headers = new StringBuffer();
	while(msList.hasMoreElements()){
	    String curKey = msList.nextElement();
	    headers.append(this.hashOfHeaders.get(curKey));
	    try{
		readers.add(new BufferedReader(new FileReader(hashOfStats.get(curKey).fileName)));
		System.err.println("\t" + hashOfStats.get(curKey).fileName);
	    }catch(IOException ioe){
		ioe.printStackTrace();
	    }
	}
	
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(projectName+"_"+representativeGenomeName+".qc.PE"));
	    bw.write(headers.toString());
	    bw.flush();
	    for(int i=0;i<readers.size();i++){
		BufferedReader curReader = readers.get(i);
		String curline = "";
		int count = 0;
		StringBuffer bf = new StringBuffer();
		while((curline=curReader.readLine())!=null){
		    count++;
		    if(curline.charAt(0) != '@')
			bf.append(curline+"\n");
		    if(count % 3000 == 0){
			bw.write(bf.toString());
			bw.flush();
			bf = new StringBuffer();
		    }
		}
		curReader.close();
		
		if(count % 3000 != 0){
		    bw.write(bf.toString());
		    bw.flush();
		    bf = new StringBuffer();
		}
	    }
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

}
