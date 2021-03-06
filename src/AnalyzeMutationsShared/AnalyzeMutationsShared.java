import java.io.*;
import java.util.*;

public class AnalyzeMutationsShared{
    /* length : 6, the order is  AT > GC , GC > AT , AT > TA , GC > TA , AT > CG , GC > CG */
    /* Notice that first two are transitions and the rest are transversions */
    private int[][] mutationSpectrumCounts; 
    private int[][] mutationCounts;
    private int[][] synonymousCounts;//count for synonymous and non-synonymous and non-coding
    private int[][] conservativeCounts;
    
    private HashMap<String, Integer> seqNameToIndex;

    private Hashtable<String, Seq> seqHash;
    private Hashtable<String, ArrayList<Gene>> geneHash;
    private Codons codons;
    private int conservativeThresholdVal;
    private boolean greaterForConservative;
    private SubstitutionMatrix sm;
    private Hashtable<String, Consensus> pos2Consensus;
    private String jBrowseTrackName;
    
    private int numLines;
    private Hashtable<Integer, String> lineNum2SampleNameHash;

    /*    public static void main(String[] args){
	if(args.length > 8){
	    String[] ptts = new String[args.length-8];
	    for(int i=8;i<args.length;i++){
		ptts[i-8]=args[i];
	    }
	    new AnalyzeMutationsShared(args[1], args[2], Integer.parseInt(args[3]), args[4], args[5], args[6], args[7], ptts).process(args[0]);
	}else
	    System.err.println("USAGE: java AnalyzeMutationsShared <putationFile> <codonfile> <genomeFasta> <numberOfLines> <substitutionMatrix> <thresholdVal for conservative mutations> < 'Yy'(default)/'Nn' : greater(Yy) or less(Nn) than for conservative mutation> <lineNumToSampleName file - tab-delimited> <ptt.1> ... <ptt.n>\n If you choose 'Y/y' mutation is conservative if lookup value in substitutionMatrix is greater than give threshold value");
	    }*/

    /*
    public static void main(String[] args){
	if(args.length > 9){
	    String[] ptts = new String[args.length-9];
	    for(int i=9;i<args.length;i++){
		ptts[i-9]=args[i];
	    }
	    new AnalyzeMutationsShared(args[1], args[2], args[3] , args[4], args[5], args[6], args[7], args[8], ptts).process(args[0]);
	}else
	    System.err.println("USAGE: java AnalyzeMutationsShared <putationFile> <codonfile> <genomeFasta> <substitutionMatrix> <thresholdVal for conservative mutations> < 'Yy'(default)/'Nn' : greater(Yy) or less(Nn) than for conservative mutation> <lineNumToSampleName file - tab-delimited> <additionMutationFile : type 'N' if no additionalMutations> <jBrose TrackName: ex) mutL_SNP> <ptt.1> ... <ptt.n>\n If you choose 'Y/y' mutation is conservative if lookup value in substitutionMatrix is greater than give threshold value");
    }
    */
    public static void main(String[] args){
	if(args.length == 10){
	    new AnalyzeMutationsShared(args[1], args[2], args[3] , args[4], args[5], args[6], args[7], args[8], getPttFiles(args[9])).process(args[0]);
	}else
	    System.err.println("USAGE: java AnalyzeMutationsShared <putationFile> <codonfile> <genomeFasta> <substitutionMatrix> <thresholdVal for conservative mutations> < 'Yy'(default)/'Nn' : greater(Yy) or less(Nn) than for conservative mutation> <lineNumToSampleName file - tab-delimited> <additionMutationFile : type 'N' if no additionalMutations> <jBrose TrackName: ex) mutL_SNP> <pttListFile>\n If you choose 'Y/y' mutation is conservative if lookup value in substitutionMatrix is greater than give threshold value");
    }

    public static String[] getPttFiles(String pttList){
	BufferedReader br = null;
	ArrayList<String> list = new ArrayList<String>();
	try{
	    br = new BufferedReader(new FileReader(pttList));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		list.add(curline);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	String[] resultList = new String[list.size()];
	return list.toArray(resultList);
    }

    public AnalyzeMutationsShared(String codonfile, String genomeFasta, String substitutionMatrixFile, String thresholdVal, String greaterOrLess, String lineNum2SampleFile, String additionalMutationFile, String jBTrackName, String[] ptts){
	this.seqNameToIndex = new HashMap<String, Integer>(ptts.length*2);

	this.loadHash(genomeFasta,lineNum2SampleFile, additionalMutationFile, ptts);
	this.codons = new Codons(codonfile, substitutionMatrixFile);
	this.mutationSpectrumCounts = new int[ptts.length][6];
	this.mutationCounts = new int[ptts.length][numLines];
	this.synonymousCounts = new int[ptts.length][3];
	this.conservativeCounts = new int[ptts.length][2];
	
	this.conservativeThresholdVal = Integer.parseInt(thresholdVal);
	this.greaterForConservative = ( greaterOrLess.toLowerCase().equals("n") ? false : true);
	this.sm = new SubstitutionMatrixParser().parse(substitutionMatrixFile);
	this.jBrowseTrackName = jBTrackName;
    }

    private void loadHash(String genomeFasta, String lineNum2Sample, String additionalMutationFile, String[] ptts){
	this.seqHash = new FastaReader().parseFasta(genomeFasta);
	this.geneHash = new PttParser().parseAll(ptts);
	this.loadLineNumToSampleNameHash(lineNum2Sample);
	this.numLines = this.lineNum2SampleNameHash.size();
	this.loadPos2Consensus(additionalMutationFile);
	
	this.checkNumberOfSequences();
	this.loadSeqNameToIndexHash(ptts);
	//System.out.println("DONE LOADING");
    }
    
    private void loadSeqNameToIndexHash(String[] ptts){
	for(int i=0; i<ptts.length; i++){
	    this.seqNameToIndex.put(this.getSequenceNameFromPtt(ptts[i]), new Integer(i));
	}
    }
    
    private String getSequenceNameFromPtt(String pttName){
	int sepLoc = pttName.lastIndexOf(File.separator);
	if(sepLoc >= 0){
	    String tmp = pttName.substring(sepLoc+1);
	    return tmp.substring(0,tmp.indexOf(".ptt"));
	}else
	    return pttName.substring(0,pttName.indexOf(".ptt"));
    }
        
    //ADDED by heewlee 01/16/14 to handle multi-fasta reference
    private int checkNumberOfSequences(){
	int numSeq = this.seqHash.size();
	if(numSeq == this.geneHash.size())
	    return numSeq;
	else{
	    System.err.println("## FATAL ERROR in [AnalyzeMutation.checkNumberOfSequences()] #refSeqs is not consistent with #ptt files");
	    System.exit(0);
	    return -1;
	}	
    }

    /* additionMutFile format
       1st : contigName
       2nd : position
       3rd : refChar
       4th : consensus
       5th : lineIndex
       6th : mutBase
     */
    public void loadPos2Consensus(String additionalMutFile){
	this.pos2Consensus = new Hashtable<String, Consensus>();
	if(!additionalMutFile.equals("N")){
	    BufferedReader br = null;
	    try{
		String curline = "";
		br = new BufferedReader(new FileReader(additionalMutFile));
		while((curline=br.readLine())!=null){
		    String[] tokens = curline.split("\\t");
		    this.pos2Consensus.put(tokens[0]+"_"+tokens[1], new Consensus(tokens[0], Integer.parseInt(tokens[1]), tokens[2].charAt(0), tokens[3].charAt(0), Integer.parseInt(tokens[4]), tokens[5].charAt(0), this.numLines));
		}
		br.close();
	    }catch(IOException ioe){
		ioe.printStackTrace();
		System.exit(1);
	    }
	}
    }

    public StringBuffer[] finishPos2Consensus(){
	StringBuffer[] tmp = new StringBuffer[this.mutationSpectrumCounts.length];
	for(int i=0;i<tmp.length;i++){
	    tmp[i] = new StringBuffer();
	}
	

	Enumeration<Consensus> cs = this.pos2Consensus.elements();
	while(cs.hasMoreElements()){
	    Consensus c = cs.nextElement();
	    int curSeqIndex = seqNameToIndex.get(c.getContig()).intValue();
	    this.mutationCounts[curSeqIndex][c.mutationType()]++;
	    int tmpMT = c.getMutType();
	    if(tmpMT == 0)
		System.out.println(c.toString());
	    this.updateSpectrum(tmpMT, curSeqIndex);
	    MutType mt = this.updateSynonymous(c, curSeqIndex);
	    String consStr = "Non-conservative";
	    if(mt.getGene() == null || mt.getSynVal() != 0)
		consStr = "-";
	    else{
		if(mt.isConservative())
		    consStr = "Conservative";
	    }
	    tmp[curSeqIndex].append(c.simpleString(this.lineNum2SampleNameHash) + "\t" + mt.mutString() + "\t"+ mutTypeToString(tmpMT)
		       + "\t" + consStr +"\t" + synValToString(mt.getSynVal()) + "\n");
	}
	return tmp;
    }

    private String generateSharedHeaders(){
	StringBuffer bf = new StringBuffer();
	int s = this.lineNum2SampleNameHash.size();
	for(int i=0;i<s;i++){
	    bf.append("\t" + this.lineNum2SampleNameHash.get(new Integer(i+1)));
	}
	return bf.toString();
    }
    
    private String indexCharArrayToString(char[] ba){
	StringBuffer bf = new StringBuffer();
	for(int i=0; i<ba.length; i++){
	    bf.append("\t" + ba[i]);
	    /*if(ba[i])
	      bf.append(lineNum2SampleNameHash.get(new Integer(i+1)));
		else
	    bf.append("-");*/
	}
	return bf.toString();
    }

    public void process(String putationFile){
	BufferedReader br = null;
	BufferedWriter bw = null;
	String curline = "";
	StringBuffer[] mutTable = new StringBuffer[this.seqHash.size()];
	StringBuffer[] jbrowsegff = new StringBuffer[this.seqHash.size()];
	
	for(int i=0;i<this.seqHash.size();i++){
	    mutTable[i] = new StringBuffer("Index\tLine\tContig\tPosition\tRefBase\tConsensus\tMutBase\tbeforeTriplet\tafterTriplet\tbeforeAA\tafterAA\tgeneName\tgeneDirection\tMutType1\tMutType2\tConservative?\tSyn/non-syn" + generateSharedHeaders() + "\n");
	    jbrowsegff[i] = new StringBuffer();
	}

	try{
	    br = new BufferedReader(new FileReader(putationFile));
	    //System.out.println("Line\tContig\tPosition\tRefBase\tConsensus\tMutBase\tbeforeTriplet\tafterTriplet\tbeforeAA\tafterAA\tgeneName\tgeneDirection\tMutType1\tMutType2\tConservative?\tSyn/non-syn");
	    while((curline = br.readLine())!=null){
		//System.out.println(curline);
		Consensus c = new Consensus(curline.split(" "));
		int curSeqIndex = seqNameToIndex.get(c.getContig()).intValue(); //ADDED heewlee 1/16/14 multiple reference
		if(this.pos2Consensus.get(c.getSignatureString()) != null){
		    c = this.pos2Consensus.get(c.getSignatureString());
		    this.pos2Consensus.remove(c.getSignatureString());
		}
		//System.out.println(c.mutationType());
		this.mutationCounts[curSeqIndex][c.mutationType()]++;
		int tmp = c.getMutType();
		if(tmp == 0)
		    System.out.println(c.toString());
		this.updateSpectrum(tmp, curSeqIndex);
		/*int synval*/ MutType mt = this.updateSynonymous(c, curSeqIndex);
		String consStr = "Non-conservative";
		if(mt.getGene() == null || mt.getSynVal() != 0)
		    consStr = "-";
		else{
		    if(mt.isConservative())
			consStr = "Conservative";
		}
		mutTable[curSeqIndex].append(c.simpleString(this.lineNum2SampleNameHash) + "\t" + mt.mutString() + "\t"+ mutTypeToString(tmp) 
				+ "\t" + consStr +"\t" + synValToString(mt.getSynVal()) );// + "\t" + mt.dist);
		mutTable[curSeqIndex].append(indexCharArrayToString(c.getIndexArrayForLinesWithMutations()) + "\n");
		jbrowsegff[curSeqIndex].append(c.toGFFString(this.lineNum2SampleNameHash, synValToString(mt.getSynVal()), mutTypeToStringForJbrowse(tmp), this.jBrowseTrackName));
	    }
	    br.close();
	    
	    StringBuffer[] pos2ConsensBuffers = finishPos2Consensus();

	    Enumeration<String> seqNames = seqHash.keys();
	    while(seqNames.hasMoreElements()){
		String curName = seqNames.nextElement();
		int curIndex = seqNameToIndex.get(curName);
		outputResult(bw, mutTable[curIndex], jbrowsegff[curIndex], pos2ConsensBuffers[curIndex], curName, putationFile);
	    }
	    /*
	    this.writeBufferToFile(jbrowsegff.toString(), putationFile + ".jbrowse.gff");

	    mutTable.append(finishPos2Consensus());

	    bw = new BufferedWriter(new FileWriter(putationFile + ".distribution"));
            StringBuffer buffer = new StringBuffer();
            for(int i=0;i<mutationCounts.length;i++)
                buffer.append( (i+1) + "\t" + this.lineNum2SampleNameHash.get(new Integer(i+1)) + "\t" + this.mutationCounts[i] + "\n");
            bw.write(buffer.toString());
            bw.close();
	    System.err.println("Putations distribution saved to : " + putationFile + ".distribution");

	    bw = new BufferedWriter(new FileWriter(putationFile + ".detail"));
	    bw.write(mutTable.toString());
	    bw.close();
	    System.err.println("Detailed annotation of putations saved to : " + putationFile + ".detail");
	    */
	}catch(IOException ioe){
	    System.err.println("CAUGHT PROCESSING : " + curline);
	    ioe.printStackTrace();
	}
	//this.printAll(putationFile+".stat");
    }

    private void outputResult(BufferedWriter bw, StringBuffer mutTable, StringBuffer jbrowsegff, StringBuffer pos2ConsensBf, String seqName, String putationFile){

	try{
	    int curSeqIndex = seqNameToIndex.get(seqName);
	    
	    this.writeBufferToFile(jbrowsegff.toString(), putationFile + "." + seqName + ".jbrowse.gff");
	    
	    mutTable.append(pos2ConsensBf);
	    
	    bw = new BufferedWriter(new FileWriter(putationFile +"."+ seqName + ".distribution"));
	    StringBuffer buffer = new StringBuffer();
	    for(int i=0;i<mutationCounts[curSeqIndex].length;i++)
		buffer.append( (i+1) + "\t" + this.lineNum2SampleNameHash.get(new Integer(i+1)) + "\t" + this.mutationCounts[curSeqIndex][i] + "\n");
	    bw.write(buffer.toString());
	    bw.close();
	    System.err.println("Putations distribution saved to : " + putationFile + "." + seqName + ".distribution");
	    
	    bw = new BufferedWriter(new FileWriter(putationFile + "." + seqName + ".detail"));
	    bw.write(mutTable.toString());
	    bw.close();

	    bw = new BufferedWriter(new FileWriter(putationFile + "."+ seqName + ".detail.wBases"));
	    new AppendBasesAround(putationFile + "." + seqName + ".detail", this.seqHash, bw);//writes surrounding bases
	    bw.close();

	    System.err.println("Detailed annotation of putations saved to : " + putationFile +"." + seqName + ".detail");
	    
	    this.printAll(putationFile + "." + seqName + ".stat", curSeqIndex);
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }


    private void writeBufferToFile(String content, String outFileName){
	
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(outFileName));
	    bw.write(content);
	    bw.close();
	    bw = null;
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
    }
    
    private void loadLineNumToSampleNameHash(String f){
	this.lineNum2SampleNameHash = new Hashtable<Integer, String>();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(f));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		lineNum2SampleNameHash.put(new Integer(tokens[0]), tokens[1]);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public String synValToString(int synval){
	if(synval == 1)
	    return "Syn";
	else if(synval == 0)
	    return "N-Syn";
	else if(synval == -2)
	    return "NC";
	else
	    return "null";
    }
    public String mutTypeToString(int tmp){
	//AT > GC , GC > AT , AT > TA , GC > TA , AT > CG , GC > CG
	if(tmp == -2)
	    return "Transition\tAT>GC";
	else if(tmp == -1)
	    return "Transition\tGC>AT";
	else if(tmp == 1)
	    return "Transversion\tAT>TA";
	else if(tmp == 2)
	    return "Transversion\tGC>TA";
	else if(tmp == 3)
	    return "Transversion\tAT>CG";
	else if(tmp == 4)
	    return "Transversion\tGC>CG";
	else
	    return "null\tnull";
    }

    private String mutTypeToStringForJbrowse(int tmp){
	//AT > GC , GC > AT , AT > TA , GC > TA , AT > CG , GC > CG
	if(tmp == -2)
	    return "Ts_AT-GC";
	else if(tmp == -1)
	    return "Ts_GC-AT";
	else if(tmp == 1)
	    return "Tv_AT-TA";
	else if(tmp == 2)
	    return "Tv_GC-TA";
	else if(tmp == 3)
	    return "Tv_AT-CG";
	else if(tmp == 4)
	    return "Tv_GC-CG";
	else
	    return "null";
    }

    /* simply updates counts */
    private void updateSpectrum(int mutType, int curSeqIndex){
	if(mutType < 0)
	    this.mutationSpectrumCounts[curSeqIndex][mutType+2]++; //transistions are -2 and -1
	else if(mutType > 0)
	    this.mutationSpectrumCounts[curSeqIndex][mutType+1]++; //transversions are 1 through 4
	else
	    ;
	//return mutType;
    }
    
    /*
      returns synVal 1 if syn, 0 non-syn, -2: non-coding, returns something else(-1) if error
     */
    private MutType updateSynonymous(Consensus c, int curSeqIndex){
	String contig = c.getContig();
	char ref = c.getRef();
	char consBase = c.getConsensus();
	char mut = c.getMutBase();
	int pos = c.getPos();
	
	Seq tmpSeq = seqHash.get(contig);
	//System.out.println(contig);
	ArrayList<Gene> tmpGenes = geneHash.get(contig);
	Gene g = null;
	for(int i=0;i<tmpGenes.size();i++){
	    if(tmpGenes.get(i).contains(pos)){
		g = tmpGenes.get(i);
		break;
	    }
	}
	if(g == null){//this means non-coding
	    this.synonymousCounts[curSeqIndex][2]++;
	    return new MutType(-2, '-','-',"-","-",null);//-2;
	}else{
	    MutType mt = tmpSeq.isSyn((g.fwd() ? g.getStart() : g.getEnd()),g.fwd(), pos, mut, this.codons, consBase); //synVal: 1 if syn, 0 if non-syn
	    mt.setGene(g);
	    if(mt.getSynVal() == 1)
		this.synonymousCounts[curSeqIndex][0]++;
	    else if(mt.getSynVal() == 0){
		this.synonymousCounts[curSeqIndex][1]++;
		mt.isConservative(sm, conservativeThresholdVal, greaterForConservative);//this updates the boolean value of isConservative(true if conservative mut, false otherwise) WE CALL THIS ONLY if it's non-syn
		if(mt.isConservative())
		    this.conservativeCounts[curSeqIndex][0]++;
		else
		    this.conservativeCounts[curSeqIndex][1]++;
	    }else
		System.err.println("FUNKY");
	    return mt;
	}
    }

    /* AT>GC GC>AT AT>TA GC>TA AT>CG GC>CG Syn NonSyn NonCoding %Coding %NonCoding NonSyn/Syn Cons NonCons NonCons/Cons*/
    public void printAll(String statFileName, int curSeqIndex){
	System.out.println("AT > GC : \t" + this.mutationSpectrumCounts[curSeqIndex][0]);
	System.out.println("GC > AT : \t" + this.mutationSpectrumCounts[curSeqIndex][1]);
	System.out.println("AT > TA : \t" + this.mutationSpectrumCounts[curSeqIndex][2]);
	System.out.println("GC > TA : \t" + this.mutationSpectrumCounts[curSeqIndex][3]);
	System.out.println("AT > CG : \t" + this.mutationSpectrumCounts[curSeqIndex][4]);
	System.out.println("GC > CG : \t" + this.mutationSpectrumCounts[curSeqIndex][5]);
	System.out.println("Synonymous : \t" + this.synonymousCounts[curSeqIndex][0]);
	System.out.println("Non-Synonymous : \t" + this.synonymousCounts[curSeqIndex][1]);
	System.out.println("#####\t" + ((this.synonymousCounts[curSeqIndex][1]*1.0d)/(this.synonymousCounts[curSeqIndex][0]*1.0d)));
	System.out.println("Non-Coding : \t" + this.synonymousCounts[curSeqIndex][2]);
	System.out.println("Conservative : \t" + this.conservativeCounts[curSeqIndex][0]);
	System.out.println("Non-conservative : \t" + this.conservativeCounts[curSeqIndex][1]); 
	


	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(statFileName));
	    for(int i=0;i<6;i++){
		bw.write(this.mutationSpectrumCounts[curSeqIndex][i] + "\t");
	    }
	    bw.write(this.synonymousCounts[curSeqIndex][0]+"\t");
	    bw.write(this.synonymousCounts[curSeqIndex][1]+"\t");
	    bw.write(this.synonymousCounts[curSeqIndex][2]+"\t");
	    bw.write(""+ (((this.synonymousCounts[curSeqIndex][0]+this.synonymousCounts[curSeqIndex][1])*1.0d) / ((this.synonymousCounts[curSeqIndex][0]+this.synonymousCounts[curSeqIndex][1]+this.synonymousCounts[curSeqIndex][2])*1.0d))  ); // proportion of coding
	    bw.write("\t"+ ( (this.synonymousCounts[curSeqIndex][2]*1.0d) / ((this.synonymousCounts[curSeqIndex][0]+this.synonymousCounts[curSeqIndex][1]+this.synonymousCounts[curSeqIndex][2])*1.0d))  ); // proportion of non-coding
	    bw.write("\t"+(this.synonymousCounts[curSeqIndex][1]*1.0d)/(this.synonymousCounts[curSeqIndex][0]*1.0d));
	    bw.write("\t" + this.conservativeCounts[curSeqIndex][0]+"\t");
	    bw.write(this.conservativeCounts[curSeqIndex][1]+"\t");
	    bw.write((this.conservativeCounts[curSeqIndex][1]*1.0d)/(this.conservativeCounts[curSeqIndex][0]*1.0d) + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}

    }
}
