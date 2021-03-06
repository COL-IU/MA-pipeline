import java.io.*;
import java.util.*;

public class AnalyzeIndel{

    private Hashtable<String, ArrayList<Gene>> geneHash; /* KEY: contigName, VALUE: list of Genes in the contig */
    private Num2Sample num2sam;                          /* Index number to sample name hashing class           */
    private Hashtable<String, Seq> seqHash;              /* KEY: contigName, VALUE: Seq object of the contig    */
    private String outputFileHeader;
    
    //public AnalyzeIndel(String num2samFile, String fastaFile, String header, String pttList){
    //	this(num2samFile, fastaFile, header, this.getPttNames(pttList));
    //}
    
    //publice String getPttNames()

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
     * 1) It loads number 2 sample name hash then it loads Seq and gene model hash.
     * 2) processIndel is called to check if it's an polymorphism indel or not 
     * 3) Then if it's NOT a POLY-indel, it outputs result
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    public AnalyzeIndel(String num2samFile, String fastaFile, String header, String[] ptts){
	this.outputFileHeader=header;
	this.num2sam = new Num2Sample(num2samFile);
	this.loadHash(fastaFile, ptts);
    }

    private void loadHash(String fastaFile, String[] ptts){
	this.geneHash = new PttParser().parseAll(ptts);
	this.seqHash = new FastaReader().parseFasta(fastaFile);
    }

    /* returns NULL if non-coding, otherwise, it returns a Gene */
    public Gene getCodingGene(IndelConsensus ic){
	ArrayList<Gene> tmpGenes = this.geneHash.get(ic.getContig());
	Gene g = null;
	for(int i=0; i<tmpGenes.size(); i++){
	    if(tmpGenes.get(i).contains(ic.getPos())){
		g = tmpGenes.get(i);
		break;
	    }
	}
	return g;
    }


    public void process(String file){
	String headerline = "SampleName\tIndexNumber\tContig\tPosition\tIndelString\tGeneName if coding\tGene Start\tGene End\tGene Direction\t5'(including upto reported indel position)\t3'(from reportedPosition+1)";
	String sindelHeader = headerline + this.num2sam.getHeaders();
	BufferedReader br = null;
	BufferedWriter bw = null;
	String curline = "";
	try{
	    br = new BufferedReader(new FileReader(file));
	    StringBuffer pindelBuffer = new StringBuffer();
	    StringBuffer sindelBuffer = new StringBuffer();
	    StringBuffer polindelBuffer = new StringBuffer();
	    while((curline=br.readLine())!=null){
		processIndel(curline.split(" "), pindelBuffer, sindelBuffer, polindelBuffer);
	    }
	    br.close();
	    bw = new BufferedWriter(new FileWriter(this.outputFileHeader+".pindels"));
	    bw.write(headerline + "\n");
	    bw.write(pindelBuffer.toString());
	    bw.close();
		
	    bw = new BufferedWriter(new FileWriter(this.outputFileHeader+".sindels"));
	    bw.write(sindelHeader + "\n");
	    bw.write(sindelBuffer.toString());
	    bw.close();

	    bw = new BufferedWriter(new FileWriter(this.outputFileHeader+".polindels"));
	    bw.write(sindelHeader + "\n");
	    bw.write(polindelBuffer.toString());
	    bw.close();

	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public void processIndel(String[] tokens, StringBuffer pindelBuffer, StringBuffer sindelBuffer, StringBuffer polindelBuffer){
	IndelConsensus ic = new IndelConsensus(tokens);
	ic.setCodingGene(this.getCodingGene(ic));
	if(!ic.getPolyMorphismIndel()){ /* if it's a unique indel */
	    pindelBuffer.append(this.num2sam.mapNum2Sample(ic.getPindelIndex()) + "\t" + ic.toTabularString(this.seqHash.get(ic.getContig())) + "\n");
	    //System.out.println(this.num2sam.mapNum2Sample(ic.getPindelIndex()) + "\t" + ic.toTabularString(this.seqHash.get(ic.getContig())));
	}else if(ic.isSharedIndel()){
	    sindelBuffer.append(this.num2sam.mapNum2Sample(ic.getPindelIndex()) + "\t" + ic.toTabularString(this.seqHash.get(ic.getContig()), true) + "\n");
	    //System.err.println(this.num2sam.mapNum2Sample(ic.getPindelIndex()) + "\t" + ic.toTabularString(this.seqHash.get(ic.getContig()), true));
	}else if(ic.isPolyIndel()){
	    polindelBuffer.append(this.num2sam.mapNum2Sample(ic.getPindelIndex()) + "\t" + ic.toTabularString(this.seqHash.get(ic.getContig()), true) + "\n");
	}


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

    public static void main(String [] args){
	if(args.length == 5){
	    if(args[3].startsWith("-O:")){
		new AnalyzeIndel(args[1], args[2], args[3].substring(3), getPttFiles(args[4])).process(args[0]);
	    }else{
		System.err.println("outputFileHeader must be appended by \"-O:\"");
		System.err.println("USAGE: java AnalayzeIndel <.indel file> <num2sampleFile> <fastaFile> <-O:outputFileHeader ex) using output will result in output.pindel amd output,sindel> <pttListFile>");
	    }
	}else
	    System.err.println("USAGE: java AnalayzeIndel <.indel file> <num2sampleFile> <fastaFile> <-O:outputFileHeader ex) using output will result in output.pindel amd output,sindel> <pttListFile>");
	//new AnalyzeIndel().process(args[0]);
    }
    /*
    public static void main(String [] args){
	if(args.length > 4){
	    if(args[3].startsWith("-O:")){
		    String[] ptts = new String[args.length - 4];
		    for(int i=4; i< args.length;i++)
			ptts[i-4] = args[i];
		    new AnalyzeIndel(args[1], args[2], args[3].substring(3), ptts).process(args[0]);
	    }else{
		System.err.println("outputFileHeader must be appended by \"-O:\"");
		System.err.println("USAGE: java AnalayzeIndel <.indel file> <num2sampleFile> <fastaFile> <-O:outputFileHeader ex) using output will result in output.pindel amd output,sindel> <ptt.1> <ptt.2> ... <ptt.n>");
	    }
	}else
	    System.err.println("USAGE: java AnalayzeIndel <.indel file> <num2sampleFile> <fastaFile> <-O:outputFileHeader ex) using output will result in output.pindel amd output,sindel> <ptt.1> <ptt.2> ... <ptt.n>");
	//new AnalyzeIndel().process(args[0]);
	}*/

}



class IndelConsensus{
    
    private String contig;
    private int pos;
    private ArrayList<String> indels;
    private int[] indelIndices;
    private String consensus;
    private boolean isInsertion;
    private double polyRatio;
    private int size;
    private Gene g;
    private int pindelIndex;
    private boolean polyMorphismIndel;
    public static final int default_flanking_len = 20;
    public static final double POLY_CUTOFF_RATIO = 0.7;
    
    public String toTabularString(){
	return this.pindelIndex + "\t" + this.contig + "\t" + this.pos + "\t" + this.consensus + "\t" + ( (g!=null)? g.getGName() + "\t" + g.getStart() + "\t" + g.getEnd() +"\t"+ g.directionStr() : "-\t-\t-\t-" ); 
    }

    public String toTabularString(Seq seq){
	return toTabularString() + "\t" + get5PrimeSeq(default_flanking_len, seq) + "\t" + get3PrimeSeq(default_flanking_len, seq);
    }
    
    public String toTabularString(Seq seq, boolean printIndelStr){
	StringBuffer buffer = new StringBuffer();
	if(printIndelStr){
	    for(int i=0; i<indels.size();i++){
		buffer.append("\t" + indels.get(i));
	    }
	}
	return toTabularString(seq) + buffer.toString();
    }

    private String get5PrimeSeq(int n, Seq seq){
	return seq.getSubSeqWithDelim(pos+1-n,pos+1, "\t");
    }

    private String get3PrimeSeq(int n, Seq seq){
	return seq.getSubSeqWithDelim(pos+1,pos+1+n, "\t");
    }

    public void setCodingGene(Gene codingGene){
	this.g = codingGene;
    }
    
    public Gene getGene(){
	return g;
    }
    
    public String getContig(){
	return this.contig;
    }
    public int getPos(){
	return this.pos;
    }
    public ArrayList<String> getIndels(){
	return this.indels;
    }
    public int[] getIndelIndices(){
	return this.indelIndices;
    }
    public String getConsensus(){
	return this.consensus;
    }
    public boolean isInsertion(){
	return this.isInsertion;
    }
    public int getSize(){
	return this.size;
    }
    
    public IndelConsensus(String[] tokens){
	this.contig = tokens[0];
	this.pos = Integer.parseInt(tokens[1]);
	this.indels = new ArrayList<String>();
	for(int i=3;i<tokens.length;i++){
	    this.indels.add(tokens[i]);
	}
	this.consensus = determineConsensus();
	if(this.consensus.charAt(0) == '+')
	    this.isInsertion = true;
	else
	    this.isInsertion = false;
	this.size = determineSize();
	this.g = null;
	this.polyMorphismIndel = this.isPolyMorphismIndel();
    }

    public int getPindelIndex(){
	return this.pindelIndex;
    }

    public boolean getPolyMorphismIndel(){
	return this.polyMorphismIndel;
    }

    /* determine if 2 or more lines has such indels, otherwise returns line number (i+1)*/
    private boolean isPolyMorphismIndel(){
	int counter = 0;
	int lastIndex = 0;
	for(int i=0; i<indels.size(); i++){
	    if(indels.get(i).equals(this.consensus)){
		counter++;
		lastIndex = i;
	    }
	}
	
	this.polyRatio = (counter*1.0d) / (indels.size()*1.0d);

	if(counter > 1){
	    this.pindelIndex = -1;
	    return true;
	}else{
	    this.pindelIndex = (lastIndex+1);
	    return false;
	}
    }
    
    public boolean isSharedIndel(){
	if(this.polyMorphismIndel && this.polyRatio < this.POLY_CUTOFF_RATIO)
	    return true;
	return false;
    }

    public boolean isPolyIndel(){
	if(this.polyMorphismIndel && this.polyRatio >= this.POLY_CUTOFF_RATIO)
	    return true;
	return false;
    }
    
    private int determineSize(){
	int end = 0;
	for(int i=1; i<this.consensus.length(); i++){
	    if(!Character.isDigit(this.consensus.charAt(i)))
		break;
	    else
		end = (i+1);
	}
	return Integer.parseInt(consensus.substring(1,end));
    }

    private String determineConsensus(){
	Hashtable<String, Integer> countHash = new Hashtable<String, Integer>();
	for(int i=0;i<this.indels.size();i++){
	    if(!this.indels.get(i).equals("-")){
		if(countHash.get(this.indels.get(i)) == null){
		    countHash.put(this.indels.get(i), new Integer(1));
		}else{//updateCounter
		    countHash.put(this.indels.get(i), new Integer(countHash.get(this.indels.get(i)).intValue()+1)); // add 1.
		}	    
	    }
	}
	
	Enumeration<String> keys = countHash.keys();
	int curMax = 0;
	String maxKey = null;
	while(keys.hasMoreElements()){
	    String key = keys.nextElement();
	    int tmpCount = countHash.get(key);
	    if(tmpCount > curMax){
		curMax = tmpCount;
		maxKey = key;
	    }
	}
	return maxKey;
    }

}

