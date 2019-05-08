import java.io.*;
import java.util.*;

public class Codons{

    private Hashtable<String, Character> codon2amino;
    private Hashtable<Character, UsageCounter> amino2UsageCounter;
    private SubstitutionMatrix sm;
    private int[] atgcNonCoding = new int[4]; /* 0:A, 1:T, 2:G, 3:C */
    private int[] atgcCoding = new int[4];    

    public SubstitutionMatrix getSM(){
	return this.sm;
    }

    /* This calculates codon usage and */
    public static void main(String[] args){
	if(args.length >= 4){
	    String[] ptts = new String[args.length-3];
	    for(int i=3;i<args.length;i++){
		ptts[i-3] = args[i];
	    }
	    new Codons(args[0], args[2]).runCodonUsageCalculator(args[1], ptts);
	}else{
	    System.err.println("USAGE: java Codons <codon file> <multiFasta> <substitutionMatrix> <ptt.1> ... <ptt.n>  (AT LEAST 1 .ptt file requires)");
	}
    }
    
    private void runCodonUsageCalculator(String fasta, String[] ptts){
	this.calculateCodonUsage(new FastaReader().parseFasta(fasta), new PttParser().parseAll(ptts));
    }

    /*
     * codonfile is a tab-delimited file where 1st column is codon triplet in UPPER characters
     * and 2nd column is a amino acid symbol in UPPER character. 
     * symbol, * is used to denote a STOP codon
     * 
     */
    public Codons(String codonfile, String smFile){
	this.sm = new SubstitutionMatrixParser().parse(smFile);
	this.codon2amino = new Hashtable<String, Character>();
	this.amino2UsageCounter = new Hashtable<Character, UsageCounter>();
	this.loadHash(codonfile);
    }
        
    private void loadHash(String codonfile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(codonfile));
	    String curLine = "";
	    while((curLine=br.readLine())!=null){
		String[] tokens = curLine.split("\\t"); /* tokens[0] --> codon, tokens[1] --> amino */
		Character curChar = new Character(tokens[1].charAt(0));
		this.codon2amino.put(tokens[0], curChar);
		UsageCounter cur = this.amino2UsageCounter.get(curChar);
		if(cur == null){
		    //		    System.out.println("Here");
		    this.amino2UsageCounter.put(curChar, new UsageCounter(tokens[0]));
		}else{
		    //System.out.println("there");
		    this.amino2UsageCounter.get(curChar).addCodon(tokens[0]);
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public boolean isConservative(SubstitutionMatrix sm, int threshold, boolean greater, char origAmino, char aAmino){
	int dist = sm.getDistance(origAmino, aAmino);
	if(greater){
	    if(dist < threshold)
		return false;
	}else{
	    if(dist > threshold)
		return false;
	}
	return true;
    }

    public void calculateCodonUsage(Hashtable<String, Seq> seqHash, Hashtable<String, ArrayList<Gene>> geneHash){
	Enumeration<String> keys = seqHash.keys(); /* get names for all contigs */
	ArrayList<String> names = new ArrayList<String>(); /* Contig Names */
	int syn=0;
	int nonsyn=0;
	int cons=0;
	int noncons=0;
	int total=0;

	int nonsense = 0;
	int other = 0;

	while(keys.hasMoreElements())
	    names.add(keys.nextElement());
	/* for each contig */
	for(int i=0; i<names.size();i++){
	    Seq curSeq = seqHash.get(names.get(i));
	    ArrayList<Gene> genes = geneHash.get(names.get(i));
	    /* for each gene in contig i */
	    for(int j=0; j<genes.size();j++){
		//System.err.print("Processing:\t" + genes.get(j).getGName());
		String geneSeq = genes.get(j).getGeneSequence(curSeq);
		if(genes.get(j).getGName().equals("prfB")){
		    String tmp1 = new Gene("ecoli", 3033206, 3034228, false, "prfB").getGeneSequence(curSeq);
		    String tmp2 = new Gene("ecoli", 3034230, 3034304, false, "prfB").getGeneSequence(curSeq);
		    geneSeq = tmp2+tmp1;
		}
		    //3033206..3034228,3034230..3034304
		//System.err.println("\t" + geneSeq.length());
		for(int k=0;k<geneSeq.length();k=k+3){
		    //total=total+3;
		    //System.err.println("\t" + geneSeq.length());
		    String curCodon = geneSeq.substring(k,k+3).toUpperCase();
		    char curAmino = this.toAmino(curCodon);
		    this.amino2UsageCounter.get(new Character(curAmino)).update(curCodon);
		    for(int l=0;l<3;l++){
			/* current base in reference is curCodon.charAt(l) */
			char[] bases = getOthers(curCodon.charAt(l));
			for(int m=0;m<bases.length;m++){
			    total++;
			    char preAmino = this.toAmino(curCodon);
			    char afterAmino = this.toAmino(makeAlternate(curCodon, bases[m],l));
			    //System.out.println(curCodon + "\t" + makeAlternate(curCodon, bases[m], l));
			    
			    if(preAmino != '*' && afterAmino == '*'){ // if nonsense mutations
				nonsense++;
			    }else{ // else
				other++;
			    }

			    if(preAmino == afterAmino ){
				syn++;
				int mutType = calcMutType(curCodon.charAt(l), bases[m]);
				this.updateSynConSyn(mutType, true);
			    }else{
				//System.out.println(".");

				nonsyn++;
				if(Character.isUpperCase(curCodon.charAt(l)) && Character.isUpperCase(bases[m]))
				    ;
				else
				    System.out.println("#####################\n#####################\n###################");
				int mutType = calcMutType(curCodon.charAt(l), bases[m]);
				this.updateSynConSyn(mutType, false);
				if(isConservative(this.sm, 0, true, preAmino, afterAmino)){
				    cons++;
				    this.updateWithMutType(mutType, true);//spectrumForConservative[mutType]++;
				}else{
				    noncons++;
				    this.updateWithMutType(mutType, false);//spectrumForNonCons[mutType]++;
				}
			    }
			    
			    //moving this block inside nonsyn block to just count for nonsyn changes 7/6/12
			    /* 
			    if(isConservative(this.sm, 0, true, preAmino, afterAmino))
				cons++;
			    else
				noncons++;
			    */
					    
			}
		    }
		}
	    }
	}
	System.out.println("#Syn\t#NonSyn\t#Cons\t#NonCons\t#Total");
	System.out.println(syn + "\t" + nonsyn + "\t" +cons + "\t" + noncons + "\t" + total);

	System.out.println("\n#nonsense\t#others");
	System.out.println(nonsense + "\t" + other);

	System.out.println("Break down for Syns");
	System.out.println("AT > GC : \t" + this.spectrumForSyn[0]);
	System.out.println("GC > AT : \t" + this.spectrumForSyn[1]);
	System.out.println("AT > TA : \t" + this.spectrumForSyn[2]);
	System.out.println("GC > TA : \t" + this.spectrumForSyn[3]);
	System.out.println("AT > CG : \t" + this.spectrumForSyn[4]);
	System.out.println("GC > CG : \t" + this.spectrumForSyn[5]);

	System.out.println("Break down for non-syn");
	System.out.println("AT > GC : \t" + this.spectrumForNonSyn[0]);
	System.out.println("GC > AT : \t" + this.spectrumForNonSyn[1]);
	System.out.println("AT > TA : \t" + this.spectrumForNonSyn[2]);
	System.out.println("GC > TA : \t" + this.spectrumForNonSyn[3]);
	System.out.println("AT > CG : \t" + this.spectrumForNonSyn[4]);
	System.out.println("GC > CG : \t" + this.spectrumForNonSyn[5]);

	System.out.println("Break down for conservatives");
	System.out.println("AT > GC : \t" + this.spectrumForConservative[0]);
	System.out.println("GC > AT : \t" + this.spectrumForConservative[1]);
	System.out.println("AT > TA : \t" + this.spectrumForConservative[2]);
	System.out.println("GC > TA : \t" + this.spectrumForConservative[3]);
	System.out.println("AT > CG : \t" + this.spectrumForConservative[4]);
	System.out.println("GC > CG : \t" + this.spectrumForConservative[5]);

	System.out.println("Break down for non-conservative");
	System.out.println("AT > GC : \t" + this.spectrumForNonCons[0]);
	System.out.println("GC > AT : \t" + this.spectrumForNonCons[1]);
	System.out.println("AT > TA : \t" + this.spectrumForNonCons[2]);
	System.out.println("GC > TA : \t" + this.spectrumForNonCons[3]);
	System.out.println("AT > CG : \t" + this.spectrumForNonCons[4]);
	System.out.println("GC > CG : \t" + this.spectrumForNonCons[5]);

	printATGCCountForNonCoding(seqHash, names, geneHash);

	printCodonUsage();
    }
    

    public void printATGCCountForNonCoding(Hashtable<String, Seq> seqHash, ArrayList<String> names, Hashtable<String, ArrayList<Gene>> geneHash){
	Seq curSeq;
	ArrayList<Gene> curGenes;
	int[] atcgxArray = new int[5];
	int sum = 0;
	for(int i=0; i<names.size(); i++){
	    curSeq = seqHash.get(names.get(i));
	    curGenes = geneHash.get(names.get(i));
	    int start = 1;
	    int prevEnd = 0;
	    for(int j=0;j<curGenes.size();j++){
		int curSt = curGenes.get(j).getStart();
		int curEnd = curGenes.get(j).getEnd();
		if(curSt <= prevEnd)
		    System.out.println(curSt + ".." + curEnd);
		if(curSt > start){
		    //System.out.println(start + "\t" + (curSt-1));
		    atcgxArray = curSeq.getATCGXArray(start, curSt-1, atcgxArray);
		    sum += (curSt - start);
		}
		start = curEnd+1;
		prevEnd = curEnd;	    
	    }
	    if(start <= curSeq.getSeqLength()){
		//System.out.println(start + "\t" + curSeq.getSeqLength());
		atcgxArray = curSeq.getATCGXArray(start, curSeq.getSeqLength(), atcgxArray);
		sum += (curSeq.getSeqLength()+1 - start);
	    }
	    System.out.println(sum);
	    System.out.println("A\tT\tC\tG\tX");
	    System.out.println(atcgxArray[0] + "\t" + atcgxArray[1] + "\t" + atcgxArray[2] + "\t" + 	    atcgxArray[3] +"\t" + 	    atcgxArray[4] );
	    //seqHash.get(names.get(i)).printATCGCounts();

	}
    }

	
    public static String makeAlternate(String codon, char c, int pos){
	//System.out.print(codon + "\t" + c + "\t" + pos);
	char[] tmp = new char[3];
	for(int i=0; i<3;i++){
	    if(i==pos)
		tmp[i]=c;
	    else
		tmp[i]=codon.charAt(i);
	}
	//System.out.println("\t" + tmp[0]+tmp[1]+tmp[2]);
	return "" + tmp[0]+tmp[1]+tmp[2];
    }

    public static char[] getOthers(char c){
	char[] tmp2 = new char[3];
	if(c == 'A'){
	    char[] tmp = new char[] {'T','G','C'};
	    tmp2 = tmp;
	}
	else if(c == 'T'){
	    char[] tmp = new char[] {'A','G','C'};
	    tmp2 = tmp;
	}else if(c == 'G'){
	    char[] tmp = new char[] {'A','T','C'};
	    tmp2 = tmp;
	}else if(c == 'C'){
	    char[] tmp = new char[] {'A','T','G'};
	    tmp2 = tmp;
	}else
	    System.err.println("WTF");
	return tmp2;
    }

    public void printCodonUsage(){
	Enumeration<Character> aminos = amino2UsageCounter.keys();
	while(aminos.hasMoreElements()){
	    Character curChar = aminos.nextElement();
	    this.amino2UsageCounter.get(curChar).calculateUsage(curChar.charValue());
	}
    }

    
    public char toAmino(String triplet){
	
	Character ch =  this.codon2amino.get(triplet);
	if(ch == null){
	    System.err.println(triplet);
	    return '!';
	    //System.exit(1);
	}
	    
	return ch.charValue();
	    
	/*if(ch !=null)
	    return ch.charValue();
	else
	return '@';*/
    }


    

    /* RETURNS 1 if Syn, 0 if non-syn, -1, if NG due to bases other than atgc */
    /*public int isSyn(String t1, String t2){
	//System.out.println(t1 + "\t" + t2);
	char c1 = this.toAmino(t1);
	char c2 = this.toAmino(t2);
	if(c1 == '!' || c2 == '!')
	    return -1;
	if(c1 == c2)
	    return 1;
	else
	    return 0;
	    }*/

    public MutType isSyn(String t1, String t2){
	//System.out.println(t1 + "\t" + t2);
	char c1 = this.toAmino(t1);
	char c2 = this.toAmino(t2);
	if(c1 == '!' || c2 == '!')
	    return new MutType(-1,c1,c2,t1,t2);
	if(c1 == c2)
	    return new MutType(1,c1,c2,t1,t2);
	else
	    return new MutType(0,c1,c2,t1,t2);
    }




    /* ADDED to KEEP TRACK of more detailed COUNTS */
    /* length : 6, the order is  AT > GC , GC > AT , AT > TA , GC > TA , AT > CG , GC > CG */
    private int[] spectrumForConservative = new int[6]; //only 
    private int[] spectrumForNonCons = new int[6];
    private int[] spectrumForSyn = new int[6];
    private int[] spectrumForNonSyn = new int[6];


    private int calcMutType(char pre, char post){
	if(pre == 'A'){
	    if(post == 'C')
		return 3;
	    else if(post == 'G')
		return -2;
	    else if(post == 'T')
		return 1;
	    else
		return 0;
	}else if(pre == 'C'){
	    if(post == 'A')
		return 2;
	    else if(post == 'G')
		return 4;
	    else if(post == 'T')
		return -1;
	    else
		return 0;
	}else if(pre == 'G'){
	    if(post == 'A')
		return -1;
	    else if(post == 'C')
		return 4;
	    else if(post == 'T')
		return 2;
	    else
		return 0;
	}else if(pre == 'T'){
	    if(post == 'A')
		return 1;
	    else if(post == 'C')
		return -2;
	    else if(post == 'G')
		return 3;
	    else
		return 0;
	}else
	    return 0;
    }

    private void updateSynConSyn(int mutType, boolean syn){
	if(syn){
	    if(mutType < 0)
		this.spectrumForSyn[mutType+2]++; //transistions are -2 and -1
	    else if(mutType > 0)
		this.spectrumForSyn[mutType+1]++; //transversions are 1 through 4
	    else
		System.out.println("FUCKKKKKKKKKKKK");
	}else{
	    if(mutType < 0)
		this.spectrumForNonSyn[mutType+2]++; //transistions are -2 and -1
	    else if(mutType > 0)
		this.spectrumForNonSyn[mutType+1]++; //transversions are 1 through 4
	    else
		System.out.println("FUCKKKKKKKKKKKK!!!!");
	}
    }

    private void updateWithMutType(int mutType, boolean cons){
	if(cons){
	    if(mutType < 0)
		this.spectrumForConservative[mutType+2]++; //transistions are -2 and -1
	    else if(mutType > 0)
		this.spectrumForConservative[mutType+1]++; //transversions are 1 through 4
	    else
		System.out.println("FUCKKKKKKKKKKKK");
	}else{
	    if(mutType < 0)
		this.spectrumForNonCons[mutType+2]++; //transistions are -2 and -1
	    else if(mutType > 0)
		this.spectrumForNonCons[mutType+1]++; //transversions are 1 through 4
	    else
		System.out.println("FUCKKKKKKKKKKKK!!!!");
	}
	    
	//return mutType;
    }



}
