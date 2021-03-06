import java.util.*;

/**
 * Seq is the base class for holding dna sequence datas.
 * 
 * @author Heewook Lee
 * @date Jan 25 2007
 * @updataed July 03, 2007
 * @ver 1.1
 */
public class Seq{
    
    private StringBuffer seqBuffer; // global variable for holding sequence data in StringBuffer Obj
    private String seqName; // golbal variable for holding sequence name in Sting Obj

    /*
     * Constructor
     *
     * @param fastaName name of the fasta seq
     * @param seq sequence in String obj
     * 
     */
    public Seq(String fastaName, String seq){
        seqName = fastaName;
        seqBuffer = new StringBuffer(seq);
    }
    
    /*
     * Constructor
     *
     * @param fastaName name of the fasta seq
     * @param seq sequence in String obj
     * 
     */
    public Seq(String fastaName, StringBuffer seq){
        seqName = fastaName;
        seqBuffer = seq;
    }
    
    /*
     * Constructor
     *
     * @param fastaName name of the fasta seq
     *
     */
    public Seq(String fastaName){
	seqName = fastaName;
	seqBuffer = new StringBuffer();
    }

    public Seq(){

    }

    /*
     * getSeqName() method simple returns the seqName
     *
     * @return seqName String obj seqName is returned
     */
    public String getSeqName(){
        return seqName;
    }
    
    /*
     * getSeq() method simple returns the sequence buffer
     *
     * @return seqBuffer StringBuffer obj sequence is returned
     */
    public StringBuffer getSeq(){
        return new StringBuffer(this.seqBuffer.toString());
    }
    
    public Seq deepCopy(){
	Seq newer = new Seq();
	newer.seqBuffer = this.getSeq();
	newer.seqName = this.seqName;
	return newer;
    }

    /*
     * formatSeq() method returns the sequence in fasta format without the seq name tag 
     * each line containts 50 neucleotide bases.
     *   
     * @return outBuffer StringBuffer obj sequnce is returned with newline character every 50 bases
     */
    public StringBuffer formatSeq(){
	StringBuffer outBuffer = new StringBuffer();
	int lowInd = 0;
	for(int i=0; i<seqBuffer.length();i++){
            if((i % 50) == 49){
                outBuffer.append(seqBuffer.substring(lowInd,i+1));
                outBuffer.append("\n");
                lowInd += 50;
            }
        }
	if( (seqBuffer.length()%50) != 0){
            outBuffer.append(seqBuffer.substring(lowInd));
            outBuffer.append("\n");
        }
        return outBuffer;
    }

    public ArrayList<IntegerCounter> countRuns(){
	System.err.println("TOTAL LEN : " + this.seqBuffer.length());
	ArrayList<IntegerCounter> runCountArr = new ArrayList<IntegerCounter>();
	int curRunLen = 0;
	char curRunChar = '0';
	for(int i=0;i<this.seqBuffer.length();i++){
	    char curChar = Character.toUpperCase(seqBuffer.charAt(i));
	    if(curChar == curRunChar)//if it's a currentRun
		curRunLen++;
	    else{/* need to update now*/
		if(curRunLen > 0){
		    if(runCountArr.size() < curRunLen)
			updateRunArr(runCountArr,curRunLen);
		    runCountArr.get(curRunLen-1).add1(curRunChar);
		}
		//if(curRunLen == 10)
		//    System.err.println(i + "\t" + curRunChar);
		curRunLen = 1;
		curRunChar = curChar;
	    }
	}
	if(curRunLen > 0 )
	    runCountArr.get(curRunLen-1).add1(curRunChar);
	
	return runCountArr;
    }



    /* ACGT, row<ACGT>, col<ACGT>
     *
     * This method outputs a 3 dimensional array: basically 4 2-dimensional array
     * where each array is 4 by 4 <ACGT> vs <ACGT> where row is the 5' neighbor base and col is the 3' neighbor base
     * 
     */
    public int[][][] getNeighboringBaseCountMatrix(){
	int[][][] nbcMatrix = new int[4][4][4];
	
	char startBase = Character.toUpperCase(this.seqBuffer.charAt(0));
	char endBase = Character.toUpperCase(this.seqBuffer.charAt(this.seqBuffer.length()-1));
	
	char prevBase = endBase;
	char curBase = startBase;//Character.toUpperCase(this.seqBuffer.charAt(0));
	char nextBase = '0';
	for(int i=1;i<this.seqBuffer.length();i++){
	    nextBase = Character.toUpperCase(this.seqBuffer.charAt(i));
	    nbcMatrix[this.baseToIndex(curBase)][this.baseToIndex(prevBase)][this.baseToIndex(nextBase)]++;
	    
	    prevBase = curBase;
	    curBase = nextBase;
	}
	nextBase = startBase;
	nbcMatrix[this.baseToIndex(curBase)][this.baseToIndex(prevBase)][this.baseToIndex(nextBase)]++;
	
	return nbcMatrix;
    }

    private int baseToIndex(char b){
	if(b == 'A')
	    return 0;
	else if(b == 'C')
	    return 1;
	else if(b == 'G')
	    return 2;
	else if(b == 'T')
	    return 3;
	else
	    return -1;
    }
    
    
    private void updateRunArr(ArrayList<IntegerCounter> arr, int curRunLen){
	while(arr.size() < curRunLen){
	    arr.add(new IntegerCounter());
	}
    }

    public ArrayList<Integer>[] getMotifLoci(String motif, boolean bothStrands){
	int fc = 0;
	int rc = 0;
	ArrayList<Integer>[] lociArr = null;
	if(bothStrands){
	    lociArr = new ArrayList[2];
	    lociArr[1] = new ArrayList<Integer>(); /* for rev strand */
	}else
	    lociArr = new ArrayList[1];
	lociArr[0] = new ArrayList<Integer>(); /* for fwd strand */

	int i = 0;
	while(i<this.seqBuffer.length()){
	    int tmp =this.seqBuffer.indexOf(motif, i);
	    if(tmp < 0)
		i = this.seqBuffer.length();
	    else{
		lociArr[0].add(new Integer(tmp+1));
		fc++;
		i = tmp + 1;
	    }
	}
	if(bothStrands){
	    StringBuffer rev = revSeqInternal();
	    i = 0;
	    while(i<rev.length()){
		int tmp = rev.indexOf(motif, i);
		if(tmp < 0)
		    i = this.seqBuffer.length();
		else{
		    lociArr[1].add(new Integer(rev.length() - tmp));
		    rc++;
		    i = tmp + 1;
		}
	    }
	}
	
	System.err.println("Fwd:\t" + fc);
	System.err.println("Rev:\t" + rc);

	return lociArr;
    }

    /*
     * appendSeq( String ) method appends the sequence 
     * at the end of the sequence held by Seq obj.
     *
     */
    public void appendSeq(String tempSeq){
        seqBuffer.append(tempSeq);
    }

    public void insertSeq(int offset, String str){
	seqBuffer.insert(offset, str);
    }

    public void insertSeq(int offset, char c){
	seqBuffer.insert(offset, c);
    }

    public void deleteSeqAt(int offset){
	seqBuffer.deleteCharAt(offset);
    }

    public void replaceSeqAt(int offset, char c){
	seqBuffer.replace(offset, offset+1, ""+c);
    }
    
    public void replaceSeqAt(int offset, String str){
	seqBuffer.replace(offset, offset+1, str);
    }
    
    public char seqAt(int index){
	return seqBuffer.charAt(index-1);
    }

    /*
     *
     * getSubSeq(int, int) method gets the substring of current seq
     * beginIndex is inclusive and end index is exclusive.
     */
    public String getSubSeq(int begin, int end){
      if (begin < 1){
        String right = seqBuffer.substring(0, end-1);
        String left = seqBuffer.substring(seqBuffer.length()+(begin-1), seqBuffer.length());
        return (left+right);
      } else if (end > seqBuffer.length()+1){
        String left = seqBuffer.substring(begin-1, seqBuffer.length());
        String right = seqBuffer.substring(0,end-seqBuffer.length()-1);
        return (left+right);
      }else{
        return seqBuffer.substring(begin-1, end-1);
      }
    }

    public String getSubSeqWithDelim(int begin, int end, String delim){
	String subseq = this.getSubSeq(begin,end);
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<subseq.length();i++){
	    bf.append(subseq.charAt(i));
	    if(i<(subseq.length()-1))
		bf.append("\t");
	}
	return bf.toString();
    }
    
    public int countGC(){
	int count = 0;
	for(int i=0; i<this.seqBuffer.length(); i++){
	    char curChar = seqBuffer.charAt(i);
	    if(curChar == 'G' || curChar == 'C' 
	       || curChar == 'g' || curChar == 'c')
		count++;
	}
	return count;
    }

    /*
     * revSeq() method returns the reverse complement sequence of the sequens held by the Seq Obj
     *
     * @return outBuffer StringBuffer holding the reverse complement seq of the seq held by the Seq Obj
     */
    public StringBuffer revSeq(){
        StringBuffer tempBuffer = revSeqInternal();
        StringBuffer outBuffer = new StringBuffer();
        int lowInd = 0;
        for(int i=0; i<tempBuffer.length();i++){
            if((i % 50) == 49){
                outBuffer.append(tempBuffer.substring(lowInd,i+1));
                outBuffer.append("\n");
                lowInd += 50;
            }
        }
	if( (tempBuffer.length()%50) != 0){
	    outBuffer.append(tempBuffer.substring(lowInd));
	    outBuffer.append("\n");
	}
        return outBuffer;
    }
    
    
    private StringBuffer revSeqInternal(){
        StringBuffer tempBuffer = seqBuffer.reverse();
        for(int i=0; i<tempBuffer.length();i++){
            char curChar = tempBuffer.charAt(i);
	    tempBuffer.setCharAt(i, rev(curChar));
        }
        return tempBuffer;
    }

    /*public SynType isSyn(int start, boolean fwd, int position, char mutbase, Codons codons, char consensus){
	
      }*/

    public MutType isSyn(int start, boolean fwd, int position, char mutbase, Codons codons){
	return this.isSyn(start,fwd,position,mutbase,codons,Character.toUpperCase(seqBuffer.charAt(position-1)));
    }
    
    
    /*
     * Changed to carry mutType information
     *
     */
    public MutType isSyn(int start, boolean fwd, int position, char mutbase, Codons codons, char consensus){
	try{
	    //System.out.print(mutbase + " ");
	    mutbase = Character.toUpperCase(mutbase);
	    //if(mutbase == consensus)
	    //		return 1;
	    //else{
	    int codonPos = this.getStartOfGene(start,fwd,position);
	    StringBuffer triple = new StringBuffer();
	    StringBuffer refTriple = new StringBuffer();
	    if(fwd){
		if(codonPos == 0){
		    refTriple.append(consensus);
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
		    triple.append(mutbase);
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
		}else if(codonPos == 1){
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    refTriple.append(consensus);
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    triple.append(mutbase);
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		}else if(codonPos == 2){
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    refTriple.append(consensus);
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    triple.append(mutbase);
		}
	    }else{//reverse
		if(codonPos == 0){
		    refTriple.append(this.rev(consensus));
		    refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
		    triple.append(this.rev(mutbase));
		    triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
		}else if(codonPos == 1){
		    refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    refTriple.append(this.rev(consensus));
		    refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    triple.append(this.rev(mutbase));
		    triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		}else if(codonPos == 2){
		    refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
		    refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    refTriple.append(this.rev(consensus));
		    triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
		    triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    triple.append(this.rev(mutbase));
		}
	    }
	    
	    return codons.isSyn(refTriple.toString(), triple.toString());
	    
		//}
	}catch(Exception e){
	    System.err.println("|" + start + "\t" + fwd + "\t" + position + "\t" + mutbase +"\t" +codons);//t start, boolean fwd, int position, char mutbase, Codons codons);
	    //System.err.println(ct.toStringBuffer().toString());
	    e.printStackTrace();
	    System.exit(1);
	}
	return null;
    }




    /*
     * start: in case of (-) stranded gene, start here is the actual start --> so 2nd column of fraggenescan output
     * in other words, start is always the 5' end of the gene(start position of the gene)
     * 
     * RETURN 1 if syn, 0 if non-synonymous, -1 if codon contains bases other than a t g c
     */
    /*public int isSyn(int start, boolean fwd, int position, char mutbase, Codons codons, char consensus){
	try{
	    //System.out.print(mutbase + " ");
	    mutbase = Character.toUpperCase(mutbase);
	    if(mutbase == consensus)
		return 1;
	    else{
		int codonPos = this.getStartOfGene(start,fwd,position);
		StringBuffer triple = new StringBuffer();
		StringBuffer refTriple = new StringBuffer();
		if(fwd){
		    if(codonPos == 0){
			refTriple.append(consensus);
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
			triple.append(mutbase);
			triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
			triple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
		    }else if(codonPos == 1){
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			refTriple.append(consensus);
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
			triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			triple.append(mutbase);
			triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    }else if(codonPos == 2){
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			refTriple.append(consensus);
			triple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
			triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			triple.append(mutbase);
		    }
		}else{//reverse
		    if(codonPos == 0){
			refTriple.append(this.rev(consensus));
			refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
			refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
			triple.append(this.rev(mutbase));
			triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
			triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
		    }else if(codonPos == 1){
			refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
			refTriple.append(this.rev(consensus));
			refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
			triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
			triple.append(this.rev(mutbase));
			triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    }else if(codonPos == 2){
			refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
			refTriple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
			refTriple.append(this.rev(consensus));
			triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
			triple.append(this.rev(Character.toUpperCase(seqBuffer.charAt(position))));
			triple.append(this.rev(mutbase));
		    }
		}
		
		return codons.isSyn(refTriple.toString(), triple.toString());
		
	    }
	}catch(Exception e){
	    System.err.println("|" + start + "\t" + fwd + "\t" + position + "\t" + mutbase +"\t" +codons);//t start, boolean fwd, int position, char mutbase, Codons codons);
	    //System.err.println(ct.toStringBuffer().toString());
	    e.printStackTrace();
	    System.exit(1);
	}
	return 0;
    }
    */

    private char rev(char curChar){
	if(curChar == 'T')
	    return 'A';
	else if(curChar == 't')
	    return 'a';
	else if(curChar == 'A')
	    return 'T';
	else if(curChar == 'a')
	    return 't';
	else if(curChar == 'G')
	    return 'C';
	else if(curChar == 'g')
	    return 'c';
	else if(curChar == 'C')
	    return 'G';
	else if(curChar == 'c')
	    return 'g';
	else
	    return curChar;
    }

    private String getCodon(boolean fwd, int codonPos, int position, char mutbase){
	if(codonPos == 0)
	    ;
	return new String();
    }
    
    /*
     *
     * @param
     * start    --> start position of target gene
     * fwd      --> true if fwd, rev otherwise
     * position --> position of polymorphism
     *
     * @returns
     * integer value of codon position for the given polymorphic pos.
     * 0 --> 1st pos of triplet
     * 1 --> 2nd pos of triplet
     * 2 --> 3rd pos of triplet
     *
     * start = 5
     * target position = 12
     * 
     * (fwd)       * (start pos)        ** (target pos)
     * (rev)       * (target pos)       ** (start pos)
     * 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18   --> position on genome
     * 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17   --> index 
     *
     * fwd strand
     * 12th position is the 2nd position of a codon.
     * --> (targetPos-start)%3 = (12-5)%3 = 1 
     *
     * rev strand
     * 5th pos is the 2nd pos of a codon 
     * --> (start - targetPos)%3 = (12-5)%3 = 1
     *
     */
    private int getStartOfGene(int start, boolean fwd, int position){
	if(fwd){//if the gene is on the fwd strand
	    return (position-start)%3;
	}else{//rev strand
	    return (start-position)%3;
	}
    }

    private int mapToRevIndex(int pos){
	return this.seqBuffer.length()-pos-1;
    }


}

class IntegerCounter{
    
    private int count;
    private int GCCount;
    private int ATCount;
    
    public IntegerCounter(){
	this(0);
    }
    
    public IntegerCounter(int initial){
	this.count = initial;
    }
    
    public int getCount(){
	return this.count;
    }
    
    public int getGCCount(){
	return this.GCCount;
    }

    public int getATCount(){
	return this.ATCount;
    }
    
    public void add1(){
	this.count++;
    }
    
    public void add1(char c){
	this.count++;
	if(c == 'G' || c == 'g' || c == 'C' || c == 'c')
	    this.GCCount++;
	else
	    this.ATCount++;
    }
    
    public void sub1(){
	this.count--;
    }
    
    public void addN(int n){
	this.count += n;
    }
    
    
}
