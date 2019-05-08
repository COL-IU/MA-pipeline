/*
 * printBases Count class to represent a object that holds alleleCount for a single locus in AaCcGgTt 8 column data with 3 columns headers
 * 1st col: contig/chromo Name
 * 2nd col: position
 * 3rd col: reference base
 * 4th - 11th : #AaCcGgTt where Upper case for forward reads and lower case for reverse reads.
 * FROM 12th col every 8cols to add next line of data.
 * each line of #AaCcGgTt is stored in an object called ExtComposition.
 */

import java.io.*;
import java.util.*;

public class PBCount{


    private static int nextId = 0;
    private int readerid;
    private String contig;
    private int pos;
    private char ref;
    private ArrayList<ExtComposition> acgtList;

    public PBCount(String ctg, int pos, char ref){
	this.readerid = nextId;
	nextId++;
	this.contig = ctg;
	this.pos = pos;
	this.ref = ref;
	this.acgtList = new ArrayList<ExtComposition>();
    }

    public PBCount(String[] tokens, int index){
	this(tokens);
	this.readerid = index;
	nextId--;
    }

    /* input is tokenized array of .printBases file*/
    public PBCount(String[] tokens){
	this(tokens[0],Integer.parseInt(tokens[1]),tokens[2].charAt(0));
	/*this.readerid = nextId;
	nextId++;
	this.contig = new String(tokens[0]);
	this.pos = Integer.parseInt(tokens[1]);
	this.ref = tokens[2].charAt(0);
	acgtList = new ArrayList<ExtComposition>();*/
	if(tokens.length > 4){
	    acgtList.add(new ExtComposition(Integer.parseInt(tokens[4]), Integer.parseInt(tokens[5]),
					    Integer.parseInt(tokens[6]), Integer.parseInt(tokens[7]),
					    Integer.parseInt(tokens[8]), Integer.parseInt(tokens[9]),
					    Integer.parseInt(tokens[10]), Integer.parseInt(tokens[11]))
			 );
	}
    }
    
    /* input is String representation of PBCount object. Note that this has one less column compared to the printBases format. Missing depth*/
    public PBCount(String line){
	String[] tokens = line.split(" ");
	nextId++;
	this.contig = new String(tokens[0]);
	this.pos = Integer.parseInt(tokens[1]);
	this.ref = tokens[2].charAt(0);
	acgtList = new ArrayList<ExtComposition>();
	if(tokens.length > 3){
	    for(int i=3; i<tokens.length;i=i+8){
		acgtList.add(new ExtComposition(Integer.parseInt(tokens[i]), Integer.parseInt(tokens[i+1]),
						Integer.parseInt(tokens[i+2]), Integer.parseInt(tokens[i+3]),
						Integer.parseInt(tokens[i+4]), Integer.parseInt(tokens[i+5]),
						Integer.parseInt(tokens[i+6]), Integer.parseInt(tokens[i+7]))
			     );
	    }
	}
    }

    public PBCount(String line, boolean forEasyString){
	String[] tokens = line.split(" ");
	nextId++;
	this.contig = new String(tokens[0]);
	this.pos = Integer.parseInt(tokens[1]);
	this.ref = tokens[2].charAt(0);
	acgtList = new ArrayList<ExtComposition>();
	for(int i=0; i<50;i++){
	    acgtList.add(null);
	}
	
	if(tokens.length > 3)
	    for(int i=3; i<tokens.length; i=i+11){
		int lineNum = Integer.parseInt(tokens[i+1]);
		acgtList.set(lineNum-1, new ExtComposition(Integer.parseInt(tokens[i+3]), Integer.parseInt(tokens[i+4]),
							  Integer.parseInt(tokens[i+5]), Integer.parseInt(tokens[i+6]),
							  Integer.parseInt(tokens[i+7]), Integer.parseInt(tokens[i+8]),
							  Integer.parseInt(tokens[i+9]), Integer.parseInt(tokens[i+10]))
			     );
			     
	    }
    }

    

    public void removeLines(int[] removal){
	for(int i=0;i<removal.length;i++){
	    acgtList.set(removal[i]-1, null);
	}
    }

    


    /* this is to parse alignPos files --> delimeted by spae and no depth field*/
    public PBCount(String[] tokens, boolean multi){
	if(!multi)
	    ;//this(tokens);
	else{
	    this.readerid = nextId;
	    nextId++;
	    this.contig = new String(tokens[0]);
	    this.pos = Integer.parseInt(tokens[1]);
	    this.ref = tokens[2].charAt(0);
	    acgtList = new ArrayList<ExtComposition>();
	    for(int i=3; i<tokens.length;i=i+8){
		acgtList.add(new ExtComposition(Integer.parseInt(tokens[i]), Integer.parseInt(tokens[i+1]),
						Integer.parseInt(tokens[i+2]), Integer.parseInt(tokens[i+3]),
						Integer.parseInt(tokens[i+4]), Integer.parseInt(tokens[i+5]),
						Integer.parseInt(tokens[i+6]), Integer.parseInt(tokens[i+7]))
			     );
	    }
	}
    }

    public ArrayList<ExtComposition> getACGTList(){
	return this.acgtList;
    }

    public int getReaderID(){
	return this.readerid;
    }

    public String getHeader(){
	return this.contig + "\t" + this.pos;
    }
    
    public String getACGTString(){
	int sumA = 0;
	int sumC = 0;
	int sumG = 0;
	int sumT = 0;
	for(int i=0; i<acgtList.size();i++){
	    ExtComposition tmp = acgtList.get(i);
	    sumA += (tmp.A + tmp.a);
	    sumC += (tmp.C + tmp.c);
	    sumG += (tmp.G + tmp.g);
	    sumT += (tmp.T + tmp.t);
	    
	}
	return "A=" + sumA + ";C=" + sumC + ";G=" + sumG + ";T=" + sumT + ";Total=" +(sumA + sumC + sumG + sumT) + ";";
    }
    
    public int getCountForA(){
	int sum = 0;
	for(int i=0; i<acgtList.size();i=i++){
	    sum += acgtList.get(i).A;
	    sum += acgtList.get(i).a;
	}
	return sum;
    }
    
    public int getCountForC(){
        int sum = 0;
        for(int i=0; i<acgtList.size();i=i++){
            sum += acgtList.get(i).C;
            sum+= acgtList.get(i).c;
        }
	return sum;
    }
    
    public int getCountForG(){
        int sum = 0;
        for(int i=0; i<acgtList.size();i=i++){
            sum += acgtList.get(i).G;
            sum+= acgtList.get(i).g;
        }
	return sum;
    }
    
    public int getCountForT(){
        int sum = 0;
        for(int i=0; i<acgtList.size();i=i++){
            sum += acgtList.get(i).T;
            sum+= acgtList.get(i).t;
        }
	return sum;
    }
    
    public boolean isEmpty(int lineNumber){
    	ExtComposition comp = acgtList.get(lineNumber);
    	return (comp.a + comp.A + comp.c + comp.C + comp.g + comp.G + comp.t + comp.T <= 5);
    }
    
    
    public String getGFFAttribute(){
	return "ID=" + this.contig + "_" + this.pos + ";Name=" + this.contig + "_" + this.pos + ";" + this.getACGTString(); 
    }
    
    public String toGFFSNPformat(char consensus, char modBase){
	StringBuffer sb = new StringBuffer(this.contig + "\t\t" +  "SNP\t" + this.pos + "\t" + this.pos + "\t.\t.\t.\t" + this.getGFFAttribute() + "Reference_Base=" + this.ref + ";Consensus=" + consensus + ";Modified_Call=" + modBase +";snp_type=Homozygous");
	return sb.toString();
    }
    
    public String getContig(){
	return this.contig;
    }

    public int getPos(){
	return this.pos;
    }
    
    public char getRef(){
	return this.ref;
    }
    

    

    public PBCount appendToACGT(ExtComposition c){
	this.acgtList.add(c);
	return this;
    }

    public PBCount merge(PBCount other){
	if(this.contig.equals(other.getContig()) && this.pos == other.getPos())
	    appendToACGT(other.getACGTList());
	return this;
    }

    public PBCount appendToACGT(ArrayList<ExtComposition> el){
	this.acgtList.addAll(el);
	return this;
    }

    /* 
     * length of printIndicies is the number of ExtComposition we are printing
     * only print the columns that are defined by printIndicies elements.
     *
     */
    public StringBuffer toStringBuffer(int[] printIndicies){
	StringBuffer buffer = new StringBuffer(this.contig + " " + pos + " " + ref);
	for(int i=0; i<printIndicies.length;i++){
	    buffer.append(" " + this.acgtList.get(printIndicies[i]).toString());
	}
	return buffer;
    }

    
    public StringBuffer toStringBuffer(){
	StringBuffer buffer = new StringBuffer(this.contig + " " + pos + " " + ref);
	for(int i=0; i<this.acgtList.size();i++){
	    buffer.append(" " + this.acgtList.get(i).toString());
	}
	return buffer;
    }

    public String toEasyString(){
	StringBuffer buffer = new StringBuffer(this.contig + " " + pos + " " + ref);
	for(int i=0; i<this.acgtList.size();i++){
	    buffer.append(" | " + (i+1) + " | "+ this.acgtList.get(i).toString());
	}
	buffer.append("\n");
	return buffer.toString();
    }

    public String toEasyStringSpecial(){
	StringBuffer buffer = new StringBuffer(this.contig + " " + pos + " " + ref);
	for(int i=0; i<this.acgtList.size();i++){
	    if(this.acgtList.get(i) != null)
		buffer.append(" | " + (i+1) + " | "+ this.acgtList.get(i).toString());
	}
	buffer.append("\n");
	return buffer.toString();
    }

    public String toSharedMutTableString(){
	StringBuffer buffer = new StringBuffer(this.contig + "\t" + pos + "\t" + ref);
	for(int i=0; i<this.acgtList.size();i++){
	    if(this.acgtList.get(i) != null)
		buffer.append("\t" + this.acgtList.get(i).getDominantAllele());
	    else
		buffer.append("\t-" );
	}
	buffer.append("\n");
	return buffer.toString();
    }

    public String toEasyString(boolean printEmpty){
	if(printEmpty)
	    return this.toEasyString();
	else{
	    if(this.acgtList.size() > 0){
		return this.toEasyString();
	    }
	    return "";
	}
    }

    public String toEasyStringByLine(){
	StringBuffer buffer = new StringBuffer("#" + this.contig + "\t" + pos + "\t" + ref+ "\n");
	for(int i=0; i<this.acgtList.size();i++){
	    buffer.append("| " + (i+1) + " |\t"+ this.acgtList.get(i).toStringTab() + "\n");
	}
	//buffer.append("\n");
	return buffer.toString();
    }

    public String toEasyString(char consensus){
	StringBuffer buffer = new StringBuffer(this.contig + " " + pos + " " + ref);
	boolean empty  = true;
	for(int i=0; i<this.acgtList.size();i++){
	    if(this.acgtList.get(i).isMutation(consensus, 0.8, 0) != 0){
		buffer.append(" | " + (i+1) + " | "+ this.acgtList.get(i).toString());
		empty = false;
	    }
	}
	buffer.append("\n");
	if(!empty)
	    return buffer.toString();
	else
	    return "";
    }


    public String getMutationSpectrum(char consensus){
	StringBuffer buffer = new StringBuffer(this.contig + " " + pos + " " + ref);
	for(int i=0; i<this.acgtList.size();i++){
	    if(this.acgtList.get(i).isMutation(consensus) != 0)
		buffer.append(" | " + (i+1)+" | " + this.acgtList.get(i).toString());
	}
	buffer.append("\n");
	return buffer.toString();
    }

    public String getMutationSpectrum(char consensus, boolean printEmpty){
	if(printEmpty)
	    return this.getMutationSpectrum(consensus);
	else{
	    if(this.acgtList.size() > 0)
		return this.getMutationSpectrum(consensus);
	    else
		return "";
	}
    }

}



class PosSorter implements Comparator<PBCount>{
    
    private Hashtable<String, Integer> g2num; /* the key is acronym for genome e.g.) ecoli for E. coli MG1655 */
    private Comparator<PBCount> posSorter;
    
    public PosSorter(String genomeHeaderFile){
	g2num = new Hashtable<String, Integer>();
	loadG2num(genomeHeaderFile);
    }

    public PosSorter(String[] headerFiles){
	g2num = new Hashtable<String, Integer>();
	loadG2num(headerFiles);
    }
    
    public int compare(PBCount s1, PBCount s2){
	if(s1 != null && s2 == null)
            return -1;
        else if(s1 == null && s2 !=null)
            return 1;
        else if(s1 == null && s2 == null)
            return 0;
	
        return this.compare(s1.getContig(), s1.getPos(), s2.getContig(), s2.getPos());
    }
    
    /* 
     *  returns 
     * -1 : if pos1 is smaller than pos2
     *  0 : if eq
     *  1 : if pos2 is smaller than pos1
     *
     */
    public int compare(String g1, int pos1, String g2, int pos2){
	if(g1.equals(g2)){
            if(pos1 < pos2)
                return -1;
            else if(pos1 == pos2)
                return 0;
            else
                return 1;
        }else
            return compareGenome(g1, g2);
    }
    
    private int compareGenome(String g1, String g2){
        if(map(g1) < map(g2))
            return -1;
        else
            return 1;
    }

    public int map(String g){
        if(this.g2num.get(g) == null){
            System.out.println("#" + g + "#");
	}
        return this.g2num.get(g).intValue();
    }

    private void loadG2num(String file){
        BufferedReader br = null;
        int index = 0;
        try{
            br = new BufferedReader(new FileReader(file));
            String curLine = "";
            while((curLine=br.readLine())!=null){
                String[] tokens = curLine.split("\\t");
                this.g2num.put(tokens[1].substring(tokens[1].indexOf(":")+1), new Integer(index));
                index++;
            }
	    br.close();
        }catch(IOException ioe){
            ioe.printStackTrace();
        }

    }

    private void loadG2num(String[] files){
	BufferedReader br = null;
	int index =0;
	try{
	    for(int i=0;i<files.length;i++){
		System.err.println("[PosSorter.loadG2num] READING IN-->" + files[i]);
		br = new BufferedReader(new FileReader(files[i]));
		String curLine = "";
		while((curLine=br.readLine())!=null){
		    System.err.println("[PosSorter.loadG2num]-->" + curLine);
		    String[] tokens = curLine.split("\\t");
		    this.g2num.put(tokens[1].substring(tokens[1].indexOf(":")+1), new Integer(index));
		    index++;
		}
		br.close();
		br = null;
	    }
	}catch(IOException ioe){
            ioe.printStackTrace();
        }
    }

}
