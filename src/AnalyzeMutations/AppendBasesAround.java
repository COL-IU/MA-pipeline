import java.io.*;
import java.util.*;

public class AppendBasesAround{
    
    private Hashtable<String, Seq> seqHash;
    private int seqLen = 20; //length of surrounding sequence. default = 20;

    public static void main(String[] args){
	if(args.length == 0){
	    new AppendBasesAround("/data/groups/heewlee/MA/references/K12MG1655/K12MG1655.fna").test();
	    
	}

	if(args.length != 2)
	    System.err.println("USAGE: java AppendBasesAround <fastaFile> <snpTable>");
	else
	    new AppendBasesAround(args[0]).run(args[1]);
    }

    public void test(){
	for(int i=4639655;i<=4639675;i++){
	    System.out.println(this.toTabularString(this.seqLen, "K12MG1655_Chr1", i));
	}
    }

    public AppendBasesAround(){
	this.seqLen = 20;
    }
    
    public AppendBasesAround(String fastaFile){
	this(fastaFile, 20);
    }
    
    public AppendBasesAround(String fastaFile, int len){
	this.loadHash(fastaFile);
	this.seqLen = len;
    }

    public AppendBasesAround(String detailFile, Hashtable<String, Seq> hash, BufferedWriter bw){
	this.seqHash = hash;
	//this.run(detailFile);
	this.processTable(detailFile, bw);
    }

    public void run(String snpTable){
	this.processTable(snpTable);
    }

    private void loadHash(String fastaFile){
	this.seqHash =  new FastaReader().parseFasta(fastaFile);
    }
    
    public void processTable(String table){
	BufferedReader br = null;
	String curline = "";
	try{
	    br = new BufferedReader(new FileReader(table));
	    while((curline=br.readLine())!=null){
		if(curline.startsWith("Line")||curline.startsWith("Index")){
		    System.out.println(curline);
		}else{
		    String[] tokens = curline.split("\\t");
		    String contigName = tokens[2];
		    int pos = Integer.parseInt(tokens[3]);
		    System.out.println(curline + "\t" + this.toTabularString(this.seqLen, contigName, pos));
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public void processTable(String table, BufferedWriter bw){
	BufferedReader br = null;
	String curline = "";
	try{
	    br = new BufferedReader(new FileReader(table));
	    while((curline=br.readLine())!=null){
		if(curline.startsWith("Line")||curline.startsWith("Index")){
		    bw.write(curline+"\n");//System.out.println(curline);
		}else{
		    String[] tokens = curline.split("\\t");
		    String contigName = tokens[2];
		    int pos = Integer.parseInt(tokens[3]);
		    bw.write(curline + "\t" + this.toTabularString(this.seqLen, contigName, pos) + "\n");//System.out.println(curline + "\t" + this.toTabularString(20, contigName, pos));
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private String toTabularString(int n, String contigName, int pos){
	return this.get5PrimeSeq(n, this.seqHash.get(contigName), pos) + "\t"
	    + this.get3PrimeSeq(n, this.seqHash.get(contigName), pos);
    }

    /*    private String get5PrimeSeq(int n, Seq seq, int pos){
        return seq.getSubSeqWithDelim(pos+1-n,pos+1, "\t");
	}*/

    private String get5PrimeSeq(int n, Seq seq, int pos){
	//if(pos == ){
	//   return seq.getSubSeqWithDelim(seq.getSeqLength()+1-n , seq.getSeqLength()+1, "\t");
	//}else 

	/*	if(pos+1-n < 1)
	    return seq.getSubSeqWithDelim(seq.getSeqLength() - (1-(pos+1-n)-1) , seq.getSeqLength()+1 , "\t") +"\t"+ seq.getSubSeqWithDelim(1,pos+1,"\t");
	else
	    return seq.getSubSeqWithDelim(pos+1-n,pos+1, "\t");
	*/
	if(pos<n)
	    return seq.getSubSeqWithDelim(seq.getSeqLength() - (n-pos) , seq.getSeqLength()+1 , "\t") +"\t"+ seq.getSubSeqWithDelim(1,pos+1,"\t");
	else
	    return seq.getSubSeqWithDelim(pos+1-n,pos+1, "\t");
    }
   
    private String get3PrimeSeq(int n, Seq seq, int pos){
	if(pos == seq.getSeqLength()){
	    return seq.getSubSeqWithDelim(1,1+n,"\t");
	}else if(pos+n > seq.getSeqLength())
	    return seq.getSubSeqWithDelim(pos+1,seq.getSeqLength()+1, "\t") + "\t" + seq.getSubSeqWithDelim(1,pos+n-seq.getSeqLength()+1, "\t");
	else
	    return seq.getSubSeqWithDelim(pos+1,pos+1+n, "\t");
    }

}
