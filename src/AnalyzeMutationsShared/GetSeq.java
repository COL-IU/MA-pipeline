import java.util.*;

public class GetSeq{
    
    public static void main(String[] args){
	if(args.length == 5)
	    getSeq(args[0],args[1],Integer.parseInt(args[2]), Integer.parseInt(args[3]), ((args[4].toUpperCase().charAt(0)=='N')?false:true));
	else
	    getSeq(args[0],args[1],Integer.parseInt(args[2]), Integer.parseInt(args[3]));
    }
    
    /* start and end inclusive */
    public static void getSeq(String fasta, String contig, int start, int end, boolean printHeader){
	Hashtable<String, Seq> seqHash = new FastaReader().parseFasta(fasta);
	Seq seq = seqHash.get(contig);
	if(seq !=null){
	    if(printHeader)
		System.out.println(">" + contig + "_"+start+"_"+end);
	    System.out.println(seq.getSubSeq(start, end+1));
	}
    }
    
    public static void getSeq(String fasta, String contig, int start, int end){
	getSeq(fasta,contig,start,end,true);
    }
    
    
    
    
    
}
