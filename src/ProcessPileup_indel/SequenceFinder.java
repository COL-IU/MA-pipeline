import java.util.*;
import java.io.*;

public class SequenceFinder{
    
    private Hashtable<String, Seq> seqHash;

    public static void main(String[] args){
	if(args.length == 3)
	    new SequenceFinder(args[0]).process(args[1], args[2]);
	else
	    System.err.println("USAGE: java SequenceFinder <Reference Fasta> <Pattern> <'Y' to look both strands>");
    }

    public SequenceFinder(String fasta){
	this.seqHash = new FastaReader().parseFasta(fasta);
    }

    public void process(String motif, String strandOption){
	
	boolean bothStrands = (strandOption.equals("Y") ? true : false);
	this.printLoci(this.getLoci(motif, bothStrands));
    }

    public Hashtable<String, ArrayList<Integer>[]> getLoci(String motif, boolean bothStrands){
	Enumeration<String> seqNames = this.seqHash.keys();
	Hashtable<String, ArrayList<Integer>[]> motifLociHash = new Hashtable<String, ArrayList<Integer>[]>();
	while(seqNames.hasMoreElements()){
	    String tmpKey = seqNames.nextElement();
	    motifLociHash.put(tmpKey, this.seqHash.get(tmpKey).getMotifLoci(motif, bothStrands));
	}
	return motifLociHash;
    }

    private void printLoci(Hashtable<String, ArrayList<Integer>[]> hash){
	StringBuffer bf = new StringBuffer();
	Enumeration<String> loci = hash.keys();
	while(loci.hasMoreElements()){
	    String tmp = loci.nextElement();
	    ArrayList<Integer>[] lists = hash.get(tmp);
	    
	    for(int j=0;j<lists.length;j++){
		bf.append("#" + tmp + "\n");
		if(j==0){
		    bf.append("#fwd\n");
		}else
		    bf.append("#rev\n");
		for(int i=0;i<lists[j].size();i++){
		    bf.append(lists[j].get(i) + "\n");
		}
	    }
	}
	System.out.print(bf.toString());
    }

}
