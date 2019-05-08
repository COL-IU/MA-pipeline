import java.util.*;

public class Count{

    public static void main(String[] args){
	run(args[0]);
    }
    
    public static void run(String fasta){
	ArrayList<Seq> list = new FastaReader().parseFastaToList(fasta);
	for(int i=0;i<list.size();i++){
	    System.out.println(list.get(i).getSeqName() + "\tLen: " + list.get(i).getSeqLength() + "\tAT: " +  list.get(i).countAT() + "\tGC: " + list.get(i).countGC() + "\tN: " + list.get(i).countN());
	}
    }

}
