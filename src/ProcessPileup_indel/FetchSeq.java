import java.util.*;
import java.io.*;
public class FetchSeq{
    
    private Hashtable<String, Seq> seqHash;
    
    public FetchSeq(String fasta){
	this.seqHash = new FastaReader().parseFasta(fasta);
    }
	
    public static void main(String[] args){
	new FetchSeq(args[0]).process(args[1], Integer.parseInt(args[2]));
    }
    
    public void process(String file, int length){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(file));
	    String curline = "";
	    String curSeq = "";
	    boolean rev = false;
	    while((curline=br.readLine())!=null){
		if(curline.startsWith("#")){
		    if(curline.equals("#fwd"))
			rev = false;
		    else if(curline.equals("#rev"))
			rev = true;
		    else
			curSeq = curline.substring(1);
		}else{
		    int index = Integer.parseInt(curline);
		    if(!rev)
			System.out.println(seqHash.get(curSeq).getSubSeq(index,index+length));
		    else
			System.out.println(seqHash.get(curSeq).getSubSeq(index-length+1,index+1));
					   
		}
		
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

}
