import java.io.*;
import java.util.*;

public class FastaReader{
    
    //    private Hashtable<String, Seq> contigHash;
    
    public Hashtable<String, Seq> parseFasta(String fileName){
	BufferedReader br = null;
	Hashtable<String, Seq> contigHash = new Hashtable<String, Seq>();
	try{
	    br = new BufferedReader(new FileReader(fileName));
	    String curLine = "";
	    boolean inTarget = false;
	    Seq curSeq = null;
	    while((curLine=br.readLine())!=null){
		curLine = curLine.trim();
		
		if(curLine.length() == 0){
		    ;
		}else {
		    if( curLine.charAt(0) == '>'){
			if(curSeq != null){
			    //System.out.println("HERE :" + curSeq.getSeqName());
			    contigHash.put(curSeq.getSeqName(), curSeq);
			    curSeq = null;
			}
			curSeq = new Seq(curLine.substring(1));
		    }else
			curSeq.appendSeq(curLine);
		}
	    }
	    contigHash.put(curSeq.getSeqName(), curSeq);
	    br.close();
	    br = null;
	}catch(IOException ioe){
	    ioe.printStackTrace();
            return null;
	}
	return contigHash;
    }

    /*public Hashtable<String, Seq> getContigHash(){
	return this.contigHash;
	}*/

}

