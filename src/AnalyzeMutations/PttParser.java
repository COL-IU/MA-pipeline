import java.io.*;
import java.util.*;

public class PttParser{
    
    /* key : contigname, value is genes in that contig */
    public Hashtable<String, ArrayList<Gene>>parseAll(String[] ptts){
	
	Hashtable<String, ArrayList<Gene>> geneHash = new Hashtable<String, ArrayList<Gene>>();
	for(int i=0;i<ptts.length;i++){
	    String key;
	    if(ptts[i].contains(File.separator))
		key = ptts[i].substring(ptts[i].lastIndexOf(File.separator)+1,ptts[i].indexOf(".ptt"));
	    else
		key = ptts[i].substring(0,ptts[i].indexOf(".ptt"));
	    //System.out.println("really inserting " + key);
	    geneHash.put(key, this.parse(ptts[i],key));
	}
	return geneHash;
    }
    
    public Hashtable<String, Gene> parseAndReturnLocusHash(String ptt, String ct){
	Hashtable<String, Gene> locus2Gene = new Hashtable<String, Gene>();
	BufferedReader br = null;
	String curline = "";
	try{
	    br = new BufferedReader(new FileReader(ptt));
	    boolean inFlag = false;
	    while((curline=br.readLine())!=null){
		if(!inFlag){
		    if(curline.startsWith("Location"))
			inFlag = true;
		}else{
		    String[] tokens = curline.split("\\t");
		    int[] stend = getStartEnd(tokens[0]);
		    locus2Gene.put(tokens[5], new Gene(ct, stend[0], stend[1], (tokens[1].equals("+") ? true : false) , tokens[4]));
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return locus2Gene;
    }

    public ArrayList<Gene> parse(String ptt, String ct){
	ArrayList<Gene> list = new ArrayList<Gene>();
	BufferedReader br = null;
	String curline = "";
	try{
	    br = new BufferedReader(new FileReader(ptt));
	    boolean inFlag = false;
	    while((curline=br.readLine())!=null){
		if(!inFlag){
		    if(curline.startsWith("Location"))
			inFlag = true;
		}else{
		    String[] tokens = curline.split("\\t");
		    int[] stend = getStartEnd(tokens[0]);
		    list.add(new Gene(ct, stend[0], stend[1], (tokens[1].equals("+") ? true : false) , tokens[4]) );
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return list;
    }

    private int[] getStartEnd(String pttIndex){
	int[] stend = new int[2];
	stend[0] = Integer.parseInt(pttIndex.substring(0,pttIndex.indexOf("..")));
	stend[1] = Integer.parseInt(pttIndex.substring(pttIndex.indexOf("..")+2));
	return stend;
    }
}
