import java.io.*;
import java.util.*;

public class PrintSharedMutationTable{

    private Hashtable<Integer, String> lineHash;
    
    public void loadHash(String file){
	this.lineHash = new Hashtable<Integer, String>();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(file));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		lineHash.put(new Integer(tokens[0]), tokens[1]);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
    }

    public PrintSharedMutationTable(String lineFile){
	this.loadHash(lineFile);
    }

    public static void main(String[] args){
	new PrintSharedMutationTable(args[1]).run(args[0]);
    }

    public void run(String inFile){
	System.out.print("contig\tposition\trefBase");
	for(int i=1;i<51;i++){
	    System.out.print("\t" + this.lineHash.get(new Integer(i)));
	}
	System.out.println();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(inFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		System.out.print(new PBCount(curline, true).toSharedMutTableString());
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
}
