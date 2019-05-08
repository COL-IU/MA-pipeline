import java.io.*;
import java.util.*;

public class Num2Sample{

    private Hashtable<Integer, String> num2sample;
    private String headers;

    public Num2Sample(String mapFile){
	this.loadHash(mapFile);
    }

    public String mapNum2Sample(int n){
	return this.num2sample.get(new Integer(n));
    }

    public String getHeaders(){
	return this.headers;
    }

    private void loadHash(String mapFile){
	StringBuffer tmp = new StringBuffer();
	this.num2sample = new Hashtable<Integer, String>();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(mapFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		this.num2sample.put(new Integer(tokens[0]), tokens[1]);
		tmp.append("\t" + tokens[1]);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	this.headers = tmp.toString();
    }
}
