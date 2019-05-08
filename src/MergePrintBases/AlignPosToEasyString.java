import java.io.*;

public class AlignPosToEasyString{

    public static void main(String[] args){
	convert(args[0]);
    }
    
    public static void convert(String file){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(file));
	    String curline = "";
	    try{
		while((curline=br.readLine())!=null){
		    PBCount tmp = new PBCount(curline);
		    //System.out.println(tmp.toEasyString());
		    System.out.println(tmp.toEasyStringByLine());
		}
		br.close();
	    }catch(IOException ioe){
		ioe.printStackTrace();
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

}
