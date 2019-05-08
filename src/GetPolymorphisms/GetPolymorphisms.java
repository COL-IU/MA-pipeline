import java.io.*;

public class GetPolymorphisms{
    
    public static void main(String[] args){
	if(args.length == 2)
	    getPolymorphisms(args[0], args[1]);
	else
	    System.err.println("USAGE: java GetPolymorphisms <.consensus file> <outFile>");
    }

    public static void getPolymorphisms(String consensusFile, String outFile){
	BufferedReader br = null;
	String curline = "";
	BufferedWriter bw = null;
	
	try{
	    br = new BufferedReader(new FileReader(consensusFile));
	    bw = new BufferedWriter(new FileWriter(outFile));
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\s");
		Consensus c = new Consensus(tokens);
		if(c.isPolymorphism()){
		    //System.err.println(c.refConsString());
		    System.out.println(c.toString());
		    bw.write(c.fakePolyString(5)+"\n");
		}
	    }
	    br.close();
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

}
