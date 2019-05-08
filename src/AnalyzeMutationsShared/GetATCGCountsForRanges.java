import java.io.*;
public class GetATCGCountsForRanges{

    public static void main(String[] args){
	new GetATCGCountsForRanges().run(args[0], args[1]);
    }
    
    public void run(String ranges, String genomeFasta){
	
	String s = "ecoli";
	Seq seq = new FastaReader().parseFasta(genomeFasta).get(s);

	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(ranges));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		this.processEach(Integer.parseInt(tokens[0]), Integer.parseInt(tokens[1]), seq);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public void processEach(int s, int e, Seq seq){
	System.out.print(s + "\t" + e);
	seq.printATCGCounts(s, e);
    }
    
}
