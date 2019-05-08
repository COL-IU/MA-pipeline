import java.io.*;

public class ConvertToGFF3{

    public static void main(String[] args){
	if(args.length == 2)
	    new ConvertToGFF3(args[0], args[1]);
	else
	    System.err.println("java ConvertToGFF3 <.putation file> <.putation.alignPos file>");
    }
    
    /*  putFile only contains the putation information and 
     *  alignPosFile is a truncated file that contains the same entry with putFile but in baseSpectrum array 
     */
    public ConvertToGFF3(String putFile, String alignPosFile){
	System.out.print(this.process(putFile, alignPosFile));
    }

    public String process(String putFile, String alignPosFile){
	BufferedReader pr = null;
	BufferedReader ar = null;
	String curlineP;
	String curlineA;
	StringBuffer bf = new StringBuffer();
	try{
	    pr = new BufferedReader(new FileReader(putFile));
	    ar = new BufferedReader(new FileReader(alignPosFile));
	    while((curlineP=pr.readLine()) != null 
		  && (curlineA=ar.readLine()) != null){
		Consensus consensus = new Consensus(curlineP.split(" "));
		PBCount pbCount = new PBCount(curlineA.split(" ") , true); /* we are loading multi bits*/
		if(consensus.getHeader().equals(pbCount.getHeader())){//looking at same entry in two different files
		    char mutBase = consensus.getMutBase();
		    char conBase = consensus.getConsensus();
		    bf.append(pbCount.toGFFSNPformat(conBase, mutBase)+"\n");
		}
	    }
	    pr.close();
	    ar.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return bf.toString();
    }

}
