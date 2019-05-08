import java.io.*;

public class ExtractAlignPosFromConsensus{

    public static void main(String[] args){
	if(args.length == 2)
	    new ExtractAlignPosFromConsensus().extract(args[0], args[1]);
	else
	    System.err.println("USAGE: java ExtractAlignPosFromConsensus <.consensus file> <.alignPos file>");
    
    }

    public void extract(String consensus, String alignpos){
	BufferedReader brc = null;
	BufferedReader bra = null;
	BufferedWriter bw = null;
	StringBuffer output = new StringBuffer();
	StringBuffer apOut = new StringBuffer();//.alignPos formatted output
	try{
	    bra = new BufferedReader(new FileReader(alignpos));
	    brc = new BufferedReader(new FileReader(consensus));
	    String consenLine = brc.readLine();
	    String curkey = this.getKey(consenLine);
	    char curConsen = consenLine.charAt(consenLine.length()-1);
	    String curline = "";
	    while((curline=bra.readLine())!=null){
		if(curline.startsWith(curkey)){
		    PBCount pbc = new PBCount(curline);
		    apOut.append(pbc.toEasyString(curConsen));
		    output.append(pbc.getMutationSpectrum(curConsen));
		    consenLine = brc.readLine();
		    if(consenLine == null)
			break;
		    curkey = this.getKey(consenLine);
		    curConsen = consenLine.charAt(consenLine.length()-1);
		    
		}
	    }
	    bra.close();
	    brc.close();
	    bw = new BufferedWriter(new FileWriter(consensus + ".EalignPos"));
	    bw.write(apOut.toString());
	    bw.close();
	    System.out.println(output.toString());
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }



    public String getKey(String cLine){
	if(cLine == null)
	    return null;
	String[] tokens = cLine.split(" ");
	return tokens[0] + " " + tokens[1] + " ";
    }
}
