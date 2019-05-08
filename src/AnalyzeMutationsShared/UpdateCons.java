import java.io.*;

public class UpdateCons{

    private Codons codons;
    
    public UpdateCons(String codonfile, String smFile){
	this.codons = new Codons(codonfile, smFile);
    }
    

    public static void main(String[] args){
	UpdateCons obj = new UpdateCons(args[0], args[1]);
	obj.process(args[2]);
    }

    public void process(String table){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(table));
	    String curline = "";
	    while( (curline=br.readLine()) !=null){
		String[] tokens = curline.split("\\t");
		if(tokens[9].charAt(0) == '-' || tokens[9].equals(tokens[10]))
		    System.out.println("-");
		else if(this.codons.isConservative(codons.getSM(), 0, true, tokens[9].charAt(0), tokens[10].charAt(0)))
		    System.out.println("Conservative");
		else
		    System.out.println("Non-conservative");
		
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
}
