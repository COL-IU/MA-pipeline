import java.io.*;

public class ProcessPileup{

    private int polyCount = 0;
    private int pMutCount = 0;
    private int sameCount = 0;
    private StringBuffer simple = new StringBuffer();

    /* args[0] --> pileupFile name */
    public static void main(String[] args){
	//new ProcessPileup().parse(args[0]);
	new ProcessPileup().parseIndel(args[0]);
    }


    public void parseIndel(String pileupFile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(pileupFile));
	    String curLine = "";
	    while((curLine=br.readLine())!=null){
		PileupRecord tmp = new PileupRecord(curLine);
		if(tmp.isThereIndels()){
		    String tmptmp = "";
		    if( (tmptmp=tmp.isThereSignificantIndels()) != null){
			System.out.println(tmp.toString() + "\t" + tmptmp);
		    //this.count(curResult);
		    }
		}
	    }
	    br.close();
	    //this.outToFile(pileupFile + ".mut");
	    //System.out.println("Poly      :\t" + this.polyCount + "\tpointMutation :\t" + this.pMutCount
	    //		       + "\tsamCount    :\t" + this.sameCount  );
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void parse(String pileupFile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(pileupFile));
	    String curLine = "";
	    while((curLine=br.readLine())!=null){
		int[] curResult = new PileupRecord(curLine).run2(simple); // this returns A T G C count spectrum for a single position
		this.count(curResult);
	    }
	    br.close();
	    this.outToFile(pileupFile + ".mut");
	    System.out.println("Poly      :\t" + this.polyCount + "\tpointMutation :\t" + this.pMutCount
			       + "\tsamCount    :\t" + this.sameCount  );
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void count(int[] result){
	if(result[0] == 2)
	    polyCount++;
	else if(result[0] == 1)
	    pMutCount++;
	else if(result[0] == 0)
	    sameCount++;
    }

    public void outToFile(String f){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(f));
	    bw.write(this.simple.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
    }

}
