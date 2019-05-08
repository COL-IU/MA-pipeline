import java.io.*;
import java.util.*;

public class PrintConsensus{

    //private int[] coverage;
    private int[] minNumReadsEachDirection;
    
    // [10/17/13] by heewlee ADDED to flag off the option to use different criteria for diff coverage threshold vals.
    private boolean accomodateVaryingCoverages; // default value is true: where we use 10-10 rule for lines with 100x or higher coverages and 3-3 rule otherwise.
    
    private static final int HIGH_COV_MIN_READS_EACH_DIRECTION = 10;
    private static final int LOW_COV_MIN_READS_EACH_DIRECTION = 3;
    
    public static void main(String[] args){
	if(args.length == 5)
	    new PrintConsensus(args[2], args[3], args[4]).process(args[0],args[1]);
	else if(args.length == 6)
	    new PrintConsensus(args[2], args[3], args[4], ( (args[5].equals("F") || args[5].equals("f")) ? false : true) ).process(args[0], args[1]);
	else
	    System.err.println("USAGE: java PrintConsensus <.alignPos> <outFileName .consensus file> <lineListFile> <outDir> <MAIN_CHROMOSOM_NAME> [Optional: accomodateVaryingCoverages Tt/Ff  :default True]");
    }
    
    public PrintConsensus(String listFile, String outDir, String mainChromo, boolean accomodateVaryingCoverages){
	this.accomodateVaryingCoverages = accomodateVaryingCoverages;
	this.loadCovAndMinNum(listFile, outDir, mainChromo);
    }

    public PrintConsensus(String listFile, String outDir, String mainChromo){
	//this.setCoverage(stat2File);
	//this.coverage = 90;
	//this.setMinNumReads();
	this(listFile,outDir, mainChromo, true);
    }

    public void loadCovAndMinNum(String listFile, String outDir, String mainChromo){
	
	BufferedReader br = null;
	ArrayList<String> lines = new ArrayList<String>();
	String genome = mainChromo;
	try{
	    br = new BufferedReader(new FileReader(listFile));
	    String curline ="";
	    while((curline=br.readLine())!=null){
		lines.add(outDir + File.separator + curline + File.separator + curline +"_"+genome+"_stat2.txt");
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
	//this.coverage = new int[lines.size()];
	this.minNumReadsEachDirection = new int[lines.size()];
	
	for(int i=0;i<lines.size();i++){
	    int curCov = this.getCoverage(lines.get(i));
	    this.minNumReadsEachDirection[i] = this.computeMinNumReads(curCov);
	}
    
    
    }

    private int getCoverage(String stat2File){
	BufferedReader br = null;
	String curline = "";
	int cov=0;
	try{
	    br = new BufferedReader(new FileReader(stat2File));
	    br.readLine();
	    br.readLine();
	    curline = br.readLine(); // read the third line
	    String[] tokens = curline.split("\\t");
	    cov = (int) Double.parseDouble(tokens[17]);
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return cov;
    }
    
    private int computeMinNumReads(int cov){
	if(!this.accomodateVaryingCoverages)
	    return HIGH_COV_MIN_READS_EACH_DIRECTION;
	else{
	    if(cov >= 100){
		return HIGH_COV_MIN_READS_EACH_DIRECTION;
	    }else{
		return LOW_COV_MIN_READS_EACH_DIRECTION;
	    }
	}
    }
    
    public void process(String alignPosFile, String consensusOutFile){
	BufferedWriter bw = null;
	BufferedReader br = null;
	StringBuffer bf = new StringBuffer();
	try{
	    br = new BufferedReader(new FileReader(alignPosFile));
	    bw = new BufferedWriter(new FileWriter(consensusOutFile));
	    String curline = "";
	    AlignPos curAlignPos = null;
	    int count = 0;
	    while((curline=br.readLine())!=null){
		count++;
		if(count%10000 == 0){
		    bw.write(bf.toString());
		    bf = new StringBuffer();
		}
		curAlignPos = new AlignPos(curline, " ");//tokenizes baseCounts and prints header
		bf.append(curAlignPos.getHeader());
		int[] bases = curAlignPos.getBases();
		int allTotal = 0;
		int Acons = 0;
		int Ccons = 0;
		int Gcons = 0;
		int Tcons = 0;
		int nthLine = 0;
		for(int i=0;i<bases.length;i=i+8){
		    //System.out.println("BASES AS PARSED :\t" +bases[i] + " " + bases[i+1] +" " + bases[i+2] +" " + bases[i+3] +" " + bases[i+4] +" " +bases[i+5] +" " +bases[i+6] +" " +bases[i+7]);
		    int As = bases[i]+bases[i+1];
		    int Cs = bases[i+2]+bases[i+3];
		    int Gs = bases[i+4]+bases[i+5];
		    int Ts = bases[i+6]+bases[i+7];
		    //System.out.print("#AsCsGsTs:\t" + As + " " +Cs + " " +Gs+ " " +Ts);
		    int curTotal = As + Cs + Gs + Ts;
		    //System.out.print("\t#curTotal:\t" + curTotal );
		    Acons += As;
		    Ccons += Cs;
		    Gcons += Gs;
		    Tcons += Ts;
		    //System.out.print("\t#ALLINES@A|C|G|T:\t" + Acons + " " +Ccons + " " +Gcons+ " " +Tcons);
		    allTotal += curTotal;
		    
		    if(curTotal == 0)
			bf.append(" -");
		    else{
			int cutoff = (int)(curTotal*0.8d);
			//System.out.println("CUTOFF:\t" +cutoff + "#MINNUMREADS:\t" + this.minNumReadsEachDirection[nthLine]);
			if(As>cutoff){
			    if(bases[i]>=this.minNumReadsEachDirection[nthLine] && bases[i+1]>=this.minNumReadsEachDirection[nthLine])
				bf.append(" A");
			    else
				bf.append(" -");
			}else if(Cs>cutoff){
			    if(bases[i+2]>=this.minNumReadsEachDirection[nthLine] && bases[i+3]>=this.minNumReadsEachDirection[nthLine])
				bf.append(" C");
			    else
				bf.append(" -");
			}else if(Gs>cutoff){
			    if(bases[i+4]>=this.minNumReadsEachDirection[nthLine] && bases[i+5]>=this.minNumReadsEachDirection[nthLine])
				bf.append(" G");
			    else
				bf.append(" -");
			}else if(Ts>cutoff){
			    if(bases[i+6]>=this.minNumReadsEachDirection[nthLine] && bases[i+7]>=this.minNumReadsEachDirection[nthLine])
				bf.append(" T");
			    else
				bf.append(" -");
			}else{
			    bf.append(" -");
			}
		    }
		    nthLine++;
		}
		//System.out.print(allTotal + "\t" + Acons + "\t" + Ccons+ "\t" + Gcons+ "\t" + Tcons );
		if(allTotal == 0)
		    bf.append(" -");
		else{
		    int cutoff = allTotal/2;
		    //System.out.println("\t" + cutoff);
		    if(Acons > cutoff)
			bf.append(" A");
		    else if(Ccons > cutoff)
			bf.append(" C");
		    else if(Gcons > cutoff)
			bf.append(" G");
		    else if(Tcons > cutoff)
			bf.append(" T");
		    else
			bf.append(" -");
		}
		bf.append("\n");
	    }
	    br.close();
	    bw.write(bf.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

}

class AlignPos{

    /*private String scaffold;
    private int position;
    private char refBase;*/
    private String header;
    private int[] bases;
        
    AlignPos(String curline, String delim){
	String[] tokens = curline.split(delim);
	/*this.scaffold = tokens[0];
	this.position = Integer.parseInt(tokens[1]);
	this.refBase = tokens[2].charAt(0);*/
	this.header = tokens[0] + " " + tokens[1] + " " + tokens[2];
	this.loadBases(tokens);
    }
     
    public String getHeader(){
	return this.header;
    }
    
    public int[] getBases(){
	return this.bases;
    }

    private void loadBases(String[] tokens){
	this.bases = new int[tokens.length-3];
	for(int i=3;i<tokens.length;i=i+8){
	    bases[i-3] = Integer.parseInt(tokens[i]);
	    bases[i-2] = Integer.parseInt(tokens[i+1]);
	    bases[i-1] = Integer.parseInt(tokens[i+2]);
	    bases[i] = Integer.parseInt(tokens[i+3]);
	    bases[i+1] = Integer.parseInt(tokens[i+4]);
	    bases[i+2] = Integer.parseInt(tokens[i+5]);
	    bases[i+3] = Integer.parseInt(tokens[i+6]);
	    bases[i+4] = Integer.parseInt(tokens[i+7]);
	}
    }

    
}
