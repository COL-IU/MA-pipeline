import java.io.*;
import java.util.*;

public class CalculateLineByLineStat{
    
    private Hashtable<String, MutationStat> nameToStatHash;
    private ArrayList<String> lineNames;

    public static void main(String[] args){
	if(args.length != 2)
	    System.err.println("Usage: java CalculateLineByLineStat <index2SampleFile> <snptabl - .detail>");
	else
	    new CalculateLineByLineStat(args[0], args[1]).printResult();
    }
    
    public void printResult(){
	System.out.println("");
	for(int i=0;i<this.lineNames.size();i++){
	    this.nameToStatHash.get(this.lineNames.get(i)).printResult();
	}
    }

    public CalculateLineByLineStat(String index2SampleFile, String table){
	this.lineNames = new ArrayList<String>();
	this.loadHash(index2SampleFile);
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(table));
	    String curline = "";
	    while((curline = br.readLine()) !=null){
		String[] tokens = curline.split("\\t");
		if(!tokens[0].equals("Line")){
		    if(!this.nameToStatHash.get(tokens[1]).updateStat(tokens))
			System.err.println("HUH!!!:\t" + curline);
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    private void loadHash(String file){
	this.nameToStatHash = new Hashtable<String, MutationStat>();
	BufferedReader br = null;
	try{
	    String curline = "";
	    br = new BufferedReader(new FileReader(file));
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		this.nameToStatHash.put(tokens[1], new MutationStat(Integer.parseInt(tokens[0]), tokens[1]));
		this.lineNames.add(tokens[1]);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    
    

}


class MutationStat{

    private int lineIndex;
    private String lineName;
    private int syn;
    private int nsyn;
    private int nc;
    private int[] spectrum;
    private int con;
    private int ncon;
    private int ts;
    private int tv;
    private int total;

    public MutationStat(int li, String ln){
	this.lineIndex = li;
	this.lineName = ln;
	this.syn = 0;
	this.nsyn = 0;
	this.nc = 0;
	this.spectrum = new int[6];
	this.con = 0;
	this.ncon = 0;
	this.ts = 0;
	this.tv = 0;
    }

    public void printResult(){
	double percent_syn = (syn*1.0d) / (total*1.0d);
	double percent_nsyn = (nsyn*1.0d) / (total*1.0d);
	double percent_nc = (nc*1.0d) / (total*1.0d);
	
	double[] percent_spectrum = new double[this.spectrum.length];
	double[] percent_spectrum_within = new double[this.spectrum.length];
	for(int i=0; i<percent_spectrum.length;i++){
	    percent_spectrum[i] = (this.spectrum[i]*1.0d) / (total*1.0d);
	    if(i<2)
		percent_spectrum_within[i] = (this.spectrum[i]*1.0d) / (ts*1.0d);
	    else
		percent_spectrum_within[i] = (this.spectrum[i]*1.0d) / (tv*1.0d);
	}
	
	double percent_con = (con*1.0d) / (total*1.0d);
	double percent_ncon = (ncon*1.0d) / (total*1.0d);
	
	double percent_ts = (ts*1.0d) / (total*1.0d);
	double percent_tv = (tv*1.0d) / (total*1.0d);
	

	double nsyn_syn_ratio = (nsyn*1.0d) / (syn*1.0d);
	double ncon_con_ratio = (ncon*1.0d) / (con*1.0d);

	//System.out.println(syn + " " + nsyn);
	System.out.print(this.lineName + "\t" + this.total + "\t" + percent_syn + "\t" + percent_nsyn + "\t" + percent_nc + "\t");
	for(int i=0;i<percent_spectrum.length;i++){
	    System.out.print(percent_spectrum[i] + "\t");
	}
	for(int i=0;i<percent_spectrum_within.length;i++){
	    System.out.print(percent_spectrum_within[i] + "\t");
	}
	System.out.println(percent_con + "\t" + percent_ncon + "\t" + percent_ts + "\t" + percent_tv + "\t" + nsyn_syn_ratio + "\t" + ncon_con_ratio);
    }

    public boolean updateStat(String[] tokens){
	if(Integer.parseInt(tokens[0]) == this.lineIndex && tokens[1].equals(lineName)){
	    this.total++;
	    if(tokens[13].equals("Transition"))
		this.ts++;
	    else if(tokens[13].equals("Transversion"))
		this.tv++;
	    
	    
	    handleSpectrum(tokens[14]);
	    
	    if(tokens[15].equals("Conservative"))
		this.con++;
	    else if(tokens[15].equals("Non-conservative"))
		this.ncon++;
	    
	    if(tokens[16].equals("Syn")){
		this.syn++;
		//System.err.println(" |HERERERERER| ");
	    }else if(tokens[16].equals("N-Syn"))
		this.nsyn++;
	    else if(tokens[16].equals("NC"))
		this.nc++;
	    
	    return true;
	}else
	    return false;
    }

    private void handleSpectrum(String spectrum){
	if(spectrum.equals("AT>GC"))
	    this.spectrum[0]++;
	else if(spectrum.equals("GC>AT"))
	    this.spectrum[1]++;
	else if(spectrum.equals("AT>TA"))
	    this.spectrum[2]++;
	else if(spectrum.equals("GC>TA"))
	    this.spectrum[3]++;
	else if(spectrum.equals("AT>CG"))
	    this.spectrum[4]++;
	else if(spectrum.equals("GC>CG"))
	    this.spectrum[5]++;
    }
    
}
