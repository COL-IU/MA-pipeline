import java.io.*;
import java.util.*;

public class CalculateGenerations{

    //private String parameterFile; //this is a 2-column, tab-delimited file where 1st column is colony diameter size and the 2nd column is the corresponding #generations
    //private String tableFile; // this is the table where each row is for a MA line with diameter sizes for all passages.
    private Hashtable<String, Double> size2generations;
    //min max diameter
    private double minD = 200000.0d;
    private double maxD = -200000.0d;
    
    public CalculateGenerations(String parameterFile, String table){
	this.size2generations = new Hashtable<String, Double>();
	this.loadParameters(parameterFile);
	this.process(table);
    }
    
    public static void main(String[] args){
	if(args.length == 2)
	    new CalculateGenerations(args[0], args[1]);//args[0] --> paprameterFile,  args[1] --> DiameterTableFile 
	else
	    System.err.println("java CalculateGenerations <parameterFile> <diameterTableFile>");
    }
    
    public void exitAndPrintError(String error){
	System.err.println("ERROR: " + error);
	System.exit(1);
    }
    
    public void loadParameters(String parameterFile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(parameterFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		double diameter = Double.parseDouble(tokens[0]);
		if(size2generations.get(tokens[0])!=null){
		    this.exitAndPrintError("Duplicate diameter entries in the parameter file");
		}
		if(!checkDFormatSimple(diameter))
		    this.exitAndPrintError("Diameter not a multiple of 0.5 [ " + diameter + " ]");
		this.updateMinMax(diameter);
		this.size2generations.put(tokens[0], Double.valueOf(tokens[1]));
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private void updateMinMax(double tmp){
	if(tmp < minD)
	    this.minD = tmp;
	else if(tmp > maxD)
	    this.maxD = tmp;
    }
    
    public boolean checkDFormatSimple(double val){
	if(val%0.5 != 0)
	    return false;
	return true;
    }

    public boolean checkDFormat(double val){
	System.err.println(val);
	if(val%0.5 != 0) //multiple of 0.5
	    return false;
	if(val<minD || val > maxD){ //within minMax
	    System.err.println("minDiameter : [ " + minD + " ]\t maxDiameter : [ " + maxD + " ]" );
	    return false;
	}
	return true;
    }
    
    public boolean checkDFormat(String diameter){
	return this.checkDFormatSimple(Double.parseDouble(diameter));
    }
    
    public void process(String table){
	BufferedReader br = null;
	StringBuffer output = new StringBuffer();
	try{
	    br = new BufferedReader(new FileReader(table));
	    String curline = "";
	    
	    while((curline=br.readLine())!=null){
		if(!curline.startsWith("Line")){//process non-header line
		    String[] tokens = curline.split("\\t");
		    output.append(processLine(tokens));
		}else//process header line
		    output.append(curline + "\tTOTAL_GENERATIONS\n");//header = curline;
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	this.toOutFile(table+".converted" , output);
    }
    
    public void toOutFile(String outFile, StringBuffer content){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(outFile));
	    bw.write(content.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public String processLine(String[] tokens){
	StringBuffer errorTag = new StringBuffer("");
	double sum = 0.0d;
	StringBuffer tmp = new StringBuffer(tokens[0] + "\t");
	for(int i=1; i<tokens.length;i++){
	    Double tmpD = null;
	    if((tmpD = size2generations.get(tokens[i]))!=null){
		tmp.append(tmpD.toString() + "\t");
		sum += tmpD.doubleValue();
	    }else{
		tmp.append(tokens[i] + "\t" );
		if(tokens[i].trim().equals("") || tokens[i].trim().equals("0") || tokens[i].trim().equals("X") | tokens[i].trim().equals("x") | tokens[i].trim().equals("-"))
		    ;
		else
		    errorTag.append("col["+(i+1)+"]"+tokens[i]+":");
	    }
	}
	String errStr = errorTag.toString();
	if(errStr.length() == 0)
	    tmp.append(sum + "\n");
	else
	    tmp.append(errStr + "\n");
	
	return tmp.toString();
    }
    
    
}
