import java.io.*;
import java.util.*;

/* 
 * This program simply looks up coverages for each lines in .list file and eliminate lines with less than MINCOV coverage
 * This doesn't physically eliminate these lines but they are removed from .list and .wNum.list files so that these lines
 * will not be processed further down.
 */
public class SelectLinesAndUpdateListFile{

    private int MINCOV;

    public SelectLinesAndUpdateListFile(){
	this.MINCOV = 30; //default of 30X
    }

    public SelectLinesAndUpdateListFile(int mincov){
	this.MINCOV = mincov;
    }

    public static void main(String[] args){
	if(args.length == 2)
	    new SelectLinesAndUpdateListFile().selectLines(args[0], args[1]);
	else if(args.length == 3)
	    new SelectLinesAndUpdateListFile().selectLines(args[0], args[1], args[2]);
	else
	    System.err.println("USAGE: java SelectLinesAndUpdateListFile <.list file> <outDir> [mainChromosomName (Optional, Default = ecoli) ]");
    }
    
    public void selectLines(String listFile, String outDir){
	this.selectLines(listFile, outDir, "ecoli");
    }
    
    public void selectLines(String listFile, String outDir, String mainChromosome){
	BufferedReader br = null;
	BufferedWriter listbw = null;
	BufferedWriter wNumListbw = null;
	String genome = mainChromosome;
	ArrayList<String> selectedLines = new ArrayList<String>();
	boolean needToSelect = false;//if false, we keep the original and exit.  If TRUE, we rewrite the list and wNum.list files
	try{
	    br = new BufferedReader(new FileReader(listFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		int curCov = this.getCoverage(outDir + File.separator 
					      + curline + File.separator 
					      + curline +"_"+genome+"_stat2.txt");
		if(curCov >= this.MINCOV){
		    selectedLines.add(curline);
		}else{
		    needToSelect = true;
		    System.err.println("@SelectLinesAndUpdateListFile:\tLine " + curline + " is getting droped. COVERAGE[" + curCov + "]");
		}
	    }
	    br.close();
	    
	    
	    if(needToSelect){
		if(renameOriginalFiles(listFile)){
		    System.err.println("@SelectLinesAndUpdateListFile:\tUpdating " + listFile + " and " + this.getStemFromListName(listFile) +".wNum.list");
		    listbw = new BufferedWriter(new FileWriter(listFile));
		    wNumListbw = new BufferedWriter(new FileWriter(getStemFromListName(listFile) + ".wNum.list"));
		    for(int i=0;i<selectedLines.size();i++){
			listbw.write(selectedLines.get(i)+"\n");
			wNumListbw.write((i+1) + "\t" + selectedLines.get(i)+"\n");
		    }
		    listbw.close();
		    wNumListbw.close();
		}
	    }else{
		System.err.println("@SelectLinesAndUpdateListFile:\tNO NEED TO UPDATE LIST FILES--> No Lines with less than " + this.MINCOV + "X avg coverage.");
	    }
	    
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private String getStemFromListName(String listFile){
	return listFile.substring(0,listFile.lastIndexOf(".list"));
    }

    private boolean renameOriginalFiles(String listFile){
	File lf = new File(listFile);
	String stem = getStemFromListName(listFile);
	File lwf = new File(stem + ".wNum.list");
	String newListName = stem + ".orig.list";
	String newWNumListName = stem + ".orig.wNum.list";
	return (lf.renameTo(new File(newListName)) && lwf.renameTo(new File(newWNumListName)));
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

}
