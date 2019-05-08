import java.io.*;


/*
 * SelectLines class: This is a simple class that selectively choose only certain lines from alignPos file. 
 * Since each column(except for header columns) in alignPos file is a MA line. it will take a flat 
 * file(tab-delimited) where the first column is a positive integer. 
 * 
 * e.g. given a file "lines.txt" that has following lines:
 *  
 * 4
 * 10
 * 15
 *
 *
 * and a alignPos file with 20 lines will output a truncated alignPos file with just lines 4, 10, and 15.
 */
public class SelectLines{
    
    private int[] selectionArray;

    public static void main(String[] args){
	if(args.length == 2)
	    new SelectLines(args[0]).process(args[1]);
	else
	    System.err.println("java SelectLines <selectionFile> <alignPosFile>");
    }
    
    public SelectLines(String selectionFile){
	this.loadSelectionArray(selectionFile);
    }
    

    public void process(String alignPosFile){
	BufferedReader br =null;
	String curline = "";
	try{
	    br = new BufferedReader(new FileReader(alignPosFile));
	    while((curline=br.readLine())!=null){
		System.out.println(new PBCount(curline).toStringBuffer(selectionArray).toString());
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public int[] getSelectionArray(){
	return this.selectionArray;
    }

    private void loadSelectionArray(String file){
	BufferedReader br = null;
	String curline = "";
	int count = 0;
	try{
	    br = new BufferedReader(new FileReader(file));
	    while((curline=br.readLine())!=null){
		if(Character.isDigit(curline.charAt(0)))
		    count++;
	    }
	    br.close();
	    this.selectionArray = new int[count];
	    br = new BufferedReader(new FileReader(file));
	    int i=0;
	    while((curline=br.readLine())!=null){
		this.selectionArray[i] = Integer.parseInt(curline.split("\\t")[0]) - 1;
		i++;
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
}
