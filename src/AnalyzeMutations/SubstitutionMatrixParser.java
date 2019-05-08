import java.util.*;
import java.io.*;

public class SubstitutionMatrixParser{

    private SubstitutionMatrix sm;
    
    public SubstitutionMatrixParser(){
	this.sm = new SubstitutionMatrix(24);
    }
    
    public SubstitutionMatrix parse(String file){
	BufferedReader br = null;
	String curline = "";
	try{
	    br = new BufferedReader(new FileReader(file));
	    int count = 0;
	    while((curline=br.readLine())!=null){
		if(curline.startsWith("#"))//comment
		    ;
		else{
		    count++;
		    if(count < 2)
			;
		    else{
			String[] tokens = curline.split("\\s+");
			this.update(tokens, count-2);
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return this.sm;
	//System.out.println(this.sm.toString());
    }

    public void update(String[] tokens, int n){
	this.sm.addAAs(tokens[0].charAt(0), n);//adding amino acid to the hash
	//this.sm.setAAs(n, tokens[0].charAt(0));//setting the amino acid
	for(int i=1; i<tokens.length;i++){
	    this.sm.setMatrix(n,i-1,Integer.parseInt(tokens[i]));
	}
    }

    public static void main(String[] args){
	new SubstitutionMatrixParser().parse(args[0]);
    }
    
}

class SubstitutionMatrix{

    private char[] AAs;
    private Hashtable<Character, Integer> AAtoIndex;
    private int[][] matrix;
    private int counter;

    public SubstitutionMatrix(int size){
	this.AAs = new char[size];
	this.AAtoIndex = new Hashtable<Character, Integer>();
	this.matrix = new int[size][size];
	this.counter = 0;
    }

    public int getDistance(char pre, char post){
	return this.matrix[this.getIndexForAmino(post)][this.getIndexForAmino(pre)];
    }

    public int getIndexForAmino(char c){
	return this.AAtoIndex.get(new Character(c)).intValue();
    }
    
    public char[] getAAs(){
	return this.AAs;
    }
    
    public Hashtable<Character, Integer> getAAtoIndex(){
	return this.AAtoIndex;
    }
    
    public int[][] getMatrix(){
	return matrix;
    }

    public void addAAs(char c, int n){
	this.AAs[n] = c;
	this.AAtoIndex.put(new Character(c), new Integer(n));
    }

    public void setMatrix(int row, int col, int val){
	this.matrix[row][col] = val;
    }

    public String toString(){
	StringBuffer buffer = new StringBuffer();
	for(int i=0;i<this.AAs.length;i++){
	    buffer.append("\t" + this.AAs[i]);
	}
	buffer.append("\n");
	for(int i=0;i<this.matrix.length;i++){
	    buffer.append(this.AAs[i]);
	    for(int j=0;j<this.matrix[i].length;j++){
		buffer.append("\t" + this.matrix[i][j]);
	    }
	    buffer.append("\n");
	}
	return buffer.toString();
    }

    
}
