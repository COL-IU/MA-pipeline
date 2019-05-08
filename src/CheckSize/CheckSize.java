import java.io.*;

public class CheckSize{
    
    private String dataLoc;
    private String lineRootName; //ex: M173, M144, M2, M5, ...

    public CheckSize(){
    
    }

    public static void main(String[] args){

    }
    
    public void process(){
	File loc = new File(this.dataLocl);
	File[] files = loc.listFiles();
	for(int i=0; i<files.length; i++){
	    String curName = files[i].getName();
	    if(this.isFromLine(curName)){
		System.out.printlnt(files[i].length());
	    }
	}
    
    }
    

    public void process(String wNumList){
	BufferedReader br = null;
	try{
	
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public boolean isFromLine(String f){
	if(f.startsWith(lineRootName+"-") && f.endsWith("_1.fq"))
	    return true;
	    
    }
    
    public String stripPairTag(String f){
	return f.substring(0,f.indexOf("_1.fq"));
    }

}
