import java.io.*;

public class ToEasyFormat{

    public static void main(String[] args){
	new ToEasyFormat().process(args[0]);
    }

    public void process(String file){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(file));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		System.out.print(new PBCount(curline).toEasyStringByLine());
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

}
