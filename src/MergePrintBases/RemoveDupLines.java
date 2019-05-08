import java.io.*;


public class RemoveDupLines{

    public static void main(String[] args){
	if(args.length < 2){
	    System.err.println("USAGE: java RemoveDupLines <easyStringEalignPosFile> <line num 1> < line num 2> ... <line num n>");
	    System.exit(0);
	}
	int[] removal = new int[args.length-1];
	for(int i=1; i<args.length; i++){
	    removal[i-1] = Integer.parseInt(args[i]);
	}
	new RemoveDupLines().run(args[0], removal);
    }
    
    public void run(String eAlignPos, int[] removal){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(eAlignPos));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		PBCount pb = new PBCount(curline, true);
		pb.removeLines(removal);
		System.out.println(pb.toEasyStringSpecial());
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
}
