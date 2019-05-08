import java.io.*;
import java.util.*;

public class MakeList{
    
    public static void main(String[] args){
	// args[0] --> datadir    /data/groups/heewlee/MA/data/run12/M304a
	// args[1] --> strainname  M304
	// args[2] --> jobname     M304a
	// args[3] --> utildir     ${MAHOME}/util/
	if(args.length == 4)
	    process(args[0], args[1], args[2], args[3]);
	else{
	    System.err.println("java MakeList <data dir :(ex) /data/groups/heewlee/MA/data/run12/M304a > <strain name : (ex) M304 > <job name : (ex) M304a > <util dir : (ex) /data/groups/heewlee/MA/MA-Analysis_heewlee_run12/util >");
	}
    }
    
    //batchname is ex) M304a, M304b, etc
    //datadir
    public static void process(String datadir, String strainName, String jobname, String utildir){
	File projdir = new File(datadir);
	final String sName = strainName;
	FilenameFilter filter = new FilenameFilter() {
		public boolean accept(File dir, String name){
			System.out.println(name);
		    if(name.startsWith(sName + "-")){
			return true;
		    }
		    return false;
		}
	    };
	File[] files = projdir.listFiles(filter);
	Arrays.sort(files);
	MakeList.makeList(files, utildir, jobname);
    }

    public static void makeList(File[] files, String utildir, String jobname){
	StringBuffer content = new StringBuffer();
	StringBuffer contentWNum = new StringBuffer();
	for(int i=0;i<files.length;i++){
	    content.append(files[i].getName() + "\n");
	    contentWNum.append((i+1) + "\t" + files[i].getName() + "\n");
	}
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(utildir + File.separator + jobname + ".list"));
	    bw.write(content.toString());
	    bw.close();
	    
	    bw = new BufferedWriter(new FileWriter(utildir + File.separator + jobname + ".wNum.list"));
	    bw.write(contentWNum.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	bw = null;
    }
    
}
