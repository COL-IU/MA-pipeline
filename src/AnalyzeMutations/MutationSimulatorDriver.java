import java.io.*;
import java.util.*;

public class MutationSimulatorDriver{

    private static final String MAHOME = "/data/groups/heewlee/MA/MA-Analysis_heewlee_run10/";
    private static final String OUTDIR = "/data/groups/heewlee/MA/BPS_Simulation";

    public static void main(String[] args){
	if(args.length == 1){
	    
	    BufferedReader br = null;
	    String curline = "";
	    try{
		br = new BufferedReader(new FileReader(args[0]));
		while((curline=br.readLine())!=null){
		    if(curline.charAt(0) != '#')
			driver(curline, OUTDIR);
		}
		br.close();
	    }catch(IOException ioe){
		ioe.printStackTrace();
	    }
	}else
	    System.err.println("java -jar MutationSimulator <parameterFile>");
    }


    public static void generateListFile(String strainName, int numLines){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(OUTDIR + File.separator + strainName+".simul.wNum.list"));
	    for(int i=1;i<=numLines;i++){
		bw.write(i + "\t" + strainName + "-" + i + "\n");
	    }
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    ////
    //#Strain Genotype        AT>GC   GC>AT   AT>TA   GC>TA   AT>CG   GC>CG   fastaFileName   seqName numLines        numRepetition
    public static void driver(String simulationFileLine, String outDir){
	String[] tokens = simulationFileLine.split("\\t");
	String strainName = tokens[0];
	String genomtype = tokens[1];
	int[] spectraCount = new int[6]; //AT>GC GC>AT AT>TA GC>TA AT>CG GC>CG
	int sum = 0;
	for(int i=2;i<8;i++){
	    spectraCount[i-2] = Integer.parseInt(tokens[i]);
	    sum += spectraCount[i-2];
	}
	String fastaFileName = tokens[8];
	String seqName = tokens[9];
	int numLines = Integer.parseInt(tokens[10]);
	int numRepetitions = Integer.parseInt(tokens[11]);

	generateListFile(strainName, numLines); //creates temporary list file

	File tmp = new File(outDir + File.separator + strainName);//creates the directory
	tmp.mkdirs();

	System.err.println("================================================================\n= Simulating parameters:\t" + simulationFileLine 
			   + "\n================================================================\n");
	
	BufferedWriter bw = null;
	BufferedReader br = null;
	try{
	    bw = new BufferedWriter(new FileWriter(outDir + File.separator + strainName + "." + numRepetitions +".simulated"));
	    for(int i=0;i<numRepetitions;i++){
		String curFile = outDir + File.separator + strainName + File.separator + strainName +"_"+i+".txt";
		new MutationSimulator(sum, spectraCount[0], spectraCount[1], spectraCount[2]
				      , spectraCount[3], spectraCount[4], spectraCount[5]
				      , new FastaReader().parseFasta(fastaFileName).get(seqName)
				      , numLines).selectMutation(curFile);
		//java -jar ${MAHOME}/bin/AnalyzeMutations.jar ${OUTDIR}/${JOB}.putations ${MAHOME}/util/standard_codon.txt ${MAHOME}/util/${REFBASE}.fna  ${MAHOME}/util/BLOSUM62.txt 0 Y ${NUMLIST} N ${OUTDIR}/${JOB}_snp ${PTTLIST}
		new AnalyzeMutations( MAHOME + "/util/standard_codon.txt" 
				      , fastaFileName
				      , MAHOME + "/util/BLOSUM62.txt"
				      , "0"
				      , "Y"
				      , OUTDIR + File.separator + strainName+".simul.wNum.list"
				      , "N"
				      , strainName + "_" + i
				      , AnalyzeMutations.getPttFiles(fastaFileName.substring(0,fastaFileName.lastIndexOf(".")) + ".ptts") 
				      ).process(curFile);
		br = new BufferedReader(new FileReader(outDir + File.separator + strainName + File.separator + strainName +"_"+i+".txt" +"."+seqName+".stat"));
		bw.write(br.readLine() + "\n");
		br.close();
		br = null;
	    }
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
    }

    
    //spectrumLine(tab-delimited) : StrainName AT>GC GC>AT AT>TA GC>TA AT>CG GC>CG
    public static void driver(String spectrumLine, String fastaFile, String seqName, int numLines, int numRepetitions){
	String[] tokens = spectrumLine.split("\\t");
	String strainName = tokens[0];
	int[] spectraCount = new int[6]; //AT>GC GC>AT AT>TA GC>TA AT>CG GC>CG
	int sum = 0;
	for(int i=1;i<7;i++){
	    spectraCount[i-1] = Integer.parseInt(tokens[i]);
	    sum += spectraCount[i-1];
	}
	
	for(int i=0;i<numRepetitions;i++){
	    new MutationSimulator(sum, spectraCount[0], spectraCount[1], spectraCount[2], spectraCount[3], spectraCount[4], spectraCount[5]
				  , new FastaReader().parseFasta(fastaFile).get(seqName), numLines).selectMutation(strainName+"_"+i+".txt");
	}

    }
}
