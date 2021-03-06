import java.io.*;
import java.util.*;

public class GenerateScript{


    public static void main(String[] args){
	new GenerateScript(args);
    }

    private String base;     /* -B(baseName)     : e.g. mutL11 should be the directory name of single line */
    private String refDir;   /* -R(refDir)       : e.g. /N/dc/projects/microbiome/MA/ref                   */
    private String refSeq;   /* -r(refSeq)       : e.g. ecoli.fasta--> Should be in <refDir>. There should be bwa index file for the ref in the same director */
    private String dataLoc;  /* -D(DataLoc)      : e.g. /N/dc/projects/microbiome/MA/ecoli/reads --> under <DataLoc>/<baseName>/<baseName_1/2.fq> should be there */
    private String outLoc;   /* -o(outLoc)       : e.g. /N/dc/projects/microbiome/MA/ecoli/mapping/mapping_5_error   */
    private int rLen;        /* -l(readLen)      : e.g. 90 (IF rLen is longer than 101, default algorithm is bwa-mem */ // THIS GETS OVERIDDEN by reading the first fastq sequence.
    private String jobName;  /* -N(jobName)      : e.g. ecoli_mutL_MA */
    private int numCpuBwa;   /* -n(numCpu)       : e.g. 4             */
    private int numMismatch; /* -e(errorPercent) : e.g. 5 --> 5% with readLenth 90 will translate to 4 (floor of 0.05*90) */
                             /* -m(numMismatch)  : e.g. 4 (ONLY use either -e or -m option) */
    private boolean skipMap; /* -s(skipMapping)  : e.g. -s T (to Skip) or -s F (not to Skip - default) */
    private String refSeqBase;
    private HashSet<String> flagSet;
    
    private boolean isMismatchSet;
    private boolean enableSeeding = false;   /* -S(enable seeding)    : e.g. -S T (to enable) or -s F (to disable - default)  */
    private boolean worstCaseSeeding = true; /* -W(worstcase seeding) : e.g. -W T (to enable - default or -W F (to disable) */
    private boolean uniquePairOnly = true;   /* -U(unique pair only)  : e.g. -U T (default) or -U F */
    private boolean keepSortedBam = false;

    private boolean useMEM = true;

    private String readSuffix = ".fq.gz";

    private String binDir = "/data/groups/heewlee/MA/MA-Analysis_heewlee_run10/bin";
    private String headerFileList = null;

    private int qThresh = 20;//minimum mapping quality score. /* -MQ(minMappingQuality) : e.g. 20 (default)*/

    public GenerateScript(String[] args){
	this.isMismatchSet = false;
	this.skipMap = false;
	if(args.length > 1 && args.length%2 == 0){
	    this.flagSet = new HashSet<String>();
	    this.loadClassVars(args);
	    this.outToFile(this.generate());
	}else{
	    System.err.println("USAGE: java GenerateScript -B <baseName> -R <refDir> -r <refSeq> -D <DataLoc> -o <outLoc> -l <readLen> -N <jobName> -n <numCpu> -e <errorPercent> -m <numMismatch> -s <skipMap T/F> -S <enable seeding T/F(default)> -W <wosrtcase seeding T(default)/F>");
	    System.err.println("-B(baseName)     : e.g. mutL11 should be the directory name of single line\n-R(refDir)       : e.g. /N/dc/projects/microbiome/MA/ref\n-r(refSeq)       : e.g. ecoli.fasta--> Should be in <refDir>. There should be bwa index file for the ref in the same directory\n-D(DataLoc)      : e.g. /N/dc/projects/microbiome/MA/ecoli/reads --> under <DataLoc>/<baseName>/<baseName_1/2.fq> should be there\n-o(outLoc)       : e.g. /N/dc/projects/microbiome/MA/ecoli/mapping/mapping_5_error\n-l(readLen)      : e.g. 90\n-N(jobName)      : e.g. ecoli_mutL_MA\n-n(numCpu)       : e.g. 4\n-e(errorPercent) : e.g. 5 --> 5% with readLenth 90 will translate to 4 (floor of 0.05\n-m(numMismatch)  : e.g. 4 (ONLY use either -e or -m option)\n-S(enable seeding) : e.g. -S T (to enable) or -s F (to disable - default)\n-W(worstcase seeding) : e.g. -W T (to enable - default or -W F (to disable)\n-U(unique pair only)  : e.g. -U T (default) or -U F\n-K(keep sorted bam) : e.g. -K T (to keep) -K F (default: to delete)\n-FS(fastq suffix): e.g. -FS <suffix i.e. .fq.gz - default>\n-M(MAHOME) : e.g. -M /data/groups/heewlee/MA/MA-Analysis_heewlee_run10/bin\n -H(headerFileList) : e.g. -H /data/groups/heewlee/MA/MA-Analysis_heewlee_run10/ecoli.headers");
	}
    }

    public void outToFile(StringBuffer bf){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.jobName + ".sh"));
	    bw.write(bf.toString());
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private int getReadLengthFromFirstReadsFile(String file){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(file));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.split("\\t");
		if(tokens[0].equals(this.base)){
		    br.close();
		    return tokens[1].length();
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	System.err.println("[@GenerateScript]\tREADLENGTH INFORMATION not found in " + file);
	System.err.println("[@GenerateScript]\tNOT WRITING SCRIPTS");
	return -1;
    }
    public StringBuffer generate(){
	this.rLen = this.getReadLengthFromFirstReadsFile(this.outLoc + File.separator + "firstReads.list");
	if(this.rLen == -1)
	    System.exit(0);
	//System.err.println("this.numCpuBwa is :\t"+this.numCpuBwa);
	StringBuffer mapCriteria = new StringBuffer(""+this.numCpuBwa); 
	//System.err.println("MAPCRITERIA IS :\t" + mapCriteria);
	if(!this.useMEM){
	    mapCriteria.append(" -n " + this.numMismatch);
	    if(worstCaseSeeding)
		mapCriteria.append(" -k " + this.numMismatch);
	    if(!enableSeeding)
		mapCriteria.append(" -l " + (this.rLen+1));
	}
	StringBuffer bf = new StringBuffer();
	
	bf.append("#!/bin/bash\n");
	bf.append("#" + this.base + "\n\n" );
	bf.append("#Alignment via bwa\n");
	if(this.skipMap)
	    bf.append("#");
	
	String outputBase = this.outLoc + File.separator + this.base + File.separator + this.base;
	String dataBase = this.dataLoc + File.separator + this.base + File.separator + this.base;
	
	if(this.useMEM){
	    bf.append("bwa mem -M -t " + mapCriteria.toString() + " " 
		      + this.refDir + File.separator + this.refSeq + " " 
		      + dataBase + "_1" + this.readSuffix + " "
		      + dataBase + "_2" + this.readSuffix 
		      + " > " + this.outLoc + File.separator + this.base + File.separator + this.base + "_PE\n\n");
	}else{
	    bf.append("bwa aln -t " + mapCriteria.toString() + " " + this.refDir + File.separator + this.refSeq + " " + dataBase + "_1" + this.readSuffix + " > " 
		      + outputBase + "_1.fq.sai;bwa aln -t " + mapCriteria.toString() + " " + this.refDir + File.separator 
		      + this.refSeq + " " + dataBase + "_2" + this.readSuffix +" > " + outputBase 
		      + "_2.fq.sai;bwa sampe -P -s " + this.refDir + File.separator + this.refSeq + " " + outputBase
		      + "_1.fq.sai " + outputBase + "_2.fq.sai " + dataBase + "_1" + this.readSuffix + " "
		      + dataBase + "_2" + this.readSuffix + " > " + this.outLoc + File.separator + this.base + File.separator
		      + this.base + "_PE;rm " + outputBase + "_1.fq.sai " + outputBase + "_2.fq.sai\n\n" );
	}
	bf.append("#Extracting Unique Mappings\n");
	if(!this.uniquePairOnly)
	    bf.append("#");

	bf.append("cd " + this.outLoc + File.separator + this.base + ";java -Xmx7g -Xms3g -jar " + this.binDir + File.separator + "ProcessAllSAMs.jar " 
		  + this.base + "_PE " + this.qThresh + " 0.05 " + this.headerFileList + " " 
		  + this.refDir + File.separator + this.refSeqBase + "_gList.txt\n\n");
	bf.append("#Getting pileups\n");
	String PEfile;
	if(this.uniquePairOnly)
	    PEfile = this.base + "_" + this.refSeqBase + ".qc.PE";
	else
	    PEfile = this.base + "_PE";
	bf.append("cd " + this.outLoc + File.separator + this.base + ";samtools view -Sb -t " + this.refDir + File.separator + this.refSeqBase + "_gList.txt " + PEfile
		  + " -o tmp.bam;samtools sort tmp.bam tmp.sorted;samtools mpileup -f " + this.refDir + File.separator + this.refSeq + " tmp.sorted.bam > " 
		  + PEfile + ".pileup;" + (keepSortedBam ? ("mv tmp.sorted.bam " + this.base + "_PE.bam;") : "") + "rm tmp.bam tmp.sorted.bam *PE\n");
	
	//System.out.println(bf.toString());
	return bf;
    }



    
    private void loadClassVars(String[] args){
	for(int i=0;i<args.length;i=i+2){
	    processArg(args[i], args[i+1]);
	}
	
    }
    
    private int percent2NumMismatch(int percent){
	return (int) ( (this.rLen*1.0d) * ( (percent*1.0d)/100.0d ) );
    }

    private void processArg(String flag, String val){
	if(flag.equals("-B") && !isDup(flag))
	    this.base = val;
	else if(flag.equals("-r") && !isDup(flag)){
	    this.refSeq = val;
	    this.refSeqBase = this.refSeq.substring(0,this.refSeq.lastIndexOf("."));
	}else if(flag.equals("-R") && !isDup(flag))
	    this.refDir = val;
	else if(flag.equals("-D") && !isDup(flag))
	    this.dataLoc = val;
	else if(flag.equals("-l") && !isDup(flag)){
	    this.rLen = Integer.parseInt(val);
	    if(!this.isMismatchSet){
		this.numMismatch = percent2NumMismatch(this.numMismatch);
		this.isMismatchSet = true;
	    }
	}
	else if(flag.equals("-N") && !isDup(flag)){
	    this.jobName = val;
	}else if(flag.equals("-n") && !isDup(flag)){
	    //System.err.println("####-->"+val);
	    this.numCpuBwa = Integer.parseInt(val);
	}else if(flag.equals("-e") && !isDup(flag)){
	    if(flagSet.contains("-l")){
		this.numMismatch = percent2NumMismatch(Integer.parseInt(val));
		this.isMismatchSet = true;
	    }else{
		this.numMismatch = Integer.parseInt(val);
	    }
	}else if(flag.equals("-m") && !isDup(flag)){
	    this.numMismatch = Integer.parseInt(val);
	    this.isMismatchSet = true;
	}else if(flag.equals("-o") && !isDup(flag)){
	    this.outLoc = val;
	}else if(flag.equals("-s") && !isDup(flag)){
	    if(val.equals("T"))
		this.skipMap = true;
	}else if(flag.equals("-S") && !isDup(flag)){
	    if(val.equals("T"))
		this.enableSeeding = true;
	    else
		this.enableSeeding = false;
	}else if(flag.equals("-W") && !isDup(flag)){
	    if(val.equals("F"))
		this.worstCaseSeeding = false;
	    else
		this.enableSeeding = true;
	}else if(flag.equals("-U") && !isDup(flag)){
	    if(val.equals("F"))
		this.uniquePairOnly = false;
	    else
		this.enableSeeding = true;
	}else if(flag.equals("-K") && !isDup(flag)){
	    if(val.equals("T"))
		this.keepSortedBam = true;
	}else if(flag.equals("-FS") && !isDup(flag)){
	    this.readSuffix = val;
	}else if(flag.equals("-M") && !isDup(flag)){
	    this.binDir = val + File.separator + "bin";
	}else if(flag.equals("-H") && !isDup(flag)){
	    this.headerFileList = val;
	}else if(flag.equals("-MQ") && !isDup(flag)){
	    this.qThresh = Integer.parseInt(val);
	}else{
	    System.err.println("Unknown option flag : " + flag);
	    System.exit(1);
	}
	    
    }
    
    private boolean isDup(String flag){
	if(flag.equals("-m")) /* this is to allow ONLY 1 between -m and -e options */
	    flag = "-e";
	if(flagSet.contains(flag)){
	    System.err.println("You have duplicate flags : " + flag);
	    System.exit(1);
	    return true;
	}else{
	    flagSet.add(flag);
	    return false;
	}
    }
 	


    public static void processBwa(String inFile)throws IOException{
	processBwa(inFile, 1);
    }

    public static void processBwa(String inFile, int multi)throws IOException{
	BufferedReader br = new BufferedReader(new FileReader(inFile));
	String curLine = "";
	String seqLoc = "~/ebi_gut/";
	ArrayList<StringBuffer> bfList = new ArrayList<StringBuffer>();
	for(int i=0; i<multi; i++){
	    bfList.add(new StringBuffer());
	}

	int count = 0;
	while((curLine=br.readLine()) != null){
	    String[] tokens = curLine.split("\\t");
	    String read1 = tokens[3].substring(tokens[3].lastIndexOf("/")+1, tokens[3].lastIndexOf(".gz")) + ".fixed";
            String read2 = tokens[4].substring(tokens[4].lastIndexOf("/")+1, tokens[3].lastIndexOf(".gz")) + ".fixed";
	    String header = read1.substring(0,read1.indexOf("_"));
	    

	    bfList.get(count%multi).append("bwa aln /nobackup/heewlee/gutGenomes/Bacteroides_vulgatus_ATCC_8482/NC_009614.fna "
			  + seqLoc + read1 + " > " + header + File.separator + read1 + ".sai;"
			  + "bwa aln /nobackup/heewlee/gutGenomes/Bacteroides_vulgatus_ATCC_8482/NC_009614.fna "
			  + seqLoc + read2 + " > " + header + File.separator + read2 + ".sai;"
			  
			  + "bwa sampe -P -s /nobackup/heewlee/gutGenomes/Bacteroides_vulgatus_ATCC_8482/NC_009614.fna " 
			  + header + File.separator + read1 + ".sai " + header + File.separator + read2 + ".sai " + seqLoc + read1 + " " + seqLoc + read2 + " > " 
			  + header + File.separator + read1.substring(0,read1.indexOf(".1.fq")) + "_on_NC_000914_PE\n" 
					   
					   /*+ "bwa samse -n 100 /nobackup/heewlee/gutGenomes/Bacteroides_vulgatus_ATCC_8482/NC_009614.fna "
			  + header + File.separator + read1 + ".sai " + seqLoc + read1 + " > " 
			  + header + File.separator + read1.substring(0,read1.indexOf(".fq")) + "_on_NC_000914_SE\n"
			  
			  + "bwa samse -n 100 /nobackup/heewlee/gutGenomes/Bacteroides_vulgatus_ATCC_8482/NC_009614.fna "
			  + header + File.separator + read2 + ".sai " + seqLoc + read2 + " > "
			  + header + File.separator + read2.substring(0,read2.indexOf(".fq")) + "_on_NC_000914_SE\n"*/
			  );

	    count++;
	}
	br.close();
	outToFile(bfList, "bwa");
    }

    public static void processNovoalign(String inFile, int multi)throws IOException{
	BufferedReader br = new BufferedReader(new FileReader(inFile));
	String curLine = "";
	String seqLoc = "~/ebi_gut/";
	ArrayList<StringBuffer> bfList = new ArrayList<StringBuffer>();
	for(int i=0; i<multi; i++){
	    bfList.add(new StringBuffer());
	}

	int count = 0;
	while((curLine=br.readLine()) != null){
	    String[] tokens = curLine.split("\\t");
	    String read1 = tokens[3].substring(tokens[3].lastIndexOf("/")+1, tokens[3].lastIndexOf(".gz"));
            String read2 = tokens[4].substring(tokens[4].lastIndexOf("/")+1, tokens[3].lastIndexOf(".gz"));
	    

	    bfList.get(count%multi).append("novoalign -d ~/gutGenomes/Bacteroides_vulgatus_ATCC_8482/NC_009614.fna.nix -r Random -F STDFQ -f "
					   + seqLoc + read1 + " " + seqLoc + read2 + " > " 
					   + read1.substring(0,read1.indexOf(".1.fq")) + "_on_NC_009614_novo.out\n"
					   
					   + "novo2sam.pl " + read1.substring(0,read1.indexOf(".1.fq")) 
					   + "_on_NC_009614_novo.out " + read1.substring(0,read1.indexOf(".1.fq")) 
					   + "_on_NC_009614_novo.out.sam" //+
					   );
	    count++;
	}
	br.close();
	outToFile(bfList, "novo");
    }

    private static void outToFile(ArrayList<StringBuffer> bfList, String prefix)throws IOException{
	BufferedWriter bw = null;
	for(int i=0; i<bfList.size(); i++){
	    bw = new BufferedWriter(new FileWriter(prefix + "_cmdList_"+i));
	    bw.write(bfList.get(i).toString());
	    bw.close();
	}
    }
      
    //    public static void main(String[] args)throws IOException{
	//processSoap(args[0]);
	//processBwa(args[0], Integer.parseInt(args[1]));
	//processNovoalign(args[0], Integer.parseInt(args[1]));
    //}

}

//mrfast --search /nobackup/heewlee/gutGenomes/Bacteroides_vulgatus_ATCC_8482/NC_009614.fna --seqcomp --pe --seq1 ~/mass_storage/ebi_gut/O2.UC-11_090112.1.fq.gz --seq2 ~/mass_storage/ebi_gut/O2.UC-11_090112.2.fq.gz  -o O2.UC-11_on_NC_009614.txt --discordant-vh --min 99 --max 185
