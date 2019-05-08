import java.util.*;

public class MutationSimulator{
    
    private int[][] mutationCounts;
    private double[] mutationProb;
    private int numDesiredMutations;
    private Seq seq;
    private int numLines;
    private int GCs;
    private int ATs;


    public MutationSimulator(int numMuts, int ATGC, int GCAT, int ATTA, int GCTA, int ATCG, int GCCG, Seq s, int nL){
	this.numDesiredMutations = numMuts;
	this.mutationProb = new double[6];
	this.mutationProb[0] = (ATGC*1.0d)/(numMuts*1.0d);
	this.mutationProb[1] = (GCAT*1.0d)/(numMuts*1.0d);
	this.mutationProb[2] = (ATTA*1.0d)/(numMuts*1.0d);
	this.mutationProb[3] = (GCTA*1.0d)/(numMuts*1.0d);
	this.mutationProb[4] = (ATCG*1.0d)/(numMuts*1.0d);
	this.mutationProb[5] = (GCCG*1.0d)/(numMuts*1.0d);
	this.mutationCounts = new int[6][2];
	this.mutationCounts[0][0] = ATGC;
	this.mutationCounts[1][0] = GCAT;
	this.mutationCounts[2][0] = ATTA;
	this.mutationCounts[3][0] = GCTA;
	this.mutationCounts[4][0] = ATCG;
	this.mutationCounts[5][0] = GCCG;

	this.GCs = s.countGC();
	this.ATs = s.getSeqLength() - GCs;
	
	//System.err.println(GCs + "\t" + ATs);

	this.mutationCounts[0][1] = ATs;
	this.mutationCounts[1][1] = GCs;
	this.mutationCounts[2][1] = ATs;
	this.mutationCounts[3][1] = GCs;
	this.mutationCounts[4][1] = ATs;
	this.mutationCounts[5][1] = GCs;
	this.seq = s;
	this.numLines = nL;
    }

    public void selectMutation(){
	int mutTypes[] = new int[numDesiredMutations];
	for(int i=0; i<this.numDesiredMutations; i++){
	    mutTypes[i] = (int)(Math.random()*numDesiredMutations);
	    int mutTypeIndex = this.getMutTypeIndex(mutTypes[i]);
	    if(mutTypeIndex < 0)
		System.err.println("FUCK");
	    int chosenNthBase = (int)(Math.random()*this.mutationCounts[mutTypeIndex][1]);
	    this.printMutation(mutTypeIndex, chosenNthBase);
	}
    }
        
    public void testAllCases(){
	for(int i=0;i<this.mutationCounts.length;i++){
	    if(i%2 == 0){//AT to smething
		for(int j=ATs-5; j<ATs;j++){
		    this.printMutation(i, j);
		}
	    }else{//GCtoSomething
		for(int j=GCs-5; j<GCs;j++){
		    this.printMutation(i,j);
		}
	    }
	}
    }

    public int getMutTypeIndex(int mutType){
	int curSum = 0;
	for(int i=0;i<this.mutationCounts.length;i++){
	    curSum += this.mutationCounts[i][0];
	    if(mutType < curSum)
		return i;//return (int)(Math.random()*this.mutationCounts[i][1]);
	}
	return -1;
    }

    public void printMutation(int mutTypeIndex, int chosenNthBase){
	int tmp = this.seq.getIndexForNthBase(mutTypeIndex, chosenNthBase, numLines);
	if(tmp < 0)
	    System.err.println("FAWK\t" + mutTypeIndex + "\t" + chosenNthBase);
    }
    

    public static void main(String[] args){
	if(args.length < 9)
	    System.err.println("USAGE: java MutationSimulator <#OfMutations> <#AT->GC> <#GC->AT> <#AT->TA> <#GC->TA> <#AT->CG> <#GC->CG> <genomeFasta> <numLines>");
	else{
	    MutationSimulator sim = new MutationSimulator(Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]), Integer.parseInt(args[3])
							  ,Integer.parseInt(args[4]), Integer.parseInt(args[5]), Integer.parseInt(args[6]), new FastaReader().parseFasta(args[7]).get("ecoli"), Integer.parseInt(args[8]));
	    
	    sim.selectMutation();
	    //sim.testAllCases();
	}
    }
}
