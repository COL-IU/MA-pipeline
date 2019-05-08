import java.util.*;

public class RunFinder{
    
    public static void main(String[] args){
	new RunFinder().findRuns(args[0]);
    }

    public void findRuns(String fastaFile){
	Hashtable<String, Seq> seqHash = new FastaReader().parseFasta(fastaFile);
	
	Enumeration<String> keys = seqHash.keys();
	while(keys.hasMoreElements()){
	    Seq curSeq = seqHash.get(keys.nextElement());
	    printRunCounts(curSeq.countRuns());
	    
	    printMatrix(curSeq.getNeighboringBaseCountMatrix());
	    
	    
	    
	}
    }

    public void printRunCounts(ArrayList<IntegerCounter> runCounts){
	for(int i=0; i<runCounts.size(); i++){
	    System.out.println((i+1) + "\t" + runCounts.get(i).getCount() + "\t" + runCounts.get(i).getATCount() + "\t" + runCounts.get(i).getGCCount());
	}
    }

    public void printMatrix(int[][][] mat){
	for(int i=0;i<mat.length;i++){
	    System.out.println("Matrix <" + index2Base(i) + ">\n");
	    System.out.println("\tA\tC\tG\tT");
	    for(int j=0;j<mat[i].length;j++){
		System.out.print(index2Base(j));
		for(int k=0;k<mat[i][j].length;k++){
		    System.out.print("\t" + mat[i][j][k]);
		}
		System.out.println();
	    }
	    System.out.println();
	}
    }

    private char index2Base(int i){
	if(i==0)
	    return 'A';
	else if(i==1)
	    return 'C';
	else if(i==2)
	    return 'G';
	else if(i==3)
	    return 'T';
	else
	    return 'X';
    }
    
}
