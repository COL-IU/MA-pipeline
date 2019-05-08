
//Counts all 64 triplets in a segment.
//usage: java Triplets ecoli.fna 1604045 3923881 -   ====> revserse strand
//usage: java Triplets ecoli.fna 3923882 1604044 +   ====> fwd strand
public class Triplets{

    //ACGT
    //0123
    //
    //0 --> 000 : AAA
    //1 --> 001 : AAC
    //2 --> 002 : AAG
    //3 --> 003 : AAT
    //4 --> 010 : ACA
    //5 --> 011 : ACC
    //6 --> 012 : ACG
    public static String index2Triplet(int n){
	if(n<0 && n>63)
	    return null;
	else{
	    int tmp = n;
	    
	    int numPower2 = tmp/16;
	    tmp = tmp - numPower2*16;
	    int numPower1 = tmp/4;
	    tmp = tmp - numPower1*4;
	    int numPower0 = tmp;
	    
	    return "" 
		+ digit2char(numPower2) 
		+ digit2char(numPower1) 
		+ digit2char(numPower0);
	}
    }

    public static int triplet2Index(String trp){
	if(trp.length() == 3){
	    return (char2digit(trp.charAt(0)) * 16 )
		+ (char2digit(trp.charAt(1)) * 4 )
		+ (char2digit(trp.charAt(2)) * 1);
	}else
	    return -1;
	    
    }

    public static int char2digit(char b){
	if(b == 'A')
	    return 0;
	else if(b == 'C')
	    return 1;
	else if(b == 'G')
	    return 2;
	else if(b == 'T')
	    return 3;
	else if(b == 'a')
	    return 0;
	else if(b == 'c')
	    return 1;
	else if(b == 'g')
	    return 2;
	else if(b == 't')
	    return 3;
	else 
	    return -1;
    }

    public static char digit2char(int d){
	if(d == 0)
	    return 'A';
	else if(d==1)
	    return 'C';
	else if(d==2)
	    return 'G';
	else if(d==3)
	    return 'T';
	else
	    return 'N';
    }	   

    public static void printCounts(int[] tripletCounts){
	for(int i=0; i<tripletCounts.length;i++)
	    System.out.println(index2Triplet(i) + "\t" + tripletCounts[i]);
    }
    
    //args[0] : fasta
    //args[1] : start
    //args[2] : end 
    //args[3] : +/- for strand
    //args[4](optional) : sequence name (default:ecoli)
    public static void main(String[] args){
	/*for(int i=0; i<64;i++){
	    System.out.println("i: " + i + " --> " + index2Triplet(i) );
	}
	for(int i=0;i<4;i++){
	    for(int j=0; j<4;j++){
		for(int k=0; k<4; k++){
		    String curTrip = ""+digit2char(i)+digit2char(j)+digit2char(k);
		    System.out.println(curTrip + " --> " + triplet2Index(curTrip) );
		}
	    }
	    }*/
	Seq seq = null;
	if(args.length == 4)
	    seq = new FastaReader().parseFasta(args[0]).get("ecoli");
	else
	    seq = new FastaReader().parseFasta(args[0]).get(args[4]);
	seq.countTriplets(Integer.parseInt(args[1]), Integer.parseInt(args[2]), (args[3].equals("+") ? true : false));
    }
}
