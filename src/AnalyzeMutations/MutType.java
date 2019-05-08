public class MutType{

    private int synVal; /* 1 if Syn, 0 if non-syn, -1, if NG due to bases other than atgc */
    private char prevAmino;
    private char afterAmino;
    private String prevTriplet;
    private String afterTriplet;
    private boolean isConservative;
    private Gene gene;
    //int dist;
    public MutType(int synVal, char pA, char aA, String pT, String aT){
	this.synVal = synVal;
	this.prevAmino = pA;
	this.afterAmino = aA;
	this.prevTriplet = pT;
	this.afterTriplet = aT;
	this.isConservative = true;
	this.gene = null;
    }

    public MutType(int synVal, char pA, char aA, String pT, String aT, Gene g){
	this(synVal,pA,aA,pT,aT);
	this.gene = g;
    }

    public Gene getGene(){
	return this.gene;
    }

    public String mutString(){
	StringBuffer bf = new StringBuffer(this.prevTriplet + "\t" + this.afterTriplet + "\t" + this.prevAmino + "\t" + this.afterAmino + "\t");
	if(this.gene == null)
	    bf.append("-\t-");
	else
	    bf.append(this.gene.getGName() + "\t" + this.gene.directionStr());
	return bf.toString();
    }
    
    public void setGene(Gene g){
	this.gene = g;
    }

    public int getSynVal(){
	return this.synVal;
    }
    
    public char getPrevAmino(){
	return this.prevAmino;
    }
    
    public char getAfterAmino(){
	return this.afterAmino;
    }
    
    public boolean isConservative(){
	return this.isConservative;
    }
    
    /* with threshold = 0, greater TRUE --> 0 or greaters : conservative, otherwise nonconservative*/
    public boolean isConservative(SubstitutionMatrix sm, int threshold, boolean greater){
	//System.out.println("here");
	int dist = sm.getDistance(this.prevAmino, this.afterAmino);
	//this.dist = dist;
	if(greater){
	    if(dist < threshold)
		this.isConservative = false;
	}else{
	    if(dist > threshold)
		this.isConservative = false;
	}
	return this.isConservative;
    }
}
