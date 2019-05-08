import java.util.*;

public class UsageCounter{

    private Hashtable<String, Counter> codon2Count;
    private int total;

    public UsageCounter(){
	this.codon2Count = new Hashtable<String, Counter>();
    }

    public UsageCounter(String codon){
	this();
	this.addCodon(codon);
    }

    public void update(String codon){
	Counter tmp = this.codon2Count.get(codon);
	//if(tmp !=null){
	//   this.codon2Count.put(codon, new Counter(1));
	//}else{
	tmp.increment();
	total++;
	    //}
    }
    
    public void calculateUsage(char amino){
	Enumeration<String> codons = codon2Count.keys();
	while(codons.hasMoreElements()){
	    String curCodon = codons.nextElement();
	    System.out.println(amino + "\t" + curCodon + "\t" + ( (this.codon2Count.get(curCodon).getCount() * 1.0d)/(this.total * 1.0d) ) + "\t" + this.total + "\t" + this.codon2Count.get(curCodon).getCount());
	    
	}
    }

    public void addCodon(String codon){
	this.codon2Count.put(codon, new Counter());
    }
    

}

class Counter{
    
    private int count;
    
    Counter(){
	this(0);
    }

    Counter(int n){
	this.count = n;
    }

    public int increment(){
	//	System.out.print("BEFORE : " + this.count + "\t");
	return this.increment(1);
    }
    
    public int increment(int n){
	this.count += n;
	//System.out.println("AFTER\t" + this.count);
	return this.count;
    }

    public int getCount(){
	return this.count;
    }

}
