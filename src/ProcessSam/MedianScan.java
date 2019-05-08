import java.util.*;
import java.io.*;

public class MedianScan{
    double[] medMAD; //array of 2 elements. index 0:median, index 1: MAD
    long totalLength;
 

    ArrayList<Integer> inserts;
    
    public MedianScan(){
	this.medMAD = new double[2];
	//this.median = 0.0d;
	this.totalLength = 0L;
	this.inserts = new ArrayList<Integer>();
    }

    public boolean isEmptySet(){
	if(this.medMAD[0] == 0.0d && this.medMAD[1] == 0.0d)
	    return true;
	return false;
    }

    public double[] getMedMAD(){
	return this.medMAD;
    }

    public void toString(String samfile){
	System.out.println(samfile + "\t" + medMAD[0] + "\t" + medMAD[1] );
    }

    public String medMADToString(){
	return "Median:\t" + medMAD[0] + "\tMAD:\t" + medMAD[1];
    }

    void updateStat(SAM s1, SAM s2){
	SAMPair sp = new SAMPair(s1, s2);
	int tmpISize = sp.getInsertSize();
	if(tmpISize > 0){
	    this.inserts.add(new Integer(sp.getInsertSize()));
	    this.totalLength += sp.sumOfReadLengths();
	}
    }
    

    public void finalizeStat(String samfileName){
	int numData = this.inserts.size();
	if(numData > 0){
	Collections.sort(this.inserts);
	if(numData % 2 == 1)
	    this.medMAD[0] = this.inserts.get(numData/2)*1.0d;
	else
	    this.medMAD[0] = ( (this.inserts.get(numData/2) + this.inserts.get(numData/2 - 1))*1.0d ) / (2*1.0d);
	//return new double[] {this.median, this.calculateMAD()};
	this.medMAD[1] = this.calculateMAD();
	    
	    //BufferedWriter bw = null;
	    //try{
	    //    bw = new BufferedWriter(new FileWriter(samfileName+"_dist.txt"));
	    //    StringBuffer buffer = new StringBuffer();
	    //   for(int i=0; i<this.inserts.size();i++){
	    //	buffer.append(inserts.get(i).intValue()+"\n");
	    //   }
	    //   bw.write(buffer.toString());
	    //   bw.close();
	    //}catch(IOException ioe){
	    //   ioe.printStackTrace();
	    //   }
	}else{
	    this.medMAD[0] = 0.0d;
	    this.medMAD[1] = 0.0d;
	}
	this.inserts = null;
    }

    private double calculateMAD(){
	double [] absDevs = new double[this.inserts.size()];
	for(int i=0; i<this.inserts.size();i++){
	    absDevs[i] = Math.abs(this.inserts.get(i).doubleValue() - this.medMAD[0]);
	}
	Arrays.sort(absDevs);
	if(absDevs.length % 2 == 1)
	    return absDevs[absDevs.length/2];
	else
	    return ( (absDevs[absDevs.length/2] + absDevs[absDevs.length/2 - 1]) / 2.0d );
    }

}
