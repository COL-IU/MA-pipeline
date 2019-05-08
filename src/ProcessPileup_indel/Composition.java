class Composition{
    int A;
    int a;
    int C;
    int c;
    int G;
    int g;
    int T;
    int t;

    public Composition(int a, int t, int g, int c){
        this.a = a;
        this.t = t;
        this.g = g;
        this.c = c;
    }

    public Composition(){
        this.a = 0;
        this.t = 0;
        this.g = 0;
        this.c = 0;
    }

    public Composition(Count ct ){
        this.a = ct.a;
        this.t = ct.t;
        this.g = ct.g;
        this.c = ct.c;
    }

    public boolean isFilled(){
	return ((a+t+g+c)>0 ? true : false);
    }

    public double alleleFrequency(char base){
	base = Character.toLowerCase(base);
	int sum = a+t+g+c;
	int num = 0;
	if(base == 'a')
	    return (this.a*1.0d) / (sum*1.0d);
	else if(base == 't')
	    return (this.t*1.0d) / (sum*1.0d);
	else if(base == 'g')
	    return (this.g*1.0d) / (sum*1.0d);
	else if(base =='c')
            return (this.c*1.0d) / (sum*1.0d);
	else{
	    System.exit(1);
	    return -1.0d;
	}
    }

    public char getMutBase(char ref){
	//System.err.println("B4 : " + ref);
	if(ref == 'a' || ref == 'A'){
	    if(this.t > 0)
		return 't';
	    else if(this.g > 0)
		return 'g';
	    else if(this.c > 0)
		return 'c';
	}else if(ref == 't' || ref == 'T'){
	    if(a > 0)
		return 'a';
	    else if(g > 0)
		return 'g';
	    else if(c > 0)
		return 'c';
	}else if(ref == 'g' || ref == 'G'){
	    if(a > 0)
		return 'a';
	    else if(t > 0)
		return 't';
	    else if(c > 0)
		return 'c';
	}else if(ref == 'c' || ref == 'C'){
	    if(a > 0)
		return 'a';
	    else if(t > 0)
		return 't';
	    else if(g > 0)
		return 'g';
	}//else{
	
	//System.err.println(this.toString());
	//System.exit(1);
	return 'x';
	//}
    }
    
    public int indexToCount(int baseIndex){
	if(baseIndex == 0)
	    return this.a;
	else if(baseIndex == 1)
	    return this.t;
	else if(baseIndex == 2)
	    return this.g;
	else if(baseIndex == 3)
	    return this.c;
	else
	    System.exit(0);
	return -1;
    }
    
    /* 
     * returns 1 if a+t+g+c becomes 0 --> need to reduce numPopCount in count object
     * otherwise returns 0
     */
    public int discardAllele(int i){
	int tmpSum = a + t + g + c;
	if(i == 0)
	    this.a = 0;
	else if(i == 1)
	    this.t = 0;
	else if(i == 2)
	    this.g = 0;
	else if(i == 3)
	    this.c = 0;
	
	if( tmpSum > 0 && (a+t+g+c) == 0)
	    return 1;
	else
	    return 0;
    }
    
    /* default value for minCount is 1.*/
    public int isPureType(){
	return isPureType(1);
    }

    /*
     * this one determines whether it's puretype or not.
     * @returns : -1 if it's not puretype, otherwise it returns 0 --> a , 1 --> t, 2 --> g, 3 --> c
     *
     */
    public int isPureType(int minCount){
	/*/int sum = a+t+g+c;
	if(sum >= minCount){
	    int aVal = (a+sum-1)/sum;
	    int tVal = (t+sum-1)/sum;
	    int gVal = (g+sum-1)/sum;
	    int cVal = (c+sum-1)/sum;
	    
	    int snpVal = aVal + tVal + gVal + cVal;
	    if(snpVal == 1){
		if(a == sum)
		    return 0;
		else if(t == sum)
		    return 1;
		else if(g == sum)
		    return 2;
		else if(c == sum)
		    return 3;
	    }
	}
	return -1;*/
	int sum = a+t+g+c;
	if(sum >= minCount){
	    if(a == sum)                                                                                                                                                            
		return 0;                                                                                                                                                           
	    else if(t == sum)                                                                                                                                                       
		return 1;                                                                                                                                                           
	    else if(g == sum)                                                                                                                                                       
		return 2;                                                                                                                                                           
	    else if(c == sum)                                                                                                                                                       
		return 3;
	}
	return -1;
    }

    public String toString(){
        return a + "\t" + t + "\t" + g + "\t" + c;
    }

    public double[] getAlleleFreqArr(){
	int sum = a+t+g+c;
	if(sum > 0){
	    int snpVal = (a+sum-1)/sum + (t+sum-1)/sum + (g+sum-1)/sum + (c+sum-1)/sum; // numAlleles
	    double[] alleleFreqArr = new double[snpVal];
	    
	    int curIndex = 0;
	    if(a != 0){
		alleleFreqArr[curIndex] = (a*1.0)/(sum*1.0);
		curIndex++;
	    }
	    
	    if(t != 0){
		alleleFreqArr[curIndex] = (t*1.0)/(sum*1.0);
		curIndex++;
	    }

	    if(g != 0){
		alleleFreqArr[curIndex] = (g*1.0)/(sum*1.0);
		curIndex++;
	    }
	    
	    if(c != 0){
		alleleFreqArr[curIndex] = (c*1.0)/(sum*1.0);
		curIndex++;
	    }
	
	    return alleleFreqArr;
	}else
	    return new double[0];
    }

    /* This method is only used for biallelic sites(overall population)*/
    /* This returns -1.0d if this microbiome does not have any sequencing coverages*/
    public double getFirstAlleleFreq(int[] ai){
	int tmpTot = this.indexToCount(ai[0]) + this.indexToCount(ai[1]);
	if(tmpTot == 0)
	    return -1.0d;
	return (1.0d * this.indexToCount(ai[0])) / (1.0d * tmpTot);
    }

    public double[] getAlleleFreq(String observed){
	int sum = a + t + g + c;
	if(sum > 0){
	    double[] freq = new double[2];
	    int curIndex = 0;
	    char base1 = Character.toLowerCase(observed.charAt(0));
	    char base2 = Character.toLowerCase(observed.charAt(1));
	    if(base1 == 'a')
		freq[0] = (a*1.0)/(sum*1.0);
	    else if(base1 == 't')
		freq[0] = (t*1.0)/(sum*1.0);
	    else if(base1 == 'g')
		freq[0] = (g*1.0)/(sum*1.0);
	    else if(base1 == 'c')
		freq[0] = (c*1.0)/(sum*1.0);
	    
	    if(base2 == 'a')
                freq[1] = (a*1.0)/(sum*1.0);
            else if(base2 =='t')
		freq[1] = (t*1.0)/(sum*1.0);
            else if(base2 == 'g')
                freq[1] = (g*1.0)/(sum*1.0);
            else if(base2 == 'c')
		freq[1] = (c*1.0)/(sum*1.0);
	    
	    return freq;
	}else
	    return null;
	    
    }
    
    
}


