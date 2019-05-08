public class ExtComposition{
    int A;
    int a;
    int C;
    int c;
    int G;
    int g;
    int T;
    int t;

    public ExtComposition(int A, int a, int C, int c, int G, int g, int T, int t){
	this.A = A;
	this.a = a;
	this.C = C;
	this.c = c;
	this.G = G;
	this.g = g;
	this.T = T;
	this.t = t;
    }
    
    public ExtComposition(){
	A = 0;
	a = 0;
	C = 0;
	c = 0;
	G = 0;
	g = 0;
	T = 0;
	t = 0;
    }


    public char isMutation(char consensus){
	return this.isMutation(consensus, 0.8, 10);
    }

    public char isMutation(char consensus, double MAF, int minDepthEachWay){
	return this.isMutation(consensus, MAF, minDepthEachWay, 10);
    }

    public char isMutation(char consensus, double MAF, int minDepthEachWay, int minDepthTotal){
	int sum = A+a+C+c+G+g+T+t;
	if(sum < minDepthTotal)
	    return 0;
	int sumA = A+a;
	int sumC = C+c;
	int sumG = G+g;
	int sumT = T+t;
	if(consensus == 'A'){
	    if( (sumC*1.0)/(sum*1.0) > MAF && C >= minDepthEachWay && c >= minDepthEachWay)
		return 'C';
	    else if( (sumG*1.0)/(sum*1.0) > MAF && G >= minDepthEachWay && g >= minDepthEachWay)
		return 'G';
	    else if( (sumT*1.0)/(sum*1.0) > MAF && T >= minDepthEachWay && t >= minDepthEachWay)
		return 'T';
	    else
		return 0;
	}else if(consensus == 'C'){
	    if( (sumA*1.0)/(sum*1.0) > MAF && A >= minDepthEachWay && a >= minDepthEachWay)
		return 'A';
	    else if( (sumG*1.0)/(sum*1.0) > MAF && G >= minDepthEachWay && g >= minDepthEachWay)
		return 'G';
	    else if( (sumT*1.0)/(sum*1.0) > MAF && T >= minDepthEachWay && t >= minDepthEachWay)
		return 'T';
	    else
		return 0;
	}else if(consensus == 'G'){
	    if( (sumC*1.0)/(sum*1.0) > MAF && C >= minDepthEachWay && c >= minDepthEachWay)
		return 'C';
	    else if( (sumA*1.0)/(sum*1.0) > MAF && A >= minDepthEachWay && a >= minDepthEachWay)
		return 'A';
	    else if( (sumT*1.0)/(sum*1.0) > MAF && T >= minDepthEachWay && t >= minDepthEachWay)
		return 'T';
	    else
		return 0;
	}else if(consensus == 'T'){
	    if( (sumC*1.0)/(sum*1.0) > MAF && C >= minDepthEachWay && c >= minDepthEachWay)
		return 'C';
	    else if( (sumG*1.0)/(sum*1.0) > MAF && G >= minDepthEachWay && g >= minDepthEachWay)
		return 'G';
	    else if( (sumA*1.0)/(sum*1.0) > MAF && A >= minDepthEachWay && a >= minDepthEachWay)
		return 'A';
	    else
		return 0;
	}else
	    return 0;
    }

    public char getDominantAllele(){
	int sum = A+a+C+c+G+g+T+t;
	double rA = ((A+a)*1.0d) / (sum*1.0d);
	double rC = ((C+c)*1.0d) / (sum*1.0d);
	double rG = ((G+g)*1.0d) / (sum*1.0d);
	double rT = ((T+t)*1.0d) / (sum*1.0d);
	
	if(rA > 0.5)
	    return 'A';
	else if(rC > 0.5)
	    return 'C';
	else if(rG > 0.5)
	    return 'G';
	else if(rT > 0.5)
	    return 'T';
	else
	    return 'X';
    }
    
    public String toString(){
	return A + " " + a + " " + C + " " + c  + " " + G + " " + g + " " + T + " " + t;
    }
    public String toStringTab(){
	return A + "\t" + a + "\t" + C + "\t" + c  + "\t" + G + "\t" + g + "\t" + T + "\t" + t;
    }
}
