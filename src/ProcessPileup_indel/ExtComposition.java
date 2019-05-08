public class ExtComposition{
    int A;
    int a;
    int C;
    int c;
    int G;
    int g;
    int T;
    int t;
    int N;
    int n;

    ExtComposition(int A, int a, int C, int c, int G, int g, int T, int t, int N, int n){
	this.A = A;
	this.a = a;
	this.C = C;
	this.c = c;
	this.G = G;
	this.g = g;
	this.T = T;
	this.t = t;
	this.N = N;
	this.n = n;
    }
    
    ExtComposition(){
	A = 0;
	a = 0;
	C = 0;
	c = 0;
	G = 0;
	g = 0;
	T = 0;
	t = 0;
	N = 0;
	n = 0;
    }
    
    public void update(char c){
	if(c == 'A')
	    this.A++;
	else if(c == 'a')
	    this.a++;
	else if(c == 'C')
	    this.C++;
	else if(c == 'c')
	    this.c++;
	else if(c == 'G')
	    this.G++;
	else if(c == 'g')
	    this.g++;
	else if(c == 'T')
	    this.T++;
	else if(c == 't')
	    this.t++;
	else if(c == 'N')
	    this.N++;
	else if(c == 'n')
	    this.n++;
    }
    /*
     * Simply return ATGC spectrum in the order of ATGC
     * any base not greater than MINCOUNT will be discarded(set to 0)
     * simplified --> not distinguishing the directionalities
     */
    public int[] getSimplifiedACGTspectrum(int MINCOUNT){
	int[] out = { this.checkMin(MINCOUNT,this.A + this.a)
			  , this.checkMin(MINCOUNT,this.C + this.c)
			  , this.checkMin(MINCOUNT,this.G + this.g)
			  , this.checkMin(MINCOUNT,this.T + this.t)};
	return out;
    }
    
    private int checkMin(int MINCOUNT, int depth){
	return ( (depth >= MINCOUNT) ? depth : 0 );
    }
    
    /* this returns a int array size of 8 elements AaCcGgTt */
    public int[] getACGTspectrum(){
	int[] out = {A,a,C,c,G,g,T,t};
	return out;
    }

}
