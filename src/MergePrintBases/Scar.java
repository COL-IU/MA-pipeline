import java.util.Comparator;

public class Scar {
	private String geneName;
	private int start;
	private int end;
	private char ref;
	private char mut;
	
	public Scar(String geneName, int start, int end, char ref, char mut){
		this.geneName = geneName;
		if (start < end){
			this.start = start;
			this.end = end;
		}else{
			this.start = end;
			this.end = start;
		}
		this.ref = ref;
		this.mut = mut;
	}
	
	
	public String getGeneName() {
		return geneName;
	}
	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public char getRef() {
    return ref;
  }
  public void setRef(char ref) {
    this.ref = ref;
  }
	public char getMut() {
    return mut;
  }
  public void setMut(char mut) {
    this.mut = mut;
  }
	
	public String toString(){
		String ret = this.geneName + "-" + this.start + ":" + this.end + "|";
		return ret;
	}
}

class ScarComparator implements Comparator<Scar>{

	@Override
	public int compare(Scar arg0, Scar arg1) {
		if (arg0.getStart() < arg1.getStart()){
			return -1;
		}else if(arg0.getStart() == arg1.getStart()){
			return 0;
		}else{
			return 1;
		}
	}
	
}
