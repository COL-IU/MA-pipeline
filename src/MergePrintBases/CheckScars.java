import java.io.*;
import java.util.*;

public class CheckScars {

	public static void main(String[] args) {
		if (args.length == 2) {
			try {
				BufferedReader sc = new BufferedReader(new FileReader(args[1]));
				ArrayList<Scar> scars = new ArrayList<Scar>();
				String scarLine = "";
				while((scarLine=sc.readLine())!=null){
					String[] words = scarLine.split("\t");
					//String[] startEnd = words[1].split("\\.\\.");
					//scars.add(new Scar(words[0],Integer.parseInt(startEnd[0]),Integer.parseInt(startEnd[1])));
					if (words.length == 3){
						scars.add(new Scar(words[0],Integer.parseInt(words[1]),Integer.parseInt(words[2]),'N','N'));
					}else{
						scars.add(new Scar(words[0],Integer.parseInt(words[1]),Integer.parseInt(words[2]),words[3].charAt(0),words[4].charAt(0)));
					}
				}
				sc.close();
				
				Collections.sort(scars,new ScarComparator());
				
				Iterator<Scar> it = scars.iterator();
				Scar curScar = it.next();
				int curStart = curScar.getStart();
				int curEnd = curScar.getEnd();
				char curRef = curScar.getRef();
				char curMut = curScar.getMut();
				int curScarIndex = 0;
				
				BufferedReader ap = new BufferedReader(new FileReader(args[0]));
				String tempLine = ap.readLine();
				int numLines = (tempLine.split(" ").length-3)/8;
				ap.close();
				
				double percentMissing[][] = new double[scars.size()][numLines];
				for (int i=0; i < scars.size(); i++)
					   for (int j=0; j < numLines; j++)
						   percentMissing[i][j] = 0;
				
				ap = new BufferedReader(new FileReader(args[0]));
				String line = "";
				int counter = 0;
				int prevline = 0; //ADDED [Heewook] ]06/23/14 : to process positions that are entirely skipped by all lines
			    outer:
				while((line=ap.readLine())!=null){
					counter++;
//					if (counter%100000 == 0){
//						System.out.println("Pos: "+counter);
//					}
					String[] tokens = line.split(" ");
					PBCount pbc = new PBCount(tokens, true);
					int position = pbc.getPos();
					
					/*----------------------------------------------------------------*/
					/* ADDED [Heewook] ]06/23/14 : to process positions that are entirely skipped by all lines */
					/*----------------------------------------------------------------*/
					if(position -1 > prevline){ //need to check if we skipped any
					    for(int innerPos=prevline+1;innerPos<position;innerPos++){
						if (innerPos >= curStart && innerPos <= curEnd){
						    if (curStart == curEnd){ // SNP
							for(int i=0; i < numLines; i++){
							    percentMissing[curScarIndex][i] = -1; //-1 is used for SNP [Heewook 06/23/14]
							}
						    }else{ // Gene Deletion
							for(int i=0; i < numLines; i++){
							    percentMissing[curScarIndex][i]++; //for deletions we do exactly same [Heewook 06/23/14]
							}
						    }
						}
						if (innerPos == curEnd){
						    if (it.hasNext()){
							curScar = it.next();
							curStart = curScar.getStart();
							curEnd = curScar.getEnd();
							curRef = curScar.getRef();
        			curMut = curScar.getMut();
							curScarIndex++;
						    }else{
							break outer;
						    }
						}
					    }
					}
					/*----------------------------------------------------------------*/
					/* END OF : ADDED [Heewook] ]06/23/14 : to process positions that are entirely skipped by all lines */
					/*----------------------------------------------------------------*/
					
					if (position >= curStart && position <= curEnd){
						if (curStart == curEnd){ // SNP
							for(int i=0; i < numLines; i++){
								if (pbc.getACGTList().get(i).isMutation(curRef)==0){
									percentMissing[curScarIndex][i] = 0;
								}else{
									if (pbc.getACGTList().get(i).isMutation(curRef)==curMut || curRef == 'N'){
										percentMissing[curScarIndex][i] = 1;
									}else{
										percentMissing[curScarIndex][i] = 0;
									}
								}
							}
						}else{ // Gene Deletion
							for(int i=0; i < numLines; i++){
								if (pbc.isEmpty(i)) percentMissing[curScarIndex][i]++;
							}
						}
					}
					if (position == curEnd){
						if (it.hasNext()){
							curScar = it.next();
							curStart = curScar.getStart();
							curEnd = curScar.getEnd();
							curRef = curScar.getRef();
        			curMut = curScar.getMut();
							curScarIndex++;
						}else{
							break outer;
						}
					}
					prevline = position;
				}
				ap.close();
				
				Scar scar;
				int length;
				for (int i=0; i < scars.size(); i++){
					scar = scars.get(i);
					length = scar.getEnd()-scar.getStart()+1;
					System.out.print(scar.getGeneName() + "(" + length + ")\t");
					for (int j=0; j < numLines; j++){
						//percentMissing[i][j] = percentMissing[i][j] * 100 / length;
						//if (percentMissing[i][j] > 20){
						//	System.out.print("1\t");
						//}else{
						//	System.out.print("0\t");
						//}
						System.out.print(percentMissing[i][j]+"\t");
					}
					System.out.print("\n");
				}
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		} else {
			System.err.println("USAGE: java CheckScars <alignPos> <scars>");
		}
	}

}
