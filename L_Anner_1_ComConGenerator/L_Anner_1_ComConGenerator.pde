/* LAnner_1_ComConGenerator - Takes the consensus sequence of a given organism 
  (file name = *-SymCon.txt, located in exampleFiles in this repository) 
  aligned to the reference sequence (file name = *-PraCon.txt located in 
  exampleFiles in this repository) that is gapped in accordance with the 
  alignment. The script then outputs a sequence file encoded to display the 
  contents of both input files (file name = *-ComCon.txt) in accordance with 
  the key located in the main directory of the GitHub repo named "ComConKey.svg". 
*/

// Load query organism abbreviations.
String[] namelist=loadTable("../sharedResources/namelist.csv","header").getStringColumn("abb"); 
int ti=namelist.length;
for(int i=0;i<ti;i++){
  // Load SymCon and PraCon.
  String SymCon=loadStrings("../sharedResources/SymCon/"+namelist[i]+"-SymCon.txt")[0];
  String PraCon=loadStrings("../sharedResources/PraCon/"+namelist[i]+"-PraCon.txt")[0];
  
  /* If the documents have been exported correctly, and came from the same file,
     then they will always be the same length. If the documents are not the 
     same length, then this indicates that something is wrong and that the 
     process should not proceed.
  */ 
  if(PraCon.length()==SymCon.length()){
    int tj=SymCon.length();
    boolean same=true;
    // Output file (*-ComCon.txt) path/name is set here.
    PrintWriter ComCon=createWriter("../sharedResources/ComCon/"+namelist[i]+"-ComCon.txt");
    
    // Encode the ComCon sequence.
    for(int j=0;j<tj;j++){
      switch(PraCon.charAt(j)){
        case'A':
          switch(SymCon.charAt(j)){
            case'A': ComCon.print("A"); break;
            case'T': ComCon.print("."); same=false;break;
            case'G': ComCon.print("}"); same=false;break;
            case'C': ComCon.print(","); same=false;break;
            case'-': ComCon.print("a"); same=false;break;
            case'N': ComCon.print("N"); same=false;break;
            case'M': ComCon.print("M"); same=false;break;
            case'R': ComCon.print("R"); same=false;break;
            case'K': ComCon.print("!"); same=false;break;
            case'Y': ComCon.print("X"); same=false;break;
            case'S': ComCon.print("z"); same=false;break;
            case'V': ComCon.print("Z"); same=false;break;
            case'B': ComCon.print("B"); same=false;break;
            case'H': ComCon.print("f"); same=false;break;
            case'D': ComCon.print("F"); same=false;break;
            case'W': ComCon.print("W"); same=false;break;
          }
        break;
        case'T':
        switch(SymCon.charAt(j)){
            case'A': ComCon.print("e"); same=false;break;
            case'T': ComCon.print("T"); break;
            case'G': ComCon.print("\\");same=false;break;
            case'C': ComCon.print("/"); same=false;break;
            case'-': ComCon.print("t"); same=false;break;
            case'N': ComCon.print("~"); same=false;break;
            case'M': ComCon.print("q"); same=false;break;
            case'R': ComCon.print(":"); same=false;break;
            case'K': ComCon.print("K"); same=false;break;
            case'Y': ComCon.print("Y"); same=false;break;
            case'S': ComCon.print(";"); same=false;break;
            case'V': ComCon.print("V"); same=false;break;
            case'B': ComCon.print("L"); same=false;break;
            case'H': ComCon.print("l"); same=false;break;
            case'D': ComCon.print("J"); same=false;break;
            case'W': ComCon.print("^"); same=false;break;
          }
        break;
        case'G':
        switch(SymCon.charAt(j)){
            case'A': ComCon.print("|"); same=false;break;
            case'T': ComCon.print("'"); same=false;break;
            case'G': ComCon.print("G"); break;
            case'C': ComCon.print("&"); same=false;break;
            case'-': ComCon.print("g"); same=false;break;
            case'N': ComCon.print("U"); same=false;break;
            case'M': ComCon.print("u"); same=false;break;
            case'R': ComCon.print("I"); same=false;break;
            case'K': ComCon.print("i"); same=false;break;
            case'Y': ComCon.print("o"); same=false;break;
            case'S': ComCon.print("S"); same=false;break;
            case'V': ComCon.print("O"); same=false;break;
            case'B': ComCon.print("p"); same=false;break;
            case'H': ComCon.print("H"); same=false;break;
            case'D': ComCon.print("P"); same=false;break;
            case'W': ComCon.print("["); same=false;break;
          }
        break;
        case'C':
        switch(SymCon.charAt(j)){
            case'A': ComCon.print("\""); same=false;break;
            case'T': ComCon.print("?"); same=false;break;
            case'G': ComCon.print("*"); same=false;break;
            case'C': ComCon.print("C"); break;
            case'-': ComCon.print("c"); same=false;break;
            case'N': ComCon.print("9"); same=false;break;
            case'M': ComCon.print("8"); same=false;break;
            case'R': ComCon.print("7"); same=false;break;
            case'K': ComCon.print("6"); same=false;break;
            case'Y': ComCon.print("5"); same=false;break;
            case'S': ComCon.print("4"); same=false;break;
            case'V': ComCon.print("3"); same=false;break;
            case'B': ComCon.print("2"); same=false;break;
            case'H': ComCon.print("="); same=false;break;
            case'D': ComCon.print("D"); same=false;break;
            case'W': ComCon.print("]"); same=false;break;
          }
        break;
        case'-':
        switch(SymCon.charAt(j)){
            case'A': ComCon.print("@"); same=false;break;
            case'T': ComCon.print("+"); same=false;break;
            case'G': ComCon.print("{"); same=false;break;
            case'C': ComCon.print("("); same=false;break;
            case'-': ComCon.print("-"); break;
            case'N': ComCon.print("n"); same=false;break;
            case'M': ComCon.print("m"); same=false;break;
            case'R': ComCon.print("r"); same=false;break;
            case'K': ComCon.print("k"); same=false;break;
            case'Y': ComCon.print("y"); same=false;break;
            case'S': ComCon.print("s"); same=false;break;
            case'V': ComCon.print("v"); same=false;break;
            case'B': ComCon.print("b"); same=false;break;
            case'H': ComCon.print("h"); same=false;break;
            case'D': ComCon.print("d"); same=false;break;
            case'W': ComCon.print("w"); same=false;break;
          }
        break;
      }
    }
    
    // Save the ComCon file.
    ComCon.flush();
    ComCon.close();
    // It is unlikely that the two consensus sequences for anything but 
    // the same exact organism will be exactly the same. If two consensus sequences for 
    // a symbiont taxa are the exact same, an alert is given to the user, 
    // so they can double check the export was preformed correctly. 
    if(same){println("warning: "+namelist[i]+" pra- and sym- cons identical");}
  }else{
    println(namelist[i]+" pra- and sym- cons have unequal length"); 
  }
}
exit();
