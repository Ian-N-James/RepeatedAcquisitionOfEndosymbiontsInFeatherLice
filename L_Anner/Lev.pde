
class Lev{

  String pra;int numsym,gene;
  String[][] sym;
  
  Lev(String TEMPpra,String[] TEMPsym){
    pra=TEMPpra;numsym=TEMPsym.length;
    sym=new String[TEMPsym.length][2];
    for (int i=0;i<TEMPsym.length;i++){
      sym[i][0]=TEMPsym[i];
      sym[i][1]="";
    }
  }
  Lev(String TEMPpra,String[] TEMPsym,int templen){
    pra=TEMPpra;numsym=TEMPsym.length;
    sym=new String[TEMPsym.length][2];
    for (int i=0;i<TEMPsym.length;i++){
      sym[i][0]=TEMPsym[i];
      sym[i][1]="";
    }
    gene=templen;
  }
  String[][] Getbest(){
    Levenshtein a;
    for(int i=0;i<numsym;i++){
      a=new Levenshtein(sym[i][0],pra);
      sym[i][1]=str(a.distance());
    }
    IntList TGBWDT1=new IntList(); //There has Got to be a Better Way to Do This 1
    int[] TGBWDT2=new int[numsym]; //There has Got to be a Better Way to Do This 2
    //println(gene,numsym);
    for(int i=0;i<numsym;i++){
      TGBWDT2[i]=int(sym[i][1]);
    }
    int TGBWDT3=min(TGBWDT2); 
    TGBWDT2=new int[0];//clear TGBWDT2
    for(int i=0;i<numsym;i++){ 
      // the reason for doing things so that we take any that are at minimum, 
      // to address the posibility of two equally good/bad 
      if (int(sym[i][1])==TGBWDT3){
        TGBWDT1.append(i);
      }
    }
    String[][] out=new String[TGBWDT1.size()][2];
    for(int i=0;i<TGBWDT1.size();i++){
      out[i][0]=sym[TGBWDT1.get(i)][0];
      out[i][1]=sym[TGBWDT1.get(i)][1];
    }
    return out; //out[n][0]=the translation, out[n][1]=the score,
  }  
}

class Levenshtein {
  String a; String b;
  Levenshtein(String ta, String tb){
    a=ta;b=tb;
  }
  int distance() {
    a = a.toLowerCase();
    b = b.toLowerCase();
    // i == 0
    int[]costs=new int[b.length() + 1];
    for (int j = 0; j < costs.length; j++){
      costs[j] = j;
    }
    for (int i = 1; i <= a.length(); i++) {
      // j == 0; nw = lev(i - 1, j)
      costs[0] = i;
      int nw = i - 1;
      for (int j = 1; j <= b.length(); j++) {
        int cj = min(1+costs[j], 1+costs[j - 1],a.charAt(i - 1) == b.charAt(j - 1) ? nw : nw + 1);
        nw = costs[j];
        costs[j] = cj;
      }
    }
    return costs[b.length()];
  }
}
//class Lev{

//  String pra;int numsym,gene;
//  String[][] sym;
  
//  Lev(String TEMPpra,String[] TEMPsym){
//    pra=TEMPpra;numsym=TEMPsym.length;
//    sym=new String[numsym][2];
//    for (int i=0;i<numsym;i++){
//      sym[i][0]=TEMPsym[i];
//      sym[i][1]="";
//    }
//  }
//  Lev(String TEMPpra,String[] TEMPsym,int templen){
//    pra=TEMPpra;numsym=TEMPsym.length;
//    sym=new String[numsym][2];
//    for (int i=0;i<numsym;i++){
//      sym[i][0]=TEMPsym[i];
//      sym[i][1]="";
//    }
//    gene=templen;
//  }
//  String[][] Getbest(){
//    for(int i=0;i<numsym;i++){
      
//      sym[i][1]=str(Levenshtein_Distance(sym[i][0],pra));
//    }
//    IntList TGBWDT1=new IntList(); //There has Got to be a Better Way to Do This 1
//    int[] TGBWDT2=new int[numsym]; //There has Got to be a Better Way to Do This 2
//    //println(gene,numsym);
//    for(int i=0;i<numsym;i++){
//      TGBWDT2[i]=int(sym[i][1]);
//    }
//    int TGBWDT3=min(TGBWDT2); 
//    TGBWDT2=new int[0];//clear TGBWDT2
//    for(int i=0;i<numsym;i++){ 
//      // the reason for doing things so that we take any that are at minimum, 
//      // to address the posibility of two equally good/bad 
//      if (int(sym[i][1])==TGBWDT3){
//        TGBWDT1.append(i);
//      }
//    }
//    String[][] out=new String[TGBWDT1.size()][2];
//    int ti=TGBWDT1.size();
//    for(int i=0;i<ti;i++){
//      out[i][0]=sym[TGBWDT1.get(i)][0];
//      out[i][1]=sym[TGBWDT1.get(i)][1];
//    }
//    return out; //out[n][0]=the translation, out[n][1]=the score,
//  }  
//}

//int Levenshtein_Distance(String a,String b) {
  
//  a = a.toLowerCase();
//  b = b.toLowerCase();
//  // i == 0
//  int[]costs=new int[b.length() + 1];
//  for (int j = 0; j < costs.length; j++){
//    costs[j] = j;
//  }
//  for (int i = 1; i <= a.length(); i++) {
//    // j == 0; nw = lev(i - 1, j)
//    costs[0] = i;
//    int nw = i - 1;
//    for (int j = 1; j <= b.length(); j++) {
//      int cj = min(1+costs[j], 1+costs[j - 1],a.charAt(i - 1) == b.charAt(j - 1) ? nw : nw + 1);
//      nw = costs[j];
//      costs[j] = cj;
//    }
//  }
//  return costs[b.length()];
//}