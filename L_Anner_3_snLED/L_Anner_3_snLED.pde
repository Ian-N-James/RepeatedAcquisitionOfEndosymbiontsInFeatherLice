/* Takes the combined sequence file (*-ComCon.txt), and the reformatted annotation 
   table (*-AnnTabFin.csv). Calculates all open reading frames (ORFs) with either 
   canonical or alternative bacterial start codons within ± 25 bp of the ortholog 
   in the reference genome. The script then determines the ORF with the lowest 
   Levenshtein edit distance (LED) from the reference ortholog for each gene. 
   Exports a modified annotation table (ZFin-*-AnnTabFin.csv) that additionally 
   includes the LED, a size normalized LED (snLED), and both the predicted 
   nucleotide and protein sequences of the ORF with the lowest LED, for each gene.
*/
int onGenome=0;
int Genome;String ComCon;Table AnnTab;
int Tag;boolean rRNA,tRNA,miscRNA,psHS,PrPgO,PrPgT,Plas,Mis,RPS,foundMatch;
float percentForIntact;String gene,split1,split2,praTrans;
String[][]symTrans;String[][] Assesment;String[] rna=new String[1];
String[] mbtSeq;
int befAf=25;
String[] namelist;
Table Namelist;
float[] Cutoffs;
String[] posStarts= {"TTG","CTG","GTG","ATG","ATT","ATC","ATA","ATN","ATK","ATM","ATR","ATY","ATS","ATW","ATB","ATV","ATH","NTG","KTG","MTG","RTG","YTG","STG","WTG","BTG","VTG","HTG"};
void setup(){
  Namelist=loadTable("../sharedResources/namelist.csv","header");
  int ti=Namelist.getRowCount();
  namelist=new String[ti];
  Genome=namelist.length;
  Cutoffs=new float[ti];
  namelist=Namelist.getStringColumn("abb");
  Cutoffs=Namelist.getFloatColumn("CutOff");
  for(int i=onGenome;i<Genome;i++){
    ComCon=loadStrings("../sharedResources/ComCon/"+namelist[onGenome]+"-ComCon.txt")[0];
    AnnTab=loadTable("../sharedResources/AnnTabFin/"+namelist[onGenome]+"-AnnTab.csv","header");
    Lev DTFA;//RawSeq Split;translator geneTrans;
    AnnTab.insertColumn(5,"ann",Table.FLOAT);
    AnnTab.insertColumn(6,"score",Table.INT);
    AnnTab.insertColumn(7,"weightScore",Table.FLOAT);
    AnnTab.insertColumn(8,"flags",Table.STRING);
    int tj=AnnTab.getRowCount();
    for(int j=0;j<tj;j++){
      AnnTab.setFloat(j,"weightScore",20);
    }
    AnnTab.addColumn("trans",Table.STRING);
    AnnTab.addColumn("seq",Table.STRING);
    for(int j=0;j<tj;j++){
      boolean reverse=false;
      rRNA=false;tRNA=false;miscRNA=false;psHS=false;PrPgO=false;PrPgT=false;Plas=false;Mis=false;RPS=false;
      if(AnnTab.getString(j,"dir").charAt(0)=='r'){
        reverse=true;
      }
      //Deal with Tags
      Tag=-1;rRNA=false;tRNA=false;miscRNA=false;psHS=false;PrPgO=false;PrPgT=false;Plas=false;Mis=false;RPS=false;
      //Acount for tRNA rRNA plasmid and HS pseudo genes
      String cTag=AnnTab.getString(j,"tag");char type=cTag.charAt(cTag.length()-5);
      int tNum=int(cTag.substring(cTag.length()-4));
      if(type=='_'){
        if(tNum>=2857 && tNum<=2946){
          AnnTab.setString(j,"flags","PpIs2");PrPgT=true;
          Tag=tNum+60000;
        }else{
          if(tNum>=2513 && tNum<=2554){
            AnnTab.setString(j,"flags","PpIs1");PrPgO=true;
            Tag=tNum+50000;
          }else{
            if(tNum==773||tNum==877||tNum==888||tNum==1078||tNum==1113){
              AnnTab.setString(j,"flags","miscRNA");miscRNA=true;
              Tag=tNum+30000;
            }else{
              Tag=tNum;
            }
          }
        }
      }else{
        if(type=='r'){
          Tag=10000+tNum;
          AnnTab.setString(j,"flags","rRNA");
          rRNA=true;
        }else{
          if(type=='t'){
            Tag=20000+tNum;
            AnnTab.setString(j,"flags","tRNA");tRNA=true;
          }else{
            if(type=='s'){
              if(cTag.substring(5,6).equals("p")){
                if(tNum>=2857 && tNum<=2946){
                  Tag=65000+tNum;
                  AnnTab.setString(j,"flags","PpIs2,pseuHS");PrPgT=true;psHS=true;
                }else{
                  if(tNum>=2513 && tNum<=2554){
                    Tag=55000+tNum;
                    AnnTab.setString(j,"flags","PpIs1,pseuHS");PrPgO=true;psHS=true;
                  }else{
                    Tag=5000+tNum;
                    AnnTab.setString(j,"flags","pseuHS");psHS=true;
                  }
                }
              }else{
                if(cTag.substring(5,6).equals("P")){
                  Tag=45000+tNum;
                  AnnTab.setString(j,"flags","Plas,pseuHS");
                  Plas=true;
                  psHS=true;
                }
              }
            }else{ 
              if(type=='P'){
                Tag=40000+tNum;
                AnnTab.setString(j,"flags","Plas");Plas=true;
              }
            }
          }
        }
      }
      if(Tag==0){AnnTab.setString(j,"flags","RPstar");RPS=true;}
      if(Tag==791){AnnTab.setString(j,"flags","PrTrFr");}
      if(Tag==-1){
        println("bad tag ");
        println(cTag.substring(cTag.length()-4));
        println(type);
        exit();
      }
      // 5000-9500 = Pseudo in HS,
      // 10000s = rRNA 20000s = tRNA, 30000s = misc. RNA, 40000s plasmid, 45000s = plasmid and pseudo, 
      // 50000s = prophage island 1 & nonpseudo in HS, 55000s = prophage island 1 & pseudo in HS,
      // 60000s = prophage island 2 & nonpseudo in HS, 65000s = prophage island 2 & pseudo in HS, 
      //Start Doing Translation stuff;
      int geneMin=AnnTab.getInt(j,"min")-1;
      int geneMax=AnnTab.getInt(j,"max");
      String praSeq=getPra(ComCon,geneMin,geneMax);
      if(Tag==791){praSeq=praSeq.substring(0,75)+praSeq.substring(75+1);}
      String symSeq=getSym(ComCon,geneMin,geneMax);
      if(symSeq.length()==0){Mis=true;AnnTab.setFloat(j,"weightScore",1);AnnTab.setFloat(j,"ann",0);}//streamline the detection of missing genes
      if(geneMin>=befAf){
        symSeq=getSym(ComCon,geneMin-befAf,geneMin)+symSeq+getSym(ComCon,geneMax,geneMax+befAf);
      }else{
        symSeq=getSym(ComCon,0,geneMin)+symSeq+getSym(ComCon,geneMax,geneMax+befAf);
      }
      gene="";
      if(!rRNA&&!tRNA&&!miscRNA&&!psHS&&!Mis&&Tag!=0){
        praTrans=transPra(praSeq,reverse);
        symTrans=(Tag!=791)?transSym(symSeq,reverse):PTF_handler(symSeq,reverse,befAf+75);
        if(symTrans[0].length!=0){
          DTFA=new Lev(praTrans,symTrans[0],Tag);
          Assesment=DTFA.Getbest();
          AnnTab.setInt(j,"score",int(Assesment[0][1]));
          AnnTab.setFloat(j,"weightScore",float(Assesment[0][1])/(float(geneMax-(geneMin))/3));
          AnnTab.setFloat(j,"ann",0.25);//placeholder value
          if(Assesment.length==1){
            AnnTab.setString(j,"trans",Assesment[0][0]);
            foundMatch=false;
            int tk=symTrans[0].length;
            for(int k=0;(k<tk&&!foundMatch);k++){
              if(symTrans[0][k].equals(Assesment[0][0])){
                foundMatch=true;
                AnnTab.setString(j,"seq",symTrans[1][k]);
              }
            }
          }else{
            int tk=Assesment[0].length;
            for(int k=0;k<tk;k++){
              gene+=Assesment[k][0];
              if(k!=Assesment[0].length-1){
                gene+="\t";
              }
            }
            mbtSeq=new String[gene.split("\t").length];
            AnnTab.setString(j,"trans",gene);
            int tl=gene.split("\t").length;
            for(int l=0;l<tl;l++){
              foundMatch=false;
              for(int m=0;(m<symTrans[0].length&&!foundMatch);m++){
                if(symTrans[0][m].equals(gene.split("\t")[l])){
                  foundMatch=true;
                  mbtSeq[l]=symTrans[1][m];
                }
              }
            }
            gene=join(mbtSeq,",");
            AnnTab.setString(j,"seq",gene);
            gene="";
          }
        }else{
          AnnTab.setFloat(j,"ann",0.0);
          AnnTab.setFloat(j,"weightScore",1);
        }
      }else{
        if((rRNA||tRNA||miscRNA||psHS||Tag==0)&&!Mis){
          AnnTab.setFloat(j,"ann",-0.1);
          AnnTab.setFloat(j,"weightScore",20);
        }
      }
    }
    saveTable(AnnTab,"../sharedResources/zFin/zFin"+"-"+namelist[onGenome]+"-AnnTabFin.csv");
    println("done with "+namelist[onGenome]+" ("+(onGenome+1)+"/"+Genome+")");
    onGenome++;
    if(onGenome==Genome){println("Success!");exit();}
  }
}

void draw(){}
String[][] transSym(String in,boolean R){
  if(R){in=revComp(in);}
  int[] startPositions=FindGivenStarts(in,posStarts).toArray();
  int ti=startPositions.length;
  String[][] out=new String[2][ti];
  for(int i=0;i<ti;i++){
    out[0][i]=doTrans(in.substring(startPositions[i]));
    out[1][i]=in.substring(startPositions[i],startPositions[i]+out[0][i].length()*3);
  }
  return out;
}
String[][]PTF_handler(String in,boolean R,int pos){
  String[][] thru1=transSym(in,R);
  String in2=in.substring(0,pos)+in.substring(pos+1);
  String[][] thru2=transSym(in2,R);
  String[][] out={concat(thru1[0],thru2[0]),concat(thru1[1],thru2[1])};
  return out;
}
String transPra(String in,boolean R){
  if(R){in=revComp(in);}
  return doTrans(in);
}
IntList FindGivenStarts(String Seq,String[] startCod){
  IntList output=new IntList();
  int th=startCod.length,ti=Seq.length()-2;
  for(int h=0;h<th;h++){
    char F=startCod[h].charAt(0);char S =startCod[h].charAt(1);char T=startCod[h].charAt(2);
    for(int i=0;i<ti;i++){
      if(Seq.charAt(i)==F&&Seq.charAt(i+1)==S&&Seq.charAt(i+2)==T){
        output.append(i);
      }
    }
  }
  return output;
}
String revComp(String in){
  String out="";
  int inlen=in.length();
  for(int i=0;i<inlen;i++){
    switch(in.charAt(inlen-i-1)){
      case'A': out+="T"; break;
      case'T': out+="A"; break;
      case'G': out+="C"; break;
      case'C': out+="G"; break;
      case'-': out+="-"; break;
      case'?': out+="?"; break;
      case'N': out+="N"; break;
      case'M': out+="K"; break;
      case'R': out+="Y"; break;
      case'K': out+="M"; break;
      case'Y': out+="R"; break;
      case'S': out+="W"; break;
      case'V': out+="B"; break;
      case'B': out+="V"; break;
      case'H': out+="D"; break;
      case'D': out+="H"; break;
      case'W': out+="S"; break;
    }
  }
  return out;
}


String getPra(String in, int start, int stop){
  String out="";
  for(int i=start;i<stop;i++){
    switch(in.charAt(i)){
      case 'A': case '.': case '}': case ',': case 'a': case 'N': case 'M': case 'R': case '!': case 'X': case 'z': case 'Z': case 'B': case 'f': case 'F': case 'W': out+="A"; break;
      case 'e': case 'T': case '\\':case '/': case 't': case '~': case 'q': case ':': case 'K': case 'Y': case ';': case 'V': case 'L': case 'l': case 'J': case '^': out+="T"; break;
      case '|':case '\'': case 'G': case '&': case 'g': case 'U': case 'u': case 'I': case 'i': case 'o': case 'S': case 'O': case 'p': case 'H': case 'P': case '[': out+="G"; break;
      case '"': case '?': case '*': case 'C': case 'c': case '9': case '8': case '7': case '6': case '5': case '4': case '3': case '2': case '=': case 'D': case ']': out+="C"; break;
      default: break;
    }
  }
  return out;
}
String getSym(String in, int start, int stop){
  String out="";
  for(int i=start;i<stop;i++){
    switch(in.charAt(i)){
      case'A':case'e': case'|':case'"':case'@': out+="A";break;
      case'.':case'T':case'\'':case'?':case'+': out+="T";break;
      case'}':case'\\':case'G':case'*':case'{': out+="G";break;
      case',':case'/': case'&':case'C':case'(': out+="C";break;
      case'N':case'~':case'U':case'9':case'n': out+="N";break;
      case'M':case'q':case'u':case'8':case'm': out+="M";break;
      case'R':case':':case'I':case'7':case'r': out+="R";break;
      case'!':case'K':case'i':case'6':case'k': out+="K";break;
      case'X':case'Y':case'o':case'5':case'y': out+="Y";break;
      case'z':case';':case'S':case'4':case's': out+="S";break;
      case'Z':case'V':case'O':case'3':case'v': out+="V";break;
      case'B':case'L':case'p':case'2':case'b': out+="B";break;
      case'f':case'l':case'H':case'=':case'h': out+="H";break;
      case'F':case'J':case'P':case'D':case'd': out+="D";break;
      case'W':case'^':case'[':case']':case'w': out+="W";break;
      default: break;
    }
  }
  return out;
}

String doTrans(String nSeq){//works on a nucleotide sequence starting at a pre-found start codon (regular or alternate)
  String pSeq="M";
  boolean foundStop=false;
  for(int i=3;i<nSeq.length()-2&&!foundStop;i+=3){
    switch(nSeq.charAt(i)){
      case 'A': //A__
      switch(nSeq.charAt(i+1)){
        case 'A': //AA_
        switch(nSeq.charAt(i+2)){
          case'T': case'C': case'Y': pSeq+="N"; break;
          case'A': case'G': case'R': pSeq+="K"; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'T'://AT_
        switch(nSeq.charAt(i+2)){
          case'A': case'T': case'C': case'Y': case'W': case'M':case'H': pSeq+="I"; break;
          case'G': pSeq+="M"; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'G'://AG_
        switch(nSeq.charAt(i+2)){
          case'T': case'C': case'Y': pSeq+="S"; break;
          case'A': case'G': case'R': pSeq+="R"; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'C': pSeq+="T"; break;
        default:pSeq+="X"; break;
      }
      break;
      case 'T'://T__
      switch(nSeq.charAt(i+1)){
        case 'A'://TA_
        switch(nSeq.charAt(i+2)){
          case'T': case'C': case'Y': pSeq+="Y"; break;
          case'A': case'G': case'R': pSeq+="*";foundStop=true; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'T'://TT_
        switch(nSeq.charAt(i+2)){
          case'T': case'C': case'Y': pSeq+="F"; break;
          case'A': case'G': case'R': pSeq+="L"; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'G'://TG_
        switch(nSeq.charAt(i+2)){
          case'T':case'C':case'Y': pSeq+="C"; break;
          case'A': pSeq+="*"; foundStop=true; break;
          case'G': pSeq+="W"; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'C': pSeq+="S";break;
        case 'R': if(nSeq.charAt(i+2)=='A'){pSeq+="*";foundStop=true;}else{pSeq+="X";}break;
        default:pSeq+="X"; break;
      }
      break;
      case 'G':
      switch(nSeq.charAt(i+1)){
        case 'A'://GA_
        switch(nSeq.charAt(i+2)){
          case'T': case'C': case'Y': pSeq+="D"; break;
          case'A': case'G': case'R': pSeq+="E"; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'T': pSeq+="V"; break;//GT_
        case 'G': pSeq+="G"; break;//GG_
        case 'C': pSeq+="A"; break;//GC_
        default:pSeq+="X"; break;
      }
      break;
      case 'C'://C__
      switch(nSeq.charAt(i+1)){
        case 'A'://CA_
        switch(nSeq.charAt(i+2)){
          case'T': case'C': case'Y': pSeq+="H"; break;
          case'A': case'G': case'R': pSeq+="Q"; break;
          default:pSeq+="X"; break;
        }
        break;
        case 'T': pSeq+="L"; break;
        case 'G': pSeq+="R"; break;
        case 'C': pSeq+="P"; break;
        default:pSeq+="X"; break;
      }
      break;
      case 'M': if(nSeq.charAt(i+1)=='G'&&(nSeq.charAt(i+2)=='G'||nSeq.charAt(i+2)=='A'||nSeq.charAt(i+2)=='Y')){pSeq+="R";}else{pSeq+="X";}break;
      case 'Y': if(nSeq.charAt(i+1)=='T'&&(nSeq.charAt(i+2)=='G'||nSeq.charAt(i+2)=='A'||nSeq.charAt(i+2)=='Y')){pSeq+="L";}else{pSeq+="X";}break;
      default:pSeq+="X"; break;
    }
     
  }
  return pSeq;
}
