/* Takes the Shannon entropies (Shannon.csv), binary strings (included in the same 
   file as the Shannon entropies) and Hamming distances (HammingDistance.csv). 
   Produces a file (PanelA.tif) containing a matrix of functional statuses select 
   genes encoding respiratory chain components in each louse endosymbiont, the 
   Shannon entropies corresponding to the binary strings of these genes, the 
   Hamming distances between the binary strings of these genes and the 
   relationship strengths.
*/
void setup(){
  size(2675,4600);background(255);
  translate(0,25);
  color intact=color(166,219,160);
  color nonIntact=color(153,112,171);
  color direct=color(25.5,84.15,255);
  color reciprocal=color(229.5,43.35,0);
  Table BinStrings=loadTable("../sharedResources/Shannon.csv","header");
  Table NameList=loadTable("../sharedResources/namelist.csv","header");
  String[] removeNames={"Coang","5_11_15_4"};
  String[] modName={"Cotas"};
  NameList.removeRow(NameList.findRowIndex(removeNames[0],"abb"));
  NameList.removeRow(NameList.findRowIndex(removeNames[1],"abb"));
  String nameMod=NameList.getString(NameList.findRowIndex(modName[0],"abb"),"name");
  NameList.setString(NameList.findRowIndex(modName[0],"abb"),"name",nameMod.substring(0,nameMod.length()-2));
  PFont name=createFont("Arial-ItalicMT",75);
  PFont gName=createFont("Arial-ItalicMT",85);
  PFont other=createFont("ArialMT",75);
  PFont bin=createFont("ArialMT",75);
  fill(0);
  String atp=BinStrings.getString(BinStrings.findRowIndex("Sant_0004","tag"),"noninverted");
  String nuo=BinStrings.getString(BinStrings.findRowIndex("Sant_1352","tag"),"noninverted");
  String ndh=BinStrings.getString(BinStrings.findRowIndex("Sant_2463","tag"),"noninverted");
  String cyo=BinStrings.getString(BinStrings.findRowIndex("Sant_3005","tag"),"noninverted");
  String cyd=BinStrings.getString(BinStrings.findRowIndex("Sant_2767","tag"),"noninverted");
  int ndh_nuo=34,nuo_atp=7,atp_cyo=8,cyo_cyd=34;
  float ndhShan=0.977418,nuoShan=0.977418,atpShan=0.833765,cydShan=0.430552,cyoShan=0.430552;
  int nameEnd=1000;int mod=80;int top=125;
  int posAdd=350;
  int ndhPos=nameEnd+150,nuoPos=ndhPos+posAdd,atpPos=nuoPos+posAdd,cyoPos=atpPos+posAdd,cydPos=cyoPos+posAdd;
  String nuoAtpS=nfc(scaleShan(min(atpShan,nuoShan),0.5)+scaleHam(sqrt(sq(nuo_atp-17)),0.5),2);
  String ndhNuoS=nfc(scaleShan(min(ndhShan,nuoShan),0.5)+scaleHam(sqrt(sq(ndh_nuo-17)),0.5),2);
  String cyoCydS=nfc(scaleShan(min(cyoShan,cydShan),0.5)+scaleHam(sqrt(sq(cyo_cyd-17)),0.5),2);
  String atpCyoS=nfc(scaleShan(min(atpShan,cyoShan),0.5)+scaleHam(sqrt(sq(atp_cyo-17)),0.5),2);
  strokeWeight(5);
  int heiSE=top+34*mod-20,heiHD=top+35*mod-20,heiRS=top+36*mod-10;
  line(nameEnd+25,top-50,nameEnd+25,top+36*mod+30);
  line(25,top+33*mod+30,cydPos+100,top+33*mod+30);
  line(25,top+35*mod+30,cydPos+100,top+35*mod+30);
  line(nameEnd+25,top-50,cydPos+100,top-50);
  textFont(gName);textAlign(CENTER,CENTER);
  text("ndh",ndhPos,top-114);
  text("nuo",nuoPos,top-114);
  text("atp",atpPos,top-114);
  text("cyo",cyoPos,top-114);
  text("cyd",cydPos,top-114);
  textFont(other);
  textAlign(RIGHT,CENTER);
  text("Shannon Entropy",nameEnd,heiSE);
  text("Hamming Distance",nameEnd,heiHD);
  text("Relationship Strength",nameEnd,heiRS);
  textAlign(CENTER,CENTER);
  text(nfc(ndhShan,2),ndhPos,heiSE);
  text(nfc(nuoShan,2),nuoPos,heiSE);
  text(nfc(atpShan,2),atpPos,heiSE);
  text(nfc(cydShan,2),cyoPos,heiSE);
  text(nfc(cydShan,2),cydPos,heiSE);
  text(ndh_nuo,((nuoPos+ndhPos)/2),heiHD);
  text(nuo_atp,((nuoPos+atpPos)/2),heiHD);
  text(atp_cyo,((atpPos+cyoPos)/2),heiHD);
  text(cyo_cyd,((cyoPos+cydPos)/2),heiHD);
  fill(reciprocal);
  text(ndhNuoS,((nuoPos+ndhPos)/2),heiRS);
  fill(direct);
  text(nuoAtpS,((nuoPos+atpPos)/2),heiRS);
  fill(direct);
  text(atpCyoS,((atpPos+cyoPos)/2),heiRS);
  fill(reciprocal);
  text(cyoCydS,((cyoPos+cydPos)/2),heiRS);
  fill(0);
  strokeWeight(5);
  for(int i=0;i<34;i++){
    textFont(name);textAlign(RIGHT,CENTER);
    text(NameList.getString(i,"name"),nameEnd,top+i*mod-20);
    textFont(bin);textAlign(CENTER,CENTER);
    if(ndh.charAt(i)=='1'){fill(intact);}else{fill(nonIntact);}
    rect(ndhPos-50,top+i*mod-50,100,mod);
    fill(0);
    if(nuo.charAt(i)=='1'){fill(intact);}else{fill(nonIntact);}
    rect(nuoPos-50,top+i*mod-50,100,mod);
    if(atp.charAt(i)=='1'){fill(intact);}else{fill(nonIntact);}
    rect(atpPos-50,top+i*mod-50,100,mod);
    if(cyo.charAt(i)=='1'){fill(intact);}else{fill(nonIntact);}
    rect(cyoPos-50,top+i*mod-50,100,mod);
    if(cyd.charAt(i)=='1'){fill(intact);}else{fill(nonIntact);}
    rect(cydPos-50,top+i*mod-50,100,mod);
    fill(0);
    text(ndh.charAt(i),ndhPos,top+i*mod-20);
    text(nuo.charAt(i),nuoPos,top+i*mod-20);
    text(atp.charAt(i),atpPos,top+i*mod-20);
    text(cyo.charAt(i),cyoPos,top+i*mod-20);
    text(cyd.charAt(i),cydPos,top+i*mod-20);
  }
  println(ndhPos+100);
  save("PanelA.tif");
  exit();
}
void draw(){}
float scaleHam(int ham,float points){
  return map(ham,0,17,0,points);
}
float scaleHam(float ham,float points){
  return map(ham,0,17,0,points);
}
float scaleShan(float shan,float points){
  return map(shan,0,1,0,points);
}
