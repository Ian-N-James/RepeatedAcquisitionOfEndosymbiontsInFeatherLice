int maxHam=34;
void setup(){
  size(5000,5000);background(255);
  
  boolean labelTags=false;
  boolean ex=true;
  int prevW=4500,prevH=4500;
  float hamPoint=127.5,shanPoint=127.5;
  float LB=-((prevW*.9)/2);
  float ellipW=prevW/300;
  float ellipH=prevH/300;
  PFont nam=createFont("Arial-ItalicMT",75);
  PFont namNonGene=createFont("ArialMT",75);
  PFont tag=createFont("ArialMT",12);
  fill(0);
  ellipse(width/2,height/2,prevW*.9,prevH*.9);
  pushMatrix();
  fill(255,0,255);
  translate(width/2,height/2);
  rotate(PI/2.0);
  Table hamming=loadTable("../sharedResources/HammingDistance.csv","header");
  Table shannon=loadTable("../sharedResources/Shannon.csv","header");
  println("tablesLoaded "+millis()/1000+" seconds");
  int ti=shannon.getRowCount();
  boolean[] aboveShanCut=new boolean[ti];
  int numAbove=0;
  int sCut=3;
  IntDict tagID=new IntDict();
  for(int i=0;i<ti;i++){
    if(shanPassFail(shannon.getFloat(i,"Shannon"),sCut)){
      aboveShanCut[i]=true;
      tagID.add(shannon.getString(i,"tag"),numAbove);
      numAbove++;
    }
  }
  println("Shannon's read "+millis()/1000+" seconds");
  int hCut=0;
  IntList[] tarNums=new IntList[numAbove];
  FloatList[] ShannonEntropies=new FloatList[numAbove];
  IntList[] HammingDistances=new IntList[numAbove];
  String[] tags=new String[numAbove];
  int tj=ti+2;
  int ID;
  String[] hamTit=hamming.getColumnTitles();
  for (int i=0;i<numAbove;i++){
    tarNums[i]=new IntList();
    ShannonEntropies[i]=new FloatList();
    HammingDistances[i]=new IntList();
  }
  Table reorg=new Table();float score1=0,score2=0;;int hamHol;float shanHol;
  Table reorg2=new Table();
  reorg.addColumn("id",Table.INT);
  reorg.addColumn("score",Table.FLOAT);
  reorg2.addColumn("id1",Table.INT);
  reorg2.addColumn("id2",Table.INT);
  reorg2.addColumn("typeR",Table.INT);
  reorg2.addColumn("score",Table.FLOAT);
  int n=0;
  for (int i=0;i<ti;i++){
    if (tagID.hasKey(hamming.getString(i,0))){
      ID=tagID.get(hamming.getString(i,0));
      tags[ID]=hamming.getString(i,0);
      reorg.addRow();reorg.setInt(ID,0,ID);score1=0;
      for (int j=2;j<tj;j++){
        if (abs(17-hamming.getInt(i,j))>hCut            //DO require relation to pass Hamming cutoff
          &&tagID.hasKey(hamTit[j])                     //DO require target to pass Shannon cutoff
        && !tarNums[tagID.get(hamTit[j])].hasValue(ID)  //do NOT add relations that already exist
        && !hamming.getString(i,0).equals(hamTit[j])    //do NOT add self relations
        ){
          reorg2.addRow();
          reorg2.setInt(n,"id1",ID);
          reorg2.setInt(n,"id2",tagID.get(hamTit[j]));
          score2=0;
          tarNums[ID].append(tagID.get(hamTit[j]));
          shanHol=float(nfc(shannon.getFloat(i,"Shannon"),4))*float(nfc(shannon.getFloat(j-2,"Shannon"),4));
          if(ex){shanHol=min(float(nfc(shannon.getFloat(i,"Shannon"),4)),float(nfc(shannon.getFloat(j-2,"Shannon"),4)));}
          ShannonEntropies[ID].append(shanHol);
          hamHol=hamming.getInt(i,j);
          HammingDistances[ID].append(hamHol);
          if(hamHol<17){reorg2.setInt(n,"typeR",0);}  
          if(hamHol>17){reorg2.setInt(n,"typeR",1);}
          if(hamHol==17){reorg2.setInt(n,"typeR",2);}
          score2=scaleHam(sqrt(sq(hamHol-17)),hamPoint)+scaleShan(shanHol,shanPoint);
          score1+=score2;
          reorg2.setFloat(n,"score",score2);
          n++;
        }
      }
      reorg.setFloat(ID,"score",score1);
    }
  }
  println("Hammings Read "+millis()/1000+" seconds");
  ti=numAbove;
  float[] corX=new float[ti];
  float[] corY=new float[ti];
  float tranFac=25;
  String[] labels=loadStrings("labels.txt");
  float oti=1.0/float(ti),iota=0,sigh=0,cosigh=0;
  int tl=labels.length;String[] taghold,colHold;
  float RST=0,pls=0;//Ring Seg Thickness
  float mod=0;
  textFont(nam);
  for(int l=0;l<tl;l+=3){
    if(labelTags){RST=125;pls=5;}else{RST=75;pls=15;}
    taghold=labels[l+1].split(",");colHold=labels[l+2].split(";");
    mod=float(labels[l].split(";")[1].split(",")[0]);
    fill(float(colHold[0].split(",")[0]),float(colHold[0].split(",")[1]),float(colHold[0].split(",")[2]),float(colHold[0].split(",")[3]));
    stroke(float(colHold[1].split(",")[0]),float(colHold[1].split(",")[1]),float(colHold[0].split(",")[2]),float(colHold[1].split(",")[3]));
    ringSeg(0,0,LB,LB,LB-RST-mod,LB-RST-mod,(float(tagID.get(taghold[0]))*oti)*TWO_PI,(float(tagID.get(taghold[1]))*oti)*TWO_PI);
    fill(float(colHold[0].split(",")[0])/3,float(colHold[0].split(",")[1])/3,float(colHold[0].split(",")[2])/3,255);
    Float mid=(float(tagID.get(taghold[0]))+float(tagID.get(taghold[1])))/2;
    if(mid/ti>=0.0/8.0 && mid/ti<1.0/8.0){println(labels[l]+" in 1");textAlign(CENTER,BOTTOM);}  // in  
    if(mid/ti>=1.0/8.0 && mid/ti<2.0/8.0){println(labels[l]+" in 2");textAlign(LEFT,CENTER);}  // in 
    if(mid/ti>=2.0/8.0 && mid/ti<3.0/8.0){println(labels[l]+" in 3");textAlign(LEFT,CENTER);}   // in 
    if(mid/ti>=3.0/8.0 && mid/ti<4.0/8.0){println(labels[l]+" in 4");textAlign(RIGHT,CENTER);} // in
    if(mid/ti>=4.0/8.0 && mid/ti<5.0/8.0){println(labels[l]+" in 5");textAlign(RIGHT,TOP);}// in  
    if(mid/ti>=5.0/8.0 && mid/ti<6.0/8.0){println(labels[l]+" in 6");textAlign(RIGHT,CENTER);}  // in 
    if(mid/ti>=6.0/8.0 && mid/ti<7.0/8.0){println(labels[l]+" in 7");textAlign(RIGHT,CENTER);}   // in 
    if(mid/ti>=7.0/8.0 && mid/ti<8.0/8.0){println(labels[l]+" in 8");textAlign(CENTER,BOTTOM);} // in
    float x=((LB-(RST+pls+mod))*cos(mid*oti*TWO_PI))+float(labels[l].split(";")[1].split(",")[1]);
    float y=((LB-(RST+pls+mod))*sin(mid*oti*TWO_PI))+float(labels[l].split(";")[1].split(",")[2]);
    String tex=labels[l].split(";")[0].substring(1).replace("\\n","\n");
    tex=tex.replace("\\t","\t");
    if(containsAny(tex," \t\n")){
      String[][] div=breaker(tex);
      int tk=div[0].length;
      float xCurrent=x,yCurrent=y;
      for(int k=0;k<tk;k++){
        pushMatrix();
        translate(xCurrent,yCurrent);
        rotate(-HALF_PI);
        if(isGeneName(div[0][k])){textFont(nam);}else{textFont(namNonGene);}
        text(div[0][k],0,0);
        popMatrix();
        if(k<tk-1){
          println(div[0][k]+", "+textWidth(div[0][k]));
          if(div[1][k].equals(" ")){yCurrent-=textWidth(div[0][k])+96;}
          //if(div[1][k].equals(" ")){y-=96;}
          if(div[1][k].equals("\t")){yCurrent-=textWidth("  ");}
          if(div[1][k].equals("\n")){xCurrent+=textAscent()+textDescent();yCurrent=y;}
        }
      }
    }else{
      if(isGeneName(tex)){textFont(nam);}else{textFont(namNonGene);}
      pushMatrix();
      translate(x,y);
      rotate(-HALF_PI);
      text(tex,0,0);
      popMatrix();
    }
  }
  textFont(tag);stroke(0);textAlign(LEFT,CENTER);
  if(labelTags){
    for(int i=0;i<ti;i++){
      fill(155,155,155);
      iota=i*oti;
      sigh=sin(iota*TWO_PI);
      cosigh=cos(iota*TWO_PI);
      corX[i]=LB*cosigh;
      corY[i]=LB*sigh;
      ellipse(corX[i],corY[i],ellipW,ellipH);
      line((LB-ellipW/2)*cosigh,(LB-ellipH/2)*sigh,(LB-ellipW*1.5)*cosigh,(LB-ellipH*1.5)*sigh);
      pushMatrix();
      fill(0);
      if(float(i)/float(ti)>0.5){
        rotate(TWO_PI*(float(i)/float(ti)));
        translate((LB-tranFac),0);
        rotate(PI);textAlign(RIGHT,CENTER);
      }else{
        rotate(TWO_PI*(float(i)/float(ti)));
        translate((LB-tranFac),0);
        textAlign(LEFT,CENTER);
      }
      rotate(PI);
      text(tags[i],0,0);
      popMatrix();
    }
  }else{
    for(int i=0;i<ti;i++){
      fill(155,155,155);
      iota=i*oti;
      sigh=sin(iota*TWO_PI);
      cosigh=cos(iota*TWO_PI);
      corX[i]=LB*cosigh;
      corY[i]=LB*sigh;
      ellipse(corX[i],corY[i],ellipW,ellipH);
    }
  }
  println("coordinates generated "+millis()/1000+" seconds");
  reorg2.sort("score");
  int[] n1=reorg2.getIntColumn("id1");
  int[] n2=reorg2.getIntColumn("id2");
  int[] Rtype=reorg2.getIntColumn("typeR");
  float[] scores=reorg2.getFloatColumn("score");
  fill(0,0);
  float x1,x2,x3,y1,y2,y3;
  strokeWeight(2);
  ti=n1.length;
  for(int i=0;i<ti;i++){    
    if(Rtype[i]==0){
      stroke(0,0,scores[i],scores[i]);
      x1=corX[n1[i]];
      y1=corY[n1[i]];
      x2=corX[n2[i]];
      y2=corY[n2[i]];
      x3=(x1+x2+0)/3;
      y3=(y1+y2+0)/3;
      bezier(x1,y1,x3,y3,x3,y3,x2,y2);
    }
    if(Rtype[i]==1){
      stroke(scores[i],0,0,scores[i]);
      x1=corX[n1[i]];
      y1=corY[n1[i]];
      x2=corX[n2[i]];
      y2=corY[n2[i]];
      x3=(x1+x2+0)/3;
      y3=(y1+y2+0)/3;
      bezier(x1,y1,x3,y3,x3,y3,x2,y2);
    }
  }
  println("plotted,Time elapsed: "+millis()/1000+" seconds");
  if(labelTags){save("CirclePlot_Tags.png");}else{save("CirclePlot_NT.png");}
  if(labelTags){save("CirclePlot_Tags.tif");}else{save("CirclePlot_NT.tif");}
  exit();
}



void draw(){
}
boolean shanPassFail(float shan,int cutOff){
  float[]val={0,0.1914333,0.32275695,0.4305519,0.52255934,0.6024308,0.6722948,0.7335379,0.78712654,0.8337649,0.87398106,0.9081783,0.9366674,0.9596869,0.9774178,0.9899928,0.99750245,1.0};
  if (shan>=val[cutOff]){
    return true;
  }
  return false;
}
float scaleHam(int ham,float points){
  return map(ham,0,17,0,points);
}
float scaleHam(float ham,float points){
  return map(ham,0,17,0,points);
}
float scaleShan(float shan,float points){
  return map(shan,0,1,0,points);
}


class sections{
  float cX;
  float cY;
  float disX;
  float disY;
  int mode=0;
  int N=0;
  
  sections(){}
  sections(float tcx,float tcy,float tdx,float tdy){
    cX=tcx;cY=tcy;disX=tdx;disY=tdy;
  }
  sections(float tcx,float tcy,float tdx,float tdy,int tPC){
    cX=tcx;cY=tcy;disX=tdx;disY=tdy;N=tPC;
  }
  void update(float tcx,float tcy,float tdx,float tdy){
    cX=tcx;cY=tcy;disX=tdx;disY=tdy;
  }
  void changeMode(int tMode){mode=tMode;}
  void lable(int staGene,int stoGene,float thick,color f,color s){
    double q=(1.0/float(N))*TWO_PI;
    double w=staGene*q-(q*.3);
    double e=stoGene*q+(q*.3);
    region(w,e,thick,f,s);
  }
  void region(double start,double stop,float thick,color fill,color stroke){    
    float a=disX+.5*thick;float b=disX-.5*thick;
    float res;
    if(mode!=0){
      if(mode==1){
        a=disX+thick;b=disX;
      }
      if(mode==2){
        a=disX;b=disX-thick;
      }
    }
    res=5000.0;
    int sa=round((float)(res*(start/PI)));
    int so=round((float)(res*(stop/PI)));
    PShape reg=createShape();
    reg.beginShape();
    reg.fill(fill);
    reg.stroke(stroke);
    
    for(int i=sa;i<=so;i++){reg.vertex(a*cos(float(i)*PI/res)-cX,cY+a*sin(float(i)*PI/res));}
    for(int i=so;i>=sa;i--){reg.vertex(b*cos(float(i)*PI/res)-cX,cY+b*sin(float(i)*PI/res));}
    
    reg.endShape(CLOSE);
    shape(reg);
  }
  
}
void ringSeg(float xPos,float yPos,float inW,float inH,float outW,float outH,float startAngle,float stopAngle){
  startAngle-=0.004;stopAngle+=0.004;
  strokeJoin(ROUND);
  //float midAngle=lerp(startAngle,stopAngle,0.5);
  curveTightness(0);
  beginShape();
  curveVertex(inW*cos(startAngle)+xPos,inH*sin(startAngle)+yPos);
  curveVertex(inW*cos(startAngle)+xPos,inH*sin(startAngle)+yPos);
  for(float i=startAngle;i<stopAngle;i+=0.05){
    curveVertex(inW*cos(i)+xPos,inH*sin(i)+yPos);
  }
  curveTightness(1);
  curveVertex(inW*cos(stopAngle)+xPos,inH*sin(stopAngle)+yPos);
  curveVertex(lerp(inW,outW,0.5)*cos(stopAngle)+xPos,lerp(inW,outW,0.5)*sin(stopAngle)+yPos);
  
  curveVertex(outW*cos(stopAngle)+xPos,outH*sin(stopAngle)+yPos);
  curveTightness(1);
  curveVertex(outW*cos(stopAngle)+xPos,outH*sin(stopAngle)+yPos);
  
  for(float i=stopAngle-.05;i>startAngle;i-=0.05){
  
    curveVertex(outW*cos(i)+xPos,outH*sin(i)+yPos);
    curveTightness(0);
  }curveTightness(1);
  curveVertex(outW*cos(startAngle)+xPos,outH*sin(startAngle)+yPos);
  curveVertex(outW*cos(startAngle)+xPos,outH*sin(startAngle)+yPos);
  endShape(CLOSE);
}
boolean isGeneName(String in){
  if(!doesNotContain(in)){return false;}
  if(in.length()>=3){
    if(in.length()==3 && in.equals(in.toLowerCase())){return true;}
    if(in.substring(0,3).equals(in.substring(0,3).toLowerCase()) &&in.substring(3).equals(in.substring(3).toUpperCase())){return true;}
  }
  return false;
}
boolean doesNotContain(String test){
  String nonAllowed=" \n\t~!@#$%^&*()_+-=,./<>?;:\'\"\\[]{}|Ω≈ç√∫˜≤≥÷åß∂ƒ©˙∆˚¬…æœ∑´†¥¨ˆπ«¡™£¢∞§¶•ªº–≠¸˛Ç◊ı˜Â¯˘¿ÅÍÎÏ˝ÓÔÒÚÆ`⁄€‹›ﬁﬂ‡°·‚—±Œ„„´‰ˇÁ¨ˆØ∏»";
  int ti=nonAllowed.length();
  for(int i=0;i<ti;i++){
    if(test.contains(str(nonAllowed.charAt(i)))){return false;}
  }
  return true;
}
boolean containsAny(String test,String list){
  int ti=list.length();
  for(int i=0;i<ti;i++){
    if(test.contains(str(list.charAt(i)))){return true;}
  }
  return false;
}
String[][] breaker(String in){
  String[][] out={{},{}};
  int ti=in.length();String temp="";
  for(int i=0;i<ti;i++){
    if(in.charAt(i)==' '||in.charAt(i)=='\n'||in.charAt(i)=='\t'){
      out[0]=append(out[0],temp);
      out[1]=append(out[1],str(in.charAt(i)));
      temp="";
    }else{
      temp+=in.charAt(i);
    }
  }
  out[0]=append(out[0],temp);
  return out;
} 
