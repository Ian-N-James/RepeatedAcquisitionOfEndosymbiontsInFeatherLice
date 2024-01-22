/* C_Map - Takes the Shannon entropies (Shannon.csv), Hamming distances 
   (HammingDistance.csv), a file containing labeling instructions 
   (labels.txt), and two color maps° (linear_ternary-blue_0-44_c57_n256.csv 
   and linear_ternary-red_0-50_c52_n256.csv). Renders a circular figure 
   (out_C_BinaryColMap_NT.tif) with bezier curves depicting the strength of
   contingent relationships between all genes in the dataset. To improve 
   clarity in visualization, genes whose binary strings have Shannon 
   entropies less than 0.43 are excluded.
   °Color maps used for rendering the relationships are from CET 
   Perceptually Uniform Colour Maps (https://colorcet.com/index.html) 
   by Peter Kovesi (See also Good Colour Maps: How to Design Them 2015) and 
   are under a Creative Commons BY License 
   (https://creativecommons.org/licenses/by/4.0/).
*/
int maxHammingDistance=-1;
float centerHammingValue=-1;
boolean testColorLoading=false,revCol=false;
boolean directOnly=false,recipOnly=false;
String namMod=""; 


void setup(){
  size(5000,5000);background(255);
  String colsNam="BinaryColMap";
  // Load the color maps.
  String colsNamD="linear_ternary-blue_0-44_c57_n256";
  float[][] colsD=loadColors("colorProfiles/"+colsNamD+".csv",revCol);
  String colsNamR="linear_ternary-red_0-50_c52_n256";
  float[][] colsR=loadColors("colorProfiles/"+colsNamR+".csv",revCol);
  
  boolean labelTags=false;
  boolean ex=true;
  int prevW=4500,prevH=4500;
  float hamPoint=127.5,shanPoint=127.5;
  
  String num="C";
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
  // Load the Hamming distances.
  Table hamming=loadTable("../sharedResources/HammingDistance.csv","header");
  // Load the Shannon entropies.
  Table shannon=loadTable("../sharedResources/Shannon.csv","header");
  
  maxHammingDistance=shannon.getString(0,"noninverted").length();
  centerHammingValue=maxHammingDistance/2.0; 
  
  
  if(directOnly && recipOnly){directOnly=false;recipOnly=false;}
  println("tablesLoaded "+millis()/1000+" seconds");
  
  /* Exclude genes with a Shannon entropy less than 0.43, corresponding to 
     the exclusion of genes having less than three cases of retention or  
     loss. This value was selected for this system on the basis that genes 
     with two or fewer were over represented in the data set. 
  */   
  int numAbove=0;
  float sCut=0.43;
  IntDict tagID=new IntDict();
  int ti=shannon.getRowCount();
  boolean[] aboveShanCut=new boolean[ti];
  for(int i=0;i<ti;i++){
    if(shannon.getFloat(i,"Shannon")>=sCut){
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
  String[] HammingTitles=hamming.getColumnTitles();
  for (int i=0;i<numAbove;i++){
    tarNums[i]=new IntList();
    ShannonEntropies[i]=new FloatList();
    HammingDistances[i]=new IntList();
  }
  
  // Prepare to reorder the relationships based on their respective relationship strength scores. 
  float score=0;
  int holdHD;float holdSE;
  Table reorg=new Table();
  reorg.addColumn("id1",Table.INT);
  reorg.addColumn("id2",Table.INT);
  reorg.addColumn("typeR",Table.INT);
  reorg.addColumn("score",Table.FLOAT);
  int n=0;
  for (int i=0;i<ti;i++){
    if (tagID.hasKey(hamming.getString(i,0))){
      ID=tagID.get(hamming.getString(i,0));
      tags[ID]=hamming.getString(i,0);
      for (int j=2;j<tj;j++){ 
        if (abs(centerHammingValue-hamming.getInt(i,j))>hCut
          &&tagID.hasKey(HammingTitles[j])                     
        && !tarNums[tagID.get(HammingTitles[j])].hasValue(ID)  
        && !hamming.getString(i,0).equals(HammingTitles[j])    
        ){
          reorg.addRow();
          reorg.setInt(n,"id1",ID);
          reorg.setInt(n,"id2",tagID.get(HammingTitles[j]));
          score=0;
          tarNums[ID].append(tagID.get(HammingTitles[j]));
          holdSE=float(nfc(shannon.getFloat(i,"Shannon"),4))*float(nfc(shannon.getFloat(j-2,"Shannon"),4));
          if(ex){holdSE=min(float(nfc(shannon.getFloat(i,"Shannon"),4)),float(nfc(shannon.getFloat(j-2,"Shannon"),4)));}
          ShannonEntropies[ID].append(holdSE);
          holdHD=hamming.getInt(i,j);
          HammingDistances[ID].append(holdHD);
          if(holdHD<centerHammingValue){reorg.setInt(n,"typeR",0);}  
          if(holdHD>centerHammingValue){reorg.setInt(n,"typeR",1);}
          if(holdHD==centerHammingValue){reorg.setInt(n,"typeR",2);}
          score=scaleHam(sqrt(sq(holdHD-centerHammingValue)),hamPoint)+scaleShan(holdSE,shanPoint);
          reorg.setFloat(n,"score",score);
          n++;
        }
      }
    }
  }
  println("Hammings Read "+millis()/1000+" seconds");
  ti=numAbove;
  float[] corX=new float[ti];
  float[] corY=new float[ti];
  float tranFac=25;
  
  // Load labeling instructions.
  String[] labels=loadStrings("labels"+namMod+".txt");
  float oti=1.0/float(ti),iota=0,sigh=0,cosigh=0;
  int tl=labels.length;String[] taghold,colHold;
  float RST=0;
  float pls=0;
  float mod=0;
  textFont(nam);
  
  // Place the labels in accordance with the instructions.
  for(int l=0;l<tl;l+=3){
    if(labelTags){RST=125;pls=5;}else{RST=75;pls=15;}
    taghold=labels[l+1].split(",");colHold=labels[l+2].split(";");
    mod=float(labels[l].split(";")[1].split(",")[0]);
    fill(float(colHold[0].split(",")[0]),float(colHold[0].split(",")[1]),float(colHold[0].split(",")[2]),float(colHold[0].split(",")[3]));
    stroke(float(colHold[1].split(",")[0]),float(colHold[1].split(",")[1]),float(colHold[0].split(",")[2]),float(colHold[1].split(",")[3]));
    ringSeg(0,0,LB,LB,LB-RST-mod,LB-RST-mod,(float(tagID.get(taghold[0]))*oti)*TWO_PI,(float(tagID.get(taghold[1]))*oti)*TWO_PI);
    fill(float(colHold[0].split(",")[0])/3,float(colHold[0].split(",")[1])/3,float(colHold[0].split(",")[2])/3,255);
    float mid=(float(tagID.get(taghold[0]))+float(tagID.get(taghold[1])))/2;
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
  
  // Generate coordinates for the ends of the bezier curves in advance. 
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
  
  // Reorder the relationships based on relationship strength.
  reorg.sort("score");
  int[] n1=reorg.getIntColumn("id1");
  int[] n2=reorg.getIntColumn("id2");
  int[] Rtype=reorg.getIntColumn("typeR");
  float[] scores=reorg.getFloatColumn("score");
  // Plot the bezier curves.
  fill(0,0);
  float x1,x2,x3,y1,y2,y3;
  strokeWeight(2);
  ti=n1.length;
  int c,f;float t,m;
  // Colors for each arc are derived from the colormaps, using linear interpolation between two colors in cases where the relationship strength value does not result in an exact match.
  if(directOnly){
    for(int i=0;i<ti;i++){
      if(Rtype[i]==0){
        t=scores[i];
        c=ceil(t);f=floor(t);m=t%1;
        stroke(lerp(colsD[0][f],colsD[0][c],m),lerp(colsD[1][f],colsD[1][c],m),lerp(colsD[2][f],colsD[2][c],m),scores[i]);
        x1=corX[n1[i]];
        y1=corY[n1[i]];
        x2=corX[n2[i]];
        y2=corY[n2[i]];
        x3=(x1+x2+0)/3;
        y3=(y1+y2+0)/3;
        bezier(x1,y1,x3,y3,x3,y3,x2,y2);
      }
    }
  }else{
    if(recipOnly){
      for(int i=0;i<ti;i++){
        if(Rtype[i]==1){
          t=scores[i];
          c=ceil(t);f=floor(t);m=t%1;
          stroke(lerp(colsR[0][f],colsR[0][c],m),lerp(colsR[1][f],colsR[1][c],m),lerp(colsR[2][f],colsR[2][c],m),scores[i]);
          x1=corX[n1[i]];
          y1=corY[n1[i]];
          x2=corX[n2[i]];
          y2=corY[n2[i]];
          x3=(x1+x2+0)/3;
          y3=(y1+y2+0)/3;
          bezier(x1,y1,x3,y3,x3,y3,x2,y2);
        }
      }
    }else{
      for(int i=0;i<ti;i++){
        if(Rtype[i]==0){
          t=scores[i];
          c=ceil(t);f=floor(t);m=t%1;
          stroke(lerp(colsD[0][f],colsD[0][c],m),lerp(colsD[1][f],colsD[1][c],m),lerp(colsD[2][f],colsD[2][c],m),scores[i]);
          x1=corX[n1[i]];
          y1=corY[n1[i]];
          x2=corX[n2[i]];
          y2=corY[n2[i]];
          x3=(x1+x2+0)/3;
          y3=(y1+y2+0)/3;
          bezier(x1,y1,x3,y3,x3,y3,x2,y2);
        }
        if(Rtype[i]==1){
          t=scores[i];
          c=ceil(t);f=floor(t);m=t%1;
          stroke(lerp(colsR[0][f],colsR[0][c],m),lerp(colsR[1][f],colsR[1][c],m),lerp(colsR[2][f],colsR[2][c],m),scores[i]);
          x1=corX[n1[i]];
          y1=corY[n1[i]];
          x2=corX[n2[i]];
          y2=corY[n2[i]];
          x3=(x1+x2+0)/3;
          y3=(y1+y2+0)/3;
          bezier(x1,y1,x3,y3,x3,y3,x2,y2);
        }
      }
    }
  }
  if(revCol){colsNam+="_R";}
  namMod+="_"+num+"_"+colsNam;
  if(directOnly){namMod+="_directOnly";}
  if(recipOnly){namMod+="_reciprocalOnly";}
  println("plotted,Time elapsed: "+millis()/1000+" seconds");
  println("testout4"+namMod+"_Tags.png");
  // Save the output file as both a .tif and a .png. 
  if(labelTags){save("outputs/out"+namMod+"_Tags.png");}else{save("outputs/out"+namMod+"_NT.png");}
  if(labelTags){save("outputs/out"+namMod+"_Tags.tif");}else{save("outputs/out"+namMod+"_NT.tif");}
  exit();
}



void draw(){}
// This function returns the component of the relationship strength value corresponding to the Hamming distance.
float scaleHam(int ham,float points){
  return map(ham,0,centerHammingValue,0,points);
}
// This function returns the component of the relationship strength value corresponding to the Hamming distance.
float scaleHam(float ham,float points){
  return map(ham,0,centerHammingValue,0,points);
}
// This function returns the component of the relationship strength value corresponding to the Shannon entropy.
float scaleShan(float shan,float points){
  return map(shan,0,1,0,points);
}

// This function draws the tab component of the labels.
void ringSeg(float xPos,float yPos,float inW,float inH,float outW,float outH,float startAngle,float stopAngle){
  startAngle-=0.004;stopAngle+=0.004;
  strokeJoin(ROUND);
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
// This function determines if text in the labels is italicized.
boolean isGeneName(String in){
  if(!doesNotContain(in)){return false;}
  if(in.length()>=3){
    if(in.length()==3 && in.equals(in.toLowerCase())){return true;}
    if(in.substring(0,3).equals(in.substring(0,3).toLowerCase()) &&in.substring(3).equals(in.substring(3).toUpperCase())){return true;}
  }
  return false;
}
// This function returns true if the string does not have any of the nonAllowed characters.
boolean doesNotContain(String test){
  String nonAllowed=" \n\t~!@#$%^&*()_+-=,./<>?;:\'\"\\[]{}|Ω≈ç√∫˜≤≥÷åß∂ƒ©˙∆˚¬…æœ∑´†¥¨ˆπ«¡™£¢∞§¶•ªº–≠¸˛Ç◊ı˜Â¯˘¿ÅÍÎÏ˝ÓÔÒÚÆ`⁄€‹›ﬁﬂ‡°·‚—±Œ„„´‰ˇÁ¨ˆØ∏»";
  int ti=nonAllowed.length();
  for(int i=0;i<ti;i++){
    if(test.contains(str(nonAllowed.charAt(i)))){return false;}
  }
  return true;
}
// This function returns true if the test string inclusively contains any of the characters in the list string.
boolean containsAny(String test,String list){
  int ti=list.length();
  for(int i=0;i<ti;i++){
    if(test.contains(str(list.charAt(i)))){return true;}
  }
  return false;
}
// This function processes label names that have whitespace.
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
