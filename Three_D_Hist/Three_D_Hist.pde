float maxpd,minpd;
IntList exclude=new IntList(0,122,262,263,264,265,266,267,301,302,303,304,344,470,471,472,473,474,475,476,484,523,717,718,728,729,730,773,797,798,802,803,804,805,811,870,871,877,888,927,928,929,930,971,1078,1113,1154,1156,1230,1259,1260,1261,1262,1264,1265,1283,1285,1286,1390,1598,1630,1632,1661,1711,1825,1973,1974,2047,2048,2139,2173,2174,2333,2409,2432,2447,2448,2514,2515,2516,2518,2519,2555,2556,2557,2567,2577,2580,2595,2755,2756,2795,2796,2797,2798,2799,2800,2801,2802,2804,2805,2806,2861,2866,2890,2891,2906,2915,2936,2947,2966,3130,3131,3141,3156,3165,3166,3167,3168,3169,3207,3220,3221,3222,3223,3436,3437,3438,3456,3457,3517,3526,3538,3557,3558,3629,3641,3782,3794,3806,3807,3808,3809,3816,3817,3902,3903,3904,3905,3906,3931,3932,3933,3934,3940,3941,3942,3943,3944,3945,3946,4102,4103,4192,4288,4318,4323,4326,4328,4329,4337,4354,4364,4367,4426);
color sideCol=color(64,238);
color topCol=color(128,238);
void setup(){
  size(2125,1675);
  background(255);
  smooth(8);
  PFont txt=createFont("ArialMT",48);
  PFont bld=createFont("Arial-BoldMT",48);
  Table namelist=loadTable("../sharedResources/namelist.csv","header");
  float[] pd=namelist.getFloatColumn("patDis");
  int syms=pd.length;
  minpd=pd[0];
  maxpd=pd[syms-1];
  String[] names=namelist.getStringColumn("abb");
  float[] cutoff=namelist.getFloatColumn("CutOff");
  int[][] counts=new int[syms][100];
  for(int i=0;i<syms;i++){
    Table zFin=loadTable("../sharedResources/zFin/zFin-"+names[i]+"-AnnTabFin.csv","header");
    int tj=zFin.getRowCount();
    for(int j=0;j<tj;j++){
      if(exclude.hasValue(j)){continue;}
      if(zFin.getFloat(j,"ann")==0){continue;}
      if(zFin.getFloat(j,"weightScore")==1){continue;}
      counts[i][floor(zFin.getFloat(j,"weightScore")*100)]++;
    }
  }
  //set up plot area:
  //original:                            X     Y      
  //bottom left corner, pre Z shift:    355, 1754
  //top left corner, pre Z shift:       355,  954
  //bottom left corner, post Z shift:  1092, 1018
  //top left corner, post Z shift:     1092,  218
  //bottom right corner, pre Z shift:  1453, 1754
  //bottom right corner, post Z shift: 2190, 1018
  //top right corner, post Z shift:    2190,  218
  // plot area                    
  // X:1098, 10.98 pixels per bar. 
  // Y:800, 4 pixels per count.   
  // Z shift:                     
  // ratio 1X:1Y                  
  // 737 pixels each
  float minXprZ=200, minYprZ=754, minXpoZ=937,minYpoZ=18;
  float maxXprZ=1298,maxYprZ=1554,maxXpoZ=2035,maxYpoZ=818;
  float plotWid=maxXprZ-minXprZ, plotHei=maxYprZ-minYprZ;
  float zShiftX=minXpoZ-minXprZ-10/3, zShiftY=minYprZ-minYpoZ-10/3;
  float barWid=plotWid/100;
  strokeWeight(2);stroke(155);
  //horizontal
  line(minXprZ,maxYprZ,maxXprZ,maxYprZ);
  line(minXpoZ,maxYpoZ,maxXpoZ,maxYpoZ);
  line(minXpoZ,minYpoZ,maxXpoZ,minYpoZ);
  //vertical
  line(minXprZ,minYprZ,minXprZ,maxYprZ);
  line(minXpoZ,minYpoZ,minXpoZ,maxYpoZ);
  line(maxXpoZ,minYpoZ,maxXpoZ,maxYpoZ);
  //diagonal
  line(minXprZ,maxYprZ,minXpoZ,maxYpoZ);
  line(minXprZ,minYprZ,minXpoZ,minYpoZ);
  line(maxXprZ,maxYprZ,maxXpoZ,maxYpoZ);
  int divlin=3;
  float dvlnScl=plotHei/(divlin+1);
  for(int h=1;h<=divlin;h++){
    line(minXpoZ,maxYpoZ-h*dvlnScl,maxXpoZ,maxYpoZ-h*dvlnScl);
    line(minXprZ,maxYprZ-h*dvlnScl,minXpoZ,maxYpoZ-h*dvlnScl);
  }
  //AxLabX
  textAlign(CENTER,TOP);fill(0);textFont(txt);strokeWeight(3);
  text("Corrected Levenshtein Edit Distance",(minXprZ+maxXprZ)/2,maxYprZ+40);
  text("1",minXprZ,maxYprZ+25);line(minXprZ,maxYprZ+20,minXprZ,maxYprZ);
  text("0",maxXprZ,maxYprZ+25);line(maxXprZ,maxYprZ+20,maxXprZ,maxYprZ);
  //AxLabY
  textAlign(CENTER,BOTTOM);
  pushMatrix();
  translate(minXprZ-125,(minYprZ+maxYprZ)/2);
  rotate(-HALF_PI);
  text("Count",0,0);
  popMatrix();
  textAlign(RIGHT,CENTER);
  int th=divlin+2;
  text("0",minXprZ-25,maxYprZ-5);
  line(minXprZ-20,maxYprZ,minXprZ,maxYprZ);
  for(int h=1;h<th;h++){
    text(str(h*50),minXprZ-25,maxYprZ-h*dvlnScl-5);
    line(minXprZ-20,maxYprZ-h*dvlnScl-1,minXprZ,maxYprZ-h*dvlnScl-1);
  }
  //AxLabZ
  stroke(0);strokeWeight(8);int adjX=55,adjY=0;
  line(maxXpoZ+adjX,maxYpoZ+adjY,maxXprZ+adjX,maxYprZ+adjY);
  line(maxXprZ+adjX,maxYprZ+adjY,maxXprZ+adjX+10*1.25,maxYprZ+adjY-40*1.25);
  line(maxXprZ+adjX,maxYprZ+adjY,maxXprZ+adjX+40*1.25,maxYprZ+adjY-10*1.25);
  stroke(psF(minpd));
  int psdYmod=110;
  line(maxXpoZ-400,maxYpoZ-2*dvlnScl+10+psdYmod,minXpoZ+50,maxYpoZ-2*dvlnScl+10+psdYmod);
  line(minXpoZ+50,maxYpoZ-2*dvlnScl+10+psdYmod,minXpoZ+50+25*1.25,maxYpoZ-2*dvlnScl+10-15*1.25+psdYmod);
  line(minXpoZ+50,maxYpoZ-2*dvlnScl+10+psdYmod,minXpoZ+50+25*1.25,maxYpoZ-2*dvlnScl+10+15*1.25+psdYmod);
  textAlign(CENTER,BOTTOM);
  text("Pseudogene Degradation",25+(maxXpoZ-400+minXpoZ+50)/2,maxYpoZ-2*dvlnScl+psdYmod);
  stroke(inF(minpd));
  line(maxXpoZ-300,maxYpoZ-645,maxXpoZ-100,maxYpoZ-445);
  line(maxXpoZ-100,maxYpoZ-445,maxXpoZ-100-10,maxYpoZ-445-35);
  line(maxXpoZ-100,maxYpoZ-445,maxXpoZ-100-35,maxYpoZ-445-10);
  text("Intact Genes",maxXpoZ-300,maxYpoZ-645);
  textAlign(CENTER,TOP);
  pushMatrix();
  translate(80+(maxXprZ+maxXpoZ)/2,(maxYprZ+maxYpoZ)/2);
  rotate(radians(-47.5));
  text("Patristic Distance",0,0);
  popMatrix();
  strokeWeight(1);
  for(int i=0;i<syms;i++){
    float xMod=map(pd[i],minpd,maxpd,zShiftX,0);
    float yMod=map(pd[i],minpd,maxpd,zShiftY,0);;
    for(int j=0;j<100;j++){
      if(counts[i][j]>0){
        if((j+1)*0.01<cutoff[i]){
          bar(maxXprZ-(j+1)*barWid+xMod,maxYprZ-yMod,10,-4*counts[i][j],inF(pd[i]),inS(pd[i]));
        }else{
          bar(maxXprZ-(j+1)*barWid+xMod,maxYprZ-yMod,10,-4*counts[i][j],psF(pd[i]),psS(pd[i]));
        }
      }
    }
  }
  save("3DHIST.png");
  save("3DHIST.tif");
  exit();
}
void draw(){}

color inF(float in){
  return color(0,200+map(in,minpd,maxpd,0,55),200+map(in,minpd,maxpd,55,0),238);
}
color inS(float in){
  return color(0,150+map(in,minpd,maxpd,4,28),150+map(in,minpd,maxpd,28,4),238);
}
color psF(float in){
  return color(200+map(in,minpd,maxpd,4,55),0,200+map(in,minpd,maxpd,55,4),238);
}
color psS(float in){
  return color(150+map(in,minpd,maxpd,4,28),0,150+map(in,minpd,maxpd,28,4),238);
}
void bar(float staX,float staY,float staW,float staH,color inColF,color inColS){
  color s=getGraphics().strokeColor;float w=getGraphics().strokeWeight;
  stroke(inColS);strokeWeight(0.1);
  fill(inColF);
  rect(staX,staY,staW,staH);
  stroke(0);strokeWeight(0.1);
  float adjX=staW/3,adjY=staW/3;
  fill(sideCol);
  quad(staX+staW,staY,staX+staW+adjX,staY-adjY,staX+staW+adjX,staY+staH-adjY,staX+staW,staY+staH);
  fill(topCol);
  quad(staX,staY+staH,staX+staW,staY+staH,staX+staW+adjX,staY+staH-adjY,staX+adjX,staY+staH-adjY);
  stroke(s);strokeWeight(w);
}
