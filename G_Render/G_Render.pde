/* Takes the annotation tables (ZFin-*-AnnTabFin.csv) that have been processed 
   by LAnner_7_Predictor and a color map (Gradient.txt). Renders several 
   outputs that show gene status (with intact genes in green, pseudo genes in 
   purple and missing genes in black) across symbiont genomes. One output, 
   saved within the directory of G_Veiwer, has the Tab containing gene names 
   saved as a separate file (tab.tif) from the gene statuses (geneRender.tif). 
   The other outputs are saved within the sketch directory. Two use sub-pixel 
   rendering (GRender.tif and GRender_2.tif) and differ by the amount of 
   whitespace to the right and top of the image. The final output 
   (GRender_1px.tif) does not use sub-pixel rendering, instead rendering at 
   one pixel per gene.
*/

PGraphics tab,ren1,ren1_2,ren2,ren4,tab2;
  int startX=100;     // the initial X start position of the genomes
  int startY=105;     // the initial Y start position of the genomes
  int startY2=10;     // the initial Y start position of the genomes
  int spacer=5;     // how much space should go between genomes
  int spacer2=5;     // how much space should go between genomes
  int genomeNum=36;  // number of genomes
  int Height=10;     // how tall should the line be for a single gene should be
  int Height2=13;     // how tall should the line be for a single gene should be
  float mult=0.4;
void setup(){
  int conLev=30;      // number of
  int GenomeNum=genomeNum+(genomeNum-conLev+1);   // number of genomes

  int Fontsize=10;   // how big the font be for the Genome Names should be
  int geneNum=4462;  // number of genes in the genome including the plasmid
  PFont GenomeNamesFont = createFont("Arial-ItalicMT", 12),sup=createFont("ArialMT",
  GenomeNamesFont.getSize()/2), LabelsFont =createFont("ArialMT", GenomeNamesFont.getSize()),SWAFont=createFont("ArialMT", GenomeNamesFont.getSize());
  PFont GenomeNamesFont2 = createFont("Arial-ItalicMT", 16),sup2=createFont("ArialMT",GenomeNamesFont2.getSize()/2), LabelsFont2 =createFont("ArialMT", GenomeNamesFont2.getSize()),SWAFont2=createFont("ArialMT", GenomeNamesFont2.getSize());
  color intact=color(166,219,160);//for 
  color pseudo=color(153,112,171);//
  color MC=color(100,100,100);
  color CF=color(255,255,0);// this color should not appear.
  color[][] gradients=loadColorsHex("Gradient.txt");
  // gradients[0] is the pseudogene gradient
  // gradients[1] is the intact gene gradient
  color riboRegion=color(53,151,143);
  color phagePlasRegions=color(223,194,125);
  int[][] findCom=new int[geneNum][3];
  
  
  
  
  println((4*geneNum)+startX+10,GenomeNum*(Height+spacer)-spacer+startY+15,"\n");
  
  
  
  
  
  Table Anntab,Namelist,RetenTab;
  //input the genome names
  String[] GenomeNames,GenomeAbbs;float[] Cutoffs;
  Namelist=loadTable("../sharedResources/namelist.csv","header");
  GenomeNames=new String[Namelist.getRowCount()];
  GenomeAbbs=new String[Namelist.getRowCount()];
  Cutoffs=new float[Namelist.getRowCount()];
  for(int i=0;i<Namelist.getRowCount();i++){
    GenomeAbbs[i]=Namelist.getString(i,"abb");
    GenomeNames[i]=Namelist.getString(i,"name");
    Cutoffs[i]=Namelist.getFloat(i,"CutOff");
  }
  genomeNum=GenomeAbbs.length;
  //String[] cutsIn=loadStrings("Cuttoffs.txt");
  Table[] anntabs=new Table[genomeNum];
  int[] tabLens=new int[genomeNum];
  boolean[] remove=new boolean[geneNum];
  for(int i=0;i<genomeNum;i++){//input the gene annotations
    anntabs[i]=loadTable("../sharedResources/zFin/zFin-"+GenomeAbbs[i]+"-AnnTabFin.csv","header");
    for(int j=anntabs[i].getRowCount()-1;j>=0;j--){
      if(anntabs[i].getString(j,"flags").equals("pseuHS")||
       anntabs[i].getString(j,"flags").equals("PpIs1,pseuHS")||
       anntabs[i].getString(j,"flags").equals("PpIs2,pseuHS")||
       anntabs[i].getString(j,"flags").equals("Plas,pseuHS")||
       anntabs[i].getFloat(j,"weightScore")==20||
       anntabs[i].getString(j,"flags").equals("tRNA")||
       anntabs[i].getString(j,"flags").equals("rRNA")||
       anntabs[i].getString(j,"flags").equals("miscRNA")||
       anntabs[i].getString(j,"flags").equals("RPstar")){
         remove[j]=true;
      }
    }
  }
  for(int i=0;i<genomeNum;i++){//input the gene annotations
    for(int j=anntabs[i].getRowCount()-1;j>=0;j--){
      if(remove[j]){anntabs[i].removeRow(j);}
    }
    tabLens[i]=anntabs[i].getRowCount();
  }
  boolean test=true;
  for(int i=1;i<genomeNum;i++){
   if(tabLens[i]!=tabLens[i-1]){
     test=false;
   }
  }
  
  if(!test){println("unequal anntabs lengths");}else{
    geneNum=tabLens[0];  
  }
  int ti=genomeNum, tj=geneNum;
  println(anntabs[0].getRowCount()+", "+tj);
  int[][] geneStatus=new int[ti][tj];
  float[][] geneWScore=new float[ti][tj];
  for(int i=0;i<ti;i++){
    for(int j=0;j<tj;j++){
      
      float statusRead=anntabs[i].getFloat(j,"weightScore");//ann
      int call= -1;//just a number that isn't used elsewere, for validation
      if(statusRead==1){call=0;}                          // Missing
      if(statusRead<1 && statusRead>Cutoffs[i]){call=1;}  // Pseudo
      if(statusRead>=0 && statusRead<=Cutoffs[i]){call=2;} // Intact
      if(statusRead==20.0){call=-2;}                      // Gene not looked at.
      if(call==-1){println("Unidentified call on "+GenomeAbbs[i]+", "+anntabs[i].getString(j,"tag"));}
      if(call==-2){println("Non protien gene / HS pseudogene on "+GenomeAbbs[i]+", "+anntabs[i].getFloat(j,"tag"));}
      geneStatus[i][j]=call;
      geneWScore[i][j]=statusRead;
    }
  }
  tab=createGraphics(145,790);
  tab2=createGraphics(145+132,790);
  ren1=createGraphics(ceil(geneNum*mult)+10,790);
  ren1_2=createGraphics(ceil(geneNum*mult)+10,790);
  ren2=createGraphics(geneNum+10,790);
  ren4=createGraphics(geneNum*4+10,790);
  tab.beginDraw();ren1.beginDraw();ren1_2.beginDraw();ren2.beginDraw();ren4.beginDraw();
  tab2.beginDraw();
  tab.background(255);ren1.background(255);ren1_2.background(255);ren2.background(255);ren4.background(255);
  tab2.background(255);
  
  //render the nonplasmid genes & genome names
  //
  //phage island 1: Sant_2513 - Sant_2554 (155,55,255,55)
  //phage island 2: Sant_2857 - Sant_2946 (155,55,255,55)
  //plasmid: Sant_P0001 - Sant_P0364 (155,55,255,55)
  //ribosomal protiens: sant_0425 (rpsJ) - Sant_0452 (rplQ) (0,55,255,55)
  String[] regionStartTags={"Sant_2513","Sant_2857","Sant_P0001","Sant_0425"};
  String[] regionStopTags={"Sant_2554","Sant_2946","Sant_P0364","Sant_0452"};
  color[]  regionColors={phagePlasRegions,phagePlasRegions,phagePlasRegions,riboRegion}; 
  String[] regionNames={"Prophage Island 1","Prophage Island 2","Plasmid","Ribosomal Proteins"};
  int tk=regionStartTags.length;
  for(int k=0;k<tk;k++){
    // star and stop exist only to improve readability.
    int star=anntabs[0].findRowIndex(regionStartTags[k],"tag");
    int stop=anntabs[0].findRowIndex(regionStopTags[k],"tag");
    println(star+", "+stop+", "+regionNames[k]);
    specialRegion(star,stop,regionColors[k],regionNames[k],LabelsFont,LabelsFont2);
  }

  tab.fill(0);tab2.fill(0);tab2.textFont(GenomeNamesFont2);
  int StartY=startY;int StartY2=startY2;
  PGraphics ren1p,ren1i,ren1e,ren1_2p,ren1_2i,ren1_2e;
  ren1p=createGraphics(ren1.width,ren1.height);ren1p.beginDraw();
  ren1i=createGraphics(ren1.width,ren1.height);ren1i.beginDraw();
  ren1e=createGraphics(ren1.width,ren1.height);ren1e.beginDraw();
  ren1_2p=createGraphics(ren1_2.width,ren1_2.height);ren1_2p.beginDraw();
  ren1_2i=createGraphics(ren1_2.width,ren1_2.height);ren1_2i.beginDraw();
  ren1_2e=createGraphics(ren1_2.width,ren1_2.height);ren1_2e.beginDraw();
  tab.textAlign(RIGHT);
  tab2.textAlign(RIGHT,CENTER);
  float x,m;
  int f,c;
  color useCol;
  for(int i=0;i<ti;i++){
    ren1.strokeWeight(.65);ren4.strokeWeight(.65);ren1_2.strokeWeight(.65);
    if(GenomeNames[i].split("@").length==1){
      tab.textFont(GenomeNamesFont);
      tab.text(GenomeNames[i],tab.width-5,startY+10);
      tab2.pushMatrix();
      tab2.translate(tab2.width-5,startY2+Height/2);
      //tab2.rotate(PI/6);
      tab2.text(GenomeNames[i],0,0);
      tab2.popMatrix();
    }else{
      tab.textFont(sup);
      tab.text(GenomeNames[i].split("@")[1],tab.width-5,startY+10-2*Fontsize/3);
      tab.textFont(GenomeNamesFont);
      tab.text(GenomeNames[i].split("@")[0],tab.width-11,startY+10);
      tab2.textFont(sup2);
      tab2.pushMatrix();
      tab2.translate(tab2.width-5,startY2+Height/2-2*Fontsize/6);
      //tab2.rotate(PI/6);
      tab2.text(GenomeNames[i].split("@")[1],0,0);
      tab2.popMatrix();
      tab2.textFont(GenomeNamesFont2);
      tab2.pushMatrix();
      tab2.translate(tab2.width-11,startY2+Height/2);
      //tab2.rotate(PI/6);
      tab2.text(GenomeNames[i].split("@")[0],0,0);
      tab2.popMatrix();
    }
    ren2.noStroke();
    for(int j=0;j<tj;j++){
      switch(geneStatus[i][j]){
        case 0: 
          ren4.fill(0); ren4.stroke(0);
          ren1.fill(0);ren1.stroke(0);
          ren1_2.fill(0);ren1_2.stroke(0);
          ren2.fill(0);
          ren1.rect((j*mult),startY,mult-mult/4,Height);
          ren1_2.rect((j*mult),startY2,mult-mult/4,Height2);
          findCom[j][1]++;
        break;
        case 1:
          x=map(geneWScore[i][j],Cutoffs[i],1,gradients[0].length-1,0);
          f=floor(x);c=ceil(x);
          m=x%1;
          useCol=lerpColor(gradients[0][f],gradients[0][c],m);
          ren4.fill(useCol);ren4.stroke(useCol);
          ren1p.fill(useCol);ren1p.stroke(useCol);
          ren1_2p.fill(useCol);ren1_2p.stroke(useCol);
          ren2.fill(useCol);
          ren1p.rect((j*mult),startY,mult-mult/4,Height);
          ren1_2p.rect((j*mult),startY2,mult-mult/4,Height2);
          findCom[j][1]++;
          findCom[j][2]++;
        break;
        case 2: 
          x=map(geneWScore[i][j],0,Cutoffs[i],gradients[1].length-1,0);
          f=floor(x);c=ceil(x);
          m=x%1;
          useCol=lerpColor(gradients[1][f],gradients[1][c],m);
          
          ren4.fill(useCol);ren4.stroke(useCol);
          ren1i.fill(useCol);ren1i.stroke(useCol);
          ren1_2i.fill(useCol);ren1_2i.stroke(useCol);
          ren2.fill(useCol);
          ren1i.rect((j*mult),startY,mult-mult/4,Height);
          ren1_2i.rect((j*mult),startY2,mult-mult/4,Height2);
          findCom[j][0]++;
          findCom[j][2]++;
        break;
        case -1: 
          ren4.fill(CF);ren4.stroke(CF);
          ren1e.fill(CF);ren1e.stroke(CF);
          ren1_2e.fill(CF);ren1_2e.stroke(CF);
          ren2.fill(CF);
          ren1e.rect((j*mult),startY,mult-mult/4,Height);
          ren1_2e.rect((j*mult),startY2,mult-mult/4,Height2);
        break;
        case -2:
          ren4.fill(CF);ren4.stroke(CF);
          ren1e.rect((j*mult),startY,mult-mult/4,Height);
          ren1_2e.rect((j*mult),startY2,mult-mult/4,Height2);
          ren2.fill(CF);
        break;
        default: break;
      }
      ren4.rect(j*4,startY,3,Height);
      ren2.rect(j,startY,1,Height);
    }
    ren1i.stroke(0);ren1_2i.stroke(0);ren4.stroke(0);ren2.stroke(0);
    ren1i.strokeWeight(.25);ren1_2i.strokeWeight(.25);ren4.strokeWeight(.25);
    ren1i.line(0,startY,(geneNum*mult),startY);
    ren1_2i.line(0,startY2,(geneNum*mult),startY2);
    ren1i.line(0,startY+Height,(geneNum*mult),startY+Height);
    ren1_2i.line(0,startY2+Height2,(geneNum*mult),startY2+Height2);
    ren4.line(0,startY,geneNum,startY);
    ren4.line(0,startY+Height,geneNum*4,startY+Height);
    ren1i.strokeWeight(.65);ren1_2i.strokeWeight(.65);ren4.strokeWeight(.65);
    startY=startY+Height+spacer;
    startY2=startY2+Height2+spacer2;
    StartY=startY;
    StartY2=startY2;
  }
  ren1p.endDraw();ren1i.endDraw();ren1e.endDraw();
  ren1_2p.endDraw();ren1_2i.endDraw();ren1_2e.endDraw();
  ren1.image(ren1p,0,0);ren1_2.image(ren1_2p,0,0);
  ren1.image(ren1i,0,0);ren1_2.image(ren1_2i,0,0);
  ren1.image(ren1e,0,0);ren1_2.image(ren1e,0,0);
  startY=startY+Height+spacer;
  startY2=startY2+Height2+spacer2;
  ren2.noStroke();
  String[] uniLab={"Universally Intact","Universally Lost","Presence in Any"};
  color[]uniCol={intact,pseudo,MC};
  int[] uniMins={genomeNum,genomeNum,1};
  tab.textFont(SWAFont);
  tab2.textFont(SWAFont2);
  for(int i=0;i<3;i++){
    for(int j=0;j<tj;j++){
      if(findCom[j][i]>=uniMins[i]){
        ren1.stroke(uniCol[i]);ren1.fill(uniCol[i]);ren4.fill(uniCol[i]);
        ren1_2.stroke(uniCol[i]);ren1_2.fill(uniCol[i]);
        ren2.fill(uniCol[i]);
      }else{
        ren1.stroke(0);ren1.fill(0);ren4.fill(0);
        ren1_2.stroke(0);ren1_2.fill(0);
        ren2.fill(0);
        
      }
      ren4.rect(j*4,startY,3,Height);
      ren1.rect((j*mult),startY,1-1/4,Height);
      ren1_2.rect((j*mult),startY2,1-1/4,Height2);
      ren2.rect(j,startY,1,Height);
    }
    tab.text(uniLab[i],tab.width-5,startY+10);
    
    tab2.text(uniLab[i],tab2.width-5,startY2+Height/2);
    startY=startY+Height+spacer;
    startY2=startY2+Height2+spacer2;
  }
  ren1.endDraw();ren1_2.endDraw();ren2.endDraw();ren4.endDraw();tab.endDraw();tab2.endDraw();
  PGraphics file1=createGraphics(tab.width+ren1.width,tab.height);
  file1.beginDraw();
  file1.image(tab,0,0);file1.image(ren1,tab.width,0);
  file1.endDraw();
  file1.save("GRender.tif");
  PGraphics file1_2=createGraphics(tab2.width+ren1_2.width,tab2.height);
  file1_2.beginDraw();
  file1_2.image(tab2,0,0);file1_2.image(ren1_2,tab2.width,0);
  file1_2.endDraw();
  file1_2.save("GRender_2.tif");
  PGraphics file2=createGraphics(tab.width+ren2.width,tab.height);
  file2.beginDraw();
  file2.image(tab,0,0);file2.image(ren2,tab.width,0);
  file2.endDraw();
  file2.save("GRender_1px.tif");
  PGraphics file3=createGraphics(tab2.width+ren2.width,tab.height);
  ren4.save("GVeiwer/geneRender.tif");
  tab.save("GVeiwer/tab.tif");
  exit();
}
void draw(){}

void specialRegion(int start,int stop,color col,String name,PFont fon,PFont fon2){
  int rHeight=(genomeNum+3)*Height+(genomeNum+7)*spacer,tHeight=rHeight+startY;
  int rHeight2=(genomeNum+3)*Height2+(genomeNum+7)*spacer2,tHeight2=rHeight2+startY2;
  color c1f=ren1.fillColor,c1s=ren1.strokeColor,c4f=ren4.fillColor,c4s=ren4.strokeColor;
  //PFont c1t=ren1.textFont,c4t=ren4.textFont; 
  ren1.fill(col);ren1_2.fill(col);ren2.fill(col);ren4.fill(col);
  ren1.stroke(0);ren1_2.stroke(0);ren2.stroke(0);ren4.stroke(0);ren1.strokeWeight(.25);
  ren1.textFont(fon);ren4.textFont(fon);
  ren1_2.textFont(fon2);
  ren1.rect(start*mult,startY-spacer,(stop-start+1)*mult,rHeight);
  ren1_2.rect(start*mult,startY2-spacer2,(stop-start+1)*mult,rHeight2);
  ren4.rect(start*4,startY-5,(stop-start+1)*4,rHeight);
  ren2.rect(start,startY-5,(stop-start+1),rHeight);
  ren1.fill(0);ren1_2.fill(0);ren4.fill(0);ren2.fill(0);
  ren1.textAlign(CENTER,TOP);ren1_2.textAlign(CENTER,TOP);ren2.textAlign(CENTER,TOP);
  ren1.text(name,(start+(stop-start+1)/2)*mult,tHeight);
  ren1_2.text(name,(start+(stop-start+1)/2)*mult,tHeight2);
  println((stop-start));
  ren4.textAlign(CENTER,TOP);
  ren4.text(name,(start+(stop-start+1)/2)*4,tHeight);
  ren2.text(name,(start+(stop-start+1)/2),tHeight);
  // Reset the values of ren1 and ren4 back to what they were
  ren1.fill(c1f);ren1.stroke(c1s);ren1_2.fill(c1f);ren1_2.stroke(c1s);ren4.fill(c4f);ren4.stroke(c4s);
}
