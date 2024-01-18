/* Takes specialized output files from G_Render saved within its sketch directory 
   (tab.tif and geneRender.tif) and a table of gene annotations (Tags.csv). 
   Serves as an interactive data visualization tool that facilitates manual 
   inspection of the whole dataset and includes scrolling and selection functions.
*/
int botX=0; int botY=0; String Id="";int gene=0;String name1="[No Gene Selected]";String name2="[No Gene Selected]";Table Ann;
int tabWidth=0;
String Func1="[No Gene Selected]";String Func2="[No Gene Selected]";
boolean Shift=false; boolean Control=false; boolean Sant=false;int Genie3;int Genie4;
int genieX=-1;int genieX2=-1;
boolean[] Genomes;
PFont font1,font2;
int[]lowBoun={108,123,138,153,168,183,198,213,228,243,258,273,288,303,318,333,348,363,378,393,408,423,438,453,468,483,498,513,528,543,558,573,588,603,618,633};
int[]higBoun={118,133,148,163,178,193,208,223,238,253,268,283,298,313,328,343,358,373,388,403,418,433,448,463,478,493,508,523,538,553,568,583,598,613,628,643};

void settings(){
  

  Genomes=new boolean[36];
  for(int i=0;i<Genomes.length;i++){
    Genomes[i]=false;
  }
  
  Genie=loadImage("geneRender.tif","tif");
  tab=loadImage("tab.tif","tif");
  tabWidth=tab.width;
  Ann=loadTable("Tags.csv","header");
  if (Genie.height<displayHeight){size(8*displayWidth/9,Genie.height);}else{size(8*displayWidth/9,5*displayHeight/6);}
  
}
boolean firstLoop=true;
PImage Genie; PImage tab; String DummyText1p1="Gene last clicked on = ";String DummyText1p2="N/A";String DummyText1;
String DummyText2p1="Gene the mouse is over = ";String DummyText2p2="";String DummyText2;
void draw(){
  strokeWeight(1);
  if(firstLoop){
  font1 = createFont("Verdana",12);
  font2  = createFont("Verdana-Italic",12);
  firstLoop=false;
  }
  
  textFont(font1);
    if(mouseX>tabWidth&&((((-botX)+mouseX-tabWidth)/4))<=4461){
    Genie3=tabWidth+((mouseX-tabWidth));
    int mod=2;
    if(Genie3%4==0){
      Genie3-=0+mod;
      }else{
        if((Genie3-1)%4==0){
          Genie3-=1+mod;
        }else{
          if((Genie3-2)%4==0){
            Genie3-=2+mod;
          }else{
            if((Genie3-3)%4==0){
              Genie3-=3+mod;
            }
          }}}}
    else{Genie3=-99;}
    if(mouseX>tabWidth &&((((-botX)+mouseX-tabWidth)/4))<=4461){
    genieX2=((((-botX)+(mouseX)-tabWidth)/4));
    Func2=Ann.findRow(str(genieX2),"id").getString("product");
    name2=Ann.findRow(str(genieX2),"id").getString("name");
    DummyText2p2=Ann.findRow(str(genieX2),"id").getString("tag");
    }else{DummyText2p2="N/A";name2="[No Gene Selected]";Func2="[No Gene Selected]";}

  DummyText1=DummyText1p1+DummyText1p2;
  DummyText2=DummyText2p1+DummyText2p2;
  rect(0,0,8*displayWidth/9,5*displayHeight/6);
  image(Genie,botX+tabWidth,botY);
  image(tab,0,0);
  //if(Shift){rect(0,0,8*displayWidth/9,5*displayHeight/6);}
  stroke(0);
  fill(255);
  rect(5,3,225,17);
  rect(500+230,3,245,17);
  rect(5,23,225,17);
  rect(500+230,23,245,17);
  rect(5,43,480,17);
  rect(495,43,480,17);
  rect(235,3,695-225+20,37);
  rect(621,6,46,31);
  rect(673,6,46,31);
  rect(985,3,265,17);
  rect(985,23,265,17);
  rect(985,43,265,17);
  fill(0);
  textAlign(CENTER,CENTER);
  text("Reset",697,20);
  text("1",645,20);
  textAlign(LEFT,BOTTOM);
  text(DummyText1,10,19);
  text(DummyText2,505+230,19);
  if(!(name1.equals("-")||name1.equals("[No Gene Selected]"))){textFont(font2);}
  text(name1,10,39);
  if(!(name2.equals("-")||name2.equals("[No Gene Selected]"))){textFont(font2);}
  text(name2,505+230,39);
  textFont(font1);
  text(Func1,10,59);
  text(Func2,500,59);
  
  
  stroke(255,155,0);
  line(Genie3-1+3,104,Genie3-1+3,716);
  line(3+Genie3+3,104,3+Genie3+3,716);
  text("Control + \"<-\" = move to start of genome",990,40);
  text("Control + \"->\" = move to end of genome",990,60);
  text("Shift = toggle selection mode",990,20);
  text("Presets:",500,30);
  if(Shift){text("Selection Mode = Genome Selection",240,20);}else{text("Selection Mode = Gene Selection",240,20);}
  
  //text("Selection Mode = Gene Selection",940,62);
  //stroke(155,0,200);strokeWeight(2);
  //for(int i=0;i<36;i++){
  //  if(Genomes[i]){
  //    line(100,lowBoun[i]-3.5,width,lowBoun[i]-3.5);
  //    line(100,higBoun[i]-2,width,higBoun[i]-2);
  //  }
  //}
  boolean anySelect=false;
  for(int i=0;i<Genomes.length && !anySelect;i++){
    if(Genomes[i]){anySelect=true;}
  }
  for(int i=0;i<Genomes.length;i++){
    if(Genomes[i]){
      stroke(155,0,200);
      line(tabWidth,lowBoun[i]-3.5,width,lowBoun[i]-3.5);
      line(tabWidth,higBoun[i]-2,width,higBoun[i]-2);
    }else{
      if(anySelect){
         noStroke();
        fill(255,155);
        rect(tabWidth,lowBoun[i]-3.5,width-100,15);
      }
    }
  }
}
void keyPressed(){
  if(keyCode==LEFT){
      if(Control){
        botX=0;
        Control=false;
      }else{if(botX<(-100)){botX+=100;}else{botX=0;}}
    }
  if(keyCode==RIGHT){
      if(Control){
        botX=(8*displayWidth/9)-(Genie.width+100);
        Control=false;
      }
      else{
        if(botX>((8*displayWidth/9)-(Genie.width+100))){botX-=100;}else{botX=(8*displayWidth/9)-(Genie.width+100);}}
    }
  if(keyCode==SHIFT){if(!Shift){Shift=true;}else{Shift=false;}}
  if(keyCode==CONTROL){Control=true;}
  if(keyCode==BACKSPACE){Shift=false;Control=false;Sant=false;}
  if(key!=CODED){if(key=='s' || key=='S'){Sant=true;}
  }
}
void mousePressed(){
   /*if(mouseX>100&&((((-botX)+mouseX-100)/4))<=4461){Genie4=99+((mouseX-100));
    if(Genie4%4==0){
    }else{
      if((Genie4-1)%4==0){
        Genie4--;
      }else{
      if((Genie4-2)%4==0){Genie4-=2;}else{if((Genie4-3)%4==0){Genie4-=3;}}}}}
    else{Genie4=-99;}*/
  if(!Shift){
    if((-botX)+mouseX>tabWidth&&((((-botX)+mouseX-tabWidth)/4))<=4461){
      genieX=((((-botX)+mouseX-tabWidth)/4));
      Func1=Ann.findRow(str(genieX),"id").getString("product");
      name1=Ann.findRow(str(genieX),"id").getString("name");
      DummyText1p2=Ann.findRow(str(genieX),"id").getString("tag");
      //DummyText1p2="Sant_";
      //if(genieX>4097){
      //  genieX-=4097;
      //  if(genieX<10){
      //    DummyText1p2+="P000"+str(genieX);
      //  }else{
      //    if(genieX<100){
      //      DummyText1p2+="P00"+str(genieX);
      //    }else{
      //      DummyText1p2+="P0"+str(genieX);
      //    }
      //  }
      //}else{
      //  if (genieX<500){
      //    if(genieX==301||genieX==302||genieX==303||genieX==304||genieX==266||genieX==267||genieX==471||genieX==472||genieX==122||genieX==263||genieX==475){
      //      DummyText1p2+="t0"+str(genieX);
      //    }else{
      //      if(genieX==262||genieX==264||genieX==265||genieX==470||genieX==473||genieX==474||genieX==476){
      //        DummyText1p2+="r0"+str(genieX);
      //      }else{
      //        if(genieX<10){
      //          DummyText1p2+="000"+str(genieX);
      //        }else{
      //          if(genieX<100){
      //            DummyText1p2+="00"+str(genieX);
      //          }else{
      //            DummyText1p2+="0"+str(genieX);
      //          }
      //        }
      //      }
      //    }
      //  }else{
      //    if (genieX<1000){
      //      if(genieX==870||genieX==871||genieX==730||genieX==797||genieX==928){
      //        DummyText1p2+="t0"+str(genieX);
      //      }else{
      //        if(genieX==927||genieX==929||genieX==930){
      //          DummyText1p2+="r0"+str(genieX);
      //        }else{DummyText1p2+="0"+str(genieX);
      //        }
      //      }
      //    }else{
      //      if (genieX<1500){
      //        if(genieX==1259||genieX==1260||genieX==1261||genieX==1262||genieX==1264||genieX==1265||genieX==1286||genieX==1390){DummyText1p2+="t"+str(genieX);}else{DummyText1p2+=str(genieX);}
      //      }else{
      //        if (genieX<2000){
      //          if(genieX==1630||genieX==1632||genieX==1661||genieX==1825){DummyText1p2+="t"+str(genieX);}else{DummyText1p2+=str(genieX);}
      //        }else{
      //          if (genieX<2500){
      //            if(genieX==2047||genieX==2048||genieX==2173||genieX==2333){DummyText1p2+="t"+str(genieX);}else{DummyText1p2+=str(genieX);}
      //          }else{
      //            if (genieX<3000){
      //              if(genieX==2795||genieX==2796||genieX==2797||genieX==2567||genieX==2580||genieX==2947||genieX==2798||genieX==2555||genieX==2556||genieX==2557||genieX==2755||genieX==2756||genieX==2800||genieX==2801){DummyText1p2+="t"+str(genieX);}else{DummyText1p2+=str(genieX);}
      //            }else{
      //              if (genieX<3500){
      //                if(genieX==3220||genieX==3221||genieX==3222||genieX==3223||genieX==3436||genieX==3437||genieX==3438||genieX==3456||genieX==3457||genieX==3141||genieX==3156||genieX==3165||genieX==3168){
      //                  DummyText1p2+="t"+str(genieX);
      //                }else{
      //                  if(genieX==3166||genieX==3167||genieX==3169){
      //                    DummyText1p2+="r"+str(genieX);
      //                  }else{DummyText1p2+=str(genieX);
      //                }}
      //              }else{
      //                if(genieX==3931||genieX==3932||genieX==3933||genieX==3934||genieX==3904||genieX==3905||genieX==3944||genieX==3945||genieX==3517||genieX==3526||genieX==3538||genieX==3557||genieX==3629||genieX==3817||genieX==3941||genieX==4090){
      //                  DummyText1p2+="t"+str(genieX);
      //                }else{
      //                  if(genieX==3902||genieX==3903||genieX==3906||genieX==3940||genieX==3942||genieX==3946||genieX==4088||genieX==4089||genieX==4091){
      //                    DummyText1p2+="r"+str(genieX);
      //                  }else{DummyText1p2+=str(genieX);}
      //                }
      //              }
      //            }
      //          }
      //        }
      //      }
      //    }
      //  }
      //}
    }else{DummyText1p2="N/A";name1="[No Gene Selected]";}
  }else{
    
    if(mouseY>=108 && mouseY<=643){
      for(int i=0;i<36;i++){
        if(mouseY>=lowBoun[i]&&mouseY<=higBoun[i]){
          Genomes[i]=!Genomes[i];
        }
      }
    }else{
      
      if(mouseX>673 && mouseX<719 && mouseY>6 && mouseY<37){
        
        for(int i=0;i<36;i++){
          Genomes[i]=false;
        }
      }
      if(mouseX>=621&&mouseX<=667&&mouseY>=6&&mouseY<=37){
        for(int i=0;i<36;i++){Genomes[i]=false;}
        Genomes[14]=true;Genomes[15]=true;Genomes[25]=true;Genomes[28]=true;
      }      
    }
  }
}
