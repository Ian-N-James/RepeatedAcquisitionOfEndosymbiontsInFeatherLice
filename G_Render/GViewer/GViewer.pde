/* Takes specialized output files from G_Render saved within its sketch directory 
   (tab.tif and geneRender.tif) and a table of gene annotations (Tags.csv). 
   Serves as an interactive data visualization tool that facilitates manual 
   inspection of the whole dataset and includes scrolling and selection functions.
*/
int botX=0; int botY=0; String Id="";int gene=0;String name1="[No Gene Selected]";String name2="[No Gene Selected]";Table Ann;
int tabWidth=0;
String Func1="[No Gene Selected]";String Func2="[No Gene Selected]";
boolean Shift=false; boolean Control=false; int Genie3;int Genie4;
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
  // Load the specialized output files from G_Render.
  Genie=loadImage("geneRender.tif","tif");
  tab=loadImage("tab.tif","tif");
  tabWidth=tab.width;
  // Load the table of gene annotations.
  Ann=loadTable("Tags.csv","header");
  // Adjust the window size in accordance with the display resolution.
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
  // Monitor which gene the mouse cursor is over.
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
    }else{DummyText2p2="N/A";name2="[No Gene Selected]";Func2="[No Gene Selected]";
  }
  DummyText1=DummyText1p1+DummyText1p2;
  DummyText2=DummyText2p1+DummyText2p2;
  // Draw GUI components.
  rect(0,0,8*displayWidth/9,5*displayHeight/6);
  image(Genie,botX+tabWidth,botY);
  image(tab,0,0);
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
  
  // Writes the locus tag of the gene last clicked on to the screen.
  text(DummyText1,10,19);
  // Writes the locus tag of the gene that the mouse cursor is over to the screen.
  text(DummyText2,505+230,19);
  if(!(name1.equals("-")||name1.equals("[No Gene Selected]"))){textFont(font2);}
  // Writes the name of the gene last clicked on to the screen.
  text(name1,10,39);
  if(!(name2.equals("-")||name2.equals("[No Gene Selected]"))){textFont(font2);}
  // Writes the name of the gene that the mouse cursor is over to the screen.
  text(name2,505+230,39);
  textFont(font1);
  // Writes the product information of the gene last clicked on to the screen.
  text(Func1,10,59);
  // Writes the product information of the gene that the mouse cursor is over to the screen.
  text(Func2,500,59);
  
  // Draw a pair of orange lines (ending below the image) in acordance with mouse cursor position.
  stroke(255,155,0);
  line(Genie3-1+3,104,Genie3-1+3,716);
  line(3+Genie3+3,104,3+Genie3+3,716);
  // Render some instructions for use on the screen.
  text("Control + \"<-\" = move to start of genome",990,40);
  text("Control + \"->\" = move to end of genome",990,60);
  text("Shift = toggle selection mode",990,20);
  text("Presets:",500,30);
  // Render the current selection mode.
  if(Shift){text("Selection Mode = Genome Selection",240,20);}else{text("Selection Mode = Gene Selection",240,20);}
  // Highlight any selected genomes. 
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
  // Scroll in in the direction of the arrow keys.
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
  // Toggle selection mode.  
  if(keyCode==SHIFT){if(!Shift){Shift=true;}else{Shift=false;}}
  // On the next time either the left or right arrow keys is pressed, jump to the corresponding end.
  if(keyCode==CONTROL){Control=true;}
  // Reset to gene selection mode, and cancel any jumps to one end.
  if(keyCode==BACKSPACE){Shift=false;Control=false;}
  
}
void mousePressed(){
  // Updates the "gene last clicked on" status when in gene selection mode.
  if(!Shift){
    if((-botX)+mouseX>tabWidth&&((((-botX)+mouseX-tabWidth)/4))<=4461){
      genieX=((((-botX)+mouseX-tabWidth)/4));
      Func1=Ann.findRow(str(genieX),"id").getString("product");
      name1=Ann.findRow(str(genieX),"id").getString("name");
      DummyText1p2=Ann.findRow(str(genieX),"id").getString("tag");
    }else{DummyText1p2="N/A";name1="[No Gene Selected]";}
  }else{
    // Updates which genome(s) are selected in genome selection mode.
    if(mouseY>=108 && mouseY<=643){
      for(int i=0;i<36;i++){
        if(mouseY>=lowBoun[i]&&mouseY<=higBoun[i]){
          Genomes[i]=!Genomes[i];
        }
      }
    }else{
      // if in genome selection mode, and the reset button is pressed, clear the selection.
      if(mouseX>673 && mouseX<719 && mouseY>6 && mouseY<37){
        for(int i=0;i<36;i++){
          Genomes[i]=false;
        }
      }
      // There is one provided example preset.
      if(mouseX>=621&&mouseX<=667&&mouseY>=6&&mouseY<=37){
        for(int i=0;i<36;i++){Genomes[i]=false;}
        Genomes[14]=true;Genomes[15]=true;Genomes[25]=true;Genomes[28]=true;
      }      
    }
  }
}
