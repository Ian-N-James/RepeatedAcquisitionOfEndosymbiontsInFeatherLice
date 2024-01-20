/* Takes the Shannon entropies (Shannon.csv). Produces a plot (PanelB.tif) 
   depicting the frequencies of strings in accordance with string (Shannon) 
   entropy.
*/
void setup(){
  size(2675,1550);
  background(255);
  int marL=400,marR=50,marT=100,marB=275;
  int oriH=height-marB;
  int graH=height-marB-marT;
  int graW=width-marL-marR;
  PFont other=createFont("ArialMT",75);
  textFont(other);
  strokeWeight(5);
  // Load the Shannon entropy values.
  Table shan=loadTable("../sharedResources/Shannon.csv","header");
  int ti=shan.getRowCount();
  // These Shannon entropy values represent the ones possible from a 34 character binary string. 
  float[]shannon={0.191433,0.322757,0.430552,0.522559,0.602431,0.672295,0.733538,0.787127,0.833765,0.873981,0.908178,0.936667,0.959687,0.977418,0.989993,0.997503,1.0};
  int sMax=shannon.length;
  for(int j=0;j<sMax;j++){}
  int[] sCount=new int[sMax];
  float[] sFreq=new float[sMax];
  for(int i=0;i<ti;i++){
    float s=shan.getFloat(i,"Shannon");
    for(int j=0;j<sMax;j++){
      if(float(nfc(s,5))==float(nfc(shannon[j],5))){sCount[j]++;}
    }
  }
  for(int j=0;j<sMax;j++){
    sFreq[j]=float(sCount[j])/ti;
  }
  line(marL,marT,marL,oriH);
  line(marL,oriH,marL+graW,oriH);
  ti=7;
  // Label the y-axis.
  textAlign(RIGHT,CENTER);fill(0);
  for(int i=0;i<ti;i++){
    line(marL,oriH-map(i*0.05,0,0.3,0,graH),marL-20,oriH-map(i*0.05,0,0.3,0,graH));
    text(nfc((i*0.05),2),marL-25,oriH-map(i*0.05,0,0.3,0,graH)-5);
  }
  textAlign(CENTER,CENTER);
  pushMatrix();
  translate(0.175*marL,marT+0.5*oriH);
  rotate(-HALF_PI);
  text("String Frequency",0,0);
  popMatrix();
  float minL=0,maxL=1;
  // Label the x-axis.
  textAlign(CENTER,TOP);
  text("String (Shannon) Entropy",marL+graW/2,oriH+135);
  float[] xLab={0.0,0.2,0.4,0.6,0.8,1.0};
  int ll=xLab.length-1;
  line(marL+map(xLab[0],minL,maxL,0,graW),oriH+50,marL+map(xLab[0],minL,maxL,0,graW),oriH);
  text(nfc(xLab[0],0),marL+map(xLab[0],minL,maxL,0,graW),oriH+50);
  for(int i=1;i<ll;i++){
    line(marL+map(xLab[i],minL,maxL,0,graW),oriH+50,marL+map(xLab[i],minL,maxL,0,graW),oriH);
    text(nfc(xLab[i],1),marL+map(xLab[i],minL,maxL,0,graW),oriH+50);
  }
  line(marL+map(xLab[ll],minL,maxL,0,graW),oriH+50,marL+map(xLab[ll],minL,maxL,0,graW),oriH);
  text(nfc(xLab[ll],1),marL+map(xLab[ll],minL,maxL,0,graW),oriH+50);
  
  //Plot the points.
  for(int j=0;j<sMax;j++){
    ellipse(marL+map(shannon[j],minL,maxL,0,graW),oriH-map(sFreq[j],0,0.3,0,graH),25,25);
  }
  // Save output file.
  save("PanelB.tif");
  exit();
}
void draw(){}
// This function returns the number of zeros in a string.
int numZero(String in){
  int out=0;
  int ti=in.length();
  for(int i=0;i<ti;i++){
    if(in.charAt(i)=='0'){out++;}
  }
  return out;
}
