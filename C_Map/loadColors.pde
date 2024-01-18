float[][] loadColors(String path, boolean reverse){
  float[][] out={{},{},{}};
  String[] in=loadStrings(path);
  int ti=in.length;
  float[] temp;
  String delimit="";
  boolean CSV=false;
  if(path.substring(path.length()-4).equals(".csv")){CSV=true;}
  if(CSV){delimit=",";}else{delimit=" ";}
  for(int i=0;i<ti;i++){
    if(!in[i].equals("")){
      temp=float(in[i].split(delimit));
      out[0]=append(out[0],temp[0]*255);
      out[1]=append(out[1],temp[1]*255);
      out[2]=append(out[2],temp[2]*255);
    }
  }
  if(reverse){
    out[0]=reverse(out[0]);
    out[1]=reverse(out[1]);
    out[2]=reverse(out[2]);
  }
  
  if(testColorLoading){// THis 
    PGraphics testOut=createGraphics(256,10);
    String nam=path.split("/")[1];
    nam=nam.substring(0,nam.lastIndexOf('.'));
    testOut.beginDraw();
    for(int i=0;i<256;i++){
      testOut.stroke(out[0][i],out[1][i],out[2][i]);
      testOut.line(i,0,i,10);
    }
    testOut.endDraw();
    testOut.save("colorTests/testout-"+nam+".png");
    println("colTest complete");
  }
  
  return out;
}
