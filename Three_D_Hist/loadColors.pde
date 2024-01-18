color[][] loadColorsHex(String in){
  color[][] output={};
  String[] input=loadStrings(in);
  int ti=input.length;
  for(int i=0;i<ti;i++){
    output=(color[][])append(output, new color[0]);
    String holdinput=input[i];
    holdinput=join(split(holdinput,"["),"");
    holdinput=join(split(holdinput,"]"),"");
    holdinput=join(split(holdinput,"'"),"");
    holdinput=join(split(holdinput,"#"),"");
    String[] holdColors=split(holdinput,", ");
    int tj=holdColors.length;
    for(int j=0;j<tj;j++){
      color temp=color(unhex(holdColors[j]));
      float r=red(temp);
      float g=green(temp);
      float b=blue(temp);
      output[i]=append(output[i],color(r,g,b));
    }
  }
  return output;
}
