void setup(){
  Table strings=loadTable("../sharedResources/Strings.csv","header");
  String nameMe;int oneCount;
  strings.addColumn("Shannon");
  for(int i=0;i<strings.getRowCount();i++){
    nameMe=strings.getString(i,"noninverted");
    oneCount=0;
    for(int j=0;j<nameMe.length();j++){
      if(nameMe.charAt(j)=='1'){oneCount++;}
    }
    strings.setString(i,"Shannon",nfc(Shannon(float(oneCount)/float(nameMe.length())),6));
  }
  saveTable(strings,"../sharedResources/Shannon.csv");
  exit();
}
void draw(){}






float Shannon (float p){
  return -1*p*logX(p,2)-(1-p)*logX(1-p,2);
}

float logX(float logof,float base){
  return log(logof)/log(base);
}
