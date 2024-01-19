/* Takes the binary strings (Strings.csv) and computes Shannon entropies for each 
   string. Exports a file (Shannon.csv) that maintains the content of Strings.csv 
   and is appended to include Shannon entropy values.
*/
void setup(){
  // Load the binary strings.
  Table strings=loadTable("../sharedResources/Strings.csv","header");
  String nameMe;int oneCount;
  
  // Compute Shannon entropies.
  strings.addColumn("Shannon");
  for(int i=0;i<strings.getRowCount();i++){
    nameMe=strings.getString(i,"noninverted");
    oneCount=0;
    for(int j=0;j<nameMe.length();j++){
      if(nameMe.charAt(j)=='1'){oneCount++;}
    }
    strings.setString(i,"Shannon",nfc(Shannon(float(oneCount)/float(nameMe.length())),6));
  }
  // Save the output file.
  saveTable(strings,"../sharedResources/Shannon.csv");
  exit();
}
void draw(){}

// This function computes Shannon entropy.
float Shannon (float p){
  return -1*p*logX(p,2)-(1-p)*logX(1-p,2);
}
// This function takes the logarithm of a number at an arbitrary base.
float logX(float logof,float base){
  return log(logof)/log(base);
}
