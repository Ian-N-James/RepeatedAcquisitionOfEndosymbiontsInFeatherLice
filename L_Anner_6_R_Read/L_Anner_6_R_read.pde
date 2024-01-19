/* Takes the EMdata.csv table and appends the entries in the namelist table with 
   snLED values estimated by E-M analysis to represent cutoffs for intact genes.
*/

void setup(){
  // Load the namelist file.
  Table namelist=loadTable("../sharedResources/namelist.csv","header");
  // Load the EMdata.csv input file. 
  Table rDat=loadTable("../sharedResources/EMdata.csv","header");
  // Compute cutoff values based on the results of the E-M analysis.
  int ti=rDat.getRowCount();
  for(int i=0;i<ti;i++){
    String abb=convert(rDat.getString(i,"name"));
    if(namelist.findRowIndex(abb,"abb")!=-1){
      namelist.setFloat(namelist.findRowIndex(abb,"abb"),"CutOff",rDat.getFloat(i,"cutoff"));
    }else{
      println(abb+" not found in namelist");
    }
  }
  // Save the cutoff values to the namelist table.
  saveTable(namelist,"../sharedResources/namelist.csv");
  exit();
}
void draw(){}

// This function converts the abbreviations back from the modified form that was required in R.
String convert(String in){
  int ti=in.length();
  String[] con={"ZERO","ONE","TWO","THREE","FOUR","FIVE","SIX","SEVEN","EIGHT","NINE"};
  for(int h=0;h<9;h++){
    int conLen=con[h].length();
    for(int i=0;i<ti-conLen+1;i++){
      if(in.substring(i,i+conLen).equals(con[h])){
        in=in.substring(0,i)+h+in.substring(i+conLen);
        i=0;
        ti-=conLen-1;
      }
    }
  }
  return in;
}
