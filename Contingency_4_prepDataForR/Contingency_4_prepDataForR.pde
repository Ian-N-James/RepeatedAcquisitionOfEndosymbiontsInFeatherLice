/* Takes the binary strings (Strings.csv). Exports a file (rdat_I.csv) to the sketch 
   directory maintaining the same information presented in a format that facilitates 
   Proximus analysis in R.
*/
void setup(){
  // Load the namelist file.
  Table namelist=loadTable("../sharedResources/namelist.csv","header");
  Table output=new Table();
  int minOfEach=0;  
  output.addColumn("sym");
  output.addColumn("patDis");
  int ti;
  
  // Remove the query organisms that are not included in contingency analyses. This is system specific.
  namelist.removeRow(namelist.findRowIndex("Coang","abb"));
  namelist.removeRow(namelist.findRowIndex("5_11_15_4","abb"));
  ti=namelist.getRowCount();
  for(int i=0;i<ti;i++){
    output.addRow();
    output.setString(i,"sym",namelist.getString(i,"abb"));
    output.setFloat(i,"patDis",namelist.getFloat(i,"patDis"));
  }
  // Load the binary strings.
  Table binstr=loadTable("../sharedResources/Strings.csv","header");
  int tj=ti;String colNam,strStore,strStoreI;
  int cols=0;
  ti=binstr.getRowCount();
  // Reformat the binary strings.
  for(int i=0;i<ti;i++){
    colNam=binstr.getString(i,"tag");
    strStore=binstr.getString(i,"noninverted");
    strStoreI=binstr.getString(i,"inverted");
    if(charNum(strStore,'1')>minOfEach && charNum(strStore,'0')>minOfEach){
      output.addColumn(colNam);
      output.addColumn(colNam+"_I");
      for(int j=0;j<tj;j++){
        output.setInt(j,colNam,int(str(strStore.charAt(j))));
        output.setInt(j,colNam+"_I",int(str(strStoreI.charAt(j))));
      }
    }
  }
  //Save the reformated data file.
  saveTable(output,"rdat_I.csv");
  exit();
}
void draw(){}

// This function returns the number of a given character in a string.
int charNum(String in,char gn){
  int out=0;
  int ti = in.length();
  for(int i=0;i<ti;i++){
    if(in.charAt(i)==gn){
      out++;
    }
  }
  return out;
}
