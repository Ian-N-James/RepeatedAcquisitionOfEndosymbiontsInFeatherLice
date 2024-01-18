/* Takes the binary strings (Strings.csv). Exports a file (rdat_I.csv) that 
   contains the binary strings transposed to a format more easily imported 
   into R.
*/
void setup(){
  Table namelist=loadTable("../sharedResources/namelist.csv","header");
  Table output=new Table();
  int minOfEach=0;// the actual minimum is one greater than this value.  
  output.addColumn("sym");
  output.addColumn("patDis");
  int ti;
  
  //remove the symbionts that are not included in contigency analyses.
  namelist.removeRow(namelist.findRowIndex("Coang","abb"));
  namelist.removeRow(namelist.findRowIndex("5_11_15_4","abb"));
  
  
  ti=namelist.getRowCount();
  for(int i=0;i<ti;i++){
    output.addRow();
    
    output.setString(i,"sym",namelist.getString(i,"abb"));
    output.setFloat(i,"patDis",namelist.getFloat(i,"patDis"));
  }
  Table binstr=loadTable("../sharedResources/Strings.csv","header");
  int tj=ti;String colNam,strStore,strStoreI;
  int cols=0;
  ti=binstr.getRowCount();
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
  saveTable(output,"rdat_I.csv");
  exit();
}
void draw(){}
Table oldOrder(Table in){
  int[] order={0,1,2,3,5,7,8,6,9,4,11,14,12,13,15,22,10,21,19,16,17,18,26,24,23,25,28,27,32,30,27,31,29,34,33,35};
  Table out=new Table();
  out.addColumn("patDis");
  out.addColumn("abb");
  for(int i=0;i<36;i++){
    out.addRow();
    out.setString(i,"abb",in.getString(order[i],"abb"));
    out.setFloat(i,"patDis",in.getFloat(order[i],"patDis"));
  }
  
  return out;
}
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
