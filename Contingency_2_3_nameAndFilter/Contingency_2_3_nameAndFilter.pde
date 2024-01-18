/* Takes the results of the analysis using Proximus (dataFromR2.txt), and a 
   file (tagsAndNams.csv) that lists the names of genes by tag. Exports 
   three files. The first file (removal Flag.csv) is the entire results 
   with separated sets that are listed by tag exclusively, and has a column 
   "remove" indicating sets that are subsets of other sets and are removed 
   in the other two files. The second file (proximusResults_tagsOnly.csv) 
   contains the same information without sets that are subsets of other 
   sets or a removal flag column. The third file (proximusResults.csv) 
   lists the same information as the second, except that genes that have 
   names are listed by their names.
*/
void setup(){
  Table in=loadTable("dataFromR2.txt","tsv,header");
  //remove a blank column.
  in.removeColumn("");
  in.removeColumn("id");
  Table tagNam=loadTable("tagsAndNams.csv","header");
  Table out=new Table();
  String[] inColNam=in.getColumnTitles();
  int ti=inColNam.length-1;
  for(int i=0;i<ti;i++){
    out.addColumn(inColNam[i]);
  }
  out.addColumn("set 1");
  out.addColumn("set 2");
  ti=in.getRowCount();String[] hold;String rec,nonrec;
  IntList[] pats=new IntList[ti];
  for(int i=0;i<ti;i++){
    hold=in.getString(i,"members").split(" ");
    pats[i]=new IntList();
    out.addRow();
    int tj=hold.length;
    rec="";
    nonrec="";
    out.setInt(i,"size",in.getInt(i,"size"));
    out.setFloat(i,"relative error",in.getFloat(i,"relative error"));
    out.setFloat(i,"Frobenius norm",in.getFloat(i,"Frobenius norm"));
    out.setFloat(i,"Jaccard similarity",in.getFloat(i,"Jaccard similarity"));
    
    for(int j=0;j<tj;j++){
      //println(hold[j].split("_"));
      pats[i].append(int(hold[j].split("_")[1]));
      if(hold[j].substring(hold[j].length()-2).equals("_I")){
        rec+=","+hold[j].substring(0,hold[j].length()-2);
      }else{
        nonrec+=","+hold[j];
      }
      //.subsrting(1) removes the comma that comes before the first gene tag.
      if(nonrec.length()>1){out.setString(i,"set 1",nonrec.substring(1));
      }else{out.setString(i,"set 1","");}
      if(rec.length()>1){out.setString(i,"set 2",rec.substring(1));
      }else{out.setString(i,"set 2","");
      }
    }
    pats[i].sort();
  }
  // Filter duplicates (resulting from inclusion of reciprocal sequences)
  IntList pattern,members;
  boolean[] inPat=new boolean[ti];
  out.addColumn("remove",Table.INT);
  for(int i=0;i<ti;i++){
    if(!inPat[i]){
      StringList A,B;
      A=new StringList(concat(out.getString(i,"set 1").split(","),out.getString(i,"set 2").split(",")));
      pattern=new IntList(i);
      members=new IntList(A.size());
      for(int j=0;j<ti;j++){
        B=new StringList(concat(out.getString(j,"set 1").split(","),out.getString(j,"set 2").split(",")));
        if(isSubset(A,B)){
          pattern.append(j);
          members.append(B.size());
        }
      }
      if(pattern.size()>1){
        int tj=pattern.size();
        int max=members.max();
        IntList maximas=new IntList(); 
        for(int j=0;j<tj;j++){
          if(members.get(j)!=max){
            out.setInt(pattern.get(j),"remove",1);
          }else{
            maximas.append(j);
          }
        }
        if(maximas.size()>1){
          tj=maximas.size();
          IntList MSMs=new IntList();// set 1 Members.
          
          for(int j=0;j<tj;j++){
            MSMs.append(out.getString(pattern.get(maximas.get(j)),"set 1").split(",").length);
          }
          int keep=pattern.get(maximas.get(MSMs.maxIndex()));
          for(int j=1;j<tj;j++){
            if(j!=keep){
              out.setInt(pattern.get(maximas.get(j)),"remove",1);
            }
          }
        }
      }
    }
  }
  saveTable(out,"removal Flag.csv");
  for(int i=ti-1;i>=0;i--){
    if(out.getInt(i,"remove")==1){out.removeRow(i);}
  }
  out.removeColumn("remove");
  saveTable(out,"proximusResults_tagsOnly.csv");
  ti=out.getRowCount();
  for(int i=0;i<ti;i++){
    String[] patHold=out.getString(i,"set 1").split(",");
    int tj=patHold.length;
    for(int j=0;j<tj;j++){
      int namIndex=tagNam.findRowIndex(patHold[j],"tag");
      if(namIndex>=0){patHold[j]=tagNam.getString(namIndex,"name");}
    }
    out.setString(i,"set 1",join(patHold,", "));
    
    patHold=out.getString(i,"set 2").split(",");
    tj=patHold.length;
    for(int j=0;j<tj;j++){
      int namIndex=tagNam.findRowIndex(patHold[j],"tag");
      if(namIndex>=0){patHold[j]=tagNam.getString(namIndex,"name");}
    }
    out.setString(i,"set 2",join(patHold,", "));
    
  }
  saveTable(out,"proximusResults.csv");
  
  
  exit();
}
void draw(){}
boolean isSubset(StringList in1,StringList in2){
  if(in1.size()>in2.size()){
    int ti=in2.size(),tj=in1.size();
    boolean[] check=new boolean[ti];
    for(int i=0;i<ti;i++){
      for(int j=0;j<tj;j++){
        if(in2.get(i).equals(in1.get(j))){
          check[i]=true;
        }
      }
      if(!check[i]){return false;}
    }
    
  }else{
    int ti=in1.size(),tj=in2.size();
    boolean[] check=new boolean[ti];
    for(int i=0;i<ti;i++){
      for(int j=0;j<tj;j++){
        if(in1.get(i).equals(in2.get(j))){
          check[i]=true;
        }
      }
      if(!check[i]){return false;}
    }
  }
  return true;
}
