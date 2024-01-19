/* Takes the results of the Proximus analysis (dataFromR2.txt), and a file 
   (tagsAndNams.csv) that lists the locus tags and names of genes and 
   exports three output files. The first (removal Flag.csv) is the entire 
   results with separated sets that are listed by tag exclusively, and has 
   a column "remove" indicating sets that are subsets of other sets and 
   are eliminated in the other two files. The second file 
   (proximusResults_tagsOnly.csv) contains the same information without 
   sets that are subsets of other sets or a removal flag column. The third
   file (proximusResults.csv) lists the same information as the second 
   with gene names utilized instead of locus tags with the caveat that 
   some genes have not been assigned proper names and can therefore only 
   be identified by locus tags.
*/
void setup(){
  // Load the results from Contingency_5_proximusAnalysis.
  Table in=loadTable("dataFromR2.txt","tsv,header");
  in.removeColumn("");
  in.removeColumn("id");
  // Load the file containing locus tags and gene names.
  Table tagNam=loadTable("tagsAndNams.csv","header");
  Table out=new Table();
  String[] inColNam=in.getColumnTitles();
  int ti=inColNam.length-1;
  for(int i=0;i<ti;i++){
    out.addColumn(inColNam[i]);
  }
  /* Split the genes in each pattern into two sets, such that those that each set 
     has genes with direct contingent relationships to each other, and reciprocal 
     relationships to genes in the opposite set.    */
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
      pats[i].append(int(hold[j].split("_")[1]));
      if(hold[j].substring(hold[j].length()-2).equals("_I")){
        rec+=","+hold[j].substring(0,hold[j].length()-2);
      }else{
        nonrec+=","+hold[j];
      }
      if(nonrec.length()>1){out.setString(i,"set 1",nonrec.substring(1));
      }else{out.setString(i,"set 1","");}
      if(rec.length()>1){out.setString(i,"set 2",rec.substring(1));
      }else{out.setString(i,"set 2","");
      }
    }
    pats[i].sort();
  }
  // Filter duplicates and subsets, initially flagging them.
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
          IntList MSMs=new IntList();
          
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
  //Save the first output file.
  saveTable(out,"removal Flag.csv");
  //Remove duplicates and subsets.
  for(int i=ti-1;i>=0;i--){
    if(out.getInt(i,"remove")==1){out.removeRow(i);}
  }
  out.removeColumn("remove");
  //Save the second output file.
  saveTable(out,"proximusResults_tagsOnly.csv");
  // Where possible, replace locus tags with gene names.
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
  // Save the third output file.
  saveTable(out,"proximusResults.csv");
  exit();
}
void draw(){}
// Check if a set of genes is a subset or duplicate of another.
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
