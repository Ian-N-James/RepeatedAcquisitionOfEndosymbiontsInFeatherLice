/* Takes the exported annotation table (file name = *-AnnTab.csv, located 
   in exampleFiles in this repository.). Exports a reformatted version of 
   the annotation table (file name = *-AnnTabFin.csv, located in 
   exampleFiles in this repository) that has non-delimiting commas 
   removed and a simplified header.
*/

void setup(){
  String[] namelist=loadTable("../sharedResources/namelist.csv","header").getStringColumn("abb");
  int ti=namelist.length;
  StringDict rowName=new StringDict();
  rowName.set("Name","name");
  rowName.set("locus_tag","tag");
  rowName.set("Direction","dir");
  rowName.set("Min (with gaps)","min");
  rowName.set("Max (with gaps)","max");
  rowName.set("product","product");
  for(int i=0;i<ti;i++){
    String[] anntab=loadStrings("../sharedResources/AnnTab/"+namelist[i]+"-AnnTab.csv");
    int tj=anntab.length;
    PrintWriter anntabfin=createWriter("../sharedResources/AnnTabFin/"+namelist[i]+"-AnnTab.csv");
    String[] line=lineSpliter(anntab[0]);
    int tk=line.length;
    int prodPos=-1,tagPos=-1;
    for(int k=0;k<tk;k++){
      line[k]=rowName.get(line[k]);
      if(line[k].equals("tag")){tagPos=k;}
      if(line[k].equals("product")){prodPos=k;}
    }
    anntabfin.print(join(line,","));
    line=lineSpliter(anntab[1]);
    if(tagPos>=0){line[tagPos]="Sant_0000";}
    if(prodPos>=0){line[prodPos]="Origin of Replication";}
    anntabfin.print("\n"+join(line,","));
    for(int j=2;j<tj;j++){
      line=lineSpliter(anntab[j]);
      anntabfin.print("\n"+join(line,","));
    }
    anntabfin.flush();
    anntabfin.close();
  }
  exit();
}
void draw(){}
String[] lineSpliter(String in){
  String[] out={};
  boolean inQuote=false;
  String tem="";
  int ti=in.length();
  for(int i=0;i<ti;i++){
    if(inQuote){
      if(in.charAt(i)=='\"'){inQuote=false;}
      else{
        if(in.charAt(i)==','){
        }else{
          tem+=in.charAt(i);
        }
      }
    }else{
      if(in.charAt(i)=='\"'){inQuote=true;}
      else{
        if(in.charAt(i)==','){
          out=append(out,tem);
          tem="";
        }else{tem+=in.charAt(i);}
      }
    }
  }
  out=append(out,tem);
  return out;
}
