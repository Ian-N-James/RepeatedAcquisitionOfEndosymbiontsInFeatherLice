String[] namelist=loadStrings("../sharedResources/namelist2.txt");
Table AnnTab;
AnnTab=loadTable("../sharedResources/zFin/zFin-"+namelist[0]+"-AnnTabFin.csv","header");
int tj=AnnTab.getRowCount();
String[] Tag=new String[tj];
String[] noninverted=new String[tj];
String[] inverted=new String[tj];
boolean[] ignore=new boolean[tj];;
for(int j=0;j<tj;j++){
  Tag[j]=AnnTab.getString(j,"tag");
  if(!AnnTab.getString(j,"flags").equals("")||AnnTab.getString(j,"flags").equals("PrTrFr")){
    ignore[j]=true;
  }else{
    ignore[j]=false;
      if(AnnTab.getFloat(j,"ann")==1){
        noninverted[j]="1";inverted[j]="0";
      }else{
        noninverted[j]="0";inverted[j]="1";
      }
    }
}
int ti=namelist.length;
for(int i=1;i<ti;i++){
  AnnTab=loadTable("../sharedResources/zFin/zFin-"+namelist[i]+"-AnnTabFin.csv","header");
  for(int j=0;j<tj;j++){
    if(!ignore[j]){
      if(AnnTab.getFloat(j,"ann")==1){
        noninverted[j]+="1";
        inverted[j]+="0";
      }else{
        noninverted[j]+="0";
        inverted[j]+="1";
      }
    }
  }
}
String Con1="";
for(int i=0;i<ti;i++){
  Con1+="1";
}
PrintWriter Output=createWriter("../sharedResources/Strings.csv");
Output.println("tag,noninverted,inverted");
for(int j=0;j<tj;j++){
  if(!ignore[j]&&!noninverted[j].equals(Con1)&&!inverted[j].equals(Con1)){
    Output.println(Tag[j]+","+noninverted[j]+","+inverted[j]);
  }
}
Output.flush();
Output.close();
exit();
