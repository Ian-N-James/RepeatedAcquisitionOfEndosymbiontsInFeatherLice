void setup(){
  Table zFin;
  Table namelist=loadTable("../sharedResources/namelist.csv","header");
  int ti=namelist.getRowCount();float cutoff;
  for(int i=0;i<ti;i++){
    cutoff=namelist.getFloat(i,"CutOff");
    zFin=loadTable("../sharedResources/zFin/zFin-"+namelist.getString(i,"abb")+"-AnnTabFin.csv","header");
    int tj=zFin.getRowCount();
    for(int j=0;j<tj;j++){
      if(zFin.getFloat(j,"ann")==-0.1){
        zFin.setFloat(j,"ann",-0.1);
        zFin.setFloat(j,"weightScore",20);
      }else{
        if(zFin.getFloat(j,"weightScore")<=cutoff){
          zFin.setFloat(j,"ann",1);
        }else{
          if(zFin.getFloat(j,"weightScore")<1){
            zFin.setFloat(j,"ann",0.5);
          }else{
            if(zFin.getFloat(j,"weightScore")==1){
              zFin.setFloat(j,"ann",0);
            }
          }
        }
      }
    }
    saveTable(zFin,"../sharedResources/zFin/zFin-"+namelist.getString(i,"abb")+"-AnnTabFin.csv");
  }
  exit();
}
void draw(){}
