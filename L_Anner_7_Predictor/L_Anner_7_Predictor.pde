/* Takes the modified annotation tables (ZFin-*-AnnTabFin.csv) and appends all 
   entries with predictions of gene status (0 = missing gene, 0.5 = 
   pseudogene, 1 = intact gene) based on the cutoffs represented in 
   namelist.csv.
*/
void setup(){
  Table zFin;
  // Load the namelist file.
  Table namelist=loadTable("../sharedResources/namelist.csv","header");
  int ti=namelist.getRowCount();float cutoff;
  for(int i=0;i<ti;i++){
    cutoff=namelist.getFloat(i,"CutOff");
    // Load the zFin-*-AnnTabFin files.
    zFin=loadTable("../sharedResources/zFin/zFin-"+namelist.getString(i,"abb")+"-AnnTabFin.csv","header");
    int tj=zFin.getRowCount();
    // Assign the predicted status. A value of -0.1 flags non protein-coding genes that are ignored in subsequent analyses.
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
    // Save the output files.
    saveTable(zFin,"../sharedResources/zFin/zFin-"+namelist.getString(i,"abb")+"-AnnTabFin.csv");
  }
  exit();
}
void draw(){}
