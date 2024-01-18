/* Takes the modified annotation tables (ZFin-*-AnnTabFin.csv). Exports 
   R code to preform expectation-maximization analysis. Produces 
   LAnner_5_ExpectationMaximization.
*/
void setup(){
  String curDir=sketchPath();
  curDir=join(shorten(curDir.split("/")),"/")+"/sharedResources/";
  PrintWriter Output=createWriter(curDir+"LAnner_5_ExpectationMaximization.R");
  Table zFin;
  Output.println("#Preforms expectation-maximization (E-M) analysis on the snLEDs. Exports a table (EMdata.csv) containing the results of the E-M analysis.");
  Output.println("library(\"mixtools\")");
  String[] namelist=loadStrings("../sharedResources/namelist.txt");
  Boolean pOneAlR; // have we printed one already for this sym?
  //print the symbiont data lines
  int ti=namelist.length;
  for(int i=0;i<ti;i++){
    Output.print(substituter(namelist[i])+"<-c(");
    zFin=loadTable("../sharedResources/zFin/zFin-"+namelist[i]+"-AnnTabFin.csv","header");
    pOneAlR=false;
    int tj=zFin.getRowCount();
    for(int j=0;j<tj;j++){
      if(
      zFin.getFloat(j,"weightScore")>=0 && 
      zFin.getFloat(j,"weightScore")<=1 &&
      zFin.getFloat(j,"ann")>0 &&
      !zFin.getString(j,"flags").equals("RPstar")&&
      !zFin.getString(j,"flags").equals("miscRNA")&&
      !zFin.getString(j,"flags").equals("PpIs1,pseuHS")&&
      !zFin.getString(j,"flags").equals("PpIs2,pseuHS")&&
      !zFin.getString(j,"flags").equals("tRNA")&&
      !zFin.getString(j,"flags").equals("rRNA")&&
      !zFin.getString(j,"flags").equals("Plas,pseuHS")&&
      !zFin.getString(j,"flags").equals("pseuHS")
      ){
        if(pOneAlR){
          Output.print(","+zFin.getFloat(j,"weightScore"));
        }else{
          Output.print(zFin.getFloat(j,"weightScore"));
          pOneAlR=true;
        }
      }
    }
    Output.print(")\n");
  }
  //print the symName Line
  Output.print("symName<-c(");
  for(int i=0;i<ti;i++){
    if(i!=0){
      Output.print(",\""+substituter(namelist[i])+"\"");
    }else{
      Output.print("\""+substituter(namelist[i]+"\""));
    }
  }
  Output.print(")\n\n## setup the Lists of data to be outputted with ridiculous values so that any that do not get a value assigned can be spotted easily\n");
  Output.print("FirstMean<-c(");for(int i=0;i<namelist.length;i++){if(i!=0){Output.print(",-20");}else{Output.print("-20");}}Output.print(")\n");
  Output.print("SecondMean<-c(");for(int i=0;i<namelist.length;i++){if(i!=0){Output.print(",-20");}else{Output.print("-20");}}Output.print(")\n");  
  Output.print("FirstSigma<-c(");for(int i=0;i<namelist.length;i++){if(i!=0){Output.print(",-20");}else{Output.print("-20");}}Output.print(")\n");
  Output.print("SecondSigma<-c(");for(int i=0;i<namelist.length;i++){if(i!=0){Output.print(",-20");}else{Output.print("-20");}}Output.print(")\n");
  Output.print("loglike<-c(");for(int i=0;i<namelist.length;i++){if(i!=0){Output.print(",-20");}else{Output.print("-20");}}Output.print(")\n");
  Output.print("\n## do the EM models\n");
  for(int i=0;i<ti;i++){
    Output.print("EMof"+substituter(namelist[i])+"<-normalmixEM("+substituter(namelist[i])+",mu=c(0,1),sigma=c(1,1))\n");
  }
  Output.print("\n## SetupOutputs\n");
  for(int i=0;i<ti;i++){
    Output.print("FirstMean["+(i+1)+"]<-EMof"+substituter(namelist[i])+"$mu[1]\n");
    Output.print("SecondMean["+(i+1)+"]<-EMof"+substituter(namelist[i])+"$mu[2]\n");
    Output.print("FirstSigma["+(i+1)+"]<-EMof"+substituter(namelist[i])+"$sigma[1]\n");
    Output.print("SecondSigma["+(i+1)+"]<-EMof"+substituter(namelist[i])+"$sigma[2]\n");
    Output.print("loglike["+(i+1)+"]<-EMof"+substituter(namelist[i])+"$loglik[1]\n");
  }
  Output.print("\n## GetOutputs\nsymName\nFirstMean\nSecondMean\nFirstSigma\nSecondSigma\nloglike\n");
  Output.print("out<-data.frame(name=symName,mean1=FirstMean,mean2=SecondMean,sigma1=FirstSigma,sigma2=SecondSigma,loglike=loglike,cutoff=(FirstMean+3*FirstSigma))\n");
  Output.print("write.table(out,\""+curDir+"EMdata.csv\",sep=\",\",row.names=FALSE)");
  Output.flush();Output.close();
  exit();
}
void draw(){
}
String substituter(String in){
  String out="";
  for(int i=0;i<in.length();i++){
    switch(in.substring(i,i+1)){
      case "1":out+="ONE";break;  case "2":out+="TWO";break; case "3":out+="THREE";break; case "4":out+="FOUR";break; case "5":out+="FIVE";break; case "6":out+="SIX";break; case "7":out+="SEVEN";break; case "8":out+="EIGHT";break; case "9":out+="NINE";break; case "0":out+="ZERO";break;
      default:out+=in.substring(i,i+1);
    }}
  return out;
}
