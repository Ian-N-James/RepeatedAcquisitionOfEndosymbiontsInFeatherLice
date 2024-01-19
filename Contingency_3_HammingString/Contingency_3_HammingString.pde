/* Takes the binary strings (Strings.csv) and computes Hamming distances between all
   pairwise combinations of strings. Exports a file (HammingDistance.csv) that 
   tabulates the Hamming distance values between the binary strings.
*/
// Load the binary strings.
String[] strings=loadTable("../sharedResources/Strings.csv","header").getStringColumn("noninverted");
int ti=strings.length;
int[][] HammingDistances=new int[ti][ti];
int tk=strings[0].length();
// Calculate the Hamming distances.
for(int i=0;i<ti;i++){
  for(int j=0;j<ti;j++){
    for(int k=0;k<tk;k++){
      if(strings[i].charAt(k)!=strings[j].charAt(k)){HammingDistances[i][j]++;}
    }
  }
}
// Load locus tags (from the same file as the binary strings).
String[] tags=loadTable("../sharedResources/Strings.csv","header").getStringColumn("tag");
Table output=new Table();
output.addColumn("tag");
output.addColumn("data");
for(int i=0;i<ti;i++){output.addColumn(tags[i]);}
for(int i=0;i<ti;i++){
  output.addRow();
  output.setString(i,"tag",tags[i]);
  output.setString(i,"data",strings[i]);
  for(int j=0;j<ti;j++){
    output.setInt(i,j+2,HammingDistances[i][j]);
  }
}
//Save the table of Hamming distances.
saveTable(output,"../sharedResources/HammingDistance.csv");
exit();
