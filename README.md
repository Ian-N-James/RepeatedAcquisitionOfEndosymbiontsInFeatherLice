# RepeatedAcquisitionOfEndosymbiontsInFeatherLice

Overview:<br>
These scripts are used for analysis of gene retention and contingency in genome degeneration.<br>
<br>
Requirements:<br>
• Processing IDE<br>
• R<br>
R packages:<br>
• mixtools<br>
• cba<br>
• scales<br>

Starting Files:<br>
Files present for each symbiont (within subfolders):<br>
• SymCon - Symbiont consensus sequence from alignment.<br>
&emsp;Names = [symbiont abbreviation]-SymCon.txt<br>
• PraCon - S. Praecaptivus sequence from alignment. The need to use one from the alignment is because insertions  will throw off the locations of genes in an ungapped version of the sequence.<br>
&emsp;Names = [symbiont abbreviation]-PraCon.txt<br>
• AnnTab - Table containing the locations of each gene in a given symbiont.<br>
&emsp;Names = [symbiont abbreviation]-AnnTab.csv<br>
Files with one copy in the main sharedResources folder:<br>
&emsp;• namelist.csv  - A CSV file containing the patristic distances to S. praecaptivus strain HS1 (column "patDis"), abbreviations (column "abb"), cutoffs for calling a gene as intact* (column "CutOff"), and names (column "name") of each symbiont. Cutoffs are set to zero until R_Read.pde is run. Symbionts should be ordered from lowest to highest patristic distance.<br>
• namelist.txt  - A list of the abbreviations for each symbiont. Symbionts should be ordered from lowest to highest patristic distance.<br>
• namelist2.txt - A list of the abbreviations for each symbiont, omitting those that should not be used in the contingency analyses. Symbionts should be ordered from lowest to highest patristic distance. <br>


Files generated/modified during the run:<br>
Files present for each symbiont (within subfolders):<br>
• AnnTabFin - AnnTabs that have under gone a minor reformatting to be compatible with the pipeline.<br>
&emsp;Generated by AnnTabReformat.pde<br>
&emsp;Names = [symbiont abbreviation]-AnnTabFin.csv<br>
• ComCon - Combines the SymCon and PraCon files.<br>
&emsp;Generated by ComConGenerator.pde<br>
&emsp;Names = [symbiont abbreviation]-ComCon.txt<br>
• zFin - AnnTabFins with additional information including LEDs (column "score"), snLEDs (column "weightScore"), and annotation values.<br>
&emsp;Generated by: LAnner<br>
&emsp;Modified by:  LAnner2<br>
&emsp;Names = zFin-[symbiont abbreviation]-AnnTabFin.csv<br>
  Files with one copy in the main sharedResources folder:<br>
• RCODE.R - Code for preforming the E-M analysis. Includes data to be used in such. <br>
&emsp;Generated by: R_Write<br>
• EMdata.csv - Outputted data from the E-M analysis.<br>
&emsp;Generated by: RCODE.R<br>
• namelist.csv - Cutoffs are added after the E-M analysis is preformed.<br>
&emsp;Generated by: N/A this is one of the intital files.<br>
&emsp;Modified by: R_Read<br>
<br>
There are two main pipelines included. These are the initial annotation pipeline, and the contingency pipeline.<br>
<br>
Initial annotation pipeline:<br>
Run scripts in the following order:<br>
• ComConGenerator.pde<br>
• AnnTabReformat.pde<br>
• L_Anner.pde<br>
• R_Write.pde<br>
• sharedResources/RCODE2.R<br>
&emsp; Note(s):<br>
&emsp;&emsp;• Is generated by R_Write.pde
&emsp;&emsp;• Preforms E-M analysis.<br>
&emsp;&emsp;• Should be opened and copied to the R console.<br>
• R_Read.pde<br>
• L_Anner2.pde<br>
&emsp;Note(s):<br>
&emsp;&emsp;• After L_Anner2.pde has run, the contingency related scripts can be run.<br>
• Three_D_Hist<br>
• G_Render.pde<br>
• G_Render/G_Veiwer.pde<br>

Contingency Analysis Scripts:<br>
Run these Scripts first:<br>
• StringMaker.pde<br>
• ShannonAdder.pde<br>
• HammString.pde<br>
Then run the next sets in whichever order:<br>
Proximus (run these in order):<br>
• prepDatForR.pde<br>
• proximusAnalysis.R<br>
&emsp; Note(s):<br>
&emsp;&emsp; • Should be sourced into R.
• namAndFilt.pde<br>
C_Map (these can run in any order, so long as "C_Map_Combined.pde" is last):<br>
• C_Map.pde<br>
• C_Map_panel1.pde<br>
• C_Map_panel2.pde<br>
• C_Map_Combined.pde<br>
