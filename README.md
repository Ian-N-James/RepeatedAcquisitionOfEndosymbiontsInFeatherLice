# RepeatedAcquisitionOfEndosymbiontsInFeatherLice

Overview:
These scripts are used for analysis of gene retention and contingency in genome degeneration.

Requirements:
• Processing IDE
• R
R packages:
• mixtools
• cba
• scales

Starting Files:
Files present for each symbiont (within subfolders):
    • SymCon - Symbiont consensus sequence from alignment.
        Names = [symbiont abbreviation]-SymCon.txt
    • PraCon - S. Praecaptivus sequence from alignment. The need to use one from the alignment is because insertions  will throw off the locations of genes in an ungapped version of the sequence.
        Names = [symbiont abbreviation]-PraCon.txt
    • AnnTab - Table containing the locations of each gene in a given symbiont.
        Names = [symbiont abbreviation]-AnnTab.csv
Files with one copy in the main sharedResources folder:
    • namelist.csv  - A CSV file containing the patristic distances to S. praecaptivus strain HS1 (column "patDis"), abbreviations (column "abb"), cutoffs for calling a gene as intact* (column "CutOff"), and names (column "name") of each symbiont. Symbionts should be ordered from lowest to highest patristic distance.
    • namelist.txt  - A list of the abbreviations for each symbiont. Symbionts should be ordered from lowest to highest patristic distance.
    • namelist2.txt - A list of the abbreviations for each symbiont, omitting those that should not be used in the contingency analyses. Symbionts should be ordered from lowest to highest patristic distance. 
* Cutoffs should be set to zero before the script R_Read is run.

Files generated/modified during the run:
Files present for each symbiont (within subfolders):
    • AnnTabFin - AnnTabs that have under gone a minor reformatting to be compatible with the pipeline.
        Generated by AnnTabReformat.pde
        Names = [symbiont abbreviation]-AnnTabFin.csv
    • ComCon - Combines the SymCon and PraCon files.
        Generated by 
        Names = [symbiont abbreviation]-ComCon.txt
    • zFin - AnnTabFins with additional information including LEDs (column "score"), snLEDs (column "weightScore"), and annotation values.
        Generated by: LAnner
        Modified by:  LAnner2
        Names = zFin-[symbiont abbreviation]-AnnTabFin.csv
  Files with one copy in the main sharedResources folder:
    • RCODE.R - Code for preforming the E-M analysis. Includes data to be used in such. 
        Generated by: R_Write
    • EMdata.csv - Outputted data from the E-M analysis.
        Generated by: RCODE.R
    • namelist.csv - Cutoffs are added after the E-M analysis is preformed.
        Generated by: N/A this is one of the intital files.
        Modified by: R_Read

There are two main pipelines included. These are the initial annotation pipeline, and the contingency pipeline.

Initial annotation pipeline:
Run scripts in the following order:
• ComConGenerator.pde
• AnnTabReformat.pde
• L_Anner.pde
• R_Write.pde
• sharedResources/RCODE2.R
    Note(s):
    • Preforms E-M analysis.
    • Should be opened and copied to the R console.
• R_Read.pde
• L_Anner2.pde
    Note(s):
    • After L_Anner2.pde has run, the contingency related scripts can be run.
• Three_D_Hist
• G_Render.pde
• G_Render/G_Veiwer.pde

Contingency Analysis Scripts:
Run these Scripts first:
• StringMaker.pde
• ShannonAdder.pde
• HammString.pde
Then run the next sets in whichever order:
Proximus (run these in order):
• prepDatForR.pde
• proximusAnalysis.R
• namAndFilt.pde
C_Map (these can run in any order, so long as "C_Map_Combined.pde" is last):
• C_Map.pde
• C_Map_panel1.pde
• C_Map_panel2.pde
• C_Map_Combined.pde
