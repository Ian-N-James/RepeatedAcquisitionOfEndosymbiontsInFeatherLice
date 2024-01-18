# RepeatedAcquisitionOfEndosymbiontsInFeatherLice

Overview:<br>
These scripts were used in a study to annotate the gene inventories of feather louse endosymbionts following alignment of their sequence reads to a whole genome sequence obtained from a close free living relative, _S. praecaptivus_. Comparative analyses then facilitated the identification of intact genes within the symbiont genomes using an Expectation Maximization (E-M) algorithm. Contingencies underpinning gene losses among the symbionts were then detected and validated using the Proximus algorithm. Additional scripts facilitate visualization of resulting data.<br>
<br>
Dependencies:<br>
• Processing<br>
• R with packages mixtools, cba and scales<br>
<br>
Usage:<br>
<br>
Identification of intact genes:<br>
Sequence reads from query organisms (feather louse endosymbionts in our study) were mapped to the whole genome sequence of a reference organism (S. praecaptivus in our study) and visualized in Geneious v8.1.8. The consensus sequence, reference sequence and annotation table were exported and then processed by the LAnner pipeline, comprising seven scripts. All scripts require a namelist.csv file (located in exampleFiles in this repository) to facilitate batch processing. This namelist.csv file and any associated input files should be located in a sharedResources directory alongside directories containing the individual processing scripts. Note that filenames are typically appended with identifiers (*-). For example Cocln-SymCon.txt represents data from the endosymbiont of Columbicola claytoni (see namelist.csv for complete list).<br>
• LAnner_1_ComConGenerator - Takes the consensus sequence of a given organism (file name = *-SymCon.txt, located in exampleFiles in this repository) aligned to the reference sequence (file name = *-PraCon.txt located in exampleFiles in this repository) that is gapped in accordance with the alignment. The script then outputs a sequence file encoded to display the contents of both input files (file name = *-ComCon.txt, located in exampleFiles in this repository).<br>
• LAnner_2_AnnTabReformatter - Takes the exported annotation table (file name = *-AnnTab.csv, located in exampleFiles in this repository.). Exports a reformatted version of the annotation table (file name = *-AnnTabFin.csv, located in exampleFiles in this repository) that has non-delimiting commas removed and a simplified header.<br>
• LAnner_3_snLED - Takes the combined sequence file (*-ComCon.txt), and the reformatted annotation table (*-AnnTabFin.csv). Calculates all open reading frames (ORFs) with either canonical or alternative bacterial start codons within ± 25 bp of the ortholog in the reference genome. The script then determines the ORF with the lowest Levenshtein edit distance (LED) from the reference ortholog for each gene. Exports a modified annotation table (ZFin-*-AnnTabFin.csv) that additionally includes the LED, a size normalized LED (snLED), and both the predicted nucleotide and protein sequences of the ORF with the lowest LED, for each gene.<br>
• LAnner_4_R_Write - Takes the modified annotation tables (ZFin-*-AnnTabFin.csv). Exports R code to preform expectation-maximization analysis. Produces LAnner_5_ExpectationMaximization. <br>
• LAnner_5_ExpectationMaximization - Preforms expectation-maximization (E-M) analysis on the snLEDs. Exports a table (EMdata.csv) containing the results of the E-M analysis.<br>
• LAnner_6_R_Read - Takes the EMdata.csv table and appends the entries in the namelist table with snLED values estimated by E-M analysis to represent cutoffs for intact genes.<br>
• LAnner_7_Predictor - Takes the modified annotation tables (ZFin-*-AnnTabFin.csv) and appends all entries with predictions of gene status (0 = missing gene, 0.5 = pseudogene, 1 = intact gene) based on the cutoffs represented in namelist.csv.<br>
<br>
Contingency Analysis:<br>
Data derived from LAnner_7_Predictor (ZFin-*-AnnTabFin.csv), maintaining the E-M-derived predictions for all genes is subject to comparative analysis to identify patterns of gene inactivation or loss that are contingent in nature. Two approaches were employed for this objective: The first focuses on analysis of Hamming distances derived from binary strings encoding predictions of functional status among individual genes. The second utilizes the Proximus algorithm in the R package ‘cba’, which identifies patterns of contingency. Again, some of these scripts require namelist.csv to be present in a sharedResources directory alongside directories containing the individual processing scripts.<br>
• Contingency_B_1_StringMaker - Takes annotation tables (ZFin-*-AnnTabFin.csv) that have been processed by LAnner_7_Predictor. Exports a file (Strings.csv) that contains binary strings for each fractionally retained gene. The output of this script is utilized by both approaches. <br>
• Contingency_1_1_ShannonAdder - Takes the binary stings (Strings.csv). Exports a file (Shannon.csv) that is equivalent to Strings.csv file appended to contain Shannon entropies of the the binary strings of each fractionally retained gene. This is the first script used in the first approach.<br>
• Contingency_1_2_HammingString - Takes the binary strings (Strings.csv). Exports a file (HammingDistance.csv) that contains the Hamming distances between the binary strings of each fractionally retained gene. This is the second script used in the first approach.<br>
• Contingency_2_1_prepDataForR - Takes the binary strings (Strings.csv). Exports a file (rdat_I.csv, located within the sketch directory) that contains the binary strings transposed to a format more easily imported into R.<br>
• Contingency_2_2_proximusAnalysis - Takes the binary strings (rdat_I.csv). Preforms analysis using the Proximus algorithim. Exports a file (dataFromR2.txt) that contains the results of the analysis. <br>
• Contingency_2_3_nameAndFilter - Takes the results of the analysis using Proximus (dataFromR2.txt), and a file (tagsAndNams.csv) that lists the names of genes by tag. Exports three files. The first file (removal Flag.csv) is the entire results with separated sets that are listed by tag exclusively, and has a column “remove” indicating sets that are subsets of other sets and are removed in the other two files. The second file (proximusResults_tagsOnly.csv) contains the same information without sets that are subsets of other sets or a removal flag column. The third file (proximusResults.csv) lists the same information as the second, except that genes that have names are listed by their names.<br>
<br>
Visualization:<br>
• G_Render - Takes the annotation tables (ZFin-*-AnnTabFin.csv) that have been processed by LAnner_7_Predictor and a color map (Gradient.txt). Renders several outputs that show gene status (with intact genes in green, pseudo genes in purple and missing genes in black) across symbiont genomes. One output, saved within the directory of G_Veiwer, has the Tab containing gene names saved as a separate file (tab.tif) from the gene statuses (geneRender.tif). The other outputs are saved within the sketch directory. Two use sub-pixel rendering (GRender.tif and GRender_2.tif) and differ by the amount of whitespace to the right and top of the image. The final output (GRender_1px.tif) does not use sub-pixel rendering, instead rendering at one pixel per gene.<br>
• G_Viewer - Takes the outputs of G_Render saved within its sketch directory (tab.tif and geneRender.tif) and a table of gene annotations (Tags.csv). Acts as an interactive viewer for a larger G_Render output (rendered at 4 pixels per gene). The sketch directory for G_Viewer is located within the sketch directory of G_Render.<br>
• C_Map - Takes the Shannon entropies (Shannon.csv), Hamming distances (HammingDistance.csv), a file containing instructions on which genes to label (labels.txt), and two color maps (linear_ternary-blue_0-44_c57_n256.csv and linear_ternary-red_0-50_c52_n256.csv). Produces the Circle map (out_C_BinaryColMap_NT.tif) displaying the contingent relationships between all fractionally retained genes. <br>
The color maps used for rendering the relationships are from CET Perceptually Uniform Colour Maps (https://colorcet.com/index.html) by Peter Kovesi (See also Good Colour Maps: How to Design Them 2015) and are under a Creative Commons BY License (https://creativecommons.org/licenses/by/4.0/).<br>
• C_Map_panel1 - Takes the Shannon entropies (Shannon.csv), binary strings (included in the same file as the Shannon entropies) and Hamming distances (HammingDistance.csv).  Produces a file (PanelA.tif) containing a matrix of functional statuses select genes encoding respiratory chain components in each louse endosymbiont, the Shannon entropies corresponding to the binary strings of these genes, the Hamming distances between the binary strings of these genes and the relationship strengths.<br>
• C_Map_panel2 - Takes the Shannon entropies (Shannon.csv). Produces a plot of the frequencies of strings with different Shannon entropies.<br>
• Three_D_Hist - Takes annotation tables (ZFin-*-AnnTabFin.csv) that have been processed by LAnner_7_Predictor, and a color map (Gradient.txt). Produces a three dimensional histogram of the number of genes with given snLEDs in accordance with patristic distance from S. praecaptivus, with the symbiont having the lowest patristic distance in the back. Predicted intact genes are plotted in green, and predicted pseudogenes are plotted in purple. <br>
