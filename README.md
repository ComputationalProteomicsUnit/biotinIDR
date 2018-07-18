# biotinIDR : Contains code and input files for Biotin IDR study 

### Main code  : R and Python code for analysing data
Note: Before you run this code, you will need to download the following scripts
1. R installed on your computer  
2. Python installed on your computer 
3. GO.R from https://github.com/TomSmithCGAT/CamProt_R/blob/master/GO.R (if running from scratch you'll need to modify the code to be able to read in the column "to.id" instead of "GO.ID". Alternatively, you can change the coloumn names in your input dataframe to be "GO.ID"
4. 'd2p2.py' and 'protinfo.py' from https://github.com/TomSmithCGAT/CamProt/tree/master/camprot/proteomics (needs to be complied)

 File | Description | Language/Format
--------------------|--------------------------- |-----  
biopepFunctions.R	|	Functions required for the R code to run	|	R  
biopep.pub.Rmd 	|	R markdown file that contains the main code used in the data analysis	|	R  
nods.pub.html	|	HTML rendering of R markdown file biopep.pub.Rmd	|	Html
Get_IDRs-DM-v2.py	|	Adds IDR info for all proteins in each of the biotin studies	|	Python  
runpy.sh	|	Calls Get_IDRs-DM-v2.py using various inputs and parameters	|	Shell
printUrl.py	|	Uses Protter URLs generated by Get_IDRs-DM-v2.py and retrieved PDF files into specific folders	|	Python  
runUrl.sh	|	Calls printUrl.py on a file of interest	|	Shell  
GO.R	|	Tom's code to add ancestor terms to existing GO terms	|	R    
d2p2.py	|	Tom's code to use the d2p2 API to use 9 different callers to identify intrinsically disordered regions	|	Python     
protinfo.py	|	 Tom's code to query uniprot/swissprot database on Ensembl	|	Python  


### Input : Folder containing all input files needed for code to run

File| Contents   
----------|------------- 
Ab-APEX-with-biotins-ptms-by-gene.rds	|	Biotinylation sites from Ab-APEX study summarised by gene with IDR and PTM information added in  
BioSITe-with-biotins-ptms-by-gene.rds	|	Biotinylation sites from BioSITE study summarised by gene with IDR and PTM information added in  
DiDBIT-with-biotins-ptms-by-gene.rds	|	Biotinylation sites from DiDBIT study summarised by gene with IDR and PTM information added in  
SpotBioID-with-biotins-ptms-by-gene.rds	|	Biotinylation sites from Ab-APEX study summarised by gene with IDR and PTM information added in  
Biotin-PTM-IDR-data-for-all-studies-and-callers.rds	|	Biotin,IDR and PTM information for all studies and caller combinations  
HEK293-annotated-list-of-proteins.rds	|	HEK293 list of proteins from Geiger et al  
HEK293-GO-annotation-with-ancestors.rds	|	GO annotations with ancestor terms for the HEK293 background protein list from Geiger et al  
HEK293-with-IDR-VSL2b-and-PTMs.rds	|	HEK293 peoteins from Geiger et al with PTMs and IDRs based on the algorithm VSL2b
PTMs-by-gene-from-Phosphositeplus.rds	|	PTM data from PhosphositePlus summarised by gene and for phosphorylation, ubiquitination, sumoylation, acetylation  
go-pro-kegg.rda	|	An object that contains all GO, Interpro and KEGG annotations for proteins in each of the 4 studies  


#### RDS files can be read using the command readRDS as shown below
```
readRDS("Input/Biotin-PTM-IDR-data-for-all-studies-and-callers.rds")
```

#### RDA files can be loaded using the command "load" as shown below
```
load("Input/go-pro-kegg.rda")
```

### Output : A description of all files generated by running nods_final.Rmd.  

Date_Output : Date-tagged output folder This folder should contain the following files
Date = Date of the analysis in format yyyy-mm-dd ; Time = Time of analysis in format hh:mm:ss  
Most of these files are related to the Figures in the main manuscript with some additional information

File| Contents   
----------|-------------  
Figure3A_Number.biotins.vs.FractionIDR-VSL2b.pdf	|	Violin plots comparing fraction of IDR and biotin number in each protein  
Figure3B_Distribution-of-Biotins-in.out.of.IDRs-by-caller-and-study.pdf	|	Bar plot showing expected and observed rates of biotin for all studies and callers  
Figure3B_Binomial-test-stats.txt	|	Bionomial test values and parameters for all studies and callers looking at obs biotins in IDRs  
Figure3C_Distribution-of-FPU-by-caller-parameters-and-study.pdf	|	Number of proteins in FPU categories in each study and by caller  
Figure3D_VSL2b_Biotins-in-IDR.vs.IDR-classes.violinPlots.pdf	|	Violin plots showing biotins in IDRs separated by FPU for caller VSL2b  
Figure4A_Num.idrs.vs.num.ptms_HEK293-vs-biotin.pdf	|	Correlation plot between number of PTMs and Biotins for HEK293 proteome (Geiger) and HEK203 biotinome  
Figure4B_Distribution-of-PTM-types-in-HEK-vs-in-Biotin-studies.pdf	|	Boxplots showing PTM counts by type in HEK293 proteome and biotinome  
Figure4C_Binomial-test-stats.txt	|	Bionomial test values and parameters for all studies and callers looking at obs PTMs in IDRs 
Figure4C_Distribution-of-PTMs-in.out.of.IDRs-by-study-VSL2b.pdf	|	Bar plot showing expected and observed rates of all types of PTMs for all studies  
Figure4D_VSL2b_PTMs-in-IDR.vs.IDR-classes.violinPlots.pdf	|	Violin plots showing all PTMs in IDRs separated by FPU for caller VSL2b  
Figure4D_VSL2b_Acytelation.vs.IDR-classes.violinPlots.pdf       |       Violin plots showing acetylation marks in IDRs separated by FPU for caller VSL2b  
Figure4D_VSL2b_Phosphorylation.vs.IDR-classes.violinPlots.pdf	       |       Violin plots showing phosphorylation marks in IDRs separated by FPU for caller VSL2b  
Figure4D_VSL2b_Sumoylation.vs.IDR-classes.violinPlots.pdf       |       Violin plots showing sumoylation marks in IDRs separated by FPU for caller VSL2b  
Figure4D_VSL2b_Ubiquitination.vs.IDR-classes.violinPlots.pdf       |       Violin plots showing ubiquitination marks in IDRs separated by FPU for caller VSL2b  
FigureS1A_VennDiagrams-for-biotinylation-data.pdf	|	Venn diagram showing overlap of proteins in biotinylation data  
FigureS1C_Biotin-GO-CC-enrichment-with-anc-terms_8_enricher-dotplot.pdf	|	Bubble plots showing top 8 GO enrichment terms (including ancestors) for "cellular component" for each of the 4 studies  
Biotin-GO-enrichment-with-anc-terms_8_enricher-dotplot.pdf	|	Bubble plots showing top 8 GO enrichment terms (including ancestors) for each of the 4 studies  
Biotin-GO-enrichment_8_enricher-dotplot.pdf	|	Bubble plots showing top 8 GO enrichment terms for each of the 4 studies. No ancestor terms used.  
Biotin-GO-CC-enrichment_8_enricher-dotplot.pdf	|	Bubble plots showing top 8 GO enrichment terms for "cellular component" for each of the 4 studies. No ancestor terms used.  
GO-KEGG-enrichment-with-whole-human-genome-bg.pdf	|	GO enrichment plotting using clusterProfiler's inbuilt code and whole human genome as background  
Effect-size-for-pairwise-comparison-of-PTMs-in-IDRs-VSL2b.pdf	|	Effect size plot based on Tukey's HSD for PTMs in IDRs  
Effect-size-for-pairwise-comparison-of-biotin-in-IDRs-VSL2b.pdf   |       Effect size plot based on Tukey's HSD for PTMs in IDRs  
Date_Time_Effect-size-for-pairwise-comparison-of-IDR-classes-and-parameters-VSL2b.pdf	   |       Effect size plot based on Tukey's HSD for multiple parameters and IDR classes using VSL2b calls for IDR
Date_Time_VSL2b_Pairwise-t-test-pvals.txt
Date_Time_VSL2b_Tukey-HSD-ci-pvals.txt

### Figures : Folder containing all main and supplementary Figures from the paper

#### SupplTables.xlsx :  
Supplementary tables referred to in the manuscript. Contain data that can be queried/used. Contains 5 datasheets and an "Index" sheet describing each of the data sheets.  
