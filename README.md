# Zymo-Oxygen
Systems-level Analysis of Oxygen Exposure in Zymomonas mobilis: Implications for Isoprenoid Production

This repository contains the raw data used to generate the metabolomics data shown in Figures 2, 3, 6, 7, 8, S2, S3, S5, and S7 of Martien et al. 2019.

The folders Alongside Proteomics, Alongside Transcriptomics, Alongside Extracellular, and For IDP/DMADP contain the raw data from four seperate iterations of the two-hour oxygen exposure timecourse, as described in Materials and Methods subsection "Two-hour oxygen exposure timecourse". Each folder contains:
1) An excel file containing the growth curve for that iteration - suffix "growth_curve"
2) A tab-delimited file containing chosen OD correction factors for each sample based on the growth curves - prefix "OD"
3) A tab-delimited file with extracted signal intensites for each metabolite - suffix "metabolites"
4) A folder titled Maven Files containing the mzxml files used to visualize metabolite peaks and extract signal intensity values using MAVEN (the folder For IDP-DMAPD has two sets - one with the complete scan and one with the suffix "IPP" that only scans for HMBDP and IDP/DMADP) 

The folder Cadmium contains:
1) An excel file containing the growth curve during cadmium treatment- suffix "growth-curve"
2) A tab-delimited file containing chosen OD correction factors for each sample based on the growth curves - prefix "OD"
3) A tab-delimited file with extracted signal intensites for each metabolite - suffix "metabolites"
4) A folder titled Maven Files containing the mzxml files used to visualize metabolite peaks and extract signal intensity values using MAVEN. 
There are two metabolites files and two Maven folders: one from treatment with 0.15 mM Cd2+, and one from 0.05 mM Cd2+ treatment.

The folder 2hr Extracellular contains:
1) A tab-delimited file containing the extracted signal intensities for each metabolite from extracellular extraction.
2) Calibration curves used to quantify extracellular gluconate and MEcDP levels - suffix "calibration_curve"
3) A folder titled Maven Files containing the mzxml files used to visualize metabolite peaks and extract signal intensity values using MAVEN.
These data were collected from the same biological samples as the intracellular data found in "Alongside Extracellular"

The folder 24hr Extracellular contains:
1) An excel files with the growth curve for the 24 hour oxygen exposure time course
2) Excel files with NMR signal intensities and absolute quantitation of Ethanol and Acetate from NMR data - suffix "quantitation"
3) An excel file containing the LC-MS signal intensites and absolute quantitation of glucose, gluconate, and pyruvate by LC-MS using 13-labeled internal standards - "lebeled extracellular-metabolites"
4) An excel file containing the four-point calibration curves used to quantify lactate, succinate, and shikimate by LC-MS - "lact-succ-shik_calibration-curve"
5) A tab-delimited file containing all LC-MS signal intensites for extacellular metabolites without internal labeled standards
6) An excel file containing the quantified values for all significant detected metabolites and calculations to obtain the percent of glucose consumption represented by each extracellular metabolite. "Percent_of_glucose"



