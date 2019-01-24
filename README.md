# Zymo-Oxygen
Systems-level Analysis of Oxygen Exposure in Zymomonas mobilis: Implications for Isoprenoid Production

This repository contains the raw data used to generate the metabolomics data shown in Figures 2, 3, 6, 7, 8, S2, S3, S5, and S7 of Martien et al. 2019.

The folders Alongside Proteomics, Alongside Transcriptomics, Alongside Extracellular, and For IDP/DMADP contain the raw data from four seperate iterations of the two-hour oxygen exposure time-course, as described in Materials and Methods subsection "Two-hour oxygen exposure time-course". Each folder contains:
1) An excel file containing the growth curve for that iteration - suffix "growth_curve"
2) A tab-delimited file containing chosen OD correction factors for each sample based on the growth curves - prefix "OD"
3) A tab-delimited file with extracted signal intensites for each metabolite - suffix "metabolites"
4) A folder titled Maven Files containing the mzxml files used to visualize metabolite peaks and extract signal intensity values using MAVEN (the folder For IDP-DMAPD has two sets - one with the complete scan and one with the suffix "IPP" that only scans for HMBDP and IDP/DMADP) 
