# R-script-examples
Examples of R codes written or modified by Eunnara Cho

## Delta delta Ct method
This script was first written as an exercise in utilizing basic R packages (i.e., tidyverse) and exploring different statements (e.g. if, elseif, while). It was written to accommodate qPCR data from multiple experiments, each with a different set of controls for normalization.
The script is currently in use for quickly processing qPCR data and visualizing fold changes in gene expression. 

## Recurrent mutation analysis for Duplex Sequencing data

## ANOVA and post-hoc tests for MicroFlow assay data
This script processes the output of 96-well MicroFlow assay, a flow cytometry-based micronucleus (MN) assay. It calculates the standard deviation and error on MN frequency at each concentration, relative survival based on nuclei-to-counting bead ratio, and performs ANOVA followed by Dunnett or Tukey post-hoc analysis to determine the statistical significance of the changes in MN frequency in treated vs vehicle control samples. Relative survival and MN frequency are then plotted with statistical significance marked with asterisks. 
