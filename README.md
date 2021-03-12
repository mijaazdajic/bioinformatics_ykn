# Statistical Analyses 

This is a repository of the collections of R scritps used to anlayse and plot bioinformatics, chemical, and experimental data from lakes in the Yellowknife area in NWT, Canada.  The results, figures, and dicussion have been published in **TBD***

# Description of scripts
## [01_Import_data_DCM.R](01_Importing_data_DCM.R)

This script contains the necessary code to import multiple analyses files from ICP-MS mercury isotope anlysis, files for sediment moisture content, and isotope spike information.  The result from this script is a large dataframe with all the necessary information orgranized in a manner approprite for statistical anlayses and vizualizations.  The files that are saved in this script will be used in all subsequent scripts

## [02_Plotting_Incubations.R](02_Plotting_Incubations.R)

This script contains the theme for all figures in the manuscript as well as the code for plotting isotpe experiment data (methylation/demethyation).
The script 01 needs to be run before this so that incubation.join object is loaded.
This script contains the loop function to make multiple plots for each experiment and save each plot individually, as well as the isotope recovery calculations for each experiments using both MA300 and ICP-MS.

## [03_Rates and Figures.R](03_Rates_and_Figures.R)

This script contains the theme for all figures in the manuscript as well as the code to plot experiment data (mean) as well as figures with regression lines.
The script 01 needs to be run before this so that incubation.join object is loaded.
This script contains the loop function to make multiple plots for each experiment and save each plot individually.
Lastly, it contains the code to calculate the methylation and demethylation rate constants (Km/Kd), which are then saved in a CSV file.

## [04_Merging_Rate_Chem.R](04_Merging_Rate_Chem.R)

This script merges the calculated isoptoe methylation dataset with water chemistry dataset, the merged data is saved as an object called "mastersheet".
StepAIC analysis is used to identify which water chemistry variable best explains the %MeHg and methylation/demethylation rates.
The second part of the script is to seperate the two group of lakes (group 1 and 2), once again stepAIC analysis is used to identify which water chemistry variable best explains the methylation/demethylation rates, as well as the figure shown in Fig 3 of the manuscript.

## [05_Methylation_Rate_W_Interactions.R](05_Methylation_Rate_W_Interactions.R)

This script merges the calculated isoptoe methylation dataset with water chemistry dataset, the merged data is saved as an object called "mastersheet".
StepAIC analysis is used to identify which water chemistry variable best explains the %MeHg and methylation/demethylation rates WITH interactions between the independent variables.
The second part of the script is to seperate the two group of lakes (group 1 and 2), once again stepAIC analysis is used to identify which water chemistry variable best explains the methylation/demethylation rates WITH interactions between the independent variables, as well as the figure shown in Fig 3 of the manuscript.

## [06_Binary_regression_Kd.R](06_Binary_regression_Kd.R)

There were many lakes for which demethylation was undetectable, therefore we needed to find a way to accurately indenfity correlations/interactions with independent water chem varialbes.  The "mastersheet" lakes are seperated into two: one with zero value Kd and one with non-zero value Kd.
First, a stepwise regression is used to identify if any chem varialbes have significant correlation with the non-zero Kd subset of lakes.  
Second, a t-test was used to analyse if there are any significant differences in water chem between lakes with and without detectable Kd.

## [07_Importing_QIIME2.R](07_Importing_QIIME2.R)

This script as well as all scrip 08 and 09 use the results from the bioinformatics pipeline completed in QIIME2.  Please refer to [ykn_QIIME_analysis repository](https://github.com/mijaazdajic/ykn_QIIME_analysis) for the pipeline script.

This code is for importing QIIME2 analysis artifacs (i.e. SVs, metadata, tree, and taxonomy).  The bioinformatic pipeline results were "cleaned/filtered" and normalized before analyses and figures.  Next, microbial abundance was filtered and plotted as bar graphs and PCAs.


## [08_Constrained_Spatial_All_Samples_NoSUVA_NoPhos.R](08_Constrained_Spatial_All_Samples_NoSUVA_NoPhos.R)

This script is to analyse the microbial community with regards to chemical results.  QIIME2 artifacts were uploaded as well as the "mastersheet".
A constrained multivariate analysis is done on the microbial community with chemical results as the constraining variables.  
In addition, variability explained was calculated.

## [09_Constrained_Spatial_WITH_Km_Kd.R](09_Constrained_Spatial_WITH_Km_Kd.R)

This script is to analyse the microbial community with regards to chemical AND INCUBATION (Kd/Km) results.  QIIME2 artifacts were uploaded as well as the "mastersheet".
A constrained multivariate analysis is done on the microbial community with chemical AND INCUBATION (Kd/Km) results as the constraining variables.  
In addition, variability explained was calculated.

## [Figures](https://github.com/mijaazdajic/ykn_R_statistical_anlayses/tree/main/figures)

Folder contains the code for all figures in the MS, excluding Figure which was made in ArcMap.  Figures are all made in ggplot2.
