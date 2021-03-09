# Statistical Analyses 

This is a repository of the collections of R scritps used to anlayse bioinformatics, chemical, and experimental data from lakes in the Yellowknife area in NWT, Canada.  The results, figures, and dicussion has been published in ______

# Description of scripts
## [01_Import_data_DCM.R](01_Importing_data_DCM.R)

This script contains the necessary code to import multiple analyses files from ICP-MS mercury isotope anlysis, files for sediment moisture content, and isotope spike information.  The result from this script is a large dataframe with all the necessary information orgranized in a manner approprite for statistical anlayses and vizualizations.

  ### Files necessary
  isotope data: .csv file that contains raw ICP-MS results of mercury isotope concentrations from all experiments.  
  The file has the following headers: 
  - "Sample_Name"	
  - "Spike_Conc"	
  - "Volume"	
  - "Weight_sample"	
  - "202_Hg_ug_kg"	
  - "199_Hg_ug_Kg"	
  - "198_Hg_ug_kg"	
  - "Rep"	
  - "Day"	
  - "Month"
  - "Year"

## [02_Plotting_Incubations.R](02_Plotting_Incubations.R)


## [03_Rates and Figures.R](03_Rates and Figures.R)
