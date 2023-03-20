# ZITS: Zero-Inflated Time-Series Clustering Via Ensemble Thick-Pen Transform

This GitHub repository contains the R code for the paper titled "Zero-Inflated Time-Series Clustering Via Ensemble Thick-Pen Transform" 

by Minji Kim, Hee-Seok Oh and Yaeji Lim, Feb 2023.

The repository includes the following files:

* ``functions.R`` ; Main functions for ZITS clustering, including TPT pen shape and clustering algorithm.

* ``simulation_models.R`` ; Data generation for simulation models.

* ``simulation_analysis.R`` ; ZITS clustering and comparison methods for simulation data, which loads "functions.R" and "simulation_models.R".

* ``covid_analysis.R`` ; ZITS clustering and comparison methods for covid-19 data, which loads "functions.R".

* ``step_readData.R`` ; Data load for step count data.

* ``step_analysis.R`` ; ZITS clustering and comparison methods for step count data, which loads "functions.R" and "step_readData.R".

* ``eda_plots.R`` ; Draw figures for the paper using the step dataset.

The COVID-19 data is publicly available and therefore provided in this repository, but the step count data is not publicly available. If you have any requests regarding the dataset, please contact mkim5@unc.edu.
