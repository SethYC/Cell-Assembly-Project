# Cell-Assembly-Project

These MATLAB scripts were used with cell assembly data using the methods outlined in Russo & Durstewitz (2017) on a dataset of aged and young rats performing a eyeblink conditioning task (see Schimanski, Lipa, & Barnes (2013) for original study on this dataset). 

`Automated_Cell_Assembly_Results_Extraction_epochs_probability_ver_Part_1.m` is run first on the original cell assembly results to create a table of all relevant data, then `Automated_Cell_Assembly_Results_Extraction_Part_2.m` uses the table to create plots comparing different groupings of results. `make_epoch_ANOVA_tables_probability_ver.m` was for preperation of the results into a format usable for IBM SPSS to conduct an ANOVA on. 


References:

Russo, E., & Durstewitz, D. (2017). Cell assemblies at multiple time scales with arbitrary lag constellations. Elife, 6, e19428.
Schimanski, L. A., Lipa, P., & Barnes, C. A. (2013). Tracking the course of hippocampal representations during learning: when is the map required?. Journal of Neuroscience, 33(7), 3094-3106.
