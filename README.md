# FutureStreams
Scripts used for PCR-GLOBWB 5arcmin discharge and water temperature simulations and data paper. 

Derived_variables:
Python code used to create the ecologically relevant derived variables. 
Also scripts to create a mask for pixels with unrealistic values from derived variables (mask_unrealistic_values.py and maskFunctions.py), and scripts to mask unrealistic values from weekly Q and WT (mask_unrealistic_values_weekly.py and maskFunctions_weekly.py). A user can set their preferred threshold values in the maskFunctions* script. 

Figures:
R code used to create Figure 4. R code for the other figures is available through https://github.com/vbarbarossa/futurestreams_figures. 

PCRGLOBWB_config:
Contains an example of a configuration file used in the PCR-GLOBWB simulations. In this example, the ipsl - rcp6.0 simulation is started, using initial conditions from the end of the ipsl historic simulation. This script is for one of the 53 regions / catchments for which the model is run. 

Warming_levels:
Contains a text file with the year in which the 30-year running mean of air temperature of each GCM- RCP-combination crosses a certain warming level (1.5, 2.0, 3.2 and 4.5 degrees, as used in Barbarossa et al. 2021, see their Supplementary Information). 
