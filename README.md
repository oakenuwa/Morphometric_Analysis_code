## Contributors:
- Oghosa Honor Akenuwa
- Steve Abel
- Andreas Nebenfuhr

## Quick Start Guide:

This repository contains the analysis code to generate the morphometric parameters for the simulated and experimental actin networks. The analysis code is implemented in MATLAB and ImageJ.
We have also included folders with the actin positions and image files of sample simulated actin networks to allow users to see how the morphometric parameters are extracted.

To use and test the code:
1. Download the folder named `Simulated_actin_positions.zip`
2. Download the folder named `Simulated_actin_images.zip`
3. Run the ImageJ Macro named `Actin_Morphometrics.ijm` on the simulated actin images to get the _measured_ morphometric parameters
4. To generate the _ground-truth_ morphometric parameters, run the following codes in MATLAB in any order on the simulated_actin_positions.zip and keep the outputs in the current workspace:
	- Angle and Ordering parameters: `Morphometric_analysis_code_Angle_and_Ordering.m`
	* Density parameters: `Morphometric_analysis_code_Distance_and_Occupancy.m`
	+ Bundling: `Morphometric_analysis_code_LFB.m`
5. To perform pca analysis of the _measured_ and _ground-truth_ morphometric parameters in MATLAB, run `Morphometric_analysis_code_pca_analysis.m`.

Our simulation framework is based on [AFINES](https://github.com/Simfreed/AFINES.git) (Active FIlament NEtwork Simulation) developed by the Dinner group at the Univeristy of Chicago.
Details of the simulaton can be found in 
* Freedman, S. L., Banerjee, S., Hocky, G. M., & Dinner, A. R. (2017). [A Versatile Framework for Simulating the Dynamic Mechanical Structure of Cytoskeletal Networks.](https://doi.org/10.1016/j.bpj.2017.06.003)
* Akenuwa, O. H., & Abel, S. M. (2023). [Organization and dynamics of cross-linked actin filaments in confined environments.](https://doi.org/10.1016/j.bpj.2022.11.2944)

