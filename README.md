# maize-abiotic-stress

The repository contains the code used in the publication REFERENCE HERE. 

# Raw Image Data
Raw Image data is available at BioImage Archive: https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD2511

# file descriptions

bell_side_maize_07232022.ipynb - This jupyter notebook is used for testing and developing the PlantCV workflow used for RGB analysis of side-view images captured on the Bellwether Phenotyping Platform. 

VIS_SV_07232022.py - This python file contains the code for analyzing shape and color of side view RGB images captured on the Bellwether Phenotyping Platform. 

x-rite_color_matrix_k2.npz - This file is used as the reference set for color correction, and is used in VIS_SV_07232022.py

VIS_SV_07252022_config_run1.json - This configuration file is used for analyzing images in control temperature conditions using PlantCV. 

VIS_SV_07282022_config_run1.json - This configuration file is used for analyzing images in heat temperature conditions using PlantCV. 

maize_naive_bayes_pdfs_new.txt - This file contains the pixel classifications for naive bayes, and is used in VIS_SV_07232022.py. 
