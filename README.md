# SPLUS-Fornax-EmissionLine

This repository is dedicated to the emission line analysis in the Fornax cluster using S-PLUS images.
It also contains codes related to Legacy Survey, in order to compare results between both surveys.

Python in this repository:
* download_legacy images

  usage: python download_legacy_images.py --legacy table_name
  
  output: creates a folder "legacy_color_images" containing legacy color images in JPEG
* spatial_butterworth.py
  usage: from spatial_butterworth import butterworth
         butterworth(input_data).apply()
  input information: input_data can be a single image or a data cube. In case of a single image, 
  the user must add image_type='single', as the default is image_type='cube'
  output: input_data with butterworth filter applied

Main code (work in progress, not yet uploaded in this repository)
Input data: list of objects to be analyzed. It is necessary to have ID, RA, DEC, radius (image cutout size in pixels)
