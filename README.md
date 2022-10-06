# 22q_spatial_transcriptomic
Analyses relating spatial patterns of gene expression from Allen Human Brain Atlas (AHBA) to fMRI markers of 22q11.2 deletion syndrome (22qDel)

## Overview

## Analysis Steps
1. Generate parcellated BOLD metrics (e.g. RSFA) on hoffman
    * scripts in this repo are intended to be run locally but require inputs files that are best generated on a HPC cluster (e.g. hoffman2) due to storage and computational constraints. For code to be run on the cluster prior to these scripts, see: https://github.com/cschleifer/22q_hoffman 


2. Convert CAB-NP parcellation in CIFTI space to fsaverage5 surface and MNI volume space for input into abagen
    * [prep_cabnp_for_abagen.sh](prep_cabnp_for_abagen.sh) uses wb_command to convert CIFTI dlabel file 
    * outputs are in [CAB_NP_CONVERTED](CAB-NP/CAB_NP_CONVERTED) directory


3. Use abagen to parcellate AHBA data based on CAB-NP atlas
    * [abagen_ahba_cabnp.py](abagen_ahba_cabnp.py) uses abagen.get_expression_data to extract gene expression for each parcel
    * outputs are [CAB-NP_subcort_abagen_expression.csv](CAB-NP_surface_abagen_expression.csv) and [CAB-NP_subcort_abagen_expression.csv](CAB-NP_subcort_abagen_expression.csv)
      * surface and volume are computed separately -- surface using fsaverage5 space, volume using MNI. Expression is normalized separately within brain structures (cortex, cerebellum, subcort/brainstem) because computing volume and surface separately complicates normalizing across the whole brain 


4. Test relationship between PVALB/SST expression gradient and fMRI group difference (22qDel - healthy controls)
    * [22q_spatial_transcriptomic.R](22q_spatial_transcriptomic.R) reads parcellated AHBA data into R, along with subject level MRI measures in csv format (see Step 1)
    * requires wb_command for ciftiTools functions that read/plot MRI data. 
      * download: https://www.humanconnectome.org/software/get-connectome-workbench
      * script expects the workbench directory to be `/Applications/workbench/bin_macosx64` (either download to this location, or edit the path in the script)
 
