# 22q_spatial_transcriptomic
Analyses relating spatial patterns of gene transcription from Allen Human Brain Atlas (AHBA) to fMRI markers of 22q11.2 deletion syndrome (22qDel)

## Overview

## Analysis Steps
1. Generate parcellated BOLD metrics (e.g. RSFA) on hoffman
  * scripts in this repo are intended to be run locally but require inputs files that are best generated on a HPC cluster (e.g. hoffman2) due to storage and computational constraints. For code to be run on the cluster prior to these scripts, see: https://github.com/cschleifer/22q_hoffman 

2. Convert CAB-NP parcellation in CIFTI space to fsaverage5 surface and MNI volume space for input into abagen
  * [prep_cabnp_for_abagen.sh](22q_spatial_transcriptomic/prep_cabnp_for_abagen.sh) uses wb_command to convert CIFTI dlabel file 
