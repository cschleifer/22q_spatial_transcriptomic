wb_command -gifti-help

# separate CAB-NP CIFTI parcels to single-hemisphere GIFTI files
cifti_in="/Users/charlie/Dropbox/PhD/bearden_lab/22q/analyses/striatum_thalamus_fc/ColeAnticevicNetPartition-master/data/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii"
metric_out_left="/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_L.func.gii"
metric_out_right="/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_R.func.gii"
wb_command -cifti-separate $cifti_in COLUMN -metric CORTEX_LEFT $metric_out_left -metric CORTEX_RIGHT $metric_out_right

for hemi in L R; do
    ## set up vars for -metric-resample to fsaverage5 
    # metric file to resample
    metric_in=/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_${hemi}.func.gii
    # path to reference surfaces
    surf_path=/Users/charlie/Dropbox/PhD/bearden_lab/22q/HCPpipelines-4.1.3/global/templates/standard_mesh_atlases/resample_fsaverage/
    # appropriate hemisphere and resolution from standard_mesh_atlases/resample_fsaverage (e.g.fs_LR-deformed_to-fsaverage.?.sphere.???.k_fs_LR.surf.gii)
    # confirmed 32k mesh by reading metric_out_left as a text file, Dim0="32492"
    current_sphere=${surf_path}/fs_LR-deformed_to-fsaverage.${hemi}.sphere.32k_fs_LR.surf.gii
    # use appropriate hemi and res (e.g fsaverage5_std_sphere.?.???.k_fsavg_?.surf.gii)
    new_sphere=${surf_path}/fsaverage5_std_sphere.${hemi}.10k_fsavg_${hemi}.surf.gii
    # specify output name *func.gii
    metric_out=/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_resample_fsaverage5_${hemi}.func.gii
    # use midthickness matching current_sphere
    current_area=${surf_path}/fs_LR.${hemi}.midthickness_va_avg.32k_fs_LR.shape.gii
    # use midthickness matching new_sphere
    new_area=${surf_path}/fsaverage5.${hemi}.midthickness_va_avg.10k_fsavg_${hemi}.shape.gii
    # check files exist
    #ls $metric_in
    #ls -d $surf_path
    #ls $current_sphere
    #ls $new_sphere
    #ls $current_area
    #ls $new_area
    # command to resample surfaces to fsaverage5
    wb_command -metric-resample $metric_in $current_sphere $new_sphere ADAP_BARY_AREA $metric_out -area-metrics $current_area $new_area
done



ls /Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/

cat /Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_resample_fsaverage5_L.func.gii| grep "Dim0"

# TODO: current resample method blurs borders (visible with random colos scheme in Surfice) 
# try creating label.gii from atlas dlabel.nii first
# CAB_NP_surface_orig_L.func.gii boundaries look good, but not resampled version
label_in="/Users/charlie/Dropbox/PhD/bearden_lab/22q/analyses/striatum_thalamus_fc/ColeAnticevicNetPartition-master/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii"
label_out_left="/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_L.label.gii"
label_out_right="/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_R.label.gii"
wb_command -cifti-separate $label_in COLUMN -label CORTEX_LEFT $label_out_left -label CORTEX_RIGHT $label_out_right


# label resample 
for hemi in L R; do
    ## set up vars for -metric-resample to fsaverage5 
    # metric file to resample
    label_in=/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_${hemi}.label.gii
    # path to reference surfaces
    surf_path=/Users/charlie/Dropbox/PhD/bearden_lab/22q/HCPpipelines-4.1.3/global/templates/standard_mesh_atlases/resample_fsaverage/
    # appropriate hemisphere and resolution from standard_mesh_atlases/resample_fsaverage (e.g.fs_LR-deformed_to-fsaverage.?.sphere.???.k_fs_LR.surf.gii)
    # confirmed 32k mesh by reading metric_out_left as a text file, Dim0="32492"
    current_sphere=${surf_path}/fs_LR-deformed_to-fsaverage.${hemi}.sphere.32k_fs_LR.surf.gii
    # use appropriate hemi and res (e.g fsaverage5_std_sphere.?.???.k_fsavg_?.surf.gii)
    new_sphere=${surf_path}/fsaverage5_std_sphere.${hemi}.10k_fsavg_${hemi}.surf.gii
    # specify output name *func.gii
    label_out=/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_resample_fsaverage5_${hemi}.label.gii
    # use midthickness matching current_sphere
    current_area=${surf_path}/fs_LR.${hemi}.midthickness_va_avg.32k_fs_LR.shape.gii
    # use midthickness matching new_sphere
    new_area=${surf_path}/fsaverage5.${hemi}.midthickness_va_avg.10k_fsavg_${hemi}.shape.gii
    # command to resample surfaces to fsaverage5
    wb_command -label-resample $label_in $current_sphere $new_sphere ADAP_BARY_AREA $label_out -area-metrics $current_area $new_area
done


# get volume parcels nii from CIFTI
# same cifti_in as surface step
cifti_in="/Users/charlie/Dropbox/PhD/bearden_lab/22q/analyses/striatum_thalamus_fc/ColeAnticevicNetPartition-master/data/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii"
# volume out name
volume_out="/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_vol_separated.nii"
wb_command -cifti-separate $cifti_in COLUMN -volume-all $volume_out

# create volume nii
#wb_command -label-to-volume-mapping 



