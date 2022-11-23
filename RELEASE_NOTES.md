# RELEASE 2022.01
In this release, portions of `atmos_param` that are unused in the `AM4 model` have been pruned.  This 2022.01 version reproduces all the CHECKSUMS and restart numerical data for AM4, CM4, ESM4, and SPEAR.  Users will however see that 
1.  the variable **restart_version** no longer exists in cg_drag.res.nc 
2.  the variable **doing_edt** no longer exists in physics_driver.res.nc
3.  the variable **tke** no longer exists in physics_driver.res.tile[1-6].nc

#### DIRECTORIES
In details, the following directories, referred to here as they are referred in the code, have been pruned:
- betts_miller 
- cloud_obs
- clouds
- cloud_zonal
- CLUBB
- diag_cloud
- diag_cloud_radiation
- diffusivity
- donner_deep
- dry_adj
- edt
- grey_radiation
- lin_cloud_microphysics
- my25_turb
- qe_moist_convection
- ras
- rh_clouds
- shallow_conv
- shallow_physics
- strat_cloud
- tke_turb
- two_stream_gray_rad 

### FILES
In addition, the following files (modules) have been pruned
- bulkphys_rad_mod
- cldwat2m_micro_mod
- donner_deep_clouds_W_mod
- micro_mg_mod
- morrison_gettelman_microp_mod
- simple_pdf_mod 

### NML VARIABLES
Removal of the abovementioned directories and files required the following nml options to be pruned:
<details><summary> cg_drag_nml</summary>
        calculate_ked <br/>
        num_diag_pts_ij <br/>
        num_diag_pts_latlon <br/>
        i_coords_gl <br/>
        j_coords_gl <br/>
        lat_coords_gl <br/>
        lon_coords_gl <br/>
        Bt_eq<br/>
        Bt_eq_width<br/>
</details>

<details><summary> cloud_rad_nml</summary>
        clubb_error <br />
        prog_ccn
</details>

<details><summary> convection_driver_nml</summary>
       do_limit_donner <br/> 
       do_unified_convective_cloure <br/>
       do_donner_before_uw <br/>
       use_updated_profiles_for_uw <br/>
       use_updated_profiles_for_donner <br/>
       only_one_conv_scheme_per_column <br/>
       force_donner_moist_conserv <br/>
       do_donner_conservation_checks <br/>
       do_donner_mca <br/> 
       conv_frac_max <br/> 
       cmt_mass_flux_source = 'donner', 'donner_and_ras', 'donner_and_uw', 'ras_and_uw', 'donner_and_ras_anduw'
       remain_detrain_bug <br/> 
       keep_icenum_detrain_bug <br/>
</details>

<details><summary> ls_cloud_driver_nml </summary>
        do_legacy_strat_cloud <br/>
        microphys_scheme = 'lin','morrison_gettelman', 'mg_ncar', 'ncar'
</details>

<details><summary> ls_cloud_macrophysics_nml </summary>
        use_updated_profiles_for_clubb
</details>

<details><summary> ls_coud_microphysics_nml </summary>
        lin_microphys_top_press <br/>
        override_liq_num <br/>
        override_ice_num <br/> 
        use_Meyers <br/> 
        use_Cooper <br/> 
        micro_begin_sec <br/> 
        min_precip_needing_adjustment
</details>

<details><summary> rotstayn_klein_mp_nml</summary>
       use_inconsistent_lh
</details>

<details><summary> moist_processes_nml</summary>
        do_mca <br/> 
        do_ras <br/>
        do_donner_deep <br/> 
        do_dryadj <br/> 
        do_bm <br/> 
        do_bmmass <br/> 
        do_bmomp <br/> 
        do_simple <br/> 
        do_rh_clouds <br/> 
        include_donmca_in_cosp
</details>

<details><summary> physics_driver_nml </summary>
        do_clubb <br/>
        donner_meso_is_largescale <br/> 
        do_grey_radiation <br/> 
        R1, R2, R3, R4 <br/> 
        l_host_applies_sfc_fluxes
</details>

<details><summary> cloud_spec_nml </summary>
        ignore_donner_cells <br/>
        cloud_type_form = 'rh', 'deep', 'stratdeep', or 'stratdeepuw', 'deepuw'
</details>

<details><summary> microphys_rad_nml </summary>
        do_orig_donner_stoch <br/> 
        ignore_donner_cells
</details>
    
<details><summary> uw_conv_nml </summary>
        use_turb_tke
</details>
  
<details><summary> vert_turb_driver_nml </summary>
        do_shallow_conv <br/> 
        do_mellor_yamada <br/> 
        do_tke_turb <br/> 
        do_diffusivity <br/> 
        do_molecular_diffusion <br/>
        do_edt <br/> 
        do_simple
</details>
