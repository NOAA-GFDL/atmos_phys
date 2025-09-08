module matrix_gfdl
    use platform_mod
    use mpp_mod,               only : input_nml_file, mpp_sync, mpp_chksum
    use fms_mod,               only : file_exist, close_file,&
        open_namelist_file, check_nml_error, &
        write_version_number, &
        error_mesg, &
        warning, &
        fatal, note, &
        lowercase, &
        mpp_pe, &
        mpp_root_pe, &
        stdlog, stdout, &
        mpp_clock_id, &
        mpp_clock_begin, &
        mpp_clock_end, &
        clock_module, &
        uppercase
    use tracer_manager_mod,    only : get_tracer_index,   &
        get_number_tracers, &
        get_tracer_names,   &
        get_tracer_indices, &
        adjust_positive_def, &
        query_method, &
        no_tracer
    use  field_manager_mod, only : model_atmos, &
        parse
    !use ieee_arithmetic, only : ieee_is_nan !not support inf 
    use ieee_arithmetic
    use  diag_manager_mod, only  : send_data, register_diag_field
    use  time_manager_mod, only  : time_type
    use  constants_mod, only     : pi, grav, rdgas, dens_h2o, wtmair, avogno, &
        pstd_mks
    !use atmos_tracer_utilities_mod, only : sedimentation_velocity, &
    !    sedimentation_flux, &
    !query_caerosol_wetdep_param
    use sat_vapor_pres_mod, only : compute_qs  
    use aero_npf, only : npfrate, steady_state_h2so4
    use aero_wet, only : aero_kohler
    use aero_condens, only : setup_kci
    use aero_coag_config, only : setup_coag_tensors, &
        giklq_control, dikl_control, &
        nm, prod_index, giklq, &
        dikl, dij, ndikl,citable, &
        mode_name, nmass_spcs   
    use aero_coag,  only : setup_kij_diameters, setup_kij_tables, get_kbarnij
    implicit none
    private

    public matrix_init, matrix_source_type,set_matrix_source, matrix_run
    public query_matrix_pop, query_matrix_info
    public query_pop_number, query_pop_dg_dry, query_pop_mspcs, query_pop_sigma
    public query_pop_kappa
    public matrix_query_dynamic_wetdep_param_2d, matrix_query_dynamic_wetdep_param_2d_mass
    public grid_wetdep_fic_uw, grid_wetdep_fic_uw_mass
    public query_tracer_in_pop
    interface set_matrix_source
        module procedure set_matrix_source_2d
        module procedure set_matrix_source_3d
    end interface set_matrix_source

    logical :: matrix_module_init = .false.
    real, parameter :: pi6 = pi/6.0
    integer, parameter :: nb_tracer_max = 9 !maxium number of tracers in a population, i.e. in mxx
    integer :: id_rh,id_kc, id_dmdt_h2so4_tot_cond_npf, id_dndt_npf, id_dmdt_h2so4_npf,id_cond_sink
    integer :: id_h2so4_source, id_pwt, id_zhalf, id_so4_emis !id for budget analysis 
    !define clocks for different processes
    integer :: ini_clock = 0 !initializetion
    integer :: npf_clock = 0 !new particle formation clock
    integer :: condgrow_clock = 0 !condesational growth
    integer :: hygrow_clock = 0 !hygroscopic growth
    integer :: hygrow_sub_clock1 = 0
    integer :: hygrow_sub_clock2 = 0
    integer :: hygrow_sub_clock3 = 0
    integer :: hygrow_sub_clock4 = 0
    integer :: dry_clock = 0
    integer :: coag_clock = 0 ! coagulation growth
    integer :: coag_sub_clock1 = 0
    integer :: coag_sub_clock2 = 0
    integer :: coag_sub_clock3 = 0
    integer :: coag_sub_clock4 = 0
    integer :: coag_sub_clock5 = 0
    integer :: coag_sub_clock6 = 0
    integer :: coag_sub_clock7 = 0
    integer :: dmodal_clock = 0 !inter-modal transfer
    integer :: aqso4_clock = 0 ! partitioning clock 
    integer :: matrix_run_clock = 0 ! test total run clock
    ! matrix source type
    type :: matrix_source_struct
        integer :: p_h2so4 = 1 !gaseous phase h2so4 to be used in nucleation scheme in akk
        integer :: e_so4   = 2 !so2 direct emission: 2.5% sulfur emission to be converted to sulfate mass in acc
        !integer :: e_dust  = 3 !dust emission to be distributed in dd1 and dd2
        integer :: e_ss    = 4 !sea salt emission to be distributed in ssa and ssc
        !integer :: e_oc    = 5 !organic carbon emission to be distributed in oc1 and oc2
        !integer :: e_bc    = 6 !black carbon emission to be distributed in bc1 and bc2     
        integer :: e_soa   = 7 !extra emission based on gfdl model
        integer :: p_aqso4 = 8 !aqueous so4 production
        integer :: e_oc_phob    = 9
        integer :: e_oc_phil    = 10
        integer :: e_bc_phob    = 11
        integer :: e_bc_phil    = 12
        integer :: e_ss_acc     = 13
        integer :: e_ss_coars   = 14
        integer :: e_dust_acc   = 15
        integer :: e_dust_coars = 16
        integer :: u_kg_m2_s = 1 !unit index for conversion
        integer :: u_vmr_s   = 2 !unit index for conversion
        integer :: u_mmr_s   = 3 !unit index for conversion
    end type matrix_source_struct

    type(matrix_source_struct) :: matrix_source_type

    !----------------------------------------------------------------
    !	-------------------type matrix_tracer-------------------
    !   ...define tracer type and encorporate tracer properties
    ! 	...tracer description in field table:
    !		"matrix_parameters", "acc", "distribution=lognormal, &
    !             & distribution_index=1,
    !		& sigma=1.8, dgn=0.068e-6, kappa=0.507, &
    !		&dens=1770, type=mass,spec=sulf"
    !----------------------------------------------------------------
    type :: matrix_tracer
        logical :: is_active = .false.
        character*32 :: type = "" ! mass or number
        character*3  :: pop = "" ! population name: akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
        integer :: pop_index = -1
        character*32 :: spec = ""! mass species name: sulf/dust/seas/ocar/bcar/alwc/gsfa...
        character*32 :: distribution = "" ! lognormal or weibull
        integer      :: distribution_index = -1
        real :: sigma   = -1 ! sigma of log-normal distribution
        real :: lnsigma  = -1 
        real :: lnsigma2 = -1 
        real :: dp0 = -1 !volume mean diameter of the current mode, unit m
        real :: dgn = -1 ! unit m: geometric diameter for initial log-normal distribution, only used for emission
        real :: kappa = 0 ! hygroscopicity facor
        real :: dens  = 0! unit kg/m3: density of specific species
        character*32 :: name,units !get_tracer_names (model_atmos, tracer_index, name= mt%name,  units = mt%units)
        real, allocatable    :: source(:,:,:) !sources in unit of µg/m3/s (mass tracers) or #/m3/s
        real, allocatable    :: value_in_matrix(:,:,:) !calculate tracer values in matrix unit
        logical :: has_emission = .false. !if the field table has distribution parameters, then has emission, otherwise, no-emission
        integer :: id_tracer_emis = -1 !diagostic emission ids, 3d µg/m3/s
        integer :: id_tracer_emis_source = -1
        integer :: id_tracer_setl = -1 !sedimentation flux, 2d, kg/m2/s or mol/m2/s
        integer :: id_tracer_col = -1 !dry deposition ids, 2d kg/m2/s
        integer :: id_tracer_tsource = -1 !matrix total source for tracer, budget test 
        integer :: id_tracer_intermodal = -1 !intermodal transfer, only happens in one direction, not reverse
        integer :: id_tracer_setl_3d = -1 !tendency from rt_dt_matrix, mmr/s
        integer :: id_tendency_matrix = -1 !matrix tendency for test, 3d ug/m3/s, #/m3/s
        integer :: id_tendency_gfdl = -1 !tendency in gfdl unit, mmr/s or vmr/s
        integer :: id_value_bf_max = -1 !around line 457, before max(value,0)
        integer :: id_value_af_max = -1 !around line 457, max(value,0)

    end type matrix_tracer

    !----------------------------------------------------------------
    !	-------------------type matrix_pop-------------------
    !   ...define population type (contain multi tracers)
    !	... to record what tracers are within a specific population
    !   ... if a tracer is in this population, record its index in the tracer array
    !   ...typical components in a population: n, m_spec1, m_spec2, ....
    !----------------------------------------------------------------
    type :: matrix_pop
        integer :: nb_tracer_pop = 0 !number of tracers in current population
        character*32 :: name = "" !name of the population
        integer :: tracer_index(nb_tracer_max) = -1
        logical :: has_emission(nb_tracer_max) = .false.
        real    :: def_dens = 1000.0
        real    :: def_kappa = 0.01
        real    :: def_dg = 1.0e-7
        integer :: i_n = -1 ! index of number in matrix_tracer array
        integer :: i_msulf = -1 ! index of sulfate mass in matrix_tracer array
        integer :: i_mdust = -1
        integer :: i_mseas = -1
        integer :: i_mocar = -1
        integer :: i_mbcar = -1
        integer :: i_mwate = -1
        integer :: i_mammo = -1
        integer :: i_mnitr = -1
        integer :: i_mgsfa = -1
        real :: sigma = 1.8 !gemetric sigma for population, temperary set for 1.8
        real :: rh_deliquescence = 0
        real :: rh_crystallization = 0
        real, allocatable    :: dg_dry(:,:,:) !unit: m, geometric diameter
        real, allocatable    :: dg_wet(:,:,:) !unit: m, geometric diameter
        real, allocatable    :: mass_dry(:,:,:) !unit: ug/m3
        real, allocatable    :: dens_wet(:,:,:) !unit: kg/m3
        real, allocatable    :: dens_dry(:,:,:) !unit: kg/m3
        real, allocatable    :: vol_dry(:,:,:) !must be double precision 
        real, allocatable    :: kappa_pop(:,:,:) !average hygroscopicity of the population
        integer :: id_dg_dry = -1
        integer :: id_dg_wet = -1
        integer :: id_pop_condens = -1
        integer :: id_pop_condens_source = -1
        integer :: id_pop_coag_source = -1
        integer :: id_pop_coag_pn = -1 !coag production of number [#/m3/s]
        integer :: id_pop_coag_ln = -1 !coag loss of number, 3d [#/m3/s]
        integer :: id_pop_coag_pmsulf = -1 !coag production of sulf mass [ug/m3/s]
        integer :: id_pop_coag_lmsulf = -1 !coag loss of sulf mass [ug/m3/s]
        integer :: id_pop_coag_pmocar = -1  !coag production of ocar mass [ug/m3/s]
        integer :: id_pop_coag_lmocar = -1  !coag loss of ocar mass [ug/m3/s]
        integer :: id_pop_coag_pmbcar = -1  !coag production of bcar mass [ug/m3/s]
        integer :: id_pop_coag_lmbcar = -1  !coag loss of bcar mass [ug/m3/s]
        integer :: id_pop_coag_pmdust = -1  !coag production of dust mass [ug/m3/s]
        integer :: id_pop_coag_lmdust = -1  !coag loss of dust mass [ug/m3/s]
        integer :: id_pop_coag_pmseas = -1  !coag production of seas mass [ug/m3/s]
        integer :: id_pop_coag_lmseas = -1  !coag loss of seas mass [ug/m3/s]
        integer :: id_pop_aqso4_partition = -1 !aqueous so4 partition to population [ug/m3/s]
    end type matrix_pop
    real,allocatable :: debug_num_domain(:,:,:,:)
    real, allocatable :: debug_mass_domain(:,:,:,:,:)
    !----------------------------------------------------------
    !output coagulation related variables
    integer :: id_mjq11 = -1, id_mjq21 = -1, id_kbar0_11 = -1, id_kbar3_11 = -1, id_kbar0_12 = -1
    integer :: id_kbar3_12 = -1, id_kbar0_21 = -1, id_kbar3_21, id_kbar0_22 = -1, id_kbar3_22 = -1 
    integer :: id_ri_1 = -1, id_ri_2 = -1, id_li_1 = -1, id_li_2 = -1
    integer :: id_bi_1 = -1, id_bi_2 = -1, id_half_knn_1 = -1, id_notself_knn_1 = -1
    integer :: id_rim_11 = -1, id_rim_21 = -1, id_fi_1 = -1, id_fi_2 = -1, id_lim_11 = -1, id_lim_21 = -1
    integer :: id_half_knn_2 = -1, id_notself_knn_2 = -1
    !end of output coagulation related variables
    !-----------------------------------------------------------
    integer :: id_aqso4_rate = -1 !diagnose total aqso4 received 
    type(matrix_tracer), allocatable :: matrix_all_tracer(:) ! define tracer array for all tracers
    type(matrix_pop),    allocatable :: matrix_all_pop(:)  ! defined population array for all population
    integer :: npop = 13 !maximum number of population in matrix
    !integer :: ntracer !maximum number of population and tracers in matrix
    integer :: i_akk = -1, i_acc = -1, i_dd1 = -1, i_dd2 = -1 ! index of population in matrix
    integer :: i_ssa = -1, i_ssc = -1, i_oc1 = -1, i_oc2 = -1 ! if exit, index >= 1, otherwise = -1
    integer :: i_bc1 = -1, i_bc2 = -1, i_mxa = -1, i_mxc =-1, i_ext = -1
    integer :: pop_active = 0 ! number of active population
    integer, allocatable :: i_pop(:)
    integer :: i_mw_h2so4 = 1, i_mw_so4 = 2, i_mw_dust = 3, i_mw_ss = 4
    integer :: i_mw_oc = 5, i_mw_bc = 6, i_mw_soa = 7
    integer :: nmspcs = 5 ! number of mass species
    real, dimension(7) :: matrix_molecular_weight = [98.07848, 96.0, 135.0, 58.5, 12.0, 12.0, 12.0]
    ! assign values to each element using the indices
    !matrix_molecular_weight(i_mw_h2so4) = 98.07848
    !matrix_molecular_weight(i_mw_so4) = 96.0
    !matrix_molecular_weight(i_mw_dust) = 135.0
    !matrix_molecular_weight(i_mw_ss) = 58.5
    !matrix_molecular_weight(i_mw_oc) = 12.0
    !matrix_molecular_weight(i_mw_bc) = 12.0
    !matrix_molecular_weight(i_mw_soa) = 12.0
    real, allocatable :: aqso4(:,:,:) !passing in: aqueous so4 production
    real, allocatable :: p_h2so4_rate(:,:,:)
    integer, parameter :: dist_lognormal = 1, dist_weibull = 2 !distribution index for lognormal and weibull
    integer :: ntrace 
    !integer :: tracer_index,  n, ntrace, nt, ntt
    !    integer :: tracer_index, ntrace
    !    integer :: ierr, io, logunit, verbose, unit
    !wet_deposition_const_param
    real :: cdust_wetdep_param(6) = 0.
    real :: cssalt_wetdep_param(6) = 0.
    integer :: nsphum, nh2so4
    logical :: do_matrix = .false., do_coag = .false., do_intermodal_transfer = .false.
    logical :: do_aqso4 = .false., do_sedimentation = .false.
    character(len=32) :: matrix_configuration = 'debug'
    character(len=32) :: coag_configuration = 'full'
    character(len=32) :: intermodal_configuration = 'full'
    character(len=7), parameter :: module_name = 'matrix'
    real :: rh_cap = 0.97
    integer :: dij_flag = 0. !used for print dij table for testing
    integer :: k_grid_size
    ! matrix_configuration is the version used in the calculation
    ! curretly matrix_configuration = debug is used to test toy model
    namelist /matrix_nml / do_matrix, matrix_configuration, & 
        do_coag, coag_configuration, &
        do_intermodal_transfer, intermodal_configuration, &
        do_aqso4, do_sedimentation, rh_cap
contains

    !----------------------------------------------------------------
    !
    !                           subroutine matrix_run
    !important: this subroutine must be after atmos_sox_chem etc, where sources are properly set-up 
    !the mass emission of different matrix_tracers are linked to different sources
    !in other f90 files, e.g. in atmos_sox_chem of atmos_sulfate.f90
    ! the source will be initialized by 0, and updated at every timestep 
    !----------------------------------------------------------------
    subroutine matrix_run(r, pfull, rh, t, dt, pwt, zhalf, rdt_matrix, time, time_next, is,ie,js,je, lon, lat, kbot)
        integer(kind=i8_kind) :: chksum_val
        real, intent(in) :: r(:,:,:,:)
        real, intent(in),    dimension(:,:)           :: lon, lat
        real, intent(in) :: dt !timestep
        integer, intent(in), optional :: kbot(:,:) ! index of bottom level
        integer, intent(in) :: is, ie, js, je ! boundaries of physical window
        type(time_type),  intent(in) :: time, time_next
        real, intent(out) :: rdt_matrix(:,:,:,:) !tendency calculated from matrix
        real, intent(in) :: pfull(:,:,:), rh(:,:,:), t(:,:,:), pwt(:,:,:), zhalf(:,:,:) !t is temperature
        real :: xnh3(size(r,1),size(r,2),size(r,3)), fland(size(r,1),size(r,2),size(r,3)) !xl
        ! assign number rate to each number tracers, the mass emission has been linked to different tracers in other files
        real :: rt(size(r,1),size(r,2),size(r,3),size(r,4)) !local variable to record updated values of tracers, in gfdl unit
        real :: rt_dt_setl(size(r,1),size(r,2),size(r,3),size(r,4)) !tendency from sedimentation
        !local variable to record absolute h2so4 mass change in ug/m3, positive
        real :: dm_h2so4(size(r,1),size(r,2),size(r,3)), dmdt_npf(size(r,1),size(r,2),size(r,3)), dndt_npf(size(r,1),size(r,2),size(r,3))
        real :: dmdt_h2so4_tot_cond_npf(size(r,1),size(r,2),size(r,3)), dmdt_h2so4_npf(size(r,1),size(r,2),size(r,3)) 
        integer :: n,mw  !local variables
        real :: m_akk_emis_source(size(r,1),size(r,2),size(r,3)), n_akk_emis_source(size(r,1),size(r,2),size(r,3))
        real :: kci_coef_pop(npop, size(r,1),size(r,2),size(r,3)), kci_aeq1_pop(npop, size(r,1),size(r,2),size(r,3)) !unit: m3/s
        logical :: used
        real :: xh2so4_nucl(size(r,1),size(r,2),size(r,3))
        integer :: ih2so4_path(size(r,1),size(r,2),size(r,3))
        real:: kc(size(pfull,1),size(pfull,2),size(pfull,3)) !kc: total condensation sink (1/s) for all population in a grid
        real:: pq_growth(size(pfull,1),size(pfull,2),size(pfull,3))
        real:: number_pop(size(pfull,1),size(pfull,2),size(pfull,3))
        real:: cond_pop(npop,  size(r,1),size(r,2),size(r,3))
        real:: emis_tmp(size(pfull,1),size(pfull,2)), col_tmp(size(pfull,1),size(pfull,2)), cond_tmp(size(pfull,1),size(pfull,2))
        integer :: it,jt,kt,i,j,k, nct !for loop use
        real :: li(npop-1,size(r,1),size(r,2),size(r,3)) !loss term of number due to intermodal coagulation [#/m^3/s]
        real :: ri(npop-1,size(r,1),size(r,2),size(r,3)) !production terms due to intermodal coagulation. [#/m^3/s]
        real :: bi(npop-1,size(r,1),size(r,2),size(r,3)) !loss terms due to intermodal coagulation. [1/s]
        real :: fi(npop-1,size(r,1),size(r,2),size(r,3)) !mass loss coefficient
        real :: lim(npop-1,nmspcs,size(r,1),size(r,2),size(r,3)) !mass loss term of number due to intermodal coagulation [ug/m^3/s]
        real :: rim(npop-1,nmspcs,size(r,1),size(r,2),size(r,3)) !mass production terms due to intermodal coagulation. [ug/m^3/s]
        real :: lim_spec(nmspcs,size(r,1),size(r,2),size(r,3)), rim_spec(nmspcs,size(r,1),size(r,2),size(r,3))
        integer :: nspecc, ig, jg, kg, kd
        real :: total_aero_num(size(r,1),size(r,2),size(r,3))
        real :: kci1, kci2, fac !debug use
        real :: pop_vset(npop-1, size(r,1),size(r,2),size(r,3))!settling velocity for population
        real :: kbar0_ij(npop-1,npop-1, size(r,1), size(r,2), size(r,3)) !mode averaged coagulation coefficient [m3/s]
        real :: mjq(npop-1,nmspcs,size(r,1),size(r,2),size(r,3)) !mass of species in a single particle: mass/number = ug/particle
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        rt = 0.
        dm_h2so4 = 0.
        rdt_matrix = 0.            
        dndt_npf = 0.
        dmdt_npf = 0.
        pq_growth = 0.
        dmdt_h2so4_tot_cond_npf = 0.
        dmdt_h2so4_npf = 0.
        xnh3 = 0.
        fland = 0.
        ih2so4_path = 0
        cond_pop = 0.
        kc = 0. 
        ri = 0.
        bi = 0.
        fi = 0.
        li = 0.
        rim = 0.
        lim = 0.
        kbar0_ij = 0.
        mjq =0.
        m_akk_emis_source = 0.
        n_akk_emis_source = 0.
        kci_coef_pop = 0.
        kci_aeq1_pop = 0.
        xh2so4_nucl = 0.
        number_pop = 0.
        emis_tmp = 0.
        col_tmp = 0.
        cond_tmp = 0.
        pop_vset = 0.
        kd = size(pfull,3)

        call mpp_clock_begin(matrix_run_clock)

        if (do_matrix) then
            rt = r
            
            !------------------------------------------------------------------
            !               step-1: settling before
            !
            !------------------------------------------------------------------
            rt_dt_setl = 0.
            !because need to retrive the dens, kappa, and diameter information, hence, need to read matrix first
            call set_matrix_value(rt, pwt, zhalf,is,ie,js,je)
            call matrix_dry_diameter(pfull,is,ie,js,je) !assign dry pop properties
            call matrix_wet_diameter(rh, t,is,ie,js,je) !assign wet pop properties, t is temperature
            if (do_sedimentation) then
                do n=1,npop
                if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
                    .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk") )  then !exclude ext, akk, oc1 and bc1
                    call rt_dt_from_sedimentation(rt, n, pfull, dt, t, pwt, is, ie, js, je,  rt_dt_setl, time_next, rh, kbot)
                endif
                enddo
            endif
            !if (mpp_root_pe().eq.mpp_pe()) then
            if (any( (rt_dt_setl >huge(0.0)) .or. (rt_dt_setl < -huge(0.0)) )) then
                call error_mesg('matrix_run check 1', 'infinity exist in rt_setl', fatal)
            endif
            !endif

            !DEBUG: check rt_dt_setl relevant to rt/dt, as don't want to have negative values
            !if (mpp_root_pe().eq.mpp_pe()) then
            do n=1,ntrace
            if (matrix_all_tracer(n)%is_active) then
                if (any(rt(:,:,:,n)+rt_dt_setl(:,:,:,n)*dt < 0)) then
                    write(*,*) 'negative value exist after settling, minval is ', matrix_all_tracer(n)%name, &
                        minval(rt(:,:,:,n)+rt_dt_setl(:,:,:,n)*dt)
                endif
            endif
            enddo
            !endif
            
            !------------------------------------------------------------------
            !                   update rt into matrix_processes
            !------------------------------------------------------------------
            rt = rt + rt_dt_setl*dt

            !after here, will redo matrix all processes

            !------------------------------------------------------------------
            !               step 0: get  matrix tracer values 
            !               and convert from gfdl to matrix unit
            !------------------------------------------------------------------
            call mpp_clock_begin(ini_clock)
            !initialize matrix species values at the beginning of the current step
            !the matrix_value(is:ie, js:je,:) is reset from this subroutine, 1e-17 for minimum n. 1e-32 for maximum n
            !the new module variable value matrix_value is updated by the rt passed in values
            call set_matrix_value(rt, pwt, zhalf,is,ie,js,je) ! assign matrix_tarcer%values_in_matrix, value in matrix unit
            !-----------------------------------------------------------------------------
            !               step 1: calculate tracer number sources from direct emission (in matrix unit)
            !                       note: this step doesn't include new particle formation
            !    note: the mass sources have been already updated before matrix_run
            !------------------------------------------------------------------------------
            !matrix_value not updated, only source get asigned
            call set_matrix_emis_number(is,ie,js,je) !update tracer source of  num_rate 
            call mpp_clock_end(ini_clock)

            !!-----------------------------------------------------------------------------
            !!  test order (0813_af)  step *: hygroscopic growth: calculate dry/wet particle diameter
            !!------------------------------------------------------------------------------
            !note: the calculation of dry and wet diamter calculation is just a side-prodcut from number/mass
            !number and mass information is not modified
            call mpp_clock_begin(dry_clock)
            call matrix_dry_diameter(pfull,is,ie,js,je) !assign dry pop properties
            call mpp_clock_end(dry_clock)
            call mpp_clock_begin(hygrow_clock)
            call matrix_wet_diameter(rh, t,is,ie,js,je) !assign wet pop properties, t is temperature
            call mpp_clock_end(hygrow_clock)

            !send data of d_wet, d_dry for test: ps. dg_dry, dg_wet in matrix_pop unit is m, only for send_data change to um
            do n = 1, npop
            if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
                if (matrix_all_pop(n)%id_dg_dry > 0) then
                    used = send_data (matrix_all_pop(n)%id_dg_dry, matrix_all_pop(n)%dg_dry(is:ie,js:je,:)*1e6, time_next, &
                        is_in=is,js_in=js, ks_in = 1) !in unit µm
                endif
                if (matrix_all_pop(n)%id_dg_wet > 0) then
                    used = send_data (matrix_all_pop(n)%id_dg_wet, matrix_all_pop(n)%dg_wet(is:ie,js:je,:)*1e6, time_next, &
                        is_in=is,js_in=js, ks_in = 1) !in unit µm
                endif

            endif
            end do

            !------------------------------------------------------------------------------------
            !               step *: test sedimentation at before
            !-------------------------------------------------------------------------------------
            !rt_dt_setl = 0.
            !if (do_sedimentation) then
            !    do n=1,npop
            !    if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
            !        .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk") )  then !exclude ext, akk, oc1 and bc1
            !        call rt_dt_from_sedimentation(rt, n, pfull, dt, t, pwt, is, ie, js, je,  rt_dt_setl, time_next, rh, kbot)
            !    endif
            !    enddo
            !endif

            !!if (mpp_root_pe().eq.mpp_pe()) then
            !if (any( (rt_dt_setl >huge(0.0)) .or. (rt_dt_setl < -huge(0.0)) )) then
            !    call error_mesg('matrix_run check 1', 'infinity exist in rt_setl', fatal)
            !endif
            !!endif

            !------------------------------------------------------------------------------------
            !               step *: h2so4 condensational growth, must prior new particle formation
            !(1) calculate the condensation coefficient for all population over all grids
            !           kci_coef_pop 4-d dimensions: (npop, is, ij, ik)
            !           kci_coef_aeq1 4-d dimensions: (npop, is, ij, ik)
            !-------------------------------------------------------------------------------------
            call mpp_clock_begin(condgrow_clock)
            call set_matrix_pop_kci(pfull,t,kci_coef_pop, kci_aeq1_pop,is,ie,js,je) !calculate the condensation coefficient for all population
            !---------------------------------------------------------------------------------
            !               step *: condensational growth w/wo new particle formation of akk
            !               calculate new particle formation rate and update sources of dndt/dmdt 
            !xlxlxlxlimportant: h2so4 need to be updated by the loss of h2so4 condensational loss
            !-----------------------------------------------------------------------------------
            call mpp_clock_begin(npf_clock)
            call set_matrix_condense_npf(i_akk, xnh3, fland, pfull,rh,t, dt, r(:,:,:,nh2so4), pwt,zhalf, kci_coef_pop, kci_aeq1_pop, &
                dndt_npf, dmdt_h2so4_npf, ih2so4_path, dmdt_h2so4_tot_cond_npf, kc, xh2so4_nucl, &
                is, ie, js, je) !the i_msulf term already updated for each population
            call mpp_clock_end(npf_clock)
            !update sources
            !assumption: h2so4(g, 98) -> aerosol(p, 96), h2so4 condensed on aerosol surface would become sulfate with molecular
            !weight 96, scale condensational h2so4 mass -> aerosol sulfate mass, unit: ug/m3/s
            pq_growth = (dmdt_h2so4_tot_cond_npf - dmdt_h2so4_npf)*matrix_molecular_weight(i_mw_so4)/matrix_molecular_weight(i_mw_h2so4)
            !scale new particle formation h2so4 mass -> aerosol sulfate mass, unit: ug/m3/s
            dmdt_npf = dmdt_h2so4_npf*matrix_molecular_weight(i_mw_so4)/matrix_molecular_weight(i_mw_h2so4)
            !-----------------------------------------------------------------------------------
            !                send emis data before source get changed with different processes
            !-----------------------------------------------------------------------------------
            do n=1, ntrace
            if (matrix_all_tracer(n)%id_tracer_emis > 0) then
                emis_tmp = 0.
                if (matrix_all_tracer(n)%type .eq. 'mass') then
                    do nct = 1, kd !ug/m3/s -> kg/m2/s
                    emis_tmp = emis_tmp + matrix_all_tracer(n)%source(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))*1e-9
                    enddo
                elseif(matrix_all_tracer(n)%type .eq. 'number') then
                    do nct = 1, kd !#/m3/s -> mol/m2/s
                    emis_tmp = emis_tmp + &
                        matrix_all_tracer(n)%source(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))/avogno
                    enddo
                endif
                used = send_data(matrix_all_tracer(n)%id_tracer_emis, emis_tmp, time_next, &
                    is_in=is,js_in=js)
                used = send_data(matrix_all_tracer(n)%id_tracer_emis_source, matrix_all_tracer(n)%source(is:ie,js:je,:), time_next, &
                    is_in=is,js_in=js, ks_in=1)
            endif
            enddo

            !-----------------------------------------------------------------------------------------------------------
            !               update matrix sources due to new particle formation and condensational growth
            !-----------------------------------------------------------------------------------------------------------            
            do n=1,npop
            cond_pop (n,:,:,:) = 0.0
            if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
                if ((i_akk > 0) .and. (n .eq. i_akk)) then
                    n_akk_emis_source = matrix_all_tracer(matrix_all_pop(i_akk)%i_n)%source(is:ie,js:je,:)
                    m_akk_emis_source = matrix_all_tracer(matrix_all_pop(i_akk)%i_msulf)%source(is:ie,js:je,:) !record ini source from emission
                    matrix_all_tracer(matrix_all_pop(i_akk)%i_n)%source(is:ie,js:je,:) = n_akk_emis_source+dndt_npf
                    matrix_all_tracer(matrix_all_pop(i_akk)%i_msulf)%source(is:ie,js:je,:) = m_akk_emis_source+dmdt_npf
                endif
                number_pop = matrix_all_tracer(matrix_all_pop(n)%i_n)%value_in_matrix(is:ie,js:je,:) !number concentration of a population
                do i=1,it
                do j=1,jt
                do k=1,kt
                if (kc(i,j,k) > 0) then
                    cond_pop (n,i,j,k) = (kci_coef_pop(n,i,j,k)*number_pop(i,j,k)/kc(i,j,k) ) * pq_growth(i,j,k)
                    matrix_all_tracer(matrix_all_pop(n)%i_msulf)%source(i+is-1,j+js-1,k) = &
                        matrix_all_tracer(matrix_all_pop(n)%i_msulf)%source(i+is-1,j+js-1,k) + cond_pop (n,i,j,k) 
                endif

                enddo
                enddo
                enddo

            endif                         
            enddo
            !update h2so4 loss
            dm_h2so4 = dm_h2so4 + dmdt_h2so4_tot_cond_npf * dt !in matrix unit ug/m3
            call mpp_clock_end(condgrow_clock)

            !-----------------------------------------------------------------------------------
            !                send condensation data before source get changed with different processes
            !-----------------------------------------------------------------------------------
            do n=1, npop
            if (matrix_all_pop(n)%id_pop_condens > 0) then
                cond_tmp = 0.
                do nct = 1, kd !ug/m3/s -> kg/m2/s
                cond_tmp = cond_tmp + cond_pop(n,:,:,nct)* (zhalf(:,:,nct)-zhalf(:,:,nct+1))*1e-9
                enddo
                used = send_data(matrix_all_pop(n)%id_pop_condens, cond_tmp, time_next, &
                    is_in=is,js_in=js)

                used = send_data(matrix_all_pop(n)%id_pop_condens_source, &
                    matrix_all_tracer(matrix_all_pop(n)%i_msulf)%source(is:ie, js:je, :), time_next, is_in=is,js_in=js, ks_in=1)
            endif
            enddo

            !-------------------------------------------------------------------------------
            !               step *: partition aqso4 to matrix: add to matrix sources
            !   currently, don't calculate activation, just distribute aqso4 based on numbers
            !---------------------------------------------------------------------------------
            if (do_aqso4) then
                call mpp_clock_begin(aqso4_clock)
                !get total number of all pops (exclude ext and akk)
                total_aero_num = 0.
                do n=1,npop
                if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
                    .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk") &
                    .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "oc1") &
                    .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "bc1") ) then !exclude ext, akk, oc1 and bc1
                    total_aero_num = total_aero_num + matrix_all_tracer(matrix_all_pop(n)%i_n)%value_in_matrix(is:ie,js:je,:) ![#/m3]
                endif
                enddo
                total_aero_num = max(total_aero_num, 1.0e-17)
                !add aqso4 source
                do n=1,npop
                if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
                    .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk") &
                    .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "oc1") &
                    .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "bc1") ) then !exclude akk, oc1 and bc1
                    matrix_all_tracer(matrix_all_pop(n)%i_msulf)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(n)%i_msulf)%source(is:ie,js:je,:) + &
                        aqso4(is:ie,js:je,:)*matrix_all_tracer(matrix_all_pop(n)%i_n)%value_in_matrix(is:ie,js:je,:)/total_aero_num ![ug/m3/s]

                    used = send_data(matrix_all_pop(n)%id_pop_aqso4_partition, &
                        aqso4(is:ie,js:je,:)*matrix_all_tracer(matrix_all_pop(n)%i_n)%value_in_matrix(is:ie,js:je,:)/total_aero_num, time_next, &
                        is_in=is,js_in=js, ks_in=1)
                endif
                enddo
                call mpp_clock_end(aqso4_clock)
            endif


            !updated coagulation: do mass solver and number solver inside do_coag, make sure conservation of number and mass
            if (do_coag) then
                call mpp_clock_begin(coag_clock)
                !     !-----------------------------------------------------------------------------
                !     !               step *: coagulation
                !     ! ri/li: (npop-1, size(r,1),size(r,2),size(r,3)), production/loss of number in mode i [#/m3/s]
                !     ! rim/lim: (npop-1,nmspcs,size(r,1),size(r,2),size(r,3)), production/loss of species qq mass in mode i [ug/m3/s]
                !     ! bi: intermodal number loss coefficient, li = bi*ni + 1/2*kbar0_ii*ni^2 is the total loss
                !     ! fi: intermodal mass loss coefficient, lim = fi*mjq
                !     !------------------------------------------------------------------------------
                !steps: (1) calculate total loss of mass: m(t)-m(t=0); assign the loss tendency to source;  
                !       (2) attribute the mass loss to the production of mass in each pop; assign it to the source
                call matrix_coag_cap_n(pfull, t, is,ie,js,je, time, time_next, dt, ri, li, rim, lim, lon, lat)
                call mpp_clock_end(coag_clock)
            endif

            !-----------------------------------------------------------------------------------------------------------
            !               update matrix sources due to coagulation and send the coagulation source
            !-----------------------------------------------------------------------------------------------------------
            do n=1,npop
            if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
                used = send_data(matrix_all_pop(n)%id_pop_coag_pn, ri(n,:,:,:), time_next, &
                    is_in=is,js_in=js, ks_in=1)
                used = send_data(matrix_all_pop(n)%id_pop_coag_ln, li(n,:,:,:), time_next, &
                    is_in=is,js_in=js, ks_in=1)
                if (matrix_all_pop(n)%i_msulf > 0) then !mass source [ug/m3/s]
                    !call module_mass_solver(n, 1, rim(n,1,:,:,:), fi(n,:,:,:), is, ie, js, je, dt)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_pmsulf, rim(n,1,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_lmsulf, lim(n,1,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                endif
                if (matrix_all_pop(n)%i_mdust > 0) then
                    !call module_mass_solver(n, 4, rim(n,4,:,:,:), fi(n,:,:,:), is, ie, js, je, dt)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_pmdust, rim(n,4,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_lmdust, lim(n,4,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                endif
                if (matrix_all_pop(n)%i_mseas > 0) then
                    !call module_mass_solver(n, 5, rim(n,5,:,:,:), fi(n,:,:,:), is, ie, js, je, dt)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_pmseas, rim(n,5,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_lmseas, lim(n,5,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                endif        
                if (matrix_all_pop(n)%i_mocar > 0) then
                    !call module_mass_solver(n, 3, rim(n,3,:,:,:), fi(n,:,:,:), is, ie, js, je, dt)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_pmocar, rim(n,3,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_lmocar, lim(n,3,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                endif        
                if (matrix_all_pop(n)%i_mbcar > 0) then
                    !call module_mass_solver(n, 2, rim(n,2,:,:,:), fi(n,:,:,:), is, ie, js, je, dt)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_pmbcar, rim(n,2,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                    used = send_data(matrix_all_pop(n)%id_pop_coag_lmbcar, lim(n,2,:,:,:), time_next, &
                        is_in=is,js_in=js, ks_in=1)
                endif
            endif
            enddo

            !---------------------debugging check: if nan exist and do value updates-----------------------
            do n=1,ntrace
            if (matrix_all_tracer(n)%is_active) then
                if (any(ieee_is_nan(matrix_all_tracer(n)%source(is:ie, js:je, :)))) then
                    if(mpp_pe() == mpp_root_pe()) then
                        write (*, *), "loc 0: nan found in matrix_source", n, matrix_all_tracer(n)%pop, matrix_all_tracer(n)%spec, ")"
                    endif
                endif

                if (any(ieee_is_nan(matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :)))) then 
                    if(mpp_pe() == mpp_root_pe()) then
                        write (*, *), "loc 1: nan found in matrix_all_tracer(", n, matrix_all_tracer(n)%pop, matrix_all_tracer(n)%spec, ")"
                    endif
                endif
            endif
            enddo


            !-----------------------------------------------------------------------------
            !               step *: matrix values get updated
            !-----------------------------------------------------------------------------
            do n=1,ntrace
            if (matrix_all_tracer(n)%is_active) then
                if (matrix_all_tracer(n)%type .eq. 'mass') then
                    matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :) = &
                        matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :) + &
                        matrix_all_tracer(n)%source(is:ie, js:je, :)*dt     
                    !matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:) &
                    !    = max(matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:), 1.e-32) !value can't < 0
                elseif (matrix_all_tracer(n)%type .eq. 'number') then
                    matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :) = &
                        matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :) + &
                        matrix_all_tracer(n)%source(is:ie, js:je, :)*dt

                    !matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:) &
                    !    = max(matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:), 1.e-17) !value can't < 0

                endif
            endif
            enddo

            !---------------------debugging check: if negative exist-----------------------
            do n=1,ntrace
            if (matrix_all_tracer(n)%is_active) then
                if (any(matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :) < 0)) then
                    if(mpp_pe() == mpp_root_pe()) then
                        write(*,*), "loc 1-: negative found before intermodal(", n, matrix_all_tracer(n)%pop, matrix_all_tracer(n)%spec, ")"
                    endif
                endif
            endif
            enddo


            !-----------------------------------------------------------------------------
            !               step *: intermodal transfer + update matrix value
            !   pseudo update matrix_tracer values -> intermodal subroutine: update value again
            !---------------------------------------------------------------------------------------------------- 
            !call matrix_dry_diameter(pfull,is,ie,js,je) !assign dry pop properties
            !call matrix_wet_diameter(rh, t,is,ie,js,je) !assign wet pop properties, t is temperature

            if (do_intermodal_transfer) then
                call mpp_clock_begin(dmodal_clock)
                !note: matrix value got updated in this matrix_intermodal_transfer subroutine
                call matrix_intermodal_transfer(intermodal_configuration, rt, is,ie,js,je, time_next)
                !do n=1,ntrace
                !if (matrix_all_tracer(n)%is_active) then
                !    if (matrix_all_tracer(n)%type .eq. 'mass') then
                !        matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:) &
                !            = max(matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:), 1.e-32) !value can't < 0
                !    elseif (matrix_all_tracer(n)%type .eq. 'number') then
                !        matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:) &
                !            = max(matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:), 1.e-17) !value can't < 0

                !    endif
                !endif
                !enddo
                call mpp_clock_end(dmodal_clock)
            endif

            !---------------------debugging check: if nan exist-----------------------
            do n=1,ntrace
            if (matrix_all_tracer(n)%is_active) then
                if (any(ieee_is_nan(matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :)))) then
                    if(mpp_pe() == mpp_root_pe()) then
                        write(*,*), "loc 2: nan found in matrix_all_tracer(", n, matrix_all_tracer(n)%pop, matrix_all_tracer(n)%spec, ")"
                    endif
                endif
                if (any(matrix_all_tracer(n)%value_in_matrix(is:ie, js:je, :) < 0)) then
                    if(mpp_pe() == mpp_root_pe()) then
                        write(*,*), "loc 2-: negative found after intermodal(", n, matrix_all_tracer(n)%pop, matrix_all_tracer(n)%spec, ")"
                    endif
                endif
            endif
            enddo
           ! !------------------------------------------------------------------------ 
           ! call matrix_dry_diameter(pfull,is,ie,js,je) !assign dry pop properties
           ! call matrix_wet_diameter(rh, t,is,ie,js,je) !assign wet pop properties, t is temperature


            !-----------------------------------------------------------------------------
            !               step *: sedimentation settling by gravity
            !  ->  calculate gravitational settling and update value again in the subroutine
            !----------------------------------------------------------------------------------------------------
            do n=1,ntrace
            if (matrix_all_tracer(n)%is_active) then
                !rt is the tracer concentration in gfdl unit
                call update_rt_from_matrix(n, rt(:,:,:,n), pwt, zhalf,is,ie,js,je)
            endif
            enddo

            if (mpp_root_pe().eq.mpp_pe()) then
                if (any( (rt >huge(0.0)) .or. (rt < -huge(0.0)) )) then
                    call error_mesg('matrix_run check 1', 'infinity exist in rt', fatal)
                endif
            endif

            ! rt_dt_setl = 0.
            !if (do_sedimentation) then
            !    do n=1,npop
            !    if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
            !        .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk") )  then !exclude ext, akk, oc1 and bc1
            !        call rt_dt_from_sedimentation(rt, n, pfull, dt, t, pwt, is, ie, js, je,  rt_dt_setl, time_next, rh, kbot)
            !    endif
            !    enddo
            !endif

            !!if (mpp_root_pe().eq.mpp_pe()) then
            !if (any( (rt_dt_setl >huge(0.0)) .or. (rt_dt_setl < -huge(0.0)) )) then
            !    call error_mesg('matrix_run check 1', 'infinity exist in rt_setl', fatal)
            !endif
            !!endif

            !!DEBUG: check rt_dt_setl relevant to rt/dt, as don't want to have negative values
            !!if (mpp_root_pe().eq.mpp_pe()) then
            !do n=1,ntrace
            !if (matrix_all_tracer(n)%is_active) then
            !    if (any(rt(:,:,:,n)+rt_dt_setl(:,:,:,n)*dt < 0)) then
            !        write(*,*) 'negative value exist after settling, minval is ', matrix_all_tracer(n)%name, &
            !            minval(rt(:,:,:,n)+rt_dt_setl(:,:,:,n)*dt)
            !    endif
            !endif
            !enddo
            !!endif


            !---------------------------------------------------------------------------------
            !               final step: pass matrix_tracer values and h2so4 concentration to gfdl tendency
            !---------------------------------------------------------------------------------
            do n=1,ntrace
            if (matrix_all_tracer(n)%is_active) then
                !call update_rt_from_matrix(n, rt(:,:,:,n), pwt, zhalf,is,ie,js,je) !convert tendency from matrix to gfdl unit, ug/m3-> mmr or vmr
                !rdt_matrix(:,:,:,n) = (rt(:,:,:,n)-r(:,:,:,n))/dt + rt_dt_setl(:,:,:,n)
                rdt_matrix(:,:,:,n) = (rt(:,:,:,n)-r(:,:,:,n))/dt !here the rt already included the change from sedimentation
            elseif (n .eq. nh2so4) then
                !do unit from matrix unit to gfdl unit: ug/m3 -> vmr
                mw = matrix_molecular_weight(i_mw_h2so4)
                dm_h2so4 = dm_h2so4/1.0e9/(pwt(:,:,:) * mw / wtmair / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))) 
                !dm_h2so4 can't exceed previous h2so4 concentration
                do i=1,it
                do j=1,jt
                do k=1,kt
                dm_h2so4(i, j, k) = min(dm_h2so4(i, j, k), r(i, j, k, nh2so4))
                dm_h2so4(i, j, k) = max(0.0, dm_h2so4(i, j, k))
                enddo
                enddo
                enddo
                rdt_matrix(:,:,:,n) = -dm_h2so4/dt !in gfdl unit
            endif
            enddo

            if (mpp_root_pe().eq.mpp_pe()) then
                if (any( (rdt_matrix >huge(0.0)) .or. (rdt_matrix < -huge(0.0)) )) then
                    call error_mesg('matrix_run check 1', 'infinity exist in rdt_matrix', fatal)
                endif
                if (any( (r>huge(0.0)) .or. (r < -huge(0.0)) )) then
                    call error_mesg('matrix_run check 1', 'infinity exist in r', fatal)
                endif
                if (any( (rt>huge(0.0)) .or. (rt < -huge(0.0)) )) then
                    call error_mesg('matrix_run check 1', 'infinity exist in rt', fatal)
                endif
            endif
            !---------------------------------------------------------------------------------------------------------------
            !---------------------------------budget analysis test: send out diagnostic of tendency, source, and values
            !-----------------------------------------------------------------------------------------------------------------
            do n=1, ntrace
            if (matrix_all_tracer(n)%id_tendency_gfdl > 0) then
                if (matrix_all_tracer(n)%type .eq. 'mass') then
                    used = send_data(matrix_all_tracer(n)%id_tendency_gfdl, rdt_matrix(:,:,:,n), time_next, &
                        is_in=is,js_in=js,ks_in=1)
                elseif(matrix_all_tracer(n)%type .eq. 'number') then
                    used = send_data(matrix_all_tracer(n)%id_tendency_gfdl, rdt_matrix(:,:,:,n), time_next, &
                        is_in=is,js_in=js,ks_in=1)
                endif
            endif
            enddo
            !-----------------------------------------------------------end-------------------------------------------------------

            !--------------------------------------------------------------------------------------------------------------------
            ! register and send data
            !--------------------------------------------------------------------------------------------------------------------
            !register and send data of rh
            used = send_data (id_rh, rh*100, time_next, &
                is_in=is,js_in=js, ks_in = 1) !in unit %
            used = send_data (id_cond_sink, pq_growth*dt, time_next, &
                is_in=is,js_in=js, ks_in = 1) !in ug/m3
            used = send_data (id_kc, kc, time_next, &
                is_in=is,js_in=js, ks_in = 1) 
            used = send_data (id_dmdt_h2so4_tot_cond_npf, dmdt_h2so4_tot_cond_npf, time_next, &
                is_in=is,js_in=js, ks_in = 1) 
            used = send_data (id_dndt_npf, dndt_npf, time_next, &
                is_in=is,js_in=js, ks_in = 1) 
            used = send_data (id_dmdt_h2so4_npf, dmdt_h2so4_npf, time_next, &
                is_in=is,js_in=js, ks_in = 1) 
            used = send_data (id_pwt, pwt, time_next, &
                is_in=is,js_in=js, ks_in = 1)
            used = send_data (id_zhalf, zhalf(:,:,1:kd)-zhalf(:,:,2:kd+1), time_next, &
                is_in=is,js_in=js, ks_in = 1)


            ! !send data of d_wet, d_dry for test: ps. dg_dry, dg_wet in matrix_pop unit is m, only for send_data change to um
            ! do n = 1, npop
            ! if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
            !     if (matrix_all_pop(n)%id_dg_dry > 0) then
            !         used = send_data (matrix_all_pop(n)%id_dg_dry, matrix_all_pop(n)%dg_dry(is:ie,js:je,:)*1e6, time_next, &
            !             is_in=is,js_in=js, ks_in = 1) !in unit µm
            !     endif
            !     if (matrix_all_pop(n)%id_dg_wet > 0) then
            !         used = send_data (matrix_all_pop(n)%id_dg_wet, matrix_all_pop(n)%dg_wet(is:ie,js:je,:)*1e6, time_next, &
            !             is_in=is,js_in=js, ks_in = 1) !in unit µm
            !     endif

            ! endif
            ! end do

            !-----------------------------------------------------------------------------------
            !                send col data for matrix tracers
            !-----------------------------------------------------------------------------------
            do n=1, ntrace
            if (matrix_all_tracer(n)%id_tracer_col > 0) then
                col_tmp = 0.
                if (matrix_all_tracer(n)%type .eq. 'mass') then
                    do nct = 1, kd !ug/m3/s -> kg/m2/s
                    col_tmp = col_tmp + matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))*1e-9
                    enddo
                elseif(matrix_all_tracer(n)%type .eq. 'number') then
                    do nct = 1, kd !#/m3/s -> mol/m2/s
                    col_tmp = col_tmp + &
                        matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))/avogno
                    enddo
                endif
                used = send_data(matrix_all_tracer(n)%id_tracer_col, col_tmp, time_next, &
                    is_in=is,js_in=js)
            endif
            enddo

        endif

        do n=1, ntrace
        if (matrix_all_tracer(n)%is_active) then !active matrix tracer in current configuration
            matrix_all_tracer(n)%source(is:ie,js:je,:) = 0.
            !matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:) = 1.e-32
        endif
        enddo

        call mpp_clock_end(matrix_run_clock)

    end subroutine matrix_run



    !----------------------------------------------------------------
    !
    !                           subroutine matrix_init(r)
    !
    !-------------------initialization-------------------
    !   1. determine configuration and number of populations: npop, allocate(matrix_pop(npop))
    !   2. determine number of tracers and get parameters
    !----------------------------------------------------------------
    subroutine matrix_init(phalf, axes, time)
        real, intent(in), dimension(:,:,:) :: phalf    
        integer, intent(in) :: axes(4)              ! diagnostic axes
        type(time_type),  intent(in) :: time       ! model time
        character(len=512) :: text_in_scheme, control
        character*32 :: name_tr_wetdep,units_tr_wetdep
        integer :: i, ip, klq, k, l, qq,  ikl, ipop, kpop, lpop, n, kd
        integer :: io, logunit, verbose, ierr
        integer :: nt, ntt
        logical :: flag
        verbose = 3
        !----------------------------------------------------------------
        !	------------------- read nml -------------------
        !----------------------------------------------------------------
        if (matrix_module_init) return

        if(file_exist('input.nml')) then
            !            #ifdef internal_file_nml
            read (input_nml_file, nml = matrix_nml, iostat = io)
            ierr = check_nml_error(io,'matrix_nml')
            ! #else
            ! unit = open_namelist_file('input.nml')
            ! ierr=1; do while (ierr /= 0)
            ! read(unit, nml = matrix_nml, iostat = io, end = 10)
            ! ierr = check_nml_error (io, 'matrix_nml')
            ! end do
            ! 10  call close_file(unit)
            ! #endif
        endif

        logunit = stdlog()
        if(mpp_pe() == mpp_root_pe()) then
            write(logunit, nml=matrix_nml)
            verbose = verbose + 1
        end if

        if ( do_matrix) then
            !save ids for important tracers
            nsphum = get_tracer_index(model_atmos,'sphum') 
            nh2so4 = get_tracer_index(model_atmos,'simpleh2so4')
            if (nh2so4.le.0) then
                nh2so4 = get_tracer_index(model_atmos,'h2so4')
            end if

            if (nh2so4.le.0) &
                call error_mesg ('matrix_gfdl','h2so4 needs to be defined!!!!', fatal)
            call set_config_pop(matrix_configuration)
            !----------------------------------------------------------------
            !	--- assign configuration population information ----
            ! akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
            ! if population exist, assign index >= 1, otherwise -1
            !----------------------------------------------------------------
            ! if (matrix_configuration .eq. "debug") then
            !     !debug version: 
            !     !pop: acc <- (n,m), ext <- (m_h2o, m_h2so4)
            !     !npop = 2 !number of poplation in the configuration
            !     !assign population index under current configuration
            !     !integer :: i_akk = 1, i_acc = 2, i_dd1 = 3, i_dd2 = 4 ! index of population in matrix
            !     !integer :: i_ssa = 5, i_ssc = 6, i_oc1 = 7, i_oc2 = 8 ! if exit, index >= 1, otherwise = -1
            !     !integer :: i_bc1 = 9, i_bc2 = 10, i_mxa = 11, i_mxc = 12, i_ext = 13
            !     i_acc = 2
            !     i_ext = 13
            !     pop_active = 1 !neglect ext
            !     allocate(i_pop(pop_active))
            !     i_pop(:)=(/i_acc/)
            ! elseif (matrix_configuration .eq. "debug_npf") then
            !     i_akk = 1
            !     pop_active = 1
            !     allocate(i_pop(pop_active))
            !     i_pop(:)=(/i_akk/)
            ! else
            !     call error_mesg('get matrix configuration','configuration not defined '//trim(matrix_configuration), fatal)
            !     !call error_mesg ('tracer_driver', 'mw needs to be defined for tracer: '//trim(tracer_name), fatal)
            ! end if

            call get_number_tracers(model_atmos, num_tracers = ntrace)

            allocate(matrix_all_tracer(ntrace))
            allocate(matrix_all_pop(npop))
            do n=1,npop
            matrix_all_pop(n)%nb_tracer_pop = 0
            !                matrix_all_pop(n)%tracer_index(:) = -1
            matrix_all_pop(n)%tracer_index = -1
            enddo

            !read in deposition fixed const from field table
            !call query_caerosol_wetdep_param(1, cdust_wetdep_param)
            !call query_caerosol_wetdep_param(2, cssalt_wetdep_param)

            !xl debug1-------------------------------------------
            if (ntrace > 0) then
                do n = 1, ntrace
                flag = query_method ('matrix_parameter',model_atmos,n, &
                    text_in_scheme,control)
                if (flag) then
                    !important: get_matrix_tracer_param only 
                    !(1) update has_emission for mass species
                    !(2) allocate source array for mass species with lognormal, dens etc parameters in field table
                    !(3) important: number tracer: has_emission & source_array allocation need to be set somewhere

                    call get_matrix_tracer_param(n,text_in_scheme, control, matrix_all_tracer(n), matrix_all_pop)!, r)
                    if (matrix_all_tracer(n)%is_active) then !active matrix tracer in current configuration
                        allocate(matrix_all_tracer(n)%source(size(phalf,1),size(phalf,2),size(phalf,3)-1))!mass tracers has sources !!!!xl!!!!pop_index the number need
                        matrix_all_tracer(n)%source(:,:,:) = 0.
                        allocate(matrix_all_tracer(n)%value_in_matrix(size(phalf,1),size(phalf,2),size(phalf,3)-1))
                        matrix_all_tracer(n)%value_in_matrix = 1.e-32
                    endif
                end if
                end do

            else
                call error_mesg('matrix_gfdl/matrix_init', 'no atmos_tracer found ', fatal)
            endif

            !set number tracer: has_emission & allocate source_array allocation
            do n=1,npop
            if (matrix_all_pop(n)%nb_tracer_pop > 0) then
                do nt = 1, matrix_all_pop(n)%nb_tracer_pop
                ntt = matrix_all_pop(n)%tracer_index(nt)
                if (matrix_all_tracer(ntt)%has_emission) then
                    matrix_all_tracer(matrix_all_pop(n)%i_n)%has_emission = .true.
                    matrix_all_pop(n)%def_dens = matrix_all_tracer(ntt)%dens
                    matrix_all_pop(n)%def_kappa = matrix_all_tracer(ntt)%kappa
                    matrix_all_pop(n)%def_dg = matrix_all_tracer(ntt)%dgn
                endif 
                enddo 
            endif    
            enddo

            allocate(debug_num_domain(npop-1, size(phalf,1),size(phalf,2),size(phalf,3)-1))
            allocate(debug_mass_domain(npop-1, 5, size(phalf,1),size(phalf,2),size(phalf,3)-1))
            debug_num_domain = 0.
            debug_mass_domain = 0.

            do n=1,npop
            allocate(matrix_all_pop(n)%dg_dry(size(phalf,1),size(phalf,2),size(phalf,3)-1))
            allocate(matrix_all_pop(n)%dg_wet(size(phalf,1),size(phalf,2),size(phalf,3)-1))
            allocate(matrix_all_pop(n)%mass_dry(size(phalf,1),size(phalf,2),size(phalf,3)-1))
            allocate(matrix_all_pop(n)%vol_dry(size(phalf,1),size(phalf,2),size(phalf,3)-1))
            allocate(matrix_all_pop(n)%kappa_pop(size(phalf,1),size(phalf,2),size(phalf,3)-1))
            allocate(matrix_all_pop(n)%dens_wet(size(phalf,1),size(phalf,2),size(phalf,3)-1))
            allocate(matrix_all_pop(n)%dens_dry(size(phalf,1),size(phalf,2),size(phalf,3)-1))
            matrix_all_pop(n)%dg_dry = 1.0e-10
            matrix_all_pop(n)%dg_wet = 1.0e-10
            matrix_all_pop(n)%mass_dry = 1.0e-26 !unit: ug/m3 (1e-32 mmr to ug/m3)
            matrix_all_pop(n)%vol_dry = 1.0e-38 !unit: m3_aerosol / m3, must be real(8), because single precision is only 1e-38
            matrix_all_pop(n)%kappa_pop = 1.0e-40
            matrix_all_pop(n)%dens_wet = 1000. !initialize the density
            matrix_all_pop(n)%dens_dry = 1000. !initilize the density
            end do

            allocate(p_h2so4_rate(size(phalf,1),size(phalf,2),size(phalf,3)-1)) !h2so4 production
            allocate(aqso4(size(phalf,1),size(phalf,2),size(phalf,3)-1)) !aqueous phase so4 production
            p_h2so4_rate = 0.0
            aqso4 = 0.0

            !assign rh_deliquescence and rh_crystalization to population
            do n = 1, npop
            if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
                select case (lowercase(trim(matrix_all_pop(n)%name)))
                case ('akk')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1770
                    matrix_all_pop(n)%sigma = 1.6
                case ('acc')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1770
                case ('dd1')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 2600
                case ('dd2')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 2600
                case ('ssa')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 2200
                case ('ssc')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 2200
                    matrix_all_pop(n)%sigma = 2.0
                case ('oc1')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                case ('oc2')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                case ('bc1')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                case ('bc2')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                case ('mxa')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                    if (i_acc > 0) then
                        matrix_all_pop(n)%def_dens = matrix_all_pop(i_acc)%def_dens
                        matrix_all_pop(n)%def_kappa = matrix_all_pop(i_acc)%def_kappa
                        matrix_all_pop(n)%def_dg = matrix_all_pop(i_acc)%def_dg
                    else
                        matrix_all_pop(n)%def_dens = 1770
                        matrix_all_pop(n)%def_kappa = 0.507
                        matrix_all_pop(n)%def_dg = 0.068e-6
                    endif
                case ('mxc')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                    matrix_all_pop(n)%sigma = 2.0
                    matrix_all_pop(n)%def_dens = 2200
                    matrix_all_pop(n)%def_kappa = 0.068
                    matrix_all_pop(n)%def_dg = 1e-6
                case ('ext')
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                case default
                    matrix_all_pop(n)%rh_deliquescence = 0.8
                    matrix_all_pop(n)%rh_crystallization = 0.35
                    matrix_all_pop(n)%dens_dry = 1000
                end select
            endif
            end do 
            !set up coagulation configuration and arrays
            if (do_coag) then 
                call setup_coag_tensors(coag_configuration)
                call setup_kij_diameters
                call setup_kij_tables
            endif

            kd = size(phalf,3)-1
            k_grid_size = kd
            !sanility check
            do n=1,npop
            if (matrix_all_pop(n)%nb_tracer_pop > 0) then !if this population exist in current configuration
                do nt = 1,matrix_all_pop(n)%nb_tracer_pop
                if (matrix_all_pop(n)%has_emission(nt)) then
                    ntt = matrix_all_pop(n)%tracer_index(nt)
                    if ((matrix_all_tracer(ntt)%sigma > 0) .and. (matrix_all_tracer(ntt)%distribution_index == dist_lognormal) &
                        .and. (matrix_all_tracer(ntt)%dgn > 0) .and. (matrix_all_tracer(ntt)%dens > 0) ) then !suggest to move this check to initialization
                    else
                        write(*,*) 'n=', n, 'pop name', matrix_all_pop(n)%name, nt, ntt, matrix_all_tracer(ntt)%name
                        call error_mesg ('matrix_gfdl','emission species not properly defined', fatal)
                    endif
                endif
                enddo
            endif
            enddo


            do n = 1, npop
            if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
                matrix_all_pop(n)%id_dg_dry = register_diag_field ( module_name,     &
                    trim(matrix_all_pop(n)%name)//'_dg_dry', axes(1:3),time,  &
                    trim(matrix_all_pop(n)%name)//'_dg_dry', 'µm',       &
                    missing_value=-999.  )

                matrix_all_pop(n)%id_dg_wet = register_diag_field ( module_name,     &
                    trim(matrix_all_pop(n)%name)//'_dg_wet', axes(1:3),time,  &
                    trim(matrix_all_pop(n)%name)//'_dg_wet', 'µm',       &
                    missing_value=-999.  )

            endif
            end do

            id_rh = register_diag_field ( module_name,'rh_matrix', axes(1:3),time, 'rh_matrix', '%',       &
                missing_value=-999.  ) 
            id_cond_sink = register_diag_field ( module_name,'cond_sink', axes(1:3),time, 'cond_sink', 'ug/m3/(30 min)',       &
                missing_value=-999.  )
            id_kc = register_diag_field ( module_name,'kc', axes(1:3),time, 'kc', '1/s',       &
                missing_value=-999.  )
            id_dmdt_h2so4_tot_cond_npf = register_diag_field ( module_name,'dmdt_h2so4_tot_cond_npf', axes(1:3),time, 'dmdt_h2so4_tot_cond_npf', 'ug/m3/s',    &
                missing_value=-999.  )
            id_dndt_npf = register_diag_field ( module_name,'dndt_npf', axes(1:3),time, 'dndt_npf', '#/m3/s',    &
                missing_value=-999.  )
            id_dmdt_h2so4_npf = register_diag_field ( module_name,'dmdt_h2so4_npf', axes(1:3),time, 'dmdt_h2so4_npf', 'ug/m3/s',    &
                missing_value=-999.  )
            id_pwt = register_diag_field ( module_name,'pwt', axes(1:3),time, 'pwt', 'kg air/m2',    &
                missing_value=-999.  )
            id_zhalf = register_diag_field ( module_name,'zhalf', axes(1:3),time, 'zhalf', 'm',    &
                missing_value=-999.  )

            !id for coagulation debug
            id_mjq11 = register_diag_field ( module_name, 'mjq11', axes(1:3),time, 'mjq11', 'ug/particle', missing_value=-999.  )
            id_mjq21 = register_diag_field ( module_name, 'mjq21', axes(1:3),time, 'mjq21', 'ug/particle', missing_value=-999.  )
            id_kbar0_11 = register_diag_field(module_name, 'kbar0_11', axes(1:3), time, 'kbar0_11', 'm3/s', missing_value=-999.)
            id_kbar3_11 = register_diag_field(module_name, 'kbar3_11', axes(1:3), time, 'kbar3_11', 'm3/s', missing_value=-999.)
            id_kbar0_12 = register_diag_field(module_name, 'kbar0_12', axes(1:3), time, 'kbar0_12', 'm3/s', missing_value=-999.)
            id_kbar3_12 = register_diag_field(module_name, 'kbar3_12', axes(1:3), time, 'kbar3_12', 'm3/s', missing_value=-999.)
            id_kbar0_21 = register_diag_field(module_name, 'kbar0_21', axes(1:3), time, 'kbar0_21', 'm3/s', missing_value=-999.)
            id_kbar3_21 = register_diag_field(module_name, 'kbar3_21', axes(1:3), time, 'kbar3_21', 'm3/s', missing_value=-999.)
            id_kbar0_22 = register_diag_field(module_name, 'kbar0_22', axes(1:3), time, 'kbar0_22', 'm3/s', missing_value=-999.)
            id_kbar3_22 = register_diag_field(module_name, 'kbar3_22', axes(1:3), time, 'kbar3_22', 'm3/s', missing_value=-999.)

            id_ri_1 = register_diag_field(module_name, 'ri_1', axes(1:3), time, 'ri_1', '#/m3/s', missing_value=-999.)
            id_ri_2 = register_diag_field(module_name, 'ri_2', axes(1:3), time, 'ri_2', '#/m3/s', missing_value=-999.)
            id_li_1 = register_diag_field(module_name, 'li_1', axes(1:3), time, 'li_1', '#/m3/s', missing_value=-999.)
            id_li_2 = register_diag_field(module_name, 'li_2', axes(1:3), time, 'li_2', '#/m3/s', missing_value=-999.)
            id_bi_1 = register_diag_field(module_name, 'bi_1', axes(1:3), time, 'bi_1', '#/m3/s', missing_value=-999.)
            id_bi_2 = register_diag_field(module_name, 'bi_2', axes(1:3), time, 'bi_2', '#/m3/s', missing_value=-999.)
            id_fi_1 = register_diag_field(module_name, 'fi_1', axes(1:3), time, 'fi_1', '#/m3/s', missing_value=-999.)
            id_fi_2 = register_diag_field(module_name, 'fi_2', axes(1:3), time, 'fi_2', '#/m3/s', missing_value=-999.)

            id_rim_11 = register_diag_field(module_name, 'rim_11', axes(1:3), time, 'rim_11', 'ug/m3/s', missing_value=-999.)
            id_rim_21 = register_diag_field(module_name, 'rim_21', axes(1:3), time, 'rim_21', 'ug/m3/s', missing_value=-999.)
            id_lim_11 = register_diag_field(module_name, 'lim_11', axes(1:3), time, 'lim_11', 'ug/m3/s', missing_value=-999.)
            id_lim_21 = register_diag_field(module_name, 'lim_21', axes(1:3), time, 'lim_21', 'ug/m3/s', missing_value=-999.)

            id_half_knn_1 = register_diag_field(module_name, 'half_knn_1', axes(1:3), time, 'half_knn_1', '#/m3/s', missing_value=-999.)
            id_notself_knn_1 = register_diag_field(module_name, 'notself_knn_1', axes(1:3), time, 'notself_knn_1', '#/m3/s', missing_value=-999.)
            id_half_knn_2 = register_diag_field(module_name, 'half_knn_2', axes(1:3), time, 'half_knn_2', '#/m3/s', missing_value=-999.)
            id_notself_knn_2 = register_diag_field(module_name, 'notself_knn_2', axes(1:3), time, 'notself_knn_2', '#/m3/s', missing_value=-999.)


            !xl matrix tracer register
            id_so4_emis = register_diag_field ( module_name,'so4_emis', axes(1:2),time, 'so4_emis', 'kg/m2/s',    &
                missing_value=-999.  )
            do n=1, ntrace
            if (matrix_all_tracer(n)%is_active) then
                if (matrix_all_tracer(n)%type .eq. 'mass') then
                    matrix_all_tracer(n)%id_tracer_emis = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_emis', axes(1:2),time,  &
                        trim(matrix_all_tracer(n)%name)//'_emis', 'kg/m2/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_emis_source = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_emis_source', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_emis_source', 'ug/m3/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_col = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_col', axes(1:2),time,  &
                        trim(matrix_all_tracer(n)%name)//'_col', 'kg/m2',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_setl = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_setl', axes(1:2),time,  &
                        trim(matrix_all_tracer(n)%name)//'_setl', 'kg/m2/s',  &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_setl_3d = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_setl_3d', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_setl_3d', 'mmr/s',  &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_tsource = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_tsource', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_tsource', 'ug/m3/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_intermodal = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_intermodal', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_intermodal', 'ug/m3',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tendency_matrix = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_tendency_matrix', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_tendency_matrix', 'ug/m3/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tendency_gfdl = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_tendency_gfdl', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_tendency_gfdl', 'mmr/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_value_bf_max = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_value_bf_max', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_value_bf_max', 'ug/m3',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_value_af_max = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_value_af_max', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_value_af_max', 'ug/m3',       &
                        missing_value=-999.  )
                elseif (matrix_all_tracer(n)%type .eq. 'number') then
                    matrix_all_tracer(n)%id_tracer_emis = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_emis', axes(1:2),time,  &
                        trim(matrix_all_tracer(n)%name)//'_emis', 'mol/m2/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_emis_source = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_emis_source', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_emis_source', '#/m3/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_col = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_col', axes(1:2),time,  &
                        trim(matrix_all_tracer(n)%name)//'_col', 'mol/m2',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_setl = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_setl', axes(1:2),time,  &
                        trim(matrix_all_tracer(n)%name)//'_setl', 'mol/m2/s',  &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_setl_3d = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_setl_3d', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_setl_3d', 'mol/m3/s',  &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_intermodal = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_intermodal', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_intermodal', '#/m3',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tracer_tsource = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_tsource', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_tsource', '#/m3/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tendency_matrix = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_tendency_matrix', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_tendency_matrix', '#/m3/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_tendency_gfdl = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_tendency_gfdl', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_tendency_gfdl', '#/kg/s',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_value_bf_max = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_value_bf_max', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_value_bf_max', '#/m3',       &
                        missing_value=-999.  )
                    matrix_all_tracer(n)%id_value_af_max = register_diag_field ( module_name,     &
                        trim(matrix_all_tracer(n)%name)//'_value_af_max', axes(1:3),time,  &
                        trim(matrix_all_tracer(n)%name)//'_value_af_max', '#/m3',       &
                        missing_value=-999.  )
                endif

            endif
            if (n .eq. nh2so4) then
                id_h2so4_source = register_diag_field ( module_name,     &
                    'h2so4_source', axes(1:2),time,  &
                    'h2so4_source', 'kg/m2/s', missing_value=-999.  )
            endif
            enddo

            do n=1,npop
            if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
                matrix_all_pop(n)%id_pop_condens = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_condens', axes(1:2),time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_condens', 'kg/m2/s',       &
                    missing_value=-999.  )
                matrix_all_pop(n)%id_pop_coag_pn = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pn', axes(1:3),time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pn', '#/m3/s',       &
                    missing_value=-999.  )
                matrix_all_pop(n)%id_pop_coag_ln = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_ln', axes(1:3),time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_ln', '#/m3/s',       &
                    missing_value=-999.  )                
                matrix_all_pop(n)%id_pop_coag_pmsulf = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmsulf', axes(1:3),time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmsulf', 'ug/m3/s',       &
                    missing_value=-999.  )  
                matrix_all_pop(n)%id_pop_coag_lmsulf = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmsulf', axes(1:3),time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmsulf', 'ug/m3/s',       &
                    missing_value=-999.  )
                matrix_all_pop(n)%id_pop_coag_pmocar = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmocar', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmocar', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_coag_lmocar = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmocar', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmocar', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_coag_pmbcar = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmbcar', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmbcar', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_coag_lmbcar = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmbcar', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmbcar', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_coag_pmdust = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmdust', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmdust', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_coag_lmdust = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmdust', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmdust', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_coag_pmseas = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmseas', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pmseas', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_coag_lmseas = register_diag_field(module_name, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmseas', axes(1:3), time, &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lmseas', 'ug/m3/s', &
                    missing_value=-999.)

                matrix_all_pop(n)%id_pop_aqso4_partition = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_aqso4_partition', axes(1:3),time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_aqso4_partition', 'ug/m3/s',       &
                    missing_value=-999. )

            endif
            enddo

            id_aqso4_rate = register_diag_field ( module_name, 'tot_aqso4_rate', axes(1:3),time,    &
                'tot_aqso4_rate', 'ug/m3/s', missing_value=-999. )
        end if      
        !initialize clock id
        ini_clock = mpp_clock_id ('matrix: initialization', grain=clock_module)
        npf_clock = mpp_clock_id ('matrix: npf', grain=clock_module)
        condgrow_clock = mpp_clock_id ('matrix: condensation_grow', grain=clock_module)
        hygrow_clock = mpp_clock_id ('matrix: hygroscopic_growth', grain=clock_module)
        hygrow_sub_clock1 = mpp_clock_id ('matrix: hygro_subclock1', grain=clock_module)
        hygrow_sub_clock2 = mpp_clock_id ('matrix: hygro_subclock2', grain=clock_module)
        hygrow_sub_clock3 = mpp_clock_id ('matrix: hygro_subclock3', grain=clock_module)
        hygrow_sub_clock4 = mpp_clock_id ('matrix: hygro_subclock4', grain=clock_module)

        dry_clock = mpp_clock_id ('matrix: dry diameter', grain=clock_module)
        coag_clock = mpp_clock_id ('matrix: coagulation', grain=clock_module)
        coag_sub_clock1 = mpp_clock_id ('matrix: coag_subclock1', grain=clock_module)
        coag_sub_clock2 = mpp_clock_id ('matrix: coag_subclock2', grain=clock_module)
        coag_sub_clock3 = mpp_clock_id ('matrix: coag_subclock3', grain=clock_module)
        coag_sub_clock4 = mpp_clock_id ('matrix: coag_subclock4', grain=clock_module)
        coag_sub_clock5 = mpp_clock_id ('matrix: coag_subclock5', grain=clock_module)
        coag_sub_clock6 = mpp_clock_id ('matrix: coag_subclock6', grain=clock_module)
        coag_sub_clock7 = mpp_clock_id ('matrix: coag_subclock7', grain=clock_module)

        dmodal_clock = mpp_clock_id ('matrix: inter-modal transfer', grain=clock_module)
        aqso4_clock = mpp_clock_id ('matrix: aqso4 partition', grain=clock_module)
        matrix_run_clock = mpp_clock_id ('matrix: total_run_clock', grain=clock_module)
        matrix_module_init = .true.
    end subroutine matrix_init


    !----------------------------------------------------------------
    !				subroutine get_matrix_tracer_param
    !	...assign tracer properties in tracer type; 
    ! 	...assign index into population 
    !----------------------------------------------------------------
    subroutine get_matrix_tracer_param(tracer_index, population,control,mt,mp)!,r)
        integer, intent(in) :: tracer_index !!tracer_index is the index in all atmos tracers from get_number_tracers
        character(len=512), intent(in)  :: population
        character(len=512), intent(in)  :: control
        type(matrix_tracer), intent(inout) :: mt !note: this is one tracer
        type(matrix_pop),    intent(inout) :: mp(npop) !note: this is all populations in the configuration
        !real, intent(in) :: r(:,:,:,:)
        integer :: pop_index
        integer :: iflag
        pop_index = 0
        ! trim(population) name: akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
        if (lowercase(trim(population))=='akk') then
            mt%pop = 'akk'
            pop_index = i_akk
            mt%pop_index = pop_index

        elseif (lowercase(trim(population))=='acc') then
            mt%pop = 'acc'
            pop_index = i_acc
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='dd1') then
            mt%pop = 'dd1'
            pop_index = i_dd1
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='dd2') then
            mt%pop = 'dd2'
            pop_index = i_dd2
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='ssa') then
            mt%pop = 'ssa'
            pop_index = i_ssa
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='ssc') then
            mt%pop = 'ssc'
            pop_index = i_ssc
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='oc1') then
            mt%pop = 'oc1'
            pop_index = i_oc1
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='oc2') then
            mt%pop = 'oc2'
            pop_index = i_oc2
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='bc1') then
            mt%pop = 'bc1'
            pop_index = i_bc1
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='bc2') then
            mt%pop = 'bc2'
            pop_index = i_bc2
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='mxa') then
            mt%pop = 'mxa'
            pop_index = i_mxa
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='mxc') then
            mt%pop = 'mxc'
            pop_index = i_mxc
            mt%pop_index = pop_index
        elseif (lowercase(trim(population))=='ext') then
            mt%pop = 'ext'
            pop_index = i_ext
            mt%pop_index = pop_index
        else
            call error_mesg('get_matrix_tracer_param', 'trim(population) not found '//trim(trim(population)), fatal )        
        endif


        if (pop_index > 0) then !if the current configuration have this population and this tracer
            !    mt%index_in_atmos = tracer_index
            iflag=parse(control,'type',mt%type)
            if (iflag>0) then
                mt%is_active=.true.
            end if
            iflag=parse(control,'spec',mt%spec)
            iflag=parse(control,'sigma',mt%sigma)
            if (iflag>0) then
                mt%lnsigma  = log(mt%sigma)
                mt%lnsigma2 = mt%lnsigma**2
            end if
            iflag=parse(control,'dgn',mt%dgn)
            if (iflag>0) then
                mt%dp0 = mt%dgn * exp(1.5*mt%lnsigma2)
            end if
            iflag=parse(control,'distribution',mt%distribution) 
            ! distribution_index !=1, lognormal; =2, weibull
            if (trim(mt%distribution).eq."lognormal") then
                mt%distribution_index = dist_lognormal
            elseif (trim(mt%distribution).eq."weibull") then
                mt%distribution_index = dist_weibull       
            end if
            !only mass species with emission have distribution parameters in the field table 
            if (iflag>0) then
                mt%has_emission = .true. !this line only updated mass tracers with emission
            else
                mt%has_emission = .false.
            end if
            iflag=parse(control,'kappa',mt%kappa)
            iflag=parse(control,'dens',mt%dens)
            call get_tracer_names (model_atmos, tracer_index, name = mt%name,  &
                units = mt%units)

            mp(pop_index)%nb_tracer_pop = mp(pop_index)%nb_tracer_pop + 1   ! count the number of tracers in the specific populations 


            mp(pop_index)%tracer_index((mp(pop_index)%nb_tracer_pop)) = tracer_index


            mp(pop_index)%has_emission(mp(pop_index)%nb_tracer_pop) = mt%has_emission
            if (trim(mt%type).eq."number") then
                mp(pop_index)%i_n = tracer_index !record the index in the matrix_tracer_array
                mp(pop_index)%name = lowercase(trim(population))
            elseif (trim(mt%type).eq."mass") then
                if  (trim(mt%spec).eq. "sulf") then
                    mp(pop_index)%i_msulf = tracer_index           
                elseif (trim(mt%spec).eq. "dust") then
                    mp(pop_index)%i_mdust = tracer_index
                elseif (trim(mt%spec).eq. "seas") then
                    mp(pop_index)%i_mseas = tracer_index
                elseif (trim(mt%spec).eq. "ocar") then
                    mp(pop_index)%i_mocar = tracer_index
                elseif (trim(mt%spec).eq. "bcar") then
                    mp(pop_index)%i_mbcar = tracer_index
                elseif (trim(mt%spec).eq. "alwc") then
                    mp(pop_index)%i_mwate = tracer_index 
                    mp(pop_index)%name = lowercase(trim(population)) 
                elseif (trim(mt%spec).eq. "ammo") then
                    mp(pop_index)%i_mammo = tracer_index
                elseif (trim(mt%spec).eq. "nitr") then
                    mp(pop_index)%i_mnitr = tracer_index    
                elseif (trim(mt%spec).eq. "gsfa") then !gaseous sulfuric acid
                    mp(pop_index)%i_mgsfa = tracer_index
                else
                    call error_mesg ('matrix_gfdl','matrix tracer: field table property not properly defined', fatal)
                endif
            endif
        endif
    end subroutine get_matrix_tracer_param

    !----------------------------------------------------------------
    !                subroutine set_matrix_value
    ! 1. given the tracer array rt, set up matrix values in matrix unit
    !----------------------------------------------------------------
    subroutine set_matrix_value(rt, pwt, zhalf,is,ie,js,je)
        real, intent(in) :: rt(:,:,:,:)
        integer, intent(in) :: is,ie,js,je
        real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
        integer :: n !local variables
        do n=1, ntrace !loop all tracers in rt array
        if (matrix_all_tracer(n)%is_active) then
            call set_unit_gfdl_to_matrix(rt(:,:,:,n),matrix_all_tracer(n)%units, &
                matrix_all_tracer(n)%type, matrix_all_tracer(n)%spec, matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,:), pwt, zhalf)
        endif
        enddo
    end subroutine set_matrix_value

    !----------------------------------------------------------------
    !                subroutine set_unit_gfdl_to_matrix
    ! given a 3d array and its units in gfdl: vmr, mmr, #/kg
    ! convert to matrix units: ug/m3 for mass, #/m3 for number
    !----------------------------------------------------------------
    subroutine set_unit_gfdl_to_matrix(monotracer, tr_unit, tr_type, tr_spec, value_in_matrix, pwt, zhalf) ! note: avoid using unit/type, because &
        !they are used previously for other purposes
        real, intent(in) :: monotracer(:,:,:)!single tracer concentration in gfdl unit
        character*32, intent(in) :: tr_unit !tracer unit defined in field table, vmr/mmr/#/kg
        character*32, intent(in) :: tr_type, tr_spec !matrix tracer type: number/mass; tracer spec: sulf/dust/seas/ocar/bcar/alwc/gsfa
        real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
        real, intent(out) :: value_in_matrix(size(pwt,1),size(pwt,2),size(pwt,3)) !value in matrix unit
        real :: mw !local variable
        integer :: kd
        kd = size(pwt, 3)
        value_in_matrix = 0.
        if (lowercase(trim(tr_type)) .eq. "number") then !number type: convert into #/m3
            if (lowercase(trim(tr_unit)) .eq. "vmr") then
                value_in_matrix = monotracer * pwt(:,:,:)  / wtmair / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
            elseif (lowercase(trim(tr_unit)) .eq. "#/kg") then
                value_in_matrix = monotracer * pwt(:,:,:) / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
            else
                call error_mesg ('matrix_gfdl','matrix tracer: number type unit not properly defined', fatal)
            endif
            value_in_matrix = max(1.e-17, value_in_matrix)
        elseif (lowercase(trim(tr_type)) .eq. "mass") then !mass type: convert into ug/m3
            if (lowercase(trim(tr_unit)) .eq. "vmr") then !only sulfate can possible to be vmr unit
                if (lowercase(trim(tr_spec)) .eq. "sulf") then
                    mw = 96
                    value_in_matrix = 1e9 * monotracer  * pwt(:,:,:) * mw / wtmair / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
                else
                    call error_mesg ('matrix_gfdl','matrix tracer: molecular weight not properly defined', fatal)
                endif
            elseif (lowercase(trim(tr_unit)) .eq. "mmr") then
                value_in_matrix  = 1e9 * monotracer * pwt(:,:,:)  / (zhalf(:,:,1:kd) - zhalf(:,:,2:kd+1))
            else
                call error_mesg ('matrix_gfdl','matrix tracer: number type unit not properly defined', fatal)
            endif
            value_in_matrix = max(1.e-32, value_in_matrix)
        endif        
    end subroutine set_unit_gfdl_to_matrix


    !---------------------subroutine sediment flux calculation-------------------------------
    subroutine rt_dt_from_sedimentation(rt, ipop, pfull, dt, t, pwt, is, ie, js, je,  rt_dt, time_next, rh,kbot)
        real, intent(in) :: rt(:,:,:,:)
        real, intent(in) :: rh(:,:,:)
        real, intent(in) :: dt !timestep
        type(time_type), intent(in) :: time_next 
        integer, intent(in) :: ipop !matrix_pop index
        integer, intent(in) :: is, ie, js, je
        integer :: is_in, js_in, ks_in
        logical :: used
        real, intent(inout) :: rt_dt(size(rt,1),size(rt,2),size(rt,3),size(rt,4))
        real, allocatable:: setl(:,:,:,:)
        integer :: nb_tracer_pop !number of tracers in ipop
        real, intent(in) :: pfull(:,:,:), t(:,:,:), pwt(:,:,:) !t is temperature in k
        integer, intent(in), optional :: kbot(:,:) ! index of bottom level
        integer ::  i, id, j, jd, k, kd, kb, icount, n_rt
        real :: dg_wet_pop(size(rt,1), size(rt,2), size(rt,3))
        real :: d_wet_pop(size(rt,1), size(rt,2), size(rt,3))
        real :: dens_wet_pop(size(rt,1), size(rt,2), size(rt,3))
        real :: vdep(size(rt,3)), air_dens(size(rt,3)), dz(size(rt,3))
        logical :: use_sj_sedimentation_solver 
        real :: mtv 
        real :: dens_dry_pop(size(rt,1),size(rt,2),size(rt,3)), kappa_pop(size(rt,1),size(rt,2),size(rt,3))
        real :: dg_dry_pop(size(rt,1),size(rt,2),size(rt,3)), pop_num(size(rt,1),size(rt,2),size(rt,3))
        use_sj_sedimentation_solver = .true.
        mtv = 1.0
        id = size(pfull, 1)
        jd = size(pfull, 2)
        kd = size(pfull, 3)
        !        rt_dt = 0.
        nb_tracer_pop = matrix_all_pop(ipop)%nb_tracer_pop
        allocate(setl(nb_tracer_pop, id, jd, kd))
        setl = 0.
        dg_wet_pop = 0.
        icount = 0
        n_rt = 0
        !------for debug purpose
        dens_dry_pop = 0.
        kappa_pop = 0.
        dg_dry_pop = 0.
        pop_num = 0.
        !------end debug variables

        nb_tracer_pop = matrix_all_pop(ipop)%nb_tracer_pop
        dg_wet_pop = matrix_all_pop(ipop)%dg_wet(is:ie,js:je,:) !unit in m
        d_wet_pop = dg_wet_pop*exp(1.5*(log(matrix_all_pop(ipop)%sigma))**2) !unit in m
        dens_wet_pop = matrix_all_pop(ipop)%dens_wet(is:ie,js:je,:) !unit kg/m3

        !------for debug purpose
        dens_dry_pop = matrix_all_pop(ipop)%dens_dry(is:ie,js:je,:)
        kappa_pop = matrix_all_pop(ipop)%kappa_pop(is:ie,js:je,:)
        dg_dry_pop =  matrix_all_pop(ipop)%dg_dry(is:ie,js:je,:)
        pop_num = matrix_all_tracer(matrix_all_pop(ipop)%i_n)%value_in_matrix(is:ie,js:je,:)
        !------end debug purpose

        do j=1,jd
        do i=1,id
        air_dens(:)=pfull(i,j,:)/t(i,j,:)/rdgas
        dz(:) = pwt(i,j,:)/air_dens(:)
        if (present(kbot)) then
            kb=kbot(i,j)
        else
            kb=kd
        endif
        !calculate deposition velocity
        !calculate flux based on mass/number tracer
        do icount = 1, nb_tracer_pop
        n_rt = matrix_all_pop(ipop)%tracer_index(icount) !tracer index in the population
        if (lowercase(matrix_all_tracer(n_rt)%type) .eq. "mass") then
            !calculate deposition velocity
            do k=1,kb
            vdep(k) = sedimentation_velocity(t(i,j,k),pfull(i,j,k),d_wet_pop(i,j,k)/2,dens_wet_pop(i,j,k)) ! settling velocity [m/s]
            !add constrain for vdep as 30um particle sediment more than 1 grid, making sj_scheme not conserve
            vdep(k) = min(dz(k)/dt, vdep(k))
            enddo
        elseif (lowercase(matrix_all_tracer(n_rt)%type) .eq. "number") then
            do k=1,kb
            vdep(k) = sedimentation_velocity(t(i,j,k),pfull(i,j,k),dg_wet_pop(i,j,k)/2,dens_wet_pop(i,j,k)) ! settling velocity [m/s]
            !add constrain for vdep as 30um particle sediment more than 1 grid, making sj_scheme not conserve
            vdep(k) = min(dz(k)/dt, vdep(k))
            enddo
        endif



        !sedimentation_flux is designed for mmr mass tracer
        !for mass: setl unit: kg/m2/s; for number: setl unit: #/m2/s
        call sedimentation_flux(use_sj_sedimentation_solver,kb, &
            dt,mtv,dz,vdep,air_dens,&
            pwt(i,j,:), rt(i,j,:,n_rt), rt_dt(i,j,:, n_rt), setl(icount, i,j,:)) !sedimentation_flux is designed for mmr mass tracer


        !    if (mpp_root_pe().eq.mpp_pe()) then
        !            if ((abs(sum(rt_dt(i,j,:, n_rt)*pwt(i, j, :) + setl(icount,i,j,:))) / sum(setl(icount,i,j,:)) > 0.10) ) then
        !            !        .and. &
        !            !(matrix_all_tracer(n_rt)%name .eq. 'ssc_mseas')) then
        !                write(*,*) '=================== DEBUG: Sedimentation Flux ==================='
        !                write(*,*) 'pop =', ipop , matrix_all_pop(ipop)%name, 'nb_tracer_pop', nb_tracer_pop
        !                write(*,*) 'icount =', icount, matrix_all_tracer(n_rt)%name
        !                write(*,*) 'Imbalance per grid:', sum(rt_dt(i,j,:,n_rt)*pwt(i,j,:) + setl(icount,i,j,:))


        !                write(*,*) '--- Column sums ---'
        !                write(*,*) 'Tracer_dt column sum: ', sum(rt_dt(i,j,:,n_rt)*pwt(i,j,:))
        !                write(*,*) 'Tracer_dt column sum2: ',  sum(rt_dt(i,j,:,n_rt)*air_dens*dz)
        !                write(*,*) 'Setl column sum:      ', sum(setl(icount,i,j,:))
        !                write(*,*) 'Setl at kb:           ', setl(icount,i,j,kb)
        !                write(*,*) ''

        !                write(*,*) '--- Vertical Profile Table (1 to kb =', kb, ') ---'
        !                write(*,'(A)') 'z   pop_num       vdep      t         pfull     rh   kappa_pop dg_dry(um)   dg_wet(um)    dens_dry      dens_wet      dz       air_dens      pwt       rt        rt_dt     setl'
        !                do k = 1, kb
        !                    write(*,'(I2,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3)') &
        !                        k, pop_num(i,j,k), vdep(k), t(i,j,k), pfull(i,j,k), rh(i,j,k), kappa_pop(i,j,k), dg_dry_pop(i,j,k)*1e6, &
        !                        dg_wet_pop(i,j,k)*1e6, dens_dry_pop(i,j,k), dens_wet_pop(i,j,k), dz(k), air_dens(k), pwt(i,j,k), &
        !                        rt(i,j,k,n_rt), rt_dt(i,j,k,n_rt), setl(icount,i,j,k)
        !                end do
        !                write(*,*) '==============================================================='
        !            endif
        !             
        !    endif
        !
        enddo
        enddo
        enddo

        ! write(*,*), 'tracer_dt colume sum', sum(rt_dt*air_dens*dz, dim=1) !mmr to kg/m2
        ! write(*,*), 'setl values:', sum(setl, dim=1)
        ! write(*,*), 'setl kb:', setl(kb)

        !send out settling diagnostic
        do icount = 1, nb_tracer_pop
        n_rt = matrix_all_pop(ipop)%tracer_index(icount)
        ! if (mpp_root_pe().eq.mpp_pe()) then
        !        do j=1,jd
        !        do i=1,id 
        !         if (abs(sum(rt_dt(i,j,:, n_rt)*pwt(i, j, :)+setl(icount,i,j,:)))/sum(setl(icount,i,j,:)) > 0.50) then !if imbalance is larger than 1%
        !                 write(*,*), 'DEBUG Sedimentation_flux: tracer_dt and setl values'
        !                 write(*,*), 'pop = ', ipop, 'icount = ', icount, &
        !                        'imbalance per grid:', sum(rt_dt(i,j,:, n_rt)*pwt(i, j, :)+setl(icount,i,j,:))
        !                 write(*,*), 'tracer_dt colume sum', sum(rt_dt(i,j,:, n_rt)*pwt(i,j, :))
        !                 write(*,*), 'setl values', sum(setl(icount,i,j,:)), setl(icount,i,j,kb)
        !                 write(*, *), 'vdep', vdep
        !                 write(*, *), '' 

        !         endif
        !         enddo
        !         enddo
        !  endif


        if (lowercase(matrix_all_tracer(n_rt)%type) .eq. "mass") then
            if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "mmr") then
                !send out setlling as kg/m2/s unit
                used = send_data(matrix_all_tracer(n_rt)%id_tracer_setl, setl(icount, :,:,kb), time_next, &
                    is_in=is,js_in=js)
                used = send_data(matrix_all_tracer(n_rt)%id_tracer_setl_3d, rt_dt(:,:,:,n_rt),  time_next, &
                    is_in=is,js_in=js, ks_in=1)
            else
                call error_mesg ('matrix_gfdl','current sedimentation scheme is not able to treat this mass unit', fatal )
            endif
        elseif (lowercase(matrix_all_tracer(n_rt)%type) .eq. "number") then
            if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "#/kg") then
                !send out setlling as mol/m2/s unit
                used = send_data(matrix_all_tracer(n_rt)%id_tracer_setl, setl(icount, :,:,kb)/avogno, time_next, &
                    is_in=is,js_in=js)
                used = send_data(matrix_all_tracer(n_rt)%id_tracer_setl_3d, rt_dt(:,:,:,n_rt)/avogno,  time_next, &
                    is_in=is,js_in=js, ks_in=1)
            else
                call error_mesg ('matrix_gfdl','current sedimentation scheme is not able to treat this number unit', fatal )
            endif
        endif
        enddo

        deallocate(setl)

    end subroutine rt_dt_from_sedimentation





    !----------------------------------------------------------------
    !                subroutine update_rt_from_matrix(n_rt, mono_rt_tracer)
    ! given the index of tracers, update gfdl tracer values from matrix%value
    !----------------------------------------------------------------
    subroutine update_rt_from_matrix(n_rt, mono_rt_tracer,pwt, zhalf,is,ie,js,je)
        integer, intent(in) :: n_rt !index in rt(:,:,:,n_rt)
        real, intent(in), optional :: pwt(:,:,:), zhalf(:,:,:)
        integer, intent(in) :: is,ie,js,je
        real, intent(out) :: mono_rt_tracer(:,:,:) !rt(:,:,:,n_rt) to be updated
        real :: mw !molecular weight, local variable
        integer :: kd
        kd = size(mono_rt_tracer, 3)
        !number or mass tracers
        mono_rt_tracer(:,:,:) = 1.0e-32
        if (lowercase(matrix_all_tracer(n_rt)%type) .eq. "mass") then
            if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "mmr") then
                mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:)/1e9/(pwt(:,:,:)/ (zhalf(:,:,1:kd) - zhalf(:,:,2:kd+1)))
            elseif (lowercase(matrix_all_tracer(n_rt)%units) .eq. "vmr") then
                if (lowercase(trim(matrix_all_tracer(n_rt)%spec)) .eq. "sulf") then
                    mw = 96
                    mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:)/1e9/(pwt(:,:,:) * mw / wtmair / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
                else
                    call error_mesg ('matrix_gfdl','matrix tracer: molecular weight not properly defined', fatal)
                endif
            endif
        elseif (lowercase(matrix_all_tracer(n_rt)%type) .eq. "number") then
            if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "#/kg") then
                mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:) / (pwt(:,:,:) / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
            elseif (lowercase(matrix_all_tracer(n_rt)%units) .eq. "vmr") then
                mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:) / (pwt(:,:,:)  / wtmair / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
            else
                call error_mesg ('matrix_gfdl','matrix tracer: number unit not properly defined', fatal)
            endif
        endif
    end subroutine update_rt_from_matrix


    !----------------------------------------------------------------
    !                subroutine set_matrix_source_2d
    ! 1. broadcast 2d source -> 3d source
    ! 2. convert unit from vmr/s, mmr/s, kg/m2/s to µg/m3/s
    ! conversion: 
    ! 1. vmr/s -> kg/m2/s: vmr/s * pwt * mw / wtmair
    ! 2. mmr/s -> kg/m2/s: mmr/s * pwt
    ! 3. kg/m2/s -> µg/m3/s: kg/m2/s / h_zgrid * 1e9
    !----------------------------------------------------------------
    subroutine set_matrix_source_2d(source_type,source,source_units,pwt,zhalf, time, time_next, is, js)

        integer, intent(in) :: source_type, source_units
        real ,   intent(in) :: source(:,:)
        real,    intent(in) :: pwt(:,:,:), zhalf(:,:,:)
        type(time_type), intent(in) :: time, time_next ! current model time
        integer, intent(in) :: is, js ! boundaries of physical window
        real :: source_processed(size(source,1),size(source,2),k_grid_size) !kd = size(r,3)
        real :: mw ! molecular weight of the specific source type
        integer :: kd
        kd = k_grid_size
        !determine molecular weight mw
        if (source_type .eq. matrix_source_type%p_h2so4) then
            mw = matrix_molecular_weight(i_mw_h2so4)
        elseif (source_type .eq. matrix_source_type%e_so4) then
            mw = matrix_molecular_weight(i_mw_so4)
        elseif (source_type .eq. matrix_source_type%e_ss_acc) then
            mw = matrix_molecular_weight(i_mw_ss)
        elseif (source_type .eq. matrix_source_type%e_ss_coars) then
            mw = matrix_molecular_weight(i_mw_ss)
        elseif (source_type .eq. matrix_source_type%e_dust_acc) then
            mw = matrix_molecular_weight(i_mw_dust)
        elseif (source_type .eq. matrix_source_type%e_dust_coars) then
            mw = matrix_molecular_weight(i_mw_dust)
        elseif(source_type .eq. matrix_source_type%e_oc_phob) then
            mw = matrix_molecular_weight(i_mw_oc)
        elseif(source_type .eq. matrix_source_type%e_oc_phil) then
            mw = matrix_molecular_weight(i_mw_oc)
        elseif(source_type .eq. matrix_source_type%e_bc_phob) then
            mw = matrix_molecular_weight(i_mw_bc)
        elseif(source_type .eq. matrix_source_type%e_bc_phil) then
            mw = matrix_molecular_weight(i_mw_bc)
        elseif(source_type .eq. matrix_source_type%e_soa) then
            mw = matrix_molecular_weight(i_mw_soa)
        elseif(source_type .eq. matrix_source_type%p_aqso4) then
            mw = matrix_molecular_weight(i_mw_so4)
        endif

        source_processed(:,:,:) = 0.

        if (do_matrix) then
            if (source_units .eq.  matrix_source_type%u_kg_m2_s) then
                source_processed(:,:,kd) = 1e9 * source / (zhalf(:,:,kd) - zhalf(:,:,kd+1)) !convert µg/m3/s
            else
                call error_mesg ('matrix_gfdl','emission source dimension is not compatible with its unit', fatal)
            end if
            call set_matrix_source_generic(source_type,source_processed, zhalf, time, time_next, is, js)  
        end if
    end subroutine set_matrix_source_2d


    subroutine set_matrix_source_3d(source_type,source,source_units,pwt,zhalf, time, time_next, is,js)

        integer, intent(in) :: source_type, source_units
        real ,   intent(in) :: source(:,:,:)
        real,    intent(in), optional :: pwt(:,:,:), zhalf(:,:,:) 
        type(time_type), intent(in) :: time, time_next ! current model time
        integer, intent(in) :: is, js ! boundaries of physical window
        real :: mw
        real ::    source_processed(size(source,1),size(source,2),size(source,3))
        integer :: kd
        kd = k_grid_size
        !determine molecular weight mw
        if (source_type .eq. matrix_source_type%p_h2so4) then
            mw = matrix_molecular_weight(i_mw_h2so4)
        elseif (source_type .eq. matrix_source_type%e_so4) then
            mw = matrix_molecular_weight(i_mw_so4)
        elseif (source_type .eq. matrix_source_type%e_ss_acc) then
            mw = matrix_molecular_weight(i_mw_ss)
        elseif (source_type .eq. matrix_source_type%e_ss_coars) then
            mw = matrix_molecular_weight(i_mw_ss)
        elseif (source_type .eq. matrix_source_type%e_dust_acc) then
            mw = matrix_molecular_weight(i_mw_dust)
        elseif (source_type .eq. matrix_source_type%e_dust_coars) then
            mw = matrix_molecular_weight(i_mw_dust)
        elseif(source_type .eq. matrix_source_type%e_oc_phob) then
            mw = matrix_molecular_weight(i_mw_oc)
        elseif(source_type .eq. matrix_source_type%e_oc_phil) then
            mw = matrix_molecular_weight(i_mw_oc)
        elseif(source_type .eq. matrix_source_type%e_bc_phob) then
            mw = matrix_molecular_weight(i_mw_bc)
        elseif(source_type .eq. matrix_source_type%e_bc_phil) then
            mw = matrix_molecular_weight(i_mw_bc)
        elseif(source_type .eq. matrix_source_type%e_soa) then
            mw = matrix_molecular_weight(i_mw_soa)
        elseif(source_type .eq. matrix_source_type%p_aqso4) then
            mw = matrix_molecular_weight(i_mw_so4)
        endif

        source_processed(:,:,:) = 0.

        if (do_matrix) then
            if (source_units .eq.  matrix_source_type%u_vmr_s) then
                source_processed(:,:,:) = 1e9 * source(:,:,:) * pwt(:,:,:) * mw / wtmair / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
            elseif (source_units .eq.  matrix_source_type%u_mmr_s) then
                source_processed(:,:,:) = 1e9 * source(:,:,:) * pwt(:,:,:)  / (zhalf(:,:,1:kd) - zhalf(:,:,2:kd+1))
            else
                call error_mesg ('matrix_gfdl','emission source dimension is not compatible with its unit', fatal)
            end if
            call set_matrix_source_generic(source_type,source_processed, zhalf, time, time_next, is,js)
        end if

    end subroutine set_matrix_source_3d

    subroutine set_matrix_source_generic(source_type,source_processed, zhalf, time, time_next, is,js)

        real, intent(in)    :: source_processed(:,:,:), zhalf(:,:,:)
        integer, intent(in) :: source_type
        type(time_type), intent(in) :: time, time_next ! current model time
        integer, intent(in) :: is, js ! boundaries of physical window
        integer :: ie, je
        integer :: nk
        real :: h2so4_source_2d(size(source_processed,1),size(source_processed,2)), so4_emis(size(source_processed,1),size(source_processed,2))
        logical :: used 
        real :: dust_emis_scale 
        integer :: kd
        kd = k_grid_size
        ie = is+size(source_processed,1)-1
        je = js+size(source_processed,2)-1
        h2so4_source_2d = 0.
        so4_emis = 0.
        dust_emis_scale = 1.0
        if (do_matrix) then
            !akk mode h2so4(g): vmr/s, 3d array
            if (source_type .eq. matrix_source_type%p_h2so4) then   
                p_h2so4_rate(is:ie,js:je,:) = source_processed
                do nk=1,kd
                h2so4_source_2d = h2so4_source_2d + source_processed(:,:,nk)*(zhalf(:,:,nk)-zhalf(:,:,nk+1))*1e-9 !convert ug/m3/s -> kg/m2/s
                enddo
                used = send_data (id_h2so4_source, h2so4_source_2d, time_next, &
                    is_in=is,js_in=js)
                !acc mode so4 emission: vmr/s, 3d array
            elseif (source_type .eq. matrix_source_type%e_so4) then
                if (i_akk > 0 .and. i_acc < 0) then
                    matrix_all_tracer(matrix_all_pop(i_akk)%i_msulf)%source(is:ie,js:je,:) =0.01* source_processed
                elseif (i_akk > 0 .and. i_acc > 0) then
                    matrix_all_tracer(matrix_all_pop(i_acc)%i_msulf)%source(is:ie,js:je,:) =0.99* source_processed
                    matrix_all_tracer(matrix_all_pop(i_akk)%i_msulf)%source(is:ie,js:je,:) =0.01* source_processed
                elseif (i_akk < 0 .and. i_acc > 0) then
                    matrix_all_tracer(matrix_all_pop(i_acc)%i_msulf)%source(is:ie,js:je,:) =source_processed
                endif
                do nk=1,kd
                so4_emis = so4_emis + source_processed(:,:,nk)*(zhalf(:,:,nk)-zhalf(:,:,nk+1))*1e-9 !convert ug/m3/s -> kg/m2/s
                enddo
                used = send_data (id_so4_emis, so4_emis, time_next, &
                    is_in=is,js_in=js)
            elseif (source_type .eq. matrix_source_type%p_aqso4) then
                aqso4(is:ie,js:je,:) =  source_processed
                used = send_data (id_aqso4_rate, aqso4(is:ie,js:je,:), time_next, is_in=is,js_in=js, ks_in=1)
                !dust emission: kg/m2/s, 2d array
            elseif (source_type .eq. matrix_source_type%e_dust_acc) then
                !test scale dust emission source as 0.5

                if (i_dd1>0 .and. i_dd2<0) then !only have 1 dust mode dd1
                    matrix_all_tracer(matrix_all_pop(i_dd1)%i_mdust)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_dd1)%i_mdust)%source(is:ie,js:je,:) + dust_emis_scale*source_processed
                elseif (i_dd1<0 .and. i_dd2>0) then !only have 1 dust mode dd2
                    matrix_all_tracer(matrix_all_pop(i_dd2)%i_mdust)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_dd2)%i_mdust)%source(is:ie,js:je,:) + dust_emis_scale*source_processed
                elseif (i_dd1>0 .and. i_dd2>0) then !have 2 dust mode: dd1 and dd2
                    matrix_all_tracer(matrix_all_pop(i_dd1)%i_mdust)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_dd1)%i_mdust)%source(is:ie,js:je,:) + dust_emis_scale*source_processed
                endif
            elseif (source_type .eq. matrix_source_type%e_dust_coars) then
                if (i_dd1>0 .and. i_dd2<0) then !only have 1 dust mode dd1
                    matrix_all_tracer(matrix_all_pop(i_dd1)%i_mdust)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_dd1)%i_mdust)%source(is:ie,js:je,:) + dust_emis_scale*source_processed
                elseif (i_dd1<0 .and. i_dd2>0) then !only have 1 dust mode dd2
                    matrix_all_tracer(matrix_all_pop(i_dd2)%i_mdust)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_dd2)%i_mdust)%source(is:ie,js:je,:) + dust_emis_scale*source_processed
                elseif (i_dd1>0 .and. i_dd2>0) then !have 2 dust mode: dd1 and dd2
                    matrix_all_tracer(matrix_all_pop(i_dd2)%i_mdust)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_dd2)%i_mdust)%source(is:ie,js:je,:) + dust_emis_scale*source_processed
                endif
                !sea salt emission: kg/m2/s, 2d array
            elseif (source_type .eq. matrix_source_type%e_ss_acc) then
                if (i_ssa>0 .and. i_ssc<0) then !only have 1 sea salt mode ssa
                    matrix_all_tracer(matrix_all_pop(i_ssa)%i_mseas)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_ssa)%i_mseas)%source(is:ie,js:je,:) + source_processed
                elseif (i_ssa<0 .and. i_ssc>0) then !only have 1 sea salt mode ssc
                    matrix_all_tracer(matrix_all_pop(i_ssc)%i_mseas)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_ssc)%i_mseas)%source(is:ie,js:je,:) + source_processed
                elseif (i_ssa>0 .and. i_ssc>0) then !have 2 sea salt mode: ssa and ssc
                    matrix_all_tracer(matrix_all_pop(i_ssa)%i_mseas)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_ssa)%i_mseas)%source(is:ie,js:je,:) + source_processed
                endif
            elseif (source_type .eq. matrix_source_type%e_ss_coars) then
                if (i_ssa>0 .and. i_ssc<0) then !only have 1 sea salt mode ssa
                    matrix_all_tracer(matrix_all_pop(i_ssa)%i_mseas)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_ssa)%i_mseas)%source(is:ie,js:je,:) + source_processed
                elseif (i_ssa<0 .and. i_ssc>0) then !only have 1 sea salt mode ssc
                    matrix_all_tracer(matrix_all_pop(i_ssc)%i_mseas)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_ssc)%i_mseas)%source(is:ie,js:je,:) + source_processed
                elseif (i_ssa>0 .and. i_ssc>0) then !have 2 sea salt mode: ssa and ssc
                    matrix_all_tracer(matrix_all_pop(i_ssc)%i_mseas)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_ssc)%i_mseas)%source(is:ie,js:je,:) + source_processed
                endif
                !soa source from gfdl model, add it to oc: kg/m2/s, 2d array
            elseif (source_type .eq. matrix_source_type%e_soa) then
                if (i_oc1>0 .and. i_oc2<0) then
                    matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) + source_processed
                elseif (i_oc1<0 .and. i_oc2>0) then
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) + source_processed
                elseif (i_oc1>0 .and. i_oc2>0) then !have 2 organic carbon mode
                    !matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) = &
                    !    & matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) + 0.2*source_processed
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) = & 
                        & matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) + source_processed
                endif
                !organic carbon: vmr, 3d array
            elseif (source_type .eq. matrix_source_type%e_oc_phob) then
                if (i_oc1>0 .and. i_oc2<0) then
                    matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) + source_processed
                elseif (i_oc1<0 .and. i_oc2>0) then
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) + source_processed
                elseif (i_oc1>0 .and. i_oc2>0) then !have 2 organic carbon mode
                    matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) + source_processed
                endif
            elseif (source_type .eq. matrix_source_type%e_oc_phil) then
                if (i_oc1>0 .and. i_oc2<0) then
                    matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%source(is:ie,js:je,:) + source_processed
                elseif (i_oc1<0 .and. i_oc2>0) then
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) + source_processed
                elseif (i_oc1>0 .and. i_oc2>0) then
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) = &
                        & matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%source(is:ie,js:je,:) + source_processed
                endif
                !black carbon: vmr, 3d array
            elseif (source_type .eq. matrix_source_type%e_bc_phob) then
                if (i_bc1>0 .and. i_bc2<0) then
                    matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%source(is:ie,js:je,:)+source_processed
                elseif (i_bc1<0 .and. i_bc2>0) then
                    matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%source(is:ie,js:je,:) + source_processed
                elseif (i_bc1>0 .and. i_bc2>0) then !have 2 organic carbon mode
                    matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%source(is:ie,js:je,:) + source_processed
                endif

            elseif (source_type .eq. matrix_source_type%e_bc_phil) then
                if (i_bc1>0 .and. i_bc2<0) then
                    matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%source(is:ie,js:je,:)+source_processed
                elseif (i_bc1<0 .and. i_bc2>0) then
                    matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%source(is:ie,js:je,:) + source_processed
                elseif (i_bc1>0 .and. i_bc2>0) then !have 2 organic carbon mode
                    matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%source(is:ie,js:je,:) = &
                        matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%source(is:ie,js:je,:) + source_processed
                endif
            endif
        endif

    end subroutine set_matrix_source_generic

    !-----------------------------------------------------------------------
    !                 set_matrix_emis_number
    ! 1. convert mass concentration into number concentration
    ! 2. assign number information to tracers with number type
    !-----------------------------------------------------------------------
    subroutine set_matrix_emis_number(is,ie,js,je)
        integer :: n,nt,ntt 
        integer, intent(in) :: is,ie,js,je
        do n=1,npop
        if (matrix_all_pop(n)%nb_tracer_pop > 0) then !if this population exist in current configuration
            do nt = 1,matrix_all_pop(n)%nb_tracer_pop
            if (matrix_all_pop(n)%has_emission(nt)) then
                ntt = matrix_all_pop(n)%tracer_index(nt)
                matrix_all_tracer(matrix_all_pop(n)%i_n)%source(is:ie,js:je,:) = &
                    matrix_all_tracer(ntt)%source(is:ie,js:je,:)/(pi6*matrix_all_tracer(ntt)%dens*matrix_all_tracer(ntt)%dp0**3)*1e-9 !emission rate converted to number rate

            endif
            end do
        endif
        end do
    end subroutine set_matrix_emis_number


    !-----------------------------------------------------------------------
    !                 set_matrix_pop_kci
    ! calculate condensational coefficient for different populations 
    !-----------------------------------------------------------------------
    subroutine set_matrix_pop_kci(pfull,t,kci_coef_pop, kci_aeq1_pop,is,ie,js,je) !calculate the condensation coefficient for all population
        real, intent(in) :: pfull(:,:,:) ! pressure on layers, pa
        real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degk
        real, intent(inout) :: kci_coef_pop(:,:,:,:)
        real, intent(inout) :: kci_aeq1_pop(:,:,:,:) !kci units: m^3/s
        integer :: n,it,jt,kt,i,j,k
        real :: fac,kci1,kci2
        real :: number_pop(size(pfull,1),size(pfull,2),size(pfull,3))
        integer, intent(in) :: is,ie,js,je
        kci_coef_pop = 0.
        kci_aeq1_pop = 0.
        number_pop = 1.0e-17
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        do n=1, npop
        if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
            fac = exp(1.5*log(matrix_all_pop(n)%sigma)**2)
            number_pop = matrix_all_tracer(matrix_all_pop(n)%i_n)%value_in_matrix(is:ie,js:je,:)
            do i=1,it
            do j=1,jt
            do k=1,kt
            !if ((number_pop(i,j,k) > 1e-17) .and. (m_dry_all(i,j,k) > 1e-32)) then
            if (number_pop(i,j,k) > 1.0e-14) then
                call setup_kci(pfull(i,j,k), t(i,j,k), matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k)*fac, matrix_all_pop(n)%sigma, &
                    kci_coef_pop(n,i,j,k), kci_aeq1_pop(n,i,j,k))
            endif
            enddo
            enddo
            enddo

        endif
        enddo

    end subroutine


    subroutine set_matrix_condense_npf(i_akk, xnh3, fland, pfull,rh,t, tstep, h2so4_vmr_0,pwt,zhalf, kci_coef_pop, kci_aeq1_pop, &
            dndt_npf, dmdt_h2so4_npf, ih2so4_path, dmdt_h2so4_tot,kc, xh2so4_nucl,is,ie,js,je)
        real, intent(in) :: pfull(:,:,:) ! pressure on layers, pa
        real, intent(in) :: rh(:,:,:) ! relative humidity
        real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
        real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degk
        real, intent(in) :: tstep !timestep, in s
        real, intent(in) :: h2so4_vmr_0(:,:,:) !gaseous h2so4 in gfdl unit, vmr
        real, intent(in) :: kci_coef_pop(:,:,:,:), kci_aeq1_pop(:,:,:,:) !m^3/s
        real, intent(in) :: xnh3(:,:,:), fland(:,:,:)
        integer, intent(in) :: i_akk
        integer, intent(in) :: is,ie,js,je
        integer, intent(inout) :: ih2so4_path(:,:,:)
        real, intent(inout):: dndt_npf(:,:,:), dmdt_h2so4_npf(:,:,:), dmdt_h2so4_tot(:,:,:) !m-3 s-1; ugso4 m-3 s-1
        integer :: n,it, jt, kt, i,j,k ! number of layers in different direction
        real :: h2so4_vmr(size(h2so4_vmr_0,1), size(h2so4_vmr_0,2), size(h2so4_vmr_0,3))
        real :: mw_so4
        real :: mw_h2so4
        real, parameter :: tinydenom = 1.0d-30
        real, intent(out):: xh2so4_nucl(size(pfull,1),size(pfull,2),size(pfull,3))  ! h2so4 (as so4) conc. used in nucleation and gr calculation [ugso4/m^3]
        real, parameter :: xntau =2.0 !number of time consants in the current time step
        real, parameter :: kcmin = 1.0d-08 ! [1/s] minimum condensational sink - see notes of 10-18-06
        real, parameter :: xh2so4_nucl_min_ncm3 = 1.00d+03 ! min. [h2so4] to enter nucleation calculations [#/cm^3]
        real :: xh2so4_nucl_min ! convert to [ugh2so4/m^3]
        real :: xh2so4_init(size(pfull,1),size(pfull,2),size(pfull,3)) !xh2so4_init: h2so4 concentration in ug/m3
        real,intent(inout) :: kc(:,:,:) !kc: total condensation sink (1/s) for all population in a grid
        real :: kc_aeq1(size(pfull,1),size(pfull,2),size(pfull,3)) !kc_aeq1: total condensation sink (1/s) using accomadation coefficient=1
        real :: number_pop(size(pfull,1),size(pfull,2),size(pfull,3)) !number concentration #/m3
        real :: h2so4rate(size(pfull,1),size(pfull,2),size(pfull,3)) ! average h2so4 production rate [ugso4/m^3/s]
        real :: xh2so4_ss(size(pfull,1),size(pfull,2),size(pfull,3)), xh2so4_ss_wnpf(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: pq_growth(size(pfull,1),size(pfull,2),size(pfull,3)), tot_h2so4_loss(size(pfull,1),size(pfull,2),size(pfull,3))
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        h2so4_vmr = 0.0
        mw_so4 = matrix_molecular_weight(i_mw_so4)
        mw_h2so4 = matrix_molecular_weight(i_mw_h2so4)
        xh2so4_nucl_min = xh2so4_nucl_min_ncm3 * mw_h2so4 * 1.0d+12 / avogno  ! convert to [ugh2so4/m^3]
        xh2so4_nucl =  xh2so4_nucl_min !tinynumer ! xh2so4_nucl_min              ! for the case  xh2so4_init .lt. xh2so4_nucl_min
        xh2so4_ss_wnpf = 0.
        pq_growth = 0.
        h2so4_vmr = h2so4_vmr_0
        if (minval(h2so4_vmr_0) < 0) then
            call error_mesg('set_matrix_condense_npf','h2so4_vmr input is negative', warning)
            h2so4_vmr = max(h2so4_vmr, 0.0)
        endif
        xh2so4_init = 1e9 * h2so4_vmr * pwt * mw_h2so4 / wtmair / (zhalf(:,:,1:kt) - zhalf(:,:,2:kt+1)) !convert gaseous h2so4 to matrix unit: ug h2so4/m3

        kc=0. !3d, total condensation sink at each grid, sum-up the population
        kc_aeq1 =0. !3d, total condensation sink at each grid, sum-up the population
        h2so4rate =0. !average so4 production rate [ugso4/m^3/s]
        dmdt_h2so4_tot = 0.
        dmdt_h2so4_npf = 0.
        dndt_npf = 0.
        do n=1,npop
        if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
            number_pop = matrix_all_tracer(matrix_all_pop(n)%i_n)%value_in_matrix(is:ie,js:je,:) !number concentration of a population
            kc = kc + kci_coef_pop(n,:,:,:)*number_pop !m3/s * #/m3 = 1/s, 3d for each grid
            kc_aeq1 = kc_aeq1 + kci_aeq1_pop(n,:,:,:)*number_pop !m3/s * #/m3 = 1/s, 3d for each grid
            !-----------test if inf exist
            if (mpp_root_pe().eq.mpp_pe()) then
                if (any( (number_pop>huge(0.0)) .or. (number_pop < -huge(0.0)) )) then
                    call error_mesg('set_matrix_condense_npf infinity exist in number_pop', matrix_all_pop(n)%name,fatal)
                endif
            endif
            !-----------------------------------
        endif
        enddo
        kc = max(1e-32, kc)
        kc_aeq1 = max(1e-32, kc)

        !test if inf exist
        if (mpp_root_pe().eq.mpp_pe()) then
            if (any((kci_coef_pop > huge(0.0)) .or. (kci_coef_pop < -huge(0.0)))) then
                call error_mesg('set_matrix_condense_npf', 'infinity exist in kci_coef_pop', fatal)
            endif
            if (any((kci_aeq1_pop > huge(0.0)) .or. (kci_aeq1_pop < -huge(0.0)))) then
                call error_mesg('set_matrix_condense_npf', 'infinity exist in kci_aeq1_pop', fatal)
            endif
            if (any((kc > huge(0.0)) .or. (kc < -huge(0.0)))) then
                call error_mesg('set_matrix_condense_npf', 'infinity exist in kc', fatal)
            endif
            if (any((kc_aeq1 > huge(0.0)) .or. (kc_aeq1 < -huge(0.0)))) then
                call error_mesg('set_matrix_condense_npf', 'infinity exist in kc', fatal)
            endif
        endif

        if (i_akk > 0) then
            h2so4rate = xh2so4_init / tstep                        ! average h2so4 production rate [ugso4/m^3/s]
            do i=1,it
            do j=1,jt
            do k=1,kt
            if ((kc(i,j,k)*tstep) .ge. xntau) then ! invoke steady-state assumption
                ih2so4_path(i,j,k) =1
                xh2so4_ss(i,j,k) = min( h2so4rate(i,j,k)/kc(i,j,k), xh2so4_init(i,j,k) )  !steady-state h2so4 [ugh2so4/m^3]

                call steady_state_h2so4(pfull(i,j,k),rh(i,j,k),t(i,j,k),fland(i,j,k),xh2so4_ss(i,j,k), &
                    h2so4rate(i,j,k),xnh3(i,j,k),kc(i,j,k),tstep,xh2so4_ss_wnpf(i,j,k))
                xh2so4_nucl(i,j,k) = xh2so4_ss_wnpf(i,j,k)                       ! [h2so4] for nucl., gr, and cond. calculation [ugso4/m^3]
            else
                ih2so4_path(i,j,k) =2
                xh2so4_nucl(i,j,k) = h2so4rate(i,j,k) / ( (2.0d+00/tstep) + kc(i,j,k) )   ! use [h2so4] at mid-time step [ugso4/m^3]
            endif


            enddo
            enddo
            enddo
            !xl: place to call matrix_npf
            call matrix_npfrate(pfull, rh, t, fland, xh2so4_nucl, h2so4rate, kc_aeq1,dndt_npf,dmdt_h2so4_npf)

            if (mpp_root_pe().eq.mpp_pe()) then
                if (any((xh2so4_nucl > huge(0.0)) .or. (xh2so4_nucl < -huge(0.0)))) then
                    call error_mesg('set_matrix_condense_npf', 'infinity exist in xh2so4_nucl', fatal)
                endif
                if (any((dmdt_h2so4_npf > huge(0.0)) .or. (dmdt_h2so4_npf < -huge(0.0)))) then
                    call error_mesg('set_matrix_condense_npf', 'infinity exist in dmdt_h2so4_npf', fatal)
                endif
            endif

        else 
            dndt_npf = 0.
            dmdt_h2so4_npf = 0.
            xh2so4_nucl = xh2so4_init
        endif

        do i=1,it
        do j=1,jt
        do k=1,kt
        pq_growth(i,j,k) = xh2so4_nucl(i,j,k) * ( 1.0d+00 - exp(-kc(i,j,k)*tstep) ) / tstep            ! [ugh2so4/m^3/s]
        tot_h2so4_loss(i,j,k) = ( dmdt_h2so4_npf(i,j,k) + pq_growth(i,j,k) ) * tstep                         ! [ugh2so4/m^3]
        if ( tot_h2so4_loss(i,j,k) .gt. xh2so4_init(i,j,k) ) then
            dmdt_h2so4_npf(i,j,k)  = dmdt_h2so4_npf(i,j,k)  * ( xh2so4_init(i,j,k) / ( tot_h2so4_loss(i,j,k)+tinydenom ) )! [ugh2so4/m^3/s]
            dndt_npf(i,j,k)  = dndt_npf(i,j,k)  * ( xh2so4_init(i,j,k) / ( tot_h2so4_loss(i,j,k) + tinydenom ))! [  #  /m^3/s]
            pq_growth(i,j,k) = pq_growth(i,j,k) * ( xh2so4_init(i,j,k) / ( tot_h2so4_loss(i,j,k) + tinydenom ))! [ugh2so4/m^3/s]
        endif
        dmdt_h2so4_tot(i,j,k) = dmdt_h2so4_npf(i,j,k)+pq_growth(i,j,k)
        enddo
        enddo
        enddo

        if (mpp_root_pe().eq.mpp_pe()) then
            if (any((dmdt_h2so4_tot > huge(0.0)) .or. (dmdt_h2so4_tot < -huge(0.0)))) then
                call error_mesg('set_matrix_condense_npf', 'infinity exist in dmdt_h2so4_tot', fatal)
            endif
            if (any((dmdt_h2so4_npf > huge(0.0)) .or. (dmdt_h2so4_npf < -huge(0.0)))) then
                call error_mesg('set_matrix_condense_npf', 'infinity exist in dmdt_h2so4_npf', fatal)
            endif
        endif

    end subroutine

    !-----------------------------------------------------------------------
    !                 matrix_npfrate
    ! 1. set-up for loop to call npfrate in matrix
    ! 2. calculate dndt and dmdt in 3d array form
    !-----------------------------------------------------------------------
    !matrix_npf(pfull, rh, t, so4rate, xh2so4_nucl,kc_aeq1,dndt,dmdt_so4)
    !npfrate(prs,rh,temp,xh2so4,so4rate,kc,dndt,dmdt_so4,icall)!xl
    !matrix_npfrate(pfull, rh, t, fland, xh2so4_nucl, so4rate, kc_aeq1,dndt_npf,dmdt_npf)
    subroutine matrix_npfrate(pfull,rh,t, fland, h2so4, so4_rate,kc_aeq1,dndt, dmdt_h2so4)
        real, intent(in) :: pfull(:,:,:) ! pressure on layers, pa
        real, intent(in) :: rh(:,:,:) ! relative humidity
        !real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
        real, intent(in) :: so4_rate(:,:,:) !h2so4 production rate in ugso4/m^3
        real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degk
        real, intent(in) :: h2so4(:,:,:) !gaseous h2so4 in ugso4/m^3
        real, intent(in) :: fland(:,:,:)
        real, intent(in) :: kc_aeq1(:,:,:) !condensation sink with accomadation coefficient =1
        real, intent(out):: dndt(:,:,:), dmdt_h2so4(:,:,:)
        integer :: it, jt, kt, i,j,k ! number of layers in different direction
        real:: pres_mt, rh_mt, temp_mt, h2so4_mt, p_h2so4rate_mt !local variable for matrix calculation
        real:: mw, dndt_tst, dmdt_tst !if directly assign dndt(i,j,k) and dmdt(i,j,k) segment error will occur
        real:: dndt_rec(size(pfull,1),size(pfull,2),size(pfull,3)), dmdt_rec(size(pfull,1),size(pfull,2),size(pfull,3)) 
        mw = matrix_molecular_weight(i_mw_h2so4)
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        do i=1,it
        do j=1,jt
        do k=1,kt
        pres_mt = pfull(i,j,k) !pressure in matrix unit: [pa]
        rh_mt = rh(i,j,k) !fractional relative humidity [1]
        !xl: check with fp for rh calculation
        temp_mt = t(i,j,k) !ambient temperature [k]
        !sulfuric acid (as so4) concentration [ugso4/m^3]
        h2so4_mt = h2so4(i,j,k)
        p_h2so4rate_mt = so4_rate(i,j,k) !gas-phase h2so4 (as so4) production rate [ugso4/m^3 s]
        !dndt: [m^-3 s^-1], dmdt: [ugso4 m^-3 s^-1]
        !npfrate(prs,rh,temp,fland,xh2so4,so4rate,xnh3,kc,dndt,dmdt_so4,icall)
        !call npfrate(pres_mt,rh_mt,temp_mt,h2so4_mt,p_h2so4rate_mt,dndt(i,j,k),dmdt(i,j,k))
        call npfrate(pres_mt,rh_mt,temp_mt,0.,h2so4_mt,p_h2so4rate_mt,0.,kc_aeq1(i,j,k),dndt_rec(i,j,k),dmdt_rec(i,j,k),0)
        enddo
        enddo
        enddo
        dndt=dndt_rec
        dmdt_h2so4=dmdt_rec
    end subroutine matrix_npfrate
    !----------------------------------------------------------------------------------------------------------------
    ! get the b_i loss       terms due to intermodal coagulation. [1/s]
    ! get the r_i production terms due to intermodal coagulation. [#/m^3/s]
    ! get the c_i terms, which include all source terms.          [#/m^3/s]
    ! for the c_i terms, the secondary particle formation term dndt must be
    !   added in after coupling to condensation below.
    ! the a_i terms for intramodal coagulation are directly computed
    !   from the coagulation coefficients when the number equations
    !   are integrated.
    ! if dikl(i,k,l) = 0, then modes k and l to not coagulate to form mode i.
    ! dij(i,j) is unity if coagulation of mode i with mode j results
    !   in the removal of particles from mode i, and zero otherwise.
    !----------------------------------------------------------------------------------------------------------------
    subroutine matrix_coag_cap_n(pfull, t, is,ie,js,je, time, time_next, tstep, ri, li, rim, lim, lon, lat)
        implicit none
        integer(kind=i8_kind) :: chksum_val
        real, intent(in),    dimension(:,:)           :: lon, lat
        logical :: used
        type(time_type),  intent(in) :: time, time_next
        integer, intent(in) :: is,ie,js,je
        real, intent(in) :: pfull(:,:,:), t(:,:,:)
        real, intent(in) :: tstep
        real :: lon_val, lat_val, max_diff
        integer:: i_loc, j_loc, k_loc, loc(3)
        integer :: index_qq, i_spec
        real, parameter :: tinydenom = 1.0d-26 !initialization 1e-32 mmr -> 1e-23, 1e-26 ug/m3 -> rate: 1e-32mmr/1800s -> 1e-26
        real, parameter :: tinydenon = 1.0d-20 !initilization 1e-17 #/m3 -> 1e-20 number rate
        real :: bi(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !loss terms due to intermodal coagulation. [1/s]
        real :: diff(size(pfull,1),size(pfull,2),size(pfull,3))
        real, intent(out) :: li(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !loss term of number due to intermodal coagulation [#/m^3/s]
        real, intent(out) :: ri(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !production terms due to intermodal coagulation. [#/m^3/s]
        real :: fi(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !mass loss coefficient
        real, intent(out) :: lim(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass loss term of number due to intermodal coagulation [ug/m^3/s]
        real, intent(out) :: rim(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass production terms due to intermodal coagulation. [ug/m^3/s]
        real :: kbar0_ij(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3)) !mode averaged coagulation coefficient [m3/s]
        real :: kbar0_ij_cap(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3))
        real :: kbar3_ij_cap(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3))
        real :: kbar3_ij(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3)) !mode averaged coagulation coefficient [m3/s]
        real :: dg_ip_um(size(pfull,1),size(pfull,2),size(pfull,3)), dg_jp_um(size(pfull,1),size(pfull,2),size(pfull,3)) !particle geometric diameter [um]
        real :: num_pop_k(size(pfull,1),size(pfull,2),size(pfull,3)), num_pop_l(size(pfull,1),size(pfull,2),size(pfull,3)), num_pop_i(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: num_pop_j(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: mjq(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass of species in a single particle: mass/number = ug/particle
        real :: sig_ip, sig_jp, pop_numi(size(pfull,1),size(pfull,2),size(pfull,3)), pop_numj(size(pfull,1),size(pfull,2),size(pfull,3))
        integer :: ip,jp,kp,it,jt,kt,i,j,k,l,ikl,ipop,kpop,lpop,q,qq,klq
        real :: sig_2pop(2), diam_2pop(2), kbar0ij_interp(2,2), kbar3ij_interp(2,2)
        real :: num_loss_max_rate(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: mass_loss_max_rate(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3))
        real :: n_scale, m_scale
        !record and updated number and mass concentration inside this subroutine
        real :: num_all_pop(npop-1, size(pfull,1),size(pfull,2),size(pfull,3)), num_all_pop_initial(npop-1, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: num_all_pop_tmp(npop-1, size(pfull,1),size(pfull,2),size(pfull,3)), dn_interm_li(npop-1, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: mass_all_pop(npop-1, 5, size(pfull,1),size(pfull,2),size(pfull,3)), mass_all_pop_initial(npop-1, 5, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: mass_all_pop_tmp(npop-1, 5, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: dni_f_kpop(size(pfull,1),size(pfull,2),size(pfull,3)), dni_f_lpop(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: diff_array(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: m_scale_qq, m_scale_used, max_bound_limit, value_kpop, value_lpop
        integer :: max_indices(3), i_mspec
        real :: expdt, pqin_thresh, expbit(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: dm(size(pfull,1),size(pfull,2),size(pfull,3)), dm_dt(size(pfull,1),size(pfull,2),size(pfull,3)), dm_thresh
        real :: dm_12_seas(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: diag_n(size(pfull,1),size(pfull,2),size(pfull,3)), diag_m(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: kbar3n(npop-1, size(pfull,1),size(pfull,2),size(pfull,3)), kbar3n_f(npop-1, size(pfull,1),size(pfull,2),size(pfull,3))
        integer :: ig, jg, kg, nspecc, np
        real :: m_spec_ini(5, size(pfull,1),size(pfull,2),size(pfull,3)), m_spec_fnl(5, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: sum_kbar3n_f_dm(npop-1, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: dm_interm_lm(npop-1, 5, size(pfull,1),size(pfull,2),size(pfull,3)), dm_interm_pm(npop-1, 5, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: dm_interm_pm_q(5, size(pfull,1),size(pfull,2),size(pfull,3)), dm_interm_lm_q(5, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: dm_scale_pm(5, size(pfull,1),size(pfull,2),size(pfull,3)),  dm_scale_lm(5, size(pfull,1),size(pfull,2),size(pfull,3))
        real :: dm_imbalance_p_l(5), mean_grid
        integer :: n
        mean_grid = 0.
        dm_imbalance_p_l = 0.
        dm_scale_pm = 1.0
        dm_scale_lm = 1.0
        dm_interm_pm_q = 0.
        dm_interm_lm_q = 0.
        sum_kbar3n_f_dm = 0.
        dm_interm_lm = 0.
        dm_interm_pm =0.
        diag_n = 0.
        diag_m = 0.
        dm_dt = 0.
        dm_12_seas = 0.
        dm_thresh = 1.0e-32
        dm = 0.
        !real, intent(out) ::
        dn_interm_li = 0.
        i_mspec = 0
        max_bound_limit = 0.5
        m_scale_qq = 1.0
        m_scale_used = 1.0
        dni_f_kpop = 0.
        dni_f_lpop = 0.
        num_all_pop = 0.
        mass_all_pop = 0.
        num_all_pop_initial = 0.
        mass_all_pop_initial = 0.
        num_all_pop_tmp = 0.
        mass_all_pop_tmp = 0.
        n_scale = 0.
        m_scale = 0.
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        num_loss_max_rate = 0.
        mass_loss_max_rate = 0.
        sig_2pop = 1.
        diam_2pop = 1e-30
        kbar0_ij = 0
        kbar3_ij = 0
        kbar0_ij_cap = 0.
        kbar3_ij_cap = 0.
        kbar0ij_interp = 0
        kbar3ij_interp = 0
        bi = 0.
        ri = 0. !production of number due to coagulation
        li = 0.
        fi = 0.
        lim = 0.
        rim = 0.
        mjq = 0. !ug/particle

        !constaruct mjq
        !mjq: mjq(j,q) is the avg. mass/particle of species q (=1-5) for mode j
        call mpp_clock_begin(coag_sub_clock1)
        do ip=1,npop-1 !exclude 'ext' pop
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_pop_i = matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            num_all_pop(ip, :,:,:) = num_pop_i
            num_all_pop_initial(ip,:,:,:) = num_pop_i 
            do q = 1, nm(ip)
            qq=prod_index(ip,q) !index of mass species: 1-5
            select case (qq)
            case (1)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_msulf)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_all_pop(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_msulf)%value_in_matrix(is:ie,js:je,:)
                mass_all_pop_initial(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_msulf)%value_in_matrix(is:ie,js:je,:)
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_msulf)%value_in_matrix(is:ie,js:je,:)/tstep*max_bound_limit
            case (4)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mdust)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_all_pop(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mdust)%value_in_matrix(is:ie,js:je,:)
                mass_all_pop_initial(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mdust)%value_in_matrix(is:ie,js:je,:)
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mdust)%value_in_matrix(is:ie,js:je,:)/tstep*max_bound_limit
            case (5)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mseas)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_all_pop(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mseas)%value_in_matrix(is:ie,js:je,:)
                mass_all_pop_initial(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mseas)%value_in_matrix(is:ie,js:je,:)
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mseas)%value_in_matrix(is:ie,js:je,:)/tstep*max_bound_limit
            case (3)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_all_pop(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(is:ie,js:je,:)
                mass_all_pop_initial(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(is:ie,js:je,:)
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(is:ie,js:je,:)/tstep*max_bound_limit
            case (2)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mbcar)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_all_pop(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mbcar)%value_in_matrix(is:ie,js:je,:)
                mass_all_pop_initial(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mbcar)%value_in_matrix(is:ie,js:je,:)
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mbcar)%value_in_matrix(is:ie,js:je,:)/tstep*max_bound_limit
            case default !call fatal error
                call error_mesg('matrix coagulation','mjq mass species not properly defined ', fatal)
            end select
            enddo
        endif
        enddo
        mjq = max(0.0, mjq)
        call mpp_clock_end(coag_sub_clock1) 


        call mpp_clock_begin(coag_sub_clock2)
        !set up kbar0_ij and kbar3_ij table (npop, npop, it, jt, kt)
        do ip = 1, npop-1 !exclude 'ext' pop
        do jp = ip, npop-1
        !do jp = 1, npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(jp)%nb_tracer_pop > 0)) then
            dg_ip_um = matrix_all_pop(ip)%dg_wet(is:ie,js:je,:)*1.0e6 !3d array [um]
            sig_ip = matrix_all_pop(ip)%sigma !real number
            dg_jp_um = matrix_all_pop(jp)%dg_wet(is:ie,js:je,:)*1.0e6 ! [um]
            sig_jp = matrix_all_pop(jp)%sigma
            pop_numi = matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            pop_numj = matrix_all_tracer(matrix_all_pop(jp)%i_n)%value_in_matrix(is:ie,js:je,:)
            do i=1,it
            do j=1,jt
            do k=1,kt
            if ((dg_ip_um(i,j,k) > 0) .and. (dg_jp_um(i,j,k) > 0)) then
                sig_2pop=(/sig_ip, sig_jp/)
                diam_2pop = (/dg_ip_um(i,j,k), dg_jp_um(i,j,k)/)
                call get_kbarnij(1, t(i,j,k), pfull(i,j,k), 2, sig_2pop, diam_2pop, kbar0ij_interp, kbar3ij_interp)
                kbar0_ij(ip, jp, i, j,k) = kbar0ij_interp(1, 2)
                kbar3_ij(ip, jp, i, j,k) = kbar3ij_interp(1, 2)
                kbar0_ij(jp, ip, i, j,k) = kbar0ij_interp(2, 1)
                kbar3_ij(jp, ip, i, j,k) = kbar3ij_interp(2, 1)
            else
                kbar0_ij(ip, jp, i, j,k) = 0.
                kbar3_ij(ip, jp, i, j,k) = 0.
                kbar0_ij(jp, ip, i, j,k) = 0.
                kbar3_ij(jp, ip, i, j,k) = 0.
            endif
            enddo
            enddo
            enddo

        endif
        enddo
        enddo
        kbar0_ij = max(kbar0_ij, 0.)
        kbar3_ij = max(kbar3_ij, 0.)
        kbar0_ij_cap = kbar0_ij
        kbar3_ij_cap = kbar3_ij
        !test kbar0 eq. kbar3
        kbar3_ij = kbar0_ij
        kbar3_ij_cap = kbar0_ij_cap
        call mpp_clock_end(coag_sub_clock2)


        !--------------------------------------------------------------------- 
        !first do self-coagulation and update the number
        !---------------------------------------------------------------------
        num_all_pop_initial = num_all_pop
        num_all_pop_tmp = num_all_pop
        li=0.
        do ip = 1, npop-1
        num_loss_max_rate = 0.
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_loss_max_rate = num_all_pop_tmp(ip,:,:,:)/tstep*max_bound_limit !max_bound_limit=0.5 
            li(ip, :,:,:) = 0.5 * kbar0_ij_cap(ip,ip,:,:,:) *num_all_pop_tmp(ip,:,:,:) *num_all_pop_tmp(ip,:,:,:)
            do i=1,it
            do j=1,jt
            do k=1,kt
            if (li(ip, i,j,k) > num_loss_max_rate(i,j,k)) then
                n_scale = (num_loss_max_rate(i,j,k)/(li(ip, i,j,k)+tinydenon))
                kbar0_ij_cap(ip,ip,i,j,k) = kbar0_ij_cap(ip,ip,i,j,k) * n_scale
            endif    
            enddo
            enddo
            enddo 
            !dN/dt = -k*N^2 --> N = N0/(1+N0*k*t)
            num_all_pop(ip,:,:,:) = num_all_pop_tmp(ip,:,:,:)/(1.0 + num_all_pop_tmp(ip,:,:,:)*kbar0_ij_cap(ip,ip,:,:,:)*tstep)  
        endif
        enddo

        if (any(ieee_is_nan(num_all_pop))) then
            if(mpp_pe() == mpp_root_pe()) then
                write (*, *), "1loc coag num self coag: nan found in num_all_pop"
            endif
        endif


        !------------------------------------------------------------------------------
        !second do intermodal number loss and production and update the number
        !------------------------------------------------------------------------------
        num_all_pop_tmp = num_all_pop !this is the initial number for intermodal transfer and will be used for mass calculation
        li=0. !hypothetical loss number
        bi=0.
        do ip =1, npop-1
        bi(ip,:,:,:) = 0 !number loss coeficient
        num_loss_max_rate = 0.
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_loss_max_rate = num_all_pop_tmp(ip,:,:,:)/tstep*max_bound_limit !max_bound_limit=0.5 
            do kp = 1,npop-1
            if ((dij(ip, kp) > 0) .and. (matrix_all_pop(kp)%nb_tracer_pop > 0)) then !dij=1 means i+j result in the number loss of i into other pops, dij not symmetric
                bi(ip,:,:,:) = bi(ip,:,:,:) + kbar0_ij_cap(ip,kp,:,:,:)*num_all_pop_tmp(kp,:,:,:)
            endif
            enddo
            li(ip, :,:,:) = bi(ip,:,:,:)*num_all_pop_tmp(ip,:,:,:)
            do i=1,it
            do j=1,jt
            do k=1,kt
            if (li(ip, i,j,k) > num_loss_max_rate(i,j,k)) then
                n_scale = (num_loss_max_rate(i,j,k)/(li(ip, i,j,k)+tinydenon))
                do kp = 1,npop-1
                if ((dij(ip, kp) > 0) .and. (matrix_all_pop(kp)%nb_tracer_pop > 0)) then
                    kbar0_ij_cap(ip,kp,i,j,k) = kbar0_ij_cap(ip,kp,i,j,k)*n_scale
                    kbar0_ij_cap(ip,kp,i,j,k) = min(kbar0_ij_cap(ip,kp,i,j,k), kbar0_ij_cap(kp,ip,i,j,k)) !kbar0 must be symetric, and must be the minimum value
                    kbar0_ij_cap(kp,ip, i,j,k) = kbar0_ij_cap(ip,kp,i,j,k)
                endif
                enddo
            endif
            enddo
            enddo
            enddo
            !recalculate bi
            bi(ip,:,:,:) = 0.
            do kp = 1,npop-1
            if ((dij(ip, kp) > 0) .and. (matrix_all_pop(kp)%nb_tracer_pop > 0)) then
                bi(ip,:,:,:) = bi(ip,:,:,:) + kbar0_ij_cap(ip,kp,:,:,:)*num_all_pop_tmp(kp,:,:,:)
            endif
            enddo
            !update num_all_pop(ip,:,:,:)
            !dN/dt = -bi*N -> N = N0*exp(-b*t)
            pqin_thresh = 1.0e-32
            do i=1,it
            do j=1,jt
            do k=1,kt
            expdt = exp(-bi(ip,i,j,k)*tstep)
            num_all_pop(ip,i,j,k) = num_all_pop_tmp(ip,i,j,k)*expdt
            enddo
            enddo
            enddo
        endif

        enddo

        if (any(ieee_is_nan(num_all_pop))) then
            if(mpp_pe() == mpp_root_pe()) then
                write (*, *), "2loc coag intermodal loss: nan found in num_all_pop"
            endif
        endif

        !------------------------------------------------------------------------------
        !third do intermodal number production and update the number
        !note: !!!!!do not update number_all_pop_tmp
        !num_all_pop_tmp record the number information right after self-coagulation, but before intermodal transfer
        !------------------------------------------------------------------------------
        dn_interm_li = num_all_pop_tmp - num_all_pop !positive value: this record the absolute change of numbers due to intermodal transfer loss
        if (minval(dn_interm_li) < 0.) then
            if(mpp_pe() == mpp_root_pe()) then
                write (*, *), "WARNING: matrix coag number loss has increased the number"
            endif
        endif


        ri=0.
        do ikl = 1, ndikl !ikl is 1-12, exclude ext already: dikl means k+l -> number in i, this means dij(k,l) and dij(l,k) both =1
        ipop = dikl_control(ikl)%i
        kpop = dikl_control(ikl)%k
        lpop = dikl_control(ikl)%l
        if ((matrix_all_pop(ipop)%nb_tracer_pop > 0) .and. (matrix_all_pop(kpop)%nb_tracer_pop > 0) &
            .and. (matrix_all_pop(lpop)%nb_tracer_pop > 0)) then
            if ((dij(kpop, lpop) > 0) .and. (dij(lpop, kpop) > 0)) then !k+l result k number loss, but i number increase
                ! this will also reult l number loss, but i number increase
                do i=1,it
                do j=1,jt
                do k=1,kt
                if (bi(kpop,i,j,k) .gt. 0) then
                    dni_f_kpop(i,j,k) = & !number change calculated from pop k
                        (kbar0_ij_cap(kpop,lpop,i,j,k)*num_all_pop_tmp(lpop,i,j,k))/bi(kpop,i,j,k) &
                        *dn_interm_li(kpop, i,j,k) !number change calculated from pop k
                endif
                if (bi(lpop,i,j,k) .gt. 0) then !!number change calculated from pop l
                    dni_f_lpop(i,j,k) = &
                        (kbar0_ij_cap(lpop,kpop,i,j,k)*num_all_pop_tmp(kpop,i,j,k))/bi(lpop,i,j,k) &
                        *dn_interm_li(lpop,i,j,k)
                endif

                enddo
                enddo
                enddo
                num_all_pop(ipop,:,:,:) = num_all_pop(ipop,:,:,:) + 1.0/2.0*(dni_f_kpop+dni_f_lpop)
            else
                call error_mesg('matrix coagulation','din and dikl not consistent ', fatal)
            endif
        endif
        enddo

        if (any(ieee_is_nan(num_all_pop))) then
            if(mpp_pe() == mpp_root_pe()) then
                write (*, *), "3loc coag intermodal prod: nan found in num_all_pop"
            endif
        endif


        !assign number tendency to the matrix_source
        do n=1,npop-1
        if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then

            matrix_all_tracer(matrix_all_pop(n)%i_n)%source(is:ie,js:je,:) = &
                matrix_all_tracer(matrix_all_pop(n)%i_n)%source(is:ie,js:je,:) + &
                (num_all_pop(n,:,:,:)-num_all_pop_initial(n,:,:,:))/tstep


            if (any(ieee_is_nan(matrix_all_tracer(matrix_all_pop(n)%i_n)%source(is:ie,js:je,:)))) then
                if(mpp_pe() == mpp_root_pe()) then
                    write (*, *), "4loc final number tendency: nan found in num_all_pop"
                endif
            endif


            !send number increase due to intermodal tranfer
            ri(n,:,:,:) = (num_all_pop(n,:,:,:)-(num_all_pop_tmp(n,:,:,:)-dn_interm_li(n,:,:,:)))/tstep
            !number decrease due to intermodal transfer and self-coagulation
            li(n,:,:,:) =  (dn_interm_li(n,:,:,:)+num_all_pop_initial(n,:,:,:)-num_all_pop_tmp(n,:,:,:))/tstep
        endif
        enddo


        !now start calculation of mass
        !-----------------------------------------------------------------------------------
        !first do intermodal coagulation mass loss
        !-----------------------------------------------------------------------------------
        mass_all_pop_tmp = mass_all_pop
        fi = 0.
        lim = 0.
        do ip = 1, npop-1
        do jp = 1, npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(jp)%nb_tracer_pop > 0)) then
            if( dij(ip,jp) > 0 ) then !dij>0 means i+j make i loss mass to other populations
                fi(ip,:,:,:) = fi(ip,:,:,:) + kbar3_ij_cap(ip,jp,:,:,:)*num_all_pop_tmp(jp,:,:,:) 
            endif
        endif
        enddo
        !--------------------------Fabien suggested solution--------------------------------
        do i=1,it
        do j=1,jt
        do k=1,kt
        m_scale = 1.0
        do q=1, nm(ip)
        qq = prod_index(ip,q)
        lim(ip,qq,i,j,k) = lim(ip,qq,i,j,k) + fi(ip,i,j,k) * mass_all_pop_tmp(ip,qq,i,j,k)
        if (lim(ip,qq, i,j,k) > mass_loss_max_rate(ip,qq,i,j,k)) then
            m_scale_qq = mass_loss_max_rate(ip,qq,i,j,k)/(lim(ip,qq,i,j,k)+tinydenom)
            m_scale = min(m_scale,  m_scale_qq)
        endif
        enddo
        kbar3_ij_cap(ip,:,i,j,k) = min(kbar3_ij_cap(ip,:,i,j,k), kbar3_ij(ip,:,i,j,k)*m_scale)
        enddo
        enddo
        enddo

        !---------------------------Uriel debuged solution
        ! m_scale = 1.0
        ! do q=1, nm(ip)
        ! qq = prod_index(ip,q)
        ! !lim_i = (SIGMA_j (kbar3_ij*N_j)) * N_i * Miq -> fi(ip,:,:,:)*mass_all_pop_tmp(ip,qq,:,:,:)
        ! lim(ip,qq,:,:,:) = lim(ip,qq,:,:,:) + fi(ip,:,:,:) * mass_all_pop_tmp(ipqq,:,:,:)
        ! do i=1,it
        ! do j=1,jt
        ! do k=1,kt
        ! if (lim(ip,qq, i,j,k) > mass_loss_max_rate(ip,qq,i,j,k)) then
        !     m_scale_qq = ( mass_loss_max_rate(ip,qq,i,j,k)/(lim(ip,qq,i,j,k)+tinydenom))
        !     !m_scale = min(m_scale, m_scale_qq) 
        !     m_scale_used = min(m_scale, m_scale_qq)
        !     !kbar3_ij_cap(ip,:,i,j,k) = min(kbar3_ij_cap(ip,:,i,j,k), kbar3_ij(ip,:,i,j,k)*m_scale)
        !     kbar3_ij_cap(ip,:,i,j,k) = min(kbar3_ij_cap(ip,:,i,j,k), kbar3_ij(ip,:,i,j,k)*m_scale_used)
        ! endif
        ! enddo
        ! enddo
        ! enddo
        ! enddo
        !---------------------------------------------------------------------------------
        ! recalculate fi
        fi(ip,:,:,:) = 0
        do jp = 1, npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(jp)%nb_tracer_pop > 0)) then
            if( dij(ip,jp) > 0 ) then
                fi(ip,:,:,:) = fi(ip,:,:,:) + kbar3_ij_cap(ip,jp,:,:,:)*num_all_pop_tmp(jp,:,:,:)
                !if (mpp_root_pe().eq.mpp_pe()) then
                !         write(*,*) 'kbar3ij_cap: 1,5,31,', ip,jp, kbar3_ij_cap(ip,jp,1,5,31),&
                !                'num_all_pop_tmp(jp,:,:,:)', num_all_pop_tmp(jp,1,5,31)
                !endif
            endif
        endif

        !if (mpp_root_pe().eq.mpp_pe()) then
        !        write(*,*) '--------SUMMARY---------ip', ip, 'fi', fi(ip,1,5,31)
        !endif
        enddo
        !update new mass due to intermodal mass transfer loss
        !dMi,q/dt= -fi*Mi,q => Mi=Mi_0*exp(-fi*t) 
        do i=1,it
        do j=1,jt
        do k=1,kt               
        expdt = exp(-fi(ip,i,j,k) * tstep)
        pqin_thresh = 1.0e-32
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0)) then
            do q=1, nm(ip)
            qq = prod_index(ip,q)
            !if ((1.0-expdt) .gt. pqin_thresh) then
            mass_all_pop(ip,qq,i,j,k) = mass_all_pop_initial(ip,qq,i,j,k)*expdt
            enddo

        endif
        enddo
        enddo
        enddo

        enddo


        if (any(ieee_is_nan(mass_all_pop))) then
            if(mpp_pe() == mpp_root_pe()) then
                write (*, *), "1_2loc coag mass interm coag: nan found in mass_all_pop"
            endif
        endif
        !----------------------------------------------------------------------------------------
        !finally, calculate the intermodal mass increase
        !in structure g_iklq means k+l will generate q in i, the contributor is pop l
        !---------------------------------------------------------------------------------------
        !test KN =FI or not
        kbar3n = 0.
        kbar3n_f = 0.
        sum_kbar3n_f_dm = 0.
        !--------------------------------------------------------------
        mass_all_pop_tmp = mass_all_pop
        dm_interm_lm = mass_all_pop_initial - mass_all_pop

        do ip = 1, npop-1
        dm_12_seas = 0.
        if (matrix_all_pop(ip)%nb_tracer_pop > 0 .and. giklq_control(ip)%n > 0) then
            do klq = 1, giklq_control(ip)%n
            k = giklq_control(ip)%k(klq)
            l = giklq_control(ip)%l(klq)
            qq = giklq_control(ip)%qq(klq)
            !the second population l is the donor
            !e.g. MXA MXA OC1; there is no MXA OC1 MXA
            !e.g. MXA OC1 BC1; MXA BC1 OC1
            if (matrix_all_pop(k)%nb_tracer_pop > 0 .and. matrix_all_pop(l)%nb_tracer_pop > 0 ) then
                !need to check if qq is a mass species in k
                select case(qq)
                case (1)
                    i_mspec = matrix_all_pop(l)%i_msulf
                case (2)
                    i_mspec = matrix_all_pop(l)%i_mbcar
                case (3)
                    i_mspec = matrix_all_pop(l)%i_mocar
                case (4)
                    i_mspec = matrix_all_pop(l)%i_mdust
                case (5)
                    i_mspec = matrix_all_pop(l)%i_mseas
                end select
                if (i_mspec > 0) then

                    !if(mpp_pe() == mpp_root_pe()) then
                    !              write(*,*) 'ipec=', i_mspec, 'ip', ip, 'k', k, 'l', l, 'qq', qq
                    !endif

                    do ig =1, it
                    do jg =1, jt
                    do kg =1, kt
                    expdt = exp(-fi(l,ig,jg,kg) * tstep)
                    pqin_thresh = 1.0e-32
                    !!if ((1.0-expdt) .gt. pqin_thresh) then
                    !if (abs(mass_all_pop_initial(l,qq,ig,jg,kg)-mass_all_pop_tmp(l,qq,ig,jg,kg)) .gt. 0.) then        
                    mass_all_pop(ip,qq,ig,jg,kg) = mass_all_pop(ip,qq,ig,jg,kg) + &
                        (kbar3_ij_cap(l,k,ig,jg,kg)*num_all_pop_tmp(k,ig,jg,kg))/fi(l,ig,jg,kg)* &
                        dm_interm_lm(l,qq,ig,jg,kg)
                    !(mass_all_pop_initial(l,qq,ig,jg,kg)-mass_all_pop_tmp(l,qq,ig,jg,kg))

                    !endif
                    !ratio is normal ourside of the if statement, with minratio and maxtario is integer
                    !kbar3n(l,ig,jg,kg) = kbar3n(l,ig,jg,kg) + (kbar3_ij_cap(l,k,ig,jg,kg)*num_all_pop_tmp(k,ig,jg,kg))
                    enddo
                    enddo
                    enddo
                endif
            endif
            enddo
        endif
        enddo

        dm_interm_pm = mass_all_pop - mass_all_pop_tmp

        !check imbalance due to computational errors
        dm_interm_pm_q = sum(dm_interm_pm, dim=1)  !calculate dm_mass in all modes 
        dm_interm_lm_q = sum(dm_interm_lm, dim=1)
        dm_scale_pm = 1.0 !dimension 5, i, j, k
        dm_scale_lm = 1.0 !dimension 5, i, j, k
        do i=1,it
        do j=1,jt
        do k=1,kt
        dm_imbalance_p_l = dm_interm_pm_q(:, i, j, k) - dm_interm_lm_q(:, i, j, k)
        do qq =1, 5
        if (abs(dm_imbalance_p_l(qq)) .gt. 0) then
            mean_grid = 0.5*(dm_interm_pm_q(qq, i, j, k) + dm_interm_lm_q(qq, i, j, k))
            if ((dm_interm_pm_q(qq, i, j, k)>0) .and. (dm_interm_lm_q(qq, i, j, k)>0)) then
                dm_scale_pm(qq, i, j, k) = mean_grid/dm_interm_pm_q(qq, i, j, k)
                dm_scale_lm(qq, i, j, k) = mean_grid/dm_interm_lm_q(qq, i, j, k)
                dm_interm_pm(:,qq, i, j, k) = dm_interm_pm(:,qq, i, j, k)*dm_scale_pm(qq, i, j, k)
                dm_interm_lm(:, qq, i, j, k) = dm_interm_lm(:, qq, i, j, k)* dm_scale_lm(qq, i, j, k)
            else
                !if (mpp_root_pe().eq.mpp_pe()) then
                !        write(*,*) 'DEBUG PM at i, j, k = ', i, j, k, '----------------------'
                !        do npop=1, 12
                !                write (*,*) 'npop = ', npop, 'qq = ', qq, &
                !                'dm_interm_pm(npop, qq, i, j, k)',dm_interm_pm(npop, qq, i, j, k)
                !        enddo
                !        write(*,*) 'DEBUG LM at i, j, k=', i, j, k, '----------------------'
                !        do npop=1, 12
                !                write (*,*) 'npop = ', npop, 'qq=', qq, &
                !                        'dm_interm_lm(npop, qq, i, j, k)',dm_interm_lm(npop, qq, i, j, k)
                !        enddo
                !        write(*,*) 'DEBUG : DENOMINATOR IS 0'
                !        write (*,*) 'qq=', qq, 'dm_interm_pm_q(qq, i, j, k)', dm_interm_pm_q(qq, i, j, k)
                !        write (*,*) 'qq=', qq, 'dm_interm_lm_q(qq, i, j, k)', dm_interm_lm_q(qq, i, j, k)
                !endif
                dm_interm_pm(:,qq, i, j, k) = 0.
                dm_interm_lm(:, qq, i, j, k) = 0.

            endif
        endif
        !check
        !if (sum(dm_interm_pm(:,qq, i, j, k)) .ne. sum(dm_interm_lm(:, qq, i, j, k))) then
        !        if (mpp_root_pe().eq.mpp_pe()) then
        !                write(*,*) 'WARNING: dm_lm and dm_pm still not balance'
        !                write(*,*) 'WARNING: dm_lm = ', sum(dm_interm_pm(:,qq, i, j, k))
        !                write(*,*) 'WARNING: dm_pm = ', sum(dm_interm_lm(:,qq, i, j, k))
        !                do npop=1, 12
        !                        write (*,*) 'npop = ', npop, 'qq = ', qq, &
        !                                'i, j, k=', i, j, k, &
        !                                'dm_interm_pm(npop, qq, i, j, k)',dm_interm_pm(npop, qq, i, j, k), &
        !                                'dm_interm_lm(npop, qq, i, j, k)',dm_interm_lm(npop, qq, i, j, k)
        !                enddo
        !       endif
        !endif
        enddo
        enddo
        enddo
        enddo

        !mass_all_pop = mass_all_pop_tmp + dm_interm_pm
        !mass_all_pop_initial = mass_all_pop_tmp + dm_interm_lm


        !if (mpp_root_pe().eq.mpp_pe()) then
        !!m_spec_ini = sum(mass_all_pop_initial, dim=1) !sum for all modes
        !!m_spec_fnl = sum(mass_all_pop, dim=1) !sum for all modes
        !m_spec_ini = sum(dm_interm_pm, dim=1)
        !m_spec_fnl = sum(dm_interm_lm, dim=1)
        !    do nspecc = 1, 5
        !    ! Compute the absolute difference
        !    diff = abs(m_spec_ini(nspecc,:,:,:) - m_spec_fnl(nspecc,:,:,:))
        !
        !    ! Find the maximum value and its location (i,j,k)
        !    max_diff = maxval(diff)
        !    loc = maxloc(diff)  ! Returns [i, j, k]
        !
        !    ! Extract indices (i,j,k)
        !    i_loc = loc(1)
        !    j_loc = loc(2)
        !    k_loc = loc(3)  ! If k is relevant (e.g., vertical level)
        !
        !    ! Get corresponding lon and lat
        !    lon_val = lon(i_loc, j_loc)
        !    lat_val = lat(i_loc, j_loc)
        !
        !    ! Print results
        !    write(*,*) 'BALANCE_SUMMARY_ nspec', nspecc, &
        !        'max imbalance per grid:', max_diff, &
        !        'at (i,j,k):', i_loc, j_loc, k_loc, dm_scale_pm(nspecc, i_loc, j_loc, k_loc), &
        !         'at (i,j,k):', dm_scale_lm(nspecc, i_loc, j_loc, k_loc), &
        !        'pm_scale: max, min = ', maxval(dm_scale_pm(nspecc, :, :,:)), minval(dm_scale_pm(nspecc, :, :,:)), &
        !        'lm_scale: max, min = ', maxval(dm_scale_lm(nspecc, :, :,:)), minval(dm_scale_lm(nspecc, :, :,:))
        !    !    'lon, lat:', lon_val, lat_val, &
        !    !    'diff_1_5_31', abs(m_spec_ini(nspecc,1,5,31) - m_spec_fnl(nspecc,1,5,31))
        !    enddo
        !endif        

        !since 
        !mass_all_pop = mass_all_pop_tmp + dm_interm_pm
        !mass_all_pop_initial = mass_all_pop_tmp + dm_interm_lm
        lim = dm_interm_lm/tstep !coagulation mass loss rate
        rim = dm_interm_pm/tstep !cogaulation mass increase rate

        !assign mass tendency to the matrix_source
        do n=1,npop-1
        if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
            if (matrix_all_pop(n)%i_msulf > 0) then
                index_qq = 1
                i_spec = matrix_all_pop(n)%i_msulf
                dm = dm_interm_pm(n,index_qq,:,:,:) - dm_interm_lm(n,index_qq,:,:,:) !mass_all_pop(n,index_qq,:,:,:)-mass_all_pop_initial(n,index_qq,:,:,:)
                dm_dt = dm/tstep
                matrix_all_tracer(i_spec)%source(is:ie,js:je,:) = &
                    matrix_all_tracer(i_spec)%source(is:ie,js:je,:) + dm_dt


                !lim(n, index_qq, :,:,:) = dm_interm_lm(n,index_qq,:,:,:)/tstep !(mass_all_pop_initial(n,1,:,:,:)-mass_all_pop_tmp(n,1,:,:,:))/tstep
                !rim(n, index_qq, :,:,:) = dm_interm_pm(n,index_qq,:,:,:)/tstep !(mass_all_pop(n,1,:,:,:)-mass_all_pop_tmp(n,1,:,:,:))/tstep
            endif

            if (matrix_all_pop(n)%i_mbcar > 0) then
                index_qq = 2
                i_spec = matrix_all_pop(n)%i_mbcar
                dm = dm_interm_pm(n,index_qq,:,:,:) - dm_interm_lm(n,index_qq,:,:,:) !mass_all_pop(n,index_qq,:,:,:)-mass_all_pop_initial(n,index_qq,:,:,:)
                dm_dt = dm/tstep

                matrix_all_tracer(i_spec)%source(is:ie,js:je,:) = &
                    matrix_all_tracer(i_spec)%source(is:ie,js:je,:) + dm_dt

                !lim(n, 2, :,:,:) = (mass_all_pop_initial(n,2,:,:,:)-mass_all_pop_tmp(n,2,:,:,:))/tstep
                !rim(n, 2, :,:,:) = (mass_all_pop(n,2,:,:,:)-mass_all_pop_tmp(n,2,:,:,:))/tstep

            endif

            if (matrix_all_pop(n)%i_mocar > 0) then
                index_qq = 3
                i_spec = matrix_all_pop(n)%i_mocar
                dm = dm_interm_pm(n,index_qq,:,:,:) - dm_interm_lm(n,index_qq,:,:,:) !mass_all_pop(n,index_qq,:,:,:)-mass_all_pop_initial(n,index_qq,:,:,:)
                dm_dt = dm/tstep
                matrix_all_tracer(i_spec)%source(is:ie,js:je,:) = &
                    matrix_all_tracer(i_spec)%source(is:ie,js:je,:) + dm_dt

                !lim(n, 3, :,:,:) = (mass_all_pop_initial(n,3,:,:,:)-mass_all_pop_tmp(n,3,:,:,:))/tstep
                !rim(n, 3, :,:,:) = (mass_all_pop(n,3,:,:,:)-mass_all_pop_tmp(n,3,:,:,:))/tstep
            endif

            if (matrix_all_pop(n)%i_mdust > 0) then
                index_qq = 4
                i_spec = matrix_all_pop(n)%i_mdust
                dm = dm_interm_pm(n,index_qq,:,:,:) - dm_interm_lm(n,index_qq,:,:,:) !mass_all_pop(n,index_qq,:,:,:)-mass_all_pop_initial(n,index_qq,:,:,:)
                dm_dt = dm/tstep
                matrix_all_tracer(i_spec)%source(is:ie,js:je,:) = &
                    matrix_all_tracer(i_spec)%source(is:ie,js:je,:) + dm_dt

                !lim(n,4, :,:,:) = (mass_all_pop_initial(n,4,:,:,:)-mass_all_pop_tmp(n,4,:,:,:))/tstep
                !rim(n,4, :,:,:) = (mass_all_pop(n,4,:,:,:)-mass_all_pop_tmp(n,4,:,:,:))/tstep

            endif

            if (matrix_all_pop(n)%i_mseas > 0) then
                index_qq = 5
                i_spec =  matrix_all_pop(n)%i_mseas                
                dm = dm_interm_pm(n,index_qq,:,:,:) - dm_interm_lm(n,index_qq,:,:,:) !mass_all_pop(n,index_qq,:,:,:)-mass_all_pop_initial(n,index_qq,:,:,:)
                dm_dt = dm/tstep
                matrix_all_tracer(i_spec)%source(is:ie,js:je,:) = &
                    matrix_all_tracer(i_spec)%source(is:ie,js:je,:) + dm_dt

                !lim(n,5, :,:,:) = (mass_all_pop_initial(n,5,:,:,:)-mass_all_pop_tmp(n,5,:,:,:))/tstep
                !rim(n,5, :,:,:) = (mass_all_pop(n,5,:,:,:)-mass_all_pop_tmp(n,5,:,:,:))/tstep

            endif

        endif
        enddo



    end subroutine





    subroutine matrix_coag_cap(pfull, t, ri, bi, li, rim, fi, lim, kbar0_ij_cap, mjq, is,ie,js,je, time, time_next, tstep)
        implicit none
        logical :: used
        type(time_type),  intent(in) :: time, time_next
        integer, intent(in) :: is,ie,js,je
        real, intent(in) :: pfull(:,:,:), t(:,:,:)
        real, intent(in) :: tstep
        real, parameter :: tinydenom = 1.0d-26 !initialization 1e-32 mmr -> 1e-23, 1e-26 ug/m3 -> rate: 1e-32mmr/1800s -> 1e-26
        real, parameter :: tinydenon = 1.0d-20 !initilization 1e-17 #/m3 -> 1e-20 number rate
        real, intent(out) :: bi(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !loss terms due to intermodal coagulation. [1/s]
        real, intent(out) :: li(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !loss term of number due to intermodal coagulation [#/m^3/s]
        real, intent(out) :: ri(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !production terms due to intermodal coagulation. [#/m^3/s]
        real, intent(out) :: fi(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !mass loss coefficient
        real, intent(out) :: lim(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass loss term of number due to intermodal coagulation [ug/m^3/s]
        real, intent(out) :: rim(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass production terms due to intermodal coagulation. [ug/m^3/s]
        real :: kbar0_ij(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3)) !mode averaged coagulation coefficient [m3/s]
        real, intent(out) :: kbar0_ij_cap(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3))
        real :: kbar3_ij_cap(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3))
        real :: kbar3_ij(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3)) !mode averaged coagulation coefficient [m3/s]
        real :: dg_ip_um(size(pfull,1),size(pfull,2),size(pfull,3)), dg_jp_um(size(pfull,1),size(pfull,2),size(pfull,3)) !particle geometric diameter [um]
        real :: num_pop_k(size(pfull,1),size(pfull,2),size(pfull,3)), num_pop_l(size(pfull,1),size(pfull,2),size(pfull,3)), num_pop_i(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: num_pop_j(size(pfull,1),size(pfull,2),size(pfull,3))
        real, intent(out) :: mjq(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass of species in a single particle: mass/number = ug/particle
        real :: sig_ip, sig_jp, pop_numi(size(pfull,1),size(pfull,2),size(pfull,3)), pop_numj(size(pfull,1),size(pfull,2),size(pfull,3))
        integer :: ip,jp,kp,it,jt,kt,i,j,k,l,ikl,ipop,kpop,lpop,q,qq,klq
        real :: sig_2pop(2), diam_2pop(2), kbar0ij_interp(2,2), kbar3ij_interp(2,2)
        real :: num_loss_max_rate(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: mass_loss_max_rate(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3))
        real :: n_scale, m_scale
        !sig_2pop, diam_2pop, kbar0ij_interp, kbar3ij_interp
        !real, intent(out) ::
        n_scale = 0.
        m_scale = 0.
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        num_loss_max_rate = 0.
        mass_loss_max_rate = 0.
        sig_2pop = 1.
        diam_2pop = 1e-30
        kbar0_ij = 0
        kbar3_ij = 0
        kbar0_ij_cap = 0.
        kbar3_ij_cap = 0.
        kbar0ij_interp = 0
        kbar3ij_interp = 0
        bi = 0.
        ri = 0. !production of number due to coagulation
        li = 0.
        fi = 0.
        lim = 0.
        rim = 0.
        mjq = 0. !ug/particle

        !constaruct mjq
        !mjq: mjq(j,q) is the avg. mass/particle of species q (=1-5) for mode j
        call mpp_clock_begin(coag_sub_clock1)
        do ip=1,npop-1 !exclude 'ext' pop
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_pop_i = matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            do q = 1, nm(ip)
            qq=prod_index(ip,q) !index of mass species: 1-5
            !if (mpp_root_pe().eq.mpp_pe()) then
            !    write(*,*) 'ip, nm(ip), q, qq', ip, nm(ip), q, qq
            !endif
            select case (qq)
            case (1)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_msulf)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_msulf)%value_in_matrix(is:ie,js:je,:)/tstep 
            case (4)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mdust)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mdust)%value_in_matrix(is:ie,js:je,:)/tstep
            case (5)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mseas)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mseas)%value_in_matrix(is:ie,js:je,:)/tstep
            case (3)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(is:ie,js:je,:)/tstep
            case (2)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mbcar)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                mass_loss_max_rate(ip, qq, :,:,:) = &
                    matrix_all_tracer(matrix_all_pop(ip)%i_mbcar)%value_in_matrix(is:ie,js:je,:)/tstep
            case default !call fatal error
                call error_mesg('matrix coagulation','mjq mass species not properly defined ', fatal)
            end select
            enddo
        endif
        enddo
        mjq = max(0.0, mjq)
        call mpp_clock_end(coag_sub_clock1)

        used = send_data(id_mjq11, mjq(1,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1) !akk: mass/num
        used = send_data(id_mjq21, mjq(2,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1) !acc: mass/num

        call mpp_clock_begin(coag_sub_clock2)
        !set up kbar0_ij and kbar3_ij table (npop, npop, it, jt, kt)
        do ip = 1, npop-1 !exclude 'ext' pop
        do jp = ip, npop-1
        !do jp = 1, npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(jp)%nb_tracer_pop > 0)) then
            dg_ip_um = matrix_all_pop(ip)%dg_wet(is:ie,js:je,:)*1.0e6 !3d array [um]
            sig_ip = matrix_all_pop(ip)%sigma !real number
            dg_jp_um = matrix_all_pop(jp)%dg_wet(is:ie,js:je,:)*1.0e6 ! [um]
            sig_jp = matrix_all_pop(jp)%sigma
            pop_numi = matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            pop_numj = matrix_all_tracer(matrix_all_pop(jp)%i_n)%value_in_matrix(is:ie,js:je,:)
            do i=1,it
            do j=1,jt
            do k=1,kt
            if ((dg_ip_um(i,j,k) > 0) .and. (dg_jp_um(i,j,k) > 0)) then
                sig_2pop=(/sig_ip, sig_jp/)
                diam_2pop = (/dg_ip_um(i,j,k), dg_jp_um(i,j,k)/)
                call get_kbarnij(1, t(i,j,k), pfull(i,j,k), 2, sig_2pop, diam_2pop, kbar0ij_interp, kbar3ij_interp)
                kbar0_ij(ip, jp, i, j,k) = kbar0ij_interp(1, 2)
                kbar3_ij(ip, jp, i, j,k) = kbar3ij_interp(1, 2)
                kbar0_ij(jp, ip, i, j,k) = kbar0ij_interp(2, 1)
                kbar3_ij(jp, ip, i, j,k) = kbar3ij_interp(2, 1)
            else
                kbar0_ij(ip, jp, i, j,k) = 0.
                kbar3_ij(ip, jp, i, j,k) = 0.
                kbar0_ij(jp, ip, i, j,k) = 0.
                kbar3_ij(jp, ip, i, j,k) = 0.
            endif
            enddo
            enddo
            enddo

        endif
        enddo
        enddo
        kbar0_ij = max(kbar0_ij, 0.)
        kbar3_ij = max(kbar3_ij, 0.)
        kbar0_ij_cap = kbar0_ij
        kbar3_ij_cap = kbar3_ij
        !test kbar0 eq. kbar3
        kbar3_ij = kbar0_ij
        kbar3_ij_cap = kbar0_ij_cap

        call mpp_clock_end(coag_sub_clock2)

        !--------------------------------------------------------
        ! calculate number loss due to coagulation
        ! li(npop-1, :, :,:) loss of number for pop-i: #/m3/s
        ! note: dij is non-symmetric
        !--------------------------------------------------------
        call mpp_clock_begin(coag_sub_clock4)
        !if ((mpp_root_pe() .eq. mpp_pe()) .and. (dij_flag .eq. 0)) then
        !if (mpp_root_pe() .eq. mpp_pe()) then
        !        dij_flag = 1
        !        write(*,*) "coag_table: dij (loss of i to other population)"
        !        write(*,*) "i+j -> loss of number i"
        !        write(*,*) "for j=1,npop-1; for i=1,npop-1: bi(i)=+kbar*dij*nj"
        !        write(*,*) "for proof, dij need to be not symmetric"
        !        do i = 1, npop-1
        !                write(*,*) 'coag_run i=', i, ':', dij(i, :)
        !        end do
        ! end if
        bi = 0.
        do kp = 1,npop-1
        do ip = 1,npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(kp)%nb_tracer_pop > 0)) then
            num_pop_k = matrix_all_tracer(matrix_all_pop(kp)%i_n)%value_in_matrix(is:ie,js:je,:)
            if (dij(ip, kp) > 0) then !i.e. dij(ip, kp) = 1: dij means i+j will result the number loss of i into other pops, dij not
                !symmetric, and the order of indices matters!
                !if (mpp_root_pe().eq.mpp_pe()) then
                !    write(*,*) 'loss_num pair for ip: dij, ip, kp ', dij(ip, kp), ip, kp
                !endif
                bi(ip,:,:,:) = bi(ip,:,:,:) + kbar0_ij(ip,kp,:,:,:)*num_pop_k !coefficient loss of mode i due to coagulation
            endif
        endif
        enddo
        enddo
        call mpp_clock_end(coag_sub_clock4)

        call mpp_clock_begin(coag_sub_clock5)
        li=0.
        do ip =1, npop-1
        num_loss_max_rate = 0.
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            !if (mpp_root_pe().eq.mpp_pe()) then
            !    if (minval(num_pop_i) < 1e-32) then
            !        write(*,*) "pop ip =", ip, 'number min value less than 1e-32'
            !        num_pop_i = max(num_pop_i, 1e-32)
            !    endif
            !endif
            num_loss_max_rate = num_pop_i/tstep*0.5 !N is changing during the coagulation process, has to be smaller than ni induced rate 
            li(ip,:,:,:) = bi(ip,:,:,:) * num_pop_i + 0.5 * kbar0_ij(ip,ip,:,:,:) * num_pop_i * num_pop_i !loss of number due to coagulation #/m3/s
            do i=1,it
            do j=1,jt
            do k=1,kt
            if (li(ip, i,j,k) > num_loss_max_rate(i,j,k)) then
                n_scale = (num_loss_max_rate(i,j,k)/(li(ip, i,j,k)+tinydenon))
                kbar0_ij_cap(ip,ip,i,j,k) = kbar0_ij(ip,ip,i,j,k) * n_scale !if loss too much, scale
                do kp = 1,npop-1
                if (matrix_all_pop(kp)%nb_tracer_pop > 0) then
                    num_pop_k = matrix_all_tracer(matrix_all_pop(kp)%i_n)%value_in_matrix(is:ie,js:je,:)
                    !if (mpp_root_pe().eq.mpp_pe()) then
                    !    if (minval(num_pop_k) < 1e-32) then
                    !        write(*,*) "pop kp =", kp, 'number min value less than 1e-32'
                    !        num_pop_k = max(num_pop_k, 1e-32)
                    !    endif
                    !endif
                    if (dij(ip, kp) > 0) then
                        kbar0_ij_cap(ip,kp,i,j,k) = min(kbar0_ij_cap(ip,kp,i,j,k), &
                            kbar0_ij(ip,kp,i,j,k)*n_scale)
                        kbar0_ij_cap(ip,kp,i,j,k) = min(kbar0_ij_cap(ip,kp,i,j,k), kbar0_ij_cap(kp,ip,i,j,k)) !kbar0 must be symetric, and must be the minimum value
                        kbar0_ij_cap(kp,ip, i,j,k) = kbar0_ij_cap(ip,kp,i,j,k)
                    endif
                endif
                enddo
                !li(ip,i,j,k) = li(ip, i,j,k)*n_scale !li has to be recalculated outside of this subroutine
            endif
            enddo
            enddo
            enddo
        endif
        enddo
        !----------------------------------------------------
        !recalculate loss based on the updated kbar0_cap
        !----------------------------------------------------
        bi = 0.
        do kp = 1,npop-1
        do ip = 1,npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(kp)%nb_tracer_pop > 0)) then
            num_pop_k = matrix_all_tracer(matrix_all_pop(kp)%i_n)%value_in_matrix(is:ie,js:je,:)
            if (dij(ip, kp) > 0) then !i.e. dij(ip, kp) = 1: dij means i+j will result the number loss of i into other pops, dij not
                !symmetric, and the order of indices matters!
                bi(ip,:,:,:) = bi(ip,:,:,:) + kbar0_ij_cap(ip,kp,:,:,:)*num_pop_k !coefficient loss of mode i due to coagulation
            endif
        endif
        enddo
        enddo

        li=0.
        do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            li(ip,:,:,:) = bi(ip,:,:,:) * num_pop_i + 0.5 * kbar0_ij_cap(ip,ip,:,:,:) * num_pop_i * num_pop_i
        endif
        enddo
        call mpp_clock_end(coag_sub_clock5)

        !--------------------------------------------------------
        ! calculate number production due to coagulation
        ! ri(npop-1, :, :,:) production of number for pop-i: #/m3/s
        !-------------------------------------------------------- 
        call mpp_clock_begin(coag_sub_clock3)
        ri=0.
        do ikl = 1, ndikl !ikl is 1-12, exclude ext already
        ipop = dikl_control(ikl)%i
        kpop = dikl_control(ikl)%k
        lpop = dikl_control(ikl)%l
        if (matrix_all_pop(ipop)%nb_tracer_pop > 0) then
            if ((matrix_all_pop(kpop)%nb_tracer_pop > 0) .and. (matrix_all_pop(lpop)%nb_tracer_pop > 0)) then

                !if (mpp_root_pe().eq.mpp_pe()) then
                !    write(*,*) 'prod_num pair: n_ikl, i, k, l = ', ikl, ipop, kpop, lpop
                !endif

                num_pop_k = matrix_all_tracer(matrix_all_pop(kpop)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_l = matrix_all_tracer(matrix_all_pop(lpop)%i_n)%value_in_matrix(is:ie,js:je,:)
                !if (mpp_root_pe().eq.mpp_pe()) then
                !    if (minval(num_pop_k) < 1e-32) then
                !        write(*,*) "pop kp =", kpop, 'number min value less than 1e-32'
                !        num_pop_k = max(num_pop_k, 1e-32)
                !    endif
                !    if (minval(num_pop_l) < 1e-32) then
                !        write(*,*) "pop lpop = ", lpop, 'number min value less than 1e-32'
                !        num_pop_l = max(num_pop_l, 1e-32)
                !    endif
                !endif
                ri(ipop,:,:,:) = ri(ipop,:,:,:) + kbar0_ij_cap(kpop,lpop,:,:,:) * num_pop_k * num_pop_l !production of number due to coagulation
            endif
        endif
        enddo
        call mpp_clock_end(coag_sub_clock3)

        !if kbar3=kbar0
        kbar3_ij_cap = kbar0_ij_cap
        !--------------------------------------------------------
        ! calculate mass loss due to coagulation
        ! lim(npop-1,nmspcs, :, :,:) loss of mass for pop-i, spec-qq :ug/m3/s
        !--------------------------------------------------------
        call mpp_clock_begin(coag_sub_clock7)
        fi = 0.
        do ip = 1,npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            do jp = 1,npop-1
            if (matrix_all_pop(jp)%nb_tracer_pop > 0) then
                num_pop_j =  matrix_all_tracer(matrix_all_pop(jp)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_j = max(num_pop_j, 1e-32)
                if( dij(ip,jp) > 0 ) then
                    !if (mpp_root_pe().eq.mpp_pe()) then
                    !    write(*,*) 'loss_mass pair for dij: ip, jp ', ip, jp
                    !endif
                    fi(ip,:,:,:) = fi(ip,:,:,:) + kbar3_ij(ip,jp,:,:,:)*num_pop_j

                endif
            endif
            enddo
        endif
        enddo

        lim = 0.
        do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            do q=1, nm(ip)
            qq = prod_index(ip,q)
            num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            num_pop_i = max(num_pop_i, 1e-32)
            lim(ip,qq,:,:,:) = lim(ip,qq,:,:,:) + fi(ip,:,:,:) * num_pop_i * mjq(ip,qq,:,:,:)
            do i=1,it
            do j=1,jt
            do k=1,kt 
            if (lim(ip,qq, i,j,k) > mass_loss_max_rate(ip,qq,i,j,k)) then
                m_scale = ( mass_loss_max_rate(ip,qq,i,j,k)/(lim(ip,qq,i,j,k)+tinydenom))
                do jp = 1,npop-1
                if (matrix_all_pop(jp)%nb_tracer_pop > 0) then
                    if( dij(ip,jp) > 0 ) then
                        kbar3_ij_cap(ip,jp,i,j,k) = kbar3_ij(ip,jp,i,j,k)* m_scale
                        !( mass_loss_max_rate(ip,qq,i,j,k)/(lim(ip,qq,i,j,k)+tinydenom)) 
                    endif
                endif
                enddo
                lim(ip,qq,i,j,k) = lim(ip,qq, i,j,k)* m_scale
            endif
            enddo
            enddo
            enddo
            !if (mpp_root_pe().eq.mpp_pe()) then
            !    write(*,*) 'loss_mass pair: ip, q, qq, i ', ip, q, qq, i
            !endif
            enddo
        endif
        enddo
        !used = send_data(id_lim_11, lim(1,1,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_lim_21, lim(2,1,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        call mpp_clock_end(coag_sub_clock7)

        !--------------------------------------------------------
        ! calculate mass production due to coagulation
        ! rim(npop-1,nmspcs, :, :,:) production of mass for pop-i, spec-qq :ug/m3/s
        !--------------------------------------------------------
        call mpp_clock_begin(coag_sub_clock6)
        rim = 0.
        do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0 .and. giklq_control(ip)%n > 0) then
            do klq = 1, giklq_control(ip)%n
            k = giklq_control(ip)%k(klq)
            l = giklq_control(ip)%l(klq)
            qq = giklq_control(ip)%qq(klq)
            if (matrix_all_pop(k)%nb_tracer_pop > 0 .and. matrix_all_pop(l)%nb_tracer_pop > 0 ) then

                !if (mpp_root_pe().eq.mpp_pe()) then
                !    write(*,*) 'prod_mass pair for ip: ip, klq, k, l, qq ', ip, klq, k, l, qq
                !endif
                num_pop_k = matrix_all_tracer(matrix_all_pop(k)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_k = max(num_pop_k, 1e-32)
                num_pop_l = matrix_all_tracer(matrix_all_pop(l)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_l = max(num_pop_l, 1e-32)
                rim(ip,qq,:,:,:) = rim(ip,qq,:,:,:) + num_pop_k*num_pop_l*kbar3_ij_cap(l,k,:,:,:)*(mjq(l,qq,:,:,:))!%+mjq(k,qq,:,:,:))
            endif
            enddo
        endif
        enddo
        call mpp_clock_end(coag_sub_clock6)

        !if (mpp_root_pe().eq.mpp_pe()) then
        !    write(*,*) "maximum ratio of kbar0/kbar0_cap to kbar3/kbar3_cap:, and max_index is: ",  &
        !        maxval(max(kbar0_ij,1e-30)/max(kbar0_ij_cap,1e-30)), &
        !        maxval(max(kbar3_ij,1e-30)/max(kbar3_ij_cap,1e-30))
        !endif
    end subroutine matrix_coag_cap        




    subroutine matrix_coag(pfull, t, ri, bi, li, rim, fi, lim, kbar0_ij, mjq, is,ie,js,je, time, time_next)
        implicit none
        logical :: used
        type(time_type),  intent(in) :: time, time_next
        integer, intent(in) :: is,ie,js,je
        real, intent(in) :: pfull(:,:,:), t(:,:,:)
        real, intent(out) :: bi(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !loss terms due to intermodal coagulation. [1/s]
        real, intent(out) :: li(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !loss term of number due to intermodal coagulation [#/m^3/s]
        real, intent(out) :: ri(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !production terms due to intermodal coagulation. [#/m^3/s]
        real, intent(out) :: fi(npop-1,size(pfull,1),size(pfull,2),size(pfull,3)) !mass loss coefficient
        real, intent(out) :: lim(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass loss term of number due to intermodal coagulation [ug/m^3/s]
        real, intent(out) :: rim(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass production terms due to intermodal coagulation. [ug/m^3/s]
        real, intent(out) :: kbar0_ij(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3)) !mode averaged coagulation coefficient [m3/s]
        real :: kbar3_ij(npop-1,npop-1, size(pfull,1), size(pfull,2), size(pfull,3)) !mode averaged coagulation coefficient [m3/s]
        real :: dg_ip_um(size(pfull,1),size(pfull,2),size(pfull,3)), dg_jp_um(size(pfull,1),size(pfull,2),size(pfull,3)) !particle geometric diameter [um]
        real :: num_pop_k(size(pfull,1),size(pfull,2),size(pfull,3)), num_pop_l(size(pfull,1),size(pfull,2),size(pfull,3)), num_pop_i(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: num_pop_j(size(pfull,1),size(pfull,2),size(pfull,3))
        real, intent(out) :: mjq(npop-1,nmspcs,size(pfull,1),size(pfull,2),size(pfull,3)) !mass of species in a single particle: mass/number = ug/particle
        real :: sig_ip, sig_jp, pop_numi(size(pfull,1),size(pfull,2),size(pfull,3)), pop_numj(size(pfull,1),size(pfull,2),size(pfull,3))
        integer :: ip,jp,kp,it,jt,kt,i,j,k,l,ikl,ipop,kpop,lpop,q,qq,klq
        real :: sig_2pop(2), diam_2pop(2), kbar0ij_interp(2,2), kbar3ij_interp(2,2)
        !sig_2pop, diam_2pop, kbar0ij_interp, kbar3ij_interp
        !real, intent(out) :: 
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        sig_2pop = 1.
        diam_2pop = 1e-30
        kbar0_ij = 0
        kbar3_ij = 0
        kbar0ij_interp = 0
        kbar3ij_interp = 0
        bi = 0.
        ri = 0. !production of number due to coagulation
        li = 0.
        fi = 0.
        lim = 0.
        rim = 0.
        mjq = 0. !ug/particle
        !constaruct mjq
        !mjq: mjq(j,q) is the avg. mass/particle of species q (=1-5) for mode j
        call mpp_clock_begin(coag_sub_clock1)
        do ip=1,npop-1 !exclude 'ext' pop
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_pop_i = matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            num_pop_i = max(num_pop_i, 1e-32)
            do q = 1, nm(ip)
            qq=prod_index(ip,q) !index of mass species: 1-5
            !if (mpp_root_pe().eq.mpp_pe()) then
            !    write(*,*) 'ip, nm(ip), q, qq', ip, nm(ip), q, qq
            !endif
            select case (qq)
            case (1)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_msulf)%value_in_matrix(is:ie,js:je,:)/num_pop_i
            case (4)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mdust)%value_in_matrix(is:ie,js:je,:)/num_pop_i
            case (5)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mseas)%value_in_matrix(is:ie,js:je,:)/num_pop_i
            case (3)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(is:ie,js:je,:)/num_pop_i
            case (2)
                mjq(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%i_mbcar)%value_in_matrix(is:ie,js:je,:)/num_pop_i
            case default !call fatal error
                call error_mesg('matrix coagulation','mjq mass species not properly defined ', fatal)
            end select
            enddo
        endif
        enddo
        mjq = max(0.0, mjq)
        call mpp_clock_end(coag_sub_clock1)

        used = send_data(id_mjq11, mjq(1,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1) !akk: mass/num
        used = send_data(id_mjq21, mjq(2,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1) !acc: mass/num

        call mpp_clock_begin(coag_sub_clock2)
        !set up kbar0_ij and kbar3_ij table (npop, npop, it, jt, kt)
        do ip = 1, npop-1 !exclude 'ext' pop
        do jp = ip, npop-1
        !do jp = 1, npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(jp)%nb_tracer_pop > 0)) then
            dg_ip_um = matrix_all_pop(ip)%dg_wet(is:ie,js:je,:)*1.0e6 !3d array [um]
            sig_ip = matrix_all_pop(ip)%sigma !real number
            dg_jp_um = matrix_all_pop(jp)%dg_wet(is:ie,js:je,:)*1.0e6 ! [um]
            sig_jp = matrix_all_pop(jp)%sigma
            pop_numi = matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            pop_numj = matrix_all_tracer(matrix_all_pop(jp)%i_n)%value_in_matrix(is:ie,js:je,:) 
            do i=1,it
            do j=1,jt
            do k=1,kt
            if ((pop_numi(i,j,k) > 1e-17) .and. (pop_numj(i,j,k) > 1e-17) .and. &
                (dg_ip_um(i,j,k) > 0) .and. (dg_jp_um(i,j,k) > 0)) then
                sig_2pop=(/sig_ip, sig_jp/)
                diam_2pop = (/dg_ip_um(i,j,k), dg_jp_um(i,j,k)/)
                call get_kbarnij(1, t(i,j,k), pfull(i,j,k), 2, sig_2pop, diam_2pop, kbar0ij_interp, kbar3ij_interp) 
                ! call get_knij(t(i,j,k), pfull(i,j,k), dg_ip_um(i,j,k), sig_ip, dg_jp_um(i,j,k), sig_jp, kbar0_ij(ip, jp, i,j,k), &
                !    kbar3_ij(ip, jp, i,j,k)) !kbar0_ij/kbar3_ij mode averaged coagulation coefficient, [m3/s]
                !if (mpp_root_pe().eq.mpp_pe()) then
                !if ( any(kbar0ij_interp > 1e-7) .or. any(kbar3ij_interp > 1e-7)) then
                !        write(*,*) 'kbar0ij_interp=', kbar0ij_interp(1, 1), kbar0ij_interp(1, 2), kbar0ij_interp(2, 1), &
                !        kbar0ij_interp(2, 2), 't=', t(i,j,k), 'pfull=', pfull(i,j,k), 'dg_ip_um=', dg_ip_um(i,j,k), &
                !        'dg_jp_um = ', dg_jp_um(i,j,k), 'sigmai=', sig_ip, 'sig_jp=', sig_jp, 'pops_ip, jp=', &
                !        matrix_all_pop(ip)%name, matrix_all_pop(jp)%name, &
                !        'num_conci, number_concj=', matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(i+is-1,j+js-1,k), &
                !        matrix_all_tracer(matrix_all_pop(jp)%i_n)%value_in_matrix(i+is-1,j+js-1,k), &
                !        'mass_conci, mass_concj=', matrix_all_tracer(matrix_all_pop(ip)%i_mocar)%value_in_matrix(i+is-1,j+js-1,k), &
                !        matrix_all_tracer(matrix_all_pop(jp)%i_mocar)%value_in_matrix(i+is-1,j+js-1,k)
                !endif
                !endif
                kbar0_ij(ip, jp, i, j,k) = kbar0ij_interp(1, 2)
                kbar3_ij(ip, jp, i, j,k) = kbar3ij_interp(1, 2)
                kbar0_ij(jp, ip, i, j,k) = kbar0ij_interp(2, 1)
                kbar3_ij(jp, ip, i, j,k) = kbar3ij_interp(2, 1)
            else
                kbar0_ij(ip, jp, i, j,k) = 0.
                kbar3_ij(ip, jp, i, j,k) = 0.
                kbar0_ij(jp, ip, i, j,k) = 0.
                kbar3_ij(jp, ip, i, j,k) = 0.
            endif
            enddo
            enddo
            enddo
            !if (mpp_root_pe().eq.mpp_pe()) then !xl debug

            !   if ((maxval(kbar0_ij(ip, jp,:,:,:))>1e-7) .or. (maxval(kbar3_ij(ip, jp,:,:,:))>1e-7)) then
            !       write(*,*) 'ip, jp=', matrix_all_pop(ip)%name, matrix_all_pop(jp)%name, &
            !           'max_kbar0=', maxval(kbar0_ij(ip, jp,:,:,:)), & 
            !           'max_kbar3=', maxval(kbar3_ij(ip, jp,:,:,:)), &  
            !           'min_kbar0=', minval(kbar0_ij(ip, jp,:,:,:)), &
            !           'min_kbar3=', minval(kbar3_ij(ip, jp,:,:,:)), &
            !           'max_dgi, min_dgi, max_dgj, min_dgj=', maxval(dg_ip_um), minval(dg_ip_um), maxval(dg_jp_um), minval(dg_jp_um), &
            !           'minp, maxp, mint, maxt=', minval(pfull), maxval(pfull), minval(t), maxval(t)
            !   endif
            !endif
        endif
        enddo
        enddo
        kbar0_ij = max(kbar0_ij, 0.)
        kbar3_ij = max(kbar3_ij, 0.)
        !if (mpp_root_pe().eq.mpp_pe()) then
        !    write(*,*) "maximum ratio of kbar0 to kbar3:, and max_index is: ",  maxval(max(kbar3_ij, 1e-30)/max(kbar0_ij, 1e-30)), &
        !        maxloc(max(kbar3_ij, 1e-30)/max(kbar0_ij, 1e-30))
        !endif
        kbar3_ij = kbar0_ij
        call mpp_clock_end(coag_sub_clock2)

        !used = send_data(id_kbar0_11, kbar0_ij(1,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_kbar3_11, kbar3_ij(1,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_kbar0_12, kbar0_ij(1,2,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_kbar3_12, kbar3_ij(1,2,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_kbar0_21, kbar0_ij(2,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_kbar3_21, kbar3_ij(2,1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_kbar0_22, kbar0_ij(2,2,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_kbar3_22, kbar3_ij(2,2,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)

        !--------------------------------------------------------
        ! calculate number production due to coagulation
        ! ri(npop-1, :, :,:) production of number for pop-i: #/m3/s
        !--------------------------------------------------------
        !if (mpp_root_pe().eq.mpp_pe()) then !xl debug 
        !        write(*,*) 'pass_test'
        !        write(*,*) 'ndikl', ndikl
        !endif
        call mpp_clock_begin(coag_sub_clock3)
        ri=0.
        do ikl = 1, ndikl !ikl is 1-12, exclude ext already
        ipop = dikl_control(ikl)%i
        kpop = dikl_control(ikl)%k
        lpop = dikl_control(ikl)%l
        if (matrix_all_pop(ipop)%nb_tracer_pop > 0) then
            if ((matrix_all_pop(kpop)%nb_tracer_pop > 0) .and. (matrix_all_pop(lpop)%nb_tracer_pop > 0)) then

                !if (mpp_root_pe().eq.mpp_pe()) then
                !    write(*,*) 'prod_num pair: n_ikl, i, k, l = ', ikl, ipop, kpop, lpop
                !endif

                num_pop_k = matrix_all_tracer(matrix_all_pop(kpop)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_l = matrix_all_tracer(matrix_all_pop(lpop)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_k = max(num_pop_k, 1e-32)
                num_pop_l = max(num_pop_l, 1e-32)
                ri(ipop,:,:,:) = ri(ipop,:,:,:) + kbar0_ij(kpop,lpop,:,:,:) * num_pop_k * num_pop_l !production of number due to coagulation
            endif
        endif
        enddo
        call mpp_clock_end(coag_sub_clock3)

        !used = send_data(id_ri_1, ri(1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_ri_2, ri(2,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)

        !--------------------------------------------------------
        ! calculate number loss due to coagulation
        ! li(npop-1, :, :,:) loss of number for pop-i: #/m3/s
        ! note: dij is non-symmetric
        !--------------------------------------------------------
        call mpp_clock_begin(coag_sub_clock4)
        bi=0.
        do kp = 1,npop-1
        do ip = 1,npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(kp)%nb_tracer_pop > 0)) then
            num_pop_k = matrix_all_tracer(matrix_all_pop(kp)%i_n)%value_in_matrix(is:ie,js:je,:)
            num_pop_k = max(num_pop_k, 1e-32)
            if (dij(ip, kp) > 0) then !i.e. dij(ip, kp) = 1

                !if (mpp_root_pe().eq.mpp_pe()) then
                !    write(*,*) 'loss_num pair for ip: dij, ip, kp ', dij(ip, kp), ip, kp
                !endif

                bi(ip,:,:,:) = bi(ip,:,:,:) + kbar0_ij(ip,kp,:,:,:)*num_pop_k !coefficient loss of mode i due to coagulation
            endif
        endif
        enddo
        enddo
        call mpp_clock_end(coag_sub_clock4)
        !used = send_data(id_bi_1, bi(1,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_bi_2, bi(2,:,:,:), time_next, is_in=is,js_in=js, ks_in=1)

        call mpp_clock_begin(coag_sub_clock5)
        li=0.
        do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            num_pop_i = max(num_pop_i, 1e-32)
            li(ip,:,:,:) = bi(ip,:,:,:) * num_pop_i + 0.5 * kbar0_ij(ip,ip,:,:,:) * num_pop_i * num_pop_i !loss of number due to coagulation #/m3/s
        endif
        enddo
        call mpp_clock_end(coag_sub_clock5)
        !ip = 1 
        !num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
        !num_pop_i = max(num_pop_i, 1e-32)
        !used = send_data(id_li_1, li(1,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_half_knn_1, 0.5 * kbar0_ij(1,1,:, :,:)*num_pop_i * num_pop_i, time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_notself_knn_1, bi(1,:,:,:) * num_pop_i, time, is_in=is,js_in=js, ks_in=1)

        !ip = 2 
        !num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
        !num_pop_i = max(num_pop_i, 1e-32)
        !used = send_data(id_li_2, li(2,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_half_knn_2, 0.5 * kbar0_ij(2,2,:, :,:)*num_pop_i * num_pop_i, time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_notself_knn_2, bi(2,:,:,:) * num_pop_i, time, is_in=is,js_in=js, ks_in=1)        



        !--------------------------------------------------------
        ! calculate mass production due to coagulation
        ! rim(npop-1,nmspcs, :, :,:) production of mass for pop-i, spec-qq :ug/m3/s
        !--------------------------------------------------------
        call mpp_clock_begin(coag_sub_clock6)
        rim = 0.
        do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0 .and. giklq_control(ip)%n > 0) then
            do klq = 1, giklq_control(ip)%n
            k = giklq_control(ip)%k(klq)
            l = giklq_control(ip)%l(klq)
            qq = giklq_control(ip)%qq(klq)
            if (matrix_all_pop(k)%nb_tracer_pop > 0 .and. matrix_all_pop(l)%nb_tracer_pop > 0 ) then

                !if (mpp_root_pe().eq.mpp_pe()) then
                !    write(*,*) 'prod_mass pair for ip: ip, klq, k, l, qq ', ip, klq, k, l, qq
                !endif
                num_pop_k = matrix_all_tracer(matrix_all_pop(k)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_k = max(num_pop_k, 1e-32)
                num_pop_l = matrix_all_tracer(matrix_all_pop(l)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_l = max(num_pop_l, 1e-32)
                rim(ip,qq,:,:,:) = rim(ip,qq,:,:,:) + num_pop_k*num_pop_l*kbar3_ij(l,k,:,:,:)*(mjq(l,qq,:,:,:))!%+mjq(k,qq,:,:,:))
            endif
            enddo
        endif
        enddo
        call mpp_clock_end(coag_sub_clock6)
        !used = send_data(id_rim_11, rim(1,1,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_rim_21,  rim(2,1,:,:,:), time, is_in=is,js_in=js, ks_in=1)


        !--------------------------------------------------------
        ! calculate mass loss due to coagulation
        ! lim(npop-1,nmspcs, :, :,:) loss of mass for pop-i, spec-qq :ug/m3/s
        !--------------------------------------------------------
        call mpp_clock_begin(coag_sub_clock7)
        fi = 0.
        do ip = 1,npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then 
            do jp = 1,npop-1
            if (matrix_all_pop(jp)%nb_tracer_pop > 0) then
                num_pop_j =  matrix_all_tracer(matrix_all_pop(jp)%i_n)%value_in_matrix(is:ie,js:je,:)
                num_pop_j = max(num_pop_j, 1e-32)
                if( dij(ip,jp) > 0 ) then
                    !if (mpp_root_pe().eq.mpp_pe()) then
                    !    write(*,*) 'loss_mass pair for dij: ip, jp ', ip, jp
                    !endif
                    fi(ip,:,:,:) = fi(ip,:,:,:) + kbar3_ij(ip,jp,:,:,:)*num_pop_j

                endif
            endif
            enddo
        endif
        enddo

        !used = send_data(id_fi_1, fi(1,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_fi_2, fi(2,:,:,:), time, is_in=is,js_in=js, ks_in=1)

        lim = 0.
        do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
            do q=1, nm(ip)
            qq = prod_index(ip,q)
            num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%i_n)%value_in_matrix(is:ie,js:je,:)
            num_pop_i = max(num_pop_i, 1e-32)
            lim(ip,qq,:,:,:) = lim(ip,qq,:,:,:) + fi(ip,:,:,:) * num_pop_i * mjq(ip,qq,:,:,:)

            !if (mpp_root_pe().eq.mpp_pe()) then
            !    write(*,*) 'loss_mass pair: ip, q, qq, i ', ip, q, qq, i
            !endif
            enddo
        endif
        enddo
        !used = send_data(id_lim_11, lim(1,1,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        !used = send_data(id_lim_21, lim(2,1,:,:,:), time, is_in=is,js_in=js, ks_in=1)
        call mpp_clock_end(coag_sub_clock7)
    end subroutine


    !----------------------------------------------------------------
    !              subroutine matrix_dry_diameter
    !       ...calculate 
    !          (1) dry density: kg/m3 (2) dry diamter: m
    !          (3) dry mass: ug/m3
    !          (4) kappa_pop: volume average of kappa for each species
    !
    !       !calc average density - sum(m_spec,i)/ sum(volume_spec,i)
    !        !calculate average hygroscopicity: (kappa* (volume_spec,i))/sum(volume_spec,i)
    !        ! 1/6*pi*d_dry^3 = sum(volume_spec,i)       
    !----------------------------------------------------------------
    subroutine matrix_dry_diameter(pfull,is,ie,js,je)
        real, intent(in) :: pfull(:,:,:)
        integer :: n, nt, ntt, tr_index !local variables
        real :: m_dry_all(size(pfull,1),size(pfull,2),size(pfull,3)), vol_dry_all(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: kappa_vol_all(size(pfull,1),size(pfull,2),size(pfull,3)), vol_spec(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: pop_num(size(pfull,1),size(pfull,2),size(pfull,3)), dens_dry(size(pfull,1),size(pfull,2),size(pfull,3))
        real :: exp_fac
        integer, intent(in) :: is,ie,js,je
        integer :: flag,it,jt,kt, i, j, k !flag for whether this population has dry mass/perform calculation
        it = size(pfull,1)
        jt = size(pfull,2)
        kt = size(pfull,3)
        do n = 1, npop
        m_dry_all = 0.
        vol_dry_all = 0.
        kappa_vol_all = 0.
        flag = 0
        m_dry_all = 0.
        vol_dry_all = 0.
        kappa_vol_all = 0.
        vol_spec = 0.
        pop_num = 1.e-17
        dens_dry = 1000.
        exp_fac = 0
        if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext"))  then
            nt = matrix_all_pop(n)%nb_tracer_pop
            do ntt = 1, nt
            tr_index = matrix_all_pop(n)%tracer_index(ntt)
            if ((lowercase(trim(matrix_all_tracer(tr_index)%spec)) .ne. "alwc") .and. &
                (lowercase(trim(matrix_all_tracer(tr_index)%type)) .ne. "number") ) then
                m_dry_all = m_dry_all + matrix_all_tracer(tr_index)%value_in_matrix(is:ie,js:je,:) !m_dry unit: ug/m3
                !vol_spec unit: m3 aerosol_volume/ m3 air
                vol_spec = (matrix_all_tracer(tr_index)%value_in_matrix(is:ie,js:je,:))/(matrix_all_tracer(tr_index)%dens)*1e-9 
                vol_dry_all = vol_dry_all + vol_spec
                kappa_vol_all = kappa_vol_all + matrix_all_tracer(tr_index)%kappa*vol_spec !sum of kappa*vol
                flag = 1 !record if this population has dry mass
            endif
            enddo
            if (flag > 0) then
                matrix_all_pop(n)%mass_dry(is:ie,js:je,:) = m_dry_all !mass unit: ug/m3
                matrix_all_pop(n)%vol_dry(is:ie,js:je,:) = vol_dry_all !volume unit: m3 aerosol/m3
                pop_num = matrix_all_tracer(matrix_all_pop(n)%i_n)%value_in_matrix(is:ie,js:je,:) !#/m3
                exp_fac = exp(1.5*(log(matrix_all_pop(n)%sigma))**2) 
                matrix_all_pop(n)%dens_dry(is:ie, js:je, :) = m_dry_all/vol_dry_all*1e-9
                matrix_all_pop(n)%kappa_pop(is:ie, js:je, :) = kappa_vol_all/vol_dry_all 

                do i=1,it
                do j=1,jt
                do k=1,kt
                if ((pop_num(i,j,k) > 1e-17) .and. (m_dry_all(i,j,k) > 1e-32)) then !not zero number on the grid
                    matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k) = &
                        max(3.0e-9, (m_dry_all(i,j,k)/pop_num(i,j,k)/dens_dry(i,j,k)*1e-9/pi6)**(1.0/3)/exp_fac ) ! dg_dry unit: m
                    !if (matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k) > 1.0e-5) then
                    !    write(*,*) 'steps, dg_dry > 10 um: dg_dry, dry_mass, num, dens_dry, exp_fac', &
                    !        step, matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k)*1.0e6, m_dry_all(i,j,k), pop_num(i,j,k), &
                    !        dens_dry(i,j,k), exp_fac
                    !endif
                    matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k) = min(3.0e-5, matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k))
                else
                    dens_dry(i,j,k) = 1000
                    matrix_all_pop(n)%dens_dry(i+is-1,j+js-1,k) = matrix_all_pop(n)%def_dens
                    matrix_all_pop(n)%kappa_pop(i+is-1,j+js-1,k)= matrix_all_pop(n)%def_kappa
                    matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k) = matrix_all_pop(n)%def_dg
                endif
                enddo
                enddo
                enddo
            endif
        endif 
        enddo
    end subroutine



    subroutine matrix_wet_diameter( rh, t, is,ie,js,je) !calculate and assgin values for pop%dg_wet
        real, intent(in) :: rh(:,:,:), t(:,:,:)
        real:: ddry_in_3d(size(rh,1),size(rh,2),size(rh,3)), hygro_3d(size(rh,1),size(rh,2),size(rh,3))
        integer :: n, nt, ntt, i, j, k, it, jt, kt, tr_index
        integer,intent(in) :: is,ie,js,je
        real :: ddry_in, hygro_in, s_in, tair_in, dwet_out, gf, rh_deliquescence, rh_crystallization
        real :: particle_vol_dry, particle_vol_water, particle_vol_wet, f_hysteresis, gf3 
        integer :: flag != 0 ! check if the population get calculated 
        flag = 0
        it = size(rh,1)
        jt = size(rh,2)
        kt = size(rh,3)
        ddry_in_3d = 0.
        hygro_3d = 0.
        dwet_out = 0.
        do n = 1, npop
        flag = 0
        call mpp_clock_begin(hygrow_sub_clock1)
        if (matrix_all_pop(n)%nb_tracer_pop > 0) then
            nt = matrix_all_pop(n)%nb_tracer_pop
            do ntt = 1, nt
            tr_index = matrix_all_pop(n)%tracer_index(ntt)
            if ((lowercase(trim(matrix_all_tracer(tr_index)%spec)) .ne. "alwc") .and. &
                (lowercase(trim(matrix_all_tracer(tr_index)%type)) .ne. "number") ) then
                flag = 1
            endif
            enddo
        endif
        call mpp_clock_end(hygrow_sub_clock1)
        call mpp_clock_begin(hygrow_sub_clock2)
        if (flag > 0) then
            ddry_in_3d = matrix_all_pop(n)%dg_dry(is:ie,js:je,:)*exp(1.5*(log(matrix_all_pop(n)%sigma))**2) !volume mean diameter
            ddry_in_3d = max(ddry_in_3d, 3.0e-9)
            hygro_3d = matrix_all_pop(n)%kappa_pop(is:ie,js:je,:)
            rh_deliquescence = matrix_all_pop(n)%rh_deliquescence
            rh_crystallization = matrix_all_pop(n)%rh_crystallization
            call mpp_clock_begin(hygrow_sub_clock3)
            do i=1,it
            do j=1,jt
            do k=1,kt
            ddry_in = ddry_in_3d(i,j,k)
            particle_vol_dry = pi6*ddry_in**3
            hygro_in = hygro_3d(i,j,k)
            s_in = max(0.0, rh(i,j,k))
            s_in = min(rh_cap, s_in) !cap the rh values for hygroscopic growth, default is 0.97
            tair_in = t(i,j,k)
            gf = 1.0 !default value
            if (ddry_in > 10.0e-9) then !only consider hygroscopic growth larger than 10 nm
                call mpp_clock_begin(hygrow_sub_clock4)
                call aero_kohler(ddry_in, hygro_in, s_in, tair_in, dwet_out)
                call mpp_clock_end(hygrow_sub_clock4)
                dwet_out = max(ddry_in, dwet_out) !dwet_out and ddry_in are volume mean diameter
                particle_vol_water = pi6*(dwet_out**3 - ddry_in**3)
                !check current rh with pop rh crystalization and deliquescence
                if (s_in < rh_crystallization) then
                    dwet_out = ddry_in
                    particle_vol_water = 0
                elseif (s_in < rh_deliquescence) then
                    f_hysteresis = 1.0 / max(1.0e-5, (rh_deliquescence - rh_crystallization))
                    particle_vol_water = f_hysteresis * (s_in - rh_crystallization) * particle_vol_water
                    particle_vol_water = max(0.0, particle_vol_water)
                    particle_vol_wet = particle_vol_dry + particle_vol_water
                    dwet_out = (particle_vol_wet / pi6)**(1.0/3)
                end if
                gf=max(dwet_out/ddry_in,1.0) !growth factor
            endif
            matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k) = gf*matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k)

            !if (matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k) > 1.0e-5) then
            !                write(*,*) 'dg_wet > 10 um: dg_dry, dg_wet, gf, hygro_in, s', &
            !                        matrix_all_pop(n)%dg_dry(i+is-1,j+js-1,k)*1.0e6, &
            !                        matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k)*1.0e6, gf, &
            !                        hygro_in, s_in
            !endif
            !matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k) = min(50.0e-6, matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k))
            !matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k) = max(3.0e-9, matrix_all_pop(n)%dg_wet(i+is-1,j+js-1,k))
            gf3 = gf**3
            matrix_all_pop(n)%dens_wet(i+is-1,j+js-1,k) = (1000*(gf3-1) + matrix_all_pop(n)%dens_dry(i+is-1,j+js-1,k))/gf3 !wet density
            enddo
            enddo
            enddo
            call mpp_clock_end(hygrow_sub_clock3)
        endif
        call mpp_clock_end(hygrow_sub_clock2)
        enddo
    end subroutine

    subroutine module_number_solver(ipop, ri_ipop, bi_ipop, kbar0_ij_iipop, is, ie, js, je, dt)
        integer, intent(in) :: ipop, is, ie, js, je
        real, intent(in) :: dt
        real, intent(in) :: kbar0_ij_iipop(:,:,:)
        real, intent(in) :: ri_ipop(:,:,:), bi_ipop(:,:,:) !ri: production from coagulation, bi_ipop: loss coefficient from coagulation
        real :: ci(size(ri_ipop,1),size(ri_ipop,2),size(ri_ipop,3)), bi(size(ri_ipop,1),size(ri_ipop,2),size(ri_ipop,3))
        real ::  kbar0ij_ii(size(ri_ipop,1),size(ri_ipop,2),size(ri_ipop,3))
        integer :: i, j, k, it, jt, kt
        it = size(ri_ipop,1)
        jt = size(ri_ipop,2)
        kt = size(ri_ipop,3)
        !dn/dt = (ri+source) -bi*n -0.5*kbar0*n^2
        ci = matrix_all_tracer(matrix_all_pop(ipop)%i_n)%source(is:ie, js:je, :) + ri_ipop 
        bi =  bi_ipop
        kbar0ij_ii = kbar0_ij_iipop
        if (any(kbar0ij_ii < 0)) then
            write(*,*) 'coagulatio coefficient negative exitst for popi=', ipop
        endif
        do i=1,it
        do j=1,jt
        do k=1,kt
        !write(*,*) 'ipop, number, kabrii, bi, ci, dt', ipop, &
        !        matrix_all_tracer(matrix_all_pop(ipop)%i_n)%value_in_matrix(i+is-1, j+js-1, k), kbar0ij_ii(i,j,k), bi(i,j,k), ci(i,j,k), dt
        call number_solver(ipop, matrix_all_tracer(matrix_all_pop(ipop)%i_n)%value_in_matrix(i+is-1, j+js-1, k), &
            kbar0ij_ii(i,j,k), bi(i,j,k), ci(i,j,k), dt)
        enddo
        enddo
        enddo
    end subroutine
    !---------------------------------------------------------------
    ! subroutine number_solver
    ! for number tracers: dni/dt = c - b*ni -a* ni^2
    ! c: total source, from emision + npf + coagulation generation
    ! b*ni: intermodal coagulation loss
    ! a*ni^2: self-coagulation loss
    !---------------------------------------------------------------
    subroutine number_solver(ipop, ni, kbar0ij_ii, bi, ci, tstep)
        integer, intent(in) :: ipop
        real, intent(inout) :: ni
        real, intent(in) :: kbar0ij_ii, bi, ci, tstep
        real :: y0, a, b, c, delta, r1, r2, gamma_i, gexpdt, y, piq_thresh, minconc, expdt
        minconc = 1.0d-17 ! [ug/m^3] and [#/m^3]
        piq_thresh = 1.0d-08  ! [1] threshold in number/mass conc. solver
        y0 = ni
        a = 0.5d+00*kbar0ij_ii
        b = bi
        c = ci
        if( c .gt. 1.0d-20 ) then
            if ((b * b + 4.0d+00 * a * c) < 0) then
                write(*,*) "ipop, ni, kbar0ij_ii, bi, ci, tstep", ipop, ni, kbar0ij_ii, bi, ci, tstep
            endif
            delta = sqrt( b * b + 4.0d+00 * a * c )
            if ((a .gt. 1.0d-30)) then
                r1 = 2.0d+00 * a * c / ( b + delta )
                r2 = - 0.5d0 * ( b + delta )
                gamma_i =  - ( r1 - a * y0 ) / ( r2 - a * y0 )
                gexpdt = gamma_i * exp( - delta * tstep )
                !xl tets --------------------------------
                if (a * ( 1.0d+00 + gexpdt ) < 1.0e-32) then
                    write(*,*) "number_solver note gexpdt, r1, r2, gamma_i, ni, a, b, c", gexpdt, r1, r2, gamma_i, ni, a, b, c
                endif
                !gexpt=gamma_i*exp, the exp can be 1, gammai can be -1, which causes  a * ( 1.0d+00 + gexpdt )=0, then explode
                !e.g, one output i got is: 
                ! -1.00000000000000       1.569014248053667e-022 -1.569014377459124e-022
                ! -1.00000000000000        1919339290.35171       1.193395344430849e-015
                ! 1.294054570989738e-029  2.062858653775375e-029
                !to resolve this problem: gexpdt must has a range condition, e.g. using piq_rhresh
                !gamma_i =-1 commonly happens as long as an0 is orders larger than b or c
                !hence, to make sure the solution valid, we only need exp(-delta*tstep) not equals 1, delta=sqrt( b * b + 4.0d+00 * a * c )
                !----------------------------------------------
                if( 1.0d+00-exp( - delta * tstep ) .gt. piq_thresh ) then
                    y = (  r1 + r2 * gexpdt ) / ( a * ( 1.0d+00 + gexpdt ) )
                else !delta is super small approching zero, which also indicates c and b are negligible
                    !since self-coagulation always happens, so solve dn/dt = -an^2
                    y = y0 / ( 1.0d+00 + a * y0 * tstep )
                endif
            elseif ((a .lt. 1.0d-30)) then
                y = y0
                call mass_solver(y, c, b,  tstep) !mass solver solving equation dy/dt=c-b*y 
            endif
        else                                        ! when c = 0.0d+00, as we assume c is not negative.
            expdt = exp( - b * tstep )
            if( 1.0d+00-expdt .gt. piq_thresh ) then                     ! if( expdt .lt. 1.0d+00 ) then
                y = b * y0 * expdt / ( b + a * y0 * ( 1.0d+00 - expdt ) )
            else
                y = y0 / ( 1.0d+00 + a * y0 * tstep )
            endif
        endif
        ni = max(y, minconc)
    end subroutine 

    subroutine module_mass_solver(ipop, ispec, rim_ipop_ispec, fi_ipop, is, ie, js, je, dt)
        integer, intent(in) :: ipop, ispec, is, ie, js, je
        real, intent(in) :: dt
        real, intent(in) :: rim_ipop_ispec(:,:,:), fi_ipop(:,:,:) !rim production from coagulation, fi_ipop: loss coefficient from coagulation
        real :: piq(size(rim_ipop_ispec,1),  size(rim_ipop_ispec,2), size(rim_ipop_ispec,3))
        integer :: i, j, k, it, jt, kt
        real :: v_bf
        it = size(rim_ipop_ispec,1)
        jt = size(rim_ipop_ispec,2)
        kt = size(rim_ipop_ispec,3)
        piq = 0.
        !dm_ipop_ispec/dt = (rim_ipop_ispec+source) -fi*m_ipop_ispec
        do i=1,it
        do j=1,jt
        do k=1,kt
        select case(ispec)
        case (1)
            piq(i,j,k) = rim_ipop_ispec(i,j,k) + matrix_all_tracer(matrix_all_pop(ipop)%i_msulf)%source(i+is-1, j+js-1, k)
            call mass_solver(matrix_all_tracer(matrix_all_pop(ipop)%i_msulf)%value_in_matrix(i+is-1, j+js-1, k), &
                piq(i,j,k), fi_ipop(i,j,k), dt)
        case (4)
            !write(*,*) 'dust_source, i,j, k, is, js:', matrix_all_tracer(matrix_all_pop(ipop)%i_mdust)%source(i+is-1, j+js-1, k), &
            !        i, j, k, is, js
            !write(*, *) 'input for mass solver: value_before, piq, fi_ipop, dt', &
            !        matrix_all_tracer(matrix_all_pop(ipop)%i_mdust)%value_in_matrix(i+is-1, j+js-1, k), &
            !        piq(i,j,k), fi_ipop(i,j,k), dt
            piq(i,j,k) = rim_ipop_ispec(i,j,k) + matrix_all_tracer(matrix_all_pop(ipop)%i_mdust)%source(i+is-1, j+js-1, k)
            call mass_solver(matrix_all_tracer(matrix_all_pop(ipop)%i_mdust)%value_in_matrix(i+is-1, j+js-1, k), &
                piq(i,j,k), fi_ipop(i,j,k), dt)
            !write(*, *) 'output for mass solver: value_after', &
            !        matrix_all_tracer(matrix_all_pop(ipop)%i_mdust)%value_in_matrix(i+is-1, j+js-1, k)
        case (5)
            piq(i,j,k) = rim_ipop_ispec(i,j,k) + matrix_all_tracer(matrix_all_pop(ipop)%i_mseas)%source(i+is-1, j+js-1, k)
            !write(*,*) 'ssalt_source, ipop, i,j, k, is, js:', matrix_all_tracer(matrix_all_pop(ipop)%i_mseas)%source(i+is-1, j+js-1, k), &
            !        ipop, i, j, k, is, js
            !write(*, *) 'input for mass solver: value_before, piq, fi_ipop, dt', &
            !        matrix_all_tracer(matrix_all_pop(ipop)%i_mseas)%value_in_matrix(i+is-1, j+js-1, k), &
            !        piq(i,j,k), fi_ipop(i,j,k), dt
            v_bf = matrix_all_tracer(matrix_all_pop(ipop)%i_mseas)%value_in_matrix(i+is-1, j+js-1, k)
            call mass_solver(matrix_all_tracer(matrix_all_pop(ipop)%i_mseas)%value_in_matrix(i+is-1, j+js-1, k), &
                piq(i,j,k), fi_ipop(i,j,k), dt)
            !write(*, *) 'output for mass solver: v_bf, v_bf+ source*dt, value_after', v_bf, v_bf+piq(i,j,k)*dt, &
            !               matrix_all_tracer(matrix_all_pop(ipop)%i_mseas)%value_in_matrix(i+is-1, j+js-1, k)

        case (3)
            piq(i,j,k) = rim_ipop_ispec(i,j,k) + matrix_all_tracer(matrix_all_pop(ipop)%i_mocar)%source(i+is-1, j+js-1, k)
            call mass_solver(matrix_all_tracer(matrix_all_pop(ipop)%i_mocar)%value_in_matrix(i+is-1, j+js-1, k), &
                piq(i,j,k), fi_ipop(i,j,k), dt)
        case (2)
            piq(i,j,k) = rim_ipop_ispec(i,j,k) + matrix_all_tracer(matrix_all_pop(ipop)%i_mbcar)%source(i+is-1, j+js-1, k)
            call mass_solver(matrix_all_tracer(matrix_all_pop(ipop)%i_mbcar)%value_in_matrix(i+is-1, j+js-1, k), &
                piq(i,j,k), fi_ipop(i,j,k), dt)
        end select
        enddo
        enddo
        enddo
    end subroutine

    !---------------------------------------------------------------
    ! subroutine mass_solver
    ! for number tracers: dmi/dt = piq - fi*mi
    ! piq: total source, from emision + npf + coagulation generation
    ! fi: intermodal coagulation loss
    !---------------------------------------------------------------
    subroutine mass_solver(miq, piq, fi_i,  tstep)
        real, intent(in) ::  piq, fi_i, tstep
        real, intent(inout) :: miq
        real :: expdt, piq_thresh, minconc, factor
        real :: miq_0
        miq_0 = miq
        minconc = 1.0e-32
        piq_thresh = 1.0d-08  ! [1] threshold in number/mass conc. solver
        expdt = exp( - fi_i * tstep )
        if ( 1.0d+00-expdt .gt. piq_thresh ) then
            factor = ( 1.0d+00 - expdt ) / fi_i
            miq = miq_0 * expdt + piq * factor
        else
            miq = miq_0 + piq * tstep
        endif
        miq = max(miq, minconc)

    end subroutine


    !----------------------------------------------------------------
    !        subroutine matrix_intermodal_transfer
    ! perform intermodal transfer : 
    !       (1) akk -> acc: aitken mode sulfate -> accumulation mode sulfate
    !       (2) oc1 -> oc2: hydrophobic organic carbon -> hydrophilic organic carbon
    !       (3) bc1 -> bc2: hydrophobic black carbon -> hydrophilic carbon
    !!!!!!!!!!note: matrix value got updated in this subroutine
    !--------------------------------------------------------------------
    subroutine matrix_intermodal_transfer(config, r, is,ie,js,je, time_next)
        character(len=32), intent(in) :: config
        real, intent(in) :: r(:,:,:,:)
        type(time_type), intent(in) :: time_next
        logical :: used
        integer, intent(in) :: is,ie,js,je
        ! these variables are used in intermodal transfer.
        real(8), parameter :: dg_akk_param      = 0.026d+00      !giss set-up
        real(8), parameter :: dg_acc_param      = 0.110d+00      ! e04, table 2, accumulation mode 
        integer, parameter :: imtr_exp          = 4              ! exponent for akk --> acc intermodal transfer
        real, parameter    :: imtr_method       = 1     ! =1 no cut of pdf, =2 fixed-dp cut, =3 variable-dp cut as in cmaq
        real(8)            :: xnum, x3                  ! error function complement arguments [1]
        real(8)            :: dgn_akk_imtr              ! geo. mean diam. of the akk mode number distribution [m]
        real(8)            :: dgn_acc_imtr              ! geo. mean diam. of the acc mode number distribution [m]
        real(8)            :: del_numb                  ! number conc. transferred from akk to acc [ #/m^3]
        real(8)            :: del_mass                  ! mass   conc. transferred from akk to acc [ug/m^3]
        real(8)      :: dpcut_imtr     ! = 0.0d+00 ! fixed diameter for intermodal transfer [um]
        real(8)      :: xnum_factor    ! = 0.0d+00 ! factor in xnum expression [1]
        real(8)      :: x3_term        ! = 0.0d+00 ! term   in x3   expression [1]
        real(8)      :: lnakk_sigma    ! = 0.0d+00 ! ln(sg_akk) [1]
        real(8)      :: lnacc_sigma    ! = 0.0d+00 ! ln(sg_acc) [1]
        real, parameter    :: dnu             = 3.0d-09   ! diameter of a new particle [nm]
        real(8), parameter :: dpakk0          = dnu     ! min. diameter of average mass for mode akk [m]
        real(8), parameter :: fnum_max        = 0.5d+00 ! max. value of fnum [1]
        real(8), parameter :: akk_minnum_imtr = 1.0d+06 ! min. akk number conc. to enable imtr [#/m^3]
        real               :: dpakk_all(size(r,1),size(r,2),size(r,3)) ! diameter of average mass for mode akk [m]
        real               :: dpacc_all(size(r,1),size(r,2),size(r,3)) !diameter of average mass for mode acc [m]
        real               :: dgakk_all(size(r,1),size(r,2),size(r,3)) !geometric mean diameter [m]
        real               :: dgacc_all(size(r,1),size(r,2),size(r,3))
        real               :: acc_num_all(size(r,1),size(r,2),size(r,3))
        real               :: akk_num_all(size(r,1),size(r,2),size(r,3))
        real               :: dpakk, dpacc
        real(8)            :: fnum_all(size(r,1),size(r,2),size(r,3))! fraction of akk number transferred over the time step [1]
        real(8)            :: f3_all(size(r,1),size(r,2),size(r,3))  ! fraction of akk mass   transferred over the time step [1]
        real               :: fnum, f3
        real, parameter :: mimr_bc1 = 0.05 !threshhold of sulfate/bc mass ratio to transfer bc1 to bc2
        real, parameter :: mimr_oc1 = 0.05 !threshhold of sulfate/oc mass ratio to transfer oc1 to oc2
        real, parameter :: transfer_factor = 0.5 !how much mass and number got transfered per timestep when modal transfer happen
        real :: dmsulf_bc1_to_bc2(size(r,1),size(r,2),size(r,3)), dmbcar_bc1_to_bc2(size(r,1),size(r,2),size(r,3)), dn_bc1_to_bc2(size(r,1),size(r,2),size(r,3))
        real :: dmsulf_oc1_to_oc2(size(r,1),size(r,2),size(r,3)), dmocar_oc1_to_oc2(size(r,1),size(r,2),size(r,3)), dn_oc1_to_oc2(size(r,1),size(r,2),size(r,3))
        real :: dmsulf_akk_to_acc(size(r,1),size(r,2),size(r,3)), dn_akk_to_acc(size(r,1),size(r,2),size(r,3))
        real :: bc1_sulf_mass(size(r,1),size(r,2),size(r,3)), bc1_bcar_mass(size(r,1),size(r,2),size(r,3)), bc1_num(size(r,1),size(r,2),size(r,3))
        real :: bc2_sulf_mass(size(r,1),size(r,2),size(r,3)), bc2_bcar_mass(size(r,1),size(r,2),size(r,3)), bc2_num(size(r,1),size(r,2),size(r,3))
        real :: oc1_sulf_mass(size(r,1),size(r,2),size(r,3)), oc1_ocar_mass(size(r,1),size(r,2),size(r,3)), oc1_num(size(r,1),size(r,2),size(r,3))
        real :: oc2_sulf_mass(size(r,1),size(r,2),size(r,3)), oc2_ocar_mass(size(r,1),size(r,2),size(r,3)), oc2_num(size(r,1),size(r,2),size(r,3))
        real :: akk_sulf_mass(size(r,1),size(r,2),size(r,3)), acc_sulf_mass(size(r,1),size(r,2),size(r,3))
        real :: i,j,k,it,jt,kt

        dpcut_imtr = 0.
        xnum_factor = 0.
        x3_term = 0.
        lnakk_sigma = 0.
        lnacc_sigma = 0.        
        fnum = 0.
        f3 = 0.
        dpakk = 0.
        dpacc = 0.
        it = size(r,1)
        jt = size(r,2)
        kt = size(r,3) 
        !bc1 -> bc2: hydrophobic to hydrophilic
        dmsulf_bc1_to_bc2 = 0.
        dmbcar_bc1_to_bc2 = 0.
        dn_bc1_to_bc2 = 0.
        !oc1 -> oc2: hydrophobic to hydrophilic
        dmsulf_oc1_to_oc2 = 0
        dmocar_oc1_to_oc2 = 0.
        dn_oc1_to_oc2 = 0.
        !akk -> acc: aitken mode to accumulation mode
        dmsulf_akk_to_acc = 0.
        dn_akk_to_acc = 0.
        bc1_sulf_mass = 0.
        bc2_sulf_mass = 0.
        bc1_bcar_mass = 0.
        bc2_bcar_mass = 0.
        bc1_num = 0.
        bc2_num = 0.
        oc1_sulf_mass = 0.
        oc2_sulf_mass = 0.
        oc1_ocar_mass = 0.
        oc2_ocar_mass = 0.
        oc1_num = 0.
        oc2_num = 0.
        if (trim(lowercase(config)) .eq. "full") then
            !-------------------------------------------------------------------
            ! transfer mode bc1 to bc2
            !-------------------------------------------------------------------
            if (i_bc1 > 0 .and. i_bc2 > 0 ) then
                if (matrix_all_pop(i_bc1)%nb_tracer_pop > 0 .and. matrix_all_pop(i_bc2)%nb_tracer_pop > 0) then
                    bc1_sulf_mass = matrix_all_tracer(matrix_all_pop(i_bc1)%i_msulf)%value_in_matrix(is:ie,js:je,:) !sulf mass in bc1 [ug/m3]
                    bc1_bcar_mass = matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%value_in_matrix(is:ie,js:je,:) !bc mass in bc1 [ug/m3]
                    bc1_num = matrix_all_tracer(matrix_all_pop(i_bc1)%i_n)%value_in_matrix(is:ie,js:je,:) !number of bc1 [#/m3]
                    bc2_sulf_mass = matrix_all_tracer(matrix_all_pop(i_bc2)%i_msulf)%value_in_matrix(is:ie,js:je,:) !sulf mass in bc2 [ug/m3]
                    bc2_bcar_mass = matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%value_in_matrix(is:ie,js:je,:) !bc mass in bc2 [ug/m3]
                    bc2_num = matrix_all_tracer(matrix_all_pop(i_bc2)%i_n)%value_in_matrix(is:ie,js:je,:) !number of bc2 [#/m3]
                    do i = 1, it
                    do j = 1, jt
                    do k = 1, kt
                    if (bc1_bcar_mass(i,j,k) > 0.) then
                        if (bc1_sulf_mass(i,j,k)/bc1_bcar_mass(i,j,k) .gt. mimr_bc1) then
                            dmsulf_bc1_to_bc2(i,j,k) = transfer_factor*bc1_sulf_mass(i,j,k)
                            dmbcar_bc1_to_bc2(i,j,k) = transfer_factor*bc1_bcar_mass(i,j,k)
                            dn_bc1_to_bc2(i,j,k) = transfer_factor*bc1_num(i,j,k)
                        endif
                    endif
                    enddo
                    enddo
                    enddo
                    !update matrix value for bc1 pop
                    matrix_all_tracer(matrix_all_pop(i_bc1)%i_msulf)%value_in_matrix(is:ie,js:je,:) = bc1_sulf_mass - dmsulf_bc1_to_bc2
                    matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%value_in_matrix(is:ie,js:je,:) = bc1_bcar_mass - dmbcar_bc1_to_bc2
                    matrix_all_tracer(matrix_all_pop(i_bc1)%i_n)%value_in_matrix(is:ie,js:je,:) = bc1_num - dn_bc1_to_bc2
                    !update matrix value for bc2 pop
                    matrix_all_tracer(matrix_all_pop(i_bc2)%i_msulf)%value_in_matrix(is:ie,js:je,:) = bc2_sulf_mass + dmsulf_bc1_to_bc2
                    matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%value_in_matrix(is:ie,js:je,:) = bc2_bcar_mass + dmbcar_bc1_to_bc2
                    matrix_all_tracer(matrix_all_pop(i_bc2)%i_n)%value_in_matrix(is:ie,js:je,:) = bc2_num + dn_bc1_to_bc2

                    !send diagnostic data
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_bc1)%i_msulf)%id_tracer_intermodal, -dmsulf_bc1_to_bc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_bc1)%i_mbcar)%id_tracer_intermodal, -dmbcar_bc1_to_bc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_bc1)%i_n)%id_tracer_intermodal, -dn_bc1_to_bc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )

                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_bc2)%i_msulf)%id_tracer_intermodal, dmsulf_bc1_to_bc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_bc2)%i_mbcar)%id_tracer_intermodal, dmbcar_bc1_to_bc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_bc2)%i_n)%id_tracer_intermodal, dn_bc1_to_bc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )


                else
                    call error_mesg ('matrix_gfdl','bc1 and bc2: population and tracers not consistent', fatal)
                endif
            endif
            !-------------------------------------------------------------------
            ! transfer mode oc1 to oc2
            !-------------------------------------------------------------------
            if (i_oc1 > 0 .and. i_oc2 > 0 ) then
                if (matrix_all_pop(i_oc1)%nb_tracer_pop > 0 .and. matrix_all_pop(i_oc2)%nb_tracer_pop > 0) then
                    oc1_sulf_mass = matrix_all_tracer(matrix_all_pop(i_oc1)%i_msulf)%value_in_matrix(is:ie,js:je,:) !sulf mass in oc1 [ug/m3]
                    oc1_ocar_mass = matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%value_in_matrix(is:ie,js:je,:) !oc mass in oc1 [ug/m3]
                    oc1_num = matrix_all_tracer(matrix_all_pop(i_oc1)%i_n)%value_in_matrix(is:ie,js:je,:) !number of oc1 [#/m3]
                    oc2_sulf_mass = matrix_all_tracer(matrix_all_pop(i_oc2)%i_msulf)%value_in_matrix(is:ie,js:je,:)
                    oc2_ocar_mass = matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%value_in_matrix(is:ie,js:je,:)
                    oc2_num = matrix_all_tracer(matrix_all_pop(i_oc2)%i_n)%value_in_matrix(is:ie,js:je,:)
                    do i = 1, it
                    do j = 1, jt
                    do k = 1, kt
                    if (oc1_ocar_mass(i,j,k) > 0) then
                        if (oc1_sulf_mass(i,j,k)/oc1_ocar_mass(i,j,k) .gt. mimr_oc1) then
                            dmsulf_oc1_to_oc2(i,j,k) = transfer_factor*oc1_sulf_mass(i,j,k)
                            dmocar_oc1_to_oc2(i,j,k) = transfer_factor*oc1_ocar_mass(i,j,k)
                            dn_oc1_to_oc2(i,j,k) = transfer_factor*oc1_num(i,j,k)
                        endif
                    endif
                    enddo
                    enddo
                    enddo
                    !update matrix value for oc1 pop
                    matrix_all_tracer(matrix_all_pop(i_oc1)%i_msulf)%value_in_matrix(is:ie,js:je,:) = oc1_sulf_mass - dmsulf_oc1_to_oc2
                    matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%value_in_matrix(is:ie,js:je,:) = oc1_ocar_mass - dmocar_oc1_to_oc2
                    matrix_all_tracer(matrix_all_pop(i_oc1)%i_n)%value_in_matrix(is:ie,js:je,:) = oc1_num - dn_oc1_to_oc2
                    !update matrix value for oc2 pop
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_msulf)%value_in_matrix(is:ie,js:je,:) = oc2_sulf_mass + dmsulf_oc1_to_oc2
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%value_in_matrix(is:ie,js:je,:) = oc2_ocar_mass + dmocar_oc1_to_oc2
                    matrix_all_tracer(matrix_all_pop(i_oc2)%i_n)%value_in_matrix(is:ie,js:je,:) = oc2_num + dn_oc1_to_oc2
                    !send diagnostic data
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_oc1)%i_msulf)%id_tracer_intermodal, -dmsulf_oc1_to_oc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_oc1)%i_mocar)%id_tracer_intermodal,-dmocar_oc1_to_oc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_oc1)%i_n)%id_tracer_intermodal, -dn_oc1_to_oc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )

                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_oc2)%i_msulf)%id_tracer_intermodal, dmsulf_oc1_to_oc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_oc2)%i_mocar)%id_tracer_intermodal,dmocar_oc1_to_oc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_oc2)%i_n)%id_tracer_intermodal, dn_oc1_to_oc2, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                else
                    call error_mesg ('matrix_gfdl','oc1 and oc2: population and tracers not consistent', fatal)
                endif
            endif
            !-------------------------------------------------------------------
            ! transfer mode akk to acc
            !-------------------------------------------------------------------
            if (i_akk > 0 .and. i_acc > 0) then
                if (matrix_all_pop(i_akk)%nb_tracer_pop > 0 .and. matrix_all_pop(i_acc)%nb_tracer_pop > 0) then
                    lnakk_sigma = log(matrix_all_pop(i_akk)%sigma)
                    lnacc_sigma = log(matrix_all_pop(i_acc)%sigma)
                    xnum_factor = 1.0d+00 / ( sqrt( 2.0d+00 ) * lnakk_sigma )
                    x3_term     = 3.0d+00 * lnakk_sigma / sqrt( 2.0d+00 )
                    dpcut_imtr  = sqrt( dg_akk_param * dg_acc_param )   ! this is the formula of easter et al. 2004. (= 0.053 um)
                    dpakk_all = matrix_all_pop(i_akk)%dg_dry(is:ie,js:je,:) * exp(1.5* lnakk_sigma**2) !diameter of average mass for akk [m]
                    dpacc_all = matrix_all_pop(i_acc)%dg_dry(is:ie,js:je,:) * exp(1.5* lnacc_sigma**2) !dry diameter of average mass for acc [m]
                    dgakk_all = matrix_all_pop(i_akk)%dg_dry(is:ie,js:je,:) !geometric diameter [m]
                    dgacc_all = matrix_all_pop(i_acc)%dg_dry(is:ie,js:je,:) !geometric diameter [m]
                    acc_num_all = matrix_all_tracer(matrix_all_pop(i_acc)%i_n)%value_in_matrix(is:ie,js:je,:)
                    akk_num_all = matrix_all_tracer(matrix_all_pop(i_akk)%i_n)%value_in_matrix(is:ie,js:je,:)
                    akk_sulf_mass = matrix_all_tracer(matrix_all_pop(i_akk)%i_msulf)%value_in_matrix(is:ie,js:je,:) !sulf mass in akk [ug/m3]
                    acc_sulf_mass = matrix_all_tracer(matrix_all_pop(i_acc)%i_msulf)%value_in_matrix(is:ie,js:je,:) !sulf mass in acc [ug/m3]
                    do i = 1, it
                    do j = 1, jt
                    do k = 1, kt
                    dpakk = dpakk_all(i,j,k)
                    dpacc = dpacc_all(i,j,k)
                    if ((dpakk > 0) .and. (dpacc > 0.)) then
                        if( imtr_method .eq. 1 ) then
                            !------------------------------------------------------------------------------------------------------------
                            ! calculate the fraction transferred based on the relative difference
                            !   in mass mean diameters of the akk and acc modes.
                            !------------------------------------------------------------------------------------------------------------
                            if( (dpakk .ge. dpakk0) .and. (dpacc .ge. dpakk0)) then                ! [m], dpakk0 = 3 nm
                                fnum = ( ( dpakk - dpakk0 ) / ( dpacc - dpakk0 ) )**imtr_exp   ! fraction transferred from akk to acc
                                fnum = max( min( fnum, fnum_max ), 0.0d+00 )                   ! limit transfer in a single transfer
                            else
                                fnum = 0.0d+00
                            endif
                            f3 = fnum
                            ! write(34,'(7d12.4)')fnum,f3
                        elseif( imtr_method .eq. 2 ) then
                            !------------------------------------------------------------------------------------------------------------
                            ! calculate the fraction transferred based on a fixed
                            !   threshold diameter dpcut_imtr.
                            !------------------------------------------------------------------------------------------------------------
                            dgn_akk_imtr = 1.0d+06 * dgakk_all(i,j,k)                      ! geometric diameter in [um]
                            xnum = xnum_factor * log( dpcut_imtr / dgn_akk_imtr )          ! [1]
                            xnum = max( xnum, x3_term )                                    ! limit for stability as in bs2003
                            x3 = xnum - x3_term                                            ! [1]
                            fnum = 0.5d+00 * erfc( xnum )                                  ! number fraction transferred from akk to acc
                            f3   = 0.5d+00 * erfc( x3   )                                  ! mass   fraction transferred from akk to acc
                            ! write(34,'(9d12.4)')dgn_akk_imtr,dpcut_imtr,dpakk*1.0d+06,aero(numb_akk_1),aero(mass_akk_sulf),fnum,f3
                        elseif( imtr_method .eq. 3 ) then
                            !------------------------------------------------------------------------------------------------------------
                            ! calculate the fraction transferred based on the
                            !   diameter of intersection of the akk and acc modes.
                            !------------------------------------------------------------------------------------------------------------
                            dgn_akk_imtr = 1.0d+06 * dgakk_all(i,j,k)                      ! [um]
                            dgn_acc_imtr = 1.0d+06 * dgacc_all(i,j,k)                      ! [um]
                            if( acc_num_all(i,j,k) .gt. 1.0d+06 ) then
                                xnum = getxnum(akk_num_all(i,j,k), acc_num_all(i,j,k), &
                                    dgn_akk_imtr,dgn_acc_imtr,lnakk_sigma,lnacc_sigma)   ! [1]
                            else                                                           ! mode acc essentially empty - use method 2
                                xnum = xnum_factor * log( dpcut_imtr / dgn_akk_imtr )        ! [1]
                            endif
                            xnum = max( xnum, x3_term )                                    ! limit for stability as in bs2003
                            x3 = xnum - x3_term                                            ! [1]
                            fnum = 0.5d+00 * erfc( xnum )                                  ! number fraction transferred from akk to acc
                            f3   = 0.5d+00 * erfc( x3   )                                  ! mass   fraction transferred from akk to acc
                        endif
                    endif
                    fnum_all(i,j,k) = fnum
                    f3_all(i,j,k) = f3
                    enddo
                    enddo
                    enddo

                    dmsulf_akk_to_acc = f3_all * akk_sulf_mass  !mass concentration transferred [ug/m3] 
                    dn_akk_to_acc = fnum_all * akk_num_all
                    !update matrix value for akk and acc
                    matrix_all_tracer(matrix_all_pop(i_akk)%i_n)%value_in_matrix(is:ie,js:je,:) = akk_num_all - dn_akk_to_acc
                    matrix_all_tracer(matrix_all_pop(i_akk)%i_msulf)%value_in_matrix(is:ie,js:je,:) = akk_sulf_mass - dmsulf_akk_to_acc
                    matrix_all_tracer(matrix_all_pop(i_acc)%i_n)%value_in_matrix(is:ie,js:je,:) = acc_num_all + dn_akk_to_acc
                    matrix_all_tracer(matrix_all_pop(i_acc)%i_msulf)%value_in_matrix(is:ie,js:je,:) = acc_sulf_mass + dmsulf_akk_to_acc

                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_akk)%i_msulf)%id_tracer_intermodal,-dmsulf_akk_to_acc, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_akk)%i_n)%id_tracer_intermodal, -dn_akk_to_acc, time_next, &
                        is_in=is, js_in=js, ks_in=1 )

                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_acc)%i_msulf)%id_tracer_intermodal, dmsulf_akk_to_acc, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                    used =  send_data(matrix_all_tracer(matrix_all_pop(i_acc)%i_n)%id_tracer_intermodal, dn_akk_to_acc, time_next, &
                        is_in=is, js_in=js, ks_in=1 )
                else
                    call error_mesg ('matrix_gfdl','akk and acc: population and tracers not consistent', fatal)
                endif
            endif

        else
            call error_mesg ('matrix_gfdl','intermodal_transfer configuration not found', fatal)
        endif

    end subroutine matrix_intermodal_transfer

    !only calculate frac_in_cloud_uw_mass
    subroutine grid_wetdep_fic_uw_mass(n, mass_in, grid_frac_in_cloud_uw)
        real, intent(in) :: mass_in(5) !mass spec for 5 specs: 1:sulf; 2:bcar; 3:ocar; 4:dust; 5:seas 
        integer, intent(in) :: n ! tarcer index
        real, intent(out) :: grid_frac_in_cloud_uw
        !hygroscopicity factor of different species
        real, parameter :: kappa_sulf = 0.507, kappa_bc_phil = 0.1, kappa_oc_phil = 0.12
        real, parameter :: kappa_dust = 0.068, kappa_seas = 1.1, kappa_bc_phob = 0, kappa_oc_phob = 0
        !frac_in_cloud_uw defined in gfdl field table based on species
        real, parameter :: fic_uw_gfdl_sulf = 0.44, fic_uw_gfdl_soa = 0.44, fic_uw_gfdl_dust = 0.3
        real, parameter :: fic_uw_gfdl_ssalt_f = 0.44, fic_uw_gfdl_ssalt_c = 0.55, fic_uw_gfdl_bcphil = 0.44
        real, parameter :: fic_uw_gfdl_ocphil = 0.44
        real :: fic_uw_seas, fic_uw_dust
        real :: mtot, msulf, mbcar, mseas, mocar, mdust
        mtot = 0.
        grid_frac_in_cloud_uw = 0.
        fic_uw_seas = 0.
        fic_uw_dust = 0.
        if (.not. matrix_all_tracer(n)%is_active) then
            call error_mesg ('matrix_gfdl', 'grid_wetdep_param shoud not be assigned to non-matrix tracers', fatal)
        else
            !if matrix tracer, dynamic is only assigned to mxc
            if (matrix_all_tracer(n)%pop == 'mxc') then
                fic_uw_seas = fic_uw_gfdl_ssalt_c
                fic_uw_dust = fic_uw_gfdl_dust
            elseif (matrix_all_tracer(n)%pop == 'mxa') then
                fic_uw_seas = fic_uw_gfdl_ssalt_f
                fic_uw_dust = fic_uw_gfdl_dust
            else
                write(*,*) 'grid_wetdep calc, pop name:', matrix_all_tracer(n)%pop
                call error_mesg ('matrix_gfdl', 'grid_wetdep has not been activated for non-mixing populations', fatal)
            endif    
            msulf = max(0., mass_in(1))
            mbcar = max(0., mass_in(2))
            mocar = max(0., mass_in(3))
            mdust = max(0., mass_in(4))
            mseas = max(0., mass_in(5))
            mtot = max(1.0e-30, msulf + mocar + mbcar + mdust + mseas)
            grid_frac_in_cloud_uw = (msulf*fic_uw_gfdl_sulf + mocar*fic_uw_gfdl_ocphil + mbcar*fic_uw_gfdl_bcphil + mseas*fic_uw_seas + &
                mdust*fic_uw_dust)/mtot
        endif
    end subroutine grid_wetdep_fic_uw_mass



    !only calculate frac_in_cloud_uw
    subroutine grid_wetdep_fic_uw(n, kappa_grid, grid_frac_in_cloud_uw)
        real, intent(in) :: kappa_grid
        integer, intent(in) :: n ! tarcer index
        real, intent(out) :: grid_frac_in_cloud_uw
        !hygroscopicity factor of different species
        real, parameter :: kappa_sulf = 0.507, kappa_bc_phil = 0.1, kappa_oc_phil = 0.12
        real, parameter :: kappa_dust = 0.068, kappa_seas = 1.1, kappa_bc_phob = 0, kappa_oc_phob = 0
        !frac_in_cloud_uw defined in gfdl field table based on species
        real, parameter :: fic_uw_gfdl_sulf = 0.44, fic_uw_gfdl_soa = 0.44, fic_uw_gfdl_dust = 0.3
        real, parameter :: fic_uw_gfdl_ssalt_f = 0.44, fic_uw_gfdl_ssalt_c = 0.55, fic_uw_gfdl_bcphil = 0.44
        real, parameter :: fic_uw_gfdl_ocphil = 0.44
        if (.not. matrix_all_tracer(n)%is_active) then
            call error_mesg ('matrix_gfdl', 'grid_wetdep_param shoud not be assigned to non-matrix tracers', fatal)
        else
            !if matrix tracer, dynamic is only assigned to mxc
            if (matrix_all_tracer(n)%pop == 'mxc') then
                grid_frac_in_cloud_uw = max(kappa_grid-kappa_dust, 1e-30)/(kappa_seas-kappa_dust) &
                    * (fic_uw_gfdl_ssalt_c-fic_uw_gfdl_dust) + fic_uw_gfdl_dust
            elseif (matrix_all_tracer(n)%pop == 'mxa') then ! if mxa mode
                grid_frac_in_cloud_uw = max(kappa_grid-kappa_dust, 1e-30)/(kappa_seas-kappa_dust) &
                    * (fic_uw_gfdl_ssalt_f-fic_uw_gfdl_dust) + fic_uw_gfdl_dust
            else
                write(*,*) 'grid_wetdep calc, pop name:', matrix_all_tracer(n)%pop
                call error_mesg ('matrix_gfdl', 'grid_wetdep has not been activated for non-mixing populations', fatal)
            endif
        endif
    end subroutine grid_wetdep_fic_uw



    !only calculate frac_in_cloud
    subroutine grid_wetdep_fic(kappa_grid, grid_frac_in_cloud, size_flag)
        real, intent(in) :: kappa_grid
        integer, intent(in) :: size_flag ! if size_flag = 1, mxa; size_flag = 2, mxc
        real, intent(out) :: grid_frac_in_cloud
        !hygroscopicity factor of different species
        real, parameter :: kappa_sulf = 0.507, kappa_bc_phil = 0.1, kappa_oc_phil = 0.12
        real, parameter :: kappa_dust = 0.068, kappa_seas = 1.1, kappa_bc_phob = 0, kappa_oc_phob = 0
        !frac_in_cloud defined in gfdl field table based on species
        real, parameter :: fic_gfdl_sulf = 0.33, fic_gfdl_soa = 0.33, fic_gfdl_dust = 0.2
        real, parameter :: fic_gfdl_ssalt_f = 0.33, fic_gfdl_ssalt_c = 0.44, fic_gfdl_bcphil = 0.33
        real, parameter :: fic_gfdl_ocphil = 0.33
        if (size_flag == 2) then !if it is for mxc mode
            grid_frac_in_cloud = max(kappa_grid-kappa_dust, 1e-30)/(kappa_seas-kappa_dust) &
                * (fic_gfdl_ssalt_c-fic_gfdl_dust) + fic_gfdl_dust
        elseif (size_flag == 1) then ! if mxa mode
            grid_frac_in_cloud = max(kappa_grid-kappa_dust, 1e-30)/(kappa_seas-kappa_dust) &
                * (fic_gfdl_ssalt_f-fic_gfdl_dust) + fic_gfdl_dust
        else
            call error_mesg ('matrix_gfdl', 'grid_wetdep_param size_flag not preoperly defined', fatal) 
        endif
    end subroutine grid_wetdep_fic

    subroutine matrix_query_dynamic_wetdep_param_2d(n, pfull, is, js, k,frac_in_cloud_2d) 
        !frac_in_cloud_snow_2d, &
        !frac_in_cloud_snow_homogeneous_2d, alphar_2d, alphas_2d, frac_in_cloud_uw_2d)
        real, intent(in), dimension(:,:,:) :: pfull
        integer, intent(in) :: n
        integer, intent(in) :: is, js, k !only 1 vertical layer of z is passed in
        real, intent(out) :: frac_in_cloud_2d(:,:)
        !real, intent(out) :: frac_in_cloud_snow_2d(:,:)
        !real, intent(out) :: frac_in_cloud_snow_homogeneous_2d(:,:)
        !real, intent(out) :: alphar_2d(:,:)
        !real, intent(out) :: alphas_2d(:,:)
        !real, intent(out), optional :: frac_in_cloud_uw_2d(:,:)
        real :: kappa_dust, kappa_seas, kappa_grid
        integer :: im, jm , ie, je
        integer, parameter :: size_flag_acc = 1, size_flag_coarse = 2
        ie = is + size(pfull, 1) -1
        je = js + size(pfull, 2) -1
        frac_in_cloud_2d = 0.
        !frac_in_cloud_snow_2d = 0.
        !frac_in_cloud_snow_homogeneous_2d = 0.
        !alphar_2d = 0.
        !alphas_2d = 0.

        !see if is matrix tracer
        if (.not. matrix_all_tracer(n)%is_active) then
            return
        else
            !if matrix tracer, dynamic is only assigned to mxc
            if (matrix_all_tracer(n)%pop == 'mxc') then
                do im = is, ie
                do jm = js, je
                kappa_grid = matrix_all_pop(i_mxc)%kappa_pop(im,jm,k)
                call grid_wetdep_fic(kappa_grid, frac_in_cloud_2d(im-is+1, jm-js+1), size_flag_coarse)
                enddo
                enddo
            elseif (matrix_all_tracer(n)%pop == 'mxa') then
                do im = is, ie
                do jm = js, je
                kappa_grid = matrix_all_pop(i_mxa)%kappa_pop(im,jm,k)
                call grid_wetdep_fic(kappa_grid, frac_in_cloud_2d(im-is+1, jm-js+1), size_flag_acc)
                enddo
                enddo
            else
                call error_mesg ('matrix_gfdl','dyanmic wet_dep property should be only assigned to mxa/mxc population', fatal) 
            endif
        endif
    end subroutine matrix_query_dynamic_wetdep_param_2d 

    
    subroutine matrix_query_dynamic_wetdep_param_2d_mass(n, pfull, is, js, k,frac_in_cloud_2d)
        real, intent(in), dimension(:,:,:) :: pfull
        integer, intent(in) :: n
        integer, intent(in) :: is, js, k !only 1 vertical layer of z is passed in
        real, intent(out) :: frac_in_cloud_2d(:,:)
        real :: msulf(size(pfull,1), size(pfull,2)), mbcar(size(pfull,1), size(pfull,2)), mseas(size(pfull,1), size(pfull,2)),mdust(size(pfull,1), size(pfull,2))
        real :: mocar(size(pfull,1), size(pfull,2)), mtot(size(pfull,1), size(pfull,2))
        integer :: im, jm , ie, je, pop_in
        integer, parameter :: size_flag_acc = 1, size_flag_coarse = 2
        real, parameter :: kappa_sulf = 0.507, kappa_bc_phil = 0.1, kappa_oc_phil = 0.12
        real, parameter :: kappa_dust = 0.068, kappa_seas = 1.1, kappa_bc_phob = 0, kappa_oc_phob = 0
        !frac_in_cloud defined in gfdl field table based on species
        real, parameter :: fic_gfdl_sulf = 0.33, fic_gfdl_soa = 0.33, fic_gfdl_dust = 0.2
        real, parameter :: fic_gfdl_ssalt_f = 0.33, fic_gfdl_ssalt_c = 0.44, fic_gfdl_bcphil = 0.33
        real, parameter :: fic_gfdl_ocphil = 0.33
        real :: fic_seas, fic_dust
        fic_seas = 0.
        fic_dust = 0.
        msulf = 0.
        mdust = 0.
        mocar = 0.
        mbcar = 0.
        mseas = 0.
        mtot = 0.
        pop_in = 0
        ie = is + size(pfull, 1) -1
        je = js + size(pfull, 2) -1
        frac_in_cloud_2d = 0.

        if (do_matrix) then
        !see if is matrix tracer
        if (.not. matrix_all_tracer(n)%is_active) then
            return
        else
            !if matrix tracer, dynamic is only assigned to mxc
            if (matrix_all_tracer(n)%pop == 'mxc') then
                pop_in = i_mxc
                fic_seas = fic_gfdl_ssalt_c
                fic_dust = fic_gfdl_dust
            elseif (matrix_all_tracer(n)%pop == 'mxa') then
                pop_in = i_mxa
                fic_seas = fic_gfdl_ssalt_f
                fic_dust = fic_gfdl_dust
           else
                call error_mesg ('matrix_gfdl','dyanmic wet_dep property should be only assigned to mxa/mxc population', fatal)
           endif
           msulf = matrix_all_tracer(matrix_all_pop(pop_in)%i_msulf)%value_in_matrix(is:ie,js:je,k)
           mocar = matrix_all_tracer(matrix_all_pop(pop_in)%i_mocar)%value_in_matrix(is:ie,js:je,k)
           mbcar = matrix_all_tracer(matrix_all_pop(pop_in)%i_mbcar)%value_in_matrix(is:ie,js:je,k)
           mdust = matrix_all_tracer(matrix_all_pop(pop_in)%i_mdust)%value_in_matrix(is:ie,js:je,k)
           mseas = matrix_all_tracer(matrix_all_pop(pop_in)%i_mseas)%value_in_matrix(is:ie,js:je,k)
           mtot = max(1.0e-30, msulf + mocar + mbcar + mdust + mseas)
           frac_in_cloud_2d = (msulf*fic_gfdl_sulf + mocar*fic_gfdl_ocphil + mbcar*fic_gfdl_bcphil + mseas*fic_seas + mdust*fic_dust)/mtot
        endif
        endif

    end subroutine matrix_query_dynamic_wetdep_param_2d_mass





    subroutine query_matrix_tracer(ind_tracer, flag) !check if a tracer is matrix tracer
        integer, intent(in) :: ind_tracer !index of tracer from query method
        integer, intent(out) :: flag  !no_active -1; active_mass_tracer = 1; active_number_tracer = 2
        integer, parameter :: flag_n = 2, flag_m = 1
        flag = -1
        if (do_matrix) then
        if (matrix_all_tracer(ind_tracer)%is_active) then
            if (lowercase(trim(matrix_all_tracer(ind_tracer)%type)) .eq. 'mass') then
                flag = flag_m
            elseif (lowercase(trim(matrix_all_tracer(ind_tracer)%type)) .eq. 'number') then
                flag = flag_n
            endif
        endif
        endif
    end subroutine

    !the following subroutine is used in shallow cu
    subroutine query_tracer_in_pop(ind_tracer, ipop_active_index, ipop_abs, flag)
        integer, intent(in) :: ind_tracer
        integer, intent(out) :: ipop_active_index
        integer, intent(out) :: ipop_abs
        integer, intent(out) :: flag
        ipop_active_index=0
        ipop_abs=0
        flag = 0 !if find, then flag=1, else flag=0
        if (do_matrix) then
                if (.not. matrix_all_tracer(ind_tracer)%is_active) then
                        flag=0
                else
                        ipop_abs = matrix_all_tracer(ind_tracer)%pop_index
                        ipop_active_index = findloc(i_pop, ipop_abs, dim=1)
                        flag=1
                endif
        endif

    end subroutine query_tracer_in_pop

    !the following subroutines are used in aerosol.f90
    !------------------------------start of aerosol.f90 call-------------------------------------
    subroutine query_matrix_info(npop_active, nspcs)
        integer, intent(out) :: npop_active, nspcs
        npop_active = 0

        if (matrix_module_init) then
            if (do_matrix) then
                npop_active = pop_active
                nspcs = 5 
            else
                npop_active = 0
                nspcs = 5
            endif
        else !matrix module hasn't been initialized when calling the function
            call error_mesg('get matrix configuration','configuration not defined '//trim(matrix_configuration), fatal) 
        end if

    end subroutine query_matrix_info


    subroutine query_pop_number(is, ie, js, je, i_order, i_abs_pop, pop_number)
        integer, intent(in) :: is, ie, js, je !physical window
        integer, intent(in) :: i_order !index of activated population, e.g. debug, i_order =1
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = i_acc
        real, intent(inout) :: pop_number(:,:, :) !number concentration of the current population
        integer :: pop_index
        pop_index = i_pop(i_order)
        i_abs_pop = pop_index
        pop_number = max(0.0, matrix_all_tracer(matrix_all_pop(pop_index)%i_n)%value_in_matrix(is:ie,js:je,:)) !#/m3
    end subroutine


    subroutine query_pop_dg_dry(is, ie, js, je, i_order, i_abs_pop, pop_dg_dry)
        integer, intent(in) :: is, ie, js, je !physical window
        integer, intent(in) :: i_order !index of activated population, e.g. debug, i_order =1
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = i_acc
        real, intent(inout) :: pop_dg_dry(:,:, :) !number concentration of the current population
        integer :: pop_index
        pop_index = i_pop(i_order)
        i_abs_pop = pop_index
        pop_dg_dry = max(0.0, matrix_all_pop(pop_index)%dg_dry(is:ie,js:je,:)) !m
    end subroutine

    subroutine query_pop_kappa(is, ie, js, je, i_order, i_abs_pop, pop_kappa)
        integer, intent(in) :: is, ie, js, je !physical window
        integer, intent(in) :: i_order !index of activated population, e.g. debug, i_order =1
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = i_acc
        real, intent(out) :: pop_kappa(ie-is+1,je-js+1, k_grid_size) !number concentration of the current population
        integer :: pop_index
        integer(kind=i8_kind) :: chksum_val
        pop_index = i_pop(i_order)
        i_abs_pop = pop_index
        pop_kappa = max(0.0, matrix_all_pop(pop_index)%kappa_pop(is:ie,js:je,:)) !m
!!$omp barrier
!! Get all of the threads in sync
!!$omp single
!        if (mpp_pe() .eq. mpp_root_pe()) then
!                write(*,'(A)') 'DEBUG:: i_order  pop_index  pop_name       is   ie   js   je  k_grid_size'
!                write(*,'(A,I8,I11,A12,4I5,I12)') '       ', i_order, pop_index, trim(matrix_all_pop(pop_index)%name), &
!                                      is, ie, js, je, k_grid_size
!        end if
!
!        chksum_val = mpp_chksum(matrix_all_pop(pop_index)%kappa_pop)
!        if (mpp_pe() .eq. mpp_root_pe()) &
!                print *, "DEBUG:: kappa_pop sum", chksum_val
!        
!
!
!!$omp end single
!!$omp barrier
    end subroutine

    subroutine query_pop_mspcs(is, ie, js, je, nspcs,  i_order, i_abs_pop, pop_mspcs)
        integer, intent(in) :: is, ie, js, je, nspcs !physical window
        integer, intent(in) :: i_order !index of activated population, e.g. debug, i_order =1
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = i_acc
        real, intent(inout) :: pop_mspcs(ie-is+1,je-js+1, k_grid_size, 1, nspcs) !mass concentration of each species in the population, [ug/m3]
        integer :: pop_index
        pop_index = i_pop(i_order)
        i_abs_pop = pop_index
        pop_mspcs = 0.0
        ! nspcs = 5, for activation code, the specs in the order of sulf, bcar, ocar, dust, seas
        if (matrix_all_pop(pop_index)%i_msulf > 0) then
            pop_mspcs(:,:,:, 1, 1) = max(0.0, &
                matrix_all_tracer(matrix_all_pop(pop_index)%i_msulf)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%i_mbcar > 0) then
            pop_mspcs(:,:,:, 1, 2) = max(0.0, &
                matrix_all_tracer(matrix_all_pop(pop_index)%i_mbcar)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%i_mocar > 0) then
            pop_mspcs(:,:,:, 1, 3) = max(0.0, &
                matrix_all_tracer(matrix_all_pop(pop_index)%i_mocar)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%i_mdust > 0) then
            pop_mspcs(:,:,:, 1, 4) = max(0.0, &
                matrix_all_tracer(matrix_all_pop(pop_index)%i_mdust)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%i_mseas > 0) then
            pop_mspcs(:,:,:, 1, 5) = max(0.0, &
                matrix_all_tracer(matrix_all_pop(pop_index)%i_mseas)%value_in_matrix(is:ie,js:je,:))
        endif

    end subroutine

    subroutine query_pop_sigma(i_order, i_abs_pop, pop_sigma)
        integer, intent(in) :: i_order
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = i_acc
        real, intent(out) :: pop_sigma
        integer :: pop_index
        pop_index = i_pop(i_order)
        i_abs_pop = pop_index
        pop_sigma = matrix_all_pop(pop_index)%sigma
    end subroutine





    !---------------------------------end of aerosol.f90 call ------------------------------------

    subroutine query_matrix_pop(ind_tracer, flag, is, ie, js, je, pop_index, pop_number, pop_dg_wet, pop_dens_wet)
        integer, intent(in) :: ind_tracer !index of tracer from query method
        integer, intent(in) :: is,ie,js,je !note: kd = size(r,3)
        integer, intent(out) :: flag  !index of population, if active and belong to matrix tracer pop
        ! number tracer: return 2, mass_tarcer: return 1, otherwise return -1
        integer, intent(out), optional :: pop_index
        real, intent(out), optional :: pop_number(ie-is+1, je-js+1, k_grid_size)
        real, intent(out), optional :: pop_dg_wet(ie-is+1, je-js+1, k_grid_size)
        real, intent(out), optional :: pop_dens_wet(ie-is+1, je-js+1, k_grid_size) 
        integer, parameter :: flag_n = 2, flag_m = 1
        integer :: kd
        kd = k_grid_size
        flag = -1
        pop_index = -1
        if (matrix_all_tracer(ind_tracer)%is_active) then
            ! mass or number tracer
            if (lowercase(trim(matrix_all_tracer(ind_tracer)%type)) .eq. 'mass') then
                flag = flag_m
            elseif (lowercase(trim(matrix_all_tracer(ind_tracer)%type)) .eq. 'number') then 
                flag = flag_n
            endif
            if (lowercase(trim(matrix_all_tracer(ind_tracer)%pop)) .ne. 'ext') then
                pop_index = matrix_all_tracer(ind_tracer)%pop_index
                pop_number = matrix_all_tracer(matrix_all_pop(pop_index)%i_n)%value_in_matrix(is:ie,js:je,:) !#/m3
                pop_dg_wet = matrix_all_pop(pop_index)%dg_wet(is:ie,js:je,:) !unit m
                pop_dens_wet = matrix_all_pop(pop_index)%dens_wet(is:ie,js:je,:) !unit: kg/m3
            endif

        else
            flag= -1
        endif
    end subroutine





    subroutine set_config_pop(imatrix_configuration)
        character(len=*), intent(in):: imatrix_configuration
        !----------------------------------------------------------------
        !       --- assign configuration population information ----
        ! akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
        ! if population exist, assign index >= 1, otherwise -1
        !----------------------------------------------------------------
        if (imatrix_configuration .eq. "debug") then
            !debug version:
            !pop: acc <- (n,m), ext <- (m_h2o, m_h2so4)
            !npop = 2 !number of poplation in the configuration
            !assign population index under current configuration
            !integer :: i_akk = 1, i_acc = 2, i_dd1 = 3, i_dd2 = 4 ! index of population in matrix
            !integer :: i_ssa = 5, i_ssc = 6, i_oc1 = 7, i_oc2 = 8 ! if exit, index >= 1, otherwise = -1
            !integer :: i_bc1 = 9, i_bc2 = 10, i_mxa = 11, i_mxc = 12, i_ext = 13
            i_acc = 2
            i_ext = 13
            pop_active = 1 !neglect ext
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_acc/)
        elseif (imatrix_configuration .eq. "full_all_pop") then  
            i_akk = 1
            i_acc = 2
            i_dd1 = 3
            i_dd2 = 4
            i_ssa = 5
            i_ssc = 6
            i_oc1 = 7
            i_oc2 = 8
            i_bc1 = 9
            i_bc2 = 10
            i_mxa = 11
            i_mxc = 12
            i_ext = 13
            pop_active = 12
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)= (/i_akk, i_acc, i_dd1, i_dd2, i_ssa, i_ssc, i_oc1, i_oc2, i_bc1, i_bc2, i_mxa, i_mxc/) 
        elseif (imatrix_configuration .eq. "akk_acc") then
            i_akk = 1
            i_acc = 2
            i_ext = 13
            pop_active = 2
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_akk, i_acc/)
        elseif (imatrix_configuration .eq. "akk_acc_mxa") then
            i_akk = 1
            i_acc = 2
            i_mxa = 11
            i_ext = 13
            pop_active = 3
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_akk, i_acc, i_mxa/)
        elseif (imatrix_configuration .eq. "oc_bc") then
            i_oc1 = 7
            i_oc2 = 8
            i_bc1 = 9
            i_bc2 = 10
            i_mxa = 11
            i_ext = 13
            pop_active = 5
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_oc1, i_oc2, i_bc1, i_bc2, i_mxa/)
        elseif (imatrix_configuration .eq. "oc_coag") then
            i_oc1 = 7
            i_oc2 = 8
            i_ext = 13
            pop_active = 2
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_oc1, i_oc2/)
        elseif (imatrix_configuration .eq. "bc_coag") then
            i_bc1 = 9
            i_bc2 = 10
            i_ext = 13
            pop_active = 2
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_bc1, i_bc2/)
        elseif (imatrix_configuration .eq. "ssa_ssc") then
            i_ssa = 5
            i_ssc = 6
            i_ext = 13
            pop_active = 2
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_ssa, i_ssc/)
        elseif (imatrix_configuration .eq. "dd1_dd2") then
            i_dd1 = 3
            i_dd2 = 4
            i_ext = 13
            pop_active = 2
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_dd1, i_dd2/)
        elseif (imatrix_configuration .eq. "debug_npf") then
            i_akk = 1
            pop_active = 1
            if (.not. allocated(i_pop)) &
                allocate(i_pop(pop_active))
            i_pop(:)=(/i_akk/)
        else
            call error_mesg('get matrix configuration','configuration not defined '//trim(imatrix_configuration), fatal)
            !call error_mesg ('tracer_driver', 'mw needs to be defined for tracer: '//trim(tracer_name), fatal)
        end if
    end subroutine













    real(8) function getxnum(ni,nj,dgni,dgnj,xlsgi,xlsgj)
        !---------------------------------------------------------------------------------------------------------------------
        ! dlw, 102306: derived from function getaf of cmaq v4.4.
        !
        ! getxnum = ln( dij / dgi ) / ( sqrt(2) * ln(sgi) ), where
        !
        !      dij is the diameter of intersection,
        !      dgi is the median diameter of the smaller size mode, and
        !      sgi is the geometric standard deviation of smaller mode.
        !
        ! a quadratic equation is solved to obtain getxnum, following the method of press et al. 1992.
        !
        ! references:
        !
        !  1. binkowski, f.s. and s.j. roselle, models-3 community multiscale air quality (cmaq)
        !     model aerosol component 1: model description.  j. geophys. res., vol 108, no d6, 4183
        !     doi:10.1029/2001jd001409, 2003.
        !  2. press, w.h., s.a. teukolsky, w.t. vetterling, and b.p. flannery, numerical recipes in
        !     fortran 77 - 2nd edition. cambridge university press, 1992.
        !----------------------------------------------------------------------------------------------------------------------
        implicit none

        ! arguments.

        real(8) :: ni         ! aitken       mode number concentration [#/m^3]
        real(8) :: nj         ! accumulation mode number concentration [#/m^3]
        real(8) :: dgni       ! aitken       mode geo. mean diameter [um]
        real(8) :: dgnj       ! accumulation mode geo. mean diameter [um]
        real(8) :: xlsgi      ! aitken       mode ln(geo. std. dev.) [1]
        real(8) :: xlsgj      ! accumulation mode ln(geo. std. dev.) [1]

        ! local variables.

        real(8) :: aa, bb, cc, disc, qq, alfa, l, yji
        real(8), parameter :: sqrt2 = 1.414213562d+00

        alfa = xlsgi / xlsgj
        yji = log( dgnj / dgni ) / ( sqrt2 * xlsgi )
        l = log( alfa * nj / ni)

        ! calculate quadratic equation coefficients & discriminant.
        aa = 1.0d+00 - alfa * alfa
        bb = 2.0d+00 * yji * alfa * alfa
        cc = l - yji * yji * alfa * alfa
        disc = bb * bb - 4.0d+00 * aa * cc

        ! if roots are imaginary, return a negative getaf value so that no imtr takes place.

        if( disc .lt. 0.0d+00 ) then
            getxnum = - 5.0d+00         ! error in intersection
            return
        endif

        ! equation 5.6.4 of press et al. 1992.

        qq = -0.5d+00 * ( bb + sign( 1.0d+00, bb ) * sqrt(disc) )

        ! return solution of the quadratic equation that corresponds to a
        ! diameter of intersection lying between the median diameters of the 2 modes.

        getxnum = cc / qq       ! see equation 5.6.5 of press et al.

        ! write(*,*)'getxnum = ', getxnum
        return
    end function getxnum
    !----------------------------------------------------------------
    !              subroutine get_matrix_dens 
    ! make it public
    ! used in wet deposition, convert #/kg/s to mol/m2/s -> kg/m2/s
    ! for test reason as well: the kg/m2/s of the number tracer should
    ! equals to the total wet deposition of mass tracers in the population
    !subroutine get_matrix_mass_per_particle(n, m_particle)
    !        integer, intent(in): 
    !
    !end subroutine

    ! calculates the vertical velocity of dust settling
    elemental real function sedimentation_velocity(t,p,rwet,rho_wet_dust) result(vdep)
        real, intent(in) :: t            ! air temperature, deg k
        real, intent(in) :: p            ! pressure, pa
        real, intent(in) :: rwet         ! radius of dust particles, m
        real, intent(in) :: rho_wet_dust ! density of dust particles, kg/m3

        real :: viscosity, free_path, c_c
        viscosity = 1.458e-6 * t**1.5/(t+110.4)     ! dynamic viscosity
        free_path = 6.6e-8*t/293.15*(pstd_mks/p)
        c_c = 1.0 + free_path/rwet * &              ! slip correction [none]
            (1.257+0.4*exp(-1.1*rwet/free_path))
        vdep = 2./9.*c_c*grav*rho_wet_dust*rwet**2/viscosity  ! settling velocity [m/s]
    end function sedimentation_velocity
    subroutine sedimentation_flux(sj_scheme,kb,dt,mtv,dz,vdep,air_dens,&
            pwt,tracer,tracer_dt,setl)
        implicit none
        integer, intent(in) :: kb !< index for the buttom layer
        real, intent(in) :: dt, mtv !< model timestep , mixing ratio factor conversion
        logical, intent(in) :: sj_scheme !< .true. if using sj's scheme
        real, intent(in),     dimension(:) :: dz, air_dens, pwt !<layer thickness, density, mass per unit area
        real, intent(in),     dimension(:) :: vdep      !< settling velocity [m/s]
        real, intent(in),     dimension(:) :: tracer    !< tracer concentration
        real, intent(inout),  dimension(:) :: tracer_dt !< tracer tendency
        real, intent(inout),  dimension(:) :: setl !< tracer settling flux [kg/m2/s]
        !    local vars
        real, dimension(size(dz)) :: qn,qn1
        integer  k

        setl(:)=0.
        tracer_dt = 0.
        qn(:)=tracer(:)
        if (sj_scheme) then
            qn1(1)=qn(1)*dz(1)/(dz(1)+dt*vdep(1))
            do k=2,kb
            qn1(k)=(qn(k)*dz(k)+dt*qn1(k-1)*vdep(k-1)*air_dens(k-1)/air_dens(k))/(dz(k)+dt*vdep(k))
            enddo
            tracer_dt(:)=tracer_dt(:)+(qn1(:)-qn(:))/dt
            setl(kb) = qn1(kb)*air_dens(kb)/mtv*vdep(kb)
        else
            do k=1,kb
            if (tracer(k) > 0.0) then
                setl(k)=tracer(k)*air_dens(k)/mtv*vdep(k)    ! settling flux [kg/m2/s]
            endif
            enddo
            tracer_dt(1)=tracer_dt(1)-setl(1)/pwt(1)*mtv
            tracer_dt(2:kb)=tracer_dt(2:kb) &
                + ( setl(1:kb-1) - setl(2:kb) )/pwt(2:kb)*mtv
        endif

        !if (mpp_root_pe().eq.mpp_pe()) then
        !        write(*,*), 'DEBUG: tracer_dt and setl values'
        !        write(*,*), 'tracer_dt colume sum', sum(tracer_dt*air_dens*dz, dim=1) !mmr to kg/m2 
        !        write(*,*), 'setl values:', sum(setl, dim=1)
        !        write(*,*), 'setl kb:', setl(kb)
        !endif
    end subroutine sedimentation_flux

end module matrix_gfdl
