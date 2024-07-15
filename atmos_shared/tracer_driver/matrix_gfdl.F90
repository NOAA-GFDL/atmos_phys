module matrix_gfdl

    use mpp_mod,               only : input_nml_file
    use fms_mod,               only : file_exist, close_file,&
        open_namelist_file, check_nml_error, &
        write_version_number, &
        error_mesg, &
        FATAL, NOTE, &
        lowercase, &
        mpp_pe, &
        mpp_root_pe, &
        stdlog, stdout, &
        mpp_clock_id, &
        mpp_clock_begin, &
        mpp_clock_end, &
        CLOCK_MODULE, &
        uppercase
    use tracer_manager_mod,    only : get_tracer_index,   &
        get_number_tracers, &
        get_tracer_names,   &
        get_tracer_indices, &
        adjust_positive_def, &
        query_method, &
        NO_TRACER
    use  field_manager_mod, only : MODEL_ATMOS, &
        parse
    use  diag_manager_mod, only  : send_data, register_diag_field
    use  time_manager_mod, only  : time_type
    use  constants_mod, only     : PI, GRAV, RDGAS, DENS_H2O, WTMAIR, AVOGNO, &
                                   PSTD_MKS
    use sat_vapor_pres_mod, only : compute_qs  
    use aero_npf, only : NPFRATE, STEADY_STATE_H2SO4
    use aero_wet, only : aero_kohler
    use aero_condens, only : SETUP_KCI
    use aero_coag_config, only : SETUP_COAG_TENSORS, &
                                 GIKLQ_control, DIKL_control, &
                                 NM, PROD_INDEX, GIKLQ, &
                                 DIKL, DIJ, nDIKL,CITABLE, &
                                 MODE_NAME, NMASS_SPCS   
    use aero_coag,  only : GET_KNIJ 
    implicit none
    private

    public matrix_init, matrix_source_type,set_matrix_source, matrix_run
    public query_matrix_pop, query_matrix_info
    public query_pop_number, query_pop_Dg_dry, query_pop_MSPCS, query_pop_SIGMA

    interface set_matrix_source
        module procedure set_matrix_source_2d
        module procedure set_matrix_source_3d
    end interface set_matrix_source
     
    logical :: matrix_module_init = .false.
    real, parameter :: PI6 = PI/6.0
    integer, parameter :: nb_tracer_max = 9 !maxium number of tracers in a population, i.e. in MXX
    integer :: kd !kd = size(r,3), defined in matrix_init,3rd dimension of emission sources
    integer :: id_RH,id_kc, id_dmdt_h2so4_tot_cond_npf, id_dndt_npf, id_dmdt_h2so4_npf,id_cond_sink
    integer :: id_h2so4_emis, id_pwt, id_zhalf, id_so4_emis !id for budget analysis 
    integer :: id_h2so4_col
    !define clocks for different processes
    integer :: ini_clock = 0 !initializetion
    integer :: npf_clock = 0 !new particle formation clock
    integer :: condgrow_clock = 0 !condesational growth
    integer :: hygrow_clock = 0 !hygroscopic growth
    integer :: coag_clock = 0 ! coagulation growth
    integer :: dmodal_clock = 0 !inter-modal transfer
    integer :: aqso4_clock = 0 ! partitioning clock 
    ! matrix source type
    type :: matrix_source_struct
        integer :: P_H2SO4 = 1 !gaseous phase H2SO4 to be used in nucleation scheme in AKK
        integer :: E_SO4   = 2 !SO2 direct emission: 2.5% sulfur emission to be converted to sulfate mass in ACC
        integer :: E_DUST  = 3 !dust emission to be distributed in DD1 and DD2
        integer :: E_SS    = 4 !sea salt emission to be distributed in SSA and SSC
        integer :: E_OC    = 5 !organic carbon emission to be distributed in OC1 and OC2
        integer :: E_BC    = 6 !black carbon emission to be distributed in BC1 and BC2     
        integer :: E_SOA   = 7 !extra emission based on GFDL model
        integer :: P_AQSO4 = 8 !aqueous so4 production
        integer :: U_KG_M2_S = 1 !unit index for conversion
        integer :: U_VMR_S   = 2 !unit index for conversion
        integer :: U_MMR_S   = 3 !unit index for conversion
    end type matrix_source_struct


    type(matrix_source_struct) :: matrix_source_type

    !----------------------------------------------------------------
    !	-------------------type matrix_tracer-------------------
    !   ...define tracer type and encorporate tracer properties
    ! 	...tracer description in field table:
    !		"matrix_parameters", "acc", "distribution=lognormal, &
    !             & distribution_index=1,
    !		& sigma=1.8, Dgn=0.068e-6, kappa=0.507, &
    !		&dens=1770, type=mass,spec=sulf"
    !----------------------------------------------------------------
    type :: matrix_tracer
        logical :: is_active = .FALSE.
        character*32 :: type = "" ! mass or number
        character*3  :: pop = "" ! population name: akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
        integer :: pop_index = -1
        character*32 :: spec = ""! mass species name: sulf/dust/seas/ocar/bcar/alwc/gsfa...
        character*32 :: distribution = "" ! lognormal or weibull
        integer      :: distribution_index = -1
        real :: sigma   = -1 ! sigma of log-normal distribution
        real :: lnsigma  = -1 
        real :: lnsigma2 = -1 
        real :: Dp0 = -1 !volume mean diameter of the current mode, unit m
        real :: Dgn = -1 ! unit m: geometric diameter for initial log-normal distribution, only used for emission
        real :: kappa = 0 ! hygroscopicity facor
        real :: dens  = 0! unit Kg/m3: density of specific species
        character*32 :: name,units !get_tracer_names (MODEL_ATMOS, tracer_index, name= MT%name,  units = MT%units)
        real, allocatable    :: source(:,:,:) !sources in unit of µg/m3/s (mass tracers) or #/m3/s
        real, allocatable    :: value_in_matrix(:,:,:) !calculate tracer values in matrix unit
        logical :: has_emission = .FALSE. !if the field table has distribution parameters, then has emission, otherwise, no-emission
        integer :: id_tracer_emis = -1 !diagostic emission IDs, 3D µg/m3/s
        integer :: id_tracer_col = -1 !dry deposition IDs, 2D kg/m2/s
    end type matrix_tracer

    !----------------------------------------------------------------
    !	-------------------type matrix_pop-------------------
    !   ...define population type (contain multi tracers)
    !	... to record what tracers are within a specific population
    !   ... if a tracer is in this population, record its index in the tracer array
    !   ...typical components in a population: N, M_spec1, M_spec2, ....
    !----------------------------------------------------------------
    type :: matrix_pop
        integer :: nb_tracer_pop = 0 !number of tracers in current population
        character*32 :: name = "" !name of the population
        integer :: tracer_index(nb_tracer_max) = -1
        logical :: has_emission(nb_tracer_max) = .FALSE.
        integer :: I_N = -1 ! index of number in matrix_tracer array
        integer :: I_MSULF = -1 ! index of sulfate mass in matrix_tracer array
        integer :: I_MDUST = -1
        integer :: I_MSEAS = -1
        integer :: I_MOCAR = -1
        integer :: I_MBCAR = -1
        integer :: I_MWATE = -1
        integer :: I_MAMMO = -1
        integer :: I_MNITR = -1
        integer :: I_MGSFA = -1
        real :: sigma = 1.8 !gemetric sigma for population, temperary set for 1.8
        real :: rh_deliquescence = 0
        real :: rh_crystallization = 0
        real, allocatable    :: Dg_dry(:,:,:) !unit: m, geometric diameter
        real, allocatable    :: Dg_wet(:,:,:) !unit: m, geometric diameter
        real, allocatable    :: mass_dry(:,:,:) !unit: ug/m3
        real, allocatable    :: dens_wet(:,:,:) !unit: kg/m3
        real, allocatable    :: dens_dry(:,:,:) !unit: kg/m3
        real, allocatable    :: vol_dry(:,:,:) 
        real, allocatable    :: kappa_pop(:,:,:) !average hygroscopicity of the population
        integer :: id_Dg_dry = -1
        integer :: id_Dg_wet = -1
        integer :: id_pop_condens = -1
        integer :: id_pop_coag_pN = -1 !coag production of number [#/m3/s]
        integer :: id_pop_coag_lN = -1 !coag loss of number, 3D [#/m3/s]
        integer :: id_pop_coag_pMSULF = -1 !coag production of sulf mass [ug/m3/s]
        integer :: id_pop_coag_lMSULF = -1 !coag loss of sulf mass [ug/m3/s]

    end type matrix_pop

    type(matrix_tracer), allocatable :: matrix_all_tracer(:) ! define tracer array for all tracers
    type(matrix_pop),    allocatable :: matrix_all_pop(:)  ! defined population array for all population
    integer :: npop = 13 !maximum number of population in matrix
    !integer :: ntracer !maximum number of population and tracers in matrix
    integer :: I_AKK = -1, I_ACC = -1, I_DD1 = -1, I_DD2 = -1 ! index of population in matrix
    integer :: I_SSA = -1, I_SSC = -1, I_OC1 = -1, I_OC2 = -1 ! if exit, index >= 1, otherwise = -1
    integer :: I_BC1 = -1, I_BC2 = -1, I_MXA = -1, I_MXC =-1, I_EXT = -1
    integer :: pop_active = 0 ! number of active population
    integer, allocatable :: I_POP(:)
    integer :: I_MW_H2SO4 = 1, I_MW_SO4 = 2, I_MW_DUST = 3, I_MW_SS = 4
    integer :: I_MW_OC = 5, I_MW_BC = 6, I_MW_SOA = 7
    integer :: NMSPCS = 5 ! number of mass species
    real, dimension(7) :: matrix_molecular_weight = [98.07848, 96.0, 135.0, 58.5, 12.0, 12.0, 12.0]
    ! Assign values to each element using the indices
    !matrix_molecular_weight(I_MW_H2SO4) = 98.07848
    !matrix_molecular_weight(I_MW_SO4) = 96.0
    !matrix_molecular_weight(I_MW_DUST) = 135.0
    !matrix_molecular_weight(I_MW_SS) = 58.5
    !matrix_molecular_weight(I_MW_OC) = 12.0
    !matrix_molecular_weight(I_MW_BC) = 12.0
    !matrix_molecular_weight(I_MW_SOA) = 12.0
    real, allocatable :: AQSO4(:,:,:) !passing in: aqueous SO4 production
    real, allocatable :: P_H2SO4_RATE(:,:,:)
    integer, parameter :: DIST_LOGNORMAL = 1, DIST_WEIBULL = 2 !distribution index for lognormal and weibull
    integer :: tracer_index,  n, ntrace, nt, ntt
    integer :: ierr, io, logunit, verbose, unit
    logical :: flag
    character(len=512) :: text_in_scheme, control

    integer :: nsphum, nh2so4
    logical :: do_matrix = .FALSE., do_coag = .FALSE., do_intermodal_transfer = .FALSE.
    logical :: do_AQSO4 = .FALSE.
    character(len=32) :: matrix_configuration = 'debug'
    character(len=32) :: coag_configuration = 'full'
    character(len=32) :: intermodal_configuration = 'full'
    character(len=7), parameter :: module_name = 'matrix'
    ! matrix_configuration is the version used in the calculation
    ! curretly matrix_configuration = debug is used to test toy model
    namelist /matrix_nml / do_matrix, matrix_configuration, & 
                           do_coag, coag_configuration, &
                           do_intermodal_transfer, intermodal_configuration, &
                           do_AQSO4 
contains

    !----------------------------------------------------------------
    !
    !                           subroutine matrix_run
    !IMPORTANT: this subroutine must be after atmos_SOx_chem etc, where sources are properly set-up 
    !the mass emission of different matrix_tracers are linked to different sources
    !in other F90 files, e.g. in atmos_SOx_chem of atmos_sulfate.F90
    ! The source will be initialized by 0, and updated at every timestep 
    !----------------------------------------------------------------
    subroutine matrix_run(r, pfull, rh, t, dt, pwt, zhalf, rdt_matrix, Time, is,ie,js,je)
            real, intent(in) :: r(:,:,:,:)
            real, intent(in) :: dt !timestep
            integer, intent(in) :: is, ie, js, je ! boundaries of physical window
            type(time_type),  intent(in) :: Time
            real, intent(out) :: rdt_matrix(:,:,:,:) !tendency calculated from matrix
            real, intent(in) :: pfull(:,:,:), rh(:,:,:), t(:,:,:), pwt(:,:,:), zhalf(:,:,:) !t is temperature
            real :: XNH3(size(r,1),size(r,2),size(r,3)), FLAND(size(r,1),size(r,2),size(r,3)) !XL
            ! assign number rate to each number tracers, the mass emission has been linked to different tracers in other files
            real :: rt(size(r,1),size(r,2),size(r,3),size(r,4)) !local variable to record updated values of tracers, in gfdl unit
            !local variable to record absolute h2so4 mass change in ug/m3, positive
            real :: dm_h2so4(size(r,1),size(r,2),size(r,3)), dmdt_npf(size(r,1),size(r,2),size(r,3)), dndt_npf(size(r,1),size(r,2),size(r,3))
            real :: dmdt_h2so4_tot_cond_npf(size(r,1),size(r,2),size(r,3)), dmdt_h2so4_npf(size(r,1),size(r,2),size(r,3)) 
            integer :: n,MW, npf_flag_npf !local variables
            real :: m_akk_emis_source(size(r,1),size(r,2),size(r,3)), n_akk_emis_source(size(r,1),size(r,2),size(r,3))
            real :: kci_coef_pop(npop, size(r,1),size(r,2),size(r,3)), kci_aeq1_pop(npop, size(r,1),size(r,2),size(r,3)) !unit: m3/s
            logical :: used
            real :: XH2SO4_NUCL(size(r,1),size(r,2),size(r,3))
            integer :: IH2SO4_PATH(size(r,1),size(r,2),size(r,3))
            real:: KC(size(pfull,1),size(pfull,2),size(pfull,3)) !KC: total condensation sink (1/s) for all population in a grid
            real:: PQ_GROWTH(size(pfull,1),size(pfull,2),size(pfull,3))
            real:: number_pop(size(pfull,1),size(pfull,2),size(pfull,3))
            real:: cond_pop(npop,  size(r,1),size(r,2),size(r,3))
            real:: emis_tmp(size(pfull,1),size(pfull,2)), col_tmp(size(pfull,1),size(pfull,2)), cond_tmp(size(pfull,1),size(pfull,2))
            integer :: it,jt,kt,i,j,k, nct !for loop use
            real :: LI(npop-1,size(r,1),size(r,2),size(r,3)) !loss term of number due to intermodal coagulation [#/m^3/s]
            real :: RI(npop-1,size(r,1),size(r,2),size(r,3)) !production terms due to intermodal coagulation. [#/m^3/s]
            real :: LIM(npop-1,NMSPCS,size(r,1),size(r,2),size(r,3)) !mass loss term of number due to intermodal coagulation [ug/m^3/s]
            real :: RIM(npop-1,NMSPCS,size(r,1),size(r,2),size(r,3)) !mass production terms due to intermodal coagulation. [ug/m^3/s]
            real :: total_aero_num(size(r,1),size(r,2),size(r,3))
            real :: kci1, kci2, fac !debug use
            real :: pop_vset(npop-1, size(r,1),size(r,2),size(r,3))!settling velocity for population
            it = size(pfull,1)
            jt = size(pfull,2)
            kt = size(pfull,3)
            dm_h2so4 = 0.
            rdt_matrix = 0.            
            dndt_npf = 0.
            dmdt_npf = 0.
            XNH3 = 0.
            FLAND = 0.
            IH2SO4_PATH = 0
            KC = 0. 
            if (do_matrix) then
                    rt = r
                    !------------------------------------------------------------------
                    !               Step 0: get  matrix tracer values 
                    !               and convert from gfdl to matrix unit
                    !------------------------------------------------------------------
                    call mpp_clock_begin(ini_clock)
                    !initialize matrix species values at the beginning of the current step
                    call set_matrix_value(rt, pwt, zhalf,is,ie,js,je) ! assign matrix_tarcer%values_in_matrix, value in matrix unit
                    !-----------------------------------------------------------------------------
                    !               Step 1: calculate tracer number sources from direct emission (in matrix unit)
                    !                       Note: this step doesn't include new particle formation
                    !    Note: the mass sources have been already updated before matrix_run
                    !------------------------------------------------------------------------------
                    call set_matrix_emis_number(is,ie,js,je) !update tracer source of  num_rate 
                    call mpp_clock_end(ini_clock)
                    !------------------------------------------------------------------------------------
                    !               Step *: H2SO4 condensational growth, must prior new particle formation
                    !(1) calculate the condensation coefficient for all population over all grids
                    !           kci_coef_pop 4-D dimensions: (npop, is, ij, ik)
                    !           kci_coef_aeq1 4-D dimensions: (npop, is, ij, ik)
                    !-------------------------------------------------------------------------------------
                    call mpp_clock_begin(condgrow_clock)
                    call set_matrix_pop_kci(pfull,t,kci_coef_pop, kci_aeq1_pop,is,ie,js,je) !calculate the condensation coefficient for all population
                    !---------------------------------------------------------------------------------
                    !               Step *: condensational growth w/wo new particle formation of AKK
                    !               calculate new particle formation rate and update sources of dndt/dmdt 
                    !XLXLXLXLIMPORTANT: H2SO4 need to be updated by the loss of H2SO4 condensational loss
                    !-----------------------------------------------------------------------------------
                    call mpp_clock_begin(npf_clock)
                    call set_matrix_condense_npf(I_AKK, XNH3, FLAND, pfull,rh,t, dt, r(:,:,:,nh2so4), pwt,zhalf, kci_coef_pop, kci_aeq1_pop, &
                            dndt_npf, dmdt_h2so4_npf, IH2SO4_PATH, dmdt_h2so4_tot_cond_npf, KC, XH2SO4_NUCL, &
                            is, ie, js, je) !the I_MSULF term already updated for each population
                    call mpp_clock_end(npf_clock)
                    !update sources
                    !assumption: H2SO4(g, 98) -> aerosol(p, 96), H2SO4 condensed on aerosol surface would become sulfate with molecular
                    !weight 96, scale condensational H2SO4 mass -> aerosol sulfate mass, unit: ug/m3/s
                    PQ_GROWTH = (dmdt_h2so4_tot_cond_npf - dmdt_h2so4_npf)*96.0/98.0
                    !scale new particle formation H2SO4 mass -> aerosol sulfate mass, unit: ug/m3/s
                    dmdt_npf = dmdt_h2so4_npf*96.0/98.0
                    !-----------------------------------------------------------------------------------
                    !                send emis data before source get changed with different processes
                    !-----------------------------------------------------------------------------------
                    do n=1, ntrace
                    if (matrix_all_tracer(n)%id_tracer_emis > 0) then
                            emis_tmp = 0.
                            if (matrix_all_tracer(n)%type .eq. 'mass') then
                                    do nct = 1, kd !ug/m3/s -> kg/m2/s
                                    emis_tmp = emis_tmp + matrix_all_tracer(n)%source(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))*1E-9
                                    enddo
                            elseif(matrix_all_tracer(n)%type .eq. 'number') then
                                    do nct = 1, kd !#/m3/s -> mol/m2/s
                                    emis_tmp = emis_tmp + &
                                            matrix_all_tracer(n)%source(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))/AVOGNO
                                    enddo
                            endif
                            used = send_data(matrix_all_tracer(n)%id_tracer_emis, emis_tmp, time, &
                                    is_in=is,js_in=js)
                    endif
                    enddo

                    !-----------------------------------------------------------------------------------------------------------
                    !               update matrix sources due to new particle formation and condensational growth
                    !-----------------------------------------------------------------------------------------------------------            
                    cond_pop = 0.
                    do n=1,npop
                    if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
                            if ((I_AKK > 0) .AND. (n .EQ. I_AKK)) then
                                    n_akk_emis_source = matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source(is:ie,js:je,:)
                                    m_akk_emis_source = matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source(is:ie,js:je,:) !record ini source from emission
                                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source(is:ie,js:je,:) = n_akk_emis_source+dndt_npf
                                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source(is:ie,js:je,:) = m_akk_emis_source+dmdt_npf
                            endif
                            number_pop = matrix_all_tracer(matrix_all_pop(n)%I_N)%value_in_matrix(is:ie,js:je,:) !number concentration of a population
                            do i=1,it
                            do j=1,jt
                            do k=1,kt
                            if (KC(i,j,k) > 0) then
                                    cond_pop (n,i,j,k) = (kci_coef_pop(n,i,j,k)*number_pop(i,j,k)/KC(i,j,k) ) * PQ_GROWTH(i,j,k)
                                    matrix_all_tracer(matrix_all_pop(n)%I_MSULF)%source(i+is-1,j+js-1,k) = &
                                            matrix_all_tracer(matrix_all_pop(n)%I_MSULF)%source(i+is-1,j+js-1,k) + cond_pop (n,i,j,k)
                            endif

                            enddo
                            enddo
                            enddo

                    endif                         
                    enddo
                    !update H2SO4 loss
                    dm_h2so4 = dm_h2so4 + dmdt_h2so4_tot_cond_npf * dt !in matrix unit ug/m3
                    call mpp_clock_end(condgrow_clock)

                    !-----------------------------------------------------------------------------------
                    !                send condensation data before source get changed with different processes
                    !-----------------------------------------------------------------------------------
                    do n=1, npop
                    if (matrix_all_pop(n)%id_pop_condens > 0) then
                            cond_tmp = 0.
                            do nct = 1, kd !ug/m3/s -> kg/m2/s
                            cond_tmp = cond_tmp + cond_pop(n,:,:,nct)* (zhalf(:,:,nct)-zhalf(:,:,nct+1))*1E-9
                            enddo
                            used = send_data(matrix_all_pop(n)%id_pop_condens, cond_tmp, time, &
                                    is_in=is,js_in=js)
                    endif
                    enddo
                    
                    if (do_coag) then
                            call mpp_clock_begin(coag_clock)
                            !-----------------------------------------------------------------------------
                            !               Step *: Coagulation
                            ! RI/LI: (npop-1, size(r,1),size(r,2),size(r,3)), production/loss of number in mode I [#/m3/s]
                            ! RIM/LIM: (npop-1,NMSPCS,size(r,1),size(r,2),size(r,3)), production/loss of species qq mass in mode I [ug/m3/s]
                            !------------------------------------------------------------------------------
                            call matrix_coag(rt,pfull, t, RI, LI, RIM, LIM, is,ie,js,je) !LI/LIM positive number

                            !-----------------------------------------------------------------------------------------------------------
                            !               update matrix sources due to coagulation and send the coagulation source
                            !-----------------------------------------------------------------------------------------------------------
                            do n=1,npop
                            if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
                                    matrix_all_tracer(matrix_all_pop(n)%I_N)%source = matrix_all_tracer(matrix_all_pop(n)%I_N)%source &
                                            + RI(n,:,:,:) - LI(n,:,:,:) !number source [#/m3/s]
                                    used = send_data(matrix_all_pop(n)%id_pop_coag_pN, RI(n,:,:,:), time, &
                                            is_in=is,js_in=js, ks_in=1)
                                    used = send_data(matrix_all_pop(n)%id_pop_coag_lN, LI(n,:,:,:), time, &
                                            is_in=is,js_in=js, ks_in=1)
                                    if (matrix_all_pop(n)%I_MSULF > 0) then !mass source [ug/m3/s]
                                            matrix_all_tracer(matrix_all_pop(n)%I_MSULF)%source = &
                                                    matrix_all_tracer(matrix_all_pop(n)%I_MSULF)%source + RIM(n,1,:,:,:) - LIM(n,1,:,:,:)
                                            used = send_data(matrix_all_pop(n)%id_pop_coag_pMSULF, RIM(n,1,:,:,:), time, &
                                                    is_in=is,js_in=js, ks_in=1)
                                            used = send_data(matrix_all_pop(n)%id_pop_coag_lMSULF, LIM(n,1,:,:,:), time, &
                                                    is_in=is,js_in=js, ks_in=1)
                                    elseif (matrix_all_pop(n)%I_MDUST > 0) then
                                            matrix_all_tracer(matrix_all_pop(n)%I_MDUST)%source = &
                                                    matrix_all_tracer(matrix_all_pop(n)%I_MDUST)%source + RIM(n,2,:,:,:) - LIM(n,2,:,:,:)
                                    elseif (matrix_all_pop(n)%I_MSEAS > 0) then
                                            matrix_all_tracer(matrix_all_pop(n)%I_MSEAS)%source = &
                                                    matrix_all_tracer(matrix_all_pop(n)%I_MSEAS)%source + RIM(n,3,:,:,:) - LIM(n,3,:,:,:)
                                    elseif (matrix_all_pop(n)%I_MOCAR > 0) then
                                            matrix_all_tracer(matrix_all_pop(n)%I_MOCAR)%source = &
                                                    matrix_all_tracer(matrix_all_pop(n)%I_MOCAR)%source + RIM(n,4,:,:,:) - LIM(n,4,:,:,:)
                                    elseif (matrix_all_pop(n)%I_MBCAR > 0) then
                                            matrix_all_tracer(matrix_all_pop(n)%I_MBCAR)%source = &
                                                    matrix_all_tracer(matrix_all_pop(n)%I_MBCAR)%source + RIM(n,5,:,:,:) - LIM(n,5,:,:,:)
                                    endif
                            endif
                            enddo
                            call mpp_clock_end(coag_clock)

                    endif
                    
                    !-------------------------------------------------------------------------------
                    !               Step *: Partition AQSO4 to matrix: add to matrix sources
                    !   Currently, don't calculate activation, just distribute AQSO4 based on numbers 
                    !---------------------------------------------------------------------------------
                    if (do_AQSO4) then
                            call mpp_clock_begin(aqso4_clock)
                            !get total number of all pops (exclude EXT and AKK)
                            total_aero_num = 0.
                            do n=1,npop
                            if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
                                 .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk")) then !exclude ext and akk
                                    total_aero_num = total_aero_num + matrix_all_tracer(matrix_all_pop(n)%I_N)%value_in_matrix(is:ie,js:je,:) ![#/m3]
                            endif
                            enddo
                            total_aero_num = max(total_aero_num, 1.0E-32)
                            !add AQSO4 source
                            do n=1,npop
                            if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
                                  .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk")) then !exclude ext and akk
                                    matrix_all_tracer(matrix_all_pop(n)%I_MSULF)%source(is:ie,js:je,:) = &
                                            matrix_all_tracer(matrix_all_pop(n)%I_MSULF)%source(is:ie,js:je,:) + &
                                    AQSO4(is:ie,js:je,:)*matrix_all_tracer(matrix_all_pop(n)%I_N)%value_in_matrix(is:ie,js:je,:)/total_aero_num ![ug/m3/s]
                            endif
                            enddo
                            call mpp_clock_end(aqso4_clock)
                    endif 


                    !-----------------------------------------------------------------------------
                    !               Step *: Hygroscopic growth: calculate dry/wet particle diameter
                    !------------------------------------------------------------------------------
                    call mpp_clock_begin(hygrow_clock)
                    call matrix_dry_diameter(rt,is,ie,js,je) !assign dry pop properties 
                    call matrix_wet_diameter(rt, rh, t,is,ie,js,je) !assign wet pop properties, t is temperature
                    call mpp_clock_end(hygrow_clock)

                    !-----------------------------------------------------------------------------
                    !               Step *: intermodal transfer + update matrix value
                    !   pseudo update matrix_tracer values -> intermodal subroutine: update value again
                    !---------------------------------------------------------------------------------------------------- 
                    do n=1,ntrace
                    if (matrix_all_tracer(n)%is_active) then
                            matrix_all_tracer(n)%value_in_matrix = matrix_all_tracer(n)%value_in_matrix + matrix_all_tracer(n)%source * dt
                            matrix_all_tracer(n)%value_in_matrix = max(matrix_all_tracer(n)%value_in_matrix, 0.0) !value can't < 0
                    endif
                    enddo

                    if (do_intermodal_transfer) then
                            call mpp_clock_begin(dmodal_clock)
                            !note: matrix value got updated in this matrix_intermodal_transfer subroutine
                            call matrix_intermodal_transfer(intermodal_configuration, rt, is,ie,js,je)
                            do n=1,ntrace
                            if (matrix_all_tracer(n)%is_active) then
                                    matrix_all_tracer(n)%value_in_matrix = max(matrix_all_tracer(n)%value_in_matrix, 0.0) !value can't < 0
                            endif
                            enddo
                            call mpp_clock_end(dmodal_clock)
                    endif
        
                    !-----------------------------------------------------------------------------
                    !               Step *: sedimentation settling by gravity
                    !  ->  calculate gravitational settling and update value again in the subroutine
                    !----------------------------------------------------------------------------------------------------
                    call matrix_sedimentation(T, pfull, zhalf, dt, is, ie, js, je)

                    !---------------------------------------------------------------------------------
                    !               Final step: pass matrix_tracer values and H2SO4 concentration to gfdl tendency
                    !---------------------------------------------------------------------------------
                    do n=1,ntrace
                    if (matrix_all_tracer(n)%is_active) then
                            call update_rt_from_matrix(n, rt(:,:,:,n), pwt, zhalf,is,ie,js,je) !convert tendency from matrix to gfdl unit, ug/m3-> mmr
                            rdt_matrix(:,:,:,n) = (rt(:,:,:,n)-r(:,:,:,n))/dt
                    elseif (n .eq. nh2so4) then
                            !do unit from matrix unit to gfdl unit: ug/m3 -> vmr
                            MW = 98
                            dm_h2so4 = dm_h2so4/1E9/(pwt(:,:,:) * MW / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))) 
                            !dm_h2so4 can't exceed previous h2so4 concentration
                            do i=1,it
                            do j=1,jt
                            do k=1,kt
                                dm_h2so4(i, j, k) = min(dm_h2so4(i, j, k), r(i, j, k, nh2so4))
                            enddo
                            enddo
                            enddo
                            rdt_matrix(:,:,:,n) = -dm_h2so4/dt !in gfdl unit
                    endif
                    enddo

                    
                    !--------------------------------------------------------------------------------------------------------------------
                    ! Need to reset AKK population source to 0. as it is not done automatically if there is no other model source but NPF
                    ! (set_matrix_source is not called)
                    !--------------------------------------------------------------------------------------------------------------------         
                    if (I_AKK > 0 ) then
                            matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source     = 0
                            matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source = 0
                    endif

                    !--------------------------------------------------------------------------------------------------------------------
                    ! register and send data
                    !--------------------------------------------------------------------------------------------------------------------
                    !register and send data of RH
                    used = send_data (id_RH, rh*100, time, &
                            is_in=is,js_in=js, ks_in = 1) !in unit %
                    used = send_data (id_cond_sink, PQ_GROWTH*dt, time, &
                            is_in=is,js_in=js, ks_in = 1) !in ug/m3
                    used = send_data (id_kc, KC, time, &
                            is_in=is,js_in=js, ks_in = 1) 
                    used = send_data (id_dmdt_h2so4_tot_cond_npf, dmdt_h2so4_tot_cond_npf, time, &
                            is_in=is,js_in=js, ks_in = 1) 
                    used = send_data (id_dndt_npf, dndt_npf, time, &
                            is_in=is,js_in=js, ks_in = 1) 
                    used = send_data (id_dmdt_h2so4_npf, dmdt_h2so4_npf, time, &
                            is_in=is,js_in=js, ks_in = 1) 
                    used = send_data (id_pwt, pwt, time, &
                            is_in=is,js_in=js, ks_in = 1)
                    used = send_data (id_zhalf, zhalf(:,:,1:kd)-zhalf(:,:,2:kd+1), time, &
                            is_in=is,js_in=js, ks_in = 1)


                    !send data of D_wet, D_dry for test: ps. Dg_dry, Dg_wet in matrix_pop unit is m, only for send_data change to um
                    do n = 1, npop
                    if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
                            if (matrix_all_pop(n)%id_Dg_dry > 0) then
                                    used = send_data (matrix_all_pop(n)%id_Dg_dry, matrix_all_pop(n)%Dg_dry(is:ie,js:je,:)*1E6, time, &
                                            is_in=is,js_in=js, ks_in = 1) !in unit µm
                            endif
                            if (matrix_all_pop(n)%id_Dg_wet > 0) then
                                    used = send_data (matrix_all_pop(n)%id_Dg_wet, matrix_all_pop(n)%Dg_wet(is:ie,js:je,:)*1E6, time, &
                                            is_in=is,js_in=js, ks_in = 1) !in unit µm
                            endif

                    endif
                    end do

                    !-----------------------------------------------------------------------------------
                    !                send col data for matrix tracers
                    !-----------------------------------------------------------------------------------
                    do n=1, ntrace
                    if (matrix_all_tracer(n)%id_tracer_col > 0) then
                            col_tmp = 0.
                            if (matrix_all_tracer(n)%type .eq. 'mass') then
                                    do nct = 1, kd !ug/m3/s -> kg/m2/s
                                    col_tmp = col_tmp + matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))*1E-9
                                    enddo
                            elseif(matrix_all_tracer(n)%type .eq. 'number') then
                                    do nct = 1, kd !#/m3/s -> mol/m2/s
                                    col_tmp = col_tmp + &
                                            matrix_all_tracer(n)%value_in_matrix(is:ie,js:je,nct)*(zhalf(:,:,nct)-zhalf(:,:,nct+1))/AVOGNO
                                    enddo
                            endif
                            used = send_data(matrix_all_tracer(n)%id_tracer_col, col_tmp, time, &
                                    is_in=is,js_in=js)
                    endif
                    enddo



            endif



    end subroutine matrix_run



    !----------------------------------------------------------------
    !
    !                           subroutine matrix_init(r)
    !
    !-------------------Initialization-------------------
    !   1. determine configuration and number of populations: npop, allocate(matrix_pop(npop))
    !   2. determine number of tracers and get parameters
    !----------------------------------------------------------------
    subroutine matrix_init(r, axes, time)
        real, intent(in), dimension(:,:,:,:) :: r    
        integer, intent(in) :: axes(4)              ! diagnostic axes
        type(time_type),  intent(in) :: time       ! model time
        !----------------------------------------------------------------
        !	------------------- Read nml -------------------
        !----------------------------------------------------------------
        if (matrix_module_init) return

        if(file_exist('input.nml')) then
            #ifdef INTERNAL_FILE_NML
            read (input_nml_file, nml = matrix_nml, iostat = io)
            ierr = check_nml_error(io,'matrix_nml')
            #else
            unit = open_namelist_file('input.nml')
            ierr=1; do while (ierr /= 0)
            read(unit, nml = matrix_nml, iostat = io, end = 10)
            ierr = check_nml_error (io, 'matrix_nml')
            end do
            10  call close_file(unit)
            #endif
        endif

        logunit = stdlog()
        if(mpp_pe() == mpp_root_pe()) then
            write(logunit, nml=matrix_nml)
            verbose = verbose + 1
        end if

        if (.not. do_matrix) then
            return
        end if

        !save ids for important tracers
        nsphum = get_tracer_index(MODEL_ATMOS,'sphum') 
        nh2so4 = get_tracer_index(MODEL_ATMOS,'simpleH2SO4')
        if (nh2so4.le.0) then
            nh2so4 = get_tracer_index(MODEL_ATMOS,'H2SO4')
        end if

        if (nh2so4.le.0) &
            call error_mesg ('matrix_gfdl','H2SO4 needs to be defined!!!!', FATAL)
        call set_config_pop(matrix_configuration)
        !----------------------------------------------------------------
        !	--- assign configuration population information ----
        ! akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
        ! if population exist, assign index >= 1, otherwise -1
        !----------------------------------------------------------------
       ! if (matrix_configuration .eq. "debug") then
       !     !debug version: 
       !     !pop: ACC <- (N,M), EXT <- (M_H2O, M_H2SO4)
       !     !npop = 2 !number of poplation in the configuration
       !     !assign population index under current configuration
       !     !integer :: I_AKK = 1, I_ACC = 2, I_DD1 = 3, I_DD2 = 4 ! index of population in matrix
       !     !integer :: I_SSA = 5, I_SSC = 6, I_OC1 = 7, I_OC2 = 8 ! if exit, index >= 1, otherwise = -1
       !     !integer :: I_BC1 = 9, I_BC2 = 10, I_MXA = 11, I_MXC = 12, I_EXT = 13
       !     I_ACC = 2
       !     I_EXT = 13
       !     pop_active = 1 !neglect ext
       !     allocate(I_POP(pop_active))
       !     I_POP(:)=(/I_ACC/)
       ! elseif (matrix_configuration .eq. "debug_npf") then
       !     I_AKK = 1
       !     pop_active = 1
       !     allocate(I_POP(pop_active))
       !     I_POP(:)=(/I_AKK/)
       ! else
       !     call ERROR_MESG('get matrix configuration','configuration not defined '//trim(matrix_configuration), FATAL)
       !     !call error_mesg ('Tracer_driver', 'mw needs to be defined for tracer: '//trim(tracer_name), FATAL)
       ! end if

        call get_number_tracers(MODEL_ATMOS, num_tracers = ntrace)

        allocate(matrix_all_tracer(ntrace))
        allocate(matrix_all_pop(npop))
        do n=1,npop
                matrix_all_pop(n)%nb_tracer_pop = 0
                !                matrix_all_pop(n)%tracer_index(:) = -1
                matrix_all_pop(n)%tracer_index = -1
        enddo
        !XL DEBUG-----------------------------------------
        !if (mpp_root_pe().eq.mpp_pe()) then
        !     write(*,*) "debug1"
        !     write(*,*) "ntrace=",ntrace
        !     write(*,*) "npop=", npop
        !  endif

          !F1P: remove when you are sure you understand what's going on
        !  do n=1,npop
        !     write(*,*) mpp_pe(),"npop=",n,"matrix_all_pop(n)%nb_tracer_pop",matrix_all_pop(n)%nb_tracer_pop,"matrix_all_pop(n)%tracer_index",matrix_all_pop(n)%tracer_index              
        !  end do
          
        !XL DEBUG1-------------------------------------------
        if (ntrace > 0) then
           do n = 1, ntrace
              flag = query_method ('matrix_parameter',MODEL_ATMOS,n, &
                   text_in_scheme,control)
              !XL DEBUG-----------------------------------------
              !if (mpp_root_pe().eq.mpp_pe()) then
              !   write(*,*) "debug2"
              !   write(*,*) "ntrace index=", n
              !   write(*,*) "text_in_scheme =", text_in_scheme
              !   write(*,*) "control = ", control
              !   write(*,*) "flag = ", flag
              !endif
              !XL DEBUG1-------------------------------------------
              if (flag) then
                 !IMPORTANT: get_matrix_tracer_param only 
                 !(1) update has_emission for mass species
                 !(2) allocate source array for mass species with lognormal, dens etc parameters in field table
                 !(3) IMPORTANT: number tracer: has_emission & source_array allocation need to be set somewhere

               !  if (mpp_root_pe().eq.mpp_pe()) then
               !     write(*,*) "debug3"
               !     write(*,*) "n=", n
               !     write(*,*) "text_in_scheme =", text_in_scheme
               !     write(*,*) "control = ", control
               !     write(*,*) "flag = ", flag
               !  endif
                 call get_matrix_tracer_param(n,text_in_scheme, control, matrix_all_tracer(n), matrix_all_pop, r)
                 if (matrix_all_tracer(n)%is_active) then !active matrix tracer in current configuration
                    allocate(matrix_all_tracer(n)%source(size(r,1),size(r,2),size(r,3)))!mass tracers has sources !!!!XL!!!!pop_index the number need
                    matrix_all_tracer(n)%source(:,:,:) = 0.
                    allocate(matrix_all_tracer(n)%value_in_matrix(size(r,1),size(r,2),size(r,3)))
                    matrix_all_tracer(n)%value_in_matrix = 0.
                 endif
        end if
        end do

    else
        call ERROR_MESG('matrix_gfdl/matrix_init', 'no atmos_tracer found ', FATAL)
    endif

    !set number tracer: has_emission & allocate source_array allocation
    do n=1,npop
    if (matrix_all_pop(n)%nb_tracer_pop > 0) then
        do nt = 1, matrix_all_pop(n)%nb_tracer_pop
        ntt = matrix_all_pop(n)%tracer_index(nt)
        if (matrix_all_tracer(ntt)%has_emission) then
            matrix_all_tracer(matrix_all_pop(n)%I_N)%has_emission = .TRUE.
        endif 
        enddo 
    endif    
    enddo


    do n=1,npop
        allocate(matrix_all_pop(n)%Dg_dry(size(r,1),size(r,2),size(r,3)))
        allocate(matrix_all_pop(n)%Dg_wet(size(r,1),size(r,2),size(r,3)))
        allocate(matrix_all_pop(n)%mass_dry(size(r,1),size(r,2),size(r,3)))
        allocate(matrix_all_pop(n)%vol_dry(size(r,1),size(r,2),size(r,3)))
        allocate(matrix_all_pop(n)%kappa_pop(size(r,1),size(r,2),size(r,3)))
        allocate(matrix_all_pop(n)%dens_wet(size(r,1),size(r,2),size(r,3)))
        allocate(matrix_all_pop(n)%dens_dry(size(r,1),size(r,2),size(r,3)))
        matrix_all_pop(n)%Dg_dry = 0.
        matrix_all_pop(n)%Dg_wet = 0.
        matrix_all_pop(n)%mass_dry = 0.
        matrix_all_pop(n)%vol_dry = 0.
        matrix_all_pop(n)%kappa_pop = 0.
        matrix_all_pop(n)%dens_wet = 0.
        matrix_all_pop(n)%dens_dry = 0.
    end do

    allocate(P_H2SO4_RATE(size(r,1),size(r,2),size(r,3))) !H2SO4 production
    allocate(AQSO4(size(r,1),size(r,2),size(r,3))) !aqueous phase SO4 production
    
    !assign rh_deliquescence and rh_crystalization to population
    do n = 1, npop
        if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
                select case (lowercase(trim(matrix_all_pop(n)%name)))
                case ('akk')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('acc')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('dd1')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('dd2')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('ssa')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('ssc')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('oc1')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('oc2')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('bc1')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('bc2')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('mxx')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case ('ext')
                        matrix_all_pop(n)%rh_deliquescence = 0.8
                        matrix_all_pop(n)%rh_crystallization = 0.2
                case default
                        matrix_all_pop(n)%rh_deliquescence = 0
                        matrix_all_pop(n)%rh_crystallization = 0
                end select
        endif
    end do 
    !set up coagulation configuration and arrays
    if (do_coag) then 
        call setup_coag_tensors(coag_configuration)
    endif
   
   !------------------------------------------=
   ! CITABLE, NM, PROD_INDEX tested 
   ! !XL DEBUG
   ! if (mpp_root_pe().eq.mpp_pe()) then
   !     write(*,*) "MODE_NAME", MODE_NAME
   !     write(*,*) "CITABLE"
   !     write(*,*) size(CITABLE, 1)
   !     do n=1,12
   !             write(*,*) n,  CITABLE(n,:)
   !     enddo
   !     write(*,*) "NM"
   !     write(*, *) NM
   !     write(*,*) "PROD_INDEX"
   !     do n=1,12
   !             write(*,*) n,  PROD_INDEX(n, 1:NMASS_SPCS)
   !     enddo
   ! endif


    kd = size(r,3)
    !sanility check
    do n=1,npop
    if (matrix_all_pop(n)%nb_tracer_pop > 0) then !if this population exist in current configuration
            do nt = 1,matrix_all_pop(n)%nb_tracer_pop
            if (matrix_all_pop(n)%has_emission(nt)) then
                    ntt = matrix_all_pop(n)%tracer_index(nt)
                    if ((matrix_all_tracer(ntt)%sigma > 0) .AND. (matrix_all_tracer(ntt)%distribution_index == DIST_LOGNORMAL) &
                            .AND. (matrix_all_tracer(ntt)%Dgn > 0) .AND. (matrix_all_tracer(ntt)%dens > 0) ) then !SUGGEST TO MOVE THIS CHECK TO INITIALIZATION
                    else
                            call error_mesg ('matrix_gfdl','Emission species not properly defined', FATAL)
                    endif
            endif
            enddo
    endif
    enddo


    do n = 1, npop
      if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
        matrix_all_pop(n)%id_Dg_dry = register_diag_field ( module_name,     &
                trim(matrix_all_pop(n)%name)//'_dg_dry', axes(1:3),Time,  &
                trim(matrix_all_pop(n)%name)//'_dg_dry', 'µm',       &
                missing_value=-999.  )
                
        matrix_all_pop(n)%id_Dg_wet = register_diag_field ( module_name,     &
                trim(matrix_all_pop(n)%name)//'_dg_wet', axes(1:3),Time,  &
                trim(matrix_all_pop(n)%name)//'_dg_wet', 'µm',       &
                missing_value=-999.  )

        !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
                 !write(*,*) "id_register test"
                 !write(*,*) "pop name = ", trim(matrix_all_pop(n)%name)
                 !write(*,*) "diag 1st wet string = ", trim(matrix_all_pop(n)%name)//'_dg_wet'
                 !write(*,*) "diag 2nd wet string = ", trim(matrix_all_pop(n)%name)//'_dg_wet'
                 !write(*,*) "dry id =", matrix_all_pop(n)%id_Dg_dry
                 !write(*,*) "wet id =", matrix_all_pop(n)%id_Dg_wet
        !endif
      endif
    end do

    id_RH = register_diag_field ( module_name,'rh_matrix', axes(1:3),Time, 'rh_matrix', '%',       &
                        missing_value=-999.  ) 
    id_cond_sink = register_diag_field ( module_name,'cond_sink', axes(1:3),Time, 'cond_sink', 'ug/m3/(30 min)',       &
                        missing_value=-999.  )
    id_kc = register_diag_field ( module_name,'kc', axes(1:3),Time, 'kc', '1/s',       &
                        missing_value=-999.  )
    id_dmdt_h2so4_tot_cond_npf = register_diag_field ( module_name,'dmdt_h2so4_tot_cond_npf', axes(1:3),Time, 'dmdt_h2so4_tot_cond_npf', 'ug/m3/s',    &
                        missing_value=-999.  )
    id_dndt_npf = register_diag_field ( module_name,'dndt_npf', axes(1:3),Time, 'dndt_npf', '#/m3/s',    &
                        missing_value=-999.  )
    id_dmdt_h2so4_npf = register_diag_field ( module_name,'dmdt_h2so4_npf', axes(1:3),Time, 'dmdt_h2so4_npf', 'ug/m3/s',    &
                        missing_value=-999.  )
    id_pwt = register_diag_field ( module_name,'pwt', axes(1:3),Time, 'pwt', 'kg air/m2',    &
                        missing_value=-999.  )
    id_zhalf = register_diag_field ( module_name,'zhalf', axes(1:3),Time, 'zhalf', 'm',    &
                        missing_value=-999.  )

    !XL matrix tracer register
    id_so4_emis = register_diag_field ( module_name,'so4_emis', axes(1:2),Time, 'so4_emis', 'kg/m2/s',    &
                        missing_value=-999.  )
    do n=1, ntrace
        if (matrix_all_tracer(n)%is_active) then
                if (matrix_all_tracer(n)%type .eq. 'mass') then
                        matrix_all_tracer(n)%id_tracer_emis = register_diag_field ( module_name,     &
                                trim(matrix_all_tracer(n)%name)//'_emis', axes(1:2),Time,  &
                                trim(matrix_all_tracer(n)%name)//'_emis', 'kg/m2/s',       &
                                missing_value=-999.  )
                         matrix_all_tracer(n)%id_tracer_col = register_diag_field ( module_name,     &
                                trim(matrix_all_tracer(n)%name)//'_col', axes(1:2),Time,  &
                                trim(matrix_all_tracer(n)%name)//'_col', 'kg/m2',       &
                                missing_value=-999.  )
                elseif (matrix_all_tracer(n)%type .eq. 'number') then
                         matrix_all_tracer(n)%id_tracer_emis = register_diag_field ( module_name,     &
                                trim(matrix_all_tracer(n)%name)//'_emis', axes(1:2),Time,  &
                                trim(matrix_all_tracer(n)%name)//'_emis', 'mol/m2/s',       &
                                missing_value=-999.  )
                         matrix_all_tracer(n)%id_tracer_col = register_diag_field ( module_name,     &
                                trim(matrix_all_tracer(n)%name)//'_col', axes(1:2),Time,  &
                                trim(matrix_all_tracer(n)%name)//'_col', 'mol/m2',       &
                                missing_value=-999.  )
               endif

        endif
        if (n .eq. nh2so4) then
                id_h2so4_emis = register_diag_field ( module_name,     &
                        'h2so4_emis', axes(1:2),Time,  &
                        'h2so4_emis', 'kg/m2/s', missing_value=-999.  )
                id_h2so4_col = register_diag_field ( module_name,     &
                        'h2so4_col', axes(1:2),Time,  &
                        'h2so4_col', 'kg/m2', missing_value=-999.  )
        endif
    enddo

    do n=1,npop
    if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
            matrix_all_pop(n)%id_pop_condens = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_condens', axes(1:2),Time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_condens', 'kg/m2/s',       &
                    missing_value=-999.  )
            matrix_all_pop(n)%id_pop_coag_pN = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pN', axes(1:3),Time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pN', '#/m3/s',       &
                    missing_value=-999.  )
            matrix_all_pop(n)%id_pop_coag_lN = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lN', axes(1:3),Time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lN', '#/m3/s',       &
                    missing_value=-999.  )                
            matrix_all_pop(n)%id_pop_coag_pMSULF = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pMSULF', axes(1:3),Time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_pMSULF', 'ug/m3/s',       &
                    missing_value=-999.  )  
            matrix_all_pop(n)%id_pop_coag_lMSULF = register_diag_field ( module_name,     &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lMSULF', axes(1:3),Time,  &
                    lowercase(trim(matrix_all_pop(n)%name))//'_coag_lMSULF', 'ug/m2/s',       &
                    missing_value=-999.  )

    endif
    enddo
    !initialize clock id
    ini_clock = mpp_clock_id ('matrix: initialization', grain=CLOCK_MODULE)
    npf_clock = mpp_clock_id ('matrix: NPF', grain=CLOCK_MODULE)
    condgrow_clock = mpp_clock_id ('matrix: condensation_grow', grain=CLOCK_MODULE)
    condgrow_clock = mpp_clock_id ('matrix: hygroscopic_growth', grain=CLOCK_MODULE)
    coag_clock = mpp_clock_id ('matrix: coagulation', grain=CLOCK_MODULE)
    dmodal_clock = mpp_clock_id ('matrix: inter-modal transfer', grain=CLOCK_MODULE)
    aqso4_clock = mpp_clock_id ('matrix: aqSO4 partition', grain=CLOCK_MODULE)
    matrix_module_init = .true.
end subroutine matrix_init


!----------------------------------------------------------------
!				subroutine get_matrix_tracer_param
!	...assign tracer properties in tracer type; 
! 	...assign index into population 
!----------------------------------------------------------------
subroutine get_matrix_tracer_param(tracer_index, population,control,MT,MP,r)
        integer, intent(in) :: tracer_index !!tracer_index is the index in all atmos tracers from get_number_tracers
        character(len=512), intent(in)  :: population
        character(len=512), intent(in)  :: control
        type(matrix_tracer), intent(inout) :: MT !note: this is one tracer
        type(matrix_pop),    intent(inout) :: MP(npop) !note: this is all populations in the configuration
        real, intent(in) :: r(:,:,:,:)
        integer :: pop_index
        integer :: iflag
        pop_index = 0
        ! trim(population) name: akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
        if (lowercase(trim(population))=='akk') then
                MT%pop = 'akk'
                pop_index = I_AKK
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='acc') then
                MT%pop = 'acc'
                pop_index = I_ACC
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='dd1') then
                MT%pop = 'dd1'
                pop_index = I_DD1
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='dd2') then
                MT%pop = 'dd2'
                pop_index = I_DD2
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='ssa') then
                MT%pop = 'ssa'
                pop_index = I_SSA
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='ssc') then
                MT%pop = 'ssc'
                pop_index = I_SSC
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='oc1') then
                MT%pop = 'oc1'
                pop_index = I_OC1
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='oc2') then
                MT%pop = 'oc2'
                pop_index = I_OC2
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='bc1') then
                MT%pop = 'bc1'
                pop_index = I_BC1
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='bc2') then
                MT%pop = 'bc2'
                pop_index = I_BC2
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='mxa') then
                MT%pop = 'mxa'
                pop_index = I_MXA
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='mxc') then
                MT%pop = 'mxc'
                pop_index = I_MXC
                MT%pop_index = pop_index
        elseif (lowercase(trim(population))=='ext') then
                MT%pop = 'ext'
                pop_index = I_EXT
                MT%pop_index = pop_index
        else
                call ERROR_MESG('get_matrix_tracer_param', 'trim(population) not found '//trim(trim(population)), FATAL )        
        endif

        if (mpp_root_pe().eq.mpp_pe()) then
                write(*,*) "debug4"
                write(*,*) "ntrace index=", tracer_index
                write(*,*) "lowercase(trim(population))", lowercase(trim(population))
                write(*,*) "pop_index", pop_index
                write(*,*) "MT%pop_index ", MT%pop_index
        endif

        if (pop_index > 0) then !if the current configuration have this population and this tracer
                !    MT%index_in_atmos = tracer_index
                iflag=parse(control,'type',MT%type)
                if (iflag>0) then
                        MT%is_active=.true.
                end if
                iflag=parse(control,'spec',MT%spec)
                iflag=parse(control,'sigma',MT%sigma)
                if (iflag>0) then
                        MT%lnsigma  = log(MT%sigma)
                        MT%lnsigma2 = MT%lnsigma**2
                end if
                iflag=parse(control,'Dgn',MT%Dgn)
                if (iflag>0) then
                        MT%DP0 = MT%DGN * exp(1.5*MT%LNSIGMA2)
                end if
                iflag=parse(control,'distribution',MT%distribution) 
                ! distribution_index !=1, lognormal; =2, weibull
                if (trim(MT%distribution).eq."lognormal") then
                        MT%distribution_index = DIST_LOGNORMAL
                elseif (trim(MT%distribution).eq."weibull") then
                        MT%distribution_index = DIST_WEIBULL       
                end if
                !only mass species with emission have distribution parameters in the field table 
                if (iflag>0) then
                        MT%has_emission = .TRUE. !this line only updated mass tracers with emission
                else
                        MT%has_emission = .FALSE.
                end if
                iflag=parse(control,'kappa',MT%kappa)
                iflag=parse(control,'dens',MT%dens)
                call get_tracer_names (MODEL_ATMOS, tracer_index, name = MT%name,  &
                        units = MT%units)
!        endif


!        if (pop_index > 0) then !if the current configuration have this population
                MP(pop_index)%nb_tracer_pop = MP(pop_index)%nb_tracer_pop + 1   ! count the number of tracers in the specific populations 
                if (mpp_root_pe().eq.mpp_pe()) then
                        write(*,*) "debug5"
                        write(*,*) "pop_index=", pop_index
                        write(*,*) "MT%name", MT%name
                        write(*,*) "MT%pop", MT%pop
                        write(*,*) "MP(pop_index)%nb_tracer_pop", MP(pop_index)%nb_tracer_pop
                        write(*,*) "tracer_index", tracer_index
                endif

                if (MP(pop_index)%nb_tracer_pop.gt.nb_tracer_max) then
                       write(*,*) "debug_crash"
                       write(*,*) "MP(pop_index)%tracer_index",MP(pop_index)%tracer_index
                        write(*,*) "pop_index=", pop_index
                        write(*,*) "MT%name", MT%name
                        write(*,*) "MT%pop", MT%pop
                        write(*,*) "MP(pop_index)%nb_tracer_pop", MP(pop_index)%nb_tracer_pop
                        write(*,*) "tracer_index", tracer_index
                end if
                     
                MP(pop_index)%tracer_index((MP(pop_index)%nb_tracer_pop)) = tracer_index
                if (mpp_root_pe().eq.mpp_pe()) then
                        write(*,*) "debug5_af"
                        write(*,*) "MP(pop_index)%nb_tracer_pop=", MP(pop_index)%nb_tracer_pop
                        write(*,*) "MP(pop_index)%tracer_index=", MP(pop_index)%tracer_index
                        write(*,*) "MP(pop_index)%tracer_index(MP(pop_index))%nb_tracer_pop)", MP(pop_index)%tracer_index(MP(pop_index)%nb_tracer_pop)
                        write(*,*) "-------------end----------------------"
                endif
                MP(pop_index)%has_emission(MP(pop_index)%nb_tracer_pop) = MT%has_emission
                if (trim(MT%type).eq."number") then
                        MP(pop_index)%I_N = tracer_index !record the index in the matrix_tracer_array
                        MP(pop_index)%name = lowercase(trim(population))
                elseif (trim(MT%type).eq."mass") then
                        if  (trim(MT%spec).eq. "sulf") then
                                MP(pop_index)%I_MSULF = tracer_index           
                        elseif (trim(MT%spec).eq. "dust") then
                                MP(pop_index)%I_MDUST = tracer_index
                        elseif (trim(MT%spec).eq. "seas") then
                                MP(pop_index)%I_MSEAS = tracer_index
                        elseif (trim(MT%spec).eq. "ocar") then
                                MP(pop_index)%I_MOCAR = tracer_index
                        elseif (trim(MT%spec).eq. "bcar") then
                                MP(pop_index)%I_MBCAR = tracer_index
                        elseif (trim(MT%spec).eq. "alwc") then
                                MP(pop_index)%I_MWATE = tracer_index 
                                MP(pop_index)%name = lowercase(trim(population)) 
                        elseif (trim(MT%spec).eq. "ammo") then
                                MP(pop_index)%I_MAMMO = tracer_index
                        elseif (trim(MT%spec).eq. "nitr") then
                                MP(pop_index)%I_MNITR = tracer_index    
                        elseif (trim(MT%spec).eq. "gsfa") then !gaseous sulfuric acid
                                MP(pop_index)%I_MGSFA = tracer_index
                        else
                                call error_mesg ('matrix_gfdl','matrix tracer: field table property not properly defined', FATAL)
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
! given a 3D array and its units in gfdl: vmr, mmr, #/kg
! convert to matrix units: ug/m3 for mass, #/m3 for number
!----------------------------------------------------------------
subroutine set_unit_gfdl_to_matrix(monotracer, tr_unit, tr_type, tr_spec, value_in_matrix, pwt, zhalf) ! note: avoid using unit/type, because &
    !they are used previously for other purposes
    real, intent(in) :: monotracer(:,:,:)!single tracer concentration in gfdl unit
    character*32, intent(in) :: tr_unit !tracer unit defined in field table, vmr/mmr/#/kg
    character*32, intent(in) :: tr_type, tr_spec !matrix tracer type: number/mass; tracer spec: sulf/dust/seas/ocar/bcar/alwc/gsfa
    real, intent(in), optional :: pwt(:,:,:), zhalf(:,:,:)
    real, intent(inout) :: value_in_matrix(:,:,:) !value in matrix unit
    real :: MW !local variable
    if (lowercase(trim(tr_type)) .eq. "number") then !number type: convert into #/m3
        if (lowercase(trim(tr_unit)) .eq. "vmr") then
            value_in_matrix = monotracer * pwt(:,:,:)  / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
        elseif (lowercase(trim(tr_unit)) .eq. "#/kg") then
            value_in_matrix = monotracer * pwt(:,:,:) / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
        else
            call error_mesg ('matrix_gfdl','matrix tracer: number type unit not properly defined', FATAL)
        endif
    elseif (lowercase(trim(tr_type)) .eq. "mass") then !mass type: convert into ug/m3
        if (lowercase(trim(tr_unit)) .eq. "vmr") then !only sulfate can possible to be vmr unit
            if (lowercase(trim(tr_spec)) .eq. "sulf") then
                MW = 96
                value_in_matrix = 1E9 * monotracer  * pwt(:,:,:) * MW / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
            else
                call error_mesg ('matrix_gfdl','matrix tracer: molecular weight not properly defined', FATAL)
            endif
        elseif (lowercase(trim(tr_unit)) .eq. "mmr") then
            value_in_matrix  = 1E9 * monotracer * pwt(:,:,:)  / (zhalf(:,:,1:kd) - zhalf(:,:,2:kd+1))
        else
            call error_mesg ('matrix_gfdl','matrix tracer: number type unit not properly defined', FATAL)
        endif
    endif          
end subroutine set_unit_gfdl_to_matrix


!----------------------------------------------------------------
!                subroutine update_rt_from_matrix(n_rt, mono_rt_tracer)
! given the index of tracers, update gfdl tracer values from matrix%value
!----------------------------------------------------------------
subroutine update_rt_from_matrix(n_rt, mono_rt_tracer,pwt, zhalf,is,ie,js,je)
    integer, intent(in) :: n_rt !index in rt(:,:,:,n_rt)
    real, intent(in), optional :: pwt(:,:,:), zhalf(:,:,:)
    integer, intent(in) :: is,ie,js,je
    real, intent(out) :: mono_rt_tracer(:,:,:) !rt(:,:,:,n_rt) to be updated
    real :: MW !molecular weight, local variable
    !number or mass tracers
    if (lowercase(matrix_all_tracer(n_rt)%type) .eq. "mass") then
        if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "mmr") then
            mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:)/1E9/(pwt(:,:,:)/ (zhalf(:,:,1:kd) - zhalf(:,:,2:kd+1)))
        elseif (lowercase(matrix_all_tracer(n_rt)%units) .eq. "vmr") then
            if (lowercase(trim(matrix_all_tracer(n_rt)%spec)) .eq. "sulf") then
                MW = 96
                mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:)/1E9/(pwt(:,:,:) * MW / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
            else
                call error_mesg ('matrix_gfdl','matrix tracer: molecular weight not properly defined', FATAL)
            endif
        endif
    elseif (lowercase(matrix_all_tracer(n_rt)%type) .eq. "number") then
        if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "#/kg") then
            mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:) / (pwt(:,:,:) / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
        elseif (lowercase(matrix_all_tracer(n_rt)%units) .eq. "vmr") then
            mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix(is:ie,js:je,:) / (pwt(:,:,:)  / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
        else
            call error_mesg ('matrix_gfdl','matrix tracer: number unit not properly defined', FATAL)
        endif
    endif
end subroutine update_rt_from_matrix


!----------------------------------------------------------------
!                subroutine set_matrix_source_2d
! 1. broadcast 2D source -> 3D source
! 2. convert unit from vmr/s, mmr/s, kg/m2/s to µg/m3/s
! conversion: 
! 1. vmr/s -> kg/m2/s: vmr/s * pwt * MW / WTMAIR
! 2. mmr/s -> kg/m2/s: mmr/s * pwt
! 3. kg/m2/s -> µg/m3/s: kg/m2/s / h_zgrid * 1E9
!----------------------------------------------------------------
subroutine set_matrix_source_2d(source_type,source,source_units,pwt,zhalf, time, is, js)

    integer, intent(in) :: source_type, source_units
    real ,   intent(in) :: source(:,:)
    real,    intent(in), optional :: pwt(:,:,:), zhalf(:,:,:)
    type(time_type), intent(in) :: time ! current model time
    integer, intent(in) :: is, js ! boundaries of physical window
    real :: source_processed(size(source,1),size(source,2),kd) !kd = size(r,3)
    real :: MW ! molecular weight of the specific source type
    !determine molecular weight MW
    if (source_type .eq. matrix_source_type%P_H2SO4) then
        MW = matrix_molecular_weight(I_MW_H2SO4)
    elseif (source_type .eq. matrix_source_type%E_SO4) then
        MW = matrix_molecular_weight(I_MW_SO4)
    elseif (source_type .eq. matrix_source_type%E_SS) then
        MW = matrix_molecular_weight(I_MW_SS)
    elseif (source_type .eq. matrix_source_type%E_DUST) then
        MW = matrix_molecular_weight(I_MW_DUST)
    elseif(source_type .eq. matrix_source_type%E_OC) then
        MW = matrix_molecular_weight(I_MW_OC)
    elseif(source_type .eq. matrix_source_type%E_BC) then
        MW = matrix_molecular_weight(I_MW_BC)
    elseif(source_type .eq. matrix_source_type%E_SOA) then
        MW = matrix_molecular_weight(I_MW_SOA)
    elseif(source_type .eq. matrix_source_type%P_AQSO4) then
        MW = matrix_molecular_weight(I_MW_SO4)
    endif

    source_processed(:,:,:) = 0.

    if (do_matrix) then
        if (source_units .eq.  matrix_source_type%U_KG_M2_S) then
            source_processed(:,:,kd) = 1E9 * source / (zhalf(:,:,kd) - zhalf(:,:,kd+1)) !convert µg/m3/s
        else
            call error_mesg ('matrix_gfdl','Emission source dimension is not compatible with its unit', FATAL)
        end if
        call set_matrix_source_generic(source_type,source_processed, zhalf, time, is, js)  
    end if
end subroutine set_matrix_source_2d


subroutine set_matrix_source_3d(source_type,source,source_units,pwt,zhalf, time, is,js)

    integer, intent(in) :: source_type, source_units
    real ,   intent(in) :: source(:,:,:)
    real,    intent(in), optional :: pwt(:,:,:), zhalf(:,:,:) 
    type(time_type), intent(in) :: time ! current model time
    integer, intent(in) :: is, js ! boundaries of physical window
    real :: MW
    real ::    source_processed(size(source,1),size(source,2),size(source,3))
    !determine molecular weight MW
    if (source_type .eq. matrix_source_type%P_H2SO4) then
        MW = matrix_molecular_weight(I_MW_H2SO4)
    elseif (source_type .eq. matrix_source_type%E_SO4) then
        MW = matrix_molecular_weight(I_MW_SO4)
    elseif (source_type .eq. matrix_source_type%E_SS) then
        MW = matrix_molecular_weight(I_MW_SS)
    elseif (source_type .eq. matrix_source_type%E_DUST) then
        MW = matrix_molecular_weight(I_MW_DUST)
    elseif(source_type .eq. matrix_source_type%E_OC) then
        MW = matrix_molecular_weight(I_MW_OC)
    elseif(source_type .eq. matrix_source_type%E_BC) then
        MW = matrix_molecular_weight(I_MW_BC)
    elseif(source_type .eq. matrix_source_type%E_SOA) then
        MW = matrix_molecular_weight(I_MW_SOA)
    elseif(source_type .eq. matrix_source_type%P_AQSO4) then
        MW = matrix_molecular_weight(I_MW_SO4)
    endif

    source_processed(:,:,:) = 0.

    if (do_matrix) then
        if (source_units .eq.  matrix_source_type%U_VMR_S) then
            source_processed(:,:,:) = 1E9 * source(:,:,:) * pwt(:,:,:) * MW / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1))
        elseif (source_units .eq.  matrix_source_type%U_MMR_S) then
            source_processed(:,:,:) = 1E9 * source(:,:,:) * pwt(:,:,:)  / (zhalf(:,:,1:kd) - zhalf(:,:,2:kd+1))
        else
            call error_mesg ('matrix_gfdl','Emission source dimension is not compatible with its unit', FATAL)
        end if
        call set_matrix_source_generic(source_type,source_processed, zhalf, time, is,js)
    end if

end subroutine set_matrix_source_3d

subroutine set_matrix_source_generic(source_type,source_processed, zhalf, time, is,js)

    real, intent(in)    :: source_processed(:,:,:), zhalf(:,:,:)
    integer, intent(in) :: source_type
    type(time_type), intent(in) :: time ! current model time
    integer, intent(in) :: is, js ! boundaries of physical window
    integer :: ie, je
    integer :: nk
    real :: h2so4_emis(size(source_processed,1),size(source_processed,2)), so4_emis(size(source_processed,1),size(source_processed,2))
    logical :: used 
    ie = is+size(source_processed,1)-1
    je = js+size(source_processed,2)-1
    h2so4_emis = 0.
    so4_emis = 0.
    if (do_matrix) then
        !AKK mode H2SO4(g): vmr/s, 3D array
        if (source_type .eq. matrix_source_type%P_H2SO4) then   
            P_H2SO4_RATE(is:ie,js:je,:) = source_processed
            do nk=1,kd
                h2so4_emis = h2so4_emis + source_processed(:,:,nk)*(zhalf(:,:,nk)-zhalf(:,:,nk+1))*1E-9 !convert ug/m3/s -> kg/m2/s
            enddo
            used = send_data (id_h2so4_emis, h2so4_emis, time, &
                                        is_in=is,js_in=js)
            !ACC mode SO4 emission: vmr/s, 3D array
        elseif (source_type .eq. matrix_source_type%E_SO4) then
            if (I_AKK > 0 .AND. I_ACC < 0) then
                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source(is:ie,js:je,:) =0.01* source_processed
            elseif (I_AKK > 0 .AND. I_ACC > 0) then
                    matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%source(is:ie,js:je,:) =0.99* source_processed
                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source(is:ie,js:je,:) =0.01* source_processed
            else
                     matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%source(is:ie,js:je,:) =source_processed
            endif
            do nk=1,kd
                so4_emis = so4_emis + source_processed(:,:,nk)*(zhalf(:,:,nk)-zhalf(:,:,nk+1))*1E-9 !convert ug/m3/s -> kg/m2/s
            enddo
            used = send_data (id_so4_emis, so4_emis, time, &
                                        is_in=is,js_in=js)
        elseif (source_type .eq. matrix_source_type%P_AQSO4) then
            AQSO4(is:ie,js:je,:) =  source_processed
            !Dust emission: Kg/m2/s, 2D array
        elseif (source_type .eq. matrix_source_type%E_DUST) then
            if (I_DD1>0 .AND. I_DD2<0) then !only have 1 dust mode DD1
                matrix_all_tracer(matrix_all_pop(I_DD1)%I_MDUST)%source(is:ie,js:je,:) = source_processed
            elseif (I_DD1<0 .AND. I_DD2>0) then !only have 1 dust mode DD2
                matrix_all_tracer(matrix_all_pop(I_DD2)%I_MDUST)%source(is:ie,js:je,:) = source_processed
            elseif (I_DD1>0 .AND. I_DD2>0) then !have 2 dust mode: DD1 and DD2
                matrix_all_tracer(matrix_all_pop(I_DD1)%I_MDUST)%source(is:ie,js:je,:) = 0.25*source_processed
                matrix_all_tracer(matrix_all_pop(I_DD2)%I_MDUST)%source(is:ie,js:je,:) = 0.75*source_processed
            endif
            !Sea salt emission: Kg/m2/s, 2D array
        elseif (source_type .eq. matrix_source_type%E_SS) then
            if (I_SSA>0 .AND. I_SSC<0) then !only have 1 sea salt mode SSA
                matrix_all_tracer(matrix_all_pop(I_SSA)%I_MSEAS)%source(is:ie,js:je,:) = source_processed
            elseif (I_SSA<0 .AND. I_SSC>0) then !only have 1 sea salt mode SSC
                matrix_all_tracer(matrix_all_pop(I_SSC)%I_MSEAS)%source(is:ie,js:je,:) = source_processed
            elseif (I_SSA>0 .AND. I_SSC>0) then !have 2 sea salt mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_SSA)%I_MSEAS)%source(is:ie,js:je,:) = 0.25*source_processed
                matrix_all_tracer(matrix_all_pop(I_SSC)%I_MSEAS)%source(is:ie,js:je,:) = 0.75*source_processed
            endif
            !SOA source from GFDL model, add it to OC: kg/m2/s, 2D array
        elseif (source_type .eq. matrix_source_type%E_SOA) then
            if (I_OC1>0 .AND. I_OC2<0) then
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) + source_processed
            elseif (I_OC1<0 .AND. I_OC2>0) then
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) = &
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) + source_processed
            elseif (I_OC1>0 .AND. I_OC2>0) then !have 2 organic carbon mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) + 0.2*source_processed
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) = & 
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) + 0.8*source_processed
            endif
            !Organic carbon: vmr, 3D array
        elseif (source_type .eq. matrix_source_type%E_OC) then
            if (I_OC1>0 .AND. I_OC2<0) then
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) + source_processed
            elseif (I_OC1<0 .AND. I_OC2>0) then
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) = &
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) + source_processed
            elseif (I_OC1>0 .AND. I_OC2>0) then !have 2 organic carbon mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source(is:ie,js:je,:) + 0.5*source_processed
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) = &
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source(is:ie,js:je,:) + 0.5*source_processed
            endif
            !Black carbon: vmr, 3D array
        elseif (source_type .eq. matrix_source_type%E_BC) then
            if (I_BC1>0 .AND. I_BC2<0) then
                matrix_all_tracer(matrix_all_pop(I_BC1)%I_MBCAR)%source(is:ie,js:je,:) = source_processed
            elseif (I_BC1<0 .AND. I_BC2>0) then
                matrix_all_tracer(matrix_all_pop(I_BC2)%I_MBCAR)%source(is:ie,js:je,:) = source_processed
            elseif (I_BC1>0 .AND. I_BC2>0) then !have 2 organic carbon mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_BC1)%I_MBCAR)%source(is:ie,js:je,:) = 0.8*source_processed
                matrix_all_tracer(matrix_all_pop(I_BC2)%I_MBCAR)%source(is:ie,js:je,:) = 0.2*source_processed
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
                    matrix_all_tracer(matrix_all_pop(n)%I_N)%source(is:ie,js:je,:) = &
                            matrix_all_tracer(ntt)%source(is:ie,js:je,:)/(PI6*matrix_all_tracer(ntt)%dens*matrix_all_tracer(ntt)%DP0**3)*1E-9 !emission rate converted to number rate

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
     real, intent(in) :: pfull(:,:,:) ! pressure on layers, Pa
     real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degK
     real, intent(inout) :: kci_coef_pop(:,:,:,:)
     real, intent(inout) :: kci_aeq1_pop(:,:,:,:) !KCI units: m^3/s
     integer :: n,it,jt,kt,i,j,k
     real :: fac,kci1,kci2
     real :: number_pop(size(pfull,1),size(pfull,2),size(pfull,3))
     integer, intent(in) :: is,ie,js,je
     kci_coef_pop = 0.
     kci_aeq1_pop = 0.
     it = size(pfull,1)
     jt = size(pfull,2)
     kt = size(pfull,3)
     do n=1, npop
         if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
            fac = exp(1.5*log(matrix_all_pop(n)%sigma)**2)
            number_pop = matrix_all_tracer(matrix_all_pop(n)%I_N)%value_in_matrix(is:ie,js:je,:)
            do i=1,it
                do j=1,jt
                    do k=1,kt
                        if (number_pop(i,j,k) > 1.0E-14) then
                                call setup_kci(pfull(i,j,k), t(i,j,k), matrix_all_pop(n)%Dg_wet(i+is-1,j+js-1,k)*fac, matrix_all_pop(n)%sigma, &
                                        kci_coef_pop(n,i,j,k), kci_aeq1_pop(n,i,j,k))
                        endif
                    enddo
                 enddo
             enddo

         endif
     enddo
 
 end subroutine


 subroutine set_matrix_condense_npf(I_AKK, XNH3, FLAND, pfull,rh,t, tstep, H2SO4_vmr,pwt,zhalf, kci_coef_pop, kci_aeq1_pop, &
                    dndt_npf, dmdt_h2so4_npf, IH2SO4_PATH, dmdt_h2so4_tot,KC, XH2SO4_NUCL,is,ie,js,je)
    real, intent(in) :: pfull(:,:,:) ! pressure on layers, Pa
    real, intent(in) :: rh(:,:,:) ! relative humidity
    real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
    real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degK
    real, intent(in) :: tstep !timestep, in s
    real, intent(in) :: H2SO4_vmr(:,:,:) !gaseous H2SO4 in gfdl unit, vmr
    real, intent(in) :: kci_coef_pop(:,:,:,:), kci_aeq1_pop(:,:,:,:) !m^3/s
    real, intent(in) :: XNH3(:,:,:), FLAND(:,:,:)
    integer, intent(in) :: I_AKK
    integer, intent(in) :: is,ie,js,je
    integer, intent(out) :: IH2SO4_PATH(:,:,:)
    real, intent(out):: dndt_npf(:,:,:), dmdt_h2so4_npf(:,:,:), dmdt_h2so4_tot(:,:,:) !m-3 s-1; ugSO4 m-3 s-1
    integer :: n,it, jt, kt, i,j,k ! number of layers in different direction
    real:: MW
    real,parameter :: MW_SO4=96
    real :: MW_H2SO4
    real, parameter :: TINYDENOM = 1.0D-30
    real, intent(out):: XH2SO4_NUCL(size(pfull,1),size(pfull,2),size(pfull,3))  ! H2SO4 (as SO4) conc. used in nucleation and GR calculation [ugSO4/m^3]
    REAL, PARAMETER :: XNTAU =2.0 !number of time consants in the current time step
    REAL, PARAMETER :: KCMIN = 1.0D-08 ! [1/s] minimum condensational sink - see notes of 10-18-06
    REAL, PARAMETER :: XH2SO4_NUCL_MIN_NCM3 = 1.00D+03 ! min. [H2SO4] to enter nucleation calculations [#/cm^3]
    REAL :: XH2SO4_NUCL_MIN ! convert to [ugH2SO4/m^3]
    real :: xH2SO4_init(size(pfull,1),size(pfull,2),size(pfull,3)) !xH2SO4_init: H2SO4 concentration in ug/m3
    real,intent(out) :: KC(:,:,:) !KC: total condensation sink (1/s) for all population in a grid
    real :: KC_AEQ1(size(pfull,1),size(pfull,2),size(pfull,3)) !KC_AEQ1: total condensation sink (1/s) using accomadation coefficient=1
    real :: number_pop(size(pfull,1),size(pfull,2),size(pfull,3)) !number concentration #/m3
    real :: SO4RATE(size(pfull,1),size(pfull,2),size(pfull,3)) ! average H2SO4 production rate [ugSO4/m^3/s]
    real :: XH2SO4_SS(size(pfull,1),size(pfull,2),size(pfull,3)), XH2SO4_SS_WNPF(size(pfull,1),size(pfull,2),size(pfull,3))
    real :: PQ_GROWTH(size(pfull,1),size(pfull,2),size(pfull,3)), TOT_H2SO4_LOSS(size(pfull,1),size(pfull,2),size(pfull,3))
    it = size(pfull,1)
    jt = size(pfull,2)
    kt = size(pfull,3)
    MW_H2SO4 = matrix_molecular_weight(I_MW_H2SO4)
    MW = matrix_molecular_weight(I_MW_H2SO4) !MW=98, H2SO4
    XH2SO4_NUCL_MIN = XH2SO4_NUCL_MIN_NCM3 * MW_H2SO4 * 1.0D+12 / AVOGNO  ! convert to [ugSO4/m^3]
    XH2SO4_NUCL =  XH2SO4_NUCL_MIN !TINYNUMER ! XH2SO4_NUCL_MIN              ! for the case  XH2SO4_INIT .LT. XH2SO4_NUCL_MIN
    XH2SO4_INIT = 1E9 * h2so4_vmr * pwt * MW / WTMAIR / (zhalf(:,:,1:kt) - zhalf(:,:,2:kt+1)) !convert gaseous H2SO4 to matrix unit: ug H2SO4/m3
    KC=0. !3D, total condensation sink at each grid, sum-up the population
    KC_AEQ1 =0. !3D, total condensation sink at each grid, sum-up the population
    SO4RATE =0. !average H2SO4 production rate [ugSO4/m^3/s]
    do n=1,npop
        if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
                number_pop = matrix_all_tracer(matrix_all_pop(n)%I_N)%value_in_matrix(is:ie,js:je,:) !number concentration of a population
                KC = KC + kci_coef_pop(n,:,:,:)*number_pop !m3/s * #/m3 = 1/s, 3D for each grid
                KC_AEQ1 = KC_AEQ1 + kci_aeq1_pop(n,:,:,:)*number_pop !m3/s * #/m3 = 1/s, 3D for each grid
        endif
    enddo
    if (I_AKK > 0) then
        SO4RATE = XH2SO4_INIT / TSTEP                        ! average H2SO4 production rate [ugSO4/m^3/s]
        do i=1,it
                do j=1,jt
                        do k=1,kt
                                if ((KC(i,j,k)*TSTEP) .GE. XNTAU) then ! invoke steady-state assumption
                                        IH2SO4_PATH(i,j,k) =1
                                        XH2SO4_SS(i,j,k) = MIN( SO4RATE(i,j,k)/KC(i,j,k), XH2SO4_INIT(i,j,k) )  !steady-state H2SO4 [ugH2SO4/m^3]
                                        
                                        CALL STEADY_STATE_H2SO4(pfull(i,j,k),rh(i,j,k),t(i,j,k),FLAND(i,j,k),XH2SO4_SS(i,j,k), &
                                                SO4RATE(i,j,k),XNH3(i,j,k),KC(i,j,k),TSTEP,XH2SO4_SS_WNPF(i,j,k))
                                        XH2SO4_NUCL(i,j,k) = XH2SO4_SS_WNPF(i,j,k)                       ! [H2SO4] for nucl., GR, and cond. calculation [ugSO4/m^3]
                                 else
                                         IH2SO4_PATH(i,j,k) =2
                                         XH2SO4_NUCL(i,j,k) = SO4RATE(i,j,k) / ( (2.0D+00/TSTEP) + KC(i,j,k) )   ! use [H2SO4] at mid-time step [ugSO4/m^3]
                                 endif

                                 
                        enddo
                 enddo
         enddo
         !XL: PLACE TO call matrix_npf
         call matrix_npfrate(pfull, rh, t, FLAND, XH2SO4_NUCL, SO4RATE, KC_AEQ1,DNDT_npf,DMDT_h2so4_npf)
    else 
        DNDT_npf = 0.
        DMDT_h2so4_npf = 0.
        XH2SO4_NUCL = XH2SO4_INIT
    endif

    do i=1,it
        do j=1,jt
                do k=1,kt
                        PQ_GROWTH(i,j,k) = XH2SO4_NUCL(i,j,k) * ( 1.0D+00 - EXP(-KC(i,j,k)*TSTEP) ) / TSTEP            ! [ugH2SO4/m^3/s]
                        TOT_H2SO4_LOSS(i,j,k) = ( dmdt_h2so4_npf(i,j,k) + PQ_GROWTH(i,j,k) ) * TSTEP                         ! [ugH2SO4/m^3]
                        IF ( TOT_H2SO4_LOSS(i,j,k) .GT. XH2SO4_INIT(i,j,k) ) THEN
                                DMDT_h2so4_npf(i,j,k)  = DMDT_h2so4_npf(i,j,k)  * ( XH2SO4_INIT(i,j,k) / ( TOT_H2SO4_LOSS(i,j,k)+TINYDENOM ) )! [ugh2SO4/m^3/s]
                                DNDT_npf(i,j,k)  = DNDT_npf(i,j,k)  * ( XH2SO4_INIT(i,j,k) / ( TOT_H2SO4_LOSS(i,j,k) + TINYDENOM ))! [  #  /m^3/s]
                                PQ_GROWTH(i,j,k) = PQ_GROWTH(i,j,k) * ( XH2SO4_INIT(i,j,k) / ( TOT_H2SO4_LOSS(i,j,k) + TINYDENOM ))! [ugH2SO4/m^3/s]
                        endif
                        dmdt_h2so4_tot(i,j,k) = DMDT_h2so4_npf(i,j,k)+PQ_GROWTH(i,j,k)
                enddo
        enddo
    enddo

 end subroutine

!-----------------------------------------------------------------------
!                 matrix_npfrate
! 1. set-up for loop to call NPFRATE in matrix
! 2. calculate DNDT and DMDT in 3D array form
!-----------------------------------------------------------------------
!matrix_npf(pfull, rh, t, SO4RATE, XH2SO4_NUCL,KC_AEQ1,DNDT,DMDT_SO4)
!NPFRATE(PRS,RH,TEMP,XH2SO4,SO4RATE,KC,DNDT,DMDT_SO4,ICALL)!XL
!matrix_npfrate(pfull, rh, t, FLAND, XH2SO4_NUCL, SO4RATE, KC_AEQ1,DNDT_npf,DMDT_npf)
subroutine matrix_npfrate(pfull,rh,t, FLAND, H2SO4, SO4_RATE,KC_AEQ1,dndt, dmdt)
    real, intent(in) :: pfull(:,:,:) ! pressure on layers, Pa
    real, intent(in) :: rh(:,:,:) ! relative humidity
    !real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
    real, intent(in) :: SO4_RATE(:,:,:) !H2SO4 production rate in ugSO4/m^3
    real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degK
    real, intent(in) :: H2SO4(:,:,:) !gaseous H2SO4 in ugSO4/m^3
    real, intent(in) :: FLAND(:,:,:)
    real, intent(in) :: KC_AEQ1(:,:,:) !condensation sink with accomadation coefficient =1
    real, intent(out):: dndt(:,:,:), dmdt(:,:,:)
    integer :: it, jt, kt, i,j,k ! number of layers in different direction
    real:: pres_mt, rh_mt, temp_mt, h2so4_mt, p_h2so4rate_mt !local variable for matrix calculation
    real:: MW, dndt_tst, dmdt_tst !if directly assign dndt(i,j,k) and dmdt(i,j,k) segment error will occur
    real:: dndt_rec(size(pfull,1),size(pfull,2),size(pfull,3)), dmdt_rec(size(pfull,1),size(pfull,2),size(pfull,3)) 
    MW = matrix_molecular_weight(I_MW_H2SO4)
    it = size(pfull,1)
    jt = size(pfull,2)
    kt = size(pfull,3)
    do i=1,it
        do j=1,jt
                do k=1,kt
                        pres_mt = pfull(i,j,k) !pressure in matrix unit: [Pa]
                        rh_mt = rh(i,j,k) !fractional relative humidity [1]
                         !XL: check with FP for RH calculation
                         temp_mt = t(i,j,k) !ambient temperature [K]
                         !sulfuric acid (as SO4) concentration [ugSO4/m^3]
                         h2so4_mt = h2so4(i,j,k)
                         p_h2so4rate_mt = SO4_RATE(i,j,k) !gas-phase H2SO4 (as SO4) production rate [ugSO4/m^3 s]
                         !dndt: [m^-3 s^-1], dmdt: [ugSO4 m^-3 s^-1]
                         !NPFRATE(PRS,RH,TEMP,FLAND,XH2SO4,SO4RATE,XNH3,KC,DNDT,DMDT_SO4,ICALL)
                         !call NPFRATE(pres_mt,rh_mt,temp_mt,h2so4_mt,p_h2so4rate_mt,dndt(i,j,k),dmdt(i,j,k))
                         call NPFRATE(pres_mt,rh_mt,temp_mt,0.,h2so4_mt,p_h2so4rate_mt,0.,kc_aeq1(i,j,k),dndt_rec(i,j,k),dmdt_rec(i,j,k),0)
                enddo
        enddo
   enddo
    dndt=dndt_rec
    dmdt=dmdt_rec
end subroutine matrix_npfrate
      !----------------------------------------------------------------------------------------------------------------
      ! Get the B_i loss       terms due to intermodal coagulation. [1/s]
      ! Get the R_i production terms due to intermodal coagulation. [#/m^3/s]
      ! Get the C_i terms, which include all source terms.          [#/m^3/s]
      ! For the C_i terms, the secondary particle formation term DNDT must be
      !   added in after coupling to condensation below.
      ! The A_i terms for intramodal coagulation are directly computed
      !   from the coagulation coefficients when the number equations
      !   are integrated.
      ! If DIKL(I,K,L) = 0, then modes K and L to not coagulate to form mode I.
      ! DIJ(I,J) is unity if coagulation of mode I with mode J results
      !   in the removal of particles from mode I, and zero otherwise.
      !----------------------------------------------------------------------------------------------------------------
subroutine matrix_coag(r,pfull, t, RI, LI, RIM, LIM, is,ie,js,je)
      implicit none
      integer, intent(in) :: is,ie,js,je
      real, intent(in) :: r(:,:,:,:), pfull(:,:,:), t(:,:,:)
      integer, parameter:: NMSPCS = 5  !5 mass spcs, check coga_config_table for details
      real :: BI(npop-1,size(r,1),size(r,2),size(r,3)) !loss terms due to intermodal coagulation. [1/s]
      real, intent(out) :: LI(npop-1,size(r,1),size(r,2),size(r,3)) !loss term of number due to intermodal coagulation [#/m^3/s]
      real, intent(out) :: RI(npop-1,size(r,1),size(r,2),size(r,3)) !production terms due to intermodal coagulation. [#/m^3/s]
      real :: FI(npop-1,size(r,1),size(r,2),size(r,3)) !mass loss coefficient
      real, intent(out) :: LIM(npop-1,NMSPCS,size(r,1),size(r,2),size(r,3)) !mass loss term of number due to intermodal coagulation [ug/m^3/s]
      real, intent(out) :: RIM(npop-1,NMSPCS,size(r,1),size(r,2),size(r,3)) !mass production terms due to intermodal coagulation. [ug/m^3/s]
      real :: KBAR0_IJ(npop-1,npop-1, size(r,1), size(r,2), size(r,3)) !mode averaged coagulation coefficient [m3/s]
      real :: KBAR3_IJ(npop-1,npop-1, size(r,1), size(r,2), size(r,3)) !mode averaged coagulation coefficient [m3/s]
      real :: dg_ip_um(size(r,1),size(r,2),size(r,3)), dg_jp_um(size(r,1),size(r,2),size(r,3)) !particle geometric diameter [um]
      real :: num_pop_k(size(r,1),size(r,2),size(r,3)), num_pop_l(size(r,1),size(r,2),size(r,3)), num_pop_i(size(r,1),size(r,2),size(r,3))
      real :: num_pop_j(size(r,1),size(r,2),size(r,3))
      real :: MJQ(npop-1,NMSPCS,size(r,1),size(r,2),size(r,3)) !mass of species in a single particle: mass/number = ug/particle
      real :: sig_ip, sig_jp
      integer :: ip,jp,kp,it,jt,kt,i,j,k,l,ikl,ipop,kpop,lpop,q,qq,klq
      !real, intent(out) :: 
      it = size(r,1)
      jt = size(r,2)
      kt = size(r,3)
      KBAR0_IJ = 0.
      KBAR3_IJ = 0.
      BI = 0.
      RI = 0. !production of number due to coagulation
      LI = 0.
      FI = 0.
      LIM = 0.
      RIM = 0.
      MJQ = 0.
      !constaruct MJQ
      !MJQ: MJQ(J,Q) is the avg. mass/particle of species Q (=1-5) for mode J
      do ip=1,npop-1 !exclude 'ext' pop
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
                num_pop_i = matrix_all_tracer(matrix_all_pop(ip)%I_N)%value_in_matrix(is:ie,js:je,:)
                num_pop_i = max(num_pop_i, 1E-15)
                do q = 1, NM(ip)
                        qq=prod_index(ip,q) !index of mass species: 1-5
                        if (mpp_root_pe().eq.mpp_pe()) then
                                write(*,*) 'ip, NM(ip), q, qq', ip, NM(ip), q, qq
                        endif
                        select case (qq)
                        case (1)
                                MJQ(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%I_MSULF)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                        case (2)
                                MJQ(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%I_MDUST)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                        case (3)
                                MJQ(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%I_MSEAS)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                        case (4)
                                MJQ(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%I_MOCAR)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                        case (5)
                                MJQ(ip,qq,:,:,:) = matrix_all_tracer(matrix_all_pop(ip)%I_MBCAR)%value_in_matrix(is:ie,js:je,:)/num_pop_i
                        case default !call fatal error
                                call ERROR_MESG('matrix coagulation','MJQ mass species not properly defined ', FATAL)
                        end select
                enddo
        endif
      enddo

      !set up KBAR0_IJ and KBAR3_IJ table (npop, npop, it, jt, kt)
      do ip = 1, npop-1 !exclude 'ext' pop
      do jp = 1, npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(jp)%nb_tracer_pop > 0)) then
                dg_ip_um = matrix_all_pop(ip)%Dg_wet(is:ie,js:je,:)*1.0e6 !3D array [um]
                sig_ip = matrix_all_pop(ip)%sigma !real number
                dg_jp_um = matrix_all_pop(jp)%Dg_wet(is:ie,js:je,:)*1.0e6 ! [um]
                sig_jp = matrix_all_pop(jp)%sigma
                do i=1,it
                do j=1,jt
                do k=1,kt
                       if (dg_ip_um(i,j,k) > 0 .and. dg_jp_um(i,j,k) > 0) then
                               !if (mpp_root_pe().eq.mpp_pe()) then
                               !        write(*,*) 'prior_test'
                               !        write(*,*) 'ip, jp', ip, jp
                               !        
                               !        write(*,*) 'i, j, k,t(i,j,k),  pfull(i,j,k), dg_ip(i,j,k), sig_ip, dg_jp(i,j,k), sig_jp'
                               !        write(*,*) i,j,k,t(i,j,k), pfull(i,j,k), dg_ip(i,j,k), sig_ip, dg_jp(i,j,k), sig_jp
                               !endif
                               call GET_KNIJ(t(i,j,k), pfull(i,j,k), dg_ip_um(i,j,k), sig_ip, dg_jp_um(i,j,k), sig_jp, KBAR0_IJ(ip, jp, i,j,k), &
                                       KBAR3_IJ(ip, jp, i,j,k)) !KBAR0_IJ/KBAR3_IJ mode averaged coagulation coefficient, [m3/s]
                               !if (mpp_root_pe().eq.mpp_pe()) then
                               !        write(*,*) 'behind_function_test'
                               !        write(*,*) 'i,j,k,t(i,j,k),  pfull(i,j,k), dg_ip(i,j,k), sig_ip, dg_jp(i,j,k), sig_jp, &
                               !        & KBAR0_IJ, KBAR3_IJ'
                               !        write(*,*) t(i,j,k), pfull(i,j,k), dg_ip(i,j,k), sig_ip, dg_jp(i,j,k), sig_jp, KBAR0_IJ(ip, &
                               !        & jp, i,j,k), KBAR3_IJ(ip, jp, i,j,k)
                               !endif
                       endif
                enddo
                enddo
                enddo
        endif
      enddo
      enddo
                  
      !--------------------------------------------------------
      ! calculate number production due to coagulation
      ! RI(npop-1, :, :,:) production of number for pop-i: #/m3/s
      !--------------------------------------------------------
      !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG 
      !        write(*,*) 'pass_test'
      !        write(*,*) 'nDIKL', nDIKL
      !endif

      do ikl = 1, nDIKL !ikl is 1-12, exclude ext already
        ipop = dikl_control(ikl)%i
        kpop = dikl_control(ikl)%k
        lpop = dikl_control(ikl)%l
        if (matrix_all_pop(ipop)%nb_tracer_pop > 0) then
                if ((matrix_all_pop(kpop)%nb_tracer_pop > 0) .and. (matrix_all_pop(lpop)%nb_tracer_pop > 0)) then
                        num_pop_k = matrix_all_tracer(matrix_all_pop(kpop)%I_N)%value_in_matrix(is:ie,js:je,:)
                        num_pop_l = matrix_all_tracer(matrix_all_pop(lpop)%I_N)%value_in_matrix(is:ie,js:je,:)
                        RI(ipop,:,:,:) = RI(ipop,:,:,:) + KBAR0_IJ(kpop,lpop,:,:,:) * num_pop_k * num_pop_l !production of number due to coagulation
                endif
        endif
      enddo
      
      !if (mpp_root_pe().eq.mpp_pe()) then
      !        write(*,*) 'debug RI and KBAR0_IJ'
      !        write(*,*) 'theoretical RI = 0'
      !        write(*,*) 'real RI max = ', maxval(RI)
      !endif

      !--------------------------------------------------------
      ! calculate number loss due to coagulation
      ! LI(npop-1, :, :,:) loss of number for pop-i: #/m3/s
      ! NOTE: DIJ is non-symmetric
      !--------------------------------------------------------
      do kp = 1,npop-1
      do ip = 1,npop-1
        if ((matrix_all_pop(ip)%nb_tracer_pop > 0) .and. (matrix_all_pop(kp)%nb_tracer_pop > 0)) then
                num_pop_k = matrix_all_tracer(matrix_all_pop(kp)%I_N)%value_in_matrix(is:ie,js:je,:)
                if (DIJ(ip, kp) > 0) then !i.e. DIJ(ip, kp) = 1
                        BI(ip,:,:,:) = BI(ip,:,:,:) + KBAR0_IJ(ip,kp,:,:,:)*num_pop_k !coefficient loss of mode I due to coagulation
                endif
        endif
      enddo
      enddo
      
      !if (mpp_root_pe().eq.mpp_pe()) then
      !        write(*,*) 'debug BI'
      !        write(*,*) 'theoretical BI = 0'
      !        write(*,*) 'real BI max = ', maxval(BI)
      !endif

      do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
                num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%I_N)%value_in_matrix(is:ie,js:je,:)
                LI(ip,:,:,:) = BI(ip,:,:,:) * num_pop_i + 0.5 * KBAR0_IJ(ip,ip,:,:,:) * num_pop_i * num_pop_i !loss of number due to coagulation #/m3/s
        endif
      enddo
     
      !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
      !        write(*,*) 'ratio of LI/num_pop_i^2/KBAR0_IJ'
      !        do ip =1, npop-1
      !        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
      !          do i=1,it
      !          do j=1,jt
      !          do k=1,kt 
      !                  if (LI(ip,i,j,k) > 0) then
      !                          num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%I_N)%value_in_matrix(is:ie,js:je,:)
      !                          write(*,*) 'ip, i,j k, LI(ip,i,j,k)', ip, i,j ,k, LI(ip,i,j,k)
      !                          write(*,*) 'ratio', LI(ip, i,j,k)/(KBAR0_IJ(ip,ip,i,j,k) * num_pop_i(i,j,k) * num_pop_i(i,j,k))
      !                  endif
      !          enddo
      !          enddo
      !          enddo         

      !        endif
      !        enddo
      !endif

      !--------------------------------------------------------
      ! calculate mass production due to coagulation
      ! RIM(npop-1,nmspcs, :, :,:) production of mass for pop-i, spec-qq :ug/m3/s
      !--------------------------------------------------------
      do ip =1, npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0 .and. giklq_control(ip)%n > 0) then
                do klq = 1, giklq_control(ip)%n
                        k = giklq_control(ip)%k(klq)
                        l = giklq_control(ip)%l(klq)
                        qq = giklq_control(ip)%qq(klq)
                        if (matrix_all_pop(k)%nb_tracer_pop > 0 .and. matrix_all_pop(l)%nb_tracer_pop > 0 ) then
                                num_pop_k = matrix_all_tracer(matrix_all_pop(k)%I_N)%value_in_matrix(is:ie,js:je,:)
                                num_pop_l = matrix_all_tracer(matrix_all_pop(l)%I_N)%value_in_matrix(is:ie,js:je,:)
                                RIM(ip,qq,:,:,:) = RIM(ip,qq,:,:,:) + num_pop_k*num_pop_l*KBAR3_IJ(l,k,:,:,:)*MJQ(l,qq,:,:,:)
                        endif
                enddo
        endif
     enddo

     !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
     !           write(*,*) 'RIM THEORETICAL = 0'
     !           write(*, *) 'real value RIM', maxval(RIM)
     !endif

      !--------------------------------------------------------
      ! calculate mass loss due to coagulation
      ! LIM(npop-1,nmspcs, :, :,:) loss of mass for pop-i, spec-qq :ug/m3/s
      !--------------------------------------------------------
      do ip = 1,npop-1
        if (matrix_all_pop(ip)%nb_tracer_pop > 0) then 
        do jp = 1,npop-1
                if (matrix_all_pop(jp)%nb_tracer_pop > 0) then
                        num_pop_j =  matrix_all_tracer(matrix_all_pop(jp)%I_N)%value_in_matrix(is:ie,js:je,:)
                        IF( DIJ(ip,jp) > 0 ) FI(ip,:,:,:) = FI(ip,:,:,:) + KBAR3_IJ(ip,jp,:,:,:)*num_pop_j
                endif
        enddo
        endif
      enddo

      DO ip =1, npop-1
      if (matrix_all_pop(ip)%nb_tracer_pop > 0) then
              DO Q=1, NM(ip)
                qq = prod_index(ip,q)
                num_pop_i =  matrix_all_tracer(matrix_all_pop(ip)%I_N)%value_in_matrix(is:ie,js:je,:)
                LIM(ip,qq,:,:,:) = LIM(ip,qq,:,:,:) + FI(ip,:,:,:) * num_pop_i * MJQ(ip,qq,:,:,:)
              ENDDO
      endif
      ENDDO

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
!        ! 1/6*PI*D_DRY^3 = sum(volume_spec,i)       
!----------------------------------------------------------------
subroutine matrix_dry_diameter(r,is,ie,js,je)
       real, intent(in) :: r(:,:,:,:)
       integer :: n, nt, ntt, tr_index !local variables
       real :: m_dry_all(size(r,1),size(r,2),size(r,3)), vol_dry_all(size(r,1),size(r,2),size(r,3))
       real :: kappa_vol_all(size(r,1),size(r,2),size(r,3)), vol_spec(size(r,1),size(r,2),size(r,3))
       real :: pop_num(size(r,1),size(r,2),size(r,3)), dens_dry(size(r,1),size(r,2),size(r,3))
       real :: exp_fac
       integer, intent(in) :: is,ie,js,je
       integer :: flag,it,jt,kt, i, j, k !flag for whether this population has dry mass/perform calculation
       integer :: step = 0 !count for debug: which step crash
       it = size(r,1)
       jt = size(r,2)
       kt = size(r,3)
       do n = 1, npop
           m_dry_all = 0.
           vol_dry_all = 0.
           kappa_vol_all = 0.
           flag = 0
           if (matrix_all_pop(n)%nb_tracer_pop > 0) then
               nt = matrix_all_pop(n)%nb_tracer_pop
               do ntt = 1, nt
                  tr_index = matrix_all_pop(n)%tracer_index(ntt)
                  if ((lowercase(trim(matrix_all_tracer(tr_index)%spec)) .ne. "alwc") .AND. &
                        (lowercase(trim(matrix_all_tracer(tr_index)%type)) .ne. "number") ) then
                      m_dry_all = m_dry_all + matrix_all_tracer(tr_index)%value_in_matrix(is:ie,js:je,:) !m_dry unit: ug/m3
                      !vol_spec unit: m3 aerosol_volume/ m3 air
                      vol_spec = m_dry_all/matrix_all_tracer(tr_index)%dens*1E-9 
                      vol_dry_all = vol_dry_all + vol_spec
                      kappa_vol_all = kappa_vol_all + matrix_all_tracer(tr_index)%kappa*vol_spec !sum of kappa*vol
                      flag = 1 !record if this population has dry mass
                  endif
               enddo
               if (flag > 0) then
                   matrix_all_pop(n)%mass_dry(is:ie,js:je,:) = m_dry_all !mass unit: ug/m3
                   matrix_all_pop(n)%vol_dry(is:ie,js:je,:) = vol_dry_all !volume unit: m3 aerosol/m3
                   pop_num = matrix_all_tracer(matrix_all_pop(n)%I_N)%value_in_matrix(is:ie,js:je,:) !#/m3
                   exp_fac = exp(1.5*(log(matrix_all_pop(n)%sigma))**2) 
                   do i=1,it
                        do j=1,jt
                                do k=1,kt
                                        if ((pop_num(i,j,k) > 1E-14) .AND. (m_dry_all(i,j,k) > 1E-32)) then !not zero number on the grid
                                            dens_dry(i,j,k) = m_dry_all(i,j,k)/vol_dry_all(i,j,k)*1E-9
                                            matrix_all_pop(n)%dens_dry(i+is-1,j+js-1,k) = dens_dry(i,j,k)
                                            matrix_all_pop(n)%kappa_pop(i+is-1,j+js-1,k) = kappa_vol_all(i,j,k)/vol_dry_all(i,j,k)
                                            matrix_all_pop(n)%Dg_dry(i+is-1,j+js-1,k) = &
                                            (m_dry_all(i,j,k)/pop_num(i,j,k)/dens_dry(i,j,k)*1E-9/PI6)**(1.0/3)/exp_fac ! Dg_dry unit: m
                                        else
                                            dens_dry(i,j,k) = 0
                                            matrix_all_pop(n)%dens_dry(i+is-1,j+js-1,k) = dens_dry(i,j,k)
                                            matrix_all_pop(n)%kappa_pop(i+is-1,j+js-1,k)=0
                                            matrix_all_pop(n)%Dg_dry(i+is-1,j+js-1,k) = 0
                                        endif
                                enddo
                        enddo
                   enddo
               endif
           endif 
       end do  
       step = step + 1
end subroutine
        


subroutine matrix_wet_diameter(r, rh, t, is,ie,js,je) !calculate and assgin values for pop%Dg_wet
       real, intent(in) :: r(:,:,:,:)
       real, intent(in) :: rh(:,:,:), t(:,:,:)
       real:: ddry_in_3d(size(r,1),size(r,2),size(r,3)), hygro_3d(size(r,1),size(r,2),size(r,3))
       integer :: n, nt, ntt, i, j, k, it, jt, kt, tr_index
       integer,intent(in) :: is,ie,js,je
       real :: ddry_in, hygro_in, s_in, tair_in, dwet_out, gf, rh_deliquescence, rh_crystallization
       real :: particle_vol_dry, particle_vol_water, particle_vol_wet, f_hysteresis, gf3 
       integer :: flag = 0 ! check if the population get calculated 
       it = size(r,1)
       jt = size(r,2)
       kt = size(r,3)
       do n = 1, npop
           flag = 0
           if (matrix_all_pop(n)%nb_tracer_pop > 0) then
               nt = matrix_all_pop(n)%nb_tracer_pop
               do ntt = 1, nt
                  tr_index = matrix_all_pop(n)%tracer_index(ntt)
                  if ((lowercase(trim(matrix_all_tracer(tr_index)%spec)) .ne. "alwc") .AND. &
                        (lowercase(trim(matrix_all_tracer(tr_index)%type)) .ne. "number") ) then
                        flag = 1
                  endif
                enddo
           endif
           if (flag > 0) then
                ddry_in_3d = matrix_all_pop(n)%Dg_dry(is:ie,js:je,:)*exp(1.5*(log(matrix_all_pop(n)%sigma))**2) !volume mean diameter
                hygro_3d = matrix_all_pop(n)%kappa_pop(is:ie,js:je,:)
                rh_deliquescence = matrix_all_pop(n)%rh_deliquescence
                rh_crystallization = matrix_all_pop(n)%rh_crystallization
                do i=1,it
                        do j=1,jt
                                do k=1,kt
                                        ddry_in = max(0.0,ddry_in_3d(i,j,k))
                                        particle_vol_dry = PI6*ddry_in**3
                                        hygro_in = hygro_3d(i,j,k)
                                        s_in = max(0.0, rh(i,j,k))
                                        s_in = min(1.0, s_in)
                                        tair_in = t(i,j,k)
                                        gf = 1.0 !default value
                                        if (ddry_in > 1E-10) then
                                                call aero_kohler(ddry_in, hygro_in, s_in, tair_in, dwet_out)
                                                dwet_out = max(ddry_in, dwet_out) !dwet_out and ddry_in are volume mean diameter
                                                particle_vol_water = PI6*(dwet_out**3 - ddry_in**3)
                                                !check current RH with pop rh crystalization and deliquescence
                                                if (s_in < rh_crystallization) then
                                                        dwet_out = ddry_in
                                                        particle_vol_water = 0
                                                elseif (s_in < rh_deliquescence) then
                                                        f_hysteresis = 1.0 / max(1.0e-5, (rh_deliquescence - rh_crystallization))
                                                        particle_vol_water = f_hysteresis * (s_in - rh_crystallization) * particle_vol_water
                                                        particle_vol_water = max(0.0, particle_vol_water)
                                                        particle_vol_wet = particle_vol_dry + particle_vol_water
                                                        dwet_out = (particle_vol_wet / PI6)**(1.0/3)
                                                end if
                                                gf=max(dwet_out/ddry_in,1.0) !growth factor
                                        endif
                                        matrix_all_pop(n)%Dg_wet(i+is-1,j+js-1,k) = gf*matrix_all_pop(n)%Dg_dry(i+is-1,j+js-1,k)
                                        gf3 = gf**3
                                        matrix_all_pop(n)%dens_wet(i+is-1,j+js-1,k) = (1000*(gf3-1) + matrix_all_pop(n)%dens_dry(i+is-1,j+js-1,k))/gf3 !wet density
                                enddo
                        enddo
                enddo
          endif
       enddo
end subroutine


!----------------------------------------------------------------
!        subroutine matrix_intermodal_transfer
! perform intermodal transfer : 
!       (1) AKK -> ACC: aitken mode sulfate -> accumulation mode sulfate
!       (2) OC1 -> OC2: hydrophobic organic carbon -> hydrophilic organic carbon
!       (3) BC1 -> BC2: hydrophobic black carbon -> hydrophilic carbon
!!!!!!!!!!NOTE: matrix value got updated in this subroutine
!--------------------------------------------------------------------
subroutine matrix_intermodal_transfer(config, r, is,ie,js,je)
        character(len=32), intent(in) :: config
        real, intent(in) :: r(:,:,:,:)
        integer, intent(in) :: is,ie,js,je
        ! These variables are used in intermodal transfer.
        REAL(8), PARAMETER :: DG_AKK_param      = 0.026D+00      !Giss set-up
        REAL(8), PARAMETER :: DG_ACC_param      = 0.110D+00      ! E04, Table 2, accumulation mode 
        INTEGER, PARAMETER :: IMTR_EXP          = 4              ! exponent for AKK --> ACC intermodal transfer
        real, parameter    :: IMTR_METHOD       = 1     ! =1 no cut of pdf, =2 fixed-Dp cut, =3 variable-Dp cut as in CMAQ
        REAL(8)            :: XNUM, X3                  ! error function complement arguments [1]
        REAL(8)            :: DGN_AKK_IMTR              ! geo. mean diam. of the AKK mode number distribution [m]
        REAL(8)            :: DGN_ACC_IMTR              ! geo. mean diam. of the ACC mode number distribution [m]
        REAL(8)            :: DEL_NUMB                  ! number conc. transferred from AKK to ACC [ #/m^3]
        REAL(8)            :: DEL_MASS                  ! mass   conc. transferred from AKK to ACC [ug/m^3]
        REAL(8), SAVE      :: DPCUT_IMTR      = 0.0D+00 ! fixed diameter for intermodal transfer [um]
        REAL(8), SAVE      :: XNUM_FACTOR     = 0.0D+00 ! factor in XNUM expression [1]
        REAL(8), SAVE      :: X3_TERM         = 0.0D+00 ! term   in X3   expression [1]
        REAL(8), SAVE      :: LNAKK_SIGMA     = 0.0D+00 ! ln(SG_AKK) [1]
        REAL(8), SAVE      :: LNACC_SIGMA     = 0.0D+00 ! ln(SG_ACC) [1]
        real, PARAMETER    :: DNU             = 3.0D-09   ! diameter of a new particle [nm]
        REAL(8), PARAMETER :: DPAKK0          = DNU     ! min. diameter of average mass for mode AKK [m]
        REAL(8), PARAMETER :: FNUM_MAX        = 0.5D+00 ! max. value of FNUM [1]
        REAL(8), PARAMETER :: AKK_MINNUM_IMTR = 1.0D+06 ! min. AKK number conc. to enable IMTR [#/m^3]
        real               :: DPAKK_all(size(r,1),size(r,2),size(r,3)) ! diameter of average mass for mode AKK [m]
        real               :: DPACC_all(size(r,1),size(r,2),size(r,3)) !diameter of average mass for mode ACC [m]
        real               :: DGAKK_all(size(r,1),size(r,2),size(r,3)) !geometric mean diameter [m]
        real               :: DGACC_all(size(r,1),size(r,2),size(r,3))
        real               :: ACC_num_all(size(r,1),size(r,2),size(r,3))
        real               :: AKK_num_all(size(r,1),size(r,2),size(r,3))
        real               :: DPAKK, DPACC
        REAL(8)            :: FNUM_all(size(r,1),size(r,2),size(r,3))! fraction of AKK number transferred over the time step [1]
        REAL(8)            :: F3_all(size(r,1),size(r,2),size(r,3))  ! fraction of AKK mass   transferred over the time step [1]
        real               :: FNUM, F3
        real, parameter :: MIMR_BC1 = 0.10 !threshhold of sulfate/bc mass ratio to transfer BC1 to BC2
        real, parameter :: MIMR_OC1 = 0.10 !threshhold of sulfate/oc mass ratio to transfer OC1 to OC2
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
        FNUM = 0.
        F3 = 0.
        DPAKK = 0.
        DPACC = 0.
        it = size(r,1)
        jt = size(r,2)
        kt = size(r,3) 
        !BC1 -> BC2: hydrophobic to hydrophilic
        dmsulf_bc1_to_bc2 = 0.
        dmbcar_bc1_to_bc2 = 0.
        dn_bc1_to_bc2 = 0.
        !OC1 -> OC2: hydrophobic to hydrophilic
        dmsulf_oc1_to_oc2 = 0
        dmocar_oc1_to_oc2 = 0.
        dn_oc1_to_oc2 = 0.
        !AKK -> ACC: aitken mode to accumulation mode
        dmsulf_akk_to_acc = 0.
        dn_akk_to_acc = 0.
        bc1_sulf_mass = 0.
        oc1_sulf_mass = 0.
        bc1_bcar_mass = 0.
        bc1_num = 0.
        oc1_ocar_mass = 0.
        oc1_num =0.
        if (trim(lowercase(config)) .eq. "full") then
                !-------------------------------------------------------------------
                ! Transfer mode BC1 to BC2
                !-------------------------------------------------------------------
                if (I_BC1 > 0 .AND. I_BC2 > 0 ) then
                        if (matrix_all_pop(I_BC1)%nb_tracer_pop > 0 .and. matrix_all_pop(I_BC2)%nb_tracer_pop > 0) then
                                bc1_sulf_mass = matrix_all_tracer(matrix_all_pop(I_BC1)%I_MSULF)%value_in_matrix(is:ie,js:je,:) !sulf mass in bc1 [ug/m3]
                                bc1_bcar_mass = matrix_all_tracer(matrix_all_pop(I_BC1)%I_MBCAR)%value_in_matrix(is:ie,js:je,:) !bc mass in bc1 [ug/m3]
                                bc1_num = matrix_all_tracer(matrix_all_pop(I_BC1)%I_N)%value_in_matrix(is:ie,js:je,:) !number of bc1 [#/m3]
                                bc2_sulf_mass = matrix_all_tracer(matrix_all_pop(I_BC2)%I_MSULF)%value_in_matrix(is:ie,js:je,:) !sulf mass in bc2 [ug/m3]
                                bc2_bcar_mass = matrix_all_tracer(matrix_all_pop(I_BC2)%I_MBCAR)%value_in_matrix(is:ie,js:je,:) !bc mass in bc2 [ug/m3]
                                bc2_num = matrix_all_tracer(matrix_all_pop(I_BC2)%I_N)%value_in_matrix(is:ie,js:je,:) !number of bc2 [#/m3]
                                do i = 1, it
                                do j = 1, jt
                                do k = 1, kt
                                if (bc1_bcar_mass(i,j,k) > 0.) then
                                        if (bc1_sulf_mass(i,j,k)/bc1_bcar_mass(i,j,k) .gt. MIMR_BC1) then
                                                dmsulf_bc1_to_bc2(i,j,k) = transfer_factor*bc1_sulf_mass(i,j,k)
                                                dmbcar_bc1_to_bc2(i,j,k) = transfer_factor*bc1_bcar_mass(i,j,k)
                                                dn_bc1_to_bc2(i,j,k) = transfer_factor*bc1_num(i,j,k)
                                        endif
                                endif
                                enddo
                                enddo
                                enddo
                                !update matrix value for BC1 pop
                                matrix_all_tracer(matrix_all_pop(I_BC1)%I_MSULF)%value_in_matrix(is:ie,js:je,:) = bc1_sulf_mass - dmsulf_bc1_to_bc2
                                matrix_all_tracer(matrix_all_pop(I_BC1)%I_MBCAR)%value_in_matrix(is:ie,js:je,:) = bc1_bcar_mass - dmbcar_bc1_to_bc2
                                matrix_all_tracer(matrix_all_pop(I_BC1)%I_N)%value_in_matrix(is:ie,js:je,:) = bc1_num - dn_bc1_to_bc2
                                !update matrix value for BC2 pop
                                matrix_all_tracer(matrix_all_pop(I_BC2)%I_MSULF)%value_in_matrix(is:ie,js:je,:) = bc2_sulf_mass + dmsulf_bc1_to_bc2
                                matrix_all_tracer(matrix_all_pop(I_BC2)%I_MBCAR)%value_in_matrix(is:ie,js:je,:) = bc2_bcar_mass + dmbcar_bc1_to_bc2
                                matrix_all_tracer(matrix_all_pop(I_BC2)%I_N)%value_in_matrix(is:ie,js:je,:) = bc2_num + dn_bc1_to_bc2

                        else
                                call error_mesg ('matrix_gfdl','BC1 and BC2: population and tracers not consistent', FATAL)
                        endif
                endif
                !-------------------------------------------------------------------
                ! Transfer mode OC1 to OC2
                !-------------------------------------------------------------------
                if (I_OC1 > 0 .AND. I_OC2 > 0 ) then
                        if (matrix_all_pop(I_OC1)%nb_tracer_pop > 0 .and. matrix_all_pop(I_OC2)%nb_tracer_pop > 0) then
                                oc1_sulf_mass = matrix_all_tracer(matrix_all_pop(I_OC1)%I_MSULF)%value_in_matrix(is:ie,js:je,:) !sulf mass in oc1 [ug/m3]
                                oc1_ocar_mass = matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%value_in_matrix(is:ie,js:je,:) !oc mass in oc1 [ug/m3]
                                oc1_num = matrix_all_tracer(matrix_all_pop(I_OC1)%I_N)%value_in_matrix(is:ie,js:je,:) !number of oc1 [#/m3]
                                do i = 1, it
                                do j = 1, jt
                                do k = 1, kt
                                if (oc1_ocar_mass(i,j,k) > 0) then
                                        if (oc1_sulf_mass(i,j,k)/oc1_ocar_mass(i,j,k) .gt. MIMR_OC1) then
                                                dmsulf_oc1_to_oc2(i,j,k) = transfer_factor*oc1_sulf_mass(i,j,k)
                                                dmocar_oc1_to_oc2(i,j,k) = transfer_factor*oc1_ocar_mass(i,j,k)
                                                dn_oc1_to_oc2(i,j,k) = transfer_factor*oc1_num(i,j,k)
                                        endif
                                endif
                                enddo
                                enddo
                                enddo
                                !update matrix value for OC1 pop
                                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MSULF)%value_in_matrix(is:ie,js:je,:) = oc1_sulf_mass - dmsulf_oc1_to_oc2
                                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%value_in_matrix(is:ie,js:je,:) = oc1_ocar_mass - dmocar_oc1_to_oc2
                                matrix_all_tracer(matrix_all_pop(I_OC1)%I_N)%value_in_matrix(is:ie,js:je,:) = oc1_num - dn_oc1_to_oc2
                                !update matrix value for OC2 pop
                                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MSULF)%value_in_matrix(is:ie,js:je,:) = oc2_sulf_mass + dmsulf_oc1_to_oc2
                                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%value_in_matrix(is:ie,js:je,:) = oc2_ocar_mass + dmocar_oc1_to_oc2
                                matrix_all_tracer(matrix_all_pop(I_OC2)%I_N)%value_in_matrix(is:ie,js:je,:) = oc2_num + dn_oc1_to_oc2
                        else
                                call error_mesg ('matrix_gfdl','OC1 and OC2: population and tracers not consistent', FATAL)
                        endif
                endif
                !-------------------------------------------------------------------
                ! Transfer mode AKK to ACC
                !-------------------------------------------------------------------
                if (I_AKK > 0 .AND. I_ACC > 0) then
                        if (matrix_all_pop(I_AKK)%nb_tracer_pop > 0 .and. matrix_all_pop(I_ACC)%nb_tracer_pop > 0) then
                                LNAKK_SIGMA = log(matrix_all_pop(I_AKK)%sigma)
                                LNACC_SIGMA = log(matrix_all_pop(I_ACC)%sigma)
                                XNUM_FACTOR = 1.0D+00 / ( SQRT( 2.0D+00 ) * LNAKK_SIGMA )
                                X3_TERM     = 3.0D+00 * LNAKK_SIGMA / SQRT( 2.0D+00 )
                                DPCUT_IMTR  = SQRT( DG_AKK_param * DG_ACC_param )   ! This is the formula of Easter et al. 2004. (= 0.053 um)
                                DPAKK_all = matrix_all_pop(I_AKK)%Dg_dry(is:ie,js:je,:) * exp(1.5* LNAKK_SIGMA**2) !diameter of average mass for AKK [m]
                                DPACC_all = matrix_all_pop(I_ACC)%Dg_dry(is:ie,js:je,:) * exp(1.5* LNACC_SIGMA**2) !dry diameter of average mass for ACC [m]
                                DGAKK_all = matrix_all_pop(I_AKK)%Dg_dry(is:ie,js:je,:) !geometric diameter [m]
                                DGACC_all = matrix_all_pop(I_ACC)%Dg_dry(is:ie,js:je,:) !geometric diameter [m]
                                ACC_num_all = matrix_all_tracer(matrix_all_pop(I_ACC)%I_N)%value_in_matrix(is:ie,js:je,:)
                                AKK_num_all = matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%value_in_matrix(is:ie,js:je,:)
                                akk_sulf_mass = matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%value_in_matrix(is:ie,js:je,:) !sulf mass in akk [ug/m3]
                                acc_sulf_mass = matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%value_in_matrix(is:ie,js:je,:) !sulf mass in acc [ug/m3]
                                do i = 1, it
                                do j = 1, jt
                                do k = 1, kt
                                DPAKK = DPAKK_all(i,j,k)
                                DPACC = DPACC_all(i,j,k)
                                if ((DPAKK > 0) .and. (DPACC > DPAKK)) then
                                        IF( IMTR_METHOD .EQ. 1 ) THEN
                                                !------------------------------------------------------------------------------------------------------------
                                                ! Calculate the fraction transferred based on the relative difference
                                                !   in mass mean diameters of the AKK and ACC modes.
                                                !------------------------------------------------------------------------------------------------------------
                                                IF( (DPAKK .GE. DPAKK0) .and. (DPACC .GE. DPAKK0)) THEN                ! [m], DPAKK0 = 3 nm
                                                        FNUM = ( ( DPAKK - DPAKK0 ) / ( DPACC - DPAKK0 ) )**IMTR_EXP   ! fraction transferred from AKK to ACC
                                                        FNUM = MAX( MIN( FNUM, FNUM_MAX ), 0.0D+00 )                   ! limit transfer in a single transfer
                                                ELSE
                                                        FNUM = 0.0D+00
                                                ENDIF
                                                F3 = FNUM
                                                ! WRITE(34,'(7D12.4)')FNUM,F3
                                        ELSEIF( IMTR_METHOD .EQ. 2 ) THEN
                                                !------------------------------------------------------------------------------------------------------------
                                                ! Calculate the fraction transferred based on a fixed
                                                !   threshold diameter DPCUT_IMTR.
                                                !------------------------------------------------------------------------------------------------------------
                                                DGN_AKK_IMTR = 1.0D+06 * DGAKK_all(i,j,k)                      ! geometric diameter in [um]
                                                XNUM = XNUM_FACTOR * LOG( DPCUT_IMTR / DGN_AKK_IMTR )          ! [1]
                                                XNUM = MAX( XNUM, X3_TERM )                                    ! limit for stability as in BS2003
                                                X3 = XNUM - X3_TERM                                            ! [1]
                                                FNUM = 0.5D+00 * ERFC( XNUM )                                  ! number fraction transferred from AKK to ACC
                                                F3   = 0.5D+00 * ERFC( X3   )                                  ! mass   fraction transferred from AKK to ACC
                                                ! WRITE(34,'(9D12.4)')DGN_AKK_IMTR,DPCUT_IMTR,DPAKK*1.0D+06,AERO(NUMB_AKK_1),AERO(MASS_AKK_SULF),FNUM,F3
                                        ELSEIF( IMTR_METHOD .EQ. 3 ) THEN
                                                !------------------------------------------------------------------------------------------------------------
                                                ! Calculate the fraction transferred based on the
                                                !   diameter of intersection of the AKK and ACC modes.
                                                !------------------------------------------------------------------------------------------------------------
                                                DGN_AKK_IMTR = 1.0D+06 * DGAKK_all(i,j,k)                      ! [um]
                                                DGN_ACC_IMTR = 1.0D+06 * DGACC_all(i,j,k)                      ! [um]
                                                IF( ACC_num_all(i,j,k) .GT. 1.0D+06 ) THEN
                                                        XNUM = GETXNUM(AKK_num_all(i,j,k), ACC_num_all(i,j,k), &
                                                                DGN_AKK_IMTR,DGN_ACC_IMTR,LNAKK_SIGMA,LNACC_SIGMA)   ! [1]
                                                ELSE                                                           ! mode ACC essentially empty - use Method 2
                                                        XNUM = XNUM_FACTOR * LOG( DPCUT_IMTR / DGN_AKK_IMTR )        ! [1]
                                                ENDIF
                                                XNUM = MAX( XNUM, X3_TERM )                                    ! limit for stability as in BS2003
                                                X3 = XNUM - X3_TERM                                            ! [1]
                                                FNUM = 0.5D+00 * ERFC( XNUM )                                  ! number fraction transferred from AKK to ACC
                                                F3   = 0.5D+00 * ERFC( X3   )                                  ! mass   fraction transferred from AKK to ACC
                                        ENDIF
                                endif
                                FNUM_all(i,j,k) = FNUM
                                F3_all(i,j,k) = F3
                                enddo
                                enddo
                                enddo

                                dmsulf_akk_to_acc = F3_all * akk_sulf_mass  !mass concentration transferred [ug/m3] 
                                dn_akk_to_acc = FNUM_all * AKK_num_all
                                !update matrix value for AKK and ACC
                                matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%value_in_matrix(is:ie,js:je,:) = AKK_num_all - dn_akk_to_acc
                                matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%value_in_matrix(is:ie,js:je,:) = akk_sulf_mass - dmsulf_akk_to_acc
                                matrix_all_tracer(matrix_all_pop(I_ACC)%I_N)%value_in_matrix(is:ie,js:je,:) = ACC_num_all + dn_akk_to_acc
                                matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%value_in_matrix(is:ie,js:je,:) = acc_sulf_mass + dmsulf_akk_to_acc
                        endif
                else
                        call error_mesg ('matrix_gfdl','AKK and ACC: population and tracers not consistent', FATAL)
                endif

        else
                call error_mesg ('matrix_gfdl','intermodal_transfer configuration not found', FATAL)
        endif

end subroutine

subroutine query_matrix_tracer(ind_tracer, flag) !check if a tracer is matrix tracer
        integer, intent(in) :: ind_tracer !index of tracer from query method
        integer, intent(out) :: flag  !no_active -1; active_mass_tracer = 1; active_number_tracer = 2
        integer :: flag_n = 2, flag_m = 1
        flag = -1
        if (matrix_all_tracer(ind_tracer)%is_active) then
                if (lowercase(trim(matrix_all_tracer(ind_tracer)%type)) .eq. 'mass') then
                        flag = flag_m
                elseif (lowercase(trim(matrix_all_tracer(ind_tracer)%type)) .eq. 'number') then
                        flag = flag_n
                endif
        endif
end subroutine

!The following subroutines are used in aerosol.F90
!------------------------------start of aerosol.F90 call-------------------------------------
subroutine query_matrix_info(npop_active, nspcs)
        integer, intent(out) :: npop_active, nspcs
        integer :: n_pop_count = 0

        !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
        !        write(*,*) 'query_matrix_debug'
        !        write(*,*) 'do_matrix', do_matrix
        !        write(*,*) 'matrix_module_init', matrix_module_init
        !        write(*,*) 'matrix_configuration', matrix_configuration
        !        write(*,*) 'matrix_configuration==debug', (matrix_configuration.eq.'debug')
        !endif

        if (matrix_module_init) then
                if (do_matrix) then
                        npop_active = pop_active
                        nspcs = 5 
                else
                        npop_active = 0
                        nspcs = 0
                endif
        else !matrix module hasn't been initialized when calling the function
                !read namelist for do_matrix
                if(file_exist('input.nml')) then
                        #ifdef INTERNAL_FILE_NML
                        read (input_nml_file, nml = matrix_nml, iostat = io)
                        ierr = check_nml_error(io,'matrix_nml')
                        #else
                        unit = open_namelist_file('input.nml')
                        ierr=1; do while (ierr /= 0)
                        read(unit, nml = matrix_nml, iostat = io, end = 10)
                        ierr = check_nml_error (io, 'matrix_nml')
                        end do
                        10  call close_file(unit)
                        #endif
                endif
                logunit = stdlog()
                if(mpp_pe() == mpp_root_pe()) then
                        write(logunit, nml=matrix_nml)
                        verbose = verbose + 1
                end if
                !generic function to assign configuration population
                if (do_matrix) then
                        call set_config_pop(matrix_configuration)
                        npop_active = pop_active
                        nspcs = 5
                else
                        npop_active = 0
                        nspcs = 0
                endif
        end if

end subroutine query_matrix_info


subroutine query_pop_number(is, ie, js, je, i_order, i_abs_pop, pop_number)
        integer, intent(in) :: is, ie, js, je !physical window
        integer, intent(in) :: i_order !index of activated population, e.g. debug, i_order =1
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = I_ACC
        real, intent(out) :: pop_number(is:ie,js:je, kd) !number concentration of the current population
        integer :: pop_index
        pop_index = I_POP(i_order)
        i_abs_pop = pop_index
        pop_number = MAX(0.0, matrix_all_tracer(matrix_all_pop(pop_index)%I_N)%value_in_matrix(is:ie,js:je,:)) !#/m3
end subroutine


subroutine query_pop_Dg_dry(is, ie, js, je, i_order, i_abs_pop, pop_Dg_dry)
        integer, intent(in) :: is, ie, js, je !physical window
        integer, intent(in) :: i_order !index of activated population, e.g. debug, i_order =1
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = I_ACC
        real, intent(out) :: pop_Dg_dry(is:ie,js:je, kd) !number concentration of the current population
        integer :: pop_index
        pop_index = I_POP(i_order)
        i_abs_pop = pop_index
        pop_Dg_dry = MAX(0.0, matrix_all_pop(pop_index)%Dg_dry(is:ie,js:je,:)) !m
end subroutine

subroutine query_pop_MSPCS(is, ie, js, je, nspcs,  i_order, i_abs_pop, pop_MSPCS)
        integer, intent(in) :: is, ie, js, je, nspcs !physical window
        integer, intent(in) :: i_order !index of activated population, e.g. debug, i_order =1
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = I_ACC
        real, intent(out) :: pop_MSPCS(is:ie,js:je, kd, 1, nspcs) !mass concentration of each species in the population, [ug/m3]
        integer :: pop_index
        pop_index = I_POP(i_order)
        i_abs_pop = pop_index
        pop_MSPCS = 0.0
        ! nspcs = 5, in the order of SULF, BCAR, OCAR, DUST, SEAS
        if (matrix_all_pop(pop_index)%I_MSULF > 0) then
                pop_MSPCS(:,:,:, 1, 1) = MAX(0.0, &
                        matrix_all_tracer(matrix_all_pop(pop_index)%I_MSULF)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%I_MBCAR > 0) then
                pop_MSPCS(:,:,:, 1, 2) = MAX(0.0, &
                        matrix_all_tracer(matrix_all_pop(pop_index)%I_MBCAR)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%I_MOCAR > 0) then
                pop_MSPCS(:,:,:, 1, 3) = MAX(0.0, &
                        matrix_all_tracer(matrix_all_pop(pop_index)%I_MOCAR)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%I_MDUST > 0) then
                pop_MSPCS(:,:,:, 1, 4) = MAX(0.0, &
                        matrix_all_tracer(matrix_all_pop(pop_index)%I_MDUST)%value_in_matrix(is:ie,js:je,:))
        elseif (matrix_all_pop(pop_index)%I_MSEAS > 0) then
                pop_MSPCS(:,:,:, 1, 5) = MAX(0.0, &
                        matrix_all_tracer(matrix_all_pop(pop_index)%I_MSEAS)%value_in_matrix(is:ie,js:je,:))
        endif

end subroutine

subroutine query_pop_sigma(i_order, i_abs_pop, pop_sigma)
        integer, intent(in) :: i_order
        integer, intent(out) :: i_abs_pop !index of absolute number, e.g. debug, i_abs_pop = I_ACC
        real, intent(out) :: pop_sigma
        integer :: pop_index
        pop_index = I_POP(i_order)
        i_abs_pop = pop_index
        pop_sigma = matrix_all_pop(pop_index)%sigma
end subroutine





!---------------------------------end of aerosol.F90 call ------------------------------------

subroutine query_matrix_pop(ind_tracer, flag, is, ie, js, je, pop_index, pop_number, pop_Dg_wet, pop_dens_wet)
        integer, intent(in) :: ind_tracer !index of tracer from query method
        integer, intent(in) :: is,ie,js,je !note: kd = size(r,3)
        integer, intent(out) :: flag  !index of population, if active and belong to matrix tracer pop
        ! number tracer: return 2, mass_tarcer: return 1, otherwise return -1
        integer, intent(out), optional :: pop_index
        real, intent(out), optional :: pop_number(ie-is+1, je-js+1, kd)
        real, intent(out), optional :: pop_Dg_wet(ie-is+1, je-js+1, kd)
        real, intent(out), optional :: pop_dens_wet(ie-is+1, je-js+1, kd) 
        integer :: flag_n = 2, flag_m = 1
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
                        pop_number = matrix_all_tracer(matrix_all_pop(pop_index)%I_N)%value_in_matrix(is:ie,js:je,:) !#/m3
                        pop_Dg_wet = matrix_all_pop(pop_index)%Dg_wet(is:ie,js:je,:) !unit m
                        pop_dens_wet = matrix_all_pop(pop_index)%dens_wet(is:ie,js:je,:) !unit: Kg/m3
                endif

        else
                flag= -1
        endif
end subroutine


subroutine matrix_sedimentation(T, P, zhalf, dt, is, ie, js, je)
        real, intent(in) :: T(:,:,:), P(:,:,:), zhalf(:,:,:) !unit K, Pa, m
        real, intent(in) :: dt
        integer, intent(in) :: is, ie, js, je
        real :: pop_vsed(size(p,1), size(p,2),size(p,3)) !gravititional settling velocity
        integer :: nt, n
        real :: rwet(size(p,1), size(p,2),size(p,3)) !wet radius m
        real :: rho_wet(size(p,1), size(p,2),size(p,3)) ! wet density, kg/m3
        real :: viscosity(size(p,1), size(p,2),size(p,3)), free_path(size(p,1), size(p,2),size(p,3)), C_c(size(p,1), size(p,2),size(p,3))
        real :: vsed(size(p,1), size(p,2),size(p,3)) !setttling velocity [m/s]
        real :: conc(size(p,1), size(p,2),size(p,3))
        real :: flux(size(p,1), size(p,2),size(p,3)) !mass/number flux, ug/m2/s or #/m2/s
        flux = 0
        do n=1,npop
        if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext") &
                .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "akk")) then !exclude ext and akk
                !code copied from Paul's sedimentation velocity in atmos_tracer_unitily
                !calculate sedimentation_velocity for the population
                viscosity = 1.458E-6 * T**1.5/(T+110.4)     ! Dynamic viscosity
                free_path = 6.6e-8*T/293.15*(PSTD_MKS/p)
                rwet = matrix_all_pop(n)%dg_wet(is:ie,js:je,:)/2
                rwet = max(1E-30, rwet) !avoid zero value issue
                rho_wet = matrix_all_pop(n)%dens_wet(is:ie,js:je,:)/2
                C_c = 1.0 + free_path/rwet * &              ! Slip correction [none]
                        (1.257+0.4*exp(-1.1*rwet/free_path))
                vsed = 2./9.*C_c*GRAV*rho_wet*rwet**2/viscosity  ! Settling velocity [m/s]
                !update matrix_value based on vdep
                do nt = 1, matrix_all_pop(n)%nb_tracer_pop
                        ntt = matrix_all_pop(n)%tracer_index(nt)
                        conc = matrix_all_tracer(ntt)%value_in_matrix(is:ie, js:je, :)
                        !flux at TOA
                        flux(:,:, 1) = -vsed(:,:,1) * conc(:,:,1) !unit: ug/m2/s
                        !flux in the middle to the bottom layer !set bottom layer as 0
                        flux(:, :, 2:kd-1) = vsed(:,:,1:kd-2) * conc(:,:,1:kd-2) - vsed(:,:,2:kd-1)*conc(:,:,2:kd-1)
                        flux(:, :, kd) = vsed(:,:,kd-1)*conc(:,:,kd-1)
                        matrix_all_tracer(ntt)%value_in_matrix(is:ie,js:je,:) = &
                               matrix_all_tracer(ntt)%value_in_matrix(is:ie,js:je,:) + flux/(zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)) * dt 
                enddo

        endif
        enddo


end subroutine



subroutine set_config_pop(matrix_configuration)
        character(len=*), intent(in):: matrix_configuration
        !----------------------------------------------------------------
        !       --- assign configuration population information ----
        ! akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
        ! if population exist, assign index >= 1, otherwise -1
        !----------------------------------------------------------------
        if (matrix_configuration .eq. "debug") then
            !debug version:
            !pop: ACC <- (N,M), EXT <- (M_H2O, M_H2SO4)
            !npop = 2 !number of poplation in the configuration
            !assign population index under current configuration
            !integer :: I_AKK = 1, I_ACC = 2, I_DD1 = 3, I_DD2 = 4 ! index of population in matrix
            !integer :: I_SSA = 5, I_SSC = 6, I_OC1 = 7, I_OC2 = 8 ! if exit, index >= 1, otherwise = -1
            !integer :: I_BC1 = 9, I_BC2 = 10, I_MXA = 11, I_MXC = 12, I_EXT = 13
            I_ACC = 2
            I_EXT = 13
            pop_active = 1 !neglect ext
            if (.not. allocated(I_POP)) &
                allocate(I_POP(pop_active))
            I_POP(:)=(/I_ACC/)
        elseif (matrix_configuration .eq. "debug_npf") then
            I_AKK = 1
            pop_active = 1
            if (.not. allocated(I_POP)) &
                allocate(I_POP(pop_active))
            I_POP(:)=(/I_AKK/)
        else
            call ERROR_MESG('get matrix configuration','configuration not defined '//trim(matrix_configuration), FATAL)
            !call error_mesg ('Tracer_driver', 'mw needs to be defined for tracer: '//trim(tracer_name), FATAL)
        end if
end subroutine



















      REAL(8) FUNCTION GETXNUM(NI,NJ,DGNI,DGNJ,XLSGI,XLSGJ)
!---------------------------------------------------------------------------------------------------------------------
! DLW, 102306: derived from function GETAF of CMAQ v4.4.
!
! GETXNUM = ln( Dij / Dgi ) / ( sqrt(2) * ln(Sgi) ), where
!
!      Dij is the diameter of intersection,
!      Dgi is the median diameter of the smaller size mode, and
!      Sgi is the geometric standard deviation of smaller mode.
!
! A quadratic equation is solved to obtain GETXNUM, following the method of Press et al. 1992.
!
! REFERENCES:
!
!  1. Binkowski, F.S. and S.J. Roselle, Models-3 Community Multiscale Air Quality (CMAQ)
!     model aerosol component 1: Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
!     doi:10.1029/2001JD001409, 2003.
!  2. Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in
!     Fortran 77 - 2nd Edition. Cambridge University Press, 1992.
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      REAL(8) :: NI         ! Aitken       mode number concentration [#/m^3]
      REAL(8) :: NJ         ! accumulation mode number concentration [#/m^3]
      REAL(8) :: DGNI       ! Aitken       mode geo. mean diameter [um]
      REAL(8) :: DGNJ       ! accumulation mode geo. mean diameter [um]
      REAL(8) :: XLSGI      ! Aitken       mode ln(geo. std. dev.) [1]
      REAL(8) :: XLSGJ      ! accumulation mode ln(geo. std. dev.) [1]

      ! Local variables.

      REAL(8) :: AA, BB, CC, DISC, QQ, ALFA, L, YJI
      REAL(8), PARAMETER :: SQRT2 = 1.414213562D+00

      ALFA = XLSGI / XLSGJ
      YJI = LOG( DGNJ / DGNI ) / ( SQRT2 * XLSGI )
      L = LOG( ALFA * NJ / NI)

      ! Calculate quadratic equation coefficients & discriminant.
            AA = 1.0D+00 - ALFA * ALFA
      BB = 2.0D+00 * YJI * ALFA * ALFA
      CC = L - YJI * YJI * ALFA * ALFA
      DISC = BB * BB - 4.0D+00 * AA * CC

      ! If roots are imaginary, return a negative GETAF value so that no IMTR takes place.

      IF( DISC .LT. 0.0D+00 ) THEN
        GETXNUM = - 5.0D+00         ! ERROR IN INTERSECTION
        RETURN
      ENDIF

      ! Equation 5.6.4 of Press et al. 1992.

      QQ = -0.5D+00 * ( BB + SIGN( 1.0D+00, BB ) * SQRT(DISC) )

      ! Return solution of the quadratic equation that corresponds to a
      ! diameter of intersection lying between the median diameters of the 2 modes.

      GETXNUM = CC / QQ       ! See Equation 5.6.5 of Press et al.

      ! WRITE(*,*)'GETXNUM = ', GETXNUM
      RETURN
      END FUNCTION GETXNUM
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


end module matrix_gfdl
