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
    use  constants_mod, only     : PI, GRAV, RDGAS, DENS_H2O, WTMAIR, AVOGNO
    use sat_vapor_pres_mod, only : compute_qs  
    use aero_npf, only : NPFRATE
    use aero_wet, only : aero_kohler
    use aero_condens, only : SETUP_KCI 
    implicit none
    private

    public matrix_init, matrix_source_type,set_matrix_source, matrix_run

    interface set_matrix_source
        module procedure set_matrix_source_2d
        module procedure set_matrix_source_3d
    end interface set_matrix_source

    real, parameter :: PI6 = PI/6.0
    integer, parameter :: nb_tracer_max = 9 !maxium number of tracers in a population, i.e. in MXX
    integer :: kd !kd = size(r,3), defined in matrix_init,3rd dimension of emission sources

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
        integer :: id_tracer_emis = -1 !diagostic IDs
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
        real, allocatable    :: vol_dry(:,:,:) !density of population: in Kg/m3
        real, allocatable    :: kappa_pop(:,:,:) !average hygroscopicity of the population
        integer :: id_Dg_dry = -1
        integer :: id_Dg_wet = -1
    end type matrix_pop

    type(matrix_tracer), allocatable :: matrix_all_tracer(:) ! define tracer array for all tracers
    type(matrix_pop),    allocatable :: matrix_all_pop(:)  ! defined population array for all population
    integer :: npop = 12 !maximum number of population in matrix
    !integer :: ntracer !maximum number of population and tracers in matrix
    integer :: I_AKK = -1, I_ACC = -1, I_DD1 = -1, I_DD2 = -1 ! index of population in matrix
    integer :: I_SSA = -1, I_SSC = -1, I_OC1 = -1, I_OC2 = -1 ! if exit, index >= 1, otherwise = -1
    integer :: I_BC1 = -1, I_BC2 = -1, I_MXX = -1, I_EXT = -1
    integer :: I_MW_H2SO4 = 1, I_MW_SO4 = 2, I_MW_DUST = 3, I_MW_SS = 4
    integer :: I_MW_OC = 5, I_MW_BC = 6, I_MW_SOA = 7
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
    logical :: do_matrix = .FALSE. 
    character(len=32) :: matrix_configuration
    character(len=7), parameter :: module_name = 'matrix'
    namelist /matrix_nml / do_matrix, matrix_configuration ! matrix_configuration is the version used in the calculation
    ! curretly matrix_configuration = debug is used to test toy model
contains

    !----------------------------------------------------------------
    !
    !                           subroutine matrix_run
    !IMPORTANT: this subroutine must be after atmos_SOx_chem etc, where sources are properly set-up 
    !the mass emission of different matrix_tracers are linked to different sources
    !in other F90 files, e.g. in atmos_SOx_chem of atmos_sulfate.F90 
    !----------------------------------------------------------------
    subroutine matrix_run(r, pfull, rh, t, dt, pwt, zhalf, rdt_matrix, Time, is,ie,js,je)
        real, intent(in) :: r(:,:,:,:)
        real, intent(in) :: dt !timestep
        integer, intent(in) :: is, ie, js, je ! boundaries of physical window
        type(time_type),  intent(in) :: Time
        real, intent(out) :: rdt_matrix(:,:,:,:) !tendency calculated from matrix
        real, intent(in) :: pfull(:,:,:), rh(:,:,:), t(:,:,:), pwt(:,:,:), zhalf(:,:,:) !t is temperature
        ! assign number rate to each number tracers, the mass emission has been linked to different tracers in other files
        real :: rt(size(r,1),size(r,2),size(r,3),size(r,4)) !local variable to record updated values of tracers, in gfdl unit
        !local variable to record absolute h2so4 mass change in ug/m3, positive
        real :: dm_h2so4(size(r,1),size(r,2),size(r,3)), dmdt_npf(size(r,1),size(r,2),size(r,3)), dndt_npf(size(r,1),size(r,2),size(r,3)) 
        integer :: n,MW, npf_flag_npf !local variables
        real :: m_akk_emis_source(size(r,1),size(r,2),size(r,3)), n_akk_emis_source(size(r,1),size(r,2),size(r,3))
        real :: kci_coef_pop(npop, size(r,1),size(r,2),size(r,3)), kci_coef_aeq1(npop, size(r,1),size(r,2),size(r,3)) !unit: m3/s
        logical :: used
        dm_h2so4 = 0.
        rdt_matrix = 0.            
        dndt_npf = 0.
        dmdt_npf = 0.
        
        if (do_matrix) then
            rt = r
            !------------------------------------------------------------------
            !               Step 0: get  matrix tracer values 
            !               and convert from gfdl to matrix unit
            !------------------------------------------------------------------
            !initialize matrix species values at the beginning of the current step
            call set_matrix_value(rt, pwt, zhalf) ! assign matrix_tarcer%values_in_matrix, value in matrix unit
            !-----------------------------------------------------------------------------
            !               Step 1: calculate tracer number sources from direct emission (in matrix unit)
            !                       Note: this step doesn't include new particle formation
            !    Note: the mass sources have been already updated before matrix_run
            !------------------------------------------------------------------------------
            call set_matrix_emis_number(pfull,pwt,zhalf) !update tracer source of  num_rate
            
            !------------------------------------------------------------------------------------
            !               Step *: H2SO4 condensational growth, must prior new particle formation
            !(1) calculate the condensation coefficient for all population over all grids
            !           kci_coef_pop 4-D dimensions: (npop, is, ij, ik)
            !           kci_coef_aeq1 4-D dimensions: (npop, is, ij, ik)
            !-------------------------------------------------------------------------------------
            call set_matrix_pop_kci(pfull,t,kci_coef_pop, kci_coef_aeq1) !calculate the condensation coefficient for all population

            !---------------------------------------------------------------------------------
            !               Step *: condensational growth w/wo new particle formation of AKK
            !               calculate new particle formation rate and update sources of dndt/dmdt 
            !XLXLXLXLIMPORTANT: H2SO4 need to be updated by the loss of H2SO4 condensational loss
            !-----------------------------------------------------------------------------------
            call set_matrix_condense_npf(I_AKK, pfull,rh,t, P_H2SO4_RATE, r(:,:,:,nh2so4),pwt,zhalf, &
                    dndt_npf, dmdt_npf, IH2SO4_PATH, XH2SO4_NUCL)

            if (I_AKK > 0) then
                    n_akk_emis_source = matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source
                    m_akk_emis_source = matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source !record ini source from emission
                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source = n_akk_emis_source+dndt_npf
                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source = m_akk_emis_source+dmdt_npf

            

            if (I_AKK > 0 ) then !if new particle formation if initiated, i.e. AKK exist
           

            else

            endif
          
          
          
          
          
!            if (I_AKK > 0 ) then !if new particle formation if initiated, i.e. AKK exist
!                 n_akk_emis_source = matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source
!                 m_akk_emis_source = matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source !record ini source from emission
!
!
!                 call matrix_npfrate(pfull,rh,t, P_H2SO4_RATE, r(:,:,:,nh2so4),pwt,zhalf, &
!                    dndt_npf, dmdt_npf)  !new particle formation number/mass rate in ug/m3
!                 !note: P_H2SO4_RATE already processes previously in ug/m3; H2SO4 unprocessed, in gfdl vmr unit
!
!
!
!
!
!
!                 matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source = n_akk_emis_source+dndt_npf
!                 matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source = m_akk_emis_source+dmdt_npf 
!                 !update H2SO4 loss
!                 dm_h2so4 = dm_h2so4 + dmdt_npf * dt !in matrix unit ug/m3
!            endif
!

            !-----------------------------------------------------------------------------
            !               Step *: Hygroscopic growth: calculate dry/wet particle diameter
            !
            !------------------------------------------------------------------------------
            call matrix_dry_diameter(rt) !assign dry pop properties 
            call matrix_wet_diameter(rt, rh, t) !assign wet pop properties, t is temperature




           !---------------------------------------------------------------------------------
            !               Final step: update matrix_tracer values and H2SO4 concentration
            !--------------------------------------------------------------------------------- 
            do n=1,ntrace
                if (matrix_all_tracer(n)%is_active) then
                    matrix_all_tracer(n)%value_in_matrix = matrix_all_tracer(n)%value_in_matrix + matrix_all_tracer(n)%source * dt
                endif
            enddo

            do n=1,ntrace
                if (matrix_all_tracer(n)%is_active) then
                    call update_rt_from_matrix(n, rt(:,:,:,n), pwt, zhalf)
                    rdt_matrix(:,:,:,n) = (rt(:,:,:,n)-r(:,:,:,n))/dt
                elseif (n .eq. nh2so4) then
                    !do unit from matrix unit to gfdl unit: ug/m3 -> vmr
                    MW = 98
                    dm_h2so4 = dm_h2so4/dt/1E9/(pwt(:,:,:) * MW / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
                    !dm_h2so4 can't exceed previous h2so4 concentration
                    dm_h2so4 = max(dm_h2so4, r(:,:,:,nh2so4))
                    rdt_matrix(:,:,:,n) = -dm_h2so4/dt
                endif
            enddo
         end if
         !--------------------------------------------------------------------------------------------------------------------
         ! Need to reset AKK population source to 0. as it is not done automatically if there is no other model source but NPF
         ! (set_matrix_source is not called)
         !--------------------------------------------------------------------------------------------------------------------         
        if (I_AKK > 0 ) then
                matrix_all_tracer(matrix_all_pop(I_AKK)%I_N)%source     = 0.
                matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source = 0.
        endif



        !send data of D_wet, D_dry for test: ps. Dg_dry, Dg_wet in matrix_pop unit is m, only for send_data change to um
        do n = 1, npop
                if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
                        if (matrix_all_pop(n)%id_Dg_dry > 0) then
                                used = send_data (matrix_all_pop(n)%id_Dg_dry, matrix_all_pop(n)%Dg_dry*1E6, time, &
                                        is_in=is,js_in=js, ks_in = 1) !in unit µm
                        endif
                        if (matrix_all_pop(n)%id_Dg_wet > 0) then
                                used = send_data (matrix_all_pop(n)%id_Dg_wet, matrix_all_pop(n)%Dg_wet*1E6, time, &
                                        is_in=is,js_in=js, ks_in = 1) !in unit µm
                        endif

                endif
         end do
        
end subroutine matrix_run


            






              

!            !------------------------------------------------------------------
!            !               Step 1: emission and nucleations
!            !------------------------------------------------------------------
!            !Step 1: emission and nucleations -> get number and mass rate
!            call set_matrix_emis_number(rt(:,:,:,nh2so4),pfull,rh,t,pwt,zhalf) !assign matrix_tracer source
!            !update matrix values from emission and nucleation
!            do n=1,ntrace
!            if (matrix_all_tracer(n)%is_active) then
!                matrix_all_tracer(n)%value_in_matrix = matrix_all_tracer(n)%value_in_matrix + matrix_all_tracer(n)%source * dt 
!            endif
!            enddo
!            !IF akk exist, i.e. npf is activated, then record h2so4 that got consumed
!            if (I_AKK > 0) then
!                dm_h2so4 = dm_h2so4 + matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source * dt 
!            endif
!            
!            
!
!            !x5l: hygroscopic growth
!
!            !matrix condensational growth
!
!            !matrix coagulation
!
!            !matrix intermodal transfer
!
!            !matrix activation (diagnostics for Huan)
!
!            !partition aqso4 mass
!
!            !overwrite rt with matrix_value with unit conversion
!
!
!
!            do n=1,ntrace
!            if (matrix_all_tracer(n)%is_active) then
!                call update_rt_from_matrix(n, rt(:,:,:,n), pwt, zhalf)
!                ! if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
!                !  write(*,*) "***********"
!                ! write(*,*) "n in ntrace", n
!                ! write(*,*) "shape of r(:,:,:,n)", shape(r(:,:,:,n))
!                ! write(*,*) "shape of rt(:,:,:,n)", shape(rt(:,:,:,n))
!                ! write(*,*) "shape of rdt_matrix(:,:,:,n)", shape(rdt_matrix(:,:,:,n))
!                ! endif 
!                rdt_matrix(:,:,:,n) = (rt(:,:,:,n)-r(:,:,:,n))/dt
!            elseif (n .eq. nh2so4) then
!                !do unit from matrix unit to gfdl unit: ug/m3 -> vmr
!                MW = 98
!                dm_h2so4 = dm_h2so4/dt/1E9/(pwt(:,:,:) * MW / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
!                !dm_h2so4 can't exceed previous h2so4 concentration
!                dm_h2so4 = max(dm_h2so4, r(:,:,:,nh2so4))
!                rdt_matrix(:,:,:,n) = -dm_h2so4/dt
!            else
!                rdt_matrix(:,:,:,n) = 0.
!            endif
!            enddo
!
!        end if
!
!
!    end subroutine matrix_run

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

        !----------------------------------------------------------------
        !	--- assign configuration population information ----
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
            !integer :: I_BC1 = 9, I_BC2 = 10, I_MXX = 11, I_EXT = 12
            I_ACC = 2
            I_EXT = 12
        elseif (matrix_configuration .eq. "debug_npf") then
            I_AKK = 1
        else
            call ERROR_MESG('get matrix configuration','configuration not defined '//trim(matrix_configuration), FATAL)
            !call error_mesg ('Tracer_driver', 'mw needs to be defined for tracer: '//trim(tracer_name), FATAL)
        end if

        call get_number_tracers(MODEL_ATMOS, num_tracers = ntrace)

        allocate(matrix_all_tracer(ntrace))
        allocate(matrix_all_pop(npop))
        if (ntrace > 0) then
            do n = 1, ntrace
            flag = query_method ('matrix_parameter',MODEL_ATMOS,n, &
                text_in_scheme,control)
            if (flag) then
                !IMPORTANT: get_matrix_tracer_param only 
                !(1) update has_emission for mass species
                !(2) allocate source array for mass species with lognormal, dens etc parameters in field table
                !(3) IMPORTANT: number tracer: has_emission & source_array allocation need to be set somewhere
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
        matrix_all_pop(n)%Dg_dry = 0.
        matrix_all_pop(n)%Dg_wet = 0.
        matrix_all_pop(n)%mass_dry = 0.
        matrix_all_pop(n)%vol_dry = 0.
        matrix_all_pop(n)%kappa_pop = 0.
    end do

    allocate(P_H2SO4_RATE(size(r,1),size(r,2),size(r,3))) !aqueous phase SO4 production
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
  

    kd = size(r,3)

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

        if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
                 write(*,*) "id_register test"
                 write(*,*) "pop name = ", trim(matrix_all_pop(n)%name)
                 write(*,*) "diag 1st wet string = ", trim(matrix_all_pop(n)%name)//'_dg_wet'
                 write(*,*) "diag 2nd wet string = ", trim(matrix_all_pop(n)%name)//'_dg_wet'
                 write(*,*) "dry id =", matrix_all_pop(n)%id_Dg_dry
                 write(*,*) "wet id =", matrix_all_pop(n)%id_Dg_wet
        endif
      endif
    end do
    
    !XL debug: test spec, dens
    do n = 1, npop
      if (matrix_all_pop(n)%nb_tracer_pop > 0 ) then
        if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
                 write(*,*) "pop index", n
                 write(*,*) "nb_tracer_pop ", matrix_all_pop(n)%nb_tracer_pop
                 nt = matrix_all_pop(n)%nb_tracer_pop
                 do ntt = 1, nt
                 write(*,*) "tracer_type ", matrix_all_tracer(matrix_all_pop(n)%tracer_index(ntt))%type
                 write(*,*) "tracer_spec ", matrix_all_tracer(matrix_all_pop(n)%tracer_index(ntt))%spec
                 write(*,*) "tracer_dens ", matrix_all_tracer(matrix_all_pop(n)%tracer_index(ntt))%dens
                 end do
        endif
      endif
    end do
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
    ! trim(population) name: akk/acc/dd1/dd2/ssa/ssc/oc1/oc2/bc1/bc2/mxx/ext
    if (lowercase(trim(population))=='akk') then
        MT%pop = 'akk'
        pop_index = I_AKK
    elseif (lowercase(trim(population))=='acc') then
        MT%pop = 'acc'
        pop_index = I_ACC
    elseif (lowercase(trim(population))=='dd1') then
        MT%pop = 'dd1'
        pop_index = I_DD1
    elseif (lowercase(trim(population))=='dd2') then
        MT%pop = 'dd2'
        pop_index = I_DD2
    elseif (lowercase(trim(population))=='ssa') then
        MT%pop = 'ssa'
        pop_index = I_SSA
    elseif (lowercase(trim(population))=='ssc') then
        MT%pop = 'ssc'
        pop_index = I_SSC
    elseif (lowercase(trim(population))=='oc1') then
        MT%pop = 'oc1'
        pop_index = I_OC1
    elseif (lowercase(trim(population))=='oc2') then
        MT%pop = 'oc2'
        pop_index = I_OC2
    elseif (lowercase(trim(population))=='bc1') then
        MT%pop = 'bc1'
        pop_index = I_BC1
    elseif (lowercase(trim(population))=='bc2') then
        MT%pop = 'bc2'
        pop_index = I_BC2
    elseif (lowercase(trim(population))=='mxx') then
        MT%pop = 'mxx'
        pop_index = I_MXX
    elseif (lowercase(trim(population))=='ext') then
        MT%pop = 'ext'
        pop_index = I_EXT
    else
        call ERROR_MESG('get_matrix_tracer_param', 'trim(population) not found '//trim(trim(population)), FATAL )        
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
    endif

    if (pop_index > 0) then !if the current configuration have this population
        MP(pop_index)%nb_tracer_pop = MP(pop_index)%nb_tracer_pop + 1   ! count the number of tracers in the specific populations 
        MP(pop_index)%tracer_index(MP(pop_index)%nb_tracer_pop) = tracer_index
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
subroutine set_matrix_value(rt, pwt, zhalf)
    real, intent(in) :: rt(:,:,:,:)
    real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
    integer :: n !local variables
    do n=1, ntrace !loop all tracers in rt array
    if (matrix_all_tracer(n)%is_active) then
        call set_unit_gfdl_to_matrix(rt(:,:,:,n),matrix_all_tracer(n)%units, &
            matrix_all_tracer(n)%type, matrix_all_tracer(n)%spec, matrix_all_tracer(n)%value_in_matrix, pwt, zhalf)
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
    real, intent(out) :: value_in_matrix(:,:,:) !value in matrix unit
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
subroutine update_rt_from_matrix(n_rt, mono_rt_tracer,pwt, zhalf)
    integer, intent(in) :: n_rt !index in rt(:,:,:,n_rt)
    real, intent(in), optional :: pwt(:,:,:), zhalf(:,:,:)
    real, intent(out) :: mono_rt_tracer(:,:,:) !rt(:,:,:,n_rt) to be updated
    real :: MW !molecular weight, local variable
    !number or mass tracers
    if (lowercase(matrix_all_tracer(n_rt)%type) .eq. "mass") then
        if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "mmr") then
            mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix/1E9/(pwt(:,:,:)/ (zhalf(:,:,1:kd) - zhalf(:,:,2:kd+1)))
        elseif (lowercase(matrix_all_tracer(n_rt)%units) .eq. "vmr") then
            if (lowercase(trim(matrix_all_tracer(n_rt)%spec)) .eq. "sulf") then
                MW = 96
                mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix/1E9/(pwt(:,:,:) * MW / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
            else
                call error_mesg ('matrix_gfdl','matrix tracer: molecular weight not properly defined', FATAL)
            endif
        endif
    elseif (lowercase(matrix_all_tracer(n_rt)%type) .eq. "number") then
        if (lowercase(matrix_all_tracer(n_rt)%units) .eq. "#/kg") then
            mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix / (pwt(:,:,:) / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
        elseif (lowercase(matrix_all_tracer(n_rt)%units) .eq. "vmr") then
            mono_rt_tracer(:,:,:) = matrix_all_tracer(n_rt)%value_in_matrix / (pwt(:,:,:)  / WTMAIR / (zhalf(:, :, 1:kd) - zhalf(:, :, 2:kd+1)))
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
        call set_matrix_source_generic(source_type,source_processed, time, is, js)  
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
        call set_matrix_source_generic(source_type,source_processed, time, is,js)
    end if

end subroutine set_matrix_source_3d

subroutine set_matrix_source_generic(source_type,source_processed, time, is,js)

    real, intent(in)    :: source_processed(:,:,:)
    integer, intent(in) :: source_type
    type(time_type), intent(in) :: time ! current model time
    integer, intent(in) :: is, js ! boundaries of physical window
    logical :: used 
    if (do_matrix) then
        !AKK mode H2SO4(g): vmr/s, 3D array
        if (source_type .eq. matrix_source_type%P_H2SO4) then   
            P_H2SO4_RATE = source_processed     
            !ACC mode SO4 emission: vmr/s, 3D array
        elseif (source_type .eq. matrix_source_type%E_SO4) then
            if (I_AKK > 0 .AND. I_ACC < 0) then
                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source =0.01* source_processed
            elseif (I_AKK > 0 .AND. I_ACC > 0) then
                    matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%source =0.99* source_processed
                    matrix_all_tracer(matrix_all_pop(I_AKK)%I_MSULF)%source =0.01* source_processed
            else
                     matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%source =source_processed
            endif
        elseif (source_type .eq. matrix_source_type%P_AQSO4) then
            AQSO4 =  source_processed
            !Dust emission: Kg/m2/s, 2D array
        elseif (source_type .eq. matrix_source_type%E_DUST) then
            if (I_DD1>0 .AND. I_DD2<0) then !only have 1 dust mode DD1
                matrix_all_tracer(matrix_all_pop(I_DD1)%I_MDUST)%source = source_processed
            elseif (I_DD1<0 .AND. I_DD2>0) then !only have 1 dust mode DD2
                matrix_all_tracer(matrix_all_pop(I_DD2)%I_MDUST)%source = source_processed
            elseif (I_DD1>0 .AND. I_DD2>0) then !have 2 dust mode: DD1 and DD2
                matrix_all_tracer(matrix_all_pop(I_DD1)%I_MDUST)%source = 0.25*source_processed
                matrix_all_tracer(matrix_all_pop(I_DD2)%I_MDUST)%source = 0.75*source_processed
            endif
            !Sea salt emission: Kg/m2/s, 2D array
        elseif (source_type .eq. matrix_source_type%E_SS) then
            if (I_SSA>0 .AND. I_SSC<0) then !only have 1 sea salt mode SSA
                matrix_all_tracer(matrix_all_pop(I_SSA)%I_MSEAS)%source = source_processed
            elseif (I_SSA<0 .AND. I_SSC>0) then !only have 1 sea salt mode SSC
                matrix_all_tracer(matrix_all_pop(I_SSC)%I_MSEAS)%source = source_processed
            elseif (I_SSA>0 .AND. I_SSC>0) then !have 2 sea salt mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_SSA)%I_MSEAS)%source = 0.25*source_processed
                matrix_all_tracer(matrix_all_pop(I_SSC)%I_MSEAS)%source = 0.75*source_processed
            endif
            !SOA source from GFDL model, add it to OC: kg/m2/s, 2D array
        elseif (source_type .eq. matrix_source_type%E_SOA) then
            if (I_OC1>0 .AND. I_OC2<0) then
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source + source_processed
            elseif (I_OC1<0 .AND. I_OC2>0) then
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source = &
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source + source_processed
            elseif (I_OC1>0 .AND. I_OC2>0) then !have 2 organic carbon mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source + 0.2*source_processed
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source = & 
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source + 0.8*source_processed
            endif
            !Organic carbon: vmr, 3D array
        elseif (source_type .eq. matrix_source_type%E_OC) then
            if (I_OC1>0 .AND. I_OC2<0) then
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source + source_processed
            elseif (I_OC1<0 .AND. I_OC2>0) then
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source = &
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source + source_processed
            elseif (I_OC1>0 .AND. I_OC2>0) then !have 2 organic carbon mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source = &
                    & matrix_all_tracer(matrix_all_pop(I_OC1)%I_MOCAR)%source + 0.5*source_processed
                matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source = &
                    & matrix_all_tracer(matrix_all_pop(I_OC2)%I_MOCAR)%source + 0.5*source_processed
            endif
            !Black carbon: vmr, 3D array
        elseif (source_type .eq. matrix_source_type%E_BC) then
            if (I_BC1>0 .AND. I_BC2<0) then
                matrix_all_tracer(matrix_all_pop(I_BC1)%I_MBCAR)%source = source_processed
            elseif (I_BC1<0 .AND. I_BC2>0) then
                matrix_all_tracer(matrix_all_pop(I_BC2)%I_MBCAR)%source = source_processed
            elseif (I_BC1>0 .AND. I_BC2>0) then !have 2 organic carbon mode: SSA and SSC
                matrix_all_tracer(matrix_all_pop(I_BC1)%I_MBCAR)%source = 0.8*source_processed
                matrix_all_tracer(matrix_all_pop(I_BC2)%I_MBCAR)%source = 0.2*source_processed
            endif

        endif

    endif

end subroutine set_matrix_source_generic

!-----------------------------------------------------------------------
!                 set_matrix_emis_number
! 1. convert mass concentration into number concentration
! 2. assign number information to tracers with number type
!-----------------------------------------------------------------------
subroutine set_matrix_emis_number(pfull,pwt,zhalf)
    real, intent(in) :: pfull(:,:,:), pwt(:,:,:), zhalf(:,:,:)
    real :: NUM_RATE
    integer :: i,j,k
    do n=1,npop
    if (matrix_all_pop(n)%nb_tracer_pop > 0) then !if this population exist in current configuration
           !if (lowercase(trim(matrix_all_pop(n)%name)) =='akk') then !if h2so4 nucleation is initiated
           !call matrix_npfrate(pfull,rh,t, P_H2SO4_RATE, H2SO4,pwt,zhalf, &                    
           !     matrix_all_tracer(matrix_all_pop(n)%I_N)%source, &    !number rate, dndt in #/m3
           !     matrix_all_tracer(matrix_all_pop(n)%I_MSULF)%source)  !new particle formation mass rate in ug/m3
           ! !note: P_H2SO4_RATE already processes previously in ug/m3; H2SO4 unprocessed, in gfdl vmr unit
            !set up akk number and mass
            do nt = 1,matrix_all_pop(n)%nb_tracer_pop
            if (matrix_all_pop(n)%has_emission(nt)) then
                ntt = matrix_all_pop(n)%tracer_index(nt)
                if ((matrix_all_tracer(ntt)%sigma > 0) .AND. (matrix_all_tracer(ntt)%distribution_index == DIST_LOGNORMAL) &
                    .AND. (matrix_all_tracer(ntt)%Dgn > 0) .AND. (matrix_all_tracer(ntt)%dens > 0) ) then !SUGGEST TO MOVE THIS CHECK TO INITIALIZATION
                    do i=1,size(pfull,1)
                        do j=1,size(pfull,2)
                                do k=1,size(pfull,3)
                                        if (matrix_all_tracer(ntt)%source(i,j,k)>0) then
                                                num_rate = matrix_all_tracer(ntt)%source(i,j,k)/(PI6*matrix_all_tracer(ntt)%dens*matrix_all_tracer(ntt)%DP0**3)*1E-9!emission rate converted to number rate
                                        matrix_all_tracer(matrix_all_pop(n)%I_N)%source(i,j,k) = num_rate !update the number variable with this number rate
                                        end if
                                 end do
                        end do
                    end do
                else
                    call error_mesg ('matrix_gfdl','Emission species not properly defined', FATAL)
                endif
            endif
            end do
        endif
    end do
end subroutine set_matrix_emis_number


!-----------------------------------------------------------------------
!                 set_matrix_pop_kci
! calculate condensational coefficient for different populations 
!-----------------------------------------------------------------------
 subroutine set_matrix_pop_kci(pfull,t,kci_coef_pop, kci_coef_aeq1) !calculate the condensation coefficient for all population
     real, intent(in) :: pfull(:,:,:) ! pressure on layers, Pa
     real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degK
     real, intent(out) :: kci_coef_pop(:,:,:,:), kci_coef_aeq1(:,:,:,:) !KCI units: m^3/s
     integer :: n,is,js,ks,i,j,k
     kci_coef_pop = 0.
     kci_coef_aeq1 = 0.
     is = size(pfull,1)
     js = size(pfull,2)
     ks = size(pfull,3)
     do n=1, npop
         if ((matrix_all_pop(n)%nb_tracer_pop > 0) .and. (lowercase(trim(matrix_all_pop(n)%name)) .ne. "ext")) then
             do i=1,is
                do j=1,js
                    do k=1,js
                        call setup_kci(pfull(i,j,k), t(i,j,k), matrix_all_pop(n)%Dg_dry(i,j,k), matrix_all_pop(n)%sigma, &
                        kci_coef_pop(n,i,j,k), kci_coef_aeq1(n,i,j,k))

                    enddo
                 enddo
             enddo
         endif
     enddo
 
 end subroutine


!-----------------------------------------------------------------------
!                 set_matrix_npfrate
! 1. set-up for loop to call NPFRATE in matrix
! 2. calculate DNDT and DMDT in 3D array form
!-----------------------------------------------------------------------
subroutine matrix_npfrate(pfull,rh,t, P_H2SO4_RATE, H2SO4,pwt,zhalf,dndt, dmdt)
    real, intent(in) :: pfull(:,:,:) ! pressure on layers, Pa
    real, intent(in) :: rh(:,:,:) ! relative humidity
    real, intent(in) :: pwt(:,:,:), zhalf(:,:,:)
    real, intent(in) :: P_H2SO4_RATE(:,:,:) !H2SO4 production rate in ugSO4/m^3
    real, intent(in) :: t(:,:,:) ! temperature of atmosphere, degK
    real, intent(in) :: H2SO4(:,:,:) !gaseous H2SO4 in gfdl unit, vmr
    real, intent(out):: dndt(:,:,:), dmdt(:,:,:)
    integer :: is, js, ks, i,j,k ! number of layers in different direction
    real:: pres_mt, rh_mt, temp_mt, h2so4_mt, p_h2so4rate_mt !local variable for matrix calculation
    real:: MW, dndt_tst, dmdt_tst !if directly assign dndt(i,j,k) and dmdt(i,j,k) segment error will occur
    real:: dndt_rec(size(pfull,1),size(pfull,2),size(pfull,3)), dmdt_rec(size(pfull,1),size(pfull,2),size(pfull,3)) 
    MW = matrix_molecular_weight(I_MW_H2SO4)
    is = size(pfull,1)
    js = size(pfull,2)
    ks = size(pfull,3)
    do i=1,is
        do j=1,js
                do k=1,ks
                        pres_mt = pfull(i,j,k) !pressure in matrix unit: [Pa]
                        rh_mt = rh(i,j,k) !fractional relative humidity [1]
                         !XL: check with FP for RH calculation
                         temp_mt = t(i,j,k) !ambient temperature [K]
                         !sulfuric acid (as SO4) concentration [ugSO4/m^3]
                         h2so4_mt = 1E9 * h2so4(i,j,k) * pwt(i,j,k) * MW / WTMAIR / (zhalf(i,j,k) - zhalf(i,j,k+1)) 
                         p_h2so4rate_mt = P_H2SO4_RATE(i,j,k) !gas-phase H2SO4 (as SO4) production rate [ugSO4/m^3 s]
                         !dndt: [m^-3 s^-1], dmdt: [ugSO4 m^-3 s^-1]
                         !call NPFRATE(pres_mt,rh_mt,temp_mt,h2so4_mt,p_h2so4rate_mt,dndt(i,j,k),dmdt(i,j,k))
                         call NPFRATE(pres_mt,rh_mt,temp_mt,h2so4_mt,p_h2so4rate_mt,dndt_rec(i,j,k),dmdt_rec(i,j,k))
                enddo
        enddo
   enddo
    dndt=dndt_rec
    dmdt=dmdt_rec
end subroutine matrix_npfrate

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
subroutine matrix_dry_diameter(r)
       real, intent(in) :: r(:,:,:,:)
       integer :: n, nt, ntt, tr_index !local variables
       real :: m_dry_all(size(r,1),size(r,2),size(r,3)), vol_dry_all(size(r,1),size(r,2),size(r,3))
       real :: kappa_vol_all(size(r,1),size(r,2),size(r,3)), vol_spec(size(r,1),size(r,2),size(r,3))
       real :: pop_num(size(r,1),size(r,2),size(r,3)), dens_dry(size(r,1),size(r,2),size(r,3))
       real :: exp_fac
       integer :: flag,is,js,ks, i, j, k !flag for whether this population has dry mass/perform calculation
       integer :: step = 0 !count for debug: which step crash
       is = size(r,1)
       js = size(r,2)
       ks = size(r,3)
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
                      m_dry_all = m_dry_all + matrix_all_tracer(tr_index)%value_in_matrix !m_dry unit: ug/m3
                      !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
                      !          write(*,*) "tracer_dens_test"
                      !          write(*,*) "tracer name = ", trim(matrix_all_tracer(tr_index)%spec)
                      !          write(*,*) "tracer density = ", matrix_all_tracer(tr_index)%dens
                      !endif
                      !vol_spec unit: m3 aerosol_volume/ m3 air
                      vol_spec = matrix_all_tracer(tr_index)%value_in_matrix/matrix_all_tracer(tr_index)%dens*1E-9 
                      vol_dry_all = vol_dry_all + vol_spec
                      kappa_vol_all = kappa_vol_all + matrix_all_tracer(tr_index)%kappa*vol_spec !sum of kappa*vol
                      flag = 1 !record if this population has dry mass
                  endif
               enddo
               if (flag > 0) then
                     !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
                     !           write(*,*) "tracer_dry_test"
                     !           write(*,*) "steps= ", step
                     !           write(*,*) "pop = ", n
                     !           write(*,*) "min m_dry_all = ", minval(m_dry_all)
                     !           write(*,*) "max m_dry_all = ", maxval(m_dry_all)
                     !           write(*,*) "min vol_dry_all = ", minval(vol_dry_all)
                     !           write(*,*) "max vol_dry_all = ", maxval(vol_dry_all)
                     ! endif
                   matrix_all_pop(n)%mass_dry = m_dry_all !mass unit: ug/m3
                   matrix_all_pop(n)%vol_dry = vol_dry_all !volume unit: m3 aerosol/m3
                   pop_num = matrix_all_tracer(matrix_all_pop(n)%I_N)%value_in_matrix !#/m3
                   exp_fac = exp(1.5*(log(matrix_all_pop(n)%sigma))**2) 
                   do i=1,is
                        do j=1,js
                                do k=1,ks
                                        if ((pop_num(i,j,k) > 1E-14) .AND. (m_dry_all(i,j,k) > 1E-32)) then !not zero number on the grid
                                            dens_dry(i,j,k) = m_dry_all(i,j,k)/vol_dry_all(i,j,k)*1E-9
                                            matrix_all_pop(n)%kappa_pop(i,j,k) = kappa_vol_all(i,j,k)/vol_dry_all(i,j,k)
                                            matrix_all_pop(n)%Dg_dry(i,j,k) = &
                                            (m_dry_all(i,j,k)/pop_num(i,j,k)/dens_dry(i,j,k)*1E-9/PI6)**(1.0/3)/exp_fac ! Dg_dry unit: m
                                                if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG,output outliers
                                                        if (abs(matrix_all_pop(n)%Dg_dry(i,j,k)-6.8E-8)/6.8E-8 > 0.1) then
                                                                write(*,*) "outlier_info_check:"
                                                                write(*,*) "step=", step
                                                                write(*,*) "dens_dry_ijk= ", dens_dry(i,j,k)
                                                                write(*,*) "m_dry_all  = ", m_dry_all(i,j,k)
                                                                write(*,*) "pop_num = ", pop_num(i,j,k)
                                                                write(*,*) "exp_term ", exp_fac
                                                                write(*,*) "DG result ", matrix_all_pop(n)%Dg_dry(i,j,k)
                                                        endif
                                                endif
                                        else
                                            dens_dry(i,j,k) = 0
                                            matrix_all_pop(n)%kappa_pop(i,j,k)=0
                                            matrix_all_pop(n)%Dg_dry(i,j,k) =0
                                        endif
                                enddo
                        enddo
                   enddo
                      !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
                      !          write(*,*) "min pop_num = ", minval(pop_num)
                      !          write(*,*) "max pop_num = ", maxval(pop_num)
                      !          write(*,*) "min dens_dry_all = ", minval(dens_dry)
                      !          write(*,*) "max dens_dry_all = ", maxval(dens_dry)

                      !          write(*,*) "min D_dry_all = ", minval(matrix_all_pop(n)%Dg_dry)
                      !          write(*,*) "max D_dry_all = ", maxval(matrix_all_pop(n)%Dg_dry)
                      !endif
               endif
           endif 
       end do  
       step = step + 1
end subroutine
        


subroutine matrix_wet_diameter(r, rh, t) !calculate and assgin values for pop%Dg_wet
       real, intent(in) :: r(:,:,:,:)
       real, intent(in) :: rh(:,:,:), t(:,:,:)
       real:: ddry_in_3d(size(r,1),size(r,2),size(r,3)), hygro_3d(size(r,1),size(r,2),size(r,3))
       integer :: n, i, j, k, is, js, ks, tr_index
       real :: ddry_in, hygro_in, s_in, tair_in, dwet_out, gf, rh_deliquescence, rh_crystallization
       real :: particle_vol_dry, particle_vol_water, particle_vol_wet, f_hysteresis 
       integer :: flag = 0 ! check if the population get calculated 
       is = size(r,1)
       js = size(r,2)
       ks = size(r,3)
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
                ddry_in_3d = matrix_all_pop(n)%Dg_dry*exp(1.5*(log(matrix_all_pop(n)%sigma))**2) !volume mean diameter
                hygro_3d = matrix_all_pop(n)%kappa_pop
                rh_deliquescence = matrix_all_pop(n)%rh_deliquescence
                rh_crystallization = matrix_all_pop(n)%rh_crystallization
                do i=1,is
                        do j=1,js
                                do k=1,ks
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
                                        matrix_all_pop(n)%Dg_wet(i,j,k) = gf*matrix_all_pop(n)%Dg_dry(i,j,k)
                                        !if (mpp_root_pe().eq.mpp_pe()) then !XL DEBUG
                                        !        write(*,*) 'D_dry=', matrix_all_pop(n)%Dg_dry(i,j,k)
                                        !        write(*,*) 'D_wet=', matrix_all_pop(n)%Dg_wet(i,j,k)
                                        !endif
                                enddo
                        enddo
                enddo
          endif
       enddo
end subroutine
        !subroutine matrix_nucleation()

!end subroutine matrix_nucleation()

!subroutine rh_calc(pmid, temp, sh, rh)

!!  implicit none
!
!  real, intent(in), dimension(:) :: pmid, temp, sh
!  real, intent(out), dimension(:) :: rh
!
!  !-----------------------------------------------------------------------
!  !       Calculate RELATIVE humidity.
!  !       This is calculated according to the formula:
!  !
!  !       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!  !
!  !       Where epsilon = Rdgas/RVgas = d622
!  !
!  !       and where 1- epsilon = d378
!  !
!  !       Note that rh does not have its proper value
!  !       until all of the following code has been executed.  That
!  !       is, rh is used to store intermediary results
!  !       in forming the full solution.
!  !-----------------------------------------------------------------------
!
!  !-----------------------------------------------------------------------
!  !calculate water saturated specific humidity
!  !-----------------------------------------------------------------------
!!  call compute_qs (temp, pmid, rh, q = sh)
!
!  !-----------------------------------------------------------------------
!  !calculate rh
!  !-----------------------------------------------------------------------
!  rh(:)= sh(:) / rh(:)
!
!end subroutine rh_calc

!Note for XL
!        matrix_all_tracer(ntt)%id_tracer_emis = register_diag_field ( module_name,     &
!            trim(matrix_all_tracer(ntt)%name)//'_src', axes(1:3),Time,  &
!            trim(matrix_all_tracer(ntt)%name)//'_src', 'µg/m3/s',       &
!            missing_value=-999.  )


!            if (matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%id_tracer_emis > 0) then
!                used = send_data (matrix_all_tracer(matrix_all_pop(I_ACC)%I_MSULF)%id_tracer_emis, source_processed, time, &
!                        is_in=is,js_in=js, ks_in = 1)
!            endif



end module matrix_gfdl
