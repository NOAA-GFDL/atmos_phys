module atmos_fire_plumerise_mod
  ! <DESCRIPTION>
  !    This module provides subroutines to calculate fire injection height
  !    and the share of biomass burning emissions in each vertical layer
  !    
  ! </DESCRIPTION>

  use mpp_mod,           only: input_nml_file
  use fms_mod,           only: fms_init, &
       mpp_pe, mpp_root_pe, stdlog, &
       write_version_number, &
       check_nml_error, error_mesg, &
       FATAL
  ! <---h1g,
  use fms2_io_mod,                only : file_exists
  use            fms_mod, only : lowercase, uppercase, &
       NOTE
  use time_manager_mod,           only : time_type, &
                                       days_in_month, days_in_year, &
                                       set_date, set_time, get_date_julian, &
                                       print_date, get_date, &
                                       operator(>), operator(+), operator(-)
  use   diag_manager_mod, only : send_data, &
       register_diag_field, get_base_time
  use atmos_cmip_diag_mod,   only : register_cmip_diag_field_3d, &
       register_cmip_diag_field_2d, &
       send_cmip_data_3d, &
       cmip_diag_id_type, &
       query_cmip_diag_id
  use tracer_manager_mod, only : query_method, &
       get_tracer_names, &
       get_tracer_index, &
       get_number_tracers, &
       MAX_TRACER_FIELDS
  use      constants_mod, only : GRAV, &     ! acceleration due to gravity [m/s2]
       RDGAS, &    ! gas constant for dry air [J/kg/deg]
       vonkarm, &
       PI, &
       DENS_H2O, & ! Water density [kg/m3]
       WTMH2O, &   ! Water molecular weight [g/mole]
       WTMAIR, &   ! Air molecular weight [g/mole]
       AVOGNO, &   ! Avogadro's number
       PSTD_MKS, &
       KAPPA     
  use   interpolator_mod, only : interpolator,  &
       obtain_interpolator_time_slices, &
       unset_interpolator_time_flag, &
       interpolate_type, &
       interpolator_init, interpolator_end, & !f1p
       CONSTANT, & !f1p
       INTERP_WEIGHTED_P !f1p
  use      astronomy_mod, only : universal_time

  implicit none
  private
  !-----------------------------------------------------------------------
  !----- interfaces -------

  public  atmos_fire_plumerise_time_vary,    &
       atmos_fire_plumerise_end, &
       atmos_fire_plumerise_endts, &
       atmos_fire_plumerise_init, &
       atmos_fire_plumerise_driver, &
       atmos_fire_do_bb_emis_diurnal
  
integer :: id_injhgt,id_FRP,id_bvf2, id_fbb
integer :: id_pfull_pt, id_temp_pt, id_z_half_pt, id_z_pbl_pt, id_z_full_pt
integer :: id_fbb_FT, id_fbb_NFT
integer, dimension(6)   :: id_fbb_perc
integer, dimension(6)   :: id_frp_perc
logical :: module_is_initialized = .FALSE.
logical :: used
type(interpolate_type),save  ::frp_value_interp, fbb_value_interp
type(time_type), save :: frp_offset
type(time_type), save :: frp_entry, fbb_entry
logical, save    :: frp_negative_offset
integer, save    :: frp_time_serie_type
character(len=80) :: frp_filename  = ' '
character(len=80) :: fbb_filename  = ' '
integer :: i
character(len=80), save, dimension(99) :: frp_value_name     = (/(' ',i=1,99)/)
character(len=80), save, dimension(99) :: fbb_value_name     = (/(' ',i=1,99)/)
character(len=80), dimension(99) :: frp_input_name  = (/(' ',i=1,99)/)
real,              dimension(99) :: perc_share      = (/(0.0,i=1,99)/)
character(len=80)     :: frp_source = ' '
character(len=80)     :: fbb_source = ' '
character(len=80)     :: frp_time_dependency_type = 'constant'
integer, dimension(6) :: frp_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
logical               :: do_LM4_fire_emis = .false.
integer               :: num_percentiles = 6
logical               :: do_inj_FRP_dist = .false. !!! BB emission injection based on distribution of FRP in each grid
logical               :: do_frp_diurnal = .false.
logical               :: do_bb_emis_diurnal = .false.
logical               :: do_bb_plumerise = .false.
integer               :: nlevel_fire = 6
integer               :: nselect_perc = 6
type(time_type) :: model_init_time     
character(len=80)  :: fire_injhgt_scheme = "EVEN_PBL"
!!! Options are "EVEN_PBL", "SURFACE", "SOFIEV", "SOFIEV-2step", "VALMARTIN", "AEROCOM"

namelist /fire_plumerise_nml/ &
 do_LM4_fire_emis, fire_injhgt_scheme, num_percentiles, nselect_perc, &
 do_inj_FRP_dist, fbb_source, fbb_filename, frp_source, frp_input_name, frp_filename, &
 perc_share, frp_time_dependency_type, frp_dataset_entry, &
 do_frp_diurnal, do_bb_emis_diurnal, do_bb_plumerise

!---- version number -----             
character(len=128) :: version = '$Id$' 
character(len=128) :: tagname = '$Name$'

type(time_type)                        :: frp_time, fbb_time


contains

subroutine atmos_fire_plumerise_init(lonb, latb, axes, Time, do_bb_plumerise_out)

    ! Routine to initialize.
    ! This registers the diag fields
    real, dimension(:,:),  intent(in) :: lonb, latb
    type(time_type),       intent(in) :: Time

    integer        , intent(in)       :: axes(4)
    logical        , intent(out)      :: do_bb_plumerise_out

    integer :: ntrace
    character(len=20) :: units =''
    !
    integer :: n, logunit
    character(len=128) :: name

    logical  :: flag

    !   local variables:
    integer   :: io, ierr
   character(len=80) :: f_name, description
   character(len=7), parameter :: mod_name = 'tracers'

   if (module_is_initialized) return
!----------------------------------
!namelist files
      if ( file_exists('input.nml')) then
        read (input_nml_file, nml=fire_plumerise_nml, iostat=io)
        ierr = check_nml_error(io,'fire_plumerise_nml')
      endif
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=fire_plumerise_nml)

    if (mpp_root_pe().eq.mpp_pe()) then
       write(*,*) 'name: error'
    end if

 call write_version_number (version, tagname)


!----------------------------------------------------------------------
!    register diag field
!----------------------------------------------------------------------
     do n = 1, num_percentiles
        write(f_name, "('fracs', i1)") n
        write(description, "('Fractions (lf_', i1, ')')") n
              
        id_fbb_perc(n)  = register_diag_field ( mod_name,           &
                    f_name, axes(1:3),Time,                 &
                    description, 'No Unit' )
     end do
     do n = 1, num_percentiles
        write(f_name, "('frp', i1)") n
        write(description, "('FRP perc (lf_', i1, ')')") n

        id_frp_perc(n)  = register_diag_field ( mod_name,           &
                    f_name, axes(1:2),Time,                 &
                    description, 'MW' )
     end do
     id_fbb = register_diag_field ( mod_name,           &
                    'fire_fbb', axes(1:3), Time,          &
                    'Vertical share of injection of fire emissions', 'Unitless' )
     id_injhgt = register_diag_field ( mod_name,           &
                    'fire_injhgt', axes(1:2), Time,          &
                    'Plume top height for injection of fire emissions', 'm' )
     id_FRP = register_diag_field ( mod_name,             &
                    'fire_FRP', axes(1:2), Time,              &
                    'Fire radiative power used for plume height calculation', 'MW/m2' )
     id_bvf2 = register_diag_field ( mod_name,             &
                    'bvf2_zpblx2', axes(1:2), Time,              &
                    'Brunt-Vaisala frequency (squared) at a height of twice the PBL', '1/sec2' )
     id_pfull_pt = register_diag_field ( mod_name,           &
                    'pfull_pt', axes(1:3), Time,          &
                    'tmp', 'tmp' )
     id_temp_pt = register_diag_field ( mod_name,           &
                    'temp_pt', axes(1:3), Time,          &
                    'tmp', 'tmp' )
     id_z_half_pt = register_diag_field ( mod_name,           &
                    'z_half_pt', axes(1:3), Time,          &
                    'tmp', 'tmp' )
     id_z_pbl_pt = register_diag_field ( mod_name,           &
                    'z_pbl_pt', axes(1:2), Time,          &
                    'tmp', 'tmp' )
     id_z_full_pt = register_diag_field ( mod_name,           &
                    'z_full_pt', axes(1:3), Time,          &
                    'tmp', 'tmp' )
     id_fbb_FT = register_diag_field ( mod_name,           &
                    'fbb_FT', axes(1:2), Time,          &
                    'Free tropospheric share of injection of fire emissions', 'Unitless' )
     id_fbb_NFT = register_diag_field ( mod_name,           &
                    'fbb_NFT', axes(1:2), Time,          &
                    'PBL share of injection of fire emissions', 'Unitless' )
!----------------------------------------------------------------------
!    initialize namelist entries
!----------------------------------------------------------------------
        frp_offset  = set_time (0,0)
        frp_entry  = set_time (0,0)
        frp_negative_offset  = .false.
        frp_time_serie_type  = 1
!----------------------------------------------------------------------
!    define the model base time  (defined in diag_table)
!----------------------------------------------------------------------
        model_init_time = get_base_time()

   if ( trim(fbb_source) .ne. ' ') then

     select case (trim(fbb_source))

       case ('MISR_clim')
!  read climatology for bb injection (Val Martin 2018)
        fbb_entry  = set_date (2008, &
                                  1,1,0,0,0)
        call error_mesg ('atmos_fire_plumerise', &
           'fbb is defined from a single annual cycle &
                &- with seasonal variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'fbb correspond to year :', &
                   '2008'
        endif
        do n = 1,25
          write(fbb_value_name(n), '("perc_inj_", I2.2)') n
        end do        

        call interpolator_init (fbb_value_interp,             &
                             trim(fbb_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = fbb_value_name(1:25),        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
       case ('AEROCOM_clim')
        fbb_entry  = set_date (2000, &
                                  1,1,0,0,0)
        call error_mesg ('atmos_fire_plumerise', &
           'fbb is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'fbb correspond to year :', &
                   '2000'
        endif
        do n = 1,7
          write(fbb_value_name(n), '("perc_inj_", I2.2)') n
        end do

        call interpolator_init (fbb_value_interp,             &
                             trim(fbb_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = fbb_value_name(1:7),        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )

     end select
   endif
!---------------------------------------------------------------------

   if ( trim(frp_source) .ne. ' ') then
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(frp_time_dependency_type) == 'constant' ) then
        frp_time_serie_type = 1
        frp_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'frp are constant in module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for frp is selected.
!---------------------------------------------------------------------
      else if (trim(frp_time_dependency_type) == 'time_varying') then
        frp_time_serie_type = 3
        if (frp_dataset_entry(1) == 1 .and. &
            frp_dataset_entry(2) == 1 .and. &
            frp_dataset_entry(3) == 1 .and. &
            frp_dataset_entry(4) == 0 .and. &
            frp_dataset_entry(5) == 0 .and. &
            frp_dataset_entry(6) == 0 ) then
          frp_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to frp_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
           frp_entry  = set_date (frp_dataset_entry(1), &
                                  frp_dataset_entry(2), &
                                  frp_dataset_entry(3), &
                                  frp_dataset_entry(4), &
                                  frp_dataset_entry(5), &
                                  frp_dataset_entry(6))
        endif
        call print_date (frp_entry , str= &
          'Data from frp timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        frp_offset = frp_entry - model_init_time
        if (model_init_time > frp_entry) then
          frp_negative_offset = .true.
        else
          frp_negative_offset = .false.
        endif
      else if (trim(frp_time_dependency_type) == 'fixed_year') then
        frp_time_serie_type = 2
        if (frp_dataset_entry(1) == 1 .and. &
            frp_dataset_entry(2) == 1 .and. &
            frp_dataset_entry(3) == 1 .and. &
            frp_dataset_entry(4) == 0 .and. &
            frp_dataset_entry(5) == 0 .and. &
            frp_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_fire_plumerise', &
            'must set frp_dataset_entry when using fixed_year source', FATAL)
        endif
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to frp_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        frp_entry  = set_date (frp_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_fire_plumerise', &
           'frp is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'frp correspond to year :', &
                    frp_dataset_entry(1)
        endif
     endif
     select case (trim(frp_source))
       case ('MODIS') 
         if (trim(frp_input_name(1)) .eq. ' ') then
           frp_value_name(1)='per10'
           frp_value_name(2)='per25'
           frp_value_name(3)='per50'
           frp_value_name(4)='per75'
           frp_value_name(5)='per90'
           frp_value_name(6)='per99'
           perc_share(1) = 0.10
           perc_share(2) = 0.15
           perc_share(3) = 0.25
           perc_share(4) = 0.25
           perc_share(5) = 0.15
           perc_share(6) = 0.10
         else
           do n=1, num_percentiles
           frp_value_name(n)(:)=trim(frp_input_name(n)(:))
           end do
         endif
     end select
     call interpolator_init (frp_value_interp,             &
                             trim(frp_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = frp_value_name(1:num_percentiles),        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
   endif

   do_bb_plumerise_out = do_bb_plumerise

 module_is_initialized = .TRUE.

end subroutine atmos_fire_plumerise_init


!####################################################################

subroutine atmos_fire_plumerise_time_vary (model_time)


 type(time_type), intent(in) :: model_time
 integer ::  yr, dum, mo_yr, mo, dy, hr, mn, sc, dayspmn

!  read climatology for bb injection (Val Martin 2018)
   if ( trim(fbb_source) .ne. ' ') then
     select case (trim(fbb_source))
       case ('MISR_clim')
         call get_date (fbb_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         fbb_time = set_date (yr, mo, 1, 0, 0, 0)
         if (mo == 12 .and. dy == 31) then
           fbb_time = set_date(yr,1,1,0,0,0)
         endif
         call obtain_interpolator_time_slices (fbb_value_interp, fbb_time)
       case ('AEROCOM_clim')
         call get_date (fbb_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         fbb_time = set_date (yr, 1, 1, 0, 0, 0)
         call obtain_interpolator_time_slices   &
                       (fbb_value_interp, fbb_time)
     end select
   endif
!---------------------------------------------------------------------
 

   if ( trim(frp_source).ne.' ') then

!--------------------------------------------------------------------
!    define the time in the frp data set from which data is to be 
!    taken. if frp is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(frp_time_serie_type .eq. 3) then
       if (frp_negative_offset) then
         frp_time = model_time - frp_offset
       else
         frp_time = model_time + frp_offset
       endif
       !!!! end of the year exception
       call get_date (frp_entry, yr, dum,dum,dum,dum,dum)
       call get_date (frp_time, mo_yr, mo, dy, hr, mn, sc)

       !!!! set all day frp values to frp max
       frp_time = set_date(mo_yr, mo, dy,0,0,0)

       if (mo == 12 .and. dy == 31) then
         frp_time = set_date(yr,12,30,0,0,0)
       endif
     else
       if(frp_time_serie_type .eq. 2 ) then
         call get_date (frp_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(frp_entry)
           if (dayspmn /= 29) then
             frp_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             frp_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           frp_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         frp_time = model_time
       endif
     endif

     call obtain_interpolator_time_slices   &
                       (frp_value_interp, frp_time)
   endif



end subroutine atmos_fire_plumerise_time_vary
!#######################################################################
!</SUBROUTINE>
!#######################################################################

subroutine atmos_fire_plumerise_endts

   if ( trim(fbb_source).ne.' ') then
     call unset_interpolator_time_flag (fbb_value_interp)
   endif

   if ( trim(frp_source).ne.' ') then
     call unset_interpolator_time_flag (frp_value_interp)
   endif

end subroutine atmos_fire_plumerise_endts

!#######################################################################

subroutine atmos_fire_plumerise_end

     call interpolator_end ( fbb_value_interp)
     call interpolator_end ( frp_value_interp)

     module_is_initialized = .FALSE.

end subroutine atmos_fire_plumerise_end

!#######################################################################
!
!<SUBROUTINE NAME="fire_fbb"> 
subroutine fire_fbb(fbbl,z_plume,bvf2,f_scheme,kd,id,jd,pfull_pt,temp_pt, &
                      z_half_pt,z_pbl_pt,z_full_pt,FRP,nlevel_fire, &
                      alt_fire_max,alt_fire_min, fbb_clim)
!
!<OVERVIEW>
! A routine to calculate the vertical profile of biomass burning emissions
! using one of several schemes chosen by the user with namelist option: 
! 'fire_injhgt_scheme'.
!</OVERVIEW>
!<DESCRIPTION>
! Routine to calculate the vertical profile of biomass burning emissions
! given environmental conditions (PBL height, Brunt-Vaisala freq) and the 
! radiative power of the fire itself (optional).  The fraction of emissions
! released at each model vertical level is returned.
!</DESCRIPTION>
! 
!<TEMPLATE>
! call fire_fbb()
!</TEMPLATE>
! INTENT IN
!<IN NAME="f_scheme" TYPE="character" DIM="(1)">
!  Injection height scheme option
!</IN>
!<IN NAME="pfull_pt" TYPE="real" DIM="(:)">
!  pressure at full model levels
!</IN>
!<IN NAME="temp_pt" TYPE="real" DIM="(:)">
!  temperature at full model levels
!</IN>
!<IN NAME="z_half_pt" TYPE="real" DIM="(:)">
!  height in meters at level interfaces (half-levels)
!</IN>
!<IN NAME="z_pbl_pt" TYPE="real" DIM="(1)">
!  Depth of the planetary boundary layer [m]
!</IN>
!<IN NAME="z_full_pt" TYPE="real" DIM="(:)">
!  height in meters at full levels 
!</IN>
!<IN NAME="kd" TYPE="real" DIM="(1)">
!  Number of vertical levels
!</IN>
!<IN NAME="nlevel_fire" TYPE="integer" DIM="(1)" (optional)>
!  Number of levels to inject fire emissions for fixed height schemes (e.g. AEROCOM)
!</IN>
!<IN NAME="alt_fire_min" TYPE="real" DIM="(6)" (optional)>
!  Height [m] of lower boundary of layer in which fire emissions are input
!  for fixed height schemes
!</IN>
!<IN NAME="alt_fire_max" TYPE="real" DIM="(6)" (optional)>
!  Height [m] of upper boundary of layer in which fire emissions are input 
!  for fixed height schemes
!</IN>
!<IN NAME="FRP" TYPE="real" DIM="(1)" (optional)>
!  Fire radiative power, used to compute injection height in the Sofiev scheme
!</IN>

! 
! INTENT INOUT
!<INOUT NAME="fbbl" TYPE="real" DIM="(:,:,:)">
!  Fraction of emissions injected into each layer,
!  second dimension is fixed height scheme layers
!</INOUT>
!<INOUT NAME="z_plume" TYPE="real" DIM="()">
!  Height of plume top calculated by the Sofiev scheme
!</INOUT>
!<INOUT NAME="bvf2" TYPE="real" DIM="()">
!  Brunt-Vaisala frequency squared at twice the pbl height
!</INOUT>

   real, intent(in), dimension(:,:,:)    :: pfull_pt
   real, intent(in), dimension(:,:,:)    :: temp_pt
   real, intent(in), dimension(:,:,:)    :: z_half_pt
   real, intent(in), dimension(:,:,:)    :: z_full_pt
   real, intent(in), dimension(:,:)      :: z_pbl_pt
   integer, intent(in)                  :: kd
   integer, intent(in) :: id
   integer, intent(in) :: jd
   real, intent(in), dimension(:), optional    :: alt_fire_min
   real, intent(in), dimension(:), optional    :: alt_fire_max
   real, intent(in), dimension(:,:), optional        :: FRP
   real, intent(in), dimension(:,:,:), optional        :: fbb_clim
   integer, intent(in), optional     :: nlevel_fire
   character(len=80), intent(in)     :: f_scheme 
   real, intent(inout), dimension(:,:,:) :: fbbl
   real, intent(inout), dimension(:,:)   :: z_plume
   real, intent(inout), dimension(:,:)   :: bvf2

!!! Local variables for Sofiev plume height parameterization (two-step)
   real, parameter :: sof_a_00 = 0.15, sof_a_1 = 0.93, sof_a_2 = 0.24  !!! alpha parameter 
   real, parameter :: sof_b_00 = 102., sof_b_1 = 298., sof_b_2 = 170.  !!! beta parameter
   real, parameter :: sof_c_00 = 0.49, sof_c_1 = 0.13, sof_c_2 = 0.35  !!! gamma parameter
   real, parameter :: sof_d_00 = 0.,  sof_d_1 = 0.7,   sof_d_2 = 0.6   !!! delta parameter
   real, parameter :: sof_N00 = 0.00025  !!! Reference squared Brunt-Vaisala freq [s-2] 
   real, parameter :: sof_P00 = 1.e6   !!! Reference FRP  [W]
!!! Local variables for Sofiev plume height parameterization (one-step)
!!! MODIS Aqua
   real, parameter :: sof_a_0 = 0.24   !!! alpha parameter, for VIIRS sof_a_0 = 0.98 
   real, parameter :: sof_b_0 = 170.  !!! beta parameter, for VIIRS sof_b_0 = 100.
   real, parameter :: sof_c_0 = 0.35  !!! gamma parameter, for VIIRS sof_c_0 = 0.30
   real, parameter :: sof_d_0 = 0.6   !!! delta parameter, for VIIRS sof_d_0 = 0.50
   real, parameter :: sof_N0 = 0.00025!!! Reference squared Brunt-Vaisala freq [s-2], for VIIRS sof_N0 = 0.00025 
   real, parameter :: sof_P0 = 1.e6   !!! Reference FRP  [W], for VIIRS sof_P0 = 7.e3

   logical         :: used  !!! for diagnostics
   integer         :: l, lf, ktop, kbot, logunit, i, j
   real            :: Z0, Z1, Z2, del, fbb_tot
!   real            :: bvf2=0.   !!! SQUARED Brunt-Vaisala freq [s-2] at twice the height of the PBL

   real :: FRPi=2.e2

!!!       changed Z1 and Z2 to be the height above the ground level,
!!!       not the height above the lowest level since this should be
!!!       more consistent with the intent of the AEROCOM levels and
!!!       with the other schemes.

   
   do j = 1, jd
   do i = 1, id

     if (f_scheme.eq."SOFIEV".or.f_scheme.eq."SOFIEV-2STEP".or.f_scheme.eq."SOFIEV-SFRP") FRPi = FRP(i,j) * 1.e6  !!! MW to W

     z_plume(i,j) = 0.0
     bvf2(i,j) = 0.0
!!! First compute plume height for SOFIEV scheme before distributing emissions
!   if (f_scheme.eq."SOFIEV".and.FRPi.gt.0.) then
     if (f_scheme.eq."SOFIEV".or.f_scheme.eq."SOFIEV-CF".and.FRPi.gt.0.) then

       !! Step 0: Compute BVfreq for twice the height of the PBL
       do l = kd,2,-1
           if (z_pbl_pt(i,j)*2. .lt. z_full_pt(i,j,kd)) z_plume(i,j)=0.0
           if (z_pbl_pt(i,j)*2. .ge. z_full_pt(i,j,l) .and. z_pbl_pt(i,j)*2. .lt. z_half_pt(i,j,l-1)) then
               ktop=l-1
               kbot=l
               bvf2(i,j)=brunt_vaisala(pfull_pt(i,j,:),temp_pt(i,j,:),ktop,kbot)
               if (bvf2(i,j).lt.0.) bvf2(i,j)=0. !! Set unstable BVF to neutral for z_plume computation
               exit
           endif
       enddo

       !! Step 1: Compute plume height with initial Sofiev parameters
       z_plume(i,j) = (sof_a_0*z_pbl_pt(i,j))+(sof_b_0*((FRPi/sof_P0)**sof_c_0))* &
                  exp(-sof_d_0*bvf2(i,j)/sof_N0)

     endif

     if (f_scheme.eq."SOFIEV-2STEP".and.FRPi.gt.0.) then

       !! Step 0: Compute BVfreq for twice the height of the PBL
       do l = kd,2,-1
           if (z_pbl_pt(i,j)*2. .lt. z_full_pt(i,j,kd)) z_plume(i,j)=0.0
           if (z_pbl_pt(i,j)*2. .ge. z_full_pt(i,j,l) .and. z_pbl_pt(i,j)*2. .lt. z_half_pt(i,j,l-1)) then
               ktop=l-1
               kbot=l
               bvf2(i,j)=brunt_vaisala(pfull_pt(i,j,:),temp_pt(i,j,:),ktop,kbot)
               if (bvf2(i,j).lt.0.) bvf2(i,j)=0. !! Set unstable BVF to neutral for z_plume computation
               exit
           endif
       enddo

       !! Step 1: Compute plume height with initial Sofiev parameters
       z_plume(i,j) = (sof_a_00*z_pbl_pt(i,j))+(sof_b_00*((FRPi/sof_P00)**sof_c_00))* &
                  exp(-sof_d_00*bvf2(i,j)/sof_N00)

       !! Step 2: Check whether the plume height is above or below the PBL and
       !!         recompute the plume height for FT or PBL injection
       if (z_plume(i,j).gt.z_pbl_pt(i,j)) then !!! FT plume
           z_plume(i,j) = (sof_a_1*z_pbl_pt(i,j))+(sof_b_1*((FRPi/sof_P00)**sof_c_1))* &
                      exp(-sof_d_1*bvf2(i,j)/sof_N0)
       else if (z_plume(i,j).le.z_pbl_pt(i,j).and.z_plume(i,j).gt.0.) then !!! PBL plume
           z_plume(i,j) = (sof_a_2*z_pbl_pt(i,j))+(sof_b_2*((FRPi/sof_P00)**sof_c_2))* &
                      exp(-sof_d_2*bvf2(i,j)/sof_N0)
       endif

     endif
     
     if (f_scheme.eq."SOFIEV-SFRP".and.FRPi.gt.0.) then

       !! Step 1: Compute BVfreq for twice the height of the PBL
       do l = kd,2,-1
           if (z_pbl_pt(i,j)*2. .lt. z_full_pt(i,j,kd)) z_plume=0.0
           if (z_pbl_pt(i,j)*2. .ge. z_full_pt(i,j,l) .and. z_pbl_pt(i,j)*2. .lt. z_half_pt(i,j,l-1)) then
               ktop=l-1
               kbot=l
               bvf2(i,j)=brunt_vaisala(pfull_pt(i,j,:),temp_pt(i,j,:),ktop,kbot)
               if (bvf2(i,j).lt.0.) bvf2(i,j)=0. !! Set unstable BVF to neutral for z_plume computation
               exit
           endif
       enddo

       !! Step 2: Compute plume height with initial Sofiev parameters
       z_plume(i,j) = (sof_a_0*z_pbl_pt(i,j))+(sof_b_0*((FRPi/sof_P0)**sof_c_0))* &
                  exp(-sof_d_0*bvf2(i,j)/sof_N0)

       if (z_plume(i,j).gt.1500) then !!! 
           FRPi = FRPi * SQRT(z_plume(i,j) / 1500) !!! FRP subgrid recalculation based on Veira et al 2015
           z_plume(i,j) = (sof_a_1*z_pbl_pt(i,j))+(sof_b_1*((FRPi/sof_P00)**sof_c_1))* &
                      exp(-sof_d_1*bvf2(i,j)/sof_N0)
       endif

     endif


     do l = kd,1,-1
       Z0 = z_half_pt(i,j,kd+1)
       Z1 = z_half_pt(i,j,l+1)
       Z2 = z_half_pt(i,j,l)
      
       if (f_scheme.eq."EVEN_PBL".or.FRPi.eq.0.) then
          if (z_pbl_pt(i,j).le.Z0) fbbl(i,j,l) = 1.
          if (z_pbl_pt(i,j).le.Z1.or.z_pbl_pt(i,j).eq.0.) exit
          if (z_pbl_pt(i,j).ge.Z2) fbbl(i,j,l) = (Z2-Z1)/z_pbl_pt(i,j)
          if (z_pbl_pt(i,j).gt.Z1.and.z_pbl_pt(i,j).lt.Z2) fbbl(i,j,l) = (z_pbl_pt(i,j)-Z1)/z_pbl_pt(i,j)

       else if (f_scheme.eq."VALMARTIN") then

          do lf = 1, 25
            if (Z1 < alt_fire_max(lf) .and. Z2 > alt_fire_min(lf)) then
                if (Z1 >= alt_fire_min(lf)) then
                    if (Z2 < alt_fire_max(lf)) then
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * (Z2 - Z1) / (alt_fire_max(lf) - alt_fire_min(lf))
                    else
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * (alt_fire_max(lf) - Z1) / (alt_fire_max(lf) - alt_fire_min(lf))
                    end if
                else
                    if (Z2 <= alt_fire_max(lf)) then
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * (Z2 - alt_fire_min(lf)) / (alt_fire_max(lf) - alt_fire_min(lf))
                    else
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * 1.0
                    end if
                end if
            end if
          end do

       else if (f_scheme.eq."SURFACE") then
          fbbl(i,j,kd)=1.

       else if (f_scheme.eq."AEROCOM") then

          do lf=1,7
            if (Z1 < alt_fire_max(lf) .and. Z2 > alt_fire_min(lf)) then
                if (Z1 >= alt_fire_min(lf)) then
                    if (Z2 < alt_fire_max(lf)) then
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * (Z2 - Z1) / (alt_fire_max(lf) - alt_fire_min(lf))
                    else
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * (alt_fire_max(lf) - Z1) / (alt_fire_max(lf) - alt_fire_min(lf))
                    end if
                else
                    if (Z2 <= alt_fire_max(lf)) then
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * (Z2 - alt_fire_min(lf)) / (alt_fire_max(lf) - alt_fire_min(lf))
                    else
                        fbbl(i, j, l) = fbbl(i, j, l) + fbb_clim(lf,i, j) * 1.0
                    end if
                end if
            end if
          end do
            
       else if (f_scheme.eq."SOFIEV".or.f_scheme.eq."SOFIEV-CF".or.f_scheme.eq."SOFIEV-2STEP".or.f_scheme.eq."SOFIEV-SFRP".and.FRPi.gt.0) then
          if (z_plume(i,j).le.Z0) fbbl(i,j,l) = 1.
          if (z_plume(i,j).le.Z1.or.z_plume(i,j).eq.0.) exit
          if (z_plume(i,j).ge.Z2) fbbl(i,j,l) = (Z2-Z1)/z_plume(i,j)
          if (z_plume(i,j).gt.Z1.and.z_plume(i,j).lt.Z2) fbbl(i,j,l) = (z_plume(i,j)-Z1)/z_plume(i,j)
       else
          call ERROR_MESG('fire_emiss', 'Fire_injhgt_scheme option in namelist is not available.', FATAL )
       endif

     enddo

!!!! conservation of mass of fbbl
     fbb_tot = sum(fbbl(i,j,:))
     if (fbb_tot.gt.0 .or. fbb_tot.lt.0) then
         fbbl(i,j,:) = fbbl(i,j,:)/fbb_tot
     else
         fbbl(i,j,kd) = 1.
     endif 

   enddo   !end of i loop
   enddo   !end of j loop

end subroutine fire_fbb
!</SUBROUTINE>

!<FUNCTION> Calculate squared brunt vaisala frequency
function brunt_vaisala(pfull,temp,kt,kb) result(bvf2)
    real :: bvf2
    real :: pfull(:)
    real :: temp(:)
    integer :: kt,kb

    ! Local variables
    real :: theta_kt, theta_kb
    real, parameter :: p00 = 1.e5

    theta_kt=temp(kt)*(pfull(kt)/p00)**(-KAPPA)
    theta_kb=temp(kb)*(pfull(kb)/p00)**(-KAPPA)
    
    !!! N=-g^2*P_ave*d(theta)/d(P)*(1/R*T_ave*theta_ave)
    bvf2=-(((GRAV**2.)*((pfull(kt)+pfull(kb))/2.)* &
          (theta_kt-theta_kb)/(pfull(kb)-pfull(kt)))/ &
          (RDGAS*((theta_kb+theta_kt)/2.)*((temp(kt)+temp(kb))/2.)))

end function

!</SUBROUTINE>
!#######################################################################
!
!<SUBROUTINE NAME="fire_emiss">
subroutine atmos_fire_plumerise_driver(fbb,pfull_pt,temp_pt, &
                      z_half_pt,z_pbl_pt,z_full_pt, &
                      tr, diag_time, is, ie, js, je, local_hour_2d)
   real, intent(in), dimension(:,:,:)    :: pfull_pt
   real, intent(in), dimension(:,:,:)    :: temp_pt
   real, intent(in), dimension(:,:,:)    :: z_half_pt
   real, intent(in), dimension(:,:,:)    :: z_full_pt
   real, intent(in), dimension(:,:)      :: z_pbl_pt
   real, intent(out), dimension(:,:,:) :: fbb
   integer, intent(in)                    :: is, ie, js, je
   real, intent(in),  dimension(:,:,:) :: tr
   type(time_type), intent(in)            :: diag_time
   real, intent(in), dimension(:,:)    :: local_hour_2d   
!  real, dimension(size(tr,1),size(tr,2),size(tr,3)) :: fbb_norm
   character(len=80)     :: f_scheme
   integer :: lf, n, np, j, i, k, id, jd, kd, npercentiles
   real, dimension(num_percentiles,size(tr,1),size(tr,2)) :: fire_intensity_perc
   real, dimension(25,size(tr,1),size(tr,2)) :: fbb_clim_misr
   real, dimension(7,size(tr,1),size(tr,2)) :: fbb_clim_aerocom
   real, dimension(size(tr,1),size(tr,2),size(tr,3),6) :: fbb_perc
   real, dimension(size(tr,1),size(tr,2)) :: z_plume,bvf2, fire_intensity_diu
   real, dimension(size(tr,1),size(tr,2)) :: diurnal_scale_factor, fbb_FT, fbb_NFT
   real, dimension(6)    :: alt_fire_min_aerocom
   real, dimension(6)    :: alt_fire_max_aerocom
   real, dimension(24)    :: alt_fire_min_misr
   real, dimension(24)    :: alt_fire_max_misr
   integer :: nlevel_fire

   id=size(tr,1); jd=size(tr,2); kd=size(tr,3)
   fbb(:,:,:) = 0.0
   z_plume(:,:) = 0.0
   bvf2(:,:) = 0.0
   fbb_perc(:,:,:,:) = 0.0  
   fire_intensity_perc(:,:,:) = 0.0
   fbb_clim_misr(:,:,:) = 0.0 
   fbb_clim_aerocom(:,:,:) = 0.0
   f_scheme = fire_injhgt_scheme

   npercentiles = num_percentiles



   if (f_scheme.eq.'SOFIEV'.or.f_scheme.eq.'SOFIEV-2STEP'.or.f_scheme.eq.'SOFIEV-SFRP'.and.trim(frp_source).ne.' ') then
     select case (trim(frp_source))
       case ('MODIS')
         do lf=1, npercentiles
           call interpolator(frp_value_interp, frp_time, fire_intensity_perc(lf,:,:), &
                             trim(frp_value_name(lf)), is, js)
         end do
     end select
   !!! percentiles

   do np = 1, npercentiles !!! default 6 sections /0-10/10-25/25-50/50-75/75-90/90-99 
   !
   ! Calculate fraction of emission at every level for open fires
   !
   !!! diurnal fire_intensity_perc calculation
        if (do_frp_diurnal) then
                fire_intensity_diu = 0

                call atmos_fire_frp_diurnal(fire_intensity_perc(np,:,:), fire_intensity_diu, tr,  &
                        local_hour_2d, id, jd, diag_time)

                fire_intensity_perc(np,:,:) = fire_intensity_diu(:,:)
        endif

        call fire_fbb(fbb_perc(:,:,:,np),z_plume(:,:),bvf2(:,:),f_scheme,kd,id,jd,pfull_pt(:,:,:),temp_pt(:,:,:), &
                     z_half_pt(:,:,:), z_pbl_pt(:,:),z_full_pt(:,:,:),fire_intensity_perc(np,:,:),nlevel_fire, &
                     alt_fire_max_aerocom,alt_fire_min_aerocom, fbb_clim_aerocom(:,:,:))
   enddo


   !!! final fbb equals to summation of all fbb_perc s based on their share
   do np = 1, npercentiles
     do j = 1, jd
       do i = 1, id
          fbb(i,j,:) = fbb(i,j,:) + (perc_share(np)*fbb_perc(i,j,:,np))
       enddo
     enddo
   enddo
 

   else if (f_scheme.eq.'VALMARTIN'.and.trim(fbb_source).ne.' ') then

       do lf=1,25
          write(fbb_value_name(lf), '("perc_inj_", I2.2)') lf
          call interpolator(fbb_value_interp, fbb_time, fbb_clim_misr(lf,:,:), &
                         trim(fbb_value_name(lf)), is, js)
       end do
       f_scheme = 'VALMARTIN'
       alt_fire_min_misr=(/0.0, 250.0, 500.0, 750.0, 1000.0, 1250.0, 1500.0, &
                              1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, &
                              3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, &
                              4750.0, 5000.0, 5250.0, 5500.0, 5750.0/)
       alt_fire_max_misr=(/250.0, 500.0, 750.0, 1000.0, 1250.0, 1500.0, &
                              1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, & 
                              3250.0, 3500.0, 3750.0, 4000.0, 4250.0, 4500.0, &
                              4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0/)

       call fire_fbb(fbb(:,:,:),z_plume(:,:),bvf2(:,:),f_scheme,kd,id,jd,pfull_pt(:,:,:),temp_pt(:,:,:), &
                     z_half_pt(:,:,:), z_pbl_pt(:,:),z_full_pt(:,:,:),fire_intensity_perc(np,:,:),nlevel_fire, &
                     alt_fire_max_misr,alt_fire_min_misr, fbb_clim_misr(:,:,:)/100)
   
   else if (f_scheme.eq.'AEROCOM'.and.trim(fbb_source).ne.'') then

        do lf=1,7
           write(fbb_value_name(lf), '("perc_inj_", I2.2)') lf
           call interpolator(fbb_value_interp, fbb_time, fbb_clim_aerocom(lf,:,:), &
                          trim(fbb_value_name(lf)), is, js)
        end do
        f_scheme = 'AEROCOM'
        alt_fire_min_aerocom=(/0.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0/)
        alt_fire_max_aerocom=(/100.0, 500.0, 1000.0, 2000.0, 3000.0, 6000.0/)

       call fire_fbb(fbb(:,:,:),z_plume(:,:),bvf2(:,:),f_scheme,kd,id,jd,pfull_pt(:,:,:),temp_pt(:,:,:), &
                     z_half_pt(:,:,:), z_pbl_pt(:,:),z_full_pt(:,:,:),fire_intensity_perc(np,:,:),nlevel_fire, &
                     alt_fire_max_aerocom,alt_fire_min_aerocom, fbb_clim_aerocom(:,:,:)/100)

   else if (f_scheme.eq.'SURFACE'.or.f_scheme.eq.'EVEN_PBL') then

       call fire_fbb(fbb(:,:,:),z_plume(:,:),bvf2(:,:),f_scheme,kd,id,jd,pfull_pt(:,:,:),temp_pt(:,:,:), &
                     z_half_pt(:,:,:), z_pbl_pt(:,:),z_full_pt(:,:,:),fire_intensity_perc(np,:,:),nlevel_fire, &
                     alt_fire_max_aerocom,alt_fire_min_aerocom, fbb_clim_aerocom(:,:,:))

   endif
     

!!! fire_emis_diunal
! fbb_norm = fbb
   if (do_bb_emis_diurnal) then
     call atmos_fire_emis_diurnal(diurnal_scale_factor, tr,  &
                        local_hour_2d, id, jd, diag_time)
     do k = 1, kd
          fbb(:,:,k) = fbb(:,:,k) * diurnal_scale_factor(:,:)
     enddo
   endif 


!     if (do_bb_emis_diurnal) then
!       if (id_fbb > 0) then
!             used = send_data(id_fbb,fbb_norm, diag_time, &
!             is_in=is,js_in=js,ks_in=1)
!       endif
!     else
        if (id_fbb > 0) then
              used = send_data(id_fbb,fbb, diag_time, &
              is_in=is,js_in=js,ks_in=1)
        endif              
!     endif
      do n = 1,npercentiles
         if (id_fbb_perc(n) > 0) then
            used = send_data ( id_fbb_perc(n), fbb_perc(:,:,:,n), diag_time, &
                 is_in=is,js_in=js,ks_in=1)
         endif
         if (id_frp_perc(n) > 0) then
            used = send_data ( id_frp_perc(n), fire_intensity_perc(n,:,:), diag_time, &
                 is_in=is,js_in=js)
         endif
      end do
      if (id_injhgt > 0) then
        used = send_data(id_injhgt,z_plume, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_FRP > 0) then
        used = send_data(id_FRP,fire_intensity_perc(npercentiles,:,:), diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_bvf2 > 0) then
        used = send_data(id_bvf2,bvf2, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_pfull_pt > 0) then
        used = send_data(id_pfull_pt,pfull_pt, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_temp_pt > 0) then
        used = send_data(id_temp_pt,temp_pt, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_z_half_pt > 0) then
        used = send_data(id_z_half_pt,z_half_pt, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_z_pbl_pt > 0) then
        used = send_data(id_z_pbl_pt,z_pbl_pt, diag_time, &
              is_in=is,js_in=js)
      endif
      if (id_z_full_pt > 0) then
        used = send_data(id_z_full_pt,z_full_pt, diag_time, &
              is_in=is,js_in=js,ks_in=1)
      endif
        
      if (id_fbb_FT > 0 .and. id_fbb_NFT > 0) then
         do j = 1, jd
            do i = 1, id
               fbb_FT(i,j) = 0.0
               fbb_NFT(i,j) = 0.0
               do k = 1, kd
                  if (z_half_pt(i,j,k) > z_pbl_pt(i,j)) then
!                    fbb_FT(i,j) = fbb_FT(i,j) + fbb_norm(i,j,k)
                     fbb_FT(i,j) = fbb_FT(i,j) + fbb(i,j,k)
                  else
!                    fbb_NFT(i,j) = fbb_NFT(i,j) + fbb_norm(i,j,k)
                     fbb_NFT(i,j) = fbb_NFT(i,j) + fbb(i,j,k)
                  endif
               enddo
            enddo
         enddo
        used = send_data(id_fbb_FT,fbb_FT, diag_time, &
              is_in=is,js_in=js)
        used = send_data(id_fbb_NFT,fbb_NFT, diag_time, &
              is_in=is,js_in=js)

      endif



end subroutine atmos_fire_plumerise_driver
! ===========================================================================

subroutine atmos_fire_emis_diurnal(emisbb_scale_factor, tr, local_hour_2d, id, jd, Time)

   real, intent(out), dimension(:,:)    :: emisbb_scale_factor
   real, intent(in),  dimension(:,:,:) :: tr   
   real, intent(in), dimension(:,:)    :: local_hour_2d
   integer, intent(in) :: id
   integer, intent(in) :: jd
   type(time_type), intent(in)            :: Time 
   real    :: ru_base, ru_peak, h_peak, sigma
   integer :: i,j
   real, dimension(size(tr,1),size(tr,2))    :: share_bb_percent
   
   emisbb_scale_factor(:,:) = 0

   !!! fitted a gaussian function on fig10a Li et al. 2019 (constraining the integral of function to 1 for conservation)
   !share_bb_function = ru_base + (ru_peak - ru_base) * exp(-(local_hour - h_peak)**2 / (2*sigma**2))
!   ru_base = 0.014216161324495054
!   ru_peak = 0.10911740870082014
!   h_peak = 13.258881845262565
!   sigma = 2.830072893797371
    ru_base = 0.013667529633792457
    ru_peak = 0.10856338628471231
    h_peak = 13.256384792906122
    sigma = 2.8253815990151745

   do j=1,jd
      do i=1,id
         share_bb_percent(i,j) = ru_base + (ru_peak - ru_base) * exp(-(local_hour_2d(i,j) - h_peak)**2 / (2*sigma**2))
         emisbb_scale_factor(i,j) = share_bb_percent(i,j) * 24 !!! conservation of mass x24
      end do
   end do

end subroutine atmos_fire_emis_diurnal
! ===========================================================================
function atmos_fire_do_bb_emis_diurnal() result(shared_logic)
                        
   logical :: shared_logic

   shared_logic = do_bb_emis_diurnal   
end function atmos_fire_do_bb_emis_diurnal
! ===========================================================================

subroutine atmos_fire_frp_diurnal(fire_arr, fire_diurnal, tr, local_hour_2d, id, jd, Time)

   real, intent(in), dimension(:,:)    :: fire_arr
   real, intent(out), dimension(:,:)    :: fire_diurnal
   real, intent(in),  dimension(:,:,:) :: tr
   real, intent(in), dimension(:,:)    :: local_hour_2d
   integer, intent(in) :: id
   integer, intent(in) :: jd   
   type(time_type), intent(in)            :: Time
   real    :: ru_base, ru_peak, h_peak, sigma
   integer :: i,j
   logical :: do_frp_diurnal_land_cats = .False.
   real, dimension(size(tr,1),size(tr,2)) ::    frp_normalized
   real, dimension(size(tr,1),size(tr,2)) :: land_cat

   !!! fitted gaussian functions on fig5 Li et al. 2022
   !share_bb_function = ru_base + (ru_peak - ru_base) * exp(-(local_hour - h_peak)**2 / (2*sigma**2))
   
   land_cat(:,:) = 0

   if (do_frp_diurnal_land_cats) then

     do j=1,jd
     do i=1,id
        if (land_cat(i,j) .eq. 0) then !!! land_cat = 0 : Water
           ru_base = 0
           ru_peak = 0
           h_peak = 1
           sigma = 1

        else if (land_cat(i,j) .eq. 1) then !!! land_cat = 1 : Forest
           ru_base = 0.5437192879498844
           ru_peak = 1.0050662881155534
           h_peak = 14.525800008024998
           sigma = 3.8603934622712863

        else if (land_cat(i,j) .eq. 2) then !!! land_cat = 2 : Shrubland
           ru_base = 0.4612955690530021
           ru_peak = 1.0010467748216656
           h_peak = 15.865467674960273
           sigma = 3.5382915261007435

        else if (land_cat(i,j) .eq. 3) then !!! land_cat = 3 : Savanna
           ru_base = 0.3578797020002618
           ru_peak = 0.9940494698442988
           h_peak = 15.209920437512379
           sigma = 3.9182059269742497

        else if (land_cat(i,j) .eq. 4) then !!! land_cat = 4 : Grassland
           ru_base = 0.48662785258430097
           ru_peak = 1.0049787306837865
           h_peak = 15.6409803300218
           sigma = 4.056755251648667

        else if (land_cat(i,j) .eq. 5) then !!! land_cat = 5 : Cropland
           ru_base = 0.4872797864699511
           ru_peak = 1.0197114730579282
           h_peak = 16.042003479419357
           sigma = 4.331292604527208 
        endif

        frp_normalized(i,j) = ru_base + (ru_peak - ru_base) * exp(-(local_hour_2d(i,j) - h_peak)**2 / (2*sigma**2))
        fire_diurnal(i,j) = fire_arr(i,j) * (frp_normalized(i,j))
     end do
     end do
   
   else   !!! general: averaged for all land categories

     ru_base = 0.46709831123362405
     ru_peak = 1.0008447676795196
     h_peak = 15.474796621482293
     sigma = 3.98145297649462

     do j=1,jd
     do i=1,id
         frp_normalized(i,j) = ru_base + (ru_peak - ru_base) * exp(-(local_hour_2d(i,j) - h_peak)**2 / (2*sigma**2))
         fire_diurnal(i,j) = fire_arr(i,j) * (frp_normalized(i,j))
     end do
     end do      

   endif

end subroutine atmos_fire_frp_diurnal
! ===========================================================================
end module atmos_fire_plumerise_mod
