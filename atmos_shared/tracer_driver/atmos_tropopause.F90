module atmos_tropopause_mod
! <CONTACT EMAIL="Larry.Horowitz@noaa.gov">
!   Larry.Horowitz
! </CONTACT>

! <OVERVIEW>
!   Tropopause diagnostics
! </OVERVIEW>

! <DESCRIPTION>
!   Tropopause diagnostics
! </DESCRIPTION>


use              fms_mod, only : write_version_number,    &
                                 mpp_pe, mpp_root_pe,                 &
                                 stdlog, stdout,          &
                                 check_nml_error, error_mesg,         &
                                 FATAL, NOTE, WARNING
use          fms2_io_mod, only : file_exists
use     diag_manager_mod, only : send_data
use atmos_cmip_diag_mod,  only : register_cmip_diag_field_2d
use     time_manager_mod, only : time_type
use              mpp_mod, only : input_nml_file 


implicit none

private

!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_tropopause, &
        atmos_tropopause_init

!-----------------------------------------------------------------------
!----------- module data -------------------
!-----------------------------------------------------------------------

character(len=48), parameter :: module_name = 'tracers'

logical :: module_is_initialized = .FALSE.
integer :: logunit

real, parameter    :: ztrop_low  = 5.e3   ! lowest tropopause level allowed (m)
real, parameter    :: ztrop_high = 20.e3  ! highest tropopause level allowed (m)
real, parameter    :: max_dtdz   = 2.e-3  ! max dt/dz for tropopause level (K/m)
real, parameter    :: ptrop_low  = 450.e2 ! lowest tropopause level allowed (Pa)
real, parameter    :: ptrop_high = 85.e2  ! highest tropopause level allowed (Pa)
real, parameter    :: deltaz     = 2.e3   ! depth to check lapse rate (m)

integer            :: tropopause_scheme_code
integer, parameter :: ESM41_TROP=1, WMO_TROP=2

!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer :: id_ptp, id_tatp, id_ztp

!---------------------------------------------------------------------
!-------- namelist  ---------
!-----------------------------------------------------------------------

logical  :: do_tropopause_diagnostics = .false.
character(len=64)  :: tropopause_scheme = 'ESM41'

namelist /atmos_tropopause_nml/  &
          do_tropopause_diagnostics, tropopause_scheme

!---- version number -----
character(len=128) :: version = '$$'
character(len=128) :: tagname = '$$'

!-----------------------------------------------------------------------

contains

!#######################################################################

!<SUBROUTINE NAME ="atmos_tropopause">
!<OVERVIEW>
!  A subroutine to calculate tropopause diagnostics.
!
! do_co2_restore   = logical to turn co2_restore on/off: default = .false.
! restore_co2_dvmr = partial pressure of co2 to which to restore  (mol/mol)
! restore_klimit   = atmospheric level to which to restore starting from top
! restore_tscale   = timescale in seconds with which to restore
!
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate tropopause diagnostics.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_tropopause (Time, Time_next, t, pfull)
!</TEMPLATE>
!
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="Time_next" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="t" TYPE="real" DIM="(:,:,:)">
!     Temperature.
!   </IN>
!   <IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!     Pressures on the model full levels.
!   </IN>
!   <IN NAME="z_full" TYPE="real" DIM="(:,:,:)">
!     Height of the model full levels.
!   </IN>

subroutine atmos_tropopause(is, ie, js, je, Time, Time_next, t, pfull, z_full, &
                            tropopause_ind)

   integer, intent(in)                   :: is, ie, js, je
   type (time_type), intent(in)          :: Time, Time_next
   real,    intent(in), dimension(:,:,:) :: t            ! K
   real,    intent(in), dimension(:,:,:) :: pfull        ! Pa
   real,    intent(in), dimension(:,:,:) :: z_full       ! m
   integer, intent(out), dimension(:,:)  :: tropopause_ind ! 1
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

integer   :: i,j,k,id,jd,kd,kk
real      :: dtemp
logical   :: used, found

real, dimension(size(t,1),size(t,2)) :: ptp, tatp, ztp

!-----------------------------------------------------------------------

    if (.not. module_is_initialized)  &
       call error_mesg ('Atmos_tropopause','atmos_tropopause_init must be called first.', FATAL)

    id=size(t,1); jd=size(t,2); kd=size(t,3)

    select case (tropopause_scheme_code)
! scheme used in ESM4.1
     case (ESM41_TROP)

       do j=1,jd
       do i=1,id

       do k = kd-1,2,-1
          if (z_full(i,j,k) < ztrop_low ) then
              cycle
          else if( z_full(i,j,k) > ztrop_high ) then
              tropopause_ind(i,j)    = k
              exit
          end if
          dtemp = t(i,j,k) - t(i,j,k-1)
          if( dtemp < max_dtdz*(z_full(i,j,k-1) - z_full(i,j,k)) ) then
             tropopause_ind(i,j)    = k
             exit
          end if
       end do

       ptp(i,j) = pfull(i,j,tropopause_ind(i,j))
       tatp(i,j) = t(i,j,tropopause_ind(i,j))
       ztp(i,j) = z_full(i,j,tropopause_ind(i,j))

       end do
       end do

     case (WMO_TROP)
! WMO scheme (International Meteorological Vocabulary, WMO, 182, 1992)
! defined as the lowest level at which the lapse rate
! decreases to 2degC km-1 or less, provided that the average
! lapse rate between this level and all higher levels within
! 2 km does not exceed 2degC km-1.
       do j=1,jd
       do i=1,id

       found = .false.
       do k = kd-1,2,-1
          if (pfull(i,j,k) > ptrop_low ) then
              cycle
          else if( pfull(i,j,k) < ptrop_high ) then
              tropopause_ind(i,j)    = k
              exit
          end if
          dtemp = t(i,j,k) - t(i,j,k-1)
          if( dtemp < max_dtdz*(z_full(i,j,k-1) - z_full(i,j,k)) ) then
! check that lapse rate LT 2K/km within 2km above
             do kk = k-2,1,-1
                if (z_full(i,j,kk)-z_full(i,j,k) > deltaz) then
                   tropopause_ind(i,j) = k
                   found = .true.
                   exit
                end if
                dtemp = t(i,j,k) - t(i,j,kk)
                if( dtemp > max_dtdz*(z_full(i,j,kk) - z_full(i,j,k)) ) then
                   exit ! search for another tropopause candidate
                end if
             end do
             if (found) exit
          end if
       end do

       ptp(i,j)  = pfull(i,j,tropopause_ind(i,j))
       tatp(i,j) = t(i,j,tropopause_ind(i,j))
       ztp(i,j)  = z_full(i,j,tropopause_ind(i,j))

       enddo
       enddo

     case default
       call error_mesg ('atmos_tropopause_init', 'undefined tropopause scheme ', FATAL )
    end select

    if (id_ptp > 0) &
       used = send_data (id_ptp, ptp, Time_next, is_in=is,js_in=js)

    if (id_tatp > 0) &
       used = send_data (id_tatp, tatp, Time_next, is_in=is,js_in=js)

    if (id_ztp > 0) &
       used = send_data (id_ztp, ztp, Time_next, is_in=is,js_in=js)

end subroutine atmos_tropopause
!</SUBROUTINE >



!#######################################################################

!<SUBROUTINE NAME ="atmos_tropopause_init">

!<OVERVIEW>
! Subroutine to initialize the tropopause diagnostics module.
!</OVERVIEW>

 subroutine atmos_tropopause_init (Time)

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!
type(time_type),       intent(in)                   :: Time
!
!-----------------------------------------------------------------------
!     local variables
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!-----------------------------------------------------------------------
!
integer :: ierr, io
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

    if (module_is_initialized) return

    call write_version_number (version, tagname)

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
    if ( file_exists('input.nml')) then
        read (input_nml_file, nml=atmos_tropopause_nml, iostat=io)
        ierr = check_nml_error(io,'atmos_tropopause_nml')
    end if

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
    logunit=stdlog()
    if (mpp_pe() == mpp_root_pe() ) &
              write (logunit, nml=atmos_tropopause_nml)

    if ( TRIM(tropopause_scheme) == 'ESM41') then
       tropopause_scheme_code = ESM41_TROP
    else if ( TRIM(tropopause_scheme) == 'WMO') then
       tropopause_scheme_code = WMO_TROP
    else
       call error_mesg ('atmos_tropopause_init', 'undefined tropopause scheme '// &
                        TRIM(tropopause_scheme), FATAL )
    end if

    id_ptp  = register_cmip_diag_field_2d ( module_name, 'ptp', Time, &
                long_name='Tropopause Air Pressure', units='Pa', &
                standard_name='tropopause_air_pressure')

    id_tatp = register_cmip_diag_field_2d ( module_name, 'tatp', Time, &
                long_name='Tropopause Air Temperature', units='K', &
                standard_name='tropopause_air_temperature')

    id_ztp  = register_cmip_diag_field_2d ( module_name, 'ztp', Time, &
                long_name='Tropopause Altitude', units='m', &
                standard_name='tropopause_altitude')

    call write_version_number (version, tagname)
    module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

end subroutine atmos_tropopause_init
!</SUBROUTINE>

end module atmos_tropopause_mod
