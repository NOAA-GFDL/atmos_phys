!FDOC_TAG_GFDL

module cu_mo_trans_mod
! <CONTACT EMAIL="Isaac.Held@noaa.gov">
!  Isaac Held
! </CONTACT>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    A simple module that computes a diffusivity proportional to the 
!    convective mass flux, for use with diffusive 
!    convective momentum transport closure
! </OVERVIEW>
! <DESCRIPTION>
!   A diffusive approximation to convective momentum transport is crude but
!    has been found to be useful in improving the simulation of tropical
!     precipitation in some models.  The diffusivity computed here is
!     simply 
!<PRE>
! diffusivity = c*W*L 
! W = M/rho  (m/sec) 
! M = convective mass flux (kg/(m2 sec)) 
! rho - density of air <p>
! L = depth of convecting layer (m)
! c = normalization constant = diff_norm/g 
!   (diff_norm is a namelist parameter;
!      the factor of g = acceleration of gravity here is an historical artifact) <p>
! for further discussion see 
!     <LINK SRC="cu_mo_trans.pdf">cu_mo_trans.pdf</LINK>
!</PRE>
! </DESCRIPTION>


!=======================================================================
!
!                 DIFFUSIVE CONVECTIVE MOMENTUM TRANSPORT MODULE
!
!=======================================================================

  use   constants_mod, only:  GRAV, RDGAS, RVGAS, CP_AIR
 

  use         mpp_mod, only: input_nml_file
  use         fms_mod, only: check_nml_error,    &
                             write_version_number,           &
                             mpp_pe, mpp_root_pe, stdlog,    &
                             error_mesg, FATAL, NOTE
  use       fms2_io_mod, only: file_exists
  use  Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use  Time_Manager_Mod, ONLY: time_type

 use convection_utilities_mod, only : conv_results_type
 use moist_proc_utils_mod, only: mp_input_type, mp_output_type, mp_nml_type
implicit none
private


! public interfaces
!=======================================================================
public :: cu_mo_trans_init, &
          cu_mo_trans,      &
          cu_mo_trans_end

!=======================================================================

! form of interfaces
!=======================================================================

      
logical :: module_is_initialized = .false.


!---------------diagnostics fields------------------------------------- 

integer :: id_diff_cmt, id_massflux_cmt

character(len=11) :: mod_name = 'cu_mo_trans'

real :: missing_value = -999.
logical  ::  do_diffusive_transport = .false.
logical  ::  do_nonlocal_transport = .false.


!--------------------- namelist variables with defaults -------------

real    :: diff_norm =   1.0
logical :: limit_mass_flux = .false.  ! when true, the mass flux 
                                      ! out of a grid box is limited to
                                      ! the mass in that grid box
character(len=64) :: transport_scheme = 'diffusive'
logical :: conserve_te = .true.  ! conserve total energy ?
real    ::  gki = 0.7  ! Gregory et. al. constant for p-gradient param
real    :: amplitude = 1.0 ! Tuning parameter (1=full strength)

namelist/cu_mo_trans_nml/ diff_norm, &
                          limit_mass_flux, &
                          conserve_te, gki, &
                          amplitude,  &
                          transport_scheme


!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: cmt_uses_uw
logical :: do_uw_conv



contains

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans_init">
!  <OVERVIEW>
!   initializes module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Reads namelist and registers one diagnostic field
!     (diff_cmt:  the kinematic diffusion coefficient)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans_init( axes, Time )
!
!  </TEMPLATE>
!  <IN NAME=" axes" TYPE="integer">
!    axes identifier needed by diag manager
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!    time at initialization needed by diag manager
!  </IN>
! </SUBROUTINE>
!
subroutine cu_mo_trans_init( axes, Time, Nml_mp, cmt_mass_flux_source)

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 type(mp_nml_type), intent(in) :: Nml_mp
 character(len=64), intent(in) :: cmt_mass_flux_source

integer :: ierr, io, logunit
integer, dimension(3)  :: half =  (/1,2,4/)

      do_uw_conv = Nml_mp%do_uw_conv


!------ read namelist ------

   if ( file_exists('input.nml')) then
      read (input_nml_file, nml=cu_mo_trans_nml, iostat=io)
      ierr = check_nml_error(io,'cu_mo_trans_nml')
   endif

!--------- write version number and namelist ------------------

      call write_version_number ( version, tagname )
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
        write ( logunit, nml=cu_mo_trans_nml )

!----------------------------------------------------------------------
!    define logicals indicating momentum transport scheme to use.
!----------------------------------------------------------------------
      if (trim(transport_scheme) == 'diffusive') then
        do_diffusive_transport = .true.
      else if (trim(transport_scheme) == 'nonlocal') then
        do_nonlocal_transport = .true.
      else
        call error_mesg ('cu_mo_trans', &
         'invalid specification of transport_scheme', FATAL)
      endif

! --- initialize quantities for diagnostics output -------------

   if (do_diffusive_transport) then
     id_diff_cmt = &
      register_diag_field ( mod_name, 'diff_cmt', axes(1:3), Time,    &
                        'cu_mo_trans coeff for momentum',  'm2/s', &
                         missing_value=missing_value               )
     id_massflux_cmt = &
      register_diag_field ( mod_name, 'massflux_cmt', axes(half), Time, &
                        'cu_mo_trans mass flux',  'kg/(m2 s)', &
                         missing_value=missing_value               )
    endif

!--------------------------------------------------------------
        if (trim(cmt_mass_flux_source) == 'uw') then
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('convection_driver_init', &
                'if cmt_uses_uw = T, then do_uw_conv must be T', FATAL)
          endif
        else if (trim(cmt_mass_flux_source) == 'all') then
          if (do_uw_conv)  then
            cmt_uses_uw = .true.
          else
            cmt_uses_uw = .false.
          endif
        else
          call error_mesg ('convection_driver_init', &
             'invalid specification of cmt_mass_flux_source', FATAL)
        endif

  module_is_initialized = .true.


end subroutine cu_mo_trans_init

!#########################################################################

subroutine cu_mo_trans ( is, js, Time, dt, num_tracers, Input_mp,  &
                         Conv_results, Output_mp, ttnd_conv)

type(time_type),          intent(in)    :: Time
integer,                  intent(in)    :: is, js, num_tracers
real,                     intent(in)    :: dt
type(mp_input_type),      intent(inout)    :: Input_mp
type(mp_output_type),     intent(inout) :: Output_mp
real, dimension(:,:,:),   intent(inout) ::  ttnd_conv
type(conv_results_type), intent(inout)   :: Conv_results

      real, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2),   &
                             size(Input_mp%tin,3)) :: ttnd, utnd, vtnd, &
                                             det_cmt
      real, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2),   &
                             size(Input_mp%tin,3)+1) :: mc_cmt
      real, dimension(size(Output_mp%rdt,1), size(Output_mp%rdt,2),      &
                              size(Output_mp%rdt,3),num_tracers) :: qtr
      integer :: n
      integer :: k, kx
      integer :: im, jm, km, nq, nq_skip

      kx = size(Input_mp%tin,3)

!----------------------------------------------------------------------
!    if doing nonlocal cmt, call cu_mo_trans for each convective scheme
!    separately.
!----------------------------------------------------------------------
      if (do_nonlocal_transport) then
        im = size(Input_mp%uin,1)
        jm = size(Input_mp%uin,2)
        km = size(Input_mp%uin,3)
        nq = size(Input_mp%tracer,4)
        nq_skip = nq
        qtr    (:,:,:,1:nq_skip) = 0.0
 
        if (cmt_uses_uw) then

!----------------------------------------------------------------------
!    CURRENTLY no detrained mass flux is provided from uw_conv; should only
!    use with 'diffusive' cmt scheme, not the non-local. (attempt to
!    use non-local will cause FATAL in _init routine.)
!----------------------------------------------------------------------
        endif

      else ! (do_diffusive_transport)

!-----------------------------------------------------------------------
!    if using diffusive cmt, call cu_mo_trans once with combined mass
!    fluxes from all desired convective schemes.
!-----------------------------------------------------------------------
        mc_cmt = 0.
        det_cmt = 0.
        if (cmt_uses_uw) then
          do k=2,kx
            mc_cmt(:,:,k) = mc_cmt(:,:,k) + Conv_results%uw_mflux(:,:,k-1)
          end do
        endif

!------------------------------------------------------------------------
!    call cu_mo_trans to calculate cumulus momentum transport.
!------------------------------------------------------------------------
        call diffusive_cu_mo_trans (is, js, Time, mc_cmt, Input_mp%tin,      &
                                    Input_mp%phalf, Input_mp%pfull, Input_mp%zhalf,  &
                                                  Input_mp%zfull, Output_mp%diff_cu_mo)
      endif ! (do_nonlocal_transport)

!-----------------------------------------------------------------------


end subroutine cu_mo_trans


!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans_end">
!  <OVERVIEW>
!   terminates module
!  </OVERVIEW>
!  <DESCRIPTION>
!   This is the destructor for cu_mo_trans
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cu_mo_trans_end()

  module_is_initialized = .false.

  end subroutine cu_mo_trans_end

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans">
!  <OVERVIEW>
!   picks one of the available cumulus momentum transport parameteriz-
!    ations based on namelist-supplied information (currently diffusive 
!    or non-local options are available).  For the diffusive scheme, it
!    returns a diffusivity proportional to the convective mass 
!    flux. For the non-local scheme temperature, specific humidity and
!    momentum tendencies are returned.
!  </OVERVIEW>
!  <DESCRIPTION>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans (is, js, Time, mass_flux, t,           &
!                p_half, p_full, z_half, z_full, diff)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! ! horizontal domain on which computation to be performed is
!    (is:is+size(t,1)-1,ie+size(t,2)-1) in global coordinates
!   (used by diag_manager only)
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! current time, used by diag_manager
!  </IN>
!  <IN NAME="mass_flux" TYPE="real">
! convective mass flux (Kg/(m**2 s)), dimension(:,:,:), 3rd dimension is
!    vertical level (top down) -- defined at interfaces, so that
!    size(mass_flux,3) = size(p_half,3); entire field processed;
!   all remaining fields are 3 dimensional;
!  size of first two dimensions must confrom for all variables
!  </IN>
!  <IN NAME="t" TYPE="real">
! temperature (K) at full levels, size(t,3) = size(p_full,3)
!  </IN>
!  <IN NAME="p_half" TYPE="real">
! pressure at interfaces (Pascals) 
!  size(p_half,3) = size(p_full,3) + 1
!  </IN>
!  <IN NAME="p_full" TYPE="real">
! pressure at full levels (levels at which temperature is defined)
!  </IN>
!  <IN NAME="z_half" TYPE="real">
! height at half levels (meters); size(z_half,3) = size(p_half,3)
!  </IN>
!  <IN NAME="z_full" TYPE="real">
! height at full levels (meters); size(z_full,3) = size(p_full,3)
!  </IN>
!  <OUT NAME="diff" TYPE="real">
! kinematic diffusivity (m*2/s); defined at half levels 
!   size(diff,3) = size(p_half,3)
!  </OUT>
! </SUBROUTINE>
!



!#######################################################################

! <SUBROUTINE NAME="diffusive_cu_mo_trans">
!  <OVERVIEW>
!   returns a diffusivity proportional to the 
!    convective mass flux, for use with diffusive 
!    convective momentum transport closure
!  </OVERVIEW>
!  <DESCRIPTION>
!   A diffusive approximation to convective momentum transport is crude but
!    has been found to be useful in inproving the simulation of tropical
!    precipitation in some models.  The diffusivity computed here is
!    simply 
!<PRE>
! diffusivity = c*W*L
! W = M/rho  (m/sec)
! M = convective mass flux (kg/(m2 sec)) 
! rho - density of air
! L = depth of convecting layer (m)
! c = normalization constant = diff_norm/g
!   (diff_norm is a namelist parameter;
!      the factor of g here is an historical artifact)
! for further discussion see cu_mo_trans.ps
!</PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diffusive_cu_mo_trans (is, js, Time, mass_flux, t,      &
!                p_half, p_full, z_half, z_full, diff)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! ! horizontal domain on which computation to be performed is
!    (is:is+size(t,1)-1,ie+size(t,2)-1) in global coordinates
!   (used by diag_manager only)
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! current time, used by diag_manager
!  </IN>
!  <IN NAME="mass_flux" TYPE="real">
! convective mass flux (Kg/(m**2 s)), dimension(:,:,:), 3rd dimension is
!    vertical level (top down) -- defined at interfaces, so that
!    size(mass_flux,3) = size(p_half,3); entire field processed;
!   all remaining fields are 3 dimensional;
!  size of first two dimensions must confrom for all variables
!  </IN>
!  <IN NAME="t" TYPE="real">
! temperature (K) at full levels, size(t,3) = size(p_full,3)
!  </IN>
!  <IN NAME="p_half" TYPE="real">
! pressure at interfaces (Pascals) 
!  size(p_half,3) = size(p_full,3) + 1
!  </IN>
!  <IN NAME="p_full" TYPE="real">
! pressure at full levels (levels at which temperature is defined)
!  </IN>
!  <IN NAME="z_half" TYPE="real">
! height at half levels (meters); size(z_half,3) = size(p_half,3)
!  </IN>
!  <IN NAME="z_full" TYPE="real">
! height at full levels (meters); size(z_full,3) = size(p_full,3)
!  </IN>
!  <OUT NAME="diff" TYPE="real">
! kinematic diffusivity (m*2/s); defined at half levels 
!   size(diff,3) = size(p_half,3)
!  </OUT>
! </SUBROUTINE>
!
subroutine diffusive_cu_mo_trans (is, js, Time, mass_flux, t,     &
                        p_half, p_full, z_half, z_full, diff)

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: is, js

real, intent(in)   , dimension(:,:,:) :: mass_flux, t, &
                                         p_half, p_full, z_half, z_full
real, intent(out), dimension(:,:,:) :: diff                      

real, dimension(size(t,1),size(t,2),size(t,3)) :: rho
real, dimension(size(t,1),size(t,2))           :: zbot, ztop

integer :: k, nlev
logical :: used

!-----------------------------------------------------------------------

 if (.not.module_is_initialized) call error_mesg ('cu_mo_trans',  &
                      'cu_mo_trans_init has not been called.', FATAL)

!-----------------------------------------------------------------------

nlev = size(t,3)

zbot = z_half(:,:,nlev+1)
ztop = z_half(:,:,nlev+1)
  
do k = 2, nlev
  where(mass_flux(:,:,k) .ne. 0.0 .and. mass_flux(:,:,k+1) == 0.0) 
    zbot = z_half(:,:,k)
  endwhere
  where(mass_flux(:,:,k-1) == 0.0 .and. mass_flux(:,:,k) .ne. 0.0) 
    ztop = z_half(:,:,k)
  endwhere
end do

rho  = p_full/(RDGAS*t)  ! density 
   ! (including the virtual temperature effect here might give the 
   ! impression that this theory is accurate to 2%!)

! diffusivity = c*W*L
! W = M/rho  (m/sec)
! M = convective mass flux (kg/(m2 sec)) 
! L = ztop - zbot = depth of convecting layer (m)
! c = normalization constant = diff_norm/g
!   (the factor of g here is an historical artifact)

diff(:,:,1) = 0.0
do k = 2, nlev
  diff(:,:,k) = diff_norm*mass_flux(:,:,k)*(ztop-zbot)/(rho(:,:,k)*GRAV)
end do


! --- diagnostics
     if ( id_diff_cmt > 0 ) then
        used = send_data ( id_diff_cmt, diff, Time, is, js, 1 )
     endif
     if ( id_massflux_cmt > 0 ) then
        used = send_data ( id_massflux_cmt, mass_flux, Time, is, js, 1 )
     endif

end subroutine diffusive_cu_mo_trans

!#######################################################################

end module cu_mo_trans_mod


