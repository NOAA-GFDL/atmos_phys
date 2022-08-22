!FDOC_TAG_GFDL
               module bulkphys_rad_mod
! <CONTACT EMAIL="fei.liu@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="">
!  
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    bulkphys_rad_mod defines cloud radiative properties based on
!    bulk cloud physics values in contrast to microphysically-based
!    properties.
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!
 
!    shared modules:

use mpp_mod,                only: input_nml_file
use fms_mod,                only: mpp_pe, &
                                  fms_init, mpp_root_pe, stdlog,  &
                                  write_version_number, & 
                                  check_nml_error, error_mesg,   &
                                  FATAL

!    individual cloud modules:

use cloudrad_types_mod, only: cld_specification_type, &
                              cldrad_properties_type, &
                              cloudrad_control_type

use strat_clouds_W_mod,     only: strat_clouds_W_init

!BW use rh_based_clouds_mod,    only: rh_based_clouds_init, &
!BW                                   obtain_bulk_sw_rh,   &
!BW                                   obtain_bulk_lw_rh

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    bulkphys_rad_mod defines cloud radiative properties based on
!    bulk cloud physics values in contrast to microphysically-based
!    properties.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'



!---------------------------------------------------------------------
!-------  interfaces --------

public                                            &
          bulkphys_rad_init,                      &
          bulkphys_rad_end

!---------------------------------------------------------------------
!-------- namelist  ---------
 
integer     :: dummy = 1

namelist /bulkphys_rad_nml /  dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!    visible band reflectivity, nir band reflectivity and nir absorp-
!    tivities are given for high, middle and low clouds. two separate
!    sets of parameters are present; set 1 has been used for all cloud
!    parameterizations except for mgrp_prscr_clds, which used the
!    values in data set 2.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   crfvis: visible band reflectivity. 
!--------------------------------------------------------------------
real  :: crfvis_hi_1  = 0.21
real  :: crfvis_hi_2  = 0.21
real  :: crfvis_mid_1 = 0.45
real  :: crfvis_mid_2 = 0.48
real  :: crfvis_low_1 = 0.59
real  :: crfvis_low_2 = 0.69

!--------------------------------------------------------------------
!   crfir: near-ir band reflectivity.
!--------------------------------------------------------------------
real  :: crfir_hi_1   = 0.21
real  :: crfir_hi_2   = 0.21
real  :: crfir_mid_1  = 0.45
real  :: crfir_mid_2  = 0.48
real  :: crfir_low_1  = 0.59
real  :: crfir_low_2  = 0.69

!--------------------------------------------------------------------
!   cabir: near-ir band absorptivity.
!--------------------------------------------------------------------
real  :: cabir_hi_1   = 0.005
real  :: cabir_hi_2   = 0.005
real  :: cabir_mid_1  = 0.02
real  :: cabir_mid_2  = 0.02
real  :: cabir_low_1  = 0.035
real  :: cabir_low_2  = 0.035

!--------------------------------------------------------------------
!   cldem:  infrared emissivity.
!--------------------------------------------------------------------
real  :: cldem_hi   = 1.00
real  :: cldem_mid  = 1.00
real  :: cldem_low  = 1.00

!------------------------------------------------------------------
!    these variables hold the values that are to be used for the
!    cloud radiative properties for the cloud parameterization 
!    activated.
!------------------------------------------------------------------
real  :: crfvis_hi, crfvis_mid, crfvis_low,    &
         crfir_hi,  crfir_mid,  crfir_low,  &
         cabir_hi,  cabir_mid,  cabir_low

real  :: min_cld_drop_rad, max_cld_drop_rad, &
         min_cld_ice_size, max_cld_ice_size

!-------------------------------------------------------------------
!    logical flag.
!-------------------------------------------------------------------
logical  :: module_is_initialized = .false.

!-------------------------------------------------------------------
!-------------------------------------------------------------------

                           contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!#####################################################################

! <SUBROUTINE NAME="bulkphys_rad_init">
!  <OVERVIEW>
!    subroutine bulkphys_rad_init is the constructor for
!    bulkphys_rad_mod.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine bulkphys_rad_init is the constructor for
!    bulkphys_rad_mod.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call bulkphys_rad_init (pref, lonb, latb)
!
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!       pref      array containing two reference pressure profiles
!                 for use in defining transmission functions [ Pa ]
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!       lonb      2d array of model longitudes on cell corners
!                 [ radians ]
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
!       latb      2d array of model latitudes at cell corners [radians]
! 
!  </IN>
! </SUBROUTINE>
!
subroutine bulkphys_rad_init (min_cld_drop_rad_in, max_cld_drop_rad_in,&
                              min_cld_ice_size_in, max_cld_ice_size_in,&
                              pref, lonb, latb, Cldrad_control)

!---------------------------------------------------------------------
!    subroutine bulkphys_rad_init is the constructor for 
!    bulkphys_rad_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real,                        intent(in) :: min_cld_drop_rad_in, &
                                           max_cld_drop_rad_in,&
                                           min_cld_ice_size_in,  &
                                           max_cld_ice_size_in
real, dimension(:,:),        intent(in) :: pref
real, dimension(:,:),        intent(in) :: lonb, latb
type(cloudrad_control_type), intent(in) :: Cldrad_control

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      2d array of model longitudes on cell corners [ radians ]
!       latb      2d array of model latitudes at cell corners [radians]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      integer  ::   ierr, io, logunit
      integer  ::   idum

!---------------------------------------------------------------------
!    local variables:
!
!      ierr     error code
!      io       error status returned from io operation  
!      idum     dummy integer argument needed to satisfy 
!               diag_clouds_W_init call
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      if (Cldrad_control%do_strat_clouds) call strat_clouds_W_init(latb, lonb, Cldrad_control)
!BW   if (Cldrad_control%do_rh_clouds)    call rh_based_clouds_init 

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      read (input_nml_file, nml=bulkphys_rad_nml, iostat=io)
      ierr = check_nml_error(io,'bulkphys_rad_nml')
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                      write (logunit, nml=bulkphys_rad_nml)

      min_cld_drop_rad = min_cld_drop_rad_in
      max_cld_drop_rad = max_cld_drop_rad_in
      min_cld_ice_size = min_cld_ice_size_in
      max_cld_ice_size = max_cld_ice_size_in

!--------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------


end subroutine bulkphys_rad_init


!###################################################################
 
! <SUBROUTINE NAME="bulkphys_rad_end">
!  <OVERVIEW>
!    bulkphys_rad_end is the destructor for bulkphys_rad_mod.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    bulkphys_rad_end is the destructor for bulkphys_rad_mod.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call bulkphys_rad_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine bulkphys_rad_end

!-------------------------------------------------------------------
!    bulkphys_rad_end is the destructor for bulkphys_rad_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
     module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine bulkphys_rad_end


!###################################################################

                     end module bulkphys_rad_mod

 
