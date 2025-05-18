module atmos_fire_emis_mod
  ! <DESCRIPTION>
  !    This module provides subroutines for passing interactive fire 
  !    emissions from land to the atmospher
  !    
  ! </DESCRIPTION>

use mpp_mod, only: mpp_pe, mpp_root_pe
use mpp_mod, only: input_nml_file

use constants_mod,   only: PI
use time_manager_mod, only : time_type, get_date, days_in_month, operator(-), &
                             time_type_to_real
use fms_mod, only : check_nml_error, error_mesg, stdlog, stdout, &
      lowercase, WARNING, FATAL, NOTE
use fms2_io_mod, only: close_file, FmsNetcdfFile_t, open_file
use diag_manager_mod, only : register_diag_field, send_data
use field_manager_mod , only : MODEL_ATMOS, MODEL_LAND, parse
use field_manager_mod, only: fm_field_name_len, fm_string_len, &
     fm_type_name_len, fm_path_name_len, fm_dump_list, fm_get_length, &
     fm_get_current_list, fm_loop_over_list, fm_change_list
use fm_util_mod, only : fm_util_get_real, fm_util_get_logical, fm_util_get_string, fm_util_get_real_array
use tracer_manager_mod, only : NO_TRACER, get_number_tracers, get_tracer_names, get_tracer_index, &
                               query_method
use mpp_mod, only: stdout, stdlog, mpp_error
use constants_mod, only: AVOGNO
use fms_mod, only : uppercase


implicit none
private

! ==== public interfaces =====================================================
public  ::  atmos_fire_emis_init, atmos_fire_emis_end
public :: atmos_fire_emis
public :: fire_emis_type

public    :: get_num_fire_tr  ! number of fire emission tracers

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'atmos_fire_emis'
logical :: used

logical         :: module_is_initialized =.FALSE.
integer, parameter :: MAX_FR_TR = 99
real            :: delta_time
real            :: dt_fast_yr      ! fast time step in years
integer :: id_fire_emis(MAX_FR_TR)
integer :: n_fire_tr  ! number of fire emission tracers

!--- Fire emissions type
type fire_emis_type
  character(fm_field_name_len):: name    = ''    ! name of the tracer
  integer                     :: tr_atm  = NO_TRACER ! index of this tracer in atmos tracer array
  real                        :: fire_mw = 1.0       ! molecular weights of fire tracers
  logical                     :: do_conversion = .false.
  real                        :: scale_factor  = 1.0 
end type
type(fire_emis_type), allocatable :: frdata(:) ! fire emissions data

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#####################################################################
! initialize vegetation tracers
subroutine atmos_fire_emis_init(axes, Time, fire_emis_ind, frdata)

 type(time_type),       intent(in) :: Time
 integer        , intent(in)           :: axes(4)
 integer        , intent(inout)           :: fire_emis_ind(:)
 type(fire_emis_type), allocatable, intent(inout) :: frdata(:)

 integer :: i, j, k, m, n, o, p, nsp, tr
!integer :: trind
 character(fm_field_name_len) :: name ! name of the vegn tracer
 character(fm_type_name_len)  :: typ  ! type of the vegn tracer
 integer :: nt_atmos
 character(len=32) :: ef_sp, unit_tr
 character(len=500)  :: method
 character(len=500) :: parameters
 real    :: value ! temporary storage for parsing input

! get number of atmos_tracers
  call get_number_tracers (MODEL_ATMOS, num_prog=nt_atmos)
 
  n_fire_tr = 0
! see if any of the atmos_tracers have bb_emis is lm4
  do tr = 1, nt_atmos
     call get_tracer_names (MODEL_ATMOS, tr, name = name)
!    trind = get_tracer_index(MODEL_ATMOS,name)
     if(query_method('emissions2dbb', MODEL_ATMOS, tr, method, parameters)) then
        if (trim(method)=='land:lm4') then
        n_fire_tr=n_fire_tr+1
        endif
     endif
  enddo 

  if (n_fire_tr > 0) then
     allocate(frdata(1:n_fire_tr))
     if (mpp_pe() == mpp_root_pe()) &
     write(*,*) 'Allocated frdata with size:', size(frdata)
  else
     if (mpp_pe() == mpp_root_pe()) &
     call mpp_error(WARNING, 'n_fire_tr is zero; cannot allocate frdata')
  endif

  i = 0
! register the frdata info
  do tr = 1, nt_atmos
     call get_tracer_names (MODEL_ATMOS, tr, name = name)
     method = ''; parameters = ''
     if(query_method('emissions2dbb', MODEL_ATMOS, tr, method, parameters)) then
        if (trim(method)=='land:lm4') then
            i = i + 1
            fire_emis_ind(tr) = i
            frdata(i)%name = trim(name)
            frdata(i)%tr_atm = tr ! get_tracer_index(MODEL_ATMOS,name)
            if (mpp_pe() == mpp_root_pe())  &
               write(*,*) 'atmos_fire_emis_init: emis2dbb tracer=', TRIM(name),tr,fire_emis_ind(tr),frdata(i)%tr_atm
            if ( parse(parameters, 'mw', value) > 0 ) then
                 frdata(i)%fire_mw = value
            endif
            if ( parse(parameters, 'scale_factor', value) > 0 ) then
                 frdata(i)%scale_factor = value
            endif
            method = ''; parameters = ''
            if(query_method('units', MODEL_ATMOS, tr, method, parameters)) then
               if (lowercase(trim(method))=='mmr') then
                  frdata(i)%do_conversion = .true.
               else
                  frdata(i)%do_conversion = .false.
               endif
            endif
        endif
     endif
  enddo 

  if (n_fire_tr .gt. MAX_FR_TR) call mpp_error(FATAL, 'Number of fire emission tracers defined exceeds the maximum of MAX_FR_TR -please increase MAX_FR_TR in vegn_data.F90')


  do i = 1,n_fire_tr

      if (frdata(i)%do_conversion) then
      unit_tr = 'kg/m2/s'
      else
      unit_tr = 'molecules/cm2/s'
      end if
    
     id_fire_emis(i) = register_diag_field ( module_name,             &
                    trim(frdata(i)%name)//'_fire_emis', axes(1:2), Time,              &
                    'Fire emissions (interactive)', trim(unit_tr) )

  enddo
  module_is_initialized = .TRUE.
end subroutine atmos_fire_emis_init

subroutine atmos_fire_emis_end
      if (allocated(frdata)) then
         deallocate(frdata)
      else
         if (mpp_pe() == mpp_root_pe()) &
         call mpp_error(WARNING, 'frdata is not allocated; cannot deallocate')
      endif
end subroutine atmos_fire_emis_end

! ============================================================================
! Added this subroutine to compute fire emissions for tracers within the atmos_tracer_driver subroutine
subroutine atmos_fire_emis(fire_emis, fire_emis_flux, frdata, diag_time, is, js)
   integer, intent(in)                    :: is, js
   real, intent(in),  dimension(:,:,:) :: fire_emis
   real, intent(inout),  dimension(:,:,:) :: fire_emis_flux
   type(fire_emis_type), allocatable, intent(inout) :: frdata(:)
   type(time_type), intent(in)            :: diag_time     
   integer :: i

   if (.not. allocated(frdata)) then
       if (mpp_pe() == mpp_root_pe()) &
       call mpp_error(FATAL, 'frdata is not allocated in atmos_fire_emis')
   endif

   do i = 1,n_fire_tr
      
      if (frdata(i)%do_conversion) then
      fire_emis_flux(:,:,i) = fire_emis(:,:,i)*1.e-3*1.e4*frdata(i)%fire_mw*(1/AVOGNO)*frdata(i)%scale_factor !molecules/cm2/s to kg/m2/s
      else
      fire_emis_flux(:,:,i) = fire_emis(:,:,i)*frdata(i)%scale_factor
      end if

     ! Mask out negative values
      fire_emis_flux(:,:,i) = max(fire_emis_flux(:,:,i), 0.0) 
      
      if (id_fire_emis(i) > 0) then
        used = send_data(id_fire_emis(i),fire_emis_flux(:,:,i), diag_time, &
              is_in=is,js_in=js)
      endif
   enddo

end subroutine atmos_fire_emis 

function get_num_fire_tr()
   integer :: get_num_fire_tr
   get_num_fire_tr = n_fire_tr
end function

end module atmos_fire_emis_mod 
