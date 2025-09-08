module matrix_actv
    use fms_mod,               only : file_exist, close_file,&
        open_namelist_file, check_nml_error, &
        write_version_number, &
        error_mesg, &
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

        public getactfrac

!#include "rundeck_opts.h"
!      module aero_actv
!      use aero_param,  only: nlays, aunit1
!      use aero_config, only: nmodes
!!-------------------------------------------------------------------------------------------------------------------------
!!@auth    susanne bauer/doug wright
!!
!!-------------------------------------------------------------------------------------------------------------------------
!      
      real(8), parameter :: dens_sulf = 1.77d+03    ! [kg/m^3] nh42so4
      real(8), parameter :: dens_bcar = 1.70d+03    ! [kg/m^3] ghan et al. (2001) - mirage
      real(8), parameter :: dens_ocar = 1.00d+03    ! [kg/m^3] ghan et al. (2001) - mirage
      real(8), parameter :: dens_dust = 2.60d+03    ! [kg/m^3] ghan et al. (2001) - mirage
      real(8), parameter :: dens_seas = 2.165d+03   ! [kg/m^3] nacl, ghan et al. (2001) used 1.90d+03
!#ifdef tracers_amp_m9
!      real(8), parameter :: dens_ocm2 = 1.00d+03
!      real(8), parameter :: dens_ocm1 = 1.00d+03
!      real(8), parameter :: dens_ocm0 = 1.00d+03
!      real(8), parameter :: dens_ocp1 = 1.00d+03
!      real(8), parameter :: dens_ocp2 = 1.00d+03
!      real(8), parameter :: dens_ocp3 = 1.00d+03
!      real(8), parameter :: dens_ocp4 = 1.00d+03
!      real(8), parameter :: dens_ocp5 = 1.00d+03
!      real(8), parameter :: dens_ocp6 = 1.00d+03
!#endif
!      real(8) :: nactiv(nmodes)       ! for use in other subroutines

      contains


      subroutine getactfrac(nmodex,xnap,xmap5,rg,sigmag,tkelvin,ptot,wupdraft, nact, mact)
                           !ac,fracactn,fracactm,nact,mact)
!----------------------------------------------------------------------------------------------------------------------
!     12-12-06, dlw: routine to set up the call to subr. actfrac_mat to calculate the 
!                    activated fraction of the number and mass concentrations, 
!                    as well as the number and mass concentrations activated 
!                    for each of nmodex modes. the minimum dry radius for activation 
!                    for each mode is also returned. 
!
!     each mode is assumed to potentially contains 5 chemical species:
!         (1) sulfate 
!         (2) bc 
!         (3) oc
!         (4) mineral dust
!         (5) sea salt 
!
!     the aerosol activation parameterizations are described in 
!
!         1. abdul-razzak et al.   1998, jgr, vol.103, p.6123-6131.
!         2. abdul-razzak and ghan 2000, jgr, vol.105, p.6837-6844. 
!
!     and values for many of the required parameters were taken from 
!
!         3. ghan et al. 2001, jgr vol 106, p.5295-5316.
!
!     with the density of sea salt set to the value used in ref. 3 (1900 kg/m^3), this routine 
!     yields values for the hygroscopicity parameters bi in agreement with ref. 3. 
!----------------------------------------------------------------------------------------------------------------------
      implicit none
!
!#ifdef tracers_amp_m9
!      integer, parameter :: ncomps = 14
!#else
      integer, parameter :: ncomps = 5
!#endif

      ! arguments.
      
      integer :: nmodex               ! number of modes [1]/population (npop)      
      real(8), intent(in) :: xnap(nmodex)         ! number concentration for each mode [#/m^3]
      real(8), intent(in) :: xmap5(nmodex,ncomps) ! mass concentration of each of the 5 species for each mode [ug/m^3]
      real(8), intent(in) :: rg(nmodex)           ! geometric mean dry radius for each mode [um]
      real(8), intent(in) :: sigmag(nmodex)       ! geometric standard deviation for each mode [um]
      real(8) :: tkelvin              ! absolute temperature [k]
      real(8) :: ptot                 ! ambient pressure [pa]
      real(8) :: wupdraft             ! updraft velocity [m/s]
      real(8) :: ac(nmodex)           ! minimum dry radius for activation for each mode [um]
      real(8) :: fracactn(nmodex)     ! activating fraction of number conc. for each mode [1]
      real(8) :: fracactm(nmodex)     ! activating fraction of mass   conc. for each mode [1]
      real(8), intent(out) :: nact(nmodex)         ! activating number concentration for each mode [#/m^3]
      real(8) :: mact(nmodex)         ! activating mass   concentration for each mode [ug/m^3]
      ! local variables. 
      integer :: i, j                 ! loop counters 
      real(8) :: xmap(nmodex)         ! total mass concentration for each mode [ug/m^3]
      real(8) :: bibar(nmodex)        ! hygroscopicity parameter for each mode [1]

      real(8) :: sumnumer, sumdenom         ! scratch variables 

      real(8), parameter :: nion_sulf = 3.00d+00    ! [1]
      real(8), parameter :: nion_bcar = 1.00d+00    ! [1]
      real(8), parameter :: nion_ocar = 1.00d+00    ! [1]
      real(8), parameter :: nion_dust = 2.30d+00    ! [1]
      real(8), parameter :: nion_seas = 2.00d+00    ! [1] nacl
!#ifdef tracers_amp_m9
!      real(8), parameter :: nion_ocm2 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocm1 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocm0 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocp1 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocp2 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocp3 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocp4 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocp5 = 1.00d+00    ! [1]
!      real(8), parameter :: nion_ocp6 = 1.00d+00    ! [1]
!#endif

      real(8), parameter :: xphi_sulf = 0.70d+00    ! [1], osmotic coefficient
      real(8), parameter :: xphi_bcar = 1.00d+00    ! [1]
      real(8), parameter :: xphi_ocar = 1.00d+00    ! [1]
      real(8), parameter :: xphi_dust = 1.00d+00    ! [1]
      real(8), parameter :: xphi_seas = 1.00d+00    ! [1] nacl
!#ifdef tracers_amp_m9
!      real(8), parameter :: xphi_ocm2 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocm1 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocm0 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocp1 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocp2 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocp3 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocp4 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocp5 = 1.00d+00    ! [1]
!      real(8), parameter :: xphi_ocp6 = 1.00d+00    ! [1]
!#endif

      real(8), parameter :: molw_sulf = 132.0d-03   ! [kg/mol]
      real(8), parameter :: molw_bcar = 100.0d-03   ! [kg/mol]
      real(8), parameter :: molw_ocar = 100.0d-03   ! [kg/mol]
      real(8), parameter :: molw_dust = 100.0d-03   ! [kg/mol]
      real(8), parameter :: molw_seas = 58.44d-03   ! [kg/m^3] nacl
!#ifdef tracers_amp_m9
!      real(8), parameter :: molw_ocm2 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocm1 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocm0 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocp1 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocp2 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocp3 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocp4 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocp5 = 100.0d-03   ! [kg/mol]
!      real(8), parameter :: molw_ocp6 = 100.0d-03   ! [kg/mol]
!#endif

      real(8), parameter :: xeps_sulf = 1.00d+00    ! [1], soluble fraction
      real(8), parameter :: xeps_bcar = 1.67d-06    ! [1]
      real(8), parameter :: xeps_ocar = 0.78d+00    ! [1]
      real(8), parameter :: xeps_dust = 0.13d+00    ! [1]
      real(8), parameter :: xeps_seas = 1.00d+00    ! [1] nacl
!#ifdef tracers_amp_m9
!      real(8), parameter :: xeps_ocm2 = 1.d+00      ! [1]
!      real(8), parameter :: xeps_ocm1 = 0.875d+00   ! [1]
!      real(8), parameter :: xeps_ocm0 = 0.75d+00    ! [1]
!      real(8), parameter :: xeps_ocp1 = 0.625d+00   ! [1]
!      real(8), parameter :: xeps_ocp2 = 0.5+00      ! [1]
!      real(8), parameter :: xeps_ocp3 = 0.375d+00   ! [1]
!      real(8), parameter :: xeps_ocp4 = 0.25d+00    ! [1]
!      real(8), parameter :: xeps_ocp5 = 0.125d+00   ! [1]
!      real(8), parameter :: xeps_ocp6 = 0.d+00      ! [1]
!#endif

      real(8), parameter :: wmolmass = 18.01528d-03 ! molar mass of h2o     [kg/mol]
      real(8), parameter :: denh2o   =  1.00d+03    ! density of water [kg/m^3]

      ! variables for mode-average hygroscopicity parameters.
      real(8)       :: xr  (nmodex,ncomps)  ! mass fraction for component j in mode i [1]

!#ifdef tracers_amp_m9
!      ! # of ions formed per formula unit solute for component j in mode i [1]
!      real(8), dimension(ncomps), parameter :: xnu=(/nion_sulf,nion_bcar,
!     &                                               nion_ocar,nion_dust,
!     &                                               nion_seas,nion_ocm2,
!     &                                               nion_ocm1,nion_ocm0,
!     &                                               nion_ocp1,nion_ocp2,
!     &                                               nion_ocp3,nion_ocp4,
!     &                                               nion_ocp5,nion_ocp6/)
!      ! osmotic coefficient for component j in mode i [1]
!      real(8), dimension(ncomps), parameter :: xphi=(/xphi_sulf,xphi_bcar,
!     &                                                xphi_ocar,xphi_dust,
!     &                                                xphi_seas,xphi_ocm2,
!     &                                                xphi_ocm1,xphi_ocm0,
!     &                                                xphi_ocp1,xphi_ocp2,
!     &                                                xphi_ocp3,xphi_ocp4,
!     &                                                xphi_ocp5,xphi_ocp6/)
!      ! density of component j in mode i [kg/m^3]
!      real(8), dimension(ncomps), parameter :: xrho=(/dens_sulf,dens_bcar,
!     &                                                dens_ocar,dens_dust,
!     &                                                dens_seas,dens_ocm2,
!     &                                                dens_ocm1,dens_ocm0,
!     &                                                dens_ocp1,dens_ocp2,
!     &                                                dens_ocp3,dens_ocp4,
!     &                                                dens_ocp5,dens_ocp6/)
!      ! soluble fraction of component j in mode i [1]
!      real(8), dimension(ncomps), parameter :: xeps=(/xeps_sulf,xeps_bcar,
!     &                                                xeps_ocar,xeps_dust,
!     &                                                xeps_seas,xeps_ocm2,
!     &                                                xeps_ocm1,xeps_ocm0,
!     &                                                xeps_ocp1,xeps_ocp2,
!     &                                                xeps_ocp3,xeps_ocp4,
!     &                                                xeps_ocp5,xeps_ocp6/)
!      ! molecular weight for component j in mode i [kg/mol]
!      real(8), dimension(ncomps), parameter :: xmw=(/molw_sulf,molw_bcar,
!     &                                               molw_ocar,molw_dust,
!     &                                               molw_seas,molw_ocm2,
!     &                                               molw_ocm1,molw_ocm0,
!     &                                               molw_ocp1,molw_ocp2,
!     &                                               molw_ocp3,molw_ocp4,
!     &                                               molw_ocp5,molw_ocp6/)
!#else
      ! # of ions formed per formula unit solute for component j in mode i [1]
      real(8), dimension(ncomps), parameter :: xnu=(/nion_sulf,nion_bcar, &
                                                    nion_ocar,nion_dust, &
                                                    nion_seas/)
      ! osmotic coefficient for component j in mode i [1]
      real(8), dimension(ncomps), parameter :: xphi=(/xphi_sulf,xphi_bcar,&
                                                     xphi_ocar,xphi_dust, &
                                                     xphi_seas/)
      ! density of component j in mode i [kg/m^3]
      real(8), dimension(ncomps), parameter :: xrho=(/dens_sulf,dens_bcar, &
                                                     dens_ocar,dens_dust, &
                                                     dens_seas/)
      ! soluble fraction of component j in mode i [1]
      real(8), dimension(ncomps), parameter :: xeps=(/xeps_sulf,xeps_bcar, &
                                                     xeps_ocar,xeps_dust, &
                                                     xeps_seas/)
      ! molecular weight for component j in mode i [kg/mol]
      real(8), dimension(ncomps), parameter :: xmw=(/molw_sulf,molw_bcar, &
                                                    molw_ocar,molw_dust, &
                                                    molw_seas/)
!#endif

      !--------------------------------------------------------------------------------------------------------------
      ! calculate the mass fraction component j for each mode i. 
      !--------------------------------------------------------------------------------------------------------------
      do i=1, nmodex
        xmap(i) = 0.0d+00
        do j=1, ncomps
                  xmap(i) = xmap(i) + xmap5(i,j)
        enddo
        xr(i,:) = xmap5(i,:) / max( xmap(i), 1.0d-30 )   
        !write(*,'(i4,5f12.6)') i,xr(i,:)
      enddo

      !--------------------------------------------------------------------------------------------------------------
      ! calculate the hygroscopicity parameter for each mode. 
      !--------------------------------------------------------------------------------------------------------------
      do i=1, nmodex
        sumnumer = 0.0d+00
        sumdenom = 0.0d+00
        do j=1, ncomps
          sumnumer = sumnumer + xr(i,j)*xnu(j)*xphi(j)*xeps(j)/xmw(j)     ! [mol/kg] 
          sumdenom = sumdenom + xr(i,j)/xrho(j)                           ! [m^3/kg] 
        enddo
        !write(*,*) 'i,xr(i,:)=', i,xr(i,:), 'sumnumer, sumdenom = ', sumnumer, sumdenom
        bibar(i) = ( wmolmass*sumnumer ) / max( denh2o*sumdenom, 1.0d-30 )            ! [1] xl for zero treatment 
      enddo

      ! write(*,'(8d15.6)') bibar(:)

      !--------------------------------------------------------------------------------------------------------------
      ! calculate the droplet activation parameters for each mode. 
      !--------------------------------------------------------------------------------------------------------------
!      write(mpp_pe()+100, *) "actfrac_mat_before", "nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact", nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact
!       if (mpp_root_pe().eq.mpp_pe()) then
!                 write(*,*) "actfrac_mat_before", "nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact", nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact
!      endif
              call actfrac_mat(nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
                      ac,fracactn,fracactm,nact,mact)
!      if (mpp_root_pe().eq.mpp_pe()) then
!                                 write(*,*) "actfrac_mat_after", "nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact", nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact
!      endif
!        write(mpp_pe()+100, *) "actfrac_mat_after", "nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact", nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
!                      ac,fracactn,fracactm,nact,mact
        !write(*,*) "actfrac_mat_after", "nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
        !              ac,fracactn,fracactm,nact,mact", nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
        !              ac,fracactn,fracactm,nact,mact
      
      do i=1, nmodex
        if(xnap(i) .lt. 1.0d-06 ) fracactn(i) = 1.0d-30
      enddo

      end subroutine getactfrac


      subroutine actfrac_mat(nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
                            ac,fracactn,fracactm,nact,mact)
!----------------------------------------------------------------------------------------------------------------------
!     12-12-06, dlw: routine to calculate the activated fraction of the number 
!                    and mass concentrations, as well as the number and mass 
!                    concentrations activated for each of nmodex modes. the 
!                    minimum dry radius for activation for each mode is also returned. 
!
!     the aerosol activation parameterizations are described in 
!
!         1. abdul-razzak et al.   1998, jgr, vol.103, p.6123-6131.
!         2. abdul-razzak and ghan 2000, jgr, vol.105, p.6837-6844. 
! 
!     this routine is for the multiple-aerosol type parameterization. 
!----------------------------------------------------------------------------------------------------------------------
!      use domain_decomp_atm,only: am_i_root
      implicit none

      ! arguments.
      
      integer, intent(in) :: nmodex            ! number of modes [1]      
      real(8), intent(in) :: xnap(nmodex)      ! number concentration for each mode [#/m^3]
      real(8), intent(in) :: xmap(nmodex)      ! mass   concentration for each mode [ug/m^3]
      real(8), intent(in) :: rg(nmodex)        ! geometric mean radius for each mode [um]
      real(8), intent(in) :: sigmag(nmodex)    ! geometric standard deviation for each mode [um]
      real(8), intent(in) :: bibar(nmodex)     ! hygroscopicity parameter for each mode [1]
      real(8), intent(in) :: tkelvin           ! absolute temperature [k]
      real(8), intent(in) :: ptot              ! ambient pressure [pa]
      real(8), intent(in) :: wupdraft          ! updraft velocity [m/s]
      real(8), intent(out) :: ac(nmodex)        ! minimum dry radius for activation for each mode [um]
      real(8) :: ac_2(nmodex)        ! minimum dry radius for activation for each mode [um]
      real(8) :: ac_3(nmodex)        ! minimum dry radius for activation for each mode [um]
      real(8) :: ac_5(nmodex)        ! minimum dry radius for activation for each mode [um]
      real(8), intent(out) :: fracactn(nmodex)  ! activating fraction of number conc. for each mode [1]
      real(8) :: fracactn_2(nmodex)  ! activating fraction of number conc. for each mode [1]
      real(8) :: fracactn_3(nmodex)  ! activating fraction of number conc. for each mode [1]
      real(8) :: fracactn_5(nmodex)  ! activating fraction of number conc. for each mode [1]
      real(8), intent(out) :: fracactm(nmodex)  ! activating fraction of mass   conc. for each mode [1]
      real(8), intent(out) :: nact(nmodex)      ! activating number concentration for each mode [#/m^3]
      real(8), intent(out) :: mact(nmodex)      ! activating mass   concentration for each mode [ug/m^3]

      ! parameters.
      
      real(8), parameter :: pi            = 3.141592653589793d+00
      real(8), parameter :: twopi         = 2.0d+00 * pi
      real(8), parameter :: sqrt2         = 1.414213562d+00
      real(8), parameter :: threesqrt2by2 = 1.5d+00 * sqrt2

      real(8), parameter :: avgnum   = 6.0221367d+23       ! [1/mol]
      real(8), parameter :: rgasjmol = 8.31451d+00         ! [j/mol/k]
      real(8), parameter :: wmolmass = 18.01528d-03        ! molar mass of h2o     [kg/mol]
      real(8), parameter :: amolmass = 28.966d-03          ! molar mass of air     [kg/mol]
      real(8), parameter :: asmolmss = 132.1406d-03        ! molar mass of nh42so4 [kg/mol]
      real(8), parameter :: denh2o   = 1.00d+03            ! density of water [kg/m^3]
      real(8), parameter :: denamsul = 1.77d+03            ! density of pure ammonium sulfate [kg/m^3]
      real(8), parameter :: xnuamsul = 3.00d+00            ! # of ions formed when the salt is dissolved in water [1]
      real(8), parameter :: phiamsul = 1.000d+00           ! osmotic coefficient value in a-r 1998. [1] 
      real(8), parameter :: gravity  = 9.81d+00            ! grav. accel. at the earth's surface [m/s/s] 
      real(8), parameter :: heatvap  = 40.66d+03/wmolmass  ! latent heat of vap. for water and tnbp [j/kg] 
      real(8), parameter :: cpair    = 1006.0d+00          ! heat capacity of air [j/kg/k] 
      real(8), parameter :: t0dij    = 273.15d+00          ! reference temp. for dv [k] 
      real(8), parameter :: p0dij    = 101325.0d+00        ! reference pressure for dv [pa] 
      real(8), parameter :: dijh2o0  = 0.211d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------    
      ! real(8), parameter :: t0dij    = 283.15d+00          ! reference temp. for dv [k] 
      ! real(8), parameter :: p0dij    = 80000.0d+00         ! reference pressure for dv [pa] 
      ! real(8), parameter :: dijh2o0  = 0.300d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------
      real(8), parameter :: deltav   = 1.096d-07           ! vapor jump length [m]  
      real(8), parameter :: deltat   = 2.160d-07           ! thermal jump length [m]  
      real(8), parameter :: alphac   = 1.000d+00           ! condensation mass accommodation coefficient [1]  
      real(8), parameter :: alphat   = 0.960d+00           ! thermal accommodation coefficient [1]  

      ! local variables. 

      integer            :: i                              ! loop counter 
      real(8)            :: dv                             ! diffusion coefficient for water [m^2/s] 
      real(8)            :: dvprime                        ! modified diffusion coefficient for water [m^2/s] 
      real(8)            :: dumw, duma                     ! scratch variables [s/m] 
      real(8)            :: wpe                            ! saturation vapor pressure of water [pa]  
      real(8)            :: surten                         ! surface tension of air-water interface [j/m^2] 
      real(8)            :: xka                            ! thermal conductivity of air [j/m/s/k]  
      real(8)            :: xkaprime                       ! modified thermal conductivity of air [j/m/s/k]  
      real(8)            :: eta(nmodex)                    ! model parameter [1]  
      real(8)            :: zeta                           ! model parameter [1]  
      real(8)            :: xlogsigm(nmodex)               ! ln(sigmag) [1]   
      real(8)            :: a                              ! [m]
      real(8)            :: g                              ! [m^2/s]   
      real(8)            :: rdrp                           ! [m]   
      real(8)            :: f1                             ! [1]   
      real(8)            :: f2                             ! [1]
      real(8)            :: alpha                          ! [1/m]
      real(8)            :: gammav                          ! [m^3/kg]   
      real(8)            :: sm(nmodex)                     ! [1]   
      real(8)            :: dum                            ! [1/m]    
      real(8)            :: u                              ! argument to error function [1]
      real(8)            :: smax                           ! maximum supersaturation [1]
      real :: tmp, tmp2
      real :: bibar_i
      tmp = 0
     tmp2 = 0 
!----------------------------------------------------------------------------------------------------------------------
!     rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
!     values close to those given in a-z et al. 1998 figure 5. 
!----------------------------------------------------------------------------------------------------------------------
      rdrp = 0.105d-06   ! [m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.  
!----------------------------------------------------------------------------------------------------------------------
!     these variables are common to all modes and need only be computed once. 
!----------------------------------------------------------------------------------------------------------------------
      dv = dijh2o0*(p0dij/ptot)*(tkelvin/t0dij)**1.94d+00                 ! [m^2/s] (p&k,2nd ed., p.503)
      surten = 76.10d-03 - 0.155d-03 * (tkelvin-273.15d+00)               ! [j/m^2] 
      wpe = exp( 77.34491296d+00 - 7235.424651d+00/tkelvin - 8.2d+00*log(tkelvin) + tkelvin*5.7113d-03 )  ! [pa] 
      dumw = sqrt(twopi*wmolmass/rgasjmol/tkelvin)                        ! [s/m] 
      dvprime = dv / ( (rdrp/(rdrp+deltav)) + (dv*dumw/(rdrp*alphac)) )   ! [m^2/s] - eq. (17) 
      xka = (5.69d+00+0.017d+00*(tkelvin-273.15d+00))*418.4d-05           ! [j/m/s/k] (0.0238 j/m/s/k at 273.15 k)
      duma = sqrt(twopi*amolmass/rgasjmol/tkelvin)                        ! [s/m]
      xkaprime = xka / ( ( rdrp/(rdrp+deltat) ) + ( xka*duma/(rdrp*alphat*denh2o*cpair) ) )   ! [j/m/s/k]
      g = 1.0d+00 / ( (denh2o*rgasjmol*tkelvin) / (wpe*dvprime*wmolmass) &
                     + ( (heatvap*denh2o) / (xkaprime*tkelvin) ) &
                     * ( (heatvap*wmolmass) / (rgasjmol*tkelvin) - 1.0d+00 ) )               ! [m^2/s]
      a = (2.0d+00*surten*wmolmass)/(denh2o*rgasjmol*tkelvin)                                 ! [m] 
      alpha = (gravity/(rgasjmol*tkelvin))*((wmolmass*heatvap)/(cpair*tkelvin) - amolmass)    ! [1/m] 
      gammav = (rgasjmol*tkelvin)/(wpe*wmolmass) &
           + (wmolmass*heatvap*heatvap)/(cpair*ptot*amolmass*tkelvin)                        ! [m^3/kg]
     
             !write(*,*) "matrix_activation"
             !write(*,*) "dv=", dv
             !write(*,*) "surten=",surten
             !write(*,*) "wpe=", wpe
             !write(*,*) "dumw=",dumw
             !write(*,*) "dvprime=",dvprime
             !write(*,*) "xka=", xka
             !write(*,*) "duma=",duma
             !write(*,*) "xkaprime=",xkaprime
             !!write(*,*) "alpha=",alpha
             !write(*,*) "wupdraft=",wupdraft
             !write(*,*) "g=", g
      dum = sqrt(alpha*wupdraft/g)                  ! [1/m] 
      zeta = 2.d+00*a*dum/3.d+00                    ! [1] 
      !----------------------------------------------------------------------------------------------------------------
      ! write(1,'(a27,4d15.5)')'surten,wpe,a            =',surten,wpe,a
      ! write(1,'(a27,4d15.5)')'xka,xkaprime,dv,dvprime =',xka,xkaprime,dv,dvprime
      ! write(1,'(a27,4d15.5)')'alpha,gammav,g, zeta     =',alpha,gammav,g,zeta
!----------------------------------------------------------------------------------------------------------------------
!     these variables must be computed for each mode. 
!----------------------------------------------------------------------------------------------------------------------
      xlogsigm(:) = log(sigmag(:))                                                    ! [1] 
      smax = 0.0d+00                                                                  ! [1]
      do i=1, nmodex
        if (bibar(i) .eq. 0) then
                write(*,*) "i, bibar(i) = ", i, bibar(i)
        endif
        if (rg(i) .eq. 0) then
                write(*,*) "i, rg(i) =", i, rg(i)
        endif
        bibar_i = max(bibar(i), 1.0d-30)
        sm(i) = ( 2.0d+00/sqrt(bibar_i) ) * ( a/(3.0d-06*rg(i)) )**1.5d+00           ! [1] 
        if (gammav .eq. 0) then
                write(*,*) "gammav = ", gammav
        endif
        if (xnap(i) .eq. 0) then
                write(*,*) "i, xnap(i) = ", i, xnap(i)
        endif
        eta(i) = dum**3 / (twopi*denh2o*gammav*xnap(i))                                ! [1] 
        !--------------------------------------------------------------------------------------------------------------
        ! write(1,'(a27,i4,4d15.5)')'i,eta(i),sm(i) =',i,eta(i),sm(i)
        !--------------------------------------------------------------------------------------------------------------
        f1 = 0.5d+00 * exp(2.50d+00 * xlogsigm(i)**2)                                 ! [1] 
        f2 = 1.0d+00 +     0.25d+00 * xlogsigm(i)                                     ! [1] 
        smax = smax + (   f1*(  zeta  / eta(i)              )**1.50d+00 & 
                       + f2*(sm(i)**2/(eta(i)+3.0d+00*zeta))**0.75d+00 ) / sm(i)**2  ! [1] - eq. (6)
      enddo 
      smax = 1.0d+00 / sqrt(smax)                                                     ! [1]
      do i=1, nmodex
        ac(i)       = rg(i) * ( sm(i) / smax )**0.66666666666666667d+00               ! [um]
        u           = log(ac(i)/rg(i)) / ( sqrt2 * xlogsigm(i) )                      ! [1]
        call erff(u, tmp)
        fracactn(i) = 0.5d+00 * (1.0d+00 - tmp)        ! [1]
        call erff(u - threesqrt2by2*xlogsigm(i), tmp2)
        fracactm(i) = 0.5d+00 * (1.0d+00 - tmp2 )      ! [1]
        nact(i)     = fracactn(i) * xnap(i)                                           ! [#/m^3]
        mact(i)     = fracactm(i) * xmap(i)                                           ! [ug/m^3]
        !--------------------------------------------------------------------------------------------------------------
      enddo 

      end subroutine actfrac_mat


      subroutine gcf(gammcf,a,x,gln)

      implicit none
!-----------------------------------------------------------------------------------------------------------------------
!     see numerical recipes, w. press et al., 2nd edition.
!-----------------------------------------------------------------------------------------------------------------------
      integer, parameter :: itmax=10000
      real(8), parameter :: eps=3.0d-07
      real(8), parameter :: fpmin=1.0d-30
      real(8) :: a,gammcf,gln,x
      integer :: i
      real(8) :: an,b,c,d,del,h
      gln=gammln(a)
      b=x+1.0d+00-a
      c=1.0d+00/fpmin
      d=1.0d+00/b
      h=d
      do i=1,itmax
        an=-i*(i-a)
        b=b+2.0d+00
        d=an*d+b
        if(abs(d).lt.fpmin)d=fpmin
        c=b+an/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1.0d+00/d
        del=d*c
        h=h*del
        if(abs(del-1.0d+00).lt.eps)goto 1
      enddo
      write(*,*)'aero_actv: subroutine gcf: a too large, itmax too small', gammcf,a,x,gln
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      end subroutine gcf


      subroutine gser(gamser,a,x,gln)

      implicit none
!-----------------------------------------------------------------------------------------------------------------------
!     see numerical recipes, w. press et al., 2nd edition.
!-----------------------------------------------------------------------------------------------------------------------
      integer, parameter :: itmax=10000  ! was itmax=100   in press et al. 
      real(8), parameter :: eps=3.0d-09  ! was eps=3.0d-07 in press et al.
      real(8) :: a,gamser,gln,x
      integer :: n
      real(8) :: ap,del,sum
      gln=gammln(a)
      if(x.le.0.d+00)then
        if(x.lt.0.)stop 'aero_actv: subroutine gser: x < 0 in gser'
        gamser=0.d+00
        return
      endif
      ap=a
      sum=1.d+00/a
      del=sum
      do n=1,itmax
        ap=ap+1.d+00
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps)goto 1
      enddo
      write(*,*)'aero_actv: subroutine gser: a too large, itmax too small'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      end subroutine gser


      double precision function gammln(xx)

      implicit none
!-----------------------------------------------------------------------------------------------------------------------
!     see numerical recipes, w. press et al., 2nd edition.
!-----------------------------------------------------------------------------------------------------------------------
      real(8) :: xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
      24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
      -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      return
      end function gammln 


     subroutine erff(x,tmp)
      implicit none
!-----------------------------------------------------------------------------------------------------------------------
!     see numerical recipes, w. press et al., 2nd edition.
!-----------------------------------------------------------------------------------------------------------------------
      real(8), intent(in) :: x
      real, intent(out) :: tmp
!u    uses gammp
      tmp = 0.d0
      if(x.lt.0.0d+00)then
        tmp=-gammp(0.5d0,x**2)
      else
        tmp= gammp(0.5d0,x**2)
      endif
      
      end subroutine


      double precision function gammp(a,x)
      implicit none
!-----------------------------------------------------------------------------------------------------------------------
!     see numerical recipes, w. press et al., 2nd edition.
!-----------------------------------------------------------------------------------------------------------------------
      real(8) :: a,x
      real(8) :: gammcf,gamser,gln
      if(x.lt.0.0d+00.or.a.le.0.0d+00)then
        write(*,*)'aero_actv: function gammp: bad arguments'
      endif
      if(x.lt.a+1.0d+00)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.0d+00-gammcf
      endif
      return
      end function gammp


end module
      
