module aero_condens
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
    !-------------------------------------------------------------------------------------------------------------------
    !
    !@sum     this module contains all sub-programs to calculate condensational growth
    !@auth  susanne bauer/doug wright, modified by x5l (xiaohan.li@noaa.gov)
    !3 diameters: fac = exp(b*log(sigma)**2) (notes by x5l)
    !           (1) geometric diameter: dg
    !           (2) number mean diameter: dg*fac, b=0.5
    !           (3) surface mean diamter: dg*fac, b=1
    !           (2) volume mean diamter: dg*fac, b=1.5 
    !-------------------------------------------------------------------------------------------------------------------
    use  constants_mod, only     : pi, avogno
    
contains   


!        setup_kci(pfull(i,j,k), t(i,j,k), matrix_all_pop(n)%dg_dry(i,j,k)*fac, matrix_all_pop(n)%sigma, &
!                        kci_coef_pop(n,i,j,k), kci_aeq1_pop(n,i,j,k))

    subroutine setup_kci(p,t,dp0,sig0, kci_coef_dp, kci_coef_dp_aeq1) !xl, volume mean diamter and sigma
    !note: figure out (1) dp0, dp -> dry or wet; (2) correction factor: thetai, what is for? 
    !subroutine setup_kci
!-----------------------------------------------------------------------------------------------------------------------
!     routine to calculate the coefficients that multiply the number
!     concentrations, or the number concentrations times the particle diameters,
!     to obtain the condensational sink for each mode or quadrature point.
!-----------------------------------------------------------------------------------------------------------------------
      implicit none
      integer :: i, l     ! indices
!      real    :: sigma    ! see subr. atmosphere below.
!      real    :: delta    ! see subr. atmosphere below.
!      real    :: theta    ! see subr. atmosphere below.
      real, intent(in) :: p        ! ambient pressure [pa]
      real, intent(in) :: t        ! ambient temperature [k]
      real, intent(in) :: dp0,sig0    ! volume mean diameter [m] and sigma
      real, intent(out) :: kci_coef_dp,  kci_coef_dp_aeq1
      real(8) :: d        ! molecular diffusivity of h2so4 in air [m^2/s]
      real(8) :: c        ! mean molecular speed of h2so4 [m/s]
      real(8) :: la       ! mean free path in air [m]
                          ! 6.6328d-08 is the sea level value given in table i.2.8
                          ! on p.10 of u.s. standard atmosphere 1962
      real(8) :: lh       ! mean free path of h2so4 in air [m]
      real(8) :: beta     ! transition regime correction to the condensational flux [1]
      real(8) :: kn       ! knudsen number for h2so4 in air [1]
      real(8) :: thetai   ! monodispersity correction factor [1] (okuyama et al. 1988)
!      real(8) :: scale_dp ! scale factor for defined table of ambient particle diameters [1]

      real(8), parameter :: alpha    = 0.86d+00      ! mass accommodation coefficient [1] (hansen, 2005)
      real(8), parameter :: d_stdatm = 1.2250d+00    ! sea-level std. density [kg/m^3]
      real(8), parameter :: p_stdatm = 101325.0d+00  ! sea-level std. pressure [pa]
      real(8), parameter :: t_stdatm = 288.15d+00    ! sea-level std. temperature [k]
      real(8), parameter :: p0       = 101325.0d+00  ! reference pressure    [pa] for d
      real(8), parameter :: t0       = 273.16d+00    ! reference temperature [k]  for d
      real(8), parameter :: d0       = 9.36d-06      ! diffusivity of h2so4 in air [m^2/s]
                                                     ! calculated using eqn 11-4.4 of reid
                                                     ! et al.(1987) at 273.16 k and 101325 pa
      real, parameter :: rgas_si = 8.314 !gas constant, j/mol/k
      real, parameter :: mw_h2so4 = 98 !g/mol
!      if( write_log ) write(aunit1,90002)'i','l','zheight','p','t','d','c','la','lh',
!     &                                   'dp0','kn','thetai','beta'
      kci_coef_dp           = 0.0d+00   ! mass accommodation coefficient is arbitrary
      kci_coef_dp_aeq1      = 0.0d+00   ! mass accommodation coefficient is unity
      if (dp0 .le. 0) then
              return
      else

      !----------------------------------------------------------------------------------------------------------------
      ! setup table of ambient particle diameters.
      !----------------------------------------------------------------------------------------------------------------
!      scale_dp = ( dp_condtable_max / dp_condtable_min )**(1.0d+00/real(n_dp_condtable-1))
!      xln_scale_dp = log( scale_dp )
!      do i=1, n_dp_condtable
!        dp_condtable(i) = dp_condtable_min * scale_dp**(i-1)              ! [m]
!      enddo
!
!      do l=1, nlays
!        call atmosphere( real( zheight(l) ), sigma, delta, theta )
!        p = delta * p_stdatm                                           ! [pa]
!        t = theta * t_stdatm                                           ! [k]
        d = d0 * ( p0 / p ) * ( t / t0 )**1.75                         ! [m^2/s]
!        diffcoef_m2s(l) = d                                            ! [m^2/s]
        c = sqrt( 8.0d+00 * rgas_si * t / ( pi * mw_h2so4*1.0d-03 ) )  ! [m/s]
        la = 6.6332d-08 * ( p_stdatm / p ) * ( t / t_stdatm )          ! [m]
        lh = 3.0d+00 * d / c                                           ! [m]
        thetai = exp( - ( log(sig0) )**2 )
!        do i=1, nweights
!          thetai = exp( - ( log(sig0(i)) )**2 )                        ! [1]
!          theta_poly(i) = thetai                                       ! [1] polydispersity adjustment factor
!          !------------------------------------------------------------------------------------------------------------
!          ! for the condensation sink for general use, the mean free path is
!          ! that for the condensing vapor (h2so4), and the mass accommodation coefficient is adjustable.
!          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * lh / dp0                                ! lh and dp0 in [m]
          beta = ( 1.0d+00 + kn ) &                                     ! [1]
     &         / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/alpha )
          kci_coef_dp     = 2.0d+00 * pi * thetai * d * beta * dp0 ! [m^3/s] may be updated in subr. matrix
          !------------------------------------------------------------------------------------------------------------
          ! for the condensation sink for use in kerminen and kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity.
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * la / dp0                                   ! la and dp0 in [m]
          beta = ( 1.0d+00 + kn ) &                                     ! [1]
     &         / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/1.0d+00 )
          kci_coef_dp_aeq1 = 2.0d+00 * pi * thetai * d * beta * dp0 ! [m^3/s] may be updated in subr. matrix
  endif
      end subroutine 


end module
