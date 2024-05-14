MODULE AERO_CONDENS
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
    !-------------------------------------------------------------------------------------------------------------------
    !
    !@sum     This module contains all sub-programs to calculate condensational growth
    !@auth  Susanne Bauer/Doug Wright, modified by x5l (Xiaohan.Li@noaa.gov)
    !3 diameters: fac = exp(b*log(sigma)**2) (Notes by x5l)
    !           (1) geometric diameter: DG
    !           (2) number mean diameter: DG*fac, b=0.5
    !           (3) surface mean diamter: DG*fac, b=1
    !           (2) volume mean diamter: DG*fac, b=1.5 
    !-------------------------------------------------------------------------------------------------------------------
    use  constants_mod, only     : PI, AVOGNO
    
contains   


!        setup_kci(pfull(i,j,k), t(i,j,k), matrix_all_pop(n)%Dg_dry(i,j,k)*fac, matrix_all_pop(n)%sigma, &
!                        kci_coef_pop(n,i,j,k), kci_aeq1_pop(n,i,j,k))

    SUBROUTINE SETUP_KCI(P,T,DP0,SIG0, KCI_COEF_DP, KCI_COEF_DP_AEQ1) !XL, volume mean diamter and sigma
    !note: figure out (1) DP0, DP -> dry or wet; (2) correction factor: thetai, what is for? 
    !SUBROUTINE SETUP_KCI
!-----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the coefficients that multiply the number
!     concentrations, or the number concentrations times the particle diameters,
!     to obtain the condensational sink for each mode or quadrature point.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I, L     ! indices
!      REAL    :: SIGMA    ! See subr. ATMOSPHERE below.
!      REAL    :: DELTA    ! See subr. ATMOSPHERE below.
!      REAL    :: THETA    ! See subr. ATMOSPHERE below.
      REAL, intent(in) :: P        ! ambient pressure [Pa]
      REAL, intent(in) :: T        ! ambient temperature [K]
      REAL, intent(in) :: DP0,SIG0    ! volume mean diameter [m] and sigma
      real, intent(out) :: KCI_COEF_DP,  KCI_COEF_DP_AEQ1
      REAL(8) :: D        ! molecular diffusivity of H2SO4 in air [m^2/s]
      REAL(8) :: C        ! mean molecular speed of H2SO4 [m/s]
      REAL(8) :: LA       ! mean free path in air [m]
                          ! 6.6328D-08 is the sea level value given in Table I.2.8
                          ! on p.10 of U.S. Standard Atmosphere 1962
      REAL(8) :: LH       ! mean free path of H2SO4 in air [m]
      REAL(8) :: BETA     ! transition regime correction to the condensational flux [1]
      REAL(8) :: KN       ! Knudsen number for H2SO4 in air [1]
      REAL(8) :: THETAI   ! monodispersity correction factor [1] (Okuyama et al. 1988)
!      REAL(8) :: SCALE_DP ! scale factor for defined table of ambient particle diameters [1]

      REAL(8), PARAMETER :: ALPHA    = 0.86D+00      ! mass accommodation coefficient [1] (Hansen, 2005)
      REAL(8), PARAMETER :: D_STDATM = 1.2250D+00    ! sea-level std. density [kg/m^3]
      REAL(8), PARAMETER :: P_STDATM = 101325.0D+00  ! sea-level std. pressure [Pa]
      REAL(8), PARAMETER :: T_STDATM = 288.15D+00    ! sea-level std. temperature [K]
      REAL(8), PARAMETER :: P0       = 101325.0D+00  ! reference pressure    [Pa] for D
      REAL(8), PARAMETER :: T0       = 273.16D+00    ! reference temperature [K]  for D
      REAL(8), PARAMETER :: D0       = 9.36D-06      ! diffusivity of H2SO4 in air [m^2/s]
                                                     ! calculated using Eqn 11-4.4 of Reid
                                                     ! et al.(1987) at 273.16 K and 101325 Pa
      real, parameter :: RGAS_SI = 8.314 !gas constant, J/mol/K
      real, parameter :: MW_H2SO4 = 98 !g/mol
!      IF( WRITE_LOG ) WRITE(AUNIT1,90002)'I','L','ZHEIGHT','P','T','D','C','LA','LH',
!     &                                   'DP0','KN','THETAI','BETA'
      KCI_COEF_DP           = 0.0D+00   ! mass accommodation coefficient is arbitrary
      KCI_COEF_DP_AEQ1      = 0.0D+00   ! mass accommodation coefficient is unity
      if (DP0 .LE. 0) then
              return
      else

      !----------------------------------------------------------------------------------------------------------------
      ! Setup table of ambient particle diameters.
      !----------------------------------------------------------------------------------------------------------------
!      SCALE_DP = ( DP_CONDTABLE_MAX / DP_CONDTABLE_MIN )**(1.0D+00/REAL(N_DP_CONDTABLE-1))
!      XLN_SCALE_DP = LOG( SCALE_DP )
!      DO I=1, N_DP_CONDTABLE
!        DP_CONDTABLE(I) = DP_CONDTABLE_MIN * SCALE_DP**(I-1)              ! [m]
!      ENDDO
!
!      DO L=1, NLAYS
!        CALL ATMOSPHERE( REAL( ZHEIGHT(L) ), SIGMA, DELTA, THETA )
!        P = DELTA * P_STDATM                                           ! [Pa]
!        T = THETA * T_STDATM                                           ! [K]
        D = D0 * ( P0 / P ) * ( T / T0 )**1.75                         ! [m^2/s]
!        DIFFCOEF_M2S(L) = D                                            ! [m^2/s]
        C = SQRT( 8.0D+00 * RGAS_SI * T / ( PI * MW_H2SO4*1.0D-03 ) )  ! [m/s]
        LA = 6.6332D-08 * ( P_STDATM / P ) * ( T / T_STDATM )          ! [m]
        LH = 3.0D+00 * D / C                                           ! [m]
        THETAI = EXP( - ( LOG(SIG0) )**2 )
!        DO I=1, NWEIGHTS
!          THETAI = EXP( - ( LOG(SIG0(I)) )**2 )                        ! [1]
!          THETA_POLY(I) = THETAI                                       ! [1] polydispersity adjustment factor
!          !------------------------------------------------------------------------------------------------------------
!          ! For the condensation sink for general use, the mean free path is
!          ! that for the condensing vapor (H2SO4), and the mass accommodation coefficient is adjustable.
!          !------------------------------------------------------------------------------------------------------------
          KN = 2.0D+00 * LH / DP0                                ! LH and DP0 in [m]
          BETA = ( 1.0D+00 + KN ) &                                     ! [1]
     &         / ( 1.0D+00 + 0.377D+00*KN + 1.33D+00*KN*(1.0D+00 + KN)/ALPHA )
          KCI_COEF_DP     = 2.0D+00 * PI * THETAI * D * BETA * DP0 ! [m^3/s] may be updated in subr. matrix
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for use in Kerminen and Kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity.
          !------------------------------------------------------------------------------------------------------------
          KN = 2.0D+00 * LA / DP0                                   ! LA and DP0 in [m]
          BETA = ( 1.0D+00 + KN ) &                                     ! [1]
     &         / ( 1.0D+00 + 0.377D+00*KN + 1.33D+00*KN*(1.0D+00 + KN)/1.0D+00 )
          KCI_COEF_DP_AEQ1 = 2.0D+00 * PI * THETAI * D * BETA * DP0 ! [m^3/s] may be updated in subr. matrix
  endif
      END SUBROUTINE 


END MODULE
