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
!    subroutine condens_calc(XH2SO4_INIT)
!      !-------------------------------------------------------------------------------------------------------------------
!      ! For the calculation of DNDT and DMDT_SO4, it is assumed that all of the
!      !   current H2SO4 concentration was produced by gas-phase oxidation of
!      !   SO2 during the current time step. SO4RATE is the average production rate
!      !   of H2SO4 (represented as SO4, MW=96g/mol) over the time step based on this
!      !   assumption. DMDT_SO4 and DNDT are limited by the available H2SO4.
!      !
!      ! When appropriate, the nucleation rate is calculated using a steady-state
!      !   H2SO4 concentration. This is done to avoid the spuriously high nucleation
!      !   rates that would result if the current H2SO4 concentration, which has
!      !   accumulated without loss over the time step, were used. When the
!      !   condensational sink (KC) is small enough such that steady-state will not
!      !   be reached during the time step, an estimate of the H2SO4 concentration
!      !   at the mid-point of the time step is used in nucleation and condensation
!      !   calculations.
!      !
!      ! The condensational sink KC is the first-order rate constant for the
!      !   loss of H2SO4 due to condensation. This is summed over all modes.
!      !   For the Kerminen and Kulmala (2002) parameterization for the conversion
!      !   of the nucleation rate to the new particle formation rate, the
!      !   condensational sink KC_AEQ1 should obtained with the mass accommodation
!      !   coefficient set of unity.
!      !----------------------------------------------------------------------------------------------------------------
!      real, intent(in) :: XH2SO4_INIT,TSTEP                        ! input H2SO4 concentration [ugSO4/m^3]
!      real, intent(in) :: npop
!      real, intent(in) :: KCI_COEF_DP(npop-1), NI(npop-1)              ! condensational coefficient and number for each population at a grid
!      real, intent(in) :: KCI_COEF_DP_AEQ1(npop-1)                                                             ! KCI_COEF_DP(npop), NI(npop) 
!      real:: XH2SO4_NUCL                                     ! H2SO4 (as SO4) conc. used in nucleation and GR calculation [ugSO4/m^3]
!      REAL(8), PARAMETER :: XNTAU = 2.0D+00                 ! number of time constants in the current time step
!      real, PARAMETER :: MW_SO4 = 96
!      REAL, PARAMETER :: KCMIN = 1.0D-08 ! [1/s] minimum condensational sink - see notes of 10-18-06
!      REAL, PARAMETER :: XH2SO4_NUCL_MIN_NCM3 = 1.00D+03 ! min. [H2SO4] to enter nucleation calculations [#/cm^3]
!      REAL, PARAMETER :: XH2SO4_NUCL_MIN = XH2SO4_NUCL_MIN_NCM3 * MW_SO4 * 1.0D+12 / AVOGNO  ! convert to [ugSO4/m^3] 
!      XH2SO4_NUCL =  XH2SO4_NUCL_MIN !TINYNUMER ! XH2SO4_NUCL_MIN              ! for the case  XH2SO4_INIT .LT. XH2SO4_NUCL_MIN
!      KC = SUM( KCI_COEF_DP(:)*NI(:) )                  ! total condensational sink [1/s] for any value of the
!      KC = MAX( KC, KCMIN )                                  !   mass accommodation coefficient [1/s]
!      IF( (I_AKK > 0) .AND. (XH2SO4_INIT > XH2SO4_NUCL_MIN) ) THEN
!        KC_AEQ1 = SUM( KCI_COEF_DP_AEQ1(:)*NI(:) )      ! total condensational sink for the
!        KC_AEQ1 = MAX( KC_AEQ1, KCMIN )                      !   mass accommodation coefficient set to unity  [1/s]
!        XNH3 = CONVNH3 * GAS( GAS_NH3 ) * TK / PRES          ! NH3 concentration; from [ug/m^3] to [ppmV]
!        SO4RATE = XH2SO4_INIT / TSTEP                        ! average H2SO4 production rate [ugSO4/m^3/s]
!        IF(KC*TSTEP .GE. XNTAU ) THEN                        ! invoke steady-state assumption
!          IH2SO4_PATH = 1
!          XH2SO4_SS = MIN( SO4RATE/KC, XH2SO4_INIT )         ! steady-state H2SO4 [ugSO4/m^3]
!          CALL STEADY_STATE_H2SO4(PRES,RH,TK,FLAND,XH2SO4_SS,SO4RATE,XNH3,KC,TSTEP,XH2SO4_SS_WNPF)
!          XH2SO4_NUCL = XH2SO4_SS_WNPF                       ! [H2SO4] for nucl., GR, and cond. calculation [ugSO4/m^3]
!        ELSE
!          IH2SO4_PATH = 2
!          XH2SO4_NUCL = SO4RATE / ( (2.0D+00/TSTEP) + KC )   ! use [H2SO4] at mid-time step [ugSO4/m^3]
!        ENDIF
!        CALL NPFRATE(PRES,RH,TK,FLAND,XH2SO4_NUCL,SO4RATE,XNH3,KC_AEQ1,DNDT,DMDT_SO4,0)
!        ! WRITE(34,90010) TK,RH,UGM3_NCM3*XH2SO4_INIT,UGM3_NCM3*XH2SO4_NUCL,KC,XNH3,1.0D-06*DNDT,IH2SO4_PATH
!!        IF( WRITE_LOG ) THEN
!!          WRITE(AUNIT1,'(/A)')'NEW PARTICLE FORMATION'
!!          WRITE(AUNIT1,'(/A)')'PRES,RH,TK,XH2SO4_INIT(ug/m3),XH2SO4_NUCL(#/cm3),SO4RATE,XNH3,KC,DNDT,DMDT_SO4'
!!          WRITE(AUNIT1,90003)  PRES,RH,TK,XH2SO4_INIT,       XH2SO4_NUCL*UGM3_NCM3,SO4RATE,XNH3,KC,DNDT,DMDT_SO4
!!          WRITE(AUNIT1,'(/A)')'TK,RH,H2SO4(#.cm^3),H2SO4_NUCL(#/cm3),KC,NH3,DNDT(#/cm3/s),IH2SO4_PATH'
!!          WRITE(AUNIT1,90010)  TK,RH,UGM3_NCM3*XH2SO4_INIT,UGM3_NCM3*XH2SO4_NUCL,KC,XNH3,1.0D-06*DNDT,IH2SO4_PATH
!!        ENDIF
!      ELSE
!        DNDT     = 0.0D+00
!        DMDT_SO4 = 0.0D+00
!      ENDIF
!
!      !----------------------------------------------------------------------------------------------------------------
!      ! Get the Pgrowth_i,q terms due to condensation and gas-particle
!      !   mass transfer for each mode (or quadrature point) and add them to the PIQ array.
!      !
!      ! The net loss of H2SO4 due to both secondary particle formation
!      !   and condensation should not exceed the current H2SO4 concentration.
!      !   This is enforced by rescaling the two H2SO4 consumption rates in a way
!      !   that preserves the relative magnitudes of these two loss processes.
!      !   The net condensation rate PQ_GROWTH is calculated from XH2SO4_NUCL, not
!      !   the total accumulated H2SO4 in XH2SO4_INIT, for balance with the
!      !   treatment of new particle formation above.
!      !
!      ! The expression in parentheses on the rhs of PIQ is that for h_i,
!      !
!      !   h_i = KCI_COEF_DP(i,ILAY) * NI(i) / KC
!      !
!      !   the ratio of the condensational sink of mode or quadrature point I
!      !   to the total condensational sink.
!      !
!      ! At this point, XH2SO4_INIT = GAS( GAS_H2SO4 ) + TINYNUMER.
!      !
!      ! The H2SO4 concentration GAS( GAS_H2SO4 ) is updated here.
!      !----------------------------------------------------------------------------------------------------------------
!      PQ_GROWTH = XH2SO4_NUCL * ( 1.0D+00 - EXP(-KC*TSTEP) ) / TSTEP            ! [ugSO4/m^3/s]
!      TOT_H2SO4_LOSS = ( DMDT_SO4 + PQ_GROWTH ) * TSTEP                         ! [ugSO4/m^3]
!      IF ( TOT_H2SO4_LOSS .GT. XH2SO4_INIT ) THEN                               ! XH2SO4_INIT=GAS(GAS_H2SO4)+TINYNUMER
!        DMDT_SO4  = DMDT_SO4  * ( XH2SO4_INIT / ( TOT_H2SO4_LOSS + TINYDENOM ) )! [ugSO4/m^3/s]
!        DNDT      = DNDT      * ( XH2SO4_INIT / ( TOT_H2SO4_LOSS + TINYDENOM ) )! [  #  /m^3/s]
!        PQ_GROWTH = PQ_GROWTH * ( XH2SO4_INIT / ( TOT_H2SO4_LOSS + TINYDENOM ) )! [ugSO4/m^3/s]
!        GAS( GAS_H2SO4 ) = TINYNUMER
!      ELSE
!        GAS( GAS_H2SO4 ) = GAS( GAS_H2SO4 ) - TOT_H2SO4_LOSS + TINYNUMER
!      ENDIF
!      DIAGTMP1(2,NUMB_MAP(1))               = DNDT
!      DIAGTMP1(9,SULF_MAP(PROD_INDEX_SULF)) = DMDT_SO4
!
!      ! WRITE(35,'(8D13.5)')PIQ(:,PROD_INDEX_SULF), KCI_COEF_DP(:,ILAY),NI(:),KC
!
!      PIQTMP(:,PROD_INDEX_SULF) = ( KCI_COEF_DP(:,ILAY)*NI(:)/KC ) * PQ_GROWTH
!      PIQ   (:,PROD_INDEX_SULF) = PIQ(:,PROD_INDEX_SULF) + PIQTMP(:,PROD_INDEX_SULF)
!
!      DIAGTMP1(11,MASS_MAP(:,PROD_INDEX_SULF)) = PIQTMP(:,PROD_INDEX_SULF)
!      !-----------------------------------------------------------------------------------------------------------------
!      ! Add the secondary particle formation term DNDT [#/m^3/s] for the number concentration term.
!      ! Add the secondary particle formation term DMDT_SO4 [ug/m^3/s], calculated above,
!      !   to the total mass production rate array PIQ for the AKK mode.
!      ! If the AKK mode is absent, the secondary particle formation terms still go into mode 1.
!      !-----------------------------------------------------------------------------------------------------------------
!      CI(1) = CI(1) + DNDT                                        ! add secondary particle formation number term
!      PIQ(1,PROD_INDEX_SULF) = PIQ(1,PROD_INDEX_SULF) + DMDT_SO4  ! add secondary particle formation mass   term
!
!      IF( WRITE_LOG ) THEN
!        WRITE(AUNIT1,'(/A,5X,3D15.8)')'XH2SO4_INIT, XH2SO4_NUCL, PQ_GROWTH = ', XH2SO4_INIT, XH2SO4_NUCL, PQ_GROWTH
!        WRITE(AUNIT1,*)'PIQ(1,PROD_INDEX_SULF) = ', PIQ(1,PROD_INDEX_SULF)
!      ENDIF
!
!    end subroutine


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
      END SUBROUTINE 


END MODULE
