      MODULE AERO_COAG
!-------------------------------------------------------------------------------------------------------------------------------------
!
!@sum     This module contains all routines for coagulation coefficients.
!@auth    Susanne Bauer/Doug Wright
!
!     The routines are included in this file in the following order.
!     Only routine KBARNIJ is called after initialization is complete.
!      
!           SETUP_KIJ_DIAMETERS
!           SETUP_KIJ
!           SETUP_KIJ_TABLES
!           BUILD_KIJ_TABLES
!           GET_KNIJ
!           GET_KBARNIJ
!           BROWNIAN_COAG_COEF
!           CBDE_COAG_COEF
!           GRAVCOLL_COAG_COEF
!           TURB_COAG_COEF
!           TOTAL_COAG_COEF
!           TEST_COAG_COEF
!
!-------------------------------------------------------------------------------------------------------------------------------------
!      USE AERO_PARAM,  ONLY: AUNIT1, AUNIT2, WRITE_LOG, KIJ_NDGS_SET, DPMIN_GLOBAL
!      USE AERO_CONFIG, ONLY: NWEIGHTS
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
        USE AERO_COAG_CONFIG,  ONLY: CITABLE ! NMODES
        IMPLICIT NONE

!-------------------------------------------------------------------------------------------------------------------------------
! Constant coagulation coefficients for each mode-mode (I-J) interaction,
! based upon characteristic sizes for each mode.
! These are not mode-averaged and are no longer in use.
!-------------------------------------------------------------------------------------------------------------------------------
!      REAL, SAVE :: KIJ(NWEIGHTS,NWEIGHTS)  ! [m^3/s]
!XLLLLLLLLLLLLLLL
INTEGER, PARAMETER :: NMODES = 12
INTEGER, PARAMETER :: NWEIGHTS = NMODES
!-------------------------------------------------------------------------------------------------------------------------------
! Mode-average coagulation coefficients for mode-mode (I-J) interactions.
!-------------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: KIJ_NDGS_SET = 81
      INTEGER, PARAMETER :: KIJ_NDGS  = KIJ_NDGS_SET ! number of geo. mean diameters
      INTEGER, PARAMETER :: KIJ_NSGS  =      3       ! number of geo. std. deviations
      REAL, PARAMETER :: KIJ_TEMP1 =    325.0     ! [K]  \.
      REAL, PARAMETER :: KIJ_TEMP2 =    260.0     ! [K]  - must have T1 > T2 > T3
      REAL, PARAMETER :: KIJ_TEMP3 =    200.0     ! [K]  / 
      REAL, PARAMETER :: KIJ_PRES1 = 101325.0     ! [Pa] \.
      REAL, PARAMETER :: KIJ_PRES2 =  10132.50    ! [Pa] - must have p1 > p2 > p3
      REAL, PARAMETER :: KIJ_PRES3 =   1013.250   ! [Pa] /
      REAL, PARAMETER :: KIJ_SIGM1 =      1.6     ! [1]
      REAL, PARAMETER :: KIJ_SIGM2 =      1.8     ! [1]
      REAL, PARAMETER :: KIJ_SIGM3 =      2.0     ! [1]
      !-------------------------------------------------------------------------------------------------------------------------------
      ! KIJ_DGMIN must be smaller than DPMIN_GLOBAL / Smax, where Smax = exp[1.5(ln Sigma_max)^2] where Sigma_max is the largest
      ! lognormal geometric standard deviation occurring for any mode. Currently set for sigma_max = 2.0. 
      !-------------------------------------------------------------------------------------------------------------------------------
      REAL, PARAMETER :: DPMIN_GLOBAL =  0.001D-06   ! [m] -  1 nm
      REAL, PARAMETER :: KIJ_DGMIN = 1.0D+06 * DPMIN_GLOBAL / 2.1D+00           ! [um] If any mode has Sigma>2.0, must modify.
      REAL, PARAMETER :: KIJ_DGMAX =    100.0000                                ! [um]
      REAL      :: K0IJ_TEMP1PRES1(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP1PRES2(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP1PRES3(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP2PRES1(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP2PRES2(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP2PRES3(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP3PRES1(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP3PRES2(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K0IJ_TEMP3PRES3(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP1PRES1(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP1PRES2(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP1PRES3(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP2PRES1(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP2PRES2(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP2PRES3(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP3PRES1(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP3PRES2(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      REAL      :: K3IJ_TEMP3PRES3(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)   ! [m^3/s]
      !-----------------------------------------------------------------------------------------------------------------
      ! KIJ_DIAMETERS contains the values of Dg used in building the lookup tables.
      ! KIJ_SIGMAS    contains the values of Sigmag used in building the lookup tables.
      ! INDEX_SIGG(I) is the index of KIJ_SIGMAS to obtain the Sigmag value for mode I.
      !-----------------------------------------------------------------------------------------------------------------
      REAL                :: KIJ_DIAMETERS(KIJ_NDGS)                             ! [um]
      REAL, DIMENSION(KIJ_NSGS) :: KIJ_SIGMAS = (/ KIJ_SIGM1, KIJ_SIGM2, KIJ_SIGM3 /)  ! [1]
      INTEGER               :: INDEX_SIGG(NWEIGHTS)                                ! [1]

     CONTAINS


!      SUBROUTINE SETUP_KIJ_DIAMETERS
!!-----------------------------------------------------------------------------------------------------------------------
!!     Routine to define the geometric mean diameters Dg of the lookup tables.
!!-----------------------------------------------------------------------------------------------------------------------
!      IMPLICIT NONE
!      INTEGER :: I
!      REAL :: E, SCALE
!      E = 1.0 / REAL( KIJ_NDGS - 1 )
!      SCALE = ( KIJ_DGMAX / KIJ_DGMIN )**E
!      DO I=1, KIJ_NDGS
!        KIJ_DIAMETERS(I) = KIJ_DGMIN * SCALE**(I-1)             ! [um]
!        ! WRITE(*,'(I6,F16.6)') I, KIJ_DIAMETERS(I)
!      ENDDO
!      RETURN
!      END SUBROUTINE SETUP_KIJ_DIAMETERS
!
!
!
!      SUBROUTINE SETUP_KIJ_TABLES
!!-----------------------------------------------------------------------------------------------------------------------
!!     Routine to setup tables of mode-average coagulation coefficients [m^3/s].
!!     Several temperatures and pressures are currently used. 
!!-----------------------------------------------------------------------------------------------------------------------
!      IMPLICIT NONE
!      INTEGER :: I
!      REAL :: K0IJ_TABLE(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)        ! [m^3/s]
!      REAL :: K3IJ_TABLE(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)        ! [m^3/s]
!
!      !-------------------------------------------------------------------------
!      ! Build the table for each choice of temperature and pressure. 
!      !-------------------------------------------------------------------------
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP1,KIJ_PRES1,K0IJ_TABLE,K3IJ_TABLE)  ! T1, p1
!      K0IJ_TEMP1PRES1(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP1PRES1(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP1,KIJ_PRES2,K0IJ_TABLE,K3IJ_TABLE)  ! T1, p2
!      K0IJ_TEMP1PRES2(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP1PRES2(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP1,KIJ_PRES3,K0IJ_TABLE,K3IJ_TABLE)  ! T1, p3
!      K0IJ_TEMP1PRES3(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP1PRES3(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP2,KIJ_PRES1,K0IJ_TABLE,K3IJ_TABLE)  ! T2, p1
!      K0IJ_TEMP2PRES1(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP2PRES1(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP2,KIJ_PRES2,K0IJ_TABLE,K3IJ_TABLE)  ! T2, p2
!      K0IJ_TEMP2PRES2(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP2PRES2(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP2,KIJ_PRES3,K0IJ_TABLE,K3IJ_TABLE)  ! T2, p3
!      K0IJ_TEMP2PRES3(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP2PRES3(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP3,KIJ_PRES1,K0IJ_TABLE,K3IJ_TABLE)  ! T3, p1
!      K0IJ_TEMP3PRES1(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP3PRES1(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP3,KIJ_PRES2,K0IJ_TABLE,K3IJ_TABLE)  ! T3, p2
!      K0IJ_TEMP3PRES2(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP3PRES2(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      CALL BUILD_KIJ_TABLES(KIJ_TEMP3,KIJ_PRES3,K0IJ_TABLE,K3IJ_TABLE)  ! T3, p3
!      K0IJ_TEMP3PRES3(:,:,:,:) = K0IJ_TABLE(:,:,:,:)      ! [m^3/s]
!      K3IJ_TEMP3PRES3(:,:,:,:) = K3IJ_TABLE(:,:,:,:)      ! [m^3/s]
!
!      RETURN
!      END SUBROUTINE SETUP_KIJ_TABLES
!
!
!      SUBROUTINE BUILD_KIJ_TABLES( TEMP, PRES, K0IJ_TABLE, K3IJ_TABLE )
!!-------------------------------------------------------------------------------------------------------------------------------------
!!     Routine to setup a table of mode-average coagulation coefficients [m^3/s]
!!       for a given temperature and pressure.
!!-------------------------------------------------------------------------------------------------------------------------------------
!      IMPLICIT NONE
!      INTEGER :: I, J, K, L
!      REAL :: TEMP                                             ! [K]
!      REAL :: PRES                                             ! [Pa]
!      REAL :: K0IJ_TABLE(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)  ! [m^3/s]
!      REAL :: K3IJ_TABLE(KIJ_NDGS,KIJ_NSGS,KIJ_NDGS,KIJ_NSGS)  ! [m^3/s]
!      REAL :: K0IJ, K3IJ                                       ! [m^3/s]
!
!      DO I=1, KIJ_NDGS
!      DO J=1, KIJ_NDGS
!      DO K=1, KIJ_NSGS
!      DO L=1, KIJ_NSGS
!        CALL GET_KNIJ(TEMP,PRES,KIJ_DIAMETERS(I),KIJ_SIGMAS(K),KIJ_DIAMETERS(J),KIJ_SIGMAS(L),K0IJ,K3IJ)
!        K0IJ_TABLE(I,K,J,L) = K0IJ    ! [m^3/s]
!        K3IJ_TABLE(I,K,J,L) = K3IJ    ! [m^3/s]
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!
!      RETURN
!      END SUBROUTINE BUILD_KIJ_TABLES
!
!       !call GET_KNIJ(t(i,j,k), pfull(i,j,k), dg_ip(i,j,k), sig_ip, dg_jp(i,j,k), sig_jp, KBAR0_IJ(ip, jp, i,j,k), &
!       !                                       KBAR3_IJ(ip, jp, i,j,k))
      SUBROUTINE GET_KNIJ( TEMP, PRES, DGI, SIGGI, DGJ, SIGGJ, K0IJ, K3IJ )
!-------------------------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the mode-average coagulation coefficients 
!     K0IJ and K3IJ [m^3/s] for the coagulation of mode I, having lognormal
!     parameters DGI and SIGGI, with mode J, having lognormal parameters
!     DGJ and SIGGJ, at a given temperature and pressure.
!
!     The integrals necessary to evaluate the mode-average coagulation
!     coefficients are evaluated using n-point quadrature from 2n moments
!     calculated from the lognormal parameters for each mode.
!-------------------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Arguments.

      REAL, intent(in) :: TEMP                           ! [K]  ambient temperature
      REAL, intent(in) :: PRES                           ! [Pa] ambient pressure
      REAL, intent(in) :: DGI, DGJ                       ! [um] geometric mean diameters
      REAL, intent(in) :: SIGGI, SIGGJ                   ! [1]  geometric standard deviations
      REAL, intent(out) :: K0IJ, K3IJ                     ! [m^3/s] mode-average coagulation coefficients

      ! Local variables.

      INTEGER, PARAMETER :: NQUADR = 4         ! number of quadrature points for each mode

      INTEGER :: I, J, K, L                     ! loop indices
      REAL(8) :: SGI, SGJ                       ! [1] function of geometric standard deviation
      REAL(8) :: UKI(2*NQUADR)                 ! [um^k] normalized moments for mode I
      REAL(8) :: UKJ(2*NQUADR)                 ! [um^k] normalized moments for mode J
      REAL(8) :: XI(NQUADR), XJ(NQUADR)       ! [um] quadrature abscissas for modes I and J
      REAL(8) :: WI(NQUADR), WJ(NQUADR)       ! [1] normalized quadrature weights for modes I and J
      REAL(8) :: ZFI, ZFJ                       ! [1] flags for failed quadrature inversion
      REAL(8) :: K0IJ_TMP, K3IJ_TMP             ! [m^3/s] double precision accumulators for K0IJ and K3IJ
      REAL :: DI, DJ                         ! [um] single precision particle diameters for modes I and J
      REAL :: BETAIJ                         ! [m^3/s] single precision coagulation coefficient

      SGI = EXP( 0.5D+00 * ( LOG( DBLE(SIGGI) )**2 ) )
      SGJ = EXP( 0.5D+00 * ( LOG( DBLE(SIGGJ) )**2 ) )
      !-------------------------------------------------------------------------
      ! WRITE(*,*)         'DGI,SIGGI,SGI,DGJ,SIGGJ,SGJ'
      ! WRITE(*,'(6F13.8)') DGI,SIGGI,SGI,DGJ,SIGGJ,SGJ 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Compute the first 2*NQUADR diameter moments for each lognormal mode.
      !-------------------------------------------------------------------------
      DO L=1, 2*NQUADR
        K = L-1
        UKI(L) = DBLE(DGI)**K * SGI**(K*K)
        UKJ(L) = DBLE(DGJ)**K * SGJ**(K*K)
        !-----------------------------------------------------------------------
        ! WRITE(*,'(I6,2D15.5)') K, UKI(L), UKJ(L)
        !-----------------------------------------------------------------------
      ENDDO

      !-------------------------------------------------------------------------
      ! Get the quadrature abscissas and weights for modes I and J.
      !-------------------------------------------------------------------------
      CALL GAUSS(NQUADR,UKI,XI,WI,ZFI)
      IF( ZFI .GT. 1.0D-15 ) THEN
        WRITE(*,*)'Failed quadrature for mode I in subr. GET_KNIJ'
        WRITE(*,*)'ABSCISSAS = ', XI(:)
        WRITE(*,*)'WEIGHTS   = ', WI(:)
        STOP
      ENDIF
      CALL GAUSS(NQUADR,UKJ,XJ,WJ,ZFJ)
      IF( ZFJ .GT. 1.0D-15 ) THEN
        WRITE(*,*)'Failed quadrature for mode J in subr. GET_KNIJ'
        WRITE(*,*)'ABSCISSAS = ', XJ(:)
        WRITE(*,*)'WEIGHTS   = ', WJ(:)
        STOP
      ENDIF

      !-------------------------------------------------------------------------
      ! Write the abscissas and weights for modes I and J.
      !-------------------------------------------------------------------------
      ! DO L=1, NQUADR
      !   WRITE(*,'(I6,4D15.5)') L, XI(L), WI(L), XJ(L), WJ(L)
      ! ENDDO
      !-------------------------------------------------------------------------

      K0IJ_TMP = 0.0D+00
      K3IJ_TMP = 0.0D+00

      DO I=1, NQUADR
      DO J=1, NQUADR
        DI = REAL(XI(I))                                           ! convert to single precision
        DJ = REAL(XJ(J))                                           ! convert to single precision
        CALL TOTAL_COAG_COEF   ( DI, DJ, TEMP, PRES, BETAIJ )      ! all variables single precision
        K0IJ_TMP = K0IJ_TMP + DBLE(BETAIJ)*WI(I)*WJ(J)             ! all factors are REAL(8)
        K3IJ_TMP = K3IJ_TMP + DBLE(BETAIJ)*WI(I)*WJ(J)*XI(I)**3    ! all factors are REAL(8)
      ENDDO
      ENDDO

      K0IJ = REAL(K0IJ_TMP)                                        ! divide by U0i*U0j=   1.0*1.0
      K3IJ = REAL(K3IJ_TMP/UKI(4))                                 ! divide by U3i*U0j=UKI(4)*1.0

      !XL debug
      if ((K0IJ > 1) .or. (K3IJ > 1)) then
                write(*,*) 'PRES, TEMP, DGI, SIGGI, DGJ, SIGGJ, K0IJ, K3IJ', PRES, TEMP, DGI, SIGGI, DGJ, SIGGJ, K0IJ, K3IJ
      endif

!-------------------------------------------------------------------------------
!     WRITE(AUNIT2,'(A,4F9.4,2E13.5)') 'DGI, SIGGI, DGJ, SIGGJ, K0IJ, K3IJ = ',
!    &                                  DGI, SIGGI, DGJ, SIGGJ, K0IJ, K3IJ
!-------------------------------------------------------------------------------
      END SUBROUTINE GET_KNIJ

     SUBROUTINE GAUSS(N,X,R,W,ZF)
!----------------------------------------------------------------------------------------------------------------------
!@sum    Perform an n-point quadrature from 2n moments to yield  n abscissas and n weights.
!@auth   Susanne Bauer/Doug Wright
!----------------------------------------------------------------------------------------------------------------------
     IMPLICIT NONE

     ! Arguments.

     INTEGER, INTENT(IN   ) :: N               ! number of quadrature points
     REAL(8), INTENT(INOUT) :: X(2*N)          ! moments
     REAL(8), INTENT(  OUT) :: R(N)            ! abscissas
     REAL(8), INTENT(  OUT) :: W(N)            ! weights
     REAL(8), INTENT(  OUT) :: ZF              ! =0.0 successful quadrature, =1.0 failed quadrature

     ! Local variables.

     INTEGER :: IFAILTQL
     REAL(8) :: A(N),B(N),ANU(2*N),AMU0

     ZF = 0.0D+00                              ! successful quadrature
     AMU0 = X(1)                               ! normalizing moment
     ANU(:) = X(:)/AMU0                        ! normalize the moments
     CALL ORTHOG(N,ANU,A,B)
     CALL GAUCOF(N,A,B,AMU0,R,W,IFAILTQL)
     IF(     MINVAL(R(:)) .LT. 0.0D+00 ) THEN  ! failed quadrature
       ZF = 1.0D+00
     ELSEIF( MINVAL(W(:)) .LT. 0.0D+00 ) THEN  ! failed quadrature
       ZF = 1.0D+00
     ELSEIF( IFAILTQL .GT. 0 ) THEN            ! failed quadrature
       ZF = 1.0D+00
     ENDIF
     END SUBROUTINE GAUSS
      
     SUBROUTINE ORTHOG(N,ANU,A,B)

             IMPLICIT NONE
             !----------------------------------------------------------------------------------------------------------------------
             !     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
             !----------------------------------------------------------------------------------------------------------------------
             INTEGER :: L,N,K,NMAX
             PARAMETER (NMAX=30)
             REAL(8) :: A(N),ANU(2*N),B(N),SIG(2*NMAX+1,2*NMAX+1)

             DO 11 L=3,2*N
             SIG(1,L)=0.D+00
             11    CONTINUE
             DO 12 L=2,2*N+1
             SIG(2,L)=ANU(L-1)
             12    CONTINUE
             A(1)=ANU(2)/ANU(1)
             B(1)=0.D+00
             DO 14 K=3,N+1
             DO 13 L=K,2*N-K+3
             SIG(K,L)=SIG(K-1,L+1)-A(K-2)*SIG(K-1,L)-B(K-2)*SIG(K-2,L)
             IF(SIG(K,K).LE.0.D+00) SIG(K,K) = 1.D-20
             13      CONTINUE
             A(K-1)=SIG(K,K+1)/SIG(K,K)-SIG(K-1,K)/SIG(K-1,K-1)
             B(K-1)=SIG(K,K)/SIG(K-1,K-1)
             14    CONTINUE
             RETURN
     END SUBROUTINE ORTHOG

     SUBROUTINE GAUCOF(N,A,B,AMU0,X,W,IFAILTQL)

             IMPLICIT NONE
             !----------------------------------------------------------------------------------------------------------------------
             !     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
             !----------------------------------------------------------------------------------------------------------------------
             INTEGER :: I,J,N,NMAX, IFAILTQL
             PARAMETER (NMAX=30)
             REAL(8) :: A(N),B(N),W(N),X(N),Z(NMAX,NMAX),AMU0
             IFAILTQL = 0
             DO 12 I=1,N
             IF(I.NE.1)B(I)=SQRT(B(I))
             DO 11 J=1,N
             IF(I.EQ.J)THEN
                     Z(I,J)=1.D+00
             ELSE
                     Z(I,J)=0.D+00
             ENDIF
             11      CONTINUE
             12    CONTINUE
             CALL TQLI(A,B,N,NMAX,Z,IFAILTQL)
             IF(IFAILTQL.GT.0) RETURN
             !----------------------------------------------------------------------------------------------------------------------
             !     Ordering of the abscissas is usually not needed.
             !----------------------------------------------------------------------------------------------------------------------
             !     CALL EIGSRT(A,Z,N,NMAX)
             !----------------------------------------------------------------------------------------------------------------------
             DO 13 I=1,N
             X(I)=A(I)
             W(I)=AMU0*Z(1,I)**2
             !--------------------------------------------------------------------------------------------------------------
             ! AVOID ZERO WEIGHTS.
             !--------------------------------------------------------------------------------------------------------------
             ! IF(W(I).EQ.0.D+00) W(I) = 1.D-30
             !--------------------------------------------------------------------------------------------------------------
             13    CONTINUE
             RETURN
     END SUBROUTINE GAUCOF


     DOUBLE PRECISION FUNCTION PYTHAG(A,B)

             IMPLICIT NONE
             !----------------------------------------------------------------------------------------------------------------------
             !     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
             !----------------------------------------------------------------------------------------------------------------------
             REAL(8) :: ABSA,ABSB, A,B,PHYTAG
             ABSA=ABS(A)
             ABSB=ABS(B)
             IF(ABSA.GT.ABSB)THEN
                     PYTHAG=ABSA*SQRT(1.D+00+(ABSB/ABSA)**2)
             ELSE
                     IF(ABSB.EQ.0.D+00)THEN
                             PYTHAG=0.D+00
                     ELSE
                             PYTHAG=ABSB*SQRT(1.D+00+(ABSA/ABSB)**2)
                     ENDIF
             ENDIF
             RETURN
             END


     SUBROUTINE TQLI(D,E,N,NP,Z,IFAILTQL)
             IMPLICIT NONE
             !----------------------------------------------------------------------------------------------------------------------
             !     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
             !----------------------------------------------------------------------------------------------------------------------
             INTEGER :: I,N,NP,M,L,K,ITER,IFAILTQL
             REAL(8) :: D(NP),E(NP),Z(NP,NP),DD,G,R,S,C,P,F,B
             DO 11 I=2,N
             E(I-1)=E(I)
             11    CONTINUE
             E(N)=0.D+00
             DO 15 L=1,N
             ITER=0
             1       DO 12 M=L,N-1
             DD=ABS(D(M))+ABS(D(M+1))
             IF (ABS(E(M))+DD.EQ.DD) GOTO 2
             12      CONTINUE
             M=N
             2       IF(M.NE.L)THEN
                     IF(ITER.EQ.300) THEN
                             IFAILTQL = 1
                             RETURN
                     ENDIF
                     ITER=ITER+1
                     G=(D(L+1)-D(L))/(2.D+00*E(L))
                     R=PYTHAG(G,1.0D0)
                     G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
                     S=1.D+00
                     C=1.D+00
                     P=0.D+00
                     DO 14 I=M-1,L,-1
                     F=S*E(I)
                     B=C*E(I)
                     R=PYTHAG(F,G)
                     E(I+1)=R
                     IF(R.EQ.0.D0)THEN
                             D(I+1)=D(I+1)-P
                             E(M)=0.D+00
                             GOTO 1
                     ENDIF
                     S=F/R
                     C=G/R
                     G=D(I+1)-P
                     R=(D(I)-G)*S+2.D0*C*B
                     P=S*R
                     D(I+1)=G+P
                     G=C*R-B
                     !     OMIT LINES FROM HERE ...
                     DO 13 K=1,N
                     F=Z(K,I+1)
                     Z(K,I+1)=S*Z(K,I)+C*F
                     Z(K,I)=C*Z(K,I)-S*F
                     13          CONTINUE
                     !     ... TO HERE WHEN FINDING ONLY EIGENVALUES.
                     14        CONTINUE
                     D(L)=D(L)-P
                     E(L)=G
                     E(M)=0.D+00
                     GOTO 1
             ENDIF
             15    CONTINUE
             RETURN
     END SUBROUTINE TQLI


!      SUBROUTINE GET_KBARNIJ( IUPDATE, TK, PRES, DIAM, KBAR0IJ, KBAR3IJ )
!!-------------------------------------------------------------------------------------------------------------------------------------
!!     Routine to setup tables of mode-average coagulation coefficients
!!     KBAR0IJ and KBAR3IJ [m^3/s] for arbitrary temperature and pressure.
!!     Modes are assumed to be lognormal, and Sigmag values are set to
!!     constants for each mode, and Dg values are derived from the current
!!     value of the diameter of average mass for each mode, before being passed
!!     to this routine.
!!-------------------------------------------------------------------------------------------------------------------------------------
!      USE AERO_SETUP, ONLY: SIG0                        ! [um], [1], default lognormal parameters
!      IMPLICIT NONE                                     !            for each mode
!
!      ! Arguments.
!
!      INTEGER :: IUPDATE                                ! [1]  control flag
!      REAL(8) :: TK                                     ! [K]  ambient temperature
!      REAL(8) :: PRES                                   ! [Pa] ambient pressure
!      REAL(8) :: DIAM(NWEIGHTS)                         ! [um] geo. mean diameter for each mode
!      REAL(8) :: KBAR0IJ(NWEIGHTS,NWEIGHTS)             ! [m^3/s] 0th mode-average coag. coef.
!      REAL(8) :: KBAR3IJ(NWEIGHTS,NWEIGHTS)             ! [m^3/s] 3th mode-average coag. coef.
!
!      ! Local variables.
!
!      INTEGER       :: INDEX_DIAMI, INDEX_DIAMJ, INDEX_DIAMIP1, INDEX_DIAMJP1 
!      INTEGER       :: I, J, ITRANGE, IPRANGE                     
!      !------------------------------------------------------------------------------------------------------------
!      ! INDEX_SIGG(I) is the index of KIJ_SIGMAS to obtain Sigmag for mode I.
!      !------------------------------------------------------------------------------------------------------------
!      INTEGER, SAVE :: INDEX_SIGG(NWEIGHTS)                               ! [1]
!      REAL, SAVE :: DELTALNDG                                          ! [1] table spacing in ln(Dg)
!      REAL       :: KBAR0IJ_LL, KBAR0IJ_LU, KBAR0IJ_UL, KBAR0IJ_UU     ! [m^3/s] for bilinear interpolation
!      REAL       :: KBAR3IJ_LL, KBAR3IJ_LU, KBAR3IJ_UL, KBAR3IJ_UU     ! [m^3/s] for bilinear interpolation
!      REAL       :: XINTERPI, XINTERPJ, XINTERPT, XINTERPP             ! [1] for bilinear interpolation
!      REAL       :: TMP0, TMP3                                         ! [m^3/s] scratch variables
!      REAL       :: TPINTERP_LL, TPINTERP_LU, TPINTERP_UL, TPINTERP_UU ! [1] scratch variables
!      REAL       :: TUSE                                               ! [K]  ambient temperature local variable 
!      REAL       :: PUSE                                               ! [Pa] ambient pressure    local variable 
!      REAL       :: KBAR0IJ_LL_LL, KBAR0IJ_LU_LL, KBAR0IJ_UL_LL, KBAR0IJ_UU_LL  ! [m^3/s] for bilinear interp.
!      REAL       :: KBAR3IJ_LL_LL, KBAR3IJ_LU_LL, KBAR3IJ_UL_LL, KBAR3IJ_UU_LL  ! [m^3/s] for bilinear interp.
!      REAL       :: KBAR0IJ_LL_LU, KBAR0IJ_LU_LU, KBAR0IJ_UL_LU, KBAR0IJ_UU_LU  ! [m^3/s] for bilinear interp.
!      REAL       :: KBAR3IJ_LL_LU, KBAR3IJ_LU_LU, KBAR3IJ_UL_LU, KBAR3IJ_UU_LU  ! [m^3/s] for bilinear interp.
!      REAL       :: KBAR0IJ_LL_UL, KBAR0IJ_LU_UL, KBAR0IJ_UL_UL, KBAR0IJ_UU_UL  ! [m^3/s] for bilinear interp.
!      REAL       :: KBAR3IJ_LL_UL, KBAR3IJ_LU_UL, KBAR3IJ_UL_UL, KBAR3IJ_UU_UL  ! [m^3/s] for bilinear interp.
!      REAL       :: KBAR0IJ_LL_UU, KBAR0IJ_LU_UU, KBAR0IJ_UL_UU, KBAR0IJ_UU_UU  ! [m^3/s] for bilinear interp.
!      REAL       :: KBAR3IJ_LL_UU, KBAR3IJ_LU_UU, KBAR3IJ_UL_UU, KBAR3IJ_UU_UU  ! [m^3/s] for bilinear interp.
!      LOGICAL, SAVE :: FIRSTIME = .TRUE.
!      LOGICAL       :: FLAG
!
!      integer :: arrindex_diam(NWEIGHTS)
!      real(4) :: arrxinterp(NWEIGHTS)
!
!      IF( FIRSTIME ) THEN
!        FIRSTIME = .FALSE.
!        !------------------------------------------------------------------------------------------------------------
!        ! DELTALNG is the table spacing in ln(Dg) and needed for interpolation.
!        !------------------------------------------------------------------------------------------------------------
!        DELTALNDG = LOG( KIJ_DGMAX / KIJ_DGMIN ) / REAL( KIJ_NDGS - 1 )  ! to interpolate in Dg
!        !------------------------------------------------------------------------------------------------------------
!        ! To efficiently identify the Sigmag value assigned to each mode.
!        !------------------------------------------------------------------------------------------------------------
!        INDEX_SIGG(:) = 1
!        DO I=1, NWEIGHTS 
!          DO J=1, KIJ_NSGS
!            IF( ABS( SIG0(I) - DBLE( KIJ_SIGMAS(J) ) ) .LT. 1.0D-03 ) INDEX_SIGG(I) = J
!          ENDDO
!          ! WRITE(AUNIT2,'(I6,F12.6,I6,F12.6)') I, SIG0(I), INDEX_SIGG(I), DIAM(I)
!        ENDDO
!      ENDIF
!
!      !--------------------------------------------------------------------------------------------------------------
!      ! Temporary code to check the new value of KIJ_DGMIN. 
!      !--------------------------------------------------------------------------------------------------------------
!      ! IF( MINVAL( DIAM(:) ) .LT. KIJ_DGMIN ) THEN
!      !  WRITE(*,*)'MINVAL( DIAM(:) ) .LT. KIJ_DGMIN  in subr. GET_KBARNIJ. :', MINVAL( DIAM(:) ) 
!      !  STOP
!      ! ENDIF
!      !--------------------------------------------------------------------------------------------------------------
!
!      ! precompute common elements
!      DO I=1, NWEIGHTS
!        !----------------------------------------------------------------------------------------------------------
!        ! For mode I, get the lower and upper bounding table diameters and the interpolation variable XINTERPI.
!        !----------------------------------------------------------------------------------------------------------
!        INDEX_DIAMI = ( log( DIAM(I) / KIJ_DGMIN ) / DELTALNDG ) + 1
!        INDEX_DIAMI = MIN( MAX( INDEX_DIAMI, 1 ), KIJ_NDGS-1 )
!        arrindex_diam(I) = INDEX_DIAMI
!
!        XINTERPI = log( DIAM(I) / KIJ_DIAMETERS(INDEX_DIAMI) ) / DELTALNDG 
!        arrxinterp(I) = XINTERPI
!      END DO
!      !--------------------------------------------------------------------------------------------------------------
!      ! IUPDATE .EQ. 0: The mode-average coagulation coefficients are held constant throughout the simulation.
!      !--------------------------------------------------------------------------------------------------------------
!      IF( IUPDATE .EQ. 0 ) THEN
!        FLAG = .FALSE.
!        DO I=1, NWEIGHTS
!          !----------------------------------------------------------------------------------------------------------
!          ! For mode I, get the lower and upper bounding table diameters and the interpolation variable XINTERPI.
!          !----------------------------------------------------------------------------------------------------------
!          INDEX_DIAMI = arrindex_diam(I)
!          INDEX_DIAMIP1 = INDEX_DIAMI+1
!          IF(DIAM(I) .LT. KIJ_DIAMETERS(INDEX_DIAMI  )) FLAG = .TRUE.
!          IF(DIAM(I) .GT. KIJ_DIAMETERS(INDEX_DIAMIP1)) FLAG = .TRUE.  
!
!          IF( FLAG ) THEN
!            WRITE(*,*)'Problem in GET_KBARNIJ for IUPDATE = 0'
!            WRITE(*,'(2I6,8F11.5)')I,KIJ_DIAMETERS(INDEX_DIAMI),DIAM(I),KIJ_DIAMETERS(INDEX_DIAMIP1)
!            STOP
!          ENDIF
!!--------------------------------------------------------------------------------------------------------------------
!!           The lower and upper table diameters bounding DIAM(I)), and the interpolation variables
!!           XINTERPI were checked and found correct.
!!--------------------------------------------------------------------------------------------------------------------
!!           WRITE(AUNIT2,'(2I6,8F11.5)')I,KIJ_DIAMETERS(INDEX_DIAMI),DIAM(I),KIJ_DIAMETERS(INDEX_DIAMIP1),XINTERPI
!!--------------------------------------------------------------------------------------------------------------------
!        END DO
!
!        DO I=1, NWEIGHTS
!          !----------------------------------------------------------------------------------------------------------
!          ! For mode I, get the lower and upper bounding table diameters and the interpolation variable XINTERPI.
!          !----------------------------------------------------------------------------------------------------------
!          INDEX_DIAMI = arrindex_diam(I)
!          INDEX_DIAMIP1 = INDEX_DIAMI+1
!          XINTERPI = arrxinterp(I)
!
!          DO J=1, NWEIGHTS
!            !--------------------------------------------------------------------------------------------------------
!            ! For mode J, get the lower and upper bounding table diameters and the interpolation variable XINTERPJ.
!            !--------------------------------------------------------------------------------------------------------
!            INDEX_DIAMJ = arrindex_diam(J)
!            INDEX_DIAMJP1 = INDEX_DIAMJ+1
!            XINTERPJ = arrxinterp(J)
!            !--------------------------------------------------------------------------------------------------------
!            ! For each of the four points needed for bilinear interpolation, get the mode-average coagulation
!            ! coefficients at the selected temperature and pressure.
!            !--------------------------------------------------------------------------------------------------------
!            KBAR0IJ_LL = K0IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!            KBAR0IJ_LU = K0IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!            KBAR0IJ_UL = K0IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!            KBAR0IJ_UU = K0IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!            KBAR3IJ_LL = K3IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!            KBAR3IJ_LU = K3IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!            KBAR3IJ_UL = K3IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!            KBAR3IJ_UU = K3IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!            !--------------------------------------------------------------------------------------------------------
!            ! Interpolate in Dg(I) and Dg(J) for modes I and J.
!            !
!            ! When DIAM(I) = KIJ_DIAMETERS(INDEX_DIAMI), the lower I-mode Dg value, XINTERPI = 0.0, so
!            ! KBAR0IJ_LL and KBAR0IJ_LU should be multiplied by (1.0 - XINTERPI ) = 1.0.
!            !--------------------------------------------------------------------------------------------------------
!            TMP0         = KBAR0IJ_LL*(1.0-XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR0IJ_LU*(1.0-XINTERPI)*(    XINTERPJ)
!     &                   + KBAR0IJ_UL*(    XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR0IJ_UU*(    XINTERPI)*(    XINTERPJ)
!            TMP3         = KBAR3IJ_LL*(1.0-XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR3IJ_LU*(1.0-XINTERPI)*(    XINTERPJ)
!     &                   + KBAR3IJ_UL*(    XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR3IJ_UU*(    XINTERPI)*(    XINTERPJ)
!            KBAR0IJ(I,J) = DBLE( TMP0 )
!            KBAR3IJ(I,J) = DBLE( TMP3 )
!            !--------------------------------------------------------------------------------------------------------
!            ! For narrow distributions, KBAR0IJ and KBAR3IJ should be nearly equal, and were found to be so
!            ! with Sigmag = 1.1 for all modes.
!            !--------------------------------------------------------------------------------------------------------
!            ! WRITE(AUNIT2,'(2I6,2E15.5)')I,J,KBAR0IJ(I,J),KBAR3IJ(I,J)
!            !--------------------------------------------------------------------------------------------------------
!          ENDDO
!        ENDDO
!      !--------------------------------------------------------------------------------------------------------------
!      ! IUPDATE .EQ. 1: The mode-average coagulation coefficients are updated at each time step.
!      !--------------------------------------------------------------------------------------------------------------
!      ELSEIF( IUPDATE .EQ. 1 ) THEN
!        TUSE = MIN( MAX( REAL(TK),   KIJ_TEMP3 ), KIJ_TEMP1 )   ! Tmin=KIJ_TEMP3, Tmax=KIJ_TEMP1  
!        PUSE = MIN( MAX( REAL(PRES), KIJ_PRES3 ), KIJ_PRES1 )   ! pmin=KIJ_PRES3, pmax=KIJ_PRES1    
!        IF( TUSE .GT. KIJ_TEMP2 ) THEN                          ! Tmiddle value=KIJ_TEMP2
!          ITRANGE = 12                                          ! use temperatures 1 and 2
!          XINTERPT = ( TUSE - KIJ_TEMP2 ) / ( KIJ_TEMP1 - KIJ_TEMP2 ) 
!        ELSE
!          ITRANGE = 23                                          ! use temperatures 2 and 3
!          XINTERPT = ( TUSE - KIJ_TEMP3 ) / ( KIJ_TEMP2 - KIJ_TEMP3 ) 
!        ENDIF
!        IF( PUSE .GT. KIJ_PRES2 ) THEN                          ! pmiddle value=KIJ_PRES2
!          IPRANGE = 12                                          ! use pressures 1 and 2
!          XINTERPP = ( PUSE - KIJ_PRES2 ) / ( KIJ_PRES1 - KIJ_PRES2 ) 
!        ELSE
!          IPRANGE = 23                                          ! use pressures 2 and 3
!          XINTERPP = ( PUSE - KIJ_PRES3 ) / ( KIJ_PRES2 - KIJ_PRES3 ) 
!        ENDIF
!        TPINTERP_LL = (1.0-XINTERPT)*(1.0-XINTERPP)             ! all weight at T_lower, p_lower
!        TPINTERP_LU = (1.0-XINTERPT)*(    XINTERPP)             ! all weight at T_lower, p_upper
!        TPINTERP_UL = (    XINTERPT)*(1.0-XINTERPP)             ! all weight at T_upper, p_lower
!        TPINTERP_UU = 1.0-TPINTERP_LL-TPINTERP_LU-TPINTERP_UL   ! all weight at T_upper, p_upper
!!--------------------------------------------------------------------------------------------------------------------
!!       WRITE(AUNIT2,'(/A/)')'new step'
!!       WRITE(AUNIT2,'(A40,2F15.6    )')'TUSE, PUSE = ', TUSE, PUSE
!!       WRITE(AUNIT2,'(A40,2I4,2F13.7)')'ITRANGE, IPRANGE, XINTERPT, XINTERPP = ',
!!    &                                   ITRANGE, IPRANGE, XINTERPT, XINTERPP
!!       WRITE(AUNIT2,'(A40,4F12.6    )')'TPINTERP_LL, TPINTERP_LU, TPINTERP_UL, TPINTERP_UU = ', 
!!    &                                   TPINTERP_LL, TPINTERP_LU, TPINTERP_UL, TPINTERP_UU
!!--------------------------------------------------------------------------------------------------------------------
!        DO I=1, NWEIGHTS
!          !----------------------------------------------------------------------------------------------------------
!          ! For mode I, get the lower and upper bounding table diameters and the interpolation variable XINTERPI.
!          !----------------------------------------------------------------------------------------------------------
!          INDEX_DIAMI = arrindex_diam(I)
!          INDEX_DIAMIP1 = INDEX_DIAMI+1
!          XINTERPI = arrxinterp(I)
!          DO J=1, NWEIGHTS
!
!            if (CITABLE(I,J) == 'OFF') then ! Turn off coagulation between selected modes.
!              KBAR0IJ(I,J) = 1.0D-30
!              KBAR3IJ(I,J) = 1.0D-30
!              cycle
!            end if
!
!            !--------------------------------------------------------------------------------------------------------
!            ! For mode J, get the lower and upper bounding table diameters and the interpolation variable XINTERPJ.
!            !--------------------------------------------------------------------------------------------------------
!            INDEX_DIAMJ = arrindex_diam(J)
!            INDEX_DIAMJP1 = INDEX_DIAMJ+1
!            XINTERPJ = arrxinterp(J)
!            !--------------------------------------------------------------------------------------------------------
!            ! For each of the four points needed for bilinear interpolation in Dg(I) and Dg(J), get the 
!            ! mode-average coagulation coefficients at each of the four temperature-pressure points.
!            !
!            ! In KBAR0IJ_AB_CD, A indicates the upper or lower value of Dg(I)
!            !                   B indicates the upper or lower value of Dg(J)
!            !                   C indicates the upper or lower value of temperature: T1 > T2 > T3
!            !                   D indicates the upper or lower value of pressure:    p1 > p2 > p3
!            !--------------------------------------------------------------------------------------------------------
!            IF( ITRANGE .EQ. 12 .AND. IPRANGE .EQ. 12 ) THEN
!
!              KBAR0IJ_LL_LL = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LL = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LL = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LL = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LL = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LL = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LL = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LL = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_LU = K0IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LU = K0IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LU = K0IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LU = K0IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LU = K3IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LU = K3IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LU = K3IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LU = K3IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UL = K0IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UL = K0IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UL = K0IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UL = K0IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UL = K3IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UL = K3IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UL = K3IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UL = K3IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UU = K0IJ_TEMP1PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UU = K0IJ_TEMP1PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UU = K0IJ_TEMP1PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UU = K0IJ_TEMP1PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UU = K3IJ_TEMP1PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UU = K3IJ_TEMP1PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UU = K3IJ_TEMP1PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UU = K3IJ_TEMP1PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!            ELSEIF( ITRANGE .EQ. 12 .AND. IPRANGE .EQ. 23 ) THEN
!
!              KBAR0IJ_LL_LL = K0IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LL = K0IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LL = K0IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LL = K0IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LL = K3IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LL = K3IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LL = K3IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LL = K3IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_LU = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LU = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LU = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LU = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LU = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LU = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LU = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LU = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UL = K0IJ_TEMP1PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UL = K0IJ_TEMP1PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UL = K0IJ_TEMP1PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UL = K0IJ_TEMP1PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UL = K3IJ_TEMP1PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UL = K3IJ_TEMP1PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UL = K3IJ_TEMP1PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UL = K3IJ_TEMP1PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UU = K0IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UU = K0IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UU = K0IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UU = K0IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UU = K3IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UU = K3IJ_TEMP1PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UU = K3IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UU = K3IJ_TEMP1PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!            ELSEIF( ITRANGE .EQ. 23 .AND. IPRANGE .EQ. 12 ) THEN
!
!              KBAR0IJ_LL_LL = K0IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LL = K0IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LL = K0IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LL = K0IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LL = K3IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LL = K3IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LL = K3IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LL = K3IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_LU = K0IJ_TEMP3PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LU = K0IJ_TEMP3PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LU = K0IJ_TEMP3PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LU = K0IJ_TEMP3PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LU = K3IJ_TEMP3PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LU = K3IJ_TEMP3PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LU = K3IJ_TEMP3PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LU = K3IJ_TEMP3PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UL = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UL = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UL = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UL = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UL = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UL = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UL = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UL = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UU = K0IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UU = K0IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UU = K0IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UU = K0IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UU = K3IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UU = K3IJ_TEMP2PRES1( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UU = K3IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UU = K3IJ_TEMP2PRES1( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!            ELSEIF( ITRANGE .EQ. 23 .AND. IPRANGE .EQ. 23 ) THEN
!
!              KBAR0IJ_LL_LL = K0IJ_TEMP3PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LL = K0IJ_TEMP3PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LL = K0IJ_TEMP3PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LL = K0IJ_TEMP3PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LL = K3IJ_TEMP3PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LL = K3IJ_TEMP3PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LL = K3IJ_TEMP3PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LL = K3IJ_TEMP3PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_LU = K0IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_LU = K0IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_LU = K0IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_LU = K0IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_LU = K3IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_LU = K3IJ_TEMP3PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_LU = K3IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_LU = K3IJ_TEMP3PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UL = K0IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UL = K0IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UL = K0IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UL = K0IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UL = K3IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UL = K3IJ_TEMP2PRES3( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UL = K3IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UL = K3IJ_TEMP2PRES3( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!              KBAR0IJ_LL_UU = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_LU_UU = K0IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR0IJ_UL_UU = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR0IJ_UU_UU = K0IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_LL_UU = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_LU_UU = K3IJ_TEMP2PRES2( INDEX_DIAMI,  INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!              KBAR3IJ_UL_UU = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJ,  INDEX_SIGG(J) )
!              KBAR3IJ_UU_UU = K3IJ_TEMP2PRES2( INDEX_DIAMIP1,INDEX_SIGG(I),INDEX_DIAMJP1,INDEX_SIGG(J) )
!
!            ELSE
!
!              WRITE(*,*)'Error in GET_KBARNIJ: ITRANGE, IPRANGE = ', ITRANGE, IPRANGE
!              STOP
!
!            ENDIF
!            !--------------------------------------------------------------------------------------------------------
!            ! Interpolate in T and p for each of the four points needed for the Dg(I) and Dg(J) interpolation.
!            !--------------------------------------------------------------------------------------------------------
!            KBAR0IJ_LL = TPINTERP_LL*KBAR0IJ_LL_LL + TPINTERP_LU*KBAR0IJ_LL_LU
!     &                 + TPINTERP_UL*KBAR0IJ_LL_UL + TPINTERP_UU*KBAR0IJ_LL_UU
!            KBAR0IJ_LU = TPINTERP_LL*KBAR0IJ_LU_LL + TPINTERP_LU*KBAR0IJ_LU_LU
!     &                 + TPINTERP_UL*KBAR0IJ_LU_UL + TPINTERP_UU*KBAR0IJ_LU_UU
!            KBAR0IJ_UL = TPINTERP_LL*KBAR0IJ_UL_LL + TPINTERP_LU*KBAR0IJ_UL_LU
!     &                 + TPINTERP_UL*KBAR0IJ_UL_UL + TPINTERP_UU*KBAR0IJ_UL_UU
!            KBAR0IJ_UU = TPINTERP_LL*KBAR0IJ_UU_LL + TPINTERP_LU*KBAR0IJ_UU_LU
!     &                 + TPINTERP_UL*KBAR0IJ_UU_UL + TPINTERP_UU*KBAR0IJ_UU_UU
!
!            KBAR3IJ_LL = TPINTERP_LL*KBAR3IJ_LL_LL + TPINTERP_LU*KBAR3IJ_LL_LU
!     &                 + TPINTERP_UL*KBAR3IJ_LL_UL + TPINTERP_UU*KBAR3IJ_LL_UU
!            KBAR3IJ_LU = TPINTERP_LL*KBAR3IJ_LU_LL + TPINTERP_LU*KBAR3IJ_LU_LU
!     &                 + TPINTERP_UL*KBAR3IJ_LU_UL + TPINTERP_UU*KBAR3IJ_LU_UU
!            KBAR3IJ_UL = TPINTERP_LL*KBAR3IJ_UL_LL + TPINTERP_LU*KBAR3IJ_UL_LU
!     &                 + TPINTERP_UL*KBAR3IJ_UL_UL + TPINTERP_UU*KBAR3IJ_UL_UU
!            KBAR3IJ_UU = TPINTERP_LL*KBAR3IJ_UU_LL + TPINTERP_LU*KBAR3IJ_UU_LU
!     &                 + TPINTERP_UL*KBAR3IJ_UU_UL + TPINTERP_UU*KBAR3IJ_UU_UU
!!--------------------------------------------------------------------------------------------------------------------
!!           WRITE(AUNIT2,'(A40,4E13.5)')'KBAR0IJ_LL, KBAR0IJ_LU, KBAR0IJ_UL, KBAR0IJ_UU = ',
!!    &                                   KBAR0IJ_LL, KBAR0IJ_LU, KBAR0IJ_UL, KBAR0IJ_UU
!!           WRITE(AUNIT2,'(A40,4E13.5)')'KBAR3IJ_LL, KBAR3IJ_LU, KBAR3IJ_UL, KBAR3IJ_UU = ',
!!    &                                   KBAR3IJ_LL, KBAR3IJ_LU, KBAR3IJ_UL, KBAR3IJ_UU
!!--------------------------------------------------------------------------------------------------------------------
!            ! Interpolate in Dg(I) and Dg(J) for modes I and J.
!            !--------------------------------------------------------------------------------------------------------
!            TMP0         = KBAR0IJ_LL*(1.0-XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR0IJ_LU*(1.0-XINTERPI)*(    XINTERPJ)
!     &                   + KBAR0IJ_UL*(    XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR0IJ_UU*(    XINTERPI)*(    XINTERPJ)
!            TMP3         = KBAR3IJ_LL*(1.0-XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR3IJ_LU*(1.0-XINTERPI)*(    XINTERPJ)
!     &                   + KBAR3IJ_UL*(    XINTERPI)*(1.0-XINTERPJ)
!     &                   + KBAR3IJ_UU*(    XINTERPI)*(    XINTERPJ)
!            KBAR0IJ(I,J) = DBLE( TMP0 )
!            KBAR3IJ(I,J) = DBLE( TMP3 )
!            !--------------------------------------------------------------------------------------------------------
!            ! For narrow distributions, KBAR0IJ and KBAR3IJ should be nearly equal, and were found to be so
!            ! with Sigmag = 1.1 for all modes.
!            !--------------------------------------------------------------------------------------------------------
!            ! WRITE(AUNIT2,'(A40,2I6,4E15.5)')'I,J,KBAR0IJ,KBAR3IJ=',I,J,KBAR0IJ(I,J),KBAR3IJ(I,J)
!            !--------------------------------------------------------------------------------------------------------
!          ENDDO
!        ENDDO
!      ENDIF
!
!      RETURN
!      END SUBROUTINE GET_KBARNIJ
!

      SUBROUTINE BROWNIAN_COAG_COEF( DI, DJ, TEMPK, PRES, BETAIJ )
!-------------------------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the Brownian coagulation coefficient BETAIJ 
!     for particles of diameters DI and DJ at ambient temperature TEMPK
!     and pressure PRES, using the interpolation formula of Fuchs (1964),
!     as given in Jacobson 1999, p.446, eq.(16.28) 
!              or Jacobson 2005, p.509, eq.(15.33). 
!-------------------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      REAL,intent(in) :: DI, DJ          ! particle diameters [um]
      REAL,intent(in) :: TEMPK           ! ambient temperature [K]
      REAL,intent(in) :: PRES            ! ambient pressure [Pa]
      REAL,intent(out) :: BETAIJ          ! coagulation coefficient [m^3/s/particle]

      REAL, PARAMETER :: PI = 3.141592653589793
      REAL, PARAMETER :: FOURPI = 4.0 * PI
      REAL, PARAMETER :: SIXPI  = 6.0 * PI
      REAL, PARAMETER :: FOURPIBY3 = 4.0 * PI / 3.0
      REAL, PARAMETER :: KB = 1.38065E-23     ! Boltzmann's constant [J/K]
      REAL, PARAMETER :: NA = 6.0221367E+23   ! Avogadro's number [#/mole]
      REAL, PARAMETER :: MWAIR = 28.9628E-03  ! molecular weight of air in [kg/mole]
      REAL, PARAMETER :: M1AIR = MWAIR/NA     ! average mass of one air molecule [kg]
      REAL, PARAMETER :: VACONST = 8.0*KB/PI/M1AIR  ! used in eq.(16.21) for VABAR
      REAL, PARAMETER :: APRIME = 1.249       ! for Cunningham slip-flow correction
      REAL, PARAMETER :: BPRIME = 0.42        !   as given in eq.(16.25);
      REAL, PARAMETER :: CPRIME = 0.87        !   Kasten (1968) values used here
      REAL, PARAMETER :: VPCONST = 8.0*KB/PI  ! used in eq.(16.27) for VPIBAR [J/K]
      REAL, PARAMETER :: DENSP_KGM3 = 1.4E+03 ! particle density [kg/m^3]

      REAL :: RI, RJ          ! particle radii [m]
      REAL :: DPI, DPJ        ! particle diffusion coefficients [m^2/s]
      REAL :: GI, GJ          ! Cunningham slip-flow corrections for particles I, J [1]
      REAL :: ETAA            ! dynamic viscosity of air [kg/m/s]
      REAL :: VABAR           ! mean thermal velocity of an air molecule [m/s]
      REAL :: RHOA            ! density of air [kg/m^3]
      REAL :: LAMBDAA         ! mean free path of air [m]
      REAL :: KNAI, KNAJ      ! Knudsen numbers in air for particles I, J [1]
      REAL :: VPIBAR, VPJBAR  ! mean thermal velocity  for particles I, J [m/s]
      REAL :: MPI, MPJ        ! particle masses [kg]                       
      REAL :: VOLI, VOLJ      ! particle volumes [m^3]
      REAL :: LAMBDAPI        ! mean free path for particle I [m]
      REAL :: LAMBDAPJ        ! mean free path for particle J [m]
      REAL :: DELTAI, DELTAJ  ! see eq.(16.29), p. 446 of MZJ 1999
      REAL :: DENOM1, DENOM2  ! scratch variables
!     REAL :: DUM1, DUM2      ! scratch variables

      RI = 0.5E-06 * DI           ! convert from [um] to [m]
      RJ = 0.5E-06 * DJ           ! convert from [um] to [m]
      VOLI = FOURPIBY3 * RI**3    ! volume of particle I [m^3]
      VOLJ = FOURPIBY3 * RJ**3    ! volume of particle J [m^3]
      !-------------------------------------------------------------------------
      ! Dynamic viscosity of air [kg/m/s], Jacobson, 2005, eq.(4.54).
      ! Mean thermal velocity of an air molecule [m/s].
      ! Density of air [kg/m^3] from the ideal gas law.
      ! Mean free path of an air molecule [m], Jacobson, 2005, eq.(15.24).
      !-------------------------------------------------------------------------
      ETAA = 1.8325E-05 * (416.16/(TEMPK + 120.0)) * (TEMPK/296.16)**1.5
      VABAR = SQRT( VACONST * TEMPK ) 
      RHOA = PRES * M1AIR / ( KB * TEMPK ) 
      LAMBDAA = 2.0 * ETAA / ( RHOA * VABAR ) 
      KNAI = LAMBDAA / RI
      KNAJ = LAMBDAA / RJ
!-------------------------------------------------------------------------------
!     DUM1 = 0.0
!     DUM2 = 0.0
!     WRITE(*,     '(12D11.3)') DUM1, DUM2
!-------------------------------------------------------------------------------
!     DUM1 = KNAI * ( APRIME + BPRIME * EXP(-CPRIME/KNAI) )
!     DUM2 = KNAJ * ( APRIME + BPRIME * EXP(-CPRIME/KNAJ) )
!-------------------------------------------------------------------------------
!     WRITE(*,     '(12D11.3)') DUM1, DUM2
!-------------------------------------------------------------------------------
!     GI = 1.0 + DUM1
!     GJ = 1.0 + DUM2
!-------------------------------------------------------------------------------
!     WRITE(*,     '(12D11.3)') DI,DJ,KNAI,KNAJ,GI,GJ
!-------------------------------------------------------------------------------
      GI = 1.0 + KNAI * ( APRIME + BPRIME * EXP(-CPRIME/KNAI) )
      GJ = 1.0 + KNAJ * ( APRIME + BPRIME * EXP(-CPRIME/KNAJ) )
      DPI = KB * TEMPK * GI / ( SIXPI * RI * ETAA ) 
      DPJ = KB * TEMPK * GJ / ( SIXPI * RJ * ETAA ) 
      MPI = VOLI * DENSP_KGM3
      MPJ = VOLJ * DENSP_KGM3
      VPIBAR = SQRT ( VPCONST * TEMPK / MPI )
      VPJBAR = SQRT ( VPCONST * TEMPK / MPJ )
      LAMBDAPI = 2.0 * DPI / ( PI * VPIBAR ) 
      LAMBDAPJ = 2.0 * DPJ / ( PI * VPJBAR ) 
      DELTAI = (2.0*RI + LAMBDAPI)**3 - (4.0*RI*RI + LAMBDAPI*LAMBDAPI)**1.5
      DELTAI = DELTAI / ( 6.0*RI*LAMBDAPI ) - 2.0*RI  
      DELTAJ = (2.0*RJ + LAMBDAPJ)**3 - (4.0*RJ*RJ + LAMBDAPJ*LAMBDAPJ)**1.5
      DELTAJ = DELTAJ / ( 6.0*RJ*LAMBDAPJ ) - 2.0*RJ  
      DENOM1 = (RI+RJ) / ( (RI+RJ) + SQRT(DELTAI*DELTAI + DELTAJ*DELTAJ) )
      DENOM2 = 4.0*(DPI+DPJ) / ( (RI+RJ) * SQRT(VPIBAR*VPIBAR + VPJBAR*VPJBAR) )
      BETAIJ = FOURPI * ( RI + RJ ) * ( DPI + DPJ ) / ( DENOM1 + DENOM2 )
!-------------------------------------------------------------------------------
!     WRITE(AUNIT2,'(12D11.3)') PRES,RHOA,LAMBDAA,DI,DJ,KNAI,KNAJ,GI,GJ,DPI,DPJ,BETAIJ
!     WRITE(*,     '(12D11.3)') PRES,RHOA,LAMBDAA,DI,DJ,KNAI,KNAJ,GI,GJ,DPI,DPJ,BETAIJ
!     WRITE(AUNIT2,*)'  '
!     WRITE(AUNIT2,*)'DI,DJ,TEMPK,PRES'
!     WRITE(AUNIT2,*) DI,DJ,TEMPK,PRES
!     WRITE(AUNIT2,*)'  '
!     WRITE(AUNIT2,*)'ETAA,VABAR,RHOA,LAMBDAA'
!     WRITE(AUNIT2,*) ETAA,VABAR,RHOA,LAMBDAA 
!     WRITE(AUNIT2,*)'  '
!     WRITE(AUNIT2,*)'KNAI,KNAJ,GI,GJ,DPI,DPJ'
!     WRITE(AUNIT2,*) KNAI,KNAJ,GI,GJ,DPI,DPJ 
!     WRITE(AUNIT2,*)'  '
!     WRITE(AUNIT2,*)'MPI,MPJ,VOLI,VOLJ,VPIBAR,VPJBAR'
!     WRITE(AUNIT2,*) MPI,MPJ,VOLI,VOLJ,VPIBAR,VPJBAR 
!     WRITE(AUNIT2,*)'  '
!     WRITE(AUNIT2,*)'LAMBDAPI,LAMBDAPJ,DELTAI,DELTAJ,DENOM1,DENOM2'
!     WRITE(AUNIT2,*) LAMBDAPI,LAMBDAPJ,DELTAI,DELTAJ,DENOM1,DENOM2 
!     WRITE(AUNIT2,*)'  '
!     WRITE(AUNIT2,*)'BETAIJ=',BETAIJ
!-------------------------------------------------------------------------------
      END SUBROUTINE BROWNIAN_COAG_COEF


      SUBROUTINE CBDE_COAG_COEF( DI, DJ, TEMPK, PRES, KDEIJ )
!-------------------------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the convective Brownian diffusion enhancement coagulation coefficient.
!     See MZJ, 2005, p. 510, Eq.15.35.
!-------------------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL,intent(in) :: DI, DJ                          ! particle diameters [um]
      REAL,intent(in) :: TEMPK                           ! ambient temperature [K]
      REAL,intent(in) :: PRES                            ! ambient pressure [Pa]
      REAL,intent(out) :: KDEIJ                           ! coagulation coefficient [m^3/s/particle]
      REAL, PARAMETER :: PI = 3.141592653589793
      REAL, PARAMETER :: SIXPI  = 6.0 * PI
      REAL, PARAMETER :: KB = 1.38065E-23     ! Boltzmann's constant [J/K]
      REAL, PARAMETER :: NA = 6.0221367E+23   ! Avogadro's number [#/mole]
      REAL, PARAMETER :: MWAIR = 28.9628E-03  ! molecular weight of air in [kg/mole]
      REAL, PARAMETER :: M1AIR = MWAIR/NA     ! average mass of one air molecule [kg]
      REAL, PARAMETER :: VACONST = 8.0*KB/PI/M1AIR  ! used in eq.(16.21) for VABAR
      REAL, PARAMETER :: APRIME = 1.249       ! for Cunningham slip-flow correction
      REAL, PARAMETER :: BPRIME = 0.42        !   as given in eq.(16.25);
      REAL, PARAMETER :: CPRIME = 0.87        !   Kasten (1968) values used here
      REAL, PARAMETER :: DENSP_KGM3 = 1.4E+03 ! particle density [kg/m^3]
      REAL, PARAMETER :: GGRAV = 9.81         ! gravitational acceleration [m/s^2]
      REAL :: RI, RJ                          ! particle radii [um]
      REAL :: DPI                             ! particle diffusion coefficients [m^2/s]
      REAL :: GI, GJ                          ! Cunningham slip-flow corrections for particles I, J [1]
      REAL :: ETAA                            ! dynamic viscosity of air [kg/m/s]
      REAL :: VABAR                           ! mean thermal velocity of an air molecule [m/s]
      REAL :: RHOA                            ! density of air [kg/m^3]
      REAL :: LAMBDAA                         ! mean free path of air [m]
      REAL :: KNAI, KNAJ                      ! Knudsen numbers in air for particles I, J [1]
      REAL :: NUA                             ! kinematic viscosity of air [m^2/s]
      REAL :: VFJ                             ! terminal fall speed [m/s]
      REAL :: REJ                             ! particle Reynolds number [1]
      REAL :: SCPI                            ! particle Schmidt  number [1]
      REAL :: KBIJ                            ! Brownian diffusion coefficient [m^3/s]

      !-----------------------------------------------------------------------------------------------------------------
      ! RJ must be larger than or equal to RI for this coagulation mechanism.
      ! Make RJ greater than or equal RI for all pairs of particles. 
      !-----------------------------------------------------------------------------------------------------------------
      IF( DJ .GE. DI ) THEN 
        RI = 0.5E-06 * DI                                                     ! convert from [um] to [m]
        RJ = 0.5E-06 * DJ                                                     ! convert from [um] to [m]
      ELSE
        RJ = 0.5E-06 * DI                                                     ! convert from [um] to [m]
        RI = 0.5E-06 * DJ                                                     ! convert from [um] to [m]
      ENDIF
      !-----------------------------------------------------------------------------------------------------------------
      ! Dynamic viscosity of air [kg/m/s], Jacobson, 2005, eq.(4.54).
      ! Mean thermal velocity of an air molecule [m/s].
      ! Density of air [kg/m^3] from the ideal gas law.
      ! Mean free path of an air molecule [m], Jacobson, 2005, eq.(15.24).
      !-----------------------------------------------------------------------------------------------------------------
      ETAA = 1.8325E-05 * (416.16/(TEMPK + 120.0)) * (TEMPK/296.16)**1.5      ! [kg/m/s]
      VABAR = SQRT( VACONST * TEMPK )                                         ! [m/s]
      RHOA = PRES * M1AIR / ( KB * TEMPK )                                    ! [kg/m^3]
      LAMBDAA = 2.0 * ETAA / ( RHOA * VABAR )                                 ! [m]
      KNAI = LAMBDAA / RI                                                     ! [1]
      KNAJ = LAMBDAA / RJ                                                     ! [1]
      GI = 1.0 + KNAI * ( APRIME + BPRIME * EXP(-CPRIME/KNAI) )               ! Cunningham slip-flow correction [1]
      GJ = 1.0 + KNAJ * ( APRIME + BPRIME * EXP(-CPRIME/KNAJ) )               ! Cunningham slip-flow correction [1]
      DPI = KB * TEMPK * GI / ( SIXPI * RI * ETAA )                           ! particle diffusion coefficient  [m^2/s]
      NUA = ETAA / RHOA                                                       ! kinematic viscosity of air [m^2/s]
      VFJ = 2.0 * RJ*RJ*(DENSP_KGM3-RHOA)*GGRAV*GJ/(9.0*ETAA)                 ! fall velocity [m/s] 
      REJ = 2.0 * RJ * VFJ / NUA                                              ! particle Reynolds number [1] 
      SCPI = NUA / DPI                                                        ! particle Schmidt  number [1] 
      CALL BROWNIAN_COAG_COEF ( DI, DJ, TEMPK, PRES, KBIJ )                   ! Brownian coag. coef. [m^3/s]
      IF( REJ .LE. 1.0 ) KDEIJ = 0.45 * KBIJ * REJ**0.3333333 * SCPI**0.3333333   ! CBDE coag. coef. [m^3/s]
      IF( REJ .GT. 1.0 ) KDEIJ = 0.45 * KBIJ * REJ**0.5000000 * SCPI**0.3333333   ! CBDE coag. coef. [m^3/s]
      ! WRITE(*,*) 'VFJ = ', VFJ
      END SUBROUTINE CBDE_COAG_COEF


      SUBROUTINE GRAVCOLL_COAG_COEF( DI, DJ, TEMPK, PRES, KGCIJ )
!-------------------------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the gravitational collection coagulation coefficient. 
!     See MZJ, 2005, p. 510, Eq.15.37.
!-------------------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL,intent(in) :: DI, DJ                          ! particle diameters [um]
      REAL,intent(in) :: TEMPK                           ! ambient temperature [K]
      REAL,intent(in) :: PRES                            ! ambient pressure [Pa]
      REAL,intent(out) :: KGCIJ                           ! coagulation coefficient [m^3/s/particle]
      REAL, PARAMETER :: PI = 3.141592653589793
      REAL, PARAMETER :: SIXPI  = 6.0 * PI
      REAL, PARAMETER :: KB = 1.38065E-23     ! Boltzmann's constant [J/K]
      REAL, PARAMETER :: NA = 6.0221367E+23   ! Avogadro's number [#/mole]
      REAL, PARAMETER :: MWAIR = 28.9628E-03  ! molecular weight of air in [kg/mole]
      REAL, PARAMETER :: M1AIR = MWAIR/NA     ! average mass of one air molecule [kg]
      REAL, PARAMETER :: VACONST = 8.0*KB/PI/M1AIR  ! used in eq.(16.21) for VABAR
      REAL, PARAMETER :: APRIME = 1.249       ! for Cunningham slip-flow correction
      REAL, PARAMETER :: BPRIME = 0.42        !   as given in eq.(16.25);
      REAL, PARAMETER :: CPRIME = 0.87        !   Kasten (1968) values used here
      REAL, PARAMETER :: DENSP_KGM3 = 1.4E+03 ! particle density [kg/m^3]
      REAL, PARAMETER :: GGRAV = 9.81         ! gravitational acceleration [m/s^2]
      REAL :: RI, RJ                          ! particle radii [um]
      REAL :: GI, GJ                          ! Cunningham slip-flow corrections for particles I, J [1]
      REAL :: ETAA                            ! dynamic viscosity of air [kg/m/s]
      REAL :: VABAR                           ! mean thermal velocity of an air molecule [m/s]
      REAL :: RHOA                            ! density of air [kg/m^3]
      REAL :: LAMBDAA                         ! mean free path of air [m]
      REAL :: KNAI, KNAJ                      ! Knudsen numbers in air for particles I, J [1]
      REAL :: VFI, VFJ                        ! terminal fall speed [m/s]
      REAL :: ECOLLIJ                         ! collision efficiency [1]
      REAL :: P                               ! p = min(ri,rj)/max(ri,rj) from Pruppacher and Klett 1980, p.377.

      RI = 0.5E-06 * DI                                                       ! convert from [um] to [m]
      RJ = 0.5E-06 * DJ                                                       ! convert from [um] to [m]
      !-----------------------------------------------------------------------------------------------------------------
      ! Dynamic viscosity of air [kg/m/s], Jacobson, 2005, eq.(4.54).
      ! Mean thermal velocity of an air molecule [m/s].
      ! Density of air [kg/m^3] from the ideal gas law.
      ! Mean free path of an air molecule [m], Jacobson, 2005, eq.(15.24).
      !-----------------------------------------------------------------------------------------------------------------
      ETAA = 1.8325E-05 * (416.16/(TEMPK + 120.0)) * (TEMPK/296.16)**1.5      ! [kg/m/s]
      VABAR = SQRT( VACONST * TEMPK )                                         ! [m/s]
      RHOA = PRES * M1AIR / ( KB * TEMPK )                                    ! [kg/m^3]
      LAMBDAA = 2.0 * ETAA / ( RHOA * VABAR )                                 ! [m]
      KNAI = LAMBDAA / RI                                                     ! [1]
      KNAJ = LAMBDAA / RJ                                                     ! [1]
      GI = 1.0 + KNAI * ( APRIME + BPRIME * EXP(-CPRIME/KNAI) )               ! Cunningham slip-flow correction [1]
      GJ = 1.0 + KNAJ * ( APRIME + BPRIME * EXP(-CPRIME/KNAJ) )               ! Cunningham slip-flow correction [1]
      VFI = 2.0 * RI*RI*(DENSP_KGM3-RHOA)*GGRAV*GI/(9.0*ETAA)                 ! fall velocity [m/s] 
      VFJ = 2.0 * RJ*RJ*(DENSP_KGM3-RHOA)*GGRAV*GJ/(9.0*ETAA)                 ! fall velocity [m/s] 
!-------------------------------------------------------------------------------------------------------------------------
!     Taking the collision efficiency from Pruppacher and Klett, 1980, Eq.12-78.
!-------------------------------------------------------------------------------------------------------------------------
      IF( RJ .GT. RI ) THEN
        P = RI / RJ
      ELSE
        P = RJ / RI
      ENDIF
      ECOLLIJ = 0.5 * ( P / ( 1.0 + P ) )**2                                  ! collision efficiency [1], Eq.12-78
      KGCIJ = ECOLLIJ * PI * ( RI + RJ )**2 * ABS( VFI - VFJ )                ! grav. collection coag. coef. [m^3/s]
!-------------------------------------------------------------------------------------------------------------------------
!     WRITE(*,'(/,A    )')'KGCIJ, RI, RJ, VFI, VFJ, ECOLLIJ, STIJ, REJ, EVIJ, EAIJ'
!     WRITE(*,'(10D12.4)') KGCIJ, RI, RJ, VFI, VFJ, ECOLLIJ, STIJ, REJ, EVIJ, EAIJ
!     WRITE(*,*) 'NUA = ', NUA
!-------------------------------------------------------------------------------------------------------------------------
      END SUBROUTINE GRAVCOLL_COAG_COEF


      SUBROUTINE TURB_COAG_COEF( DI, DJ, TEMPK, PRES, KTIIJ, KTSIJ )
!-------------------------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the turbulent inertial coagulation coefficient. See MZJ, 2005, p. 511, Eq.15.40.
!     Routine to calculate the turbulent shear    coagulation coefficient. See MZJ, 2005, p. 511, Eq.15.41.
!-------------------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL,intent(in) :: DI, DJ                          ! particle diameters [um]
      REAL,intent(in) :: TEMPK                           ! ambient temperature [K]
      REAL,intent(in) :: PRES                            ! ambient pressure [Pa]
      REAL,intent(out) :: KTIIJ                           ! turbulent inertial coagulation coefficient [m^3/s/particle]
      REAL,intent(out) :: KTSIJ                           ! turbulent shear    coagulation coefficient [m^3/s/particle]
      REAL, PARAMETER :: PI = 3.141592653589793
      REAL, PARAMETER :: KB = 1.38065E-23     ! Boltzmann's constant [J/K]
      REAL, PARAMETER :: NA = 6.0221367E+23   ! Avogadro's number [#/mole]
      REAL, PARAMETER :: MWAIR = 28.9628E-03  ! molecular weight of air in [kg/mole]
      REAL, PARAMETER :: M1AIR = MWAIR/NA     ! average mass of one air molecule [kg]
      REAL, PARAMETER :: VACONST = 8.0*KB/PI/M1AIR  ! used in eq.(16.21) for VABAR
      REAL, PARAMETER :: APRIME = 1.249       ! for Cunningham slip-flow correction
      REAL, PARAMETER :: BPRIME = 0.42        !   as given in eq.(16.25);
      REAL, PARAMETER :: CPRIME = 0.87        !   Kasten (1968) values used here
      REAL, PARAMETER :: DENSP_KGM3 = 1.4E+03 ! particle density [kg/m^3]
      REAL, PARAMETER :: GGRAV = 9.81         ! gravitational acceleration [m/s^2]
      REAL, PARAMETER :: ED = 5.0E-04         ! dissipation rate of turbulent kinetic energy per gram of medium [m^2/s^3]
                                                 ! Typical clear-sky value from MZJ (2005, p.511,p.236) of 5.0E-05 [m^2/s^3]
                                                 !   from Pruppacher and Klett (1997).
      REAL :: RI, RJ                          ! particle radii [um]
      REAL :: GI, GJ                          ! Cunningham slip-flow corrections for particles I, J [1]
      REAL :: ETAA                            ! dynamic viscosity of air [kg/m/s]
      REAL :: VABAR                           ! mean thermal velocity of an air molecule [m/s]
      REAL :: RHOA                            ! density of air [kg/m^3]
      REAL :: LAMBDAA                         ! mean free path of air [m]
      REAL :: KNAI, KNAJ                      ! Knudsen numbers in air for particles I, J [1]
      REAL :: NUA                             ! kinematic viscosity of air [m^2/s]
      REAL :: VFI, VFJ                        ! terminal fall speed [m/s]

      RI = 0.5E-06 * DI                                                       ! convert from [um] to [m]
      RJ = 0.5E-06 * DJ                                                       ! convert from [um] to [m]
      !-----------------------------------------------------------------------------------------------------------------
      ! Dynamic viscosity of air [kg/m/s], Jacobson, 2005, eq.(4.54).
      ! Mean thermal velocity of an air molecule [m/s].
      ! Density of air [kg/m^3] from the ideal gas law.
      ! Mean free path of an air molecule [m], Jacobson, 2005, eq.(15.24).
      !-----------------------------------------------------------------------------------------------------------------
      ETAA = 1.8325E-05 * (416.16/(TEMPK + 120.0)) * (TEMPK/296.16)**1.5      ! [kg/m/s]
      VABAR = SQRT( VACONST * TEMPK )                                         ! [m/s]
      RHOA = PRES * M1AIR / ( KB * TEMPK )                                    ! [kg/m^3]
      LAMBDAA = 2.0 * ETAA / ( RHOA * VABAR )                                 ! [m]
      KNAI = LAMBDAA / RI                                                     ! [1]
      KNAJ = LAMBDAA / RJ                                                     ! [1]
      GI = 1.0 + KNAI * ( APRIME + BPRIME * EXP(-CPRIME/KNAI) )               ! Cunningham slip-flow correction [1]
      GJ = 1.0 + KNAJ * ( APRIME + BPRIME * EXP(-CPRIME/KNAJ) )               ! Cunningham slip-flow correction [1]
      VFI = 2.0 * RI*RI*(DENSP_KGM3-RHOA)*GGRAV*GI/(9.0*ETAA)                 ! fall velocity [m/s] 
      VFJ = 2.0 * RJ*RJ*(DENSP_KGM3-RHOA)*GGRAV*GJ/(9.0*ETAA)                 ! fall velocity [m/s] 
      NUA = ETAA / RHOA                                                       ! kinematic viscosity of air [m^2/s]
      KTIIJ = (PI*ED**0.75/(GGRAV*NUA**0.25))*(RI+RJ)**2*ABS(VFI-VFJ)         ! turbulent inertial coag. coef. [m^3/s]
      KTSIJ = SQRT(8.0*PI*ED/(15.0*NUA))*(RI+RJ)**3                           ! turbulent shear    coag. coef. [m^3/s]
      ! WRITE(*,'(10D12.4)') KTIIJ, KTSIJ, RI, RJ, VFI, VFJ
      END SUBROUTINE TURB_COAG_COEF


      SUBROUTINE TOTAL_COAG_COEF( DI, DJ, TEMPK, PRES, KTOTIJ )
!-------------------------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the total coagulation coefficient. 
!-------------------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL, intent(in) :: DI, DJ                          ! particle diameters [um]
      REAL, intent(in) :: TEMPK                           ! ambient temperature [K]
      REAL, intent(in) :: PRES                            ! ambient pressure [Pa]
      REAL :: KBIJ                            ! Brownian                                  kernel [m^3/s/particle]
      REAL :: KDEIJ                           ! convective Brownian diffusion enhancement kernel [m^3/s/particle]
      REAL :: KGCIJ                           ! gravitational collection                  kernel [m^3/s/particle]
      REAL :: KTIIJ                           ! turbulent inertial                        kernel [m^3/s/particle]
      REAL :: KTSIJ                           ! turbulent shear                           kernel [m^3/s/particle]
      REAL, intent(out) :: KTOTIJ                          ! total                                     kernel [m^3/s/particle]
      CALL BROWNIAN_COAG_COEF( DI, DJ, TEMPK, PRES, KBIJ  )
      CALL CBDE_COAG_COEF    ( DI, DJ, TEMPK, PRES, KDEIJ )
      CALL GRAVCOLL_COAG_COEF( DI, DJ, TEMPK, PRES, KGCIJ )
      CALL TURB_COAG_COEF    ( DI, DJ, TEMPK, PRES, KTIIJ, KTSIJ )
      KTOTIJ = KBIJ + KDEIJ + KGCIJ + KTIIJ + KTSIJ
      END SUBROUTINE TOTAL_COAG_COEF

      END MODULE AERO_COAG


