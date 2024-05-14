MODULE AERO_COAG_CONFIG

    !----------------------------------------------------------
    ! TWO STEPS (x5l): 
    !1. construct Production/Loss 1/0 coefficient
    !             Production -GIKLQ(I,K,L,Q): K+L -> I, increase mass specie Q in pop I, K/L symmetric
    !             Production -DIKL(I,K,L)   : K+L -> I, increase particle number in pop I, K/L symmetric
    !             Loss       -DIJ(I,J)      : Loss of I from I+J, not symmetric
    !
    !2. calculate dyanmical coagulation coefficient
    !
    !---------------------------------------------------------
    use mpp_mod,               only : input_nml_file
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
    implicit none
    private
    public  SETUP_COAG_TENSORS
    character(len=3), allocatable,public, protected :: MODE_NAME(:), CITABLE(:,:)
    integer, public, protected :: NMASS_SPCS 
    integer, allocatable, public, protected :: NM(:), PROD_INDEX(:,:), GIKLQ(:,:,:,:), DIKL(:,:,:), DIJ(:,:), nDIKL
    character(len=4), allocatable :: NM_SPC_NAME(:,:)
    integer :: NMODES, NWEIGHTS

    type GIKLQ_type
        integer :: n
        integer, allocatable :: l(:), k(:)
        integer, allocatable :: qq(:)
     end type GIKLQ_type

    type DIKL_type
        integer :: i
        integer :: k
        integer :: l
    end type DIKL_type

    type (GIKLQ_type), allocatable, public, protected :: GIKLQ_control(:)
    type (DIKL_type), allocatable, public, protected :: DIKL_control(:)

    integer :: ierr, io, logunit, verbose, unit
    character(len=7), parameter :: module_name = 'matrix'
    contains

    SUBROUTINE SETUP_COAG_TENSORS(coag_configuration)
            !NMODES=npop-1=11, NMASS=5
        !!-------------------------------------------------------------------------------------------------------------------
        !     Step 1: Production/Loss 1/0 coefficient: Production-GIKLQ(I,K,L,Q), Production-DIKL(I,K,L), Loss-DIJ(I,J)
        !      Routine to define the g_ikl,q, the d_ikl, and the d_ij.
        !      All elements GIKLQ(I,K,L,Q), DIKL(I,K,L), and DIJ(I,J) were checked through printouts available below.
        !     NM_SPC_NAME(I,:) contains the names of the mass species defined for mode I.
        !-------------------------------------------------------------------------------------------------------------------
        IMPLICIT NONE
        character(len=32), intent(in) :: coag_configuration
        INTEGER :: I,J,K,L,Q,QQ, ntot, n
        integer :: ikl, ip, klq !debug output

        !-------------------------------------------------------------
        !configuration set-up based on matrix_nml
        !------------------------------------------------------
        if (lowercase(trim(coag_configuration)) .eq. "full") then
                NMODES = 12
                NWEIGHTS = NMODES
                allocate(MODE_NAME(NMODES)) !in total 12 modes
                !integer :: I_AKK = 1, I_ACC = 2, I_DD1 = 3, I_DD2 = 4 ! index of population in matrix
                !integer :: I_SSA = 5, I_SSC = 6, I_OC1 = 7, I_OC2 = 8 ! if exit, index >= 1, otherwise = -1
                !integer :: I_BC1 = 9, I_BC2 = 10, I_MXA = 11, I_MXC = 12, I_EXT = 13
                !MXA:MXX in accumulation mode, MXC: MXX in coarse mode
                MODE_NAME(1:NMODES) = (/'AKK','ACC','DD1','DD2','SSA','SSC','OC1','OC2','BC1','BC2','MXA','MXC'/)
                allocate(CITABLE(NMODES,NMODES))
                CITABLE(1:NMODES, 1)=(/'AKK','ACC','DD1','DD2','SSA','SSC','OC1','OC2','BC1','BC2','MXA','MXC'/) ! AKK
                CITABLE(1:NMODES, 2)=(/'ACC','ACC','DD1','DD2','SSA','SSC','OC1','OC2','BC1','BC2','MXA','MXC'/) ! ACC
                CITABLE(1:NMODES, 3)=(/'DD1','DD1','DD1','DD2','MXA','MXC','MXA','MXA','MXA','MXA','MXA','MXC'/) ! DD1
                CITABLE(1:NMODES, 4)=(/'DD2','DD2','DD2','DD2','MXC','MXC','MXC','MXC','MXC','MXC','MXC','MXC'/) ! DD2
                CITABLE(1:NMODES, 5)=(/'SSA','SSA','MXA','MXC','SSA','SSC','MXA','MXA','MXA','MXA','MXA','MXC'/) ! SSA
                CITABLE(1:NMODES, 6)=(/'SSC','SSC','MXC','MXC','SSC','SSC','MXC','MXC','MXC','MXC','MXC','MXC'/) ! SSC
                CITABLE(1:NMODES, 7)=(/'OC1','OC1','MXA','MXC','MXA','MXC','OC1','OC2','MXA','MXA','MXA','MXC'/) ! OC1
                CITABLE(1:NMODES, 8)=(/'OC2','OC2','MXA','MXC','MXA','MXC','OC2','OC2','MXA','MXA','MXA','MXC'/) ! OC2
                CITABLE(1:NMODES, 9)=(/'BC1','BC1','MXA','MXC','MXA','MXC','MXA','MXA','BC1','BC2','MXA','MXC'/) ! BC1
                CITABLE(1:NMODES,10)=(/'BC2','BC2','MXA','MXC','MXA','MXC','MXA','MXA','BC2','BC2','MXA','MXC'/) ! BC2
                CITABLE(1:NMODES,11)=(/'MXA','MXA','MXA','MXC','MXA','MXC','MXA','MXA','MXA','MXA','MXA','MXC'/) ! MXA
                CITABLE(1:NMODES,12)=(/'MXC','MXC','MXC','MXC','MXC','MXC','MXC','MXC','MXC','MXC','MXC','MXC'/) ! MXC
                allocate( NM(NMODES)) !number of mass species in each mode
                NMASS_SPCS = 5 !sulf, dust, seas, ocar, bcar 
                NM(1:NMODES)=(/1,1,2,2,2,2,2,2,2,2,5,5/)
                allocate(NM_SPC_NAME(NMODES,NMASS_SPCS))
                NM_SPC_NAME(1, 1:NMASS_SPCS)=(/'SULF','    ','    ','    ','    '/) !AKK
                NM_SPC_NAME(2, 1:NMASS_SPCS)=(/'SULF','    ','    ','    ','    '/) !ACC
                NM_SPC_NAME(3, 1:NMASS_SPCS)=(/'SULF','DUST','    ','    ','    '/) !DD1
                NM_SPC_NAME(4, 1:NMASS_SPCS)=(/'SULF','DUST','    ','    ','    '/) !DD2
                NM_SPC_NAME(5, 1:NMASS_SPCS)=(/'SULF','SEAS','    ','    ','    '/) !SSA
                NM_SPC_NAME(6, 1:NMASS_SPCS)=(/'SULF','SEAS','    ','    ','    '/) !SSC
                NM_SPC_NAME(7, 1:NMASS_SPCS)=(/'SULF','OCAR','    ','    ','    '/) !OC1
                NM_SPC_NAME(8, 1:NMASS_SPCS)=(/'SULF','OCAR','    ','    ','    '/) !OC2
                NM_SPC_NAME(9, 1:NMASS_SPCS)=(/'SULF','BCAR','    ','    ','    '/) !BC1
                NM_SPC_NAME(10,1:NMASS_SPCS)=(/'SULF','BCAR','    ','    ','    '/) !BC2
                NM_SPC_NAME(11,1:NMASS_SPCS)=(/'SULF','DUST','SEAS','OCAR','BCAR'/) !MXA
                NM_SPC_NAME(12,1:NMASS_SPCS)=(/'SULF','DUST','SEAS','OCAR','BCAR'/) !MXC
                allocate(PROD_INDEX(NMODES,NMASS_SPCS))
                PROD_INDEX(1, 1:NMASS_SPCS)=(/1,0,0,0,0/) !AKK
                PROD_INDEX(2, 1:NMASS_SPCS)=(/1,0,0,0,0/) !ACC
                PROD_INDEX(3, 1:NMASS_SPCS)=(/1,2,0,0,0/) !DD1
                PROD_INDEX(4, 1:NMASS_SPCS)=(/1,2,0,0,0/) !DD2
                PROD_INDEX(5, 1:NMASS_SPCS)=(/1,3,0,0,0/) !SSA
                PROD_INDEX(6, 1:NMASS_SPCS)=(/1,3,0,0,0/) !SSC
                PROD_INDEX(7, 1:NMASS_SPCS)=(/1,4,0,0,0/) !OC1
                PROD_INDEX(8, 1:NMASS_SPCS)=(/1,4,0,0,0/) !OC2
                PROD_INDEX(9, 1:NMASS_SPCS)=(/1,5,0,0,0/) !BC1
                PROD_INDEX(10,1:NMASS_SPCS)=(/1,5,0,0,0/) !BC2
                PROD_INDEX(11,1:NMASS_SPCS)=(/1,2,3,4,5/) !MXA
                PROD_INDEX(12,1:NMASS_SPCS)=(/1,2,3,4,5/) !MXC 
        else
               call ERROR_MESG('matrix coagulation configuration not properly defined','check nml for coag_configuration ', FATAL) 
        endif

!----------------------------------------------
!               debug tested
!----------------------------------------------
!     if (mpp_root_pe().eq.mpp_pe()) then
!        write(*,*) 'inside'
!        write(*,*) "MODE_NAME", MODE_NAME
!        write(*,*) "CITABLE"
!        write(*,*) size(CITABLE, 1)
!        do n=1,12
!                write(*,*) n,  CITABLE(n,1)
!                write(*,*) n,  CITABLE(n,:)
!        enddo
!        write(*,*) "NM"
!        write(*, *) NM
!        write(*,*) "PROD_INDEX"
!        do n=1,12
!                write(*,*) n,  PROD_INDEX(n, 1:NMASS_SPCS)
!        enddo
!        write(*,*) 'end_of_inside'
!    endif        
!


         allocate(GIKLQ(NMODES,NMODES,NMODES,NMASS_SPCS))
         allocate(DIKL(NMODES,NMODES, NMODES))
         allocate(DIJ(NMODES,NMODES))
         GIKLQ(:,:,:,:) = 0
         DIKL(:,:,:) = 0
         DIJ(:,:) = 0
         !-------------------------------------------------------------------------------------------------------------
         ! The tensors g_ikl,q and d_ikl are symmetric in K and L.
         !
         ! GIKLQ is unity if coagulation of modes K and L produce mass of species Q
         !       in mode I, and zero otherwise.
         !
         ! DIKL is unity if coagulation of modes K and L produce particles
         !      in mode I, and zero otherwise.
         !      Neither mode K nor mode L can be mode I for a nonzero DIKL:
         !      all three modes I, K, L must be different modes.
         !-------------------------------------------------------------------------------------------------------------
         DO I=1, NMODES
         DO K=1, NMODES
         DO L=K+1, NMODES                              ! Mode L is the same as mode K.
         IF ( CITABLE(K,L) .EQ. MODE_NAME(I) ) THEN  ! modes K and L produce mode I
                 IF ( I .NE. K  .AND. I .NE. L ) THEN      ! omit intramodal coagulation
                         DIKL(I,K,L) = 1
                         DIKL(I,L,K) = 1
                 ENDIF
                 DO Q=1, NM(I)                             ! loop over all mass species in mode I
                 DO QQ=1, NMASS_SPCS                     ! loop over all principal mass species
                 !-----------------------------------------------------------------------------------------------------
                 ! Compare the name of mass species Q in mode I with that of mass species QQ in mode K (or L).
                 ! The inner loop is over all principal mass species since all species must be checked for
                 !   mode K (or L) for a potential match with species Q in mode I.
                 !-----------------------------------------------------------------------------------------------------
                 IF( NM_SPC_NAME(K,QQ) .EQ. NM_SPC_NAME(I,Q) ) THEN   ! mode K contains Q
                         IF( I .NE. K ) GIKLQ(I,K,L,Q) = 1   ! I and K must be different modes
                         IF( I .NE. K ) GIKLQ(I,L,K,Q) = 1   ! I and K must be different modes
                 ENDIF
                 IF( NM_SPC_NAME(L,QQ) .EQ. NM_SPC_NAME(I,Q) ) THEN   ! mode L contains Q
                         IF( I .NE. L ) GIKLQ(I,K,L,Q) = 1   ! I and L must be different modes
                         IF( I .NE. L ) GIKLQ(I,L,K,Q) = 1   ! I and L must be different modes
                 ENDIF
                 ENDDO
                 ENDDO
         ENDIF
         ENDDO
         ENDDO
         ENDDO

         if (allocated(dikl_control)) then
                 deallocate(dikl_control)
         end if
         allocate(dikl_control(count(dikl /= 0))) !XL note: count(dikl /= 0) = 2*nDIKL, control only record half od DIKL

         call initializeDiklControl(dikl_control, DIKL)

         if (allocated(giklq_control)) then
                 do i = 1, nweights
                 deallocate(giklq_control(i)%k)
                 deallocate(giklq_control(i)%l)
                 deallocate(giklq_control(i)%qq)
                 end do
         else
                 allocate(GIKLQ_control(NWEIGHTS))
         end if

         call initializeGiklqControl(GIKLQ_control, GIKLQ)  
         !-------------------------------------------------------------------------------------------------------------
         ! The tensor d_ij is not symmetric in I,J.
         !
         ! DIJ(I,J) is unity if coagulation of mode I with mode J results!
         !   in the removal of particles from mode I, and zero otherwise.
         !-------------------------------------------------------------------------------------------------------------
         DO I=1, NMODES
         DO J=1, NMODES
         DO K=1, NMODES                               ! Find the product mode of the I-J coagulation.
         IF( I .EQ. J ) CYCLE                       ! Omit intramodal interactions: --> I .NE. J .
         IF( CITABLE(I,J) .EQ. MODE_NAME(K) ) THEN  ! I-particles and J-particles are lost; K-particles are formed.
                 IF( I .NE. K ) DIJ(I,J) = 1              ! The K-particles are not I-particles (but may be J-particles),
         ENDIF                                      !   so I-particles are lost by this I-J interaction.
         ENDDO
         ENDDO
         ENDDO

        ! !----------------------------------------------
        ! !               debug tested
        ! !----------------------------------------------
        ! if (mpp_root_pe().eq.mpp_pe()) then
        !         write(*,*) 'inside'
        !         !        write(*,*) 'nDIKL', nDIKL
        !         !        do ikl = 1, nDIKL !ikl is 1-12, exclude ext already
        !         !        write(*,*) 'i,k,l', dikl_control(ikl)%i, dikl_control(ikl)%k, dikl_control(ikl)%l
        !         !        enddo
        !         write(*,*) 'giklq'
        !         DO I=1, NMODES
        !         DO K=1, NMODES
        !         DO L=K+1, NMODES                              ! Mode L is the same as mode K.
        !         DO Q=1, NMASS_SPCS                             ! loop over all mass species in mode I
        !         !-----------------------------------------------------------------------------------------------------
        !         ! Compare the name of mass species Q in mode I with that of mass species QQ in mode K (or L).
        !         ! The inner loop is over all principal mass species since all species must be checked for
        !         !   mode K (or L) for a potential match with species Q in mode I.
        !         !-----------------------------------------------------------------------------------------------------
        !         IF (GIKLQ(I,K,L,Q) > 0) THEN
        !                 write(*,*) 'I,K,L,Q,QQ, *', I,K,L,Q,QQ
        !                 write(*,*) 'GIKLQ(I,K,L,Q)', GIKLQ(I,K,L,Q)
        !         ENDIF
        !         IF (GIKLQ(I,L,K,Q) > 0) THEN
        !                 write(*,*) 'I,L,K,Q,QQ, *', I,L,K,Q, QQ
        !                 write(*,*) 'GIKLQ(I,L,K,Q)', GIKLQ(I,L,K,Q)
        !         ENDIF
        !         ENDDO
        !         ENDDO
        !         ENDDO
        !         ENDDO
        !         do ip = 1, NMODES
        !         write(*,*) 'ipop', ip
        !         do klq = 1, giklq_control(ip)%n
        !         write(*,*) 'klq', klq
        !         write(*,*) 'k,l,qq',  giklq_control(ip)%k(klq), giklq_control(ip)%l(klq), giklq_control(ip)%qq(klq)
        !         enddo
        !         enddo        
        !         write(*,*) 'DIJ'
        !         do I=1, NMODES
        !         DO J=1, NMODES
        !         if (DIJ(I,J)>0) then
        !                 write(*,*) 'i,j', I, J
        !         endif
        !         enddo
        !         enddo
        ! endif



 END SUBROUTINE




    !--------------------------------------------------------------------
    !record the pair with non-zero D_ikl
    !type DIKL_type
    !        integer :: i
    !        integer :: k
    !        integer :: l
    !end type DIKL_type
    !---------------------------------------------------------------------
    subroutine initializeDiklControl(control, mask)
            type (dikl_type) :: control(:)
            integer, intent(in) :: mask(:,:,:)
            integer :: i, k, l, n
            n = 0
            do k = 1, NWEIGHTS
            do l = k+1, NWEIGHTS
            do i = 1, NWEIGHTS
            if (mask(i,k,l) /= 0) then
                    n = n + 1
                    control(n)%i = i
                    control(n)%k = k
                    control(n)%l = l
            end if
            end do
            end do
            end do
            NDIKL = n
    end subroutine initializeDiklControl
    !--------------------------------------------------------------------
    !record the pair with non-zero Giklq
    !   type GIKLQ_type
    !    integer :: n
    !    integer, allocatable :: l(:), k(:)
    !    integer, allocatable :: qq(:)
    !  end type GIKLQ_type
    !---------------------------------------------------------------------
    subroutine initializeGiklqControl(control, mask)
            type (GIKLQ_type) :: control(:)
            integer, intent(in) :: mask(:,:,:,:)
            integer :: i, q, k, l, n, qq, ip,  nTotal

        do i = 1, NWEIGHTS
            ! 1) count contributing cases for mode i
            n = 0
            do q = 1, NM(i) !mass species in each mode
            do k = 1, nmodes
            do l = k+1, nmodes
            if (mask(i,k,l,q) /= 0) then
                    if (i /= l) then
                            n = n + 1
                    end if
                    if (i /= k) then
                            n = n + 1
                    end if
            end if
            end do
            end do
            end do
            nTotal = n
            ! 2) allocate nTotal entries
            control(i)%n = nTotal
            allocate(control(i)%k(nTotal))
            allocate(control(i)%l(nTotal))
            allocate(control(i)%qq(nTotal))
            ! 3) repeat sweep, but now assign k,l,qq
            n = 0
            do q = 1, nm(i)
            do k = 1, nmodes
            do l = k+1, nmodes
            if (mask(i,k,l,q) /= 0) then
                    if (i /= l) then
                            n = n + 1
                            control(i)%k(n) = k
                            control(i)%l(n) = l
                            qq = prod_index(i,q)
                            control(i)%qq(n) = qq
                    end if
                    if (i /= k) then
                            n = n + 1
                            control(i)%k(n) = l
                            control(i)%l(n) = k
                            qq = prod_index(i,q)
                            control(i)%qq(n) = qq
                    end if
            end if
            end do
            end do
            end do
            control(i)%n = n

            !!GIKL_control printout
            !if (mpp_root_pe().eq.mpp_pe()) then
            !     write(*,*) 'gikl_control'
            !     write(*,*) 'pop_num', i
            !     write(*,*) 'control(i)%n', control(i)%n
            !     if (i > 10) then
            !             do ip=1, control(i)%n 
            !                    write(*,*) 'k,l,qq', control(i)%k(ip),control(i)%l(ip), control(i)%qq(ip)
            !            enddo
            !    endif

            !endif

        end do

    end subroutine initializeGiklqControl




END MODULE




!
!      !-------------------------------------------------------------------------------------------------------------
!      ! The tensors g_ikl,q and d_ikl are symmetric in K and L.
!      !
!      ! GIKLQ is unity if coagulation of modes K and L produce mass of species Q
!      !       in mode I, and zero otherwise.
!      !
!      ! DIKL is unity if coagulation of modes K and L produce particles
!      !      in mode I, and zero otherwise.
!      !      Neither mode K nor mode L can be mode I for a nonzero DIKL:
!      !      all three modes I, K, L must be different modes.
!      !-------------------------------------------------------------------------------------------------------------
!      DO I=1, NMODES
!      DO K=1, NMODES
!      DO L=K+1, NMODES                              ! Mode L is the same as mode K.
!        IF ( CITABLE(K,L) .EQ. MODE_NAME(I) ) THEN  ! modes K and L produce mode I
!          ! WRITE(36,*)'MODE_NAME(I) = ', MODE_NAME(I)
!          IF ( I .NE. K  .AND. I .NE. L ) THEN      ! omit intramodal coagulation
!            DIKL(I,K,L) = 1
!            DIKL(I,L,K) = 1
!          ENDIF
!          DO Q=1, NM(I)                             ! loop over all mass species in mode I
!            DO QQ=1, NMASS_SPCS                     ! loop over all principal mass species
!              !-----------------------------------------------------------------------------------------------------
!              ! Compare the name of mass species Q in mode I with that of mass species QQ in mode K (or L).
!              ! The inner loop is over all principal mass species since all species must be checked for
!              !   mode K (or L) for a potential match with species Q in mode I.
!              !-----------------------------------------------------------------------------------------------------
!              IF( NM_SPC_NAME(K,QQ) .EQ. NM_SPC_NAME(I,Q) ) THEN   ! mode K contains Q
!                IF( I .NE. K ) GIKLQ(I,K,L,Q) = 1   ! I and K must be different modes
!                IF( I .NE. K ) GIKLQ(I,L,K,Q) = 1   ! I and K must be different modes
!              ENDIF
!              IF( NM_SPC_NAME(L,QQ) .EQ. NM_SPC_NAME(I,Q) ) THEN   ! mode L contains Q
!                IF( I .NE. L ) GIKLQ(I,K,L,Q) = 1   ! I and L must be different modes
!                IF( I .NE. L ) GIKLQ(I,L,K,Q) = 1   ! I and L must be different modes
!              ENDIF
!            ENDDO
!          ENDDO
!        ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!
!      if (allocated(dikl_control)) then
!        deallocate(dikl_control)
!      end if
!      allocate(dikl_control(count(dikl /= 0)))
!
!      call initializeDiklControl(dikl_control, DIKL)
!
!      if (allocated(giklq_control)) then
!        do i = 1, nweights
!          deallocate(giklq_control(i)%k)
!          deallocate(giklq_control(i)%l)
!          deallocate(giklq_control(i)%qq)
!        end do
!      else
!        allocate(GIKLQ_control(NWEIGHTS))
!      end if
!
!      call initializeGiklqControl(GIKLQ_control, GIKLQ)
!
!      !-------------------------------------------------------------------------------------------------------------
!      ! The tensor d_ij is not symmetric in I,J.
!      !
!      ! DIJ(I,J) is unity if coagulation of mode I with mode J results
!      !   in the removal of particles from mode I, and zero otherwise.
!      !-------------------------------------------------------------------------------------------------------------
!      DO I=1, NMODES
!      DO J=1, NMODES
!        DO K=1, NMODES                               ! Find the product mode of the I-J coagulation.
!          IF( I .EQ. J ) CYCLE                       ! Omit intramodal interactions: --> I .NE. J .
!          IF( CITABLE(I,J) .EQ. MODE_NAME(K) ) THEN  ! I-particles and J-particles are lost; K-particles are formed.
!            IF( I .NE. K ) DIJ(I,J) = 1              ! The K-particles are not I-particles (but may be J-particles),
!          ENDIF                                      !   so I-particles are lost by this I-J interaction.
!        ENDDO
!      ENDDO
!      ENDDO
!      xDIJ = DIJ
!
!      IF( .NOT. WRITE_TENSORS ) RETURN
!
!      !-------------------------------------------------------------------------
!      ! Write the g_ikl,q.
!      !-------------------------------------------------------------------------
!      DO I=1, NMODES
!        WRITE(AUNIT1,'(/2A)') 'g_iklq for MODE ', MODE_NAME(I)
!        DO Q=1, NM(I)
!          WRITE(AUNIT1,'(/A,I3,3X,3A5/)') 'Q, NM_SPC_NAME(I,Q), MODE',
!     &                        Q, NM_SPC_NAME(I,Q), '-->', MODE_NAME(I)
!          IF ( SUM( GIKLQ(I,1:NMODES,1:NMODES,Q) ) .EQ. 0 ) CYCLE
!          WRITE(AUNIT1,'(5X,16A5)') MODE_NAME(1:NMODES)
!          DO K=1, NMODES
!            WRITE(AUNIT1,'(A5,16I5)') MODE_NAME(K),GIKLQ(I,K,1:NMODES,Q)
!          ENDDO
!        ENDDO
!      ENDDO
!
!      !-------------------------------------------------------------------------
!      ! Write the d_ikl.
!      !-------------------------------------------------------------------------
!      DO I=1, NMODES
!        WRITE(AUNIT1,'(/2A)') 'd_ikl for MODE ', MODE_NAME(I)
!        IF ( SUM( DIKL(I,1:NMODES,1:NMODES) ) .EQ. 0 ) CYCLE
!        WRITE(AUNIT1,'(5X,16A5)') MODE_NAME(1:NMODES)
!        DO K=1, NMODES
!          WRITE(AUNIT1,'(A5,16I5)') MODE_NAME(K),DIKL(I,K,1:NMODES)
!        ENDDO
!      ENDDO
!
!      !-------------------------------------------------------------------------
!      ! Write the d_ij.
!      !-------------------------------------------------------------------------
!      WRITE(AUNIT1,'(/2A)') 'd_ij'
!      WRITE(AUNIT1,'(5X,16A5)') MODE_NAME(1:NMODES)
!      DO I=1, NMODES
!        WRITE(AUNIT1,'(A5,16I5)') MODE_NAME(I),DIJ(I,1:NMODES)
!      ENDDO
!
!90000 FORMAT(14A4)
!      RETURN
!
!      contains
!
!      subroutine initializeDiklControl(control, mask)
!      type (dikl_type) :: control(:)
!      integer, intent(in) :: mask(:,:,:)
!      integer :: i, k, l, n
!
!      n = 0
!      do k = 1, NWEIGHTS
!        do l = k+1, NWEIGHTS
!          do i = 1, NWEIGHTS
!            if (mask(i,k,l) /= 0) then
!              n = n + 1
!              control(n)%i = i
!              control(n)%k = k
!              control(n)%l = l
!            end if
!          end do
!        end do
!      end do
!      NDIKL = n
!      end subroutine initializeDiklControl
!
!      subroutine initializeGiklqControl(control, mask)
!      type (GIKLQ_type) :: control(:)
!      integer, intent(in) :: mask(:,:,:,:)
!
!      integer :: i, q, k, l, n, nTotal
!
!      do i = 1, NWEIGHTS
!        ! 1) count contributing cases for mode i
!        n = 0
!        do q = 1, nm(i)
!          do k = 1, nmodes
!            do l = k+1, nmodes
!              if (mask(i,k,l,q) /= 0) then
!                if (i /= l) then
!                  n = n + 1
!                end if
!                if (i /= k) then
!                  n = n + 1
!                end if
!              end if
!            end do
!          end do
!        end do
!        nTotal = n
!        ! 2) allocate nTotal entries
!        control(i)%n = nTotal
!        allocate(control(i)%k(nTotal))
!        allocate(control(i)%l(nTotal))
!        allocate(control(i)%qq(nTotal))
!
!        ! 3) repeat sweep, but now assign k,l,qq
!        n = 0
!        do q = 1, nm(i)
!          do k = 1, nmodes
!            do l = k+1, nmodes
!              if (mask(i,k,l,q) /= 0) then
!                if (i /= l) then
!                  n = n + 1
!                  control(i)%k(n) = k
!                  control(i)%l(n) = l
!                  qq = prod_index(i,q)
!                  control(i)%qq(n) = qq
!                end if
!                if (i /= k) then
!                  n = n + 1
!                  control(i)%k(n) = l
!                  control(i)%l(n) = k
!                  qq = prod_index(i,q)
!                  control(i)%qq(n) = qq
!                end if
!              end if
!            end do
!          end do
!        end do
!        control(i)%n = n
!      end do
!
!      end subroutine initializeGiklqControl
!
!
!      END SUBROUTINE SETUP_COAG_TENSORS
!
