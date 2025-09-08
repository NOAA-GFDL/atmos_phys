module aero_coag_config

    !----------------------------------------------------------
    ! two steps (x5l): 
    !1. construct production/loss 1/0 coefficient
    !             production -giklq(i,k,l,q): k+l -> i, increase mass specie q in pop i, k/l symmetric
    !             production -dikl(i,k,l)   : k+l -> i, increase particle number in pop i, k/l symmetric
    !             loss       -dij(i,j)      : loss of i from i+j, not symmetric
    !
    !2. calculate dyanmical coagulation coefficient
    !
    !---------------------------------------------------------
    use mpp_mod,               only : input_nml_file
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
    implicit none
    private
    public  setup_coag_tensors
    character(len=3), allocatable,public, protected :: mode_name(:), citable(:,:)
    integer, public, protected :: nmass_spcs 
    integer, allocatable, public, protected :: nm(:), prod_index(:,:), giklq(:,:,:,:), dikl(:,:,:), dij(:,:), ndikl
    character(len=4), allocatable :: nm_spc_name(:,:)
    integer :: nmodes, nweights

    type giklq_type
        integer :: n
        integer, allocatable :: l(:), k(:)
        integer, allocatable :: qq(:)
     end type giklq_type

    type dikl_type
        integer :: i
        integer :: k
        integer :: l
    end type dikl_type

    type (giklq_type), allocatable, public, protected :: giklq_control(:)
    type (dikl_type), allocatable, public, protected :: dikl_control(:)

    integer :: ierr, io, logunit, verbose, unit
    character(len=7), parameter :: module_name = 'matrix'
    contains

    subroutine setup_coag_tensors(coag_configuration)
            !nmodes=npop-1=11, nmass=5
        !!-------------------------------------------------------------------------------------------------------------------
        !     step 1: production/loss 1/0 coefficient: production-giklq(i,k,l,q), production-dikl(i,k,l), loss-dij(i,j)
        !      routine to define the g_ikl,q, the d_ikl, and the d_ij.
        !      all elements giklq(i,k,l,q), dikl(i,k,l), and dij(i,j) were checked through printouts available below.
        !     nm_spc_name(i,:) contains the names of the mass species defined for mode i.
        !-------------------------------------------------------------------------------------------------------------------
        implicit none
        character(len=32), intent(in) :: coag_configuration
        integer :: i,j,k,l,q,qq, ntot, n
        integer :: ikl, ip, klq !debug output

        !-------------------------------------------------------------
        !configuration set-up based on matrix_nml
        !------------------------------------------------------
        if (lowercase(trim(coag_configuration)) .eq. "full") then
                nmodes = 12
                nweights = nmodes
                allocate(mode_name(nmodes)) !in total 12 modes
                !integer :: i_akk = 1, i_acc = 2, i_dd1 = 3, i_dd2 = 4 ! index of population in matrix
                !integer :: i_ssa = 5, i_ssc = 6, i_oc1 = 7, i_oc2 = 8 ! if exit, index >= 1, otherwise = -1
                !integer :: i_bc1 = 9, i_bc2 = 10, i_mxa = 11, i_mxc = 12, i_ext = 13
                !mxa:mxx in accumulation mode, mxc: mxx in coarse mode
                mode_name(1:nmodes) = (/'akk','acc','dd1','dd2','ssa','ssc','oc1','oc2','bc1','bc2','mxa','mxc'/)
                allocate(citable(nmodes,nmodes))
                citable(1:nmodes, 1)=(/'akk','acc','dd1','dd2','ssa','ssc','oc1','oc2','bc1','bc2','mxa','mxc'/) ! akk
                citable(1:nmodes, 2)=(/'acc','acc','dd1','dd2','ssa','ssc','oc1','oc2','bc1','bc2','mxa','mxc'/) ! acc
                citable(1:nmodes, 3)=(/'dd1','dd1','dd1','dd2','mxa','mxc','mxa','mxa','mxa','mxa','mxa','mxc'/) ! dd1
                citable(1:nmodes, 4)=(/'dd2','dd2','dd2','dd2','mxc','mxc','mxc','mxc','mxc','mxc','mxc','mxc'/) ! dd2
                citable(1:nmodes, 5)=(/'ssa','ssa','mxa','mxc','ssa','ssc','mxa','mxa','mxa','mxa','mxa','mxc'/) ! ssa
                citable(1:nmodes, 6)=(/'ssc','ssc','mxc','mxc','ssc','ssc','mxc','mxc','mxc','mxc','mxc','mxc'/) ! ssc
                citable(1:nmodes, 7)=(/'oc1','oc1','mxa','mxc','mxa','mxc','oc1','oc1','mxa','mxa','mxa','mxc'/) ! oc1
                citable(1:nmodes, 8)=(/'oc2','oc2','mxa','mxc','mxa','mxc','oc1','oc2','mxa','mxa','mxa','mxc'/) ! oc2
                citable(1:nmodes, 9)=(/'bc1','bc1','mxa','mxc','mxa','mxc','mxa','mxa','bc1','bc1','mxa','mxc'/) ! bc1
                citable(1:nmodes,10)=(/'bc2','bc2','mxa','mxc','mxa','mxc','mxa','mxa','bc1','bc2','mxa','mxc'/) ! bc2
                citable(1:nmodes,11)=(/'mxa','mxa','mxa','mxc','mxa','mxc','mxa','mxa','mxa','mxa','mxa','mxc'/) ! mxa
                citable(1:nmodes,12)=(/'mxc','mxc','mxc','mxc','mxc','mxc','mxc','mxc','mxc','mxc','mxc','mxc'/) ! mxc
                allocate( nm(nmodes)) !number of mass species in each mode
                nmass_spcs = 5 !order: sulf, bcar, ocar, dust, seas 
                nm(1:nmodes)=(/1,1,2,2,2,2,2,2,2,2,5,5/)
                allocate(nm_spc_name(nmodes,nmass_spcs))
                nm_spc_name(1, 1:nmass_spcs)=(/'sulf','    ','    ','    ','    '/) !akk
                nm_spc_name(2, 1:nmass_spcs)=(/'sulf','    ','    ','    ','    '/) !acc
                nm_spc_name(3, 1:nmass_spcs)=(/'sulf','dust','    ','    ','    '/) !dd1
                nm_spc_name(4, 1:nmass_spcs)=(/'sulf','dust','    ','    ','    '/) !dd2
                nm_spc_name(5, 1:nmass_spcs)=(/'sulf','seas','    ','    ','    '/) !ssa
                nm_spc_name(6, 1:nmass_spcs)=(/'sulf','seas','    ','    ','    '/) !ssc
                nm_spc_name(7, 1:nmass_spcs)=(/'sulf','ocar','    ','    ','    '/) !oc1
                nm_spc_name(8, 1:nmass_spcs)=(/'sulf','ocar','    ','    ','    '/) !oc2
                nm_spc_name(9, 1:nmass_spcs)=(/'sulf','bcar','    ','    ','    '/) !bc1
                nm_spc_name(10,1:nmass_spcs)=(/'sulf','bcar','    ','    ','    '/) !bc2
                nm_spc_name(11,1:nmass_spcs)=(/'sulf','dust','seas','ocar','bcar'/) !mxa
                nm_spc_name(12,1:nmass_spcs)=(/'sulf','dust','seas','ocar','bcar'/) !mxc
                allocate(prod_index(nmodes,nmass_spcs))
                prod_index(1, 1:nmass_spcs)=(/1,0,0,0,0/) !akk
                prod_index(2, 1:nmass_spcs)=(/1,0,0,0,0/) !acc
                prod_index(3, 1:nmass_spcs)=(/1,4,0,0,0/) !dd1
                prod_index(4, 1:nmass_spcs)=(/1,4,0,0,0/) !dd2
                prod_index(5, 1:nmass_spcs)=(/1,5,0,0,0/) !ssa
                prod_index(6, 1:nmass_spcs)=(/1,5,0,0,0/) !ssc
                prod_index(7, 1:nmass_spcs)=(/1,3,0,0,0/) !oc1
                prod_index(8, 1:nmass_spcs)=(/1,3,0,0,0/) !oc2
                prod_index(9, 1:nmass_spcs)=(/1,2,0,0,0/) !bc1
                prod_index(10,1:nmass_spcs)=(/1,2,0,0,0/) !bc2
                prod_index(11,1:nmass_spcs)=(/1,4,5,3,2/) !mxa
                prod_index(12,1:nmass_spcs)=(/1,4,5,3,2/) !mxc 
        else
               call error_mesg('matrix coagulation configuration not properly defined','check nml for coag_configuration ', fatal) 
        endif

!----------------------------------------------
!               debug tested
!----------------------------------------------
!     if (mpp_root_pe().eq.mpp_pe()) then
!        write(*,*) 'inside'
!        write(*,*) "mode_name", mode_name
!        write(*,*) "citable"
!        write(*,*) size(citable, 1)
!        do n=1,12
!                write(*,*) n,  citable(n,1)
!                write(*,*) n,  citable(n,:)
!        enddo
!        write(*,*) "nm"
!        write(*, *) nm
!        write(*,*) "prod_index"
!        do n=1,12
!                write(*,*) n,  prod_index(n, 1:nmass_spcs)
!        enddo
!        write(*,*) 'end_of_inside'
!    endif        
!


         allocate(giklq(nmodes,nmodes,nmodes,nmass_spcs))
         allocate(dikl(nmodes,nmodes, nmodes))
         allocate(dij(nmodes,nmodes))
         giklq(:,:,:,:) = 0
         dikl(:,:,:) = 0
         dij(:,:) = 0
         !-------------------------------------------------------------------------------------------------------------
         ! the tensors g_ikl,q and d_ikl are symmetric in k and l.
         !
         ! giklq is unity if coagulation of modes k and l produce mass of species q
         !       in mode i, and zero otherwise.
         !
         ! dikl is unity if coagulation of modes k and l produce particles
         !      in mode i, and zero otherwise.
         !      neither mode k nor mode l can be mode i for a nonzero dikl:
         !      all three modes i, k, l must be different modes.
         !-------------------------------------------------------------------------------------------------------------
         do i=1, nmodes
         do k=1, nmodes
         do l=k+1, nmodes                              ! mode l is the same as mode k.
         if ( citable(k,l) .eq. mode_name(i) ) then  ! modes k and l produce mode i
                 if ( i .ne. k  .and. i .ne. l ) then      ! omit intramodal coagulation
                         dikl(i,k,l) = 1
                         dikl(i,l,k) = 1
                 endif
                 do q=1, nm(i)                             ! loop over all mass species in mode i
                 do qq=1, nmass_spcs                     ! loop over all principal mass species
                 !-----------------------------------------------------------------------------------------------------
                 ! compare the name of mass species q in mode i with that of mass species qq in mode k (or l).
                 ! the inner loop is over all principal mass species since all species must be checked for
                 !   mode k (or l) for a potential match with species q in mode i.
                 !-----------------------------------------------------------------------------------------------------
                 if( nm_spc_name(k,qq) .eq. nm_spc_name(i,q) ) then   ! mode k contains q
                         if( i .ne. k ) giklq(i,k,l,q) = 1   ! i and k must be different modes
                         if( i .ne. k ) giklq(i,l,k,q) = 1   ! i and k must be different modes
                 endif
                 if( nm_spc_name(l,qq) .eq. nm_spc_name(i,q) ) then   ! mode l contains q
                         if( i .ne. l ) giklq(i,k,l,q) = 1   ! i and l must be different modes
                         if( i .ne. l ) giklq(i,l,k,q) = 1   ! i and l must be different modes
                 endif
                 enddo
                 enddo
         endif
         enddo
         enddo
         enddo

         if (allocated(dikl_control)) then
                 deallocate(dikl_control)
         end if
         allocate(dikl_control(count(dikl /= 0))) !xl note: count(dikl /= 0) = 2*ndikl, control only record half od dikl

         call initializediklcontrol(dikl_control, dikl)

         if (allocated(giklq_control)) then
                 do i = 1, nweights
                 deallocate(giklq_control(i)%k)
                 deallocate(giklq_control(i)%l)
                 deallocate(giklq_control(i)%qq)
                 end do
         else
                 allocate(giklq_control(nweights))
         end if

         call initializegiklqcontrol(giklq_control, giklq)  
         !-------------------------------------------------------------------------------------------------------------
         ! the tensor d_ij is not symmetric in i,j.
         !
         ! dij(i,j) is unity if coagulation of mode i with mode j results!
         !   in the removal of particles from mode i, and zero otherwise.
         !-------------------------------------------------------------------------------------------------------------
         do i=1, nmodes
         do j=1, nmodes
         do k=1, nmodes                               ! find the product mode of the i-j coagulation.
         if( i .eq. j ) cycle                       ! omit intramodal interactions: --> i .ne. j .
         if( citable(i,j) .eq. mode_name(k) ) then  ! i-particles and j-particles are lost; k-particles are formed.
                 if( i .ne. k ) dij(i,j) = 1              ! the k-particles are not i-particles (but may be j-particles),
         endif                                      !   so i-particles are lost by this i-j interaction.
         enddo
         enddo
         enddo

        ! !----------------------------------------------
        ! !               debug tested
        ! !----------------------------------------------
        ! if (mpp_root_pe().eq.mpp_pe()) then
        !         write(*,*) 'inside'
        !         !        write(*,*) 'ndikl', ndikl
        !         !        do ikl = 1, ndikl !ikl is 1-12, exclude ext already
        !         !        write(*,*) 'i,k,l', dikl_control(ikl)%i, dikl_control(ikl)%k, dikl_control(ikl)%l
        !         !        enddo
        !         write(*,*) 'giklq'
        !         do i=1, nmodes
        !         do k=1, nmodes
        !         do l=k+1, nmodes                              ! mode l is the same as mode k.
        !         do q=1, nmass_spcs                             ! loop over all mass species in mode i
        !         !-----------------------------------------------------------------------------------------------------
        !         ! compare the name of mass species q in mode i with that of mass species qq in mode k (or l).
        !         ! the inner loop is over all principal mass species since all species must be checked for
        !         !   mode k (or l) for a potential match with species q in mode i.
        !         !-----------------------------------------------------------------------------------------------------
        !         if (giklq(i,k,l,q) > 0) then
        !                 write(*,*) 'i,k,l,q,qq, *', i,k,l,q,qq
        !                 write(*,*) 'giklq(i,k,l,q)', giklq(i,k,l,q)
        !         endif
        !         if (giklq(i,l,k,q) > 0) then
        !                 write(*,*) 'i,l,k,q,qq, *', i,l,k,q, qq
        !                 write(*,*) 'giklq(i,l,k,q)', giklq(i,l,k,q)
        !         endif
        !         enddo
        !         enddo
        !         enddo
        !         enddo
        !         do ip = 1, nmodes
        !         write(*,*) 'ipop', ip
        !         do klq = 1, giklq_control(ip)%n
        !         write(*,*) 'klq', klq
        !         write(*,*) 'k,l,qq',  giklq_control(ip)%k(klq), giklq_control(ip)%l(klq), giklq_control(ip)%qq(klq)
        !         enddo
        !         enddo        
        !         write(*,*) 'dij'
        !         do i=1, nmodes
        !         do j=1, nmodes
        !         if (dij(i,j)>0) then
        !                 write(*,*) 'i,j', i, j
        !         endif
        !         enddo
        !         enddo
        ! endif



 end subroutine




    !--------------------------------------------------------------------
    !record the pair with non-zero d_ikl
    !type dikl_type
    !        integer :: i
    !        integer :: k
    !        integer :: l
    !end type dikl_type
    !---------------------------------------------------------------------
    subroutine initializediklcontrol(control, mask)
            type (dikl_type) :: control(:)
            integer, intent(in) :: mask(:,:,:)
            integer :: i, k, l, n
            n = 0
            do k = 1, nweights
            do l = k+1, nweights
            do i = 1, nweights
            if (mask(i,k,l) /= 0) then
                    n = n + 1
                    control(n)%i = i
                    control(n)%k = k
                    control(n)%l = l
            end if
            end do
            end do
            end do
            ndikl = n
    end subroutine initializediklcontrol
    !--------------------------------------------------------------------
    !record the pair with non-zero giklq
    !   type giklq_type
    !    integer :: n
    !    integer, allocatable :: l(:), k(:)
    !    integer, allocatable :: qq(:)
    !  end type giklq_type
    !---------------------------------------------------------------------
    subroutine initializegiklqcontrol(control, mask)
            type (giklq_type) :: control(:)
            integer, intent(in) :: mask(:,:,:,:)
            integer :: i, q, k, l, n, qq, ip,  ntotal

        do i = 1, nweights
            ! 1) count contributing cases for mode i
            n = 0
            do q = 1, nm(i) !mass species in each mode
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
            ntotal = n
            ! 2) allocate ntotal entries
            control(i)%n = ntotal
            allocate(control(i)%k(ntotal))
            allocate(control(i)%l(ntotal))
            allocate(control(i)%qq(ntotal))
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

            !!gikl_control printout
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

    end subroutine initializegiklqcontrol




end module




!
!      !-------------------------------------------------------------------------------------------------------------
!      ! the tensors g_ikl,q and d_ikl are symmetric in k and l.
!      !
!      ! giklq is unity if coagulation of modes k and l produce mass of species q
!      !       in mode i, and zero otherwise.
!      !
!      ! dikl is unity if coagulation of modes k and l produce particles
!      !      in mode i, and zero otherwise.
!      !      neither mode k nor mode l can be mode i for a nonzero dikl:
!      !      all three modes i, k, l must be different modes.
!      !-------------------------------------------------------------------------------------------------------------
!      do i=1, nmodes
!      do k=1, nmodes
!      do l=k+1, nmodes                              ! mode l is the same as mode k.
!        if ( citable(k,l) .eq. mode_name(i) ) then  ! modes k and l produce mode i
!          ! write(36,*)'mode_name(i) = ', mode_name(i)
!          if ( i .ne. k  .and. i .ne. l ) then      ! omit intramodal coagulation
!            dikl(i,k,l) = 1
!            dikl(i,l,k) = 1
!          endif
!          do q=1, nm(i)                             ! loop over all mass species in mode i
!            do qq=1, nmass_spcs                     ! loop over all principal mass species
!              !-----------------------------------------------------------------------------------------------------
!              ! compare the name of mass species q in mode i with that of mass species qq in mode k (or l).
!              ! the inner loop is over all principal mass species since all species must be checked for
!              !   mode k (or l) for a potential match with species q in mode i.
!              !-----------------------------------------------------------------------------------------------------
!              if( nm_spc_name(k,qq) .eq. nm_spc_name(i,q) ) then   ! mode k contains q
!                if( i .ne. k ) giklq(i,k,l,q) = 1   ! i and k must be different modes
!                if( i .ne. k ) giklq(i,l,k,q) = 1   ! i and k must be different modes
!              endif
!              if( nm_spc_name(l,qq) .eq. nm_spc_name(i,q) ) then   ! mode l contains q
!                if( i .ne. l ) giklq(i,k,l,q) = 1   ! i and l must be different modes
!                if( i .ne. l ) giklq(i,l,k,q) = 1   ! i and l must be different modes
!              endif
!            enddo
!          enddo
!        endif
!      enddo
!      enddo
!      enddo
!
!      if (allocated(dikl_control)) then
!        deallocate(dikl_control)
!      end if
!      allocate(dikl_control(count(dikl /= 0)))
!
!      call initializediklcontrol(dikl_control, dikl)
!
!      if (allocated(giklq_control)) then
!        do i = 1, nweights
!          deallocate(giklq_control(i)%k)
!          deallocate(giklq_control(i)%l)
!          deallocate(giklq_control(i)%qq)
!        end do
!      else
!        allocate(giklq_control(nweights))
!      end if
!
!      call initializegiklqcontrol(giklq_control, giklq)
!
!      !-------------------------------------------------------------------------------------------------------------
!      ! the tensor d_ij is not symmetric in i,j.
!      !
!      ! dij(i,j) is unity if coagulation of mode i with mode j results
!      !   in the removal of particles from mode i, and zero otherwise.
!      !-------------------------------------------------------------------------------------------------------------
!      do i=1, nmodes
!      do j=1, nmodes
!        do k=1, nmodes                               ! find the product mode of the i-j coagulation.
!          if( i .eq. j ) cycle                       ! omit intramodal interactions: --> i .ne. j .
!          if( citable(i,j) .eq. mode_name(k) ) then  ! i-particles and j-particles are lost; k-particles are formed.
!            if( i .ne. k ) dij(i,j) = 1              ! the k-particles are not i-particles (but may be j-particles),
!          endif                                      !   so i-particles are lost by this i-j interaction.
!        enddo
!      enddo
!      enddo
!      xdij = dij
!
!      if( .not. write_tensors ) return
!
!      !-------------------------------------------------------------------------
!      ! write the g_ikl,q.
!      !-------------------------------------------------------------------------
!      do i=1, nmodes
!        write(aunit1,'(/2a)') 'g_iklq for mode ', mode_name(i)
!        do q=1, nm(i)
!          write(aunit1,'(/a,i3,3x,3a5/)') 'q, nm_spc_name(i,q), mode',
!     &                        q, nm_spc_name(i,q), '-->', mode_name(i)
!          if ( sum( giklq(i,1:nmodes,1:nmodes,q) ) .eq. 0 ) cycle
!          write(aunit1,'(5x,16a5)') mode_name(1:nmodes)
!          do k=1, nmodes
!            write(aunit1,'(a5,16i5)') mode_name(k),giklq(i,k,1:nmodes,q)
!          enddo
!        enddo
!      enddo
!
!      !-------------------------------------------------------------------------
!      ! write the d_ikl.
!      !-------------------------------------------------------------------------
!      do i=1, nmodes
!        write(aunit1,'(/2a)') 'd_ikl for mode ', mode_name(i)
!        if ( sum( dikl(i,1:nmodes,1:nmodes) ) .eq. 0 ) cycle
!        write(aunit1,'(5x,16a5)') mode_name(1:nmodes)
!        do k=1, nmodes
!          write(aunit1,'(a5,16i5)') mode_name(k),dikl(i,k,1:nmodes)
!        enddo
!      enddo
!
!      !-------------------------------------------------------------------------
!      ! write the d_ij.
!      !-------------------------------------------------------------------------
!      write(aunit1,'(/2a)') 'd_ij'
!      write(aunit1,'(5x,16a5)') mode_name(1:nmodes)
!      do i=1, nmodes
!        write(aunit1,'(a5,16i5)') mode_name(i),dij(i,1:nmodes)
!      enddo
!
!90000 format(14a4)
!      return
!
!      contains
!
!      subroutine initializediklcontrol(control, mask)
!      type (dikl_type) :: control(:)
!      integer, intent(in) :: mask(:,:,:)
!      integer :: i, k, l, n
!
!      n = 0
!      do k = 1, nweights
!        do l = k+1, nweights
!          do i = 1, nweights
!            if (mask(i,k,l) /= 0) then
!              n = n + 1
!              control(n)%i = i
!              control(n)%k = k
!              control(n)%l = l
!            end if
!          end do
!        end do
!      end do
!      ndikl = n
!      end subroutine initializediklcontrol
!
!      subroutine initializegiklqcontrol(control, mask)
!      type (giklq_type) :: control(:)
!      integer, intent(in) :: mask(:,:,:,:)
!
!      integer :: i, q, k, l, n, ntotal
!
!      do i = 1, nweights
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
!        ntotal = n
!        ! 2) allocate ntotal entries
!        control(i)%n = ntotal
!        allocate(control(i)%k(ntotal))
!        allocate(control(i)%l(ntotal))
!        allocate(control(i)%qq(ntotal))
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
!      end subroutine initializegiklqcontrol
!
!
!      end subroutine setup_coag_tensors
!
