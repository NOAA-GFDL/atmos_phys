#include "cosp_defs.H"
! version number = 1.4.3
! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 09:08:38 -0500 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/llnl/llnl_stats.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list 
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Lawrence Livermore National Security Limited Liability Corporation 
!       nor the names of its contributors may be used to endorse or promote products derived from 
!       this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
!
! Jan 2013 - G. Cesana        - Added betaperp_tot and temp_tot arguments 
!


MODULE MOD_LLNL_STATS
  USE MOD_COSP_CONSTANTS
  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- FUNCTION COSP_CFAD ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION COSP_CFAD(Npoints,Ncolumns,Nlevels,Nbins,x,xmin,xmax,bmin,bwidth)
   ! Input arguments
   integer,intent(in) :: Npoints,Ncolumns,Nlevels,Nbins
   real,dimension(Npoints,Ncolumns,Nlevels),intent(in) :: x
   real,intent(in) :: xmin,xmax 
   real,intent(in) :: bmin,bwidth
   
   real,dimension(Npoints,Nbins,Nlevels) :: cosp_cfad
   ! Local variables
   integer :: i, j, k
   integer :: ibin
   
   !--- Input arguments
   ! Npoints: Number of horizontal points
   ! Ncolumns: Number of subcolumns
   ! Nlevels: Number of levels
   ! Nbins: Number of x axis bins
   ! x: variable to process (Npoints,Ncolumns,Nlevels)
   ! xmin: minimum value allowed for x
   ! xmax: minimum value allowed for x
   ! bmin: mimumum value of first bin
   ! bwidth: bin width
   !
   ! Output: 2D histogram on each horizontal point (Npoints,Nbins,Nlevels)
   
   cosp_cfad = 0.0
   ! bwidth intervals in the range [bmin,bmax=bmin+Nbins*hwidth]
   ! Valid x values smaller than bmin and larger than bmax are set 
   ! into the smallest bin and largest bin, respectively.
   do j = 1, Nlevels, 1
      do k = 1, Ncolumns, 1
         do i = 1, Npoints, 1
            if (x(i,k,j) == R_GROUND) then
               cosp_cfad(i,:,j) = R_UNDEF
            elseif ((x(i,k,j) >= xmin) .and. (x(i,k,j) <= xmax)) then 
               ibin = ceiling((x(i,k,j) - bmin)/bwidth)
               if (ibin > Nbins) ibin = Nbins
               if (ibin < 1)     ibin = 1
               cosp_cfad(i,ibin,j) = cosp_cfad(i,ibin,j) + 1.0 
            end if
         enddo  !i
      enddo  !k
   enddo  !j
   where ((cosp_cfad /= R_UNDEF).and.(cosp_cfad /= 0.0)) cosp_cfad = cosp_cfad / Ncolumns
END FUNCTION COSP_CFAD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_LIDAR_ONLY_CLOUD -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_LIDAR_ONLY_CLOUD(Npoints,Ncolumns,Nlevels,temp_tot,beta_tot, &
                   betaperp_tot,beta_mol,Ze_tot,lidar_only_freq_cloud,tcc)
   ! Input arguments
   integer,intent(in) :: Npoints,Ncolumns,Nlevels
   real,dimension(Npoints,Nlevels),intent(in) :: beta_mol   ! Molecular backscatter
   real,dimension(Npoints,Ncolumns,Nlevels),intent(in) :: beta_tot   ! Total backscattered signal
   real,dimension(Npoints,Ncolumns,Nlevels),intent(in) :: temp_tot   ! Total backscattered signal
   real,dimension(Npoints,Ncolumns,Nlevels),intent(in) :: betaperp_tot   ! perpendicular Total backscattered signal
   real,dimension(Npoints,Ncolumns,Nlevels),intent(in) :: Ze_tot     ! Radar reflectivity
   ! Output arguments
   real,dimension(Npoints,Nlevels),intent(out) :: lidar_only_freq_cloud
   real,dimension(Npoints),intent(out) :: tcc

   ! local variables
   real :: sc_ratio
   real :: s_cld, s_att
   parameter (S_cld = 5.0)
   parameter (s_att = 0.01)
   integer :: flag_sat !first saturated level encountered from top
   integer :: flag_cld !cloudy column
   integer :: pr,i,j

   lidar_only_freq_cloud = 0.0
   tcc = 0.0
   do pr=1,Npoints
     do i=1,Ncolumns
       flag_sat = 0
       flag_cld = 0
       do j=Nlevels,1,-1 !top->surf
        sc_ratio = beta_tot(pr,i,j)/beta_mol(pr,j)
        if ((sc_ratio .le. s_att) .and. (flag_sat .eq. 0)) flag_sat = j
        if (Ze_tot(pr,i,j) .lt. -30.) then  !radar can't detect cloud
         if ( (sc_ratio .gt. s_cld) .or. (flag_sat .eq. j) ) then  !lidar sense cloud
            lidar_only_freq_cloud(pr,j)=lidar_only_freq_cloud(pr,j)+1. !top->surf
            flag_cld=1
         endif
        else  !radar sense cloud (z%Ze_tot(pr,i,j) .ge. -30.)
           flag_cld=1
        endif
       enddo !levels
       if (flag_cld .eq. 1) tcc(pr)=tcc(pr)+1.
     enddo !columns
   enddo !points
   lidar_only_freq_cloud=lidar_only_freq_cloud/Ncolumns
   tcc=tcc/Ncolumns

END SUBROUTINE COSP_LIDAR_ONLY_CLOUD
END MODULE MOD_LLNL_STATS
