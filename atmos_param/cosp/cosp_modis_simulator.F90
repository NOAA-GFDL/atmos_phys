#include "cosp_defs.H"
! version number = 1.4.3
! (c) 2009, Regents of the Unversity of Colorado
!   Author: Robert Pincus, Cooperative Institute for Research in the Environmental Sciences
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 09:08:38 -0500 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/cosp_modis_simulator.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
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

!
! History:
!   May 2009 - Robert Pincus - Initial version
!   Dec 2009 - Robert Pincus - Tiny revisions
!
MODULE MOD_COSP_Modis_Simulator
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  use mod_modis_sim, numModisTauBins      => numTauHistogramBins,      &
                     numModisPressureBins => numPressureHistogramBins, &
                     MODIS_TAU      => nominalTauHistogramCenters,     &
                     MODIS_TAU_BNDS => nominalTauHistogramBoundaries,  &
                     MODIS_PC       => nominalPressureHistogramCenters, &
                     MODIS_PC_BNDS  => nominalPressureHistogramBoundaries                     
  implicit none
  !------------------------------------------------------------------------------------------------
  ! Public type
  !
  ! Summary statistics from MODIS retrievals
  type COSP_MODIS
     ! Dimensions
     integer :: Npoints   ! Number of gridpoints
#ifdef COSP_GFDL
     integer :: Ncolumns  ! Number of columns
#endif
     
     !
     ! Grid means; dimension nPoints
     ! 
     real, dimension(:),       pointer :: & 
       Cloud_Fraction_Total_Mean,       Cloud_Fraction_Water_Mean,       Cloud_Fraction_Ice_Mean,       &
       Cloud_Fraction_High_Mean,        Cloud_Fraction_Mid_Mean,         Cloud_Fraction_Low_Mean,       &
       Optical_Thickness_Total_Mean,    Optical_Thickness_Water_Mean,    Optical_Thickness_Ice_Mean,    &
       Optical_Thickness_Total_LogMean, Optical_Thickness_Water_LogMean, Optical_Thickness_Ice_LogMean, &
                                        Cloud_Particle_Size_Water_Mean,  Cloud_Particle_Size_Ice_Mean,  &
       Cloud_Top_Pressure_Total_Mean,                                                                   &
                                        Liquid_Water_Path_Mean,          Ice_Water_Path_Mean

#ifdef COSP_GFDL
      !  Subcolumn particle sizes, optical thickness and cloud top pressure
      !  dimensions (npoints, ncolumns)
      real, dimension(:,:),       pointer :: &
        Column_Particle_Size, Column_Optical_Thickness, &
        Column_Cloud_Top_Pressure, &
        retrievedPhase
#endif
     !
     ! Also need the ISCCP-type optical thickness/cloud top pressure histogram
     !
     real, dimension(:, :, :), pointer :: Optical_Thickness_vs_Cloud_Top_Pressure
     real, dimension(:, :, :), pointer :: Optical_Thickness_vs_ReffICE
     real, dimension(:, :, :), pointer :: Optical_Thickness_vs_ReffLIQ

  end type COSP_MODIS 
  
contains
  !------------------------------------------------------------------------------------------------
#ifdef COSP_GFDL
  subroutine COSP_Modis_Simulator(gridBox, subCols, subcolHydro, isccpSim, nSunlit, modisSim)
#else
  subroutine COSP_Modis_Simulator(gridBox, subCols, subcolHydro, isccpSim, modisSim)
#endif
    ! Arguments
    type(cosp_gridbox), intent(in   ) :: gridBox     ! Gridbox info
    type(cosp_subgrid), intent(in   ) :: subCols     ! subCol indicators of convective/stratiform 
    type(cosp_sghydro), intent(in   ) :: subcolHydro ! subcol hydrometeor contens
    type(cosp_isccp),   intent(in   ) :: isccpSim    ! ISCCP simulator output
#ifdef COSP_GFDL
    integer,            intent(in   ) :: nSunlit     ! Are there any sunlit points?
#endif
    type(cosp_modis),   intent(  out) :: modisSim    ! MODIS simulator subcol output
    
    ! ------------------------------------------------------------
    ! Local variables 
    !   Leave space only for sunlit points
    
#ifdef COSP_GFDL
    integer :: nPoints, nSubCols, nLevels, i, j, k
#else
    integer :: nPoints, nSubCols, nLevels, nSunlit, i, j, k
#endif
    
#ifdef COSP_GFDL
    ! Grid-mean quanties;  dimensions nPoints, nLevels
    real, &
      dimension(nSunlit,                  gridBox%nLevels) :: &
        temperature, pressureLayers
    real, &
      dimension(nSunlit,                  gridBox%nLevels + 1) :: &
        pressureLevels

    ! Subcol quantities, dimension nPoints, nSubCols, nLevels
    real, &
      dimension(nSunlit, subCols%nColumns, gridBox%nLevels) :: &
         opticalThickness, cloudWater, cloudIce, waterSize, iceSize, &
         liquid_opticalThickness, ice_opticalThickness

    ! Vertically-integrated subcol quantities; dimensions nPoints, nSubcols
    integer, &
      dimension(nSunlit, subCols%nColumns) :: &
        retrievedPhase
    real, &
      dimension(nSunlit, subCols%nColumns) :: &
        isccpTau, isccpCloudTopPressure, retrievedCloudTopPressure, retrievedTau, retrievedSize

    ! Vertically-integrated results
    real, dimension(nSunlit) :: &
        cfTotal, cfLiquid, cfIce,                &
        cfHigh,  cfMid,    cfLow,                &
        meanTauTotal, meanTauLiquid, meanTauIce, &
        meanLogTauTotal, meanLogTauLiquid, meanLogTauIce , &
        meanSizeLiquid, meanSizeIce,             &
        meanCloudTopPressure,                    &
        meanLiquidWaterPath, meanIceWaterPath

    real, dimension(nSunlit, numModisTauBins, numModisPressureBins) :: &
       jointHistogram
    real, dimension(nSunlit, numModisTauBins, numMODISReffIceBins) :: & 
         jointHistogram2
    real, dimension(nSunlit, numModisTauBins, numMODISReffLiqBins) :: & 
         jointHistogram3

    integer, dimension(nSunlit) :: sunlit
    integer, dimension(:), allocatable :: notSunlit
#else
    ! Grid-mean quanties;  dimensions nPoints, nLevels
    real, &
      dimension(count(gridBox%sunlit(:) > 0),                  gridBox%nLevels) :: &
        temperature, pressureLayers
    real, &
      dimension(count(gridBox%sunlit(:) > 0),                  gridBox%nLevels + 1) :: &
        pressureLevels
    
    ! Subcol quantities, dimension nPoints, nSubCols, nLevels 
    real, &
      dimension(count(gridBox%sunlit(:) > 0), subCols%nColumns, gridBox%nLevels) :: & 
        opticalThickness, cloudWater, cloudIce, waterSize, iceSize
    
    ! Vertically-integrated subcol quantities; dimensions nPoints, nSubcols 
    integer, &
      dimension(count(gridBox%sunlit(:) > 0), subCols%nColumns) :: & 
        retrievedPhase
    real, &
      dimension(count(gridBox%sunlit(:) > 0), subCols%nColumns) :: & 
        isccpTau, isccpCloudTopPressure, retrievedCloudTopPressure, retrievedTau, retrievedSize  
    
    ! Vertically-integrated results
    real, dimension(count(gridBox%sunlit(:) > 0)) :: & 
        cfTotal, cfLiquid, cfIce,                &
        cfHigh,  cfMid,    cfLow,                &
        meanTauTotal, meanTauLiquid, meanTauIce, &
        meanLogTauTotal, meanLogTauLiquid, meanLogTauIce , &
        meanSizeLiquid, meanSizeIce,             &
        meanCloudTopPressure,                    &
        meanLiquidWaterPath, meanIceWaterPath
        
    real, dimension(count(gridBox%sunlit(:) > 0), numModisTauBins, numModisPressureBins) :: & 
       jointHistogram
    real, dimension(count(gridBox%sunlit(:) > 0), numModisTauBins, numMODISReffIceBins) :: & 
         jointHistogram2
    real, dimension(count(gridBox%sunlit(:) > 0), numModisTauBins, numMODISReffLiqBins) :: & 
         jointHistogram3
    
    
    integer, dimension(count(gridBox%sunlit(:) >  0)) :: sunlit
    integer, dimension(count(gridBox%sunlit(:) <= 0)) :: notSunlit
#endif
    ! ------------------------------------------------------------
    
    !
    ! Are there any sunlit points? 
    !
#ifdef COSP_GFDL
    allocate(notSunlit(count(gridBox%sunlit(:) <= 0)))
#else
    nSunlit = count(gridBox%sunlit(:) > 0)
#endif
    if(nSunlit > 0) then 
      nLevels  = gridBox%Nlevels
      nPoints  = gridBox%Npoints
      nSubCols = subCols%Ncolumns
      !
      ! This is a vector index indicating which points are sunlit
      !
      sunlit(:)    = pack((/ (i, i = 1, nPoints ) /), mask =       gridBox%sunlit(:) > 0)
      notSunlit(:) = pack((/ (i, i = 1, nPoints ) /), mask = .not. gridBox%sunlit(:) > 0)
               
      !
      ! Copy needed quantities, reversing vertical order and removing points with no sunlight 
      !
      pressureLevels(:, 1) = 0.0 ! Top of model, following ISCCP sim
      temperature(:, :)     = gridBox%T (sunlit(:), nLevels:1:-1) 
      pressureLayers(:, :)  = gridBox%p (sunlit(:), nLevels:1:-1) 
      pressureLevels(:, 2:) = gridBox%ph(sunlit(:), nLevels:1:-1) 
      
      !
      ! Subcolumn properties - first stratiform cloud...
      ! 
#ifdef COSP_GFDL
      if (subCols%cols_input_from_model) then
      ! Use stochastic column taus
      opticalThickness(:, :, :) = subCols%dtau_col(sunlit(:),:,:)
      where(subCols%frac_out(sunlit(:), :, :) == I_LSC)
        cloudWater(:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_LSCLIQ)
        waterSize (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_LSCLIQ)
        cloudIce  (:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_LSCICE)
        iceSize   (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_LSCICE)
      elsewhere
        cloudWater      (:, :, :) = 0.
        cloudIce        (:, :, :) = 0.
        waterSize       (:, :, :) = 0.
        iceSize         (:, :, :) = 0.
      end where

      !
      ! .. then add convective cloud
      !
       where(subCols%frac_out(sunlit(:), :, :) == I_CVC)
        cloudWater(:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_CVCLIQ)
        waterSize (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_CVCLIQ)
        cloudIce  (:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_CVCICE)
        iceSize   (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_CVCICE)
      end where
      else
#endif

      where(subCols%frac_out(sunlit(:), :, :) == I_LSC)
        !opticalThickness(:, :, :) = & 
        !               spread(gridBox%dtau_s      (sunlit(:),    :), dim   = 2, nCopies = nSubCols)
        cloudWater(:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_LSCLIQ)
        waterSize (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_LSCLIQ)
        cloudIce  (:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_LSCICE)
        iceSize   (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_LSCICE)
      elsewhere
!       opticalThickness(:, :, :) = 0.
        cloudWater      (:, :, :) = 0.
        cloudIce        (:, :, :) = 0.
        waterSize       (:, :, :) = 0.
        iceSize         (:, :, :) = 0.
      end where 

      ! Loop version of spread above - intrinsic doesn't work on certain platforms. 
      do k = 1, nLevels
        do j = 1, nSubCols
          do i = 1, nSunlit
            if(subCols%frac_out(sunlit(i), j, k) == I_LSC) then
              opticalThickness(i, j, k) = gridBox%dtau_s(sunlit(i), k)
            else
              opticalThickness(i, j, k) = 0.   
            end if 
          end do 
        end do
      end do

      !
      ! .. then add convective cloud...
      !
      where(subCols%frac_out(sunlit(:), :, :) == I_CVC) 
        !opticalThickness(:, :, :) = &
        !               spread(gridBox%dtau_c(      sunlit(:),    :), dim = 2, nCopies = nSubCols)
        cloudWater(:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_CVCLIQ)
        waterSize (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_CVCLIQ)
        cloudIce  (:, :, :) = subcolHydro%mr_hydro(sunlit(:), :, :, I_CVCICE)
        iceSize   (:, :, :) = subcolHydro%reff    (sunlit(:), :, :, I_CVCICE)
      end where

      ! Loop version of spread above - intrinsic doesn't work on certain platforms. 
      do k = 1, nLevels
        do j = 1, nSubCols
          do i = 1, nSunlit
            if(subCols%frac_out(sunlit(i), j, k) == I_CVC) opticalThickness(i, j, k) = gridBox%dtau_c(sunlit(i), k)
          end do 
        end do
      end do
#ifdef COSP_GFDL
     endif
#endif

      !
      ! Reverse vertical order 
      !
#ifdef COSP_GFDL
      ! Optical thickness need not be reversed when using subcolumn values
#else
      opticalThickness(:, :, :)  = opticalThickness(:, :, nLevels:1:-1)
#endif
      cloudWater      (:, :, :)  = cloudWater      (:, :, nLevels:1:-1)
      waterSize       (:, :, :)  = waterSize       (:, :, nLevels:1:-1)
      cloudIce        (:, :, :)  = cloudIce        (:, :, nLevels:1:-1)
      iceSize         (:, :, :)  = iceSize         (:, :, nLevels:1:-1)
      
      isccpTau(:, :)              = isccpSim%boxtau (sunlit(:), :)
      isccpCloudTopPressure(:, :) = isccpSim%boxptop(sunlit(:), :)
      
      do i = 1, nSunlit
        call modis_L2_simulator(temperature(i, :), pressureLayers(i, :), pressureLevels(i, :),     &
                                opticalThickness(i, :, :), cloudWater(i, :, :), cloudIce(i, :, :), &
                                waterSize(i, :, :), iceSize(i, :, :),                       &
                                isccpTau(i, :), isccpCloudTopPressure(i, :),                &
                                retrievedPhase(i, :), retrievedCloudTopPressure(i, :),      & 
                                retrievedTau(i, :), retrievedSize(i, :))
      end do

      ! DJS2015: Call L3 modis simulator used by cospv2.0
     ! call modis_L3_simulator(retrievedPhase,              &
     !                         retrievedCloudTopPressure,   &
     !                         retrievedTau, retrievedSize, &
     !                         cfTotal,         cfLiquid,         cfIce,         &
     !                         cfHigh,          cfMid,            cfLow,         &
     !                         meanTauTotal,    meanTauLiquid,    meanTauIce,    &
     !                         meanLogTauTotal, meanLogTauLiquid, meanLogTauIce, &
     !                         meanSizeLiquid,   meanSizeIce,   &
     !                         meanCloudTopPressure,                             &
     !                         meanLiquidWaterPath, meanIceWaterPath, &
     !                         jointHistogram)
     call modis_column(nSunlit,nSubcols,retrievedPhase,retrievedCloudTopPressure,   &
                        retrievedTau,retrievedSize,cfTotal,cfLiquid,cfIce,cfHigh,    &
                        cfMid,cfLow,meanTauTotal,meanTauLiquid,meanTauIce,           &
                        meanLogTauTotal,meanLogTauLiquid,meanLogTauIce,              &
                        meanSizeLiquid,meanSizeIce,meanCloudTopPressure,             &
                        meanLiquidWaterPath, meanIceWaterPath,                       &
                        jointHistogram,jointHistogram2,jointHistogram3)
      ! DJS2015: END
      !
      ! Copy results into COSP structure
      !
      modisSim%Cloud_Fraction_Total_Mean(sunlit(:)) = cfTotal(:)
      modisSim%Cloud_Fraction_Water_Mean(sunlit(:)) = cfLiquid
      modisSim%Cloud_Fraction_Ice_Mean  (sunlit(:)) = cfIce
  
      modisSim%Cloud_Fraction_High_Mean(sunlit(:)) = cfHigh
      modisSim%Cloud_Fraction_Mid_Mean (sunlit(:)) = cfMid
      modisSim%Cloud_Fraction_Low_Mean (sunlit(:)) = cfLow
  
      modisSim%Optical_Thickness_Total_Mean(sunlit(:)) = meanTauTotal
      modisSim%Optical_Thickness_Water_Mean(sunlit(:)) = meanTauLiquid
      modisSim%Optical_Thickness_Ice_Mean  (sunlit(:)) = meanTauIce
  
      modisSim%Optical_Thickness_Total_LogMean(sunlit(:)) = meanLogTauTotal
      modisSim%Optical_Thickness_Water_LogMean(sunlit(:)) = meanLogTauLiquid
      modisSim%Optical_Thickness_Ice_LogMean  (sunlit(:)) = meanLogTauIce
  
      modisSim%Cloud_Particle_Size_Water_Mean(sunlit(:)) = meanSizeLiquid
      modisSim%Cloud_Particle_Size_Ice_Mean  (sunlit(:)) = meanSizeIce
  
      modisSim%Cloud_Top_Pressure_Total_Mean(sunlit(:)) = meanCloudTopPressure
  
      modisSim%Liquid_Water_Path_Mean(sunlit(:)) = meanLiquidWaterPath
      modisSim%Ice_Water_Path_Mean   (sunlit(:)) = meanIceWaterPath
      
#ifdef COSP_GFDL
      modisSim%Column_Particle_Size(sunlit(:),:) = RetrievedSize(:,:)
      modisSim%retrievedPhase      (sunlit(:),:) = RetrievedPhase(:,:)
      modisSim%Column_Optical_Thickness(sunlit(:),:) = RetrievedTau(:,:)
      modisSim%Column_Cloud_Top_Pressure(sunlit(:),:) = RetrievedCloudTopPressure(:,:)

#endif
      modisSim%Optical_Thickness_vs_Cloud_Top_Pressure(sunlit(:), 2:numModisTauBins+1, :) = jointHistogram(:, :, :)
      modisSim%Optical_Thickness_vs_ReffICE(sunlit(:),2:numModisTauBins+1,:)              = jointHistogram2(:, :, :)
      modisSim%Optical_Thickness_vs_ReffLIQ(sunlit(:),2:numModisTauBins+1,:)              = jointHistogram3(:, :, :)
      ! 
      ! Reorder pressure bins in joint histogram to go from surface to TOA 
      !
      modisSim%Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = modisSim%Optical_Thickness_vs_Cloud_Top_Pressure(:, :, numModisPressureBins:1:-1)
      if(nSunlit < nPoints) then 
        !
        ! Where it's night and we haven't done the retrievals the values are undefined
        !
        modisSim%Cloud_Fraction_Total_Mean(notSunlit(:)) = R_UNDEF
        modisSim%Cloud_Fraction_Water_Mean(notSunlit(:)) = R_UNDEF
        modisSim%Cloud_Fraction_Ice_Mean  (notSunlit(:)) = R_UNDEF
    
        modisSim%Cloud_Fraction_High_Mean(notSunlit(:)) = R_UNDEF
        modisSim%Cloud_Fraction_Mid_Mean (notSunlit(:)) = R_UNDEF
        modisSim%Cloud_Fraction_Low_Mean (notSunlit(:)) = R_UNDEF

        modisSim%Optical_Thickness_Total_Mean(notSunlit(:)) = R_UNDEF
        modisSim%Optical_Thickness_Water_Mean(notSunlit(:)) = R_UNDEF
        modisSim%Optical_Thickness_Ice_Mean  (notSunlit(:)) = R_UNDEF
    
        modisSim%Optical_Thickness_Total_LogMean(notSunlit(:)) = R_UNDEF
        modisSim%Optical_Thickness_Water_LogMean(notSunlit(:)) = R_UNDEF
        modisSim%Optical_Thickness_Ice_LogMean  (notSunlit(:)) = R_UNDEF
    
        modisSim%Cloud_Particle_Size_Water_Mean(notSunlit(:)) = R_UNDEF
        modisSim%Cloud_Particle_Size_Ice_Mean  (notSunlit(:)) = R_UNDEF
    
        modisSim%Cloud_Top_Pressure_Total_Mean(notSunlit(:)) = R_UNDEF
    
        modisSim%Liquid_Water_Path_Mean(notSunlit(:)) = R_UNDEF
        modisSim%Ice_Water_Path_Mean   (notSunlit(:)) = R_UNDEF
#ifdef COSP_GFDL
        modisSim%Column_Particle_Size(notSunlit(:),:) = R_UNDEF
        modisSim%retrievedPhase(notSunlit(:),:) = R_UNDEF
        modisSim%Column_Optical_Thickness(notSunlit(:),:) = R_UNDEF
        modisSim%Column_Cloud_Top_Pressure(notSunlit(:),:) = R_UNDEF
#endif
  
        modisSim%Optical_Thickness_vs_Cloud_Top_Pressure(notSunlit(:), :, :) = R_UNDEF
        modisSim%Optical_Thickness_vs_ReffICE(notSunlit(:), :, :) = R_UNDEF
        modisSim%Optical_Thickness_vs_ReffLIQ(notSunlit(:), :, :) = R_UNDEF
      end if 
    else
      !
      ! It's nightime everywhere - everything is undefined
      !
      modisSim%Cloud_Fraction_Total_Mean(:) = R_UNDEF
      modisSim%Cloud_Fraction_Water_Mean(:) = R_UNDEF
      modisSim%Cloud_Fraction_Ice_Mean  (:) = R_UNDEF
  
      modisSim%Cloud_Fraction_High_Mean(:) = R_UNDEF
      modisSim%Cloud_Fraction_Mid_Mean (:) = R_UNDEF
      modisSim%Cloud_Fraction_Low_Mean (:) = R_UNDEF

      modisSim%Optical_Thickness_Total_Mean(:) = R_UNDEF
      modisSim%Optical_Thickness_Water_Mean(:) = R_UNDEF
      modisSim%Optical_Thickness_Ice_Mean  (:) = R_UNDEF
  
      modisSim%Optical_Thickness_Total_LogMean(:) = R_UNDEF
      modisSim%Optical_Thickness_Water_LogMean(:) = R_UNDEF
      modisSim%Optical_Thickness_Ice_LogMean  (:) = R_UNDEF
  
      modisSim%Cloud_Particle_Size_Water_Mean(:) = R_UNDEF
      modisSim%Cloud_Particle_Size_Ice_Mean  (:) = R_UNDEF
  
      modisSim%Cloud_Top_Pressure_Total_Mean(:) = R_UNDEF
  
      modisSim%Liquid_Water_Path_Mean(:) = R_UNDEF
      modisSim%Ice_Water_Path_Mean   (:) = R_UNDEF

#ifdef COSP_GFDL
      modisSim%Column_Particle_Size(:,:) = R_UNDEF
      modisSim%retrievedPhase      (:,:) = R_UNDEF
      modisSim%Column_Optical_Thickness(:,:) = R_UNDEF
      modisSim%Column_Cloud_Top_Pressure(:,:) = R_UNDEF
#endif

      modisSim%Optical_Thickness_vs_Cloud_Top_Pressure(:, :, :) = R_UNDEF
      modisSim%Optical_Thickness_vs_ReffICE(:, :, :) = R_UNDEF
      modisSim%Optical_Thickness_vs_ReffLIQ(:, :, :) = R_UNDEF
    end if 

#ifdef COSP_GFDL
    deallocate (notSunlit)
#endif

  end subroutine COSP_Modis_Simulator
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------- SUBROUTINE CONSTRUCT_COSP_MODIS ------------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef COSP_GFDL
  SUBROUTINE CONSTRUCT_COSP_MODIS(cfg, nPoints, nColumns, x)
#else
  SUBROUTINE CONSTRUCT_COSP_MODIS(cfg, nPoints, x)
#endif
    type(cosp_config), intent(in)  :: cfg ! Configuration options
    integer,           intent(in)  :: Npoints  ! Number of sampled points
#ifdef COSP_GFDL
    integer,           intent(in)  :: Ncolumns  ! Number of subgrid columns
#endif
    type(cosp_MODIS),  intent(out) :: x
    !
    ! Allocate minumum storage if simulator not used
    !
    if (cfg%LMODIS_sim) then
      x%nPoints  = nPoints
#ifdef COSP_GFDL
      x%nColumns  = nColumns
#endif
    else
      x%Npoints  = 1
#ifdef COSP_GFDL
      x%Ncolumns  = 1
#endif
    endif
    
    ! --- Allocate arrays ---
    allocate(x%Cloud_Fraction_Total_Mean(x%nPoints)) 
    allocate(x%Cloud_Fraction_Water_Mean(x%nPoints)) 
    allocate(x%Cloud_Fraction_Ice_Mean(x%nPoints)) 
    
    allocate(x%Cloud_Fraction_High_Mean(x%nPoints)) 
    allocate(x%Cloud_Fraction_Mid_Mean(x%nPoints)) 
    allocate(x%Cloud_Fraction_Low_Mean(x%nPoints)) 
    
    allocate(x%Optical_Thickness_Total_Mean(x%nPoints)) 
    allocate(x%Optical_Thickness_Water_Mean(x%nPoints)) 
    allocate(x%Optical_Thickness_Ice_Mean(x%nPoints)) 
    
    allocate(x%Optical_Thickness_Total_LogMean(x%nPoints)) 
    allocate(x%Optical_Thickness_Water_LogMean(x%nPoints)) 
    allocate(x%Optical_Thickness_Ice_LogMean(x%nPoints)) 
    
    allocate(x%Cloud_Particle_Size_Water_Mean(x%nPoints)) 
    allocate(x%Cloud_Particle_Size_Ice_Mean(x%nPoints)) 
    
    allocate(x%Cloud_Top_Pressure_Total_Mean(x%nPoints)) 
    
    allocate(x%Liquid_Water_Path_Mean(x%nPoints)) 
    allocate(x%Ice_Water_Path_Mean(x%nPoints)) 
#ifdef COSP_GFDL
    allocate(x%Column_Particle_Size(x%nPoints,x%nColumns))
    allocate(x%retrievedPhase      (x%nPoints,x%nColumns))
    allocate(x%Column_Optical_Thickness(x%nPoints,x%nColumns))
    allocate(x%Column_Cloud_Top_Pressure(x%nPoints,x%nColumns))
#endif
      
    allocate(x%Optical_Thickness_vs_Cloud_Top_Pressure(nPoints, numModisTauBins+1, numModisPressureBins))
    allocate(x%Optical_Thickness_vs_ReffICE(nPoints, numModisTauBins+1, numModisReffIceBins))
    allocate(x%Optical_Thickness_vs_ReffLIQ(nPoints, numModisTauBins+1, numModisReffLiqBins))
    x%Optical_Thickness_vs_Cloud_Top_Pressure(:, :, :) = R_UNDEF
    x%Optical_Thickness_vs_ReffICE(:, :, :)            = R_UNDEF
    x%Optical_Thickness_vs_ReffLIQ(:, :, :)            = R_UNDEF
  END SUBROUTINE CONSTRUCT_COSP_MODIS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------- SUBROUTINE FREE_COSP_MODIS -----------------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE FREE_COSP_MODIS(x)
    type(cosp_MODIS),intent(inout) :: x
    !
    ! Free space used by cosp_modis variable. 
    !
    
    if(associated(x%Cloud_Fraction_Total_Mean)) deallocate(x%Cloud_Fraction_Total_Mean) 
    if(associated(x%Cloud_Fraction_Water_Mean)) deallocate(x%Cloud_Fraction_Water_Mean) 
    if(associated(x%Cloud_Fraction_Ice_Mean  )) deallocate(x%Cloud_Fraction_Ice_Mean) 
    
    if(associated(x%Cloud_Fraction_High_Mean)) deallocate(x%Cloud_Fraction_High_Mean) 
    if(associated(x%Cloud_Fraction_Mid_Mean )) deallocate(x%Cloud_Fraction_Mid_Mean) 
    if(associated(x%Cloud_Fraction_Low_Mean )) deallocate(x%Cloud_Fraction_Low_Mean) 
    
    if(associated(x%Optical_Thickness_Total_Mean)) deallocate(x%Optical_Thickness_Total_Mean) 
    if(associated(x%Optical_Thickness_Water_Mean)) deallocate(x%Optical_Thickness_Water_Mean) 
    if(associated(x%Optical_Thickness_Ice_Mean  )) deallocate(x%Optical_Thickness_Ice_Mean) 
    
    if(associated(x%Optical_Thickness_Total_LogMean)) deallocate(x%Optical_Thickness_Total_LogMean) 
    if(associated(x%Optical_Thickness_Water_LogMean)) deallocate(x%Optical_Thickness_Water_LogMean) 
    if(associated(x%Optical_Thickness_Ice_LogMean  )) deallocate(x%Optical_Thickness_Ice_LogMean) 
    
    if(associated(x%Cloud_Particle_Size_Water_Mean)) deallocate(x%Cloud_Particle_Size_Water_Mean) 
    if(associated(x%Cloud_Particle_Size_Ice_Mean  )) deallocate(x%Cloud_Particle_Size_Ice_Mean) 
    
    if(associated(x%Cloud_Top_Pressure_Total_Mean )) deallocate(x%Cloud_Top_Pressure_Total_Mean   ) 
    
    if(associated(x%Liquid_Water_Path_Mean)) deallocate(x%Liquid_Water_Path_Mean   ) 
    if(associated(x%Ice_Water_Path_Mean   )) deallocate(x%Ice_Water_Path_Mean   ) 
    
#ifdef COSP_GFDL
    if(associated(x%Column_Particle_Size  )) deallocate(x%Column_Particle_Size  )
    if(associated(x%retrievedPhase        )) deallocate(x%retrievedPhase )
    if(associated(x%Column_Optical_Thickness)) deallocate(x%Column_Optical_Thickness)
    if(associated(x%Column_Cloud_Top_Pressure)) deallocate(x%Column_Cloud_Top_Pressure  )
#endif

    if(associated(x%Optical_Thickness_vs_Cloud_Top_Pressure)) deallocate(x%Optical_Thickness_vs_Cloud_Top_Pressure   ) 
    if(associated(x%Optical_Thickness_vs_ReffIce)) deallocate(x%Optical_Thickness_vs_ReffIce) 
    if(associated(x%Optical_Thickness_vs_ReffLiq)) deallocate(x%Optical_Thickness_vs_ReffLiq) 
  END SUBROUTINE FREE_COSP_MODIS
  ! -----------------------------------------------------

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------- SUBROUTINE COSP_MODIS_CPSECTION -----------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_MODIS_CPSECTION(ix, iy, orig, copy)
    integer, dimension(2), intent(in) :: ix, iy
    type(cosp_modis),      intent(in   ) :: orig
    type(cosp_modis),      intent(  out) :: copy
    !
    ! Copy a set of grid points from one cosp_modis variable to another.
    !   Should test to be sure ix and iy refer to the same number of grid points 
    !
    integer :: orig_start, orig_end, copy_start, copy_end
    
    orig_start = ix(1); orig_end = ix(2)
    copy_start = iy(1); copy_end = iy(2) 
    
    copy%Cloud_Fraction_Total_Mean(copy_start:copy_end) = orig%Cloud_Fraction_Total_Mean(orig_start:orig_end)
    copy%Cloud_Fraction_Water_Mean(copy_start:copy_end) = orig%Cloud_Fraction_Water_Mean(orig_start:orig_end)
    copy%Cloud_Fraction_Ice_Mean  (copy_start:copy_end) = orig%Cloud_Fraction_Ice_Mean  (orig_start:orig_end)
    
    copy%Cloud_Fraction_High_Mean(copy_start:copy_end) = orig%Cloud_Fraction_High_Mean(orig_start:orig_end)
    copy%Cloud_Fraction_Mid_Mean (copy_start:copy_end) = orig%Cloud_Fraction_Mid_Mean (orig_start:orig_end)
    copy%Cloud_Fraction_Low_Mean (copy_start:copy_end) = orig%Cloud_Fraction_Low_Mean (orig_start:orig_end)
    
    copy%Optical_Thickness_Total_Mean(copy_start:copy_end) = orig%Optical_Thickness_Total_Mean(orig_start:orig_end)
    copy%Optical_Thickness_Water_Mean(copy_start:copy_end) = orig%Optical_Thickness_Water_Mean(orig_start:orig_end)
    copy%Optical_Thickness_Ice_Mean  (copy_start:copy_end) = orig%Optical_Thickness_Ice_Mean  (orig_start:orig_end)
    
    copy%Optical_Thickness_Total_LogMean(copy_start:copy_end) = &
                                                          orig%Optical_Thickness_Total_LogMean(orig_start:orig_end)
    copy%Optical_Thickness_Water_LogMean(copy_start:copy_end) = &
                                                          orig%Optical_Thickness_Water_LogMean(orig_start:orig_end)
    copy%Optical_Thickness_Ice_LogMean  (copy_start:copy_end) = &
                                                          orig%Optical_Thickness_Ice_LogMean  (orig_start:orig_end)

    copy%Cloud_Particle_Size_Water_Mean(copy_start:copy_end) = orig%Cloud_Particle_Size_Water_Mean(orig_start:orig_end)
    copy%Cloud_Particle_Size_Ice_Mean  (copy_start:copy_end) = orig%Cloud_Particle_Size_Ice_Mean  (orig_start:orig_end)

    copy%Cloud_Top_Pressure_Total_Mean(copy_start:copy_end) = orig%Cloud_Top_Pressure_Total_Mean(orig_start:orig_end)
    
    copy%Liquid_Water_Path_Mean(copy_start:copy_end) = orig%Liquid_Water_Path_Mean(orig_start:orig_end)
    copy%Ice_Water_Path_Mean   (copy_start:copy_end) = orig%Ice_Water_Path_Mean  (orig_start:orig_end)
    
#ifdef COSP_GFDL
    copy%Column_Particle_Size  (copy_start:copy_end,:) = orig%Column_Particle_Size (orig_start:orig_end,:)
    copy%retrievedPhase        (copy_start:copy_end,:) = orig%retrievedPhase       (orig_start:orig_end,:)
    copy%Column_Optical_Thickness  (copy_start:copy_end,:) = orig%Column_Optical_Thickness (orig_start:orig_end,:)
    copy%Column_Cloud_Top_Pressure  (copy_start:copy_end,:) = orig%Column_Cloud_Top_Pressure (orig_start:orig_end,:)
#endif

    copy%Optical_Thickness_vs_Cloud_Top_Pressure(copy_start:copy_end, :, :) = &
                          orig%Optical_Thickness_vs_Cloud_Top_Pressure(orig_start:orig_end, :, :)
    copy%Optical_Thickness_vs_ReffIce(copy_start:copy_end, :, :) = &
         orig%Optical_Thickness_vs_ReffIce(orig_start:orig_end, :, :)
    copy%Optical_Thickness_vs_ReffLiq(copy_start:copy_end, :, :) = &
         orig%Optical_Thickness_vs_ReffLiq(orig_start:orig_end, :, :)

  END SUBROUTINE COSP_MODIS_CPSECTION
  ! -----------------------------------------------------

END MODULE MOD_COSP_Modis_Simulator
