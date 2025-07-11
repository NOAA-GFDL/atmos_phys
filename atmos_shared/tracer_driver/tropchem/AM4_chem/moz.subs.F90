      module mo_setrxt_mod
      private
      public :: setrxt
      contains
      subroutine setrxt( rate, temp, m, plonl, plev, plnplv )
      use CHEM_MODS_MOD, only : rxntot
      use mo_jpl_mod, only : jpl
      implicit none
!-------------------------------------------------------
! ... Dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: plonl, plev, plnplv
      real, intent(in) :: temp(plonl,plev), m(plonl,plev)
      real, intent(inout) :: rate(plonl,plev,rxntot)
!-------------------------------------------------------
! ... Local variables
!-------------------------------------------------------
      real :: itemp(plonl,plev), exp_fac(plonl,plev)
      real, dimension(plonl,plev) :: ko, kinf
      rate(:,:,55) = 1.2e-10
      rate(:,:,61) = 1.8e-12
      rate(:,:,63) = 1.8e-12
      rate(:,:,74) = 3.5e-12
      rate(:,:,77) = 0.
      rate(:,:,85) = 1.5e-10
      rate(:,:,95) = 1.e-14
      rate(:,:,116) = 3.e-13
      rate(:,:,117) = 4.1e-14
      rate(:,:,133) = 8.37E-14
      rate(:,:,134) = 1.54E-13
      rate(:,:,138) = 3.7E-19
      rate(:,:,148) = 8.37E-14
      rate(:,:,150) = 1.6E-12
      rate(:,:,154) = 3.4E-15
      rate(:,:,157) = 8.37E-14
      rate(:,:,160) = 3.20E-12
      rate(:,:,172) = 8.37E-14
      rate(:,:,179) = 2.90E-11
      rate(:,:,180) = 1.00E-11
      rate(:,:,182) = 4.0E-16
      rate(:,:,183) = 1.50E-11
      rate(:,:,188) = 2.30E-12
      rate(:,:,193) = 4.00E-12
      rate(:,:,198) = 1.6e-12
      rate(:,:,207) = 6.5e-15
      rate(:,:,211) = 9.2e-13
      rate(:,:,215) = 4e-13
      rate(:,:,216) = 1.2e-11
      rate(:,:,217) = 1.1e-11
      rate(:,:,218) = 2.2e-11
      rate(:,:,221) = 3.0e-13
      rate(:,:,223) = 5e-11
      rate(:,:,236) = 2e-12
      rate(:,:,237) = 1e-11
      rate(:,:,238) = 1e-11
      rate(:,:,241) = 1.11e-11
      rate(:,:,304) = 1.16e-5
      rate(:,:,305) = 1.16e-5
      itemp(:,:) = 1. / temp(:,:)
      rate(:,:,50) = 8e-12 * exp( -2060. * itemp(:,:) )
      rate(:,:,51) = 1.5e-11 * exp( -3600. * itemp(:,:) )
      rate(:,:,52) = 2.1e-11 * exp( 100. * itemp(:,:) )
      exp_fac(:,:) = exp( 180. * itemp(:,:) )
      rate(:,:,56) = 1.8e-11 * exp_fac(:,:)
      rate(:,:,94) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,121) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,203) = 4.2e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 200. * itemp(:,:) )
      rate(:,:,57) = 3e-11 * exp_fac(:,:)
      rate(:,:,90) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,103) = 8.78E-12 * exp_fac(:,:)
      rate(:,:,110) = 3.80e-12 * exp_fac(:,:)
      rate(:,:,118) = 5.18e-12 * exp_fac(:,:)
      rate(:,:,124) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,139) = 4.75E-12 * exp_fac(:,:)
      rate(:,:,149) = 8.78E-12 * exp_fac(:,:)
      rate(:,:,158) = 1.84E-12 * exp_fac(:,:)
      rate(:,:,169) = 6.13E-13 * exp_fac(:,:)
      rate(:,:,177) = 2.66e-12 * exp_fac(:,:)
      rate(:,:,178) = 1.14e-12 * exp_fac(:,:)
      rate(:,:,191) = 5.18E-12 * exp_fac(:,:)
      rate(:,:,205) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,288) = 5.5e-12 * exp_fac(:,:)
      rate(:,:,58) = 1.7e-12 * exp( -940. * itemp(:,:) )
      rate(:,:,59) = 1e-14 * exp( -490. * itemp(:,:) )
      rate(:,:,62) = 4.8e-11 * exp( 250. * itemp(:,:) )
      rate(:,:,64) = 2.8e-12 * exp( -1800. * itemp(:,:) )
      rate(:,:,65) = 2.15e-11 * exp( 110. * itemp(:,:) )
      rate(:,:,66) = 3.3e-11 * exp( 55. * itemp(:,:) )
      rate(:,:,67) = 1.63e-10 * exp( 60. * itemp(:,:) )
      exp_fac(:,:) = exp( 20. * itemp(:,:) )
      rate(:,:,68) = 7.25e-11 * exp_fac(:,:)
      rate(:,:,69) = 4.63e-11 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270. * itemp(:,:) )
      rate(:,:,70) = 3.3e-12 * exp_fac(:,:)
      rate(:,:,106) = 8.1e-12 * exp_fac(:,:)
      rate(:,:,275) = 7.4e-12 * exp_fac(:,:)
      rate(:,:,71) = 3e-12 * exp( -1500. * itemp(:,:) )
      rate(:,:,72) = 5.1e-12 * exp( 210. * itemp(:,:) )
      exp_fac(:,:) = exp( -2450. * itemp(:,:) )
      rate(:,:,73) = 1.2e-13 * exp_fac(:,:)
      rate(:,:,292) = 8.5e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170. * itemp(:,:) )
      rate(:,:,80) = 1.5e-11 * exp_fac(:,:)
      rate(:,:,273) = 1.8e-11 * exp_fac(:,:)
      exp_fac(:,:) = exp( 380. * itemp(:,:) )
      rate(:,:,82) = 1.3e-12 * exp_fac(:,:)
      rate(:,:,135) = 3.61E-12 * exp_fac(:,:)
      rate(:,:,151) = 8.00E-12 * exp_fac(:,:)
      rate(:,:,159) = 4.40E-12 * exp_fac(:,:)
      rate(:,:,170) = 3.60E-12 * exp_fac(:,:)
      rate(:,:,84) = 2.45e-12 * exp( -1775. * itemp(:,:) )
      exp_fac(:,:) = exp( 300. * itemp(:,:) )
      rate(:,:,86) = 2.8e-12 * exp_fac(:,:)
      rate(:,:,175) = 2.80E-12 * exp_fac(:,:)
      rate(:,:,87) = 6.03e-13 * exp( -453. * itemp(:,:) )
      rate(:,:,88) = 2.30e-14 * exp( 677. * itemp(:,:) )
      rate(:,:,89) = 4.1e-13 * exp( 750. * itemp(:,:) )
      exp_fac(:,:) = exp( -1900. * itemp(:,:) )
      rate(:,:,91) = 3.4e-13 * exp_fac(:,:)
      rate(:,:,105) = 1.4e-12 * exp_fac(:,:)
      rate(:,:,92) = 5.5e-12 * exp( 125. * itemp(:,:) )
      rate(:,:,96) = 1.6e11 * exp( -4150. * itemp(:,:) )
      rate(:,:,97) = 1.2e-14 * exp( -2630. * itemp(:,:) )
      rate(:,:,99) = 5.50E-15 * exp( -1880. * itemp(:,:) )
      rate(:,:,100) = 4.6e-13 * exp( -1156. * itemp(:,:) )
      exp_fac(:,:) = exp( 350. * itemp(:,:) )
      rate(:,:,101) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,104) = 4.63E-12 * exp_fac(:,:)
      rate(:,:,128) = 3.10E-11 * exp_fac(:,:)
      rate(:,:,131) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,143) = 2.70E-12 * exp_fac(:,:)
      rate(:,:,146) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,155) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,161) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,173) = 2.35E-12 * exp_fac(:,:)
      rate(:,:,174) = 0.35E-12 * exp_fac(:,:)
      rate(:,:,187) = 2.70E-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 700. * itemp(:,:) )
      rate(:,:,102) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,115) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,122) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,176) = 8.60E-13 * exp_fac(:,:)
      rate(:,:,204) = 7.5e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 980. * itemp(:,:) )
      rate(:,:,108) = 5.2e-13 * exp_fac(:,:)
      rate(:,:,195) = 5.20E-13 * exp_fac(:,:)
      rate(:,:,209) = 5.2e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500. * itemp(:,:) )
      rate(:,:,109) = 2.0e-12 * exp_fac(:,:)
      rate(:,:,112) = 2.5e-12 * exp_fac(:,:)
      rate(:,:,165) = 1.68E-12 * exp_fac(:,:)
      rate(:,:,166) = 1.87E-13 * exp_fac(:,:)
      rate(:,:,113) = 7.66e-12 * exp( -1020. * itemp(:,:) )
      rate(:,:,114) = 2.6e-12 * exp( 365. * itemp(:,:) )
      rate(:,:,119) = 1.55e-11 * exp( -540. * itemp(:,:) )
      rate(:,:,120) = 8.7e-12 * exp( -615. * itemp(:,:) )
      rate(:,:,123) = 3.75e-13 * exp( -40. * itemp(:,:) )
      rate(:,:,126) = 2.9e-12 * exp( -345. * itemp(:,:) )
      rate(:,:,127) = 6.9e-12 * exp( -230. * itemp(:,:) )
      rate(:,:,129) = 4.07E+08 * exp( -7694. * itemp(:,:) )
      rate(:,:,130) = 1.00e-14 * exp( -1970. * itemp(:,:) )
      exp_fac(:,:) = exp( 1300. * itemp(:,:) )
      rate(:,:,132) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,142) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,147) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,156) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,162) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,171) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,189) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,234) = 1.13e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 360. * itemp(:,:) )
      rate(:,:,136) = 2.4e-12 * exp_fac(:,:)
      rate(:,:,210) = 9.9e-12 * exp_fac(:,:)
      rate(:,:,137) = 8.7e-14 * exp( 1650. * itemp(:,:) )
      exp_fac(:,:) = exp( 390. * itemp(:,:) )
      rate(:,:,140) = 1.90E-11 * exp_fac(:,:)
      rate(:,:,190) = 1.90E-11 * exp_fac(:,:)
      rate(:,:,141) = 5.78E-11 * exp( -400. * itemp(:,:) )
      rate(:,:,144) = 2.6E-12 * exp( 610. * itemp(:,:) )
      exp_fac(:,:) = exp( -1520. * itemp(:,:) )
      rate(:,:,145) = 8.5E-16 * exp_fac(:,:)
      rate(:,:,196) = 4.15E-15 * exp_fac(:,:)
      rate(:,:,152) = 2.90E+07 * exp( -5297. * itemp(:,:) )
      rate(:,:,153) = 1.40E-15 * exp( -2100. * itemp(:,:) )
      exp_fac(:,:) = exp( 340. * itemp(:,:) )
      rate(:,:,163) = 6.7E-12 * exp_fac(:,:)
      rate(:,:,181) = 3.1E-12 * exp_fac(:,:)
      rate(:,:,194) = 6.70E-12 * exp_fac(:,:)
      rate(:,:,208) = 6.7e-12 * exp_fac(:,:)
      rate(:,:,164) = 4.3E-13 * exp( 1040. * itemp(:,:) )
      rate(:,:,184) = 1.40E-12 * exp( -1860. * itemp(:,:) )
      rate(:,:,185) = 1.60E-12 * exp( 305. * itemp(:,:) )
      rate(:,:,186) = 3.30E-12 * exp( -450. * itemp(:,:) )
      rate(:,:,192) = 3.15E-13 * exp( -448. * itemp(:,:) )
      rate(:,:,197) = 7.48E-12 * exp( 410. * itemp(:,:) )
      rate(:,:,199) = 1.2e-11 * exp( 440. * itemp(:,:) )
      rate(:,:,200) = 5.3e-16 * exp( -530. * itemp(:,:) )
      rate(:,:,201) = 1.2e-12 * exp( 490. * itemp(:,:) )
      rate(:,:,202) = 2.3e-12 * exp( -170. * itemp(:,:) )
      rate(:,:,206) = 1.5e-12 * exp( -90. * itemp(:,:) )
      rate(:,:,212) = 6.0e-11 * exp( 240. * itemp(:,:) )
      rate(:,:,213) = 1.15e-12 * exp( 430. * itemp(:,:) )
      rate(:,:,214) = 1.2e-16 * exp( 1580. * itemp(:,:) )
      rate(:,:,219) = 5.6e16 * exp( -10870. * itemp(:,:) )
      rate(:,:,220) = 3.5e10 * exp( -3560. * itemp(:,:) )
      rate(:,:,222) = 5e13 * exp( -9673. * itemp(:,:) )
      rate(:,:,224) = 5e13 * exp( -9946. * itemp(:,:) )
      rate(:,:,231) = 1.2e-11 * exp( -280. * itemp(:,:) )
      rate(:,:,233) = 1.90e-13 * exp( 520. * itemp(:,:) )
      exp_fac(:,:) = exp( 260. * itemp(:,:) )
      rate(:,:,235) = 4.9e-12 * exp_fac(:,:)
      rate(:,:,282) = 2.3e-12 * exp_fac(:,:)
      rate(:,:,284) = 8.8e-12 * exp_fac(:,:)
      rate(:,:,263) = 1.7e-12 * exp( -710. * itemp(:,:) )
      rate(:,:,264) = 1.4e-10 * exp( -470. * itemp(:,:) )
      rate(:,:,266) = 2.3e-11 * exp( -200. * itemp(:,:) )
      rate(:,:,267) = 2.8e-11 * exp( 85. * itemp(:,:) )
      exp_fac(:,:) = exp( 290. * itemp(:,:) )
      rate(:,:,268) = 6.4e-12 * exp_fac(:,:)
      rate(:,:,289) = 4.1e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800. * itemp(:,:) )
      rate(:,:,270) = 2.9e-12 * exp_fac(:,:)
      rate(:,:,280) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,287) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,271) = 7.3e-12 * exp( -1280. * itemp(:,:) )
      rate(:,:,272) = 2.6e-12 * exp( -350. * itemp(:,:) )
      exp_fac(:,:) = exp( 220. * itemp(:,:) )
      rate(:,:,274) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,294) = 5.8e-12 * exp_fac(:,:)
      rate(:,:,276) = 8.1e-11 * exp( -30. * itemp(:,:) )
      rate(:,:,283) = 4.5e-12 * exp( 460. * itemp(:,:) )
      rate(:,:,285) = 1.2e-10 * exp( -430. * itemp(:,:) )
      rate(:,:,286) = 4.8e-12 * exp( -310. * itemp(:,:) )
      rate(:,:,290) = 6.0e-13 * exp( 230. * itemp(:,:) )
      rate(:,:,291) = 4.5e-14 * exp( -1260. * itemp(:,:) )
      itemp(:,:) = 300. * itemp(:,:)
      ko(:,:) = 5.9e-33 * itemp(:,:)**1.4
      kinf(:,:) = 1.1e-12 * itemp(:,:)**-1.3
      call jpl( rate(1,1,53), m, .6, ko, kinf, plnplv )
      ko(:,:) = 1.5e-13 * itemp(:,:)**-0.6
      kinf(:,:) = 2.1e9 * itemp(:,:)**-6.1
      call jpl( rate(1,1,54), m, .6, ko, kinf, plnplv )
      ko(:,:) = 2.e-30 * itemp(:,:)**4.4
      kinf(:,:) = 1.4e-12 * itemp(:,:)**.7
      call jpl( rate(1,1,75), m, .6, ko, kinf, plnplv )
      ko(:,:) = 1.8e-30 * itemp(:,:)**3.0
      kinf(:,:) = 2.8e-11
      call jpl( rate(1,1,78), m, .6, ko, kinf, plnplv )
      ko(:,:) = 2.0e-31 * itemp(:,:)**3.4
      kinf(:,:) = 2.9e-12 * itemp(:,:)**1.1
      call jpl( rate(1,1,81), m, .6, ko, kinf, plnplv )
      ko(:,:) = 1.e-28 * itemp(:,:)**4.5
      kinf(:,:) = 7.5e-12 * itemp(:,:)**0.85
      call jpl( rate(1,1,93), m, .6, ko, kinf, plnplv )
      ko(:,:) = 8.e-27 * itemp(:,:)**3.5
      kinf(:,:) = 3.e-11
      call jpl( rate(1,1,98), m, .5, ko, kinf, plnplv )
      ko(:,:) = 9.7e-29 * itemp(:,:)**5.6
      kinf(:,:) = 9.3e-12 * itemp(:,:)**1.5
      call jpl( rate(1,1,107), m, .6, ko, kinf, plnplv )
      ko(:,:) = 9.0e-28 * itemp(:,:)**8.9
      kinf(:,:) = 7.7e-12 * itemp(:,:)**.2
      call jpl( rate(1,1,167), m, .6, ko, kinf, plnplv )
      ko(:,:) = 3.3e-31 * itemp(:,:)**4.3
      kinf(:,:) = 1.6e-12
      call jpl( rate(1,1,230), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 4.4e-32 * itemp(:,:)**1.3
      kinf(:,:) = 4.7e-11 * itemp(:,:)**0.2
      call jpl( rate(1,1,265), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 1.8e-31 * itemp(:,:)**3.4
      kinf(:,:) = 1.5e-11 * itemp(:,:)**1.9
      call jpl( rate(1,1,269), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 6.9e-31 * itemp(:,:)**1.0
      kinf(:,:) = 2.6e-11
      call jpl( rate(1,1,277), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 1.6e-32 * itemp(:,:)**4.5
      kinf(:,:) = 2.0e-12 * itemp(:,:)**2.4
      call jpl( rate(1,1,278), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 5.2e-31 * itemp(:,:)**3.2
      kinf(:,:) = 6.9e-12 * itemp(:,:)**2.9
      call jpl( rate(1,1,281), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 9.0e-32 * itemp(:,:)**1.5
      kinf(:,:) = 3.0e-11
      call jpl( rate(1,1,293), m, 0.6, ko, kinf, plnplv )
      end subroutine setrxt
      end module mo_setrxt_mod
      module mo_adjrxt_mod
      private
      public :: adjrxt
      contains
      subroutine adjrxt( rate, inv, m, plnplv )
      use CHEM_MODS_MOD, only : nfs, rxntot
      implicit none
!--------------------------------------------------------------------
! ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in) :: inv(plnplv,nfs)
      real, intent(in) :: m(plnplv)
      real, intent(inout) :: rate(plnplv,rxntot)
!--------------------------------------------------------------------
! ... Local variables
!--------------------------------------------------------------------
      real :: im(plnplv)
      rate(:, 51) = rate(:, 51) * inv(:, 3)
      rate(:, 53) = rate(:, 53) * inv(:, 1)
      rate(:, 65) = rate(:, 65) * inv(:, 2)
      rate(:, 66) = rate(:, 66) * inv(:, 3)
      rate(:, 75) = rate(:, 75) * inv(:, 1)
      rate(:, 76) = rate(:, 76) * inv(:, 1)
      rate(:, 78) = rate(:, 78) * inv(:, 1)
      rate(:, 81) = rate(:, 81) * inv(:, 1)
      rate(:, 83) = rate(:, 83) * inv(:, 1)
      rate(:, 93) = rate(:, 93) * inv(:, 1)
      rate(:, 95) = rate(:, 95) * inv(:, 3)
      rate(:, 98) = rate(:, 98) * inv(:, 1)
      rate(:,107) = rate(:,107) * inv(:, 1)
      rate(:,111) = rate(:,111) * inv(:, 1)
      rate(:,167) = rate(:,167) * inv(:, 1)
      rate(:,214) = rate(:,214) * inv(:, 3)
      rate(:,230) = rate(:,230) * inv(:, 1)
      rate(:,269) = rate(:,269) * inv(:, 1)
      rate(:,277) = rate(:,277) * inv(:, 1)
      rate(:,278) = rate(:,278) * inv(:, 1)
      rate(:,279) = rate(:,279) * inv(:, 1)
      rate(:,281) = rate(:,281) * inv(:, 1)
      rate(:,293) = rate(:,293) * inv(:, 1)
      rate(:, 49) = rate(:, 49) * inv(:, 3) * inv(:, 1)
      rate(:,265) = rate(:,265) * inv(:, 3) * inv(:, 1)
      rate(:, 50) = rate(:, 50) * m(:)
      rate(:, 52) = rate(:, 52) * m(:)
      rate(:, 53) = rate(:, 53) * m(:)
      rate(:, 54) = rate(:, 54) * m(:)
      rate(:, 55) = rate(:, 55) * m(:)
      rate(:, 56) = rate(:, 56) * m(:)
      rate(:, 57) = rate(:, 57) * m(:)
      rate(:, 58) = rate(:, 58) * m(:)
      rate(:, 59) = rate(:, 59) * m(:)
      rate(:, 60) = rate(:, 60) * m(:)
      rate(:, 61) = rate(:, 61) * m(:)
      rate(:, 62) = rate(:, 62) * m(:)
      rate(:, 63) = rate(:, 63) * m(:)
      rate(:, 64) = rate(:, 64) * m(:)
      rate(:, 67) = rate(:, 67) * m(:)
      rate(:, 68) = rate(:, 68) * m(:)
      rate(:, 69) = rate(:, 69) * m(:)
      rate(:, 70) = rate(:, 70) * m(:)
      rate(:, 71) = rate(:, 71) * m(:)
      rate(:, 72) = rate(:, 72) * m(:)
      rate(:, 73) = rate(:, 73) * m(:)
      rate(:, 74) = rate(:, 74) * m(:)
      rate(:, 75) = rate(:, 75) * m(:)
      rate(:, 77) = rate(:, 77) * m(:)
      rate(:, 78) = rate(:, 78) * m(:)
      rate(:, 79) = rate(:, 79) * m(:)
      rate(:, 80) = rate(:, 80) * m(:)
      rate(:, 81) = rate(:, 81) * m(:)
      rate(:, 82) = rate(:, 82) * m(:)
      rate(:, 84) = rate(:, 84) * m(:)
      rate(:, 85) = rate(:, 85) * m(:)
      rate(:, 86) = rate(:, 86) * m(:)
      rate(:, 87) = rate(:, 87) * m(:)
      rate(:, 88) = rate(:, 88) * m(:)
      rate(:, 89) = rate(:, 89) * m(:)
      rate(:, 90) = rate(:, 90) * m(:)
      rate(:, 91) = rate(:, 91) * m(:)
      rate(:, 92) = rate(:, 92) * m(:)
      rate(:, 93) = rate(:, 93) * m(:)
      rate(:, 94) = rate(:, 94) * m(:)
      rate(:, 97) = rate(:, 97) * m(:)
      rate(:, 98) = rate(:, 98) * m(:)
      rate(:, 99) = rate(:, 99) * m(:)
      rate(:,100) = rate(:,100) * m(:)
      rate(:,101) = rate(:,101) * m(:)
      rate(:,102) = rate(:,102) * m(:)
      rate(:,103) = rate(:,103) * m(:)
      rate(:,104) = rate(:,104) * m(:)
      rate(:,105) = rate(:,105) * m(:)
      rate(:,106) = rate(:,106) * m(:)
      rate(:,107) = rate(:,107) * m(:)
      rate(:,108) = rate(:,108) * m(:)
      rate(:,109) = rate(:,109) * m(:)
      rate(:,110) = rate(:,110) * m(:)
      rate(:,112) = rate(:,112) * m(:)
      rate(:,113) = rate(:,113) * m(:)
      rate(:,114) = rate(:,114) * m(:)
      rate(:,115) = rate(:,115) * m(:)
      rate(:,116) = rate(:,116) * m(:)
      rate(:,117) = rate(:,117) * m(:)
      rate(:,118) = rate(:,118) * m(:)
      rate(:,119) = rate(:,119) * m(:)
      rate(:,120) = rate(:,120) * m(:)
      rate(:,121) = rate(:,121) * m(:)
      rate(:,122) = rate(:,122) * m(:)
      rate(:,123) = rate(:,123) * m(:)
      rate(:,124) = rate(:,124) * m(:)
      rate(:,125) = rate(:,125) * m(:)
      rate(:,126) = rate(:,126) * m(:)
      rate(:,127) = rate(:,127) * m(:)
      rate(:,128) = rate(:,128) * m(:)
      rate(:,130) = rate(:,130) * m(:)
      rate(:,131) = rate(:,131) * m(:)
      rate(:,132) = rate(:,132) * m(:)
      rate(:,133) = rate(:,133) * m(:)
      rate(:,134) = rate(:,134) * m(:)
      rate(:,135) = rate(:,135) * m(:)
      rate(:,136) = rate(:,136) * m(:)
      rate(:,137) = rate(:,137) * m(:)
      rate(:,138) = rate(:,138) * m(:)
      rate(:,139) = rate(:,139) * m(:)
      rate(:,140) = rate(:,140) * m(:)
      rate(:,141) = rate(:,141) * m(:)
      rate(:,142) = rate(:,142) * m(:)
      rate(:,143) = rate(:,143) * m(:)
      rate(:,144) = rate(:,144) * m(:)
      rate(:,145) = rate(:,145) * m(:)
      rate(:,146) = rate(:,146) * m(:)
      rate(:,147) = rate(:,147) * m(:)
      rate(:,148) = rate(:,148) * m(:)
      rate(:,149) = rate(:,149) * m(:)
      rate(:,150) = rate(:,150) * m(:)
      rate(:,151) = rate(:,151) * m(:)
      rate(:,153) = rate(:,153) * m(:)
      rate(:,154) = rate(:,154) * m(:)
      rate(:,155) = rate(:,155) * m(:)
      rate(:,156) = rate(:,156) * m(:)
      rate(:,157) = rate(:,157) * m(:)
      rate(:,158) = rate(:,158) * m(:)
      rate(:,159) = rate(:,159) * m(:)
      rate(:,160) = rate(:,160) * m(:)
      rate(:,161) = rate(:,161) * m(:)
      rate(:,162) = rate(:,162) * m(:)
      rate(:,163) = rate(:,163) * m(:)
      rate(:,164) = rate(:,164) * m(:)
      rate(:,165) = rate(:,165) * m(:)
      rate(:,166) = rate(:,166) * m(:)
      rate(:,167) = rate(:,167) * m(:)
      rate(:,169) = rate(:,169) * m(:)
      rate(:,170) = rate(:,170) * m(:)
      rate(:,171) = rate(:,171) * m(:)
      rate(:,172) = rate(:,172) * m(:)
      rate(:,173) = rate(:,173) * m(:)
      rate(:,174) = rate(:,174) * m(:)
      rate(:,175) = rate(:,175) * m(:)
      rate(:,176) = rate(:,176) * m(:)
      rate(:,177) = rate(:,177) * m(:)
      rate(:,178) = rate(:,178) * m(:)
      rate(:,179) = rate(:,179) * m(:)
      rate(:,180) = rate(:,180) * m(:)
      rate(:,181) = rate(:,181) * m(:)
      rate(:,182) = rate(:,182) * m(:)
      rate(:,183) = rate(:,183) * m(:)
      rate(:,184) = rate(:,184) * m(:)
      rate(:,185) = rate(:,185) * m(:)
      rate(:,186) = rate(:,186) * m(:)
      rate(:,187) = rate(:,187) * m(:)
      rate(:,188) = rate(:,188) * m(:)
      rate(:,189) = rate(:,189) * m(:)
      rate(:,190) = rate(:,190) * m(:)
      rate(:,191) = rate(:,191) * m(:)
      rate(:,192) = rate(:,192) * m(:)
      rate(:,193) = rate(:,193) * m(:)
      rate(:,194) = rate(:,194) * m(:)
      rate(:,195) = rate(:,195) * m(:)
      rate(:,196) = rate(:,196) * m(:)
      rate(:,197) = rate(:,197) * m(:)
      rate(:,198) = rate(:,198) * m(:)
      rate(:,199) = rate(:,199) * m(:)
      rate(:,200) = rate(:,200) * m(:)
      rate(:,201) = rate(:,201) * m(:)
      rate(:,202) = rate(:,202) * m(:)
      rate(:,203) = rate(:,203) * m(:)
      rate(:,204) = rate(:,204) * m(:)
      rate(:,205) = rate(:,205) * m(:)
      rate(:,206) = rate(:,206) * m(:)
      rate(:,207) = rate(:,207) * m(:)
      rate(:,208) = rate(:,208) * m(:)
      rate(:,209) = rate(:,209) * m(:)
      rate(:,210) = rate(:,210) * m(:)
      rate(:,211) = rate(:,211) * m(:)
      rate(:,212) = rate(:,212) * m(:)
      rate(:,213) = rate(:,213) * m(:)
      rate(:,215) = rate(:,215) * m(:)
      rate(:,216) = rate(:,216) * m(:)
      rate(:,217) = rate(:,217) * m(:)
      rate(:,218) = rate(:,218) * m(:)
      rate(:,221) = rate(:,221) * m(:)
      rate(:,223) = rate(:,223) * m(:)
      rate(:,230) = rate(:,230) * m(:)
      rate(:,231) = rate(:,231) * m(:)
      rate(:,232) = rate(:,232) * m(:)
      rate(:,233) = rate(:,233) * m(:)
      rate(:,234) = rate(:,234) * m(:)
      rate(:,235) = rate(:,235) * m(:)
      rate(:,236) = rate(:,236) * m(:)
      rate(:,237) = rate(:,237) * m(:)
      rate(:,238) = rate(:,238) * m(:)
      rate(:,241) = rate(:,241) * m(:)
      rate(:,263) = rate(:,263) * m(:)
      rate(:,264) = rate(:,264) * m(:)
      rate(:,266) = rate(:,266) * m(:)
      rate(:,267) = rate(:,267) * m(:)
      rate(:,268) = rate(:,268) * m(:)
      rate(:,269) = rate(:,269) * m(:)
      rate(:,270) = rate(:,270) * m(:)
      rate(:,271) = rate(:,271) * m(:)
      rate(:,272) = rate(:,272) * m(:)
      rate(:,273) = rate(:,273) * m(:)
      rate(:,274) = rate(:,274) * m(:)
      rate(:,275) = rate(:,275) * m(:)
      rate(:,276) = rate(:,276) * m(:)
      rate(:,277) = rate(:,277) * m(:)
      rate(:,278) = rate(:,278) * m(:)
      rate(:,280) = rate(:,280) * m(:)
      rate(:,281) = rate(:,281) * m(:)
      rate(:,282) = rate(:,282) * m(:)
      rate(:,283) = rate(:,283) * m(:)
      rate(:,284) = rate(:,284) * m(:)
      rate(:,285) = rate(:,285) * m(:)
      rate(:,286) = rate(:,286) * m(:)
      rate(:,287) = rate(:,287) * m(:)
      rate(:,288) = rate(:,288) * m(:)
      rate(:,289) = rate(:,289) * m(:)
      rate(:,290) = rate(:,290) * m(:)
      rate(:,291) = rate(:,291) * m(:)
      rate(:,292) = rate(:,292) * m(:)
      rate(:,293) = rate(:,293) * m(:)
      rate(:,294) = rate(:,294) * m(:)
      rate(:,295) = rate(:,295) * m(:)
      rate(:,296) = rate(:,296) * m(:)
      rate(:,297) = rate(:,297) * m(:)
      rate(:,298) = rate(:,298) * m(:)
      rate(:,299) = rate(:,299) * m(:)
      rate(:,300) = rate(:,300) * m(:)
      rate(:,301) = rate(:,301) * m(:)
      rate(:,302) = rate(:,302) * m(:)
      rate(:,303) = rate(:,303) * m(:)
      end subroutine adjrxt
      end module mo_adjrxt_mod
      module mo_phtadj_mod
      private
      public :: phtadj
      contains
      subroutine phtadj( p_rate, inv, m, plnplv )
      use CHEM_MODS_MOD, only : nfs, phtcnt
      implicit none
!--------------------------------------------------------------------
! ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in) :: inv(plnplv,nfs)
      real, intent(in) :: m(plnplv)
      real, intent(inout) :: p_rate(plnplv,phtcnt)
!--------------------------------------------------------------------
! ... Local variables
!--------------------------------------------------------------------
      real :: im(plnplv)
      im(:) = 1. / m(:)
      p_rate(:, 1) = p_rate(:, 1) * inv(:, 3) * im(:)
      end subroutine phtadj
      end module mo_phtadj_mod
      module mo_rxt_mod
      private
      public :: rxt_mod
      contains
      subroutine rxt_mod( rate, het_rates, grp_ratios, plnplv )
      use CHEM_MODS_MOD, only : rxntot, hetcnt, grpcnt
      implicit none
!---------------------------------------------------------------------------
! ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(inout) :: rate(plnplv,rxntot)
      real, intent(inout) :: het_rates(plnplv,hetcnt)
      real, intent(in) :: grp_ratios(plnplv,grpcnt)
      end subroutine rxt_mod
      end module mo_rxt_mod
      module mo_make_grp_vmr_mod
      private
      public :: mak_grp_vmr
      contains
      subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )
      use MO_GRID_MOD, only : plev, pcnstm1
      use CHEM_MODS_MOD, only : grpcnt
      implicit none
!----------------------------------------------------------------------------
! ... Dummy arguments
!----------------------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in) :: vmr(plonl,plev,pcnstm1)
      real, intent(in) :: group_ratios(plonl,plev,grpcnt)
      real, intent(out) :: group_vmrs(plonl,plev,grpcnt)
!----------------------------------------------------------------------------
! ... Local variables
!----------------------------------------------------------------------------
      integer :: k
      end subroutine mak_grp_vmr
      end module mo_make_grp_vmr_mod
