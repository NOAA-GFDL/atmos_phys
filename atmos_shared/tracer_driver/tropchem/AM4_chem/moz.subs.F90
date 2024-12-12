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
      rate(:,:,52) = 1.2e-10
      rate(:,:,58) = 1.8e-12
      rate(:,:,60) = 1.8e-12
      rate(:,:,71) = 3.5e-12
      rate(:,:,74) = 0.
      rate(:,:,82) = 1.5e-10
      rate(:,:,92) = 1.e-14
      rate(:,:,113) = 3.e-13
      rate(:,:,114) = 4.1e-14
      rate(:,:,130) = 8.37E-14
      rate(:,:,131) = 1.54E-13
      rate(:,:,135) = 3.7E-19
      rate(:,:,145) = 8.37E-14
      rate(:,:,147) = 1.6E-12
      rate(:,:,151) = 3.4E-15
      rate(:,:,154) = 8.37E-14
      rate(:,:,157) = 3.20E-12
      rate(:,:,169) = 8.37E-14
      rate(:,:,176) = 2.90E-11
      rate(:,:,177) = 1.00E-11
      rate(:,:,179) = 4.0E-16
      rate(:,:,180) = 1.50E-11
      rate(:,:,185) = 2.30E-12
      rate(:,:,190) = 4.00E-12
      rate(:,:,195) = 1.6e-12
      rate(:,:,270) = 3.17e-8
      itemp(:,:) = 1. / temp(:,:)
      rate(:,:,47) = 8e-12 * exp( -2060. * itemp(:,:) )
      rate(:,:,48) = 1.5e-11 * exp( -3600. * itemp(:,:) )
      rate(:,:,49) = 2.1e-11 * exp( 100. * itemp(:,:) )
      exp_fac(:,:) = exp( 180. * itemp(:,:) )
      rate(:,:,53) = 1.8e-11 * exp_fac(:,:)
      rate(:,:,91) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,118) = 4.2e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 200. * itemp(:,:) )
      rate(:,:,54) = 3e-11 * exp_fac(:,:)
      rate(:,:,87) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,100) = 8.78E-12 * exp_fac(:,:)
      rate(:,:,107) = 3.80e-12 * exp_fac(:,:)
      rate(:,:,115) = 5.18e-12 * exp_fac(:,:)
      rate(:,:,121) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,136) = 4.75E-12 * exp_fac(:,:)
      rate(:,:,146) = 8.78E-12 * exp_fac(:,:)
      rate(:,:,155) = 1.84E-12 * exp_fac(:,:)
      rate(:,:,166) = 6.13E-13 * exp_fac(:,:)
      rate(:,:,174) = 2.66e-12 * exp_fac(:,:)
      rate(:,:,175) = 1.14e-12 * exp_fac(:,:)
      rate(:,:,188) = 5.18E-12 * exp_fac(:,:)
      rate(:,:,254) = 5.5e-12 * exp_fac(:,:)
      rate(:,:,55) = 1.7e-12 * exp( -940. * itemp(:,:) )
      rate(:,:,56) = 1e-14 * exp( -490. * itemp(:,:) )
      rate(:,:,59) = 4.8e-11 * exp( 250. * itemp(:,:) )
      rate(:,:,61) = 2.8e-12 * exp( -1800. * itemp(:,:) )
      rate(:,:,62) = 2.15e-11 * exp( 110. * itemp(:,:) )
      rate(:,:,63) = 3.3e-11 * exp( 55. * itemp(:,:) )
      rate(:,:,64) = 1.63e-10 * exp( 60. * itemp(:,:) )
      exp_fac(:,:) = exp( 20. * itemp(:,:) )
      rate(:,:,65) = 7.25e-11 * exp_fac(:,:)
      rate(:,:,66) = 4.63e-11 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270. * itemp(:,:) )
      rate(:,:,67) = 3.3e-12 * exp_fac(:,:)
      rate(:,:,103) = 8.1e-12 * exp_fac(:,:)
      rate(:,:,241) = 7.4e-12 * exp_fac(:,:)
      rate(:,:,68) = 3e-12 * exp( -1500. * itemp(:,:) )
      rate(:,:,69) = 5.1e-12 * exp( 210. * itemp(:,:) )
      exp_fac(:,:) = exp( -2450. * itemp(:,:) )
      rate(:,:,70) = 1.2e-13 * exp_fac(:,:)
      rate(:,:,258) = 8.5e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170. * itemp(:,:) )
      rate(:,:,77) = 1.5e-11 * exp_fac(:,:)
      rate(:,:,239) = 1.8e-11 * exp_fac(:,:)
      exp_fac(:,:) = exp( 380. * itemp(:,:) )
      rate(:,:,79) = 1.3e-12 * exp_fac(:,:)
      rate(:,:,132) = 3.61E-12 * exp_fac(:,:)
      rate(:,:,148) = 8.00E-12 * exp_fac(:,:)
      rate(:,:,156) = 4.40E-12 * exp_fac(:,:)
      rate(:,:,167) = 3.60E-12 * exp_fac(:,:)
      rate(:,:,81) = 2.45e-12 * exp( -1775. * itemp(:,:) )
      exp_fac(:,:) = exp( 300. * itemp(:,:) )
      rate(:,:,83) = 2.8e-12 * exp_fac(:,:)
      rate(:,:,172) = 2.80E-12 * exp_fac(:,:)
      rate(:,:,84) = 6.03e-13 * exp( -453. * itemp(:,:) )
      rate(:,:,85) = 2.30e-14 * exp( 677. * itemp(:,:) )
      rate(:,:,86) = 4.1e-13 * exp( 750. * itemp(:,:) )
      exp_fac(:,:) = exp( -1900. * itemp(:,:) )
      rate(:,:,88) = 3.4e-13 * exp_fac(:,:)
      rate(:,:,102) = 1.4e-12 * exp_fac(:,:)
      rate(:,:,89) = 5.5e-12 * exp( 125. * itemp(:,:) )
      rate(:,:,93) = 1.6e11 * exp( -4150. * itemp(:,:) )
      rate(:,:,94) = 1.2e-14 * exp( -2630. * itemp(:,:) )
      rate(:,:,96) = 5.50E-15 * exp( -1880. * itemp(:,:) )
      rate(:,:,97) = 4.6e-13 * exp( -1156. * itemp(:,:) )
      exp_fac(:,:) = exp( 350. * itemp(:,:) )
      rate(:,:,98) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,101) = 4.63E-12 * exp_fac(:,:)
      rate(:,:,125) = 3.10E-11 * exp_fac(:,:)
      rate(:,:,128) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,140) = 2.70E-12 * exp_fac(:,:)
      rate(:,:,143) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,152) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,158) = 2.7E-12 * exp_fac(:,:)
      rate(:,:,170) = 2.35E-12 * exp_fac(:,:)
      rate(:,:,171) = 0.35E-12 * exp_fac(:,:)
      rate(:,:,184) = 2.70E-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 700. * itemp(:,:) )
      rate(:,:,99) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,112) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,119) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,173) = 8.60E-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 980. * itemp(:,:) )
      rate(:,:,105) = 5.2e-13 * exp_fac(:,:)
      rate(:,:,192) = 5.20E-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500. * itemp(:,:) )
      rate(:,:,106) = 2.0e-12 * exp_fac(:,:)
      rate(:,:,109) = 2.5e-12 * exp_fac(:,:)
      rate(:,:,162) = 1.68E-12 * exp_fac(:,:)
      rate(:,:,163) = 1.87E-13 * exp_fac(:,:)
      rate(:,:,110) = 7.66e-12 * exp( -1020. * itemp(:,:) )
      rate(:,:,111) = 2.6e-12 * exp( 365. * itemp(:,:) )
      rate(:,:,116) = 1.55e-11 * exp( -540. * itemp(:,:) )
      rate(:,:,117) = 8.7e-12 * exp( -615. * itemp(:,:) )
      rate(:,:,120) = 3.75e-13 * exp( -40. * itemp(:,:) )
      rate(:,:,123) = 2.9e-12 * exp( -345. * itemp(:,:) )
      rate(:,:,124) = 6.9e-12 * exp( -230. * itemp(:,:) )
      rate(:,:,126) = 4.07E+08 * exp( -7694. * itemp(:,:) )
      rate(:,:,127) = 1.00e-14 * exp( -1970. * itemp(:,:) )
      exp_fac(:,:) = exp( 1300. * itemp(:,:) )
      rate(:,:,129) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,139) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,144) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,153) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,159) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,168) = 1.82E-13 * exp_fac(:,:)
      rate(:,:,186) = 2.06E-13 * exp_fac(:,:)
      rate(:,:,133) = 2.4e-12 * exp( 360. * itemp(:,:) )
      rate(:,:,134) = 8.7e-14 * exp( 1650. * itemp(:,:) )
      exp_fac(:,:) = exp( 390. * itemp(:,:) )
      rate(:,:,137) = 1.90E-11 * exp_fac(:,:)
      rate(:,:,187) = 1.90E-11 * exp_fac(:,:)
      rate(:,:,138) = 5.78E-11 * exp( -400. * itemp(:,:) )
      rate(:,:,141) = 2.6E-12 * exp( 610. * itemp(:,:) )
      exp_fac(:,:) = exp( -1520. * itemp(:,:) )
      rate(:,:,142) = 8.5E-16 * exp_fac(:,:)
      rate(:,:,193) = 4.15E-15 * exp_fac(:,:)
      rate(:,:,149) = 2.90E+07 * exp( -5297. * itemp(:,:) )
      rate(:,:,150) = 1.40E-15 * exp( -2100. * itemp(:,:) )
      exp_fac(:,:) = exp( 340. * itemp(:,:) )
      rate(:,:,160) = 6.7E-12 * exp_fac(:,:)
      rate(:,:,178) = 3.1E-12 * exp_fac(:,:)
      rate(:,:,191) = 6.70E-12 * exp_fac(:,:)
      rate(:,:,161) = 4.3E-13 * exp( 1040. * itemp(:,:) )
      rate(:,:,181) = 1.40E-12 * exp( -1860. * itemp(:,:) )
      rate(:,:,182) = 1.60E-12 * exp( 305. * itemp(:,:) )
      rate(:,:,183) = 3.30E-12 * exp( -450. * itemp(:,:) )
      rate(:,:,189) = 3.15E-13 * exp( -448. * itemp(:,:) )
      rate(:,:,194) = 7.48E-12 * exp( 410. * itemp(:,:) )
      rate(:,:,196) = 1.2e-11 * exp( 440. * itemp(:,:) )
      rate(:,:,197) = 5.3e-16 * exp( -530. * itemp(:,:) )
      rate(:,:,198) = 1.2e-12 * exp( 490. * itemp(:,:) )
      rate(:,:,205) = 1.2e-11 * exp( -280. * itemp(:,:) )
      rate(:,:,207) = 1.90e-13 * exp( 530. * itemp(:,:) )
      rate(:,:,229) = 1.7e-12 * exp( -710. * itemp(:,:) )
      rate(:,:,230) = 1.4e-10 * exp( -470. * itemp(:,:) )
      rate(:,:,232) = 2.3e-11 * exp( -200. * itemp(:,:) )
      rate(:,:,233) = 2.8e-11 * exp( 85. * itemp(:,:) )
      exp_fac(:,:) = exp( 290. * itemp(:,:) )
      rate(:,:,234) = 6.4e-12 * exp_fac(:,:)
      rate(:,:,255) = 4.1e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800. * itemp(:,:) )
      rate(:,:,236) = 2.9e-12 * exp_fac(:,:)
      rate(:,:,246) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,253) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,237) = 7.3e-12 * exp( -1280. * itemp(:,:) )
      rate(:,:,238) = 2.6e-12 * exp( -350. * itemp(:,:) )
      exp_fac(:,:) = exp( 220. * itemp(:,:) )
      rate(:,:,240) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,260) = 5.8e-12 * exp_fac(:,:)
      rate(:,:,242) = 8.1e-11 * exp( -30. * itemp(:,:) )
      exp_fac(:,:) = exp( 260. * itemp(:,:) )
      rate(:,:,248) = 2.3e-12 * exp_fac(:,:)
      rate(:,:,250) = 8.8e-12 * exp_fac(:,:)
      rate(:,:,249) = 4.5e-12 * exp( 460. * itemp(:,:) )
      rate(:,:,251) = 1.2e-10 * exp( -430. * itemp(:,:) )
      rate(:,:,252) = 4.8e-12 * exp( -310. * itemp(:,:) )
      rate(:,:,256) = 6.0e-13 * exp( 230. * itemp(:,:) )
      rate(:,:,257) = 4.5e-14 * exp( -1260. * itemp(:,:) )
      itemp(:,:) = 300. * itemp(:,:)
      ko(:,:) = 5.9e-33 * itemp(:,:)**1.4
      kinf(:,:) = 1.1e-12 * itemp(:,:)**-1.3
      call jpl( rate(1,1,50), m, .6, ko, kinf, plnplv )
      ko(:,:) = 1.5e-13 * itemp(:,:)**-0.6
      kinf(:,:) = 2.1e9 * itemp(:,:)**-6.1
      call jpl( rate(1,1,51), m, .6, ko, kinf, plnplv )
      ko(:,:) = 2.e-30 * itemp(:,:)**4.4
      kinf(:,:) = 1.4e-12 * itemp(:,:)**.7
      call jpl( rate(1,1,72), m, .6, ko, kinf, plnplv )
      ko(:,:) = 1.8e-30 * itemp(:,:)**3.0
      kinf(:,:) = 2.8e-11
      call jpl( rate(1,1,75), m, .6, ko, kinf, plnplv )
      ko(:,:) = 2.0e-31 * itemp(:,:)**3.4
      kinf(:,:) = 2.9e-12 * itemp(:,:)**1.1
      call jpl( rate(1,1,78), m, .6, ko, kinf, plnplv )
      ko(:,:) = 1.e-28 * itemp(:,:)**4.5
      kinf(:,:) = 7.5e-12 * itemp(:,:)**0.85
      call jpl( rate(1,1,90), m, .6, ko, kinf, plnplv )
      ko(:,:) = 8.e-27 * itemp(:,:)**3.5
      kinf(:,:) = 3.e-11
      call jpl( rate(1,1,95), m, .5, ko, kinf, plnplv )
      ko(:,:) = 9.7e-29 * itemp(:,:)**5.6
      kinf(:,:) = 9.3e-12 * itemp(:,:)**1.5
      call jpl( rate(1,1,104), m, .6, ko, kinf, plnplv )
      ko(:,:) = 9.0e-28 * itemp(:,:)**8.9
      kinf(:,:) = 7.7e-12 * itemp(:,:)**.2
      call jpl( rate(1,1,164), m, .6, ko, kinf, plnplv )
      ko(:,:) = 3.3e-31 * itemp(:,:)**4.3
      kinf(:,:) = 1.6e-12
      call jpl( rate(1,1,204), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 4.4e-32 * itemp(:,:)**1.3
      kinf(:,:) = 4.7e-11 * itemp(:,:)**0.2
      call jpl( rate(1,1,231), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 1.8e-31 * itemp(:,:)**3.4
      kinf(:,:) = 1.5e-11 * itemp(:,:)**1.9
      call jpl( rate(1,1,235), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 6.9e-31 * itemp(:,:)**1.0
      kinf(:,:) = 2.6e-11
      call jpl( rate(1,1,243), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 1.6e-32 * itemp(:,:)**4.5
      kinf(:,:) = 2.0e-12 * itemp(:,:)**2.4
      call jpl( rate(1,1,244), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 5.2e-31 * itemp(:,:)**3.2
      kinf(:,:) = 6.9e-12 * itemp(:,:)**2.9
      call jpl( rate(1,1,247), m, 0.6, ko, kinf, plnplv )
      ko(:,:) = 9.0e-32 * itemp(:,:)**1.5
      kinf(:,:) = 3.0e-11
      call jpl( rate(1,1,259), m, 0.6, ko, kinf, plnplv )
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
      rate(:, 48) = rate(:, 48) * inv(:, 3)
      rate(:, 50) = rate(:, 50) * inv(:, 1)
      rate(:, 62) = rate(:, 62) * inv(:, 2)
      rate(:, 63) = rate(:, 63) * inv(:, 3)
      rate(:, 72) = rate(:, 72) * inv(:, 1)
      rate(:, 73) = rate(:, 73) * inv(:, 1)
      rate(:, 75) = rate(:, 75) * inv(:, 1)
      rate(:, 78) = rate(:, 78) * inv(:, 1)
      rate(:, 80) = rate(:, 80) * inv(:, 1)
      rate(:, 90) = rate(:, 90) * inv(:, 1)
      rate(:, 92) = rate(:, 92) * inv(:, 3)
      rate(:, 95) = rate(:, 95) * inv(:, 1)
      rate(:,104) = rate(:,104) * inv(:, 1)
      rate(:,108) = rate(:,108) * inv(:, 1)
      rate(:,164) = rate(:,164) * inv(:, 1)
      rate(:,204) = rate(:,204) * inv(:, 1)
      rate(:,235) = rate(:,235) * inv(:, 1)
      rate(:,243) = rate(:,243) * inv(:, 1)
      rate(:,244) = rate(:,244) * inv(:, 1)
      rate(:,245) = rate(:,245) * inv(:, 1)
      rate(:,247) = rate(:,247) * inv(:, 1)
      rate(:,259) = rate(:,259) * inv(:, 1)
      rate(:, 46) = rate(:, 46) * inv(:, 3) * inv(:, 1)
      rate(:,231) = rate(:,231) * inv(:, 3) * inv(:, 1)
      rate(:, 47) = rate(:, 47) * m(:)
      rate(:, 49) = rate(:, 49) * m(:)
      rate(:, 50) = rate(:, 50) * m(:)
      rate(:, 51) = rate(:, 51) * m(:)
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
      rate(:, 64) = rate(:, 64) * m(:)
      rate(:, 65) = rate(:, 65) * m(:)
      rate(:, 66) = rate(:, 66) * m(:)
      rate(:, 67) = rate(:, 67) * m(:)
      rate(:, 68) = rate(:, 68) * m(:)
      rate(:, 69) = rate(:, 69) * m(:)
      rate(:, 70) = rate(:, 70) * m(:)
      rate(:, 71) = rate(:, 71) * m(:)
      rate(:, 72) = rate(:, 72) * m(:)
      rate(:, 74) = rate(:, 74) * m(:)
      rate(:, 75) = rate(:, 75) * m(:)
      rate(:, 76) = rate(:, 76) * m(:)
      rate(:, 77) = rate(:, 77) * m(:)
      rate(:, 78) = rate(:, 78) * m(:)
      rate(:, 79) = rate(:, 79) * m(:)
      rate(:, 81) = rate(:, 81) * m(:)
      rate(:, 82) = rate(:, 82) * m(:)
      rate(:, 83) = rate(:, 83) * m(:)
      rate(:, 84) = rate(:, 84) * m(:)
      rate(:, 85) = rate(:, 85) * m(:)
      rate(:, 86) = rate(:, 86) * m(:)
      rate(:, 87) = rate(:, 87) * m(:)
      rate(:, 88) = rate(:, 88) * m(:)
      rate(:, 89) = rate(:, 89) * m(:)
      rate(:, 90) = rate(:, 90) * m(:)
      rate(:, 91) = rate(:, 91) * m(:)
      rate(:, 94) = rate(:, 94) * m(:)
      rate(:, 95) = rate(:, 95) * m(:)
      rate(:, 96) = rate(:, 96) * m(:)
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
      rate(:,109) = rate(:,109) * m(:)
      rate(:,110) = rate(:,110) * m(:)
      rate(:,111) = rate(:,111) * m(:)
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
      rate(:,127) = rate(:,127) * m(:)
      rate(:,128) = rate(:,128) * m(:)
      rate(:,129) = rate(:,129) * m(:)
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
      rate(:,150) = rate(:,150) * m(:)
      rate(:,151) = rate(:,151) * m(:)
      rate(:,152) = rate(:,152) * m(:)
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
      rate(:,166) = rate(:,166) * m(:)
      rate(:,167) = rate(:,167) * m(:)
      rate(:,168) = rate(:,168) * m(:)
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
      rate(:,204) = rate(:,204) * m(:)
      rate(:,205) = rate(:,205) * m(:)
      rate(:,206) = rate(:,206) * m(:)
      rate(:,207) = rate(:,207) * m(:)
      rate(:,229) = rate(:,229) * m(:)
      rate(:,230) = rate(:,230) * m(:)
      rate(:,232) = rate(:,232) * m(:)
      rate(:,233) = rate(:,233) * m(:)
      rate(:,234) = rate(:,234) * m(:)
      rate(:,235) = rate(:,235) * m(:)
      rate(:,236) = rate(:,236) * m(:)
      rate(:,237) = rate(:,237) * m(:)
      rate(:,238) = rate(:,238) * m(:)
      rate(:,239) = rate(:,239) * m(:)
      rate(:,240) = rate(:,240) * m(:)
      rate(:,241) = rate(:,241) * m(:)
      rate(:,242) = rate(:,242) * m(:)
      rate(:,243) = rate(:,243) * m(:)
      rate(:,244) = rate(:,244) * m(:)
      rate(:,246) = rate(:,246) * m(:)
      rate(:,247) = rate(:,247) * m(:)
      rate(:,248) = rate(:,248) * m(:)
      rate(:,249) = rate(:,249) * m(:)
      rate(:,250) = rate(:,250) * m(:)
      rate(:,251) = rate(:,251) * m(:)
      rate(:,252) = rate(:,252) * m(:)
      rate(:,253) = rate(:,253) * m(:)
      rate(:,254) = rate(:,254) * m(:)
      rate(:,255) = rate(:,255) * m(:)
      rate(:,256) = rate(:,256) * m(:)
      rate(:,257) = rate(:,257) * m(:)
      rate(:,258) = rate(:,258) * m(:)
      rate(:,259) = rate(:,259) * m(:)
      rate(:,260) = rate(:,260) * m(:)
      rate(:,261) = rate(:,261) * m(:)
      rate(:,262) = rate(:,262) * m(:)
      rate(:,263) = rate(:,263) * m(:)
      rate(:,264) = rate(:,264) * m(:)
      rate(:,265) = rate(:,265) * m(:)
      rate(:,266) = rate(:,266) * m(:)
      rate(:,267) = rate(:,267) * m(:)
      rate(:,268) = rate(:,268) * m(:)
      rate(:,269) = rate(:,269) * m(:)
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
