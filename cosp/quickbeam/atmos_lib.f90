
!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------
        
! $Id: atmos_lib.f90,v 1.1.2.1.2.1 2009/08/10 10:48:14 rsh Exp $
! $Name: riga_201012 $

! ATMOS_LIB: Atmospheric science procedures for F90
! Compiled/Modified:
!   07/01/06  John Haynes (haynes@atmos.colostate.edu)
!
! mcclatchey (subroutine)
  
  module atmos_lib
  implicit none
  
  contains
  
! ----------------------------------------------------------------------------
! subroutine MCCLATCHEY
! ----------------------------------------------------------------------------
  subroutine mcclatchey(stype,hgt,prs,tk,rh)
  implicit none
!
! Purpose:
!   returns a standard atmospheric profile
!
! Input:
!   [stype]   type of profile to return
!             1 = mid-latitude summer
!             2 = mid-latitude winter
!             3 = tropical
!
! Outputs:
!   [hgt]     height (m)
!   [prs]     pressure (hPa)
!   [tk]      temperature (K)
!   [rh]      relative humidity (%)
!
! Created:
!   06/01/2006  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----
  integer, intent(in) :: &
  stype

  integer, parameter :: ndat = 33

! ----- OUTPUTS -----
  real*8, intent(out), dimension(ndat) :: &
  hgt, &                        ! height (m)
  prs, &                        ! pressure (hPa)
  tk, &                         ! temperature (K)
  rh                            ! relative humidity (%)
  
  hgt = (/0.00000,1000.00,2000.00,3000.00,4000.00,5000.00, &
          6000.00,7000.00,8000.00,9000.00,10000.0,11000.0, &
          12000.0,13000.0,14000.0,15000.0,16000.0,17000.0, &
          18000.0,19000.0,20000.0,21000.0,22000.0,23000.0, &
          24000.0,25000.0,30000.0,35000.0,40000.0,45000.0, &
          50000.0,70000.0,100000./)

  select case(stype)

  case(1)
!   // mid-latitide summer  
    prs = (/1013.00, 902.000, 802.000, 710.000, 628.000, 554.000, &
            487.000, 426.000, 372.000, 324.000, 281.000, 243.000, &
            209.000, 179.000, 153.000, 130.000, 111.000, 95.0000, &
            81.2000, 69.5000, 59.5000, 51.0000, 43.7000, 37.6000, &
            32.2000, 27.7000, 13.2000, 6.52000, 3.33000, 1.76000, &
            0.951000,0.0671000,0.000300000/)
	   
    tk =  (/294.000, 290.000, 285.000, 279.000, 273.000, 267.000, &
            261.000, 255.000, 248.000, 242.000, 235.000, 229.000, &
            222.000, 216.000, 216.000, 216.000, 216.000, 216.000, &
            216.000, 217.000, 218.000, 219.000, 220.000, 222.000, &
            223.000, 224.000, 234.000, 245.000, 258.000, 270.000, &
            276.000, 218.000, 210.000/)

    rh =  (/74.8384, 63.4602, 55.0485, 45.4953, 39.3805, 31.7965, &
            30.3958, 29.5966, 30.1626, 29.3624, 30.3334, 19.0768, &
            11.0450, 6.61278, 3.67379, 2.79209, 2.35123, 2.05732, &
            1.83690, 1.59930, 1.30655, 1.31890, 1.17620,0.994076, &
            0.988566,0.989143,0.188288,0.0205613,0.00271164,0.000488798, &
            0.000107066,0.000406489,7.68645e-06/)

  case(2)
!   // mid-latitude winter
    prs = (/1018.00, 897.300, 789.700, 693.800, 608.100, 531.300, &
            462.700, 401.600, 347.300, 299.200, 256.800, 219.900, &
            188.200, 161.000, 137.800, 117.800, 100.700, 86.1000, &
            73.5000, 62.8000, 53.7000, 45.8000, 39.1000, 33.4000, &
            28.6000, 24.3000, 11.1000, 5.18000, 2.53000, 1.29000, &
            0.682000,0.0467000,0.000300000/)

    tk =  (/272.200, 268.700, 265.200, 261.700, 255.700, 249.700, &
            243.700, 237.700, 231.700, 225.700, 219.700, 219.200, &
            218.700, 218.200, 217.700, 217.200, 216.700, 216.200, &
            215.700, 215.200, 215.200, 215.200, 215.200, 215.200, &
            215.200, 215.200, 217.400, 227.800, 243.200, 258.500, &
            265.700, 230.700, 210.200/)

    rh =  (/76.6175, 70.1686, 65.2478, 56.6267, 49.8755, 47.1765, &
            44.0477, 31.0565, 23.0244, 19.6510, 17.8987, 17.4376, &
            16.0621, 5.10608, 3.00679, 2.42293, 2.16406, 2.00901, &
            1.90374, 1.98072, 1.81902, 2.06155, 2.06154, 2.18280, &
            2.42531,2.70824,1.12105,0.108119,0.00944200,0.00115201, &
            0.000221094,0.000101946,7.49350e-06/)

  case(3)
!   // tropical
    prs = (/1013.00, 904.000, 805.000, 715.000, 633.000, 559.000, &
            492.000, 432.000, 378.000, 329.000, 286.000, 247.000, &
            213.000, 182.000, 156.000, 132.000, 111.000, 93.7000, &
            78.9000, 66.6000, 56.5000, 48.0000, 40.9000, 35.0000, &
            30.0000, 25.7000, 12.2000, 6.00000, 3.05000, 1.59000, &
            0.854000,0.0579000,0.000300000/)

    tk =  (/300.000, 294.000, 288.000, 284.000, 277.000, 270.000, &
            264.000, 257.000, 250.000, 244.000, 237.000, 230.000, &
            224.000, 217.000, 210.000, 204.000, 197.000, 195.000, &
            199.000, 203.000, 207.000, 211.000, 215.000, 217.000, &
            219.000, 221.000, 232.000, 243.000, 254.000, 265.000, &
            270.000, 219.000, 210.000/)

    rh =  (/71.4334, 69.4097, 71.4488, 46.7724, 34.7129, 38.3820, &
            33.7214, 32.0122, 30.2607, 24.5059, 19.5321, 13.2966, &
            8.85795, 5.87496, 7.68644, 12.8879, 29.4976, 34.9351, &
            17.1606, 9.53422, 5.10154, 3.45407, 2.11168, 1.76247, &
            1.55162,1.37966,0.229799,0.0245943,0.00373686,0.000702138, &
            0.000162076,0.000362055,7.68645e-06/)
	    
  case default
    print *, 'Must enter a profile type'
    stop
    
  end select
  
  end subroutine mcclatchey
  
  end module atmos_lib
