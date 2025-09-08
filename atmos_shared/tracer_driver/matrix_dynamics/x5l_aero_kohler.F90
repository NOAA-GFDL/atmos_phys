module aero_wet 
use mpp_mod, only: mpp_pe, mpp_sync

        !-----------------------------------------------------------------------
    ! written by x5l (Xiaohan.Li@noaa.gov), used to calculate water uptake by k-Kohler theory
    !--------------
    ! use kappa-kohler theory: Petters & Kreidenweis, 2007
    !        S = (r_p**3 - r_d**3)/(r_p**3 - r_d**3 + kappa*r_d**3) * exp(a/r_p)
    ! here S is relative humidity, r_p is the wet-particle radius, r_d is dry radius
    ! kappa is hygroscopicity, a=2*sigma*M_w/(R*T*rho_w) is the surface tension factor
    ! ------------- 
    !        S = (1-kappa*r_d**3/(r_p**3 - r_d**3 + kappa*r_d**3)) * exp(a/r_p)
    !(1) using taylor expansion: if kappa*r_d**3/(r_p**3 - r_d**3 + kappa*r_d**3) is small
    ! then
    !     p43 = -a/logS, p42 = 0, P41 = (kappa-1)*r_d**3 + kappa/logS*r_d**3, p40 = a*r_d**3*(1-kappa)/logS
    !(2) if S is aournd 1, then 
    !    log(1-kappa*r_d**3/(r_p**3 - r_d**3 + kappa*r_d**3)) + a/r_p = 0
    !by taylor expansion, then
    !     p32 = 0 , p31 = -kappa*r_d**3/a , p30 = (kappa - 1)* r_d**3    
    !
    !-> follow MAM: if RH > 1-epsilon, then interpolate the solution from P3 and P4
    !-----------------------------------------------------------------------
CONTAINS
    subroutine aero_kohler(ddry_in, hygro, s, tair, dwet_out) !note: ddry_in and dwet_out are volume mean diamter, not Dg
        implicit none
        real, intent(in):: ddry_in !dry diameter: m, volume mean
        real, intent(in):: hygro !hygroscopicity: volume mean in pop
        real, intent(in):: s !realtive humidity
        real, intent(in):: tair ! ambient temperature: K
        real, intent(out):: dwet_out !wet diameter: m
        real :: rdry_in, rwet_out !local variables
        integer :: i, n, nsol
        real :: a, b
        real :: p40,p41,p42,p43 ! coefficients of quartic polynomial
        real :: p30,p31,p32 ! coefficients of cubic polynomial
        real :: p
        real :: r3, r4
        real :: r       ! wet radius: um
        real :: rdry    ! radius of dry particle: microns
        real :: ss      ! relative humidity (1 = saturated) in 10^-4 to 1-epsilon
        real :: slog    ! log relative humidity
        real :: vol     ! total volume of particle: um^3
        real :: xi, xr
        complex :: cx4(4),cx3(3) !solution of quartic cx4/cubic cx3
        real, parameter :: eps = 1.e-4
        real, parameter :: mw = 18.
        real, parameter :: pi = 3.14159
        real, parameter :: rhow = 1.
        real, parameter :: surften = 72.
        real, parameter :: third = 1./3.
        real, parameter :: ugascon = 8.314e7
!        call mpp_sync()
!!$omp critical
!        write(mpp_pe()+100, *) ddry_in, hygro, s, tair
!!$omp end critical

        rdry_in = ddry_in/2
        !effect of organics on surface tension is neglected
        a=2.e4*mw*surften/(ugascon*tair*rhow) !in um, note: 2.e4 is due to ugascon in 10^7
        rdry = rdry_in*1.0e6   ! convert (m) to (microns)
        vol = rdry**3          ! vol is r**3, not volume: um^3
        b = vol*hygro
        !quartic coefficient based on taylor expansion
        ss=min(s,1.-eps)
        ss=max(ss,1.e-10)
        slog=log(ss)
        p43=-a/slog
        p42=0.
        !p41=b/slog-vol+b
        p41=b/slog-vol
        !p40=a*vol*(1-hygro)/slog
        p40=a*vol/slog
        !cubic coefficient for rh=1
        p32=0.
        p31=-b/a
        !p30=-vol+b
        p30=-vol

        !if particle radius is smaller than 1E-4 um -> 0.1 nm, don't consider hygroscopic growth 
        if(vol .le. 1.e-12) then !if particle is too small, here vol in um^3
            r=rdry
        else
            p=abs(p31)/(rdry*rdry)
            if(p.lt.eps)then
                r=rdry*(1.+p*third/(1.-slog*rdry/a))
            else
                call makoh_quartic(cx4,p43,p42,p41,p40)
                !dummy r to find solutions
                r = 1000*rdry
                nsol = 0
                !find the smallest real solution
                do n=1,4 !4 solutions
                xr=real(cx4(n))
                xi=aimag(cx4(n))
                if(abs(xi).gt.abs(xr)*eps) cycle
                if(xr.gt.r) cycle
                if(xr.lt.rdry*(1.-eps)) cycle
                if(xr.ne.xr) cycle
                r=xr
                nsol=n
                end do
                if (nsol .eq. 0) then !no solution found
                    r=rdry
                endif
            endif
        endif

        !if rh is large, put constrains
        if(s .gt. 1-eps) then
            !save quartic solution at s=1-eps
            r4=r
            p=abs(p31)/(rdry*rdry)
            if (p.lt.eps) then
                r=rdry*(1.+p*third)
            else
                call makoh_cubic(cx3,p32,p31,p30)
                !find smallest real(r8) solution
                r=1000.*rdry
                nsol=0
                do n=1,3
                xr=real(cx3(n))
                xi=aimag(cx3(n))
                if(abs(xi).gt.abs(xr)*eps) cycle
                if(xr.gt.r) cycle
                if(xr.lt.rdry*(1.-eps)) cycle
                if(xr.ne.xr) cycle
                r=xr
                nsol=n
                end do
                if(nsol.eq.0)then
                    r=rdry
                endif
            endif
            r3=r
            r=(r4*(1.-s)+r3*(s-1.+eps))/eps !now interpolate between quartic, cubic solutions
        endif

        ! bound and convert from microns to m
        !r = min(r,30.) ! upper bound based on 1 day lifetime
        rwet_out = r*1.e-6
        dwet_out = 2*rwet_out

    end subroutine 



    !-----------------------------------------------------------------------
    !      subroutine makoh_cubic( cx, p2, p1, p0, im )
    !      ----solves  x**3 + p2 x**2 + p1 x + p0 = 0
    !      ----where p0, p1, p2 are real
    !      ----modified from MAM SRC, tested by x5l MATHMATICA 
    !----------------------------------------------------------------------
    subroutine makoh_cubic(cx, p2, p1, p0)
        implicit none
        real,intent(in) :: p0, p1, p2 !coefficient of the equation, real
        complex, intent(out) :: cx(3) !3 solutions of the cubic equation, complex
        integer :: i
        real :: eps, q, r, sqrt3, third !lcoal variables
        complex :: ci, cq, crad, cw, cwsq, cy, cz
        eps = 1e-20 !epsilon, very small number
        third=1./3.
        sqrt3=sqrt(3.)
        ci=cmplx(0.,1.) !wait to check
        cw=0.5*(-1+ci*sqrt3)
        cwsq=0.5*(-1-ci*sqrt3)
        if(p1 .eq. 0.)then
            !completely insoluble particle
            cx(1)=(-p0)**third
            cx(2)=cx(1)
            cx(3)=cx(1)
        else
            q=p1/3.
            r=p0/2.
            crad=r*r+q*q*q
            crad=sqrt(crad)
            cy=r-crad
            if (abs(cy).gt.eps) cy=cy**third
            cq=q
            cz=-cq/cy
            cx(1)=-cy-cz
            cx(2)=-cw*cy-cwsq*cz
            cx(3)=-cwsq*cy-cw*cz
        endif

    end subroutine

    !-----------------------------------------------------------------------
    !      subroutine makoh_quartic( cx, p3, p2, p1, p0, im )
    !      ----solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
    !      ----where p0, p1, p2, p3 are real
    !      ----modified from MAM SRC, tested by x5l MATHMATICA
    !----------------------------------------------------------------------
    subroutine makoh_quartic(cx, p3, p2, p1, p0)
        implicit none
        real, intent(in) :: p0, p1, p2, p3
        complex, intent(out) :: cx(4)
        integer :: i
        real :: third, q, r
        complex :: cb, cb0, cb1, crad, cy, czero
        czero=cmplx(0.0,0.0)
        third=1./3.
        q=-p2*p2/36.+(p3*p1-4*p0)/12.
        r=-(p2/6.)**3+p2*(p3*p1-4*p0)/48.+(4.*p0*p2-p0*p3*p3-p1*p1)/16.
        crad=r*r+q*q*q
        crad=sqrt(crad)
        cb=r-crad
        if(cb.eq.czero)then
            !        insoluble particle
            cx(1)=(-p1)**third
            cx(2)=cx(1)
            cx(3)=cx(1)
            cx(4)=cx(1)
        else
            cb=cb**third
            cy=-cb+q/cb+p2/6.
            cb0=sqrt(cy*cy-p0)
            cb1=(p3*cy-p1)/(2.*cb0)
            cb=p3/2.+cb1
            crad=cb*cb-4.*(cy+cb0)
            crad=sqrt(crad)
            cx(1)=(-cb+crad)/2.
            cx(2)=(-cb-crad)/2.
            cb=p3/2.-cb1
            crad=cb*cb-4.*(cy-cb0)
            crad=sqrt(crad)
            cx(3)=(-cb+crad)/2.
            cx(4)=(-cb-crad)/2.
        endif
    end subroutine
end module
