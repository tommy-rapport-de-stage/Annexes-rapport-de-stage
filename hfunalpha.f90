!#######################################################################
!## Hamiltonian
subroutine hfunalpha(M, deltaI, pdeltaI, alpha, pars, cont, h)
        implicit none
        double precision, intent(in)  :: M, alpha, deltaI(5), pdeltaI(5), pars(11), cont
        double precision, intent(out) :: h

        !local variables
        double precision :: mu, I(5), bSail(3), sDir(3) 
        double precision :: Fx(6), Fy(6), Fz(6), pIG(3) 
        double precision :: uSail(3), u(3), nnOrb, hSail, aaOrb

        mu          = pars(1)                  ! Gravitational constant
        I           = pars(2 : 6)              ! Orbital elements
        bSail       = pars(9 : 11)             ! Optical coefficients

        sDir        = (/0.0D0, 0.0D0, -1.0D0/) ! Direction of the Sun

        ! GVEs
        call        gveeci(M, I, mu, Fx, Fy, Fz)

        call        dot(5, pdeltaI, Fx, pIG(1))
        call        dot(5, pdeltaI, Fy, pIG(2))
        call        dot(5, pdeltaI, Fz, pIG(3))
        
        
        ! Hamiltonian of the sail
        call        pmpsail(pIG, sDir,alpha, bSail, uSail)
        call        dot(3, pIG, uSail, hSail)
        
        h = hsail
            
        ! Mean longitude as time variable
        aaOrb        = pars(5)
        nnOrb       = sqrt(mu / aaOrb**3)
        h           = h / nnOrb
                        
        !!!!!!!!!
        !h = u(1)
        !!!!!!!!!
        
end subroutine hfunalpha

!Calcul of hcone
subroutine hfuncone(M, pdeltaI, pars, hcone)
        implicit none
        double precision, intent(in)  :: pdeltaI(5), M, pars(11)
        double precision, intent(out) :: hcone

        ! local variables
        double precision :: Fx(6), Fy(6), Fz(6), pIG(3), sDir(3), mu, fCone(2), uCone(3), I(5)
        
        mu = pars(1)
        I           = pars(2 : 6)              ! Orbital elements
        fCone       = pars(7 : 8)              ! Forces on convex cone
        
        sDir        = (/0.0D0, 0.0D0, -1.0D0/)

        call        gveeci(M, I, mu, Fx, Fy, Fz)

        call        dot(5, pdeltaI, Fx, pIG(1))
        call        dot(5, pdeltaI, Fy, pIG(2))
        call        dot(5, pdeltaI, Fz, pIG(3))
        
        ! Hamiltonian of the cone
        call        pmpcone(pIG, sDir, fCone, uCone)
        call        dot(3, pIG, uCone, hCone)

end subroutine hfuncone

!#######################################################################
!## PMP of the cone
subroutine pmpcone(pIG, sDir, fCone, u)
        implicit none
        double precision, intent(in)  :: pIG(3), sDir(3), fCone(2)
        double precision, intent(out) :: u(3)

        ! local variables
        double precision :: sPerp(3)

        sPerp   = pIG - (pIG(1) * sDir(1) + pIG(2) * sDir(2) + pIG(3) &
&                 * sDir(3)) * sDir
        sPerp   = sPerp / sqrt(sPerp(1)**2 + sPerp(2)**2 + sPerp(3)**2)

        u       = fCone(1) * sDir + fCone(2) * sPerp

end subroutine pmpcone

!#######################################################################
!# PMP of the sail
subroutine pmpsail(pIG, sDir, alpha, b, u)
        implicit none
        double precision, intent(in)  :: pIG(3), alpha, sDir(3), b(3)
        double precision, intent(out) :: u(3)

        ! local variables
        double precision :: b1, b2, b3, pIs, pInorm, theta, Kpolar, &
&                           sPerp(3), cAlpha, sAlpha, fs, fperp

        b1              = b(1)
        b2              = b(2)
        b3              = b(3)

        pIs             = pIG(1) * sDir(1) + pIG(2) * sDir(2) &
&                         + pIG(3) * sDir(3)
        pInorm          = sqrt(pIG(1)**2 + pIG(2)**2 + pIG(3)**2)
        theta           = acos(pIs / pInorm)
        
        ! Definition of the unit vector sPerp
        sPerp           = pIG - pIs * sDir
        sPerp           = sPerp / sqrt(sPerp(1)**2 + sPerp(2)**2 + &
&                         sPerp(3)**2) ! Note: nan for collinear vectors

        cAlpha          = cos(alpha)
        sAlpha          = sin(alpha)

        ! Optimal force
        fs              = b1 * cAlpha + (b2 * cAlpha**2 + b3 * cAlpha) &
&                         * cAlpha
        fperp           = (b2 * cAlpha**2 + b3 * cAlpha) * sAlpha

        u               = fs * sDir + fperp * sPerp
        
end subroutine pmpsail


! ######################################################################
! ## GVE in ECI frame
subroutine gveeci(M, I, mu, Fx, Fy, Fz)
        implicit none
        double precision, intent(in)  :: M, I(5), mu
        double precision, intent(out) :: Fx(6), Fy(6), Fz(6)

        ! Local variables
        double precision :: R(6), T(6), N(6), f, Om, inc, w, a, e, &
&                           theta, sOm, cOm, sI, cI, sTh, cTh

        ! Orbital elements
        Om       = I(1)
        inc      = I(2)
        w        = I(3)
        a        = I(4)
        e        = I(5)

        ! GVE in LVLH frame
        call     gvelvlh(M, I, mu, R, T, N, f)
        theta    = w + f

        sOm      = SIN(Om)
        cOm      = COS(Om)
        sI       = SIN(inc)
        cI       = COS(inc)
        sTh      = SIN(theta)
        cTh      = COS(theta)

        ! GVE in ECI frame
        Fx         = R * ((- sOm * cI * sTh + cOm * cTh)) + &
&                    T * ((- sOm * cI * cTh - cOm * sTh)) + &
&                    N * ((sOm * sI)                    )
        Fy         = R * ((  cOm * cI * sTh + sOm * cTh)) + &
&                    T * ((  cOm * cI * cTh - sOm * sTh)) + &
&                    N * (- (cOm * sI)                  )
        Fz         = R * (sI * sTh) + T * (sI * cTh) + N * cI

end subroutine gveeci

! ######################################################################
! ## GVE in LVLH frame
subroutine gvelvlh(M, I, mu, R, T, N, f)
        implicit none
        double precision, intent(in)  :: M, I(5), mu
        double precision, intent(out) :: R(6), T(6), N(6), f

        ! Local variables
        double precision :: inc, w, a, e, cF, sF, p, rad, b, nnOrb, h, theta

        ! Orbital elements
        ! Omega  = I(1)
        inc      = I(2)
        w        = I(3)
        a        = I(4)
        e        = I(5)

        ! Solving Kepler's equation
        call     kepler(e, M, f)

        ! Useful variables
        cF       = COS(f)
        sF       = SIN(f)
        p        = a * (1.0D0 - e**2)
        rad      = p / (1.0D0 + e * cF)
        b        = a * SQRT(1.0D0 - e**2)
        nnOrb     = SQRT(mu / a**3)
        h        = nnOrb * a * b
        theta    = w + f

        ! GVEs
        R(1)     = 0.0D0
        T(1)     = 0.0D0
        N(1)     = rad * SIN(theta) / h / SIN(inc)

        R(2)     = 0.0D0
        T(2)     = 0.0D0
        N(2)     = rad * COS(theta) / h

        R(3)     = - p * cF / h / e
        T(3)     = (p + rad) * sF / h / e
        N(3)     = - rad * SIN(theta) * COS(inc) / h / SIN(inc)

        R(4)     = 2.0D0 * a**2 * e * sF / h
        T(4)     = 2.0D0 * a**2 * p / h / rad
        N(4)     = 0.0D0

        R(5)     = p * sF / h
        T(5)     = ((p + rad) * cF + rad * e) / h
        N(5)     = 0.0D0

        R(6)     = p * cF - 2.0D0 * rad * e * b / a / h / e
        T(6)     =  - (p + rad) * sF * b / a / h / e
        N(6)     = 0.0D0

end subroutine gvelvlh


! ######################################################################
! ## Kepler's equation
subroutine kepler(ecc, M, f)
        implicit none
        double precision, intent(in)  :: ecc, M
        double precision, intent(out) :: f

        ! Local variables
        integer          :: Nmax, ii
        double precision :: E, K, dK

        ! Set parameters
        Nmax     = 20

        ! Newton-Rapshon
        E        = M
        do ii    = 1, Nmax
            K    = E - ecc * SIN(E) - M
            dK   = 1.0D0 - ecc * COS(E)
            E    = E - K / dK
        end do

        ! True Anomaly
        f        = 2.0D0 * ATAN2(SQRT(1.0D0 + ecc) * SIN(E / 2.0D0), &
&                                SQRT(1.0D0 - ecc) * COS(E / 2.0D0))

end subroutine kepler

! ######################################################################
! ## Dot product
subroutine dot(n,u,v,res)
    implicit none
    integer, intent(in)             :: n
    double precision, intent(in)    :: u(n), v(n)
    double precision, intent(out)   :: res

    ! local variables
    integer :: i

    res = 0.0D0
    do i=1,n
        res = res + u(i)*v(i)
    end do

end subroutine dot

!#######################################################################
!## control
subroutine control(M, deltaI, pdeltaI, alpha, pars, cont, u)
        implicit none
        double precision, intent(in)  :: M, deltaI(5), pdeltaI(5), alpha, pars(11), cont
        double precision, intent(out) :: u(3)

        !local variables
        double precision :: mu, I(5), fCone(2), bSail(3), sDir(3), &
&                           Fx(6), Fy(6), Fz(6), pIG(3), uCone(3), &
&                           uSail(3), hCone

        mu          = pars(1)                  ! Gravitational constant
        I           = pars(2 : 6)              ! Orbital elements
        fCone       = pars(7 : 8)              ! Forces on convex cone
        bSail       = pars(9 : 11)             ! Optical coefficients

        sDir        = (/0.0D0, 0.0D0, -1.0D0/) ! Direction of the Sun

        ! GVEs
        call        gveeci(M, I, mu, Fx, Fy, Fz)

        call        dot(5, pdeltaI, Fx, pIG(1))
        call        dot(5, pdeltaI, Fy, pIG(2))
        call        dot(5, pdeltaI, Fz, pIG(3))

        ! Hamiltonian of the cone
        call        pmpcone(pIG, sDir, fCone, uCone)
        call        dot(3, pIG, uCone, hCone)
        
        ! Hamiltonian of the sail
        call        pmpsail(pIG, sDir, alpha, bSail, uSail)
        
        ! Hamiltonian (note that cont is the continuation parameter)
        u           = (/0.0D0, 0.0D0, -0.0D0/)
        if (hCone .GE. 0.0D0) then
            u       = uCone * (1.0D0 - cont) + uSail * cont
        end if                        
        
end subroutine control


