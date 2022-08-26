subroutine hfunimplicit(M, q, p, pars, cont, h)
        !implicit none
        double precision, intent(in)  :: M, q(6), p(6), pars(11), cont
        double precision, intent(out) :: h

        !local variables
        double precision :: hhat, g, nnOrb, mu, palpha, alpha, deltaI(5), pdeltaI(5), dum, h_u, h_p, h_cp, h_cI
        double precision :: h_I, h_m, h_uu, h_pu, h_Iu, h_mu, ddeltaI(5), dpdeltaI(5)
        double precision :: I(5), bSail(3), sDir(3) 
        double precision :: Fx(6), Fy(6), Fz(6), pIG(3) 
        double precision ::  hSail
        double precision ::  fCone(2), uCone(3), hcone
 
        integer :: idx
        
        mu = pars(1)
        I           = pars(2 : 6)              ! Orbital elements
        fCone       = pars(7 : 8)              ! Forces on convex cone
        bSail       = pars(9 : 11)             ! Optical coefficients

        deltaI = q(1:5)
        alpha = q(6)
        pdeltaI = p(1:5)
        palpha = p(6) 
        
        sDir        = (/0.0D0, 0.0D0, -1.0D0/) ! Direction of the Sun
        
        ! GVEs
        call        gveeci(M, I, mu, Fx, Fy, Fz)

        call        dot(5, pdeltaI, Fx, pIG(1))
        call        dot(5, pdeltaI, Fy, pIG(2))
        call        dot(5, pdeltaI, Fz, pIG(3))
        
        ! Hamiltonian of the cone
        call        pmpcone(pIG, sDir, fCone, uCone)
        call        dot(3, pIG, uCone, hCone)
        
        
        !-------------------------------------------------
        
        ddeltaI = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
        dpdeltaI = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
        
        !call hfunalpha_d(M, deltaI, ddeltaI, pdeltaI, dpdeltaI, alpha, 1.0D0, dum, h_u)
        call hfunalpha_da_da(M, 0.0D0, deltaI, ddeltaI, pdeltaI, dpdeltaI, alpha, 1.0D0, 1.0D0, pars, cont, hhat, h_u, h_uu)
        
        !call hfunalpha_da(M, 1.0D0, deltaI, ddeltaI, pdeltaI, dpdeltaI, alpha, 0.0D0, pars, cont, dum, h_m)
        call hfunalpha_da_da(M, 1.0D0, deltaI, ddeltaI, pdeltaI, dpdeltaI, alpha, 1.0D0, 0.0D0, pars, cont, dum, h_m, h_mu)
        
        g = h_mu
        
        do idx   = 1,5
            ddeltaI = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
            dpdeltaI = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
            ddeltaI(idx) = 1.0D0
            call hfunalpha_da_da(M, 0.0D0,deltaI, ddeltaI, pdeltaI, dpdeltaI, alpha, 1.0D0, 0.0D0, pars, cont, dum, h_I, h_Iu)
           
            ddeltaI = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
            dpdeltaI = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
            dpdeltaI(idx) = 1.0D0
            call hfunalpha_da_da(M, 0.0D0, deltaI, ddeltaI, pdeltaI, dpdeltaI, alpha, 1.0D0, 0.0D0, pars, cont, dum, h_p, h_pu) 
            call hfuncone_dc(M, 0.0, pdeltaI, dpdeltaI, pars, dum, h_cp)
           
            g = g + (h_Iu * (h_cp * (1.0D0 - cont) + h_p * cont) - h_pu * h_I)       
        end do      
        
        g = - g / h_uu
           
        !if (hCone .GE. 0.0D0) then
        h       = hCone * (1.0D0 - cont) + hhat * cont + palpha * g 
        !else
         !   h       = 0.0D0
        !end if

end subroutine hfunimplicit

