! ==================================================================
      SUBROUTINE get_abcd(rhoin,s,sumMs,Cs,tkin, ssfc,a,b,c,d, massflux, dosfcbcneuman, kinflx)

      use grid
      use params
      use vars
      implicit none


!! to solve system of linear eqs on tridiagonal matrix n times n
!! after Peaceman and Rachford, 1955
!! we need a,b,c,d - which all are vectors of order n 
!! a,b,c - are coefficients on the LHS
!! d - is initially RHS on the output becomes a solution vector

!-------------------------------------------------------------------

        ! input
        REAL, DIMENSION(nzm),    INTENT(in)  :: rhoin        ! density, S at cell center
        REAL, DIMENSION(nzm),    INTENT(in)  :: s, tkin    ! S and eddy viscosity at cell center
        REAL, DIMENSION(nz    ), INTENT(in)  :: sumMs      ! sum(MiSi) on faces
        REAL,                    INTENT(in)  :: Cs,ssfc    ! drag coefficient, surface value of S
        LOGICAL,                 INTENT(in)  :: dosfcbcneuman ! If true, then neuman bc conditions, dirichlet otherwise
        LOGICAL,                 INTENT(in)  :: massflux   ! If true, then sumM is used otherwise sumM=0; note that sumMs needs to be zero then
        REAL, INTENT(in)   ,     OPTIONAL    :: kinflx
        ! output
        REAL, DIMENSION(nzm),    INTENT(out) :: a,b,c,d
        

        ! local
        INTEGER :: k
        REAL, DIMENSION(nz) :: tkf, rhof, sumMloc                   ! tk and rho on face
 
        tkf(1)  = Cs
        tkf(nz) = 0.
        do k=2,nzm
           tkf(k)  = 0.5*(tkin(k-1)+tkin(k))
           rhof(k) = 0.5*(rhoin(k-1)+rhoin(k))
        end do
        rhof(1) = 2.*rhof(2) - rhof(3)
        rhof(nz)= 2.*rhof(nzm) - rhof(nzm-1)

        if ((dosfcbcneuman).and..not.present(kinflx)) then
          write(*,*) 'ERROR: Surface fluxes need to be passed in case' 
          write (*,*)'of Neuman BC.'
          write (*,*)'stopping'
          stop
        end if

        if (massflux) then
          sumMloc = sumM
        else
          sumMloc = 0.0
        end if


        DO k=1,nzm
           
           
           IF (k.gt.1.and.k.lt.nzm) THEN
              a(k) = rhof(k)/rhoin(k)/adz(k)/dz * betap *           &
                    (tkf(k)/adzw(k)/dz)
              !a(k) = rhof(k)/rhoin(k)/adz(k)/dz *                   &
              !      (betap*tkf(k)/adzw(k)/dz - 0.5 * sumMloc(k))
              b(k) = -1./dt + betap* rhof(k+1)/rhoin(k)/adz(k)/dz * & 
                    (-tkf(k+1)/adzw(k+1)/dz ) - & 
                    betap*rhof(k)/rhoin(k)/adz(k)/dz*               &
                    (tkf(k)/adzw(k)/dz + sumMloc(k))
              !b(k) = -1./dt +  rhof(k+1)/rhoin(k)/adz(k)/dz * & 
              !      (-betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1) ) - & 
              !       rhof(k)/rhoin(k)/adz(k)/dz*               &
              !      (betap*tkf(k)/adzw(k)/dz + 0.5 * sumMloc(k))
              c(k) = rhof(k+1)/rhoin(k)/adz(k)/dz * betap *         &
                     (tkf(k+1)/adzw(k+1)/dz + sumMloc(k+1))
              !c(k) = rhof(k+1)/rhoin(k)/adz(k)/dz *              &
              !       (betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1))
              d(k) = -s(k)/dt - 1./rhoin(k)/adz(k)/dz * (                       &
                    tkf(k+1)*rhof(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
                   -tkf(k)*rhof(k)/adzw(k)/dz*betam*(s(k)-s(k-1))         &
                   +(rhof(k)*(sumMs(k) - betam*s(k)*sumMloc(k)))&
                   -(rhof(k+1)*(sumMs(k+1) - betam*s(k+1)*sumMloc(k+1))) )
              !d(k) = -s(k)/dt - 1./rhoin(k)/adz(k)/dz * (                       &
              !      tkf(k+1)*rhof(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
              !     -tkf(k)*rhof(k)/adzw(k)/dz*betam*(s(k)-s(k-1))         &
              !     +(rhof(k)* sumMs(k))&
              !     -(rhof(k+1)*sumMs(k+1)) )
           ELSEIF (k.eq.1) THEN
              a(k) = 0.0
              b(k) = -1./dt + betap* rhof(k+1)/rhoin(k)/adz(k)/dz * &
                    (-tkf(k+1)/adzw(k+1)/dz ) 
              !b(k) = -1./dt +  rhof(k+1)/rhoin(k)/adz(k)/dz * &
              !      (-betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1) ) 
              if (.not.dosfcbcneuman) then
                  b(k) = b(k) - betap*vmag*tkf(1)/adz(k)/dz
              end if
              c(k) = rhof(k+1)/rhoin(k)/adz(k)/dz * betap *         &
                     (tkf(k+1)/adzw(k+1)/dz + sumMloc(k+1))
              !c(k) = rhof(k+1)/rhoin(k)/adz(k)/dz *                  &
              !       (betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1))
              d(k) = -s(k)/dt - 1./rhoin(k)/adz(k)/dz * (                       &
                    tkf(k+1)*rhof(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
                   -rhof(k+1)*(sumMs(k+1) - betam*s(k+1)*sumMloc(k+1)))
              !d(k) = -s(k)/dt - 1./rhoin(k)/adz(k)/dz * (                       &
              !      tkf(k+1)*rhof(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
              !     -rhof(k+1)*sumMs(k+1) )
              if (.not.dosfcbcneuman) then
                  d(k) = d(k) + rhof(k)/rhoin(k)*vmag*tkf(1)/adz(k)/dz*(betam*s(k)-ssfc)
              else ! neuman bc, always fully explicit right now
                  d(k) = d(k) - rhof(k)/rhoin(k)*kinflx/adz(k)/dz
              end if
           ELSE
              a(k) = rhof(k)/rhoin(k)/adz(k)/dz * betap *                  &
                    (tkf(k)/adzw(k)/dz)
              !a(k) = rhof(k)/rhoin(k)/adz(k)/dz *                           &
              !      (betap*tkf(k)/adzw(k)/dz - 0.5 * sumMloc(k))
              b(k) = -1./dt - betap* rhof(k)/rhoin(k)/adz(k)/dz *          &
                    (tkf(k)/adzw(k)/dz + sumMloc(k) )
              !b(k) = -1./dt - rhof(k)/rhoin(k)/adz(k)/dz *          &
              !      (betap * tkf(k)/adzw(k)/dz + 0.5 * sumMloc(k) )
              c(k) = 0.0
              d(k) = -s(k)/dt + 1./rhoin(k)/adz(k)/dz * (              &
                    tkf(k)*rhof(k)/adzw(k)/dz*betam*(s(k)-s(k-1))&
                   -rhof(k)*(sumMs(k) - betam*s(k)*sumMloc(k)))
              !d(k) = -s(k)/dt + 1./rhoin(k)/adz(k)/dz * (              &
              !      tkf(k)*rhof(k)/adzw(k)/dz*betam*(s(k)-s(k-1))&
              !     -rhof(k)*sumMs(k) )
           END IF
        ENDDO


      END SUBROUTINE get_abcd

