module advect_mod

   implicit none


   contains

   
      !==============================================!
      ! Function advect_rho                          !
      !                                              !
      ! Advect_rho returns the tendency of rho from  !
      ! from transport (advection and convergence).  !
      ! This function must be called before calls to !
      ! advect_mass_fraction and advect_mixing_ratio !
      ! because it sets rho on the faces.            !
      !==============================================!
      function advect_rho(rho)

         use grid
         use vars, only: w, rho_i
       
         implicit none

         ! In
         real, intent(in), dimension(:) :: rho

         ! Out
         real, dimension(nzm) :: advect_rho

         integer :: k, kp
         real :: one_over_dz
            
         ! Find the values of rho on the faces
         call thirdorder(rho,rho_i)

         advect_rho(:) = 0.

         ! Calculate the change in rho
         advect_rho(1) = - w(2)*rho_i(2) * (1./adz(1)/dz)
         do k = 2,nzm-1
            kp = k + 1
            one_over_dz = 1./adz(k)/dz
                  advect_rho(k) = advect_rho(k) + &
                    ( w(k ) * rho_i(k ) - &
                      w(kp) * rho_i(kp) ) * one_over_dz 
         end do
         advect_rho(nzm) = w(nzm)*rho_i(nzm) * (1./adz(nzm)/dz)
         advect_rho = advect_rho / rho

      end function advect_rho


      !==============================================!
      ! Function advect_mass_fraction                !
      !                                              !
      ! Advect_density returns the tendency of q*rho !
      ! from transport (advection and convergence).  !
      !==============================================!
      function advect_mass_fraction(q,rho,limiter,zflux)

         use grid
         use vars, only:  w, rho_i
         use thuburn_mod
         use params

         implicit none

         ! In
         real, intent(in), dimension(:) :: q, rho
         logical, intent(in), optional :: limiter   ! Defaults to true

         ! Out
         real, dimension(nzm) :: advect_mass_fraction
         real, dimension(nz), optional :: zflux

         real, dimension(2:nzm) :: d_i, q_i
         integer :: k, kp
         real :: one_over_dz
         logical :: limiter2

         limiter2 = .true.
         if (present(limiter)) then
            limiter2 = limiter
         end if

         call thirdorder(q,q_i)

         ! Apply a flux limiter
         if (limiter2) then
            call thuburn_vertical(q,q_i)
            d_i = rho_i*q_i
         else
            d_i = rho_i*q_i
         end if

         if (present(zflux)) then
            zflux(1    ) = 0.
            zflux(2:nzm) = d_i(2:nzm)*w(2:nzm)
            zflux(nz ) = 0.
         end if
         advect_mass_fraction(:) = 0.

         ! Calculate the change in d
         advect_mass_fraction(1) = -  w(  2     )*d_i(  2     ) * (1./adz(1)/dz)
         do k = 2,nzm-1
            kp = k + 1
            one_over_dz = (1./adz(k)/dz)
                  advect_mass_fraction(k) = advect_mass_fraction(k) + &
                    ( w(k ) * d_i(k ) - &
                      w(kp) * d_i(kp) ) * one_over_dz 
         end do
         advect_mass_fraction(nzm) =  w( nzm)*d_i(nzm) * (1./adz(nzm)/dz)

         advect_mass_fraction = advect_mass_fraction / rho

      end function advect_mass_fraction

      !==============================================!
      ! Function advect_mass_fraction2                !
      !                                              !
      ! Advect_density returns the tendency of q*rho !
      ! from transport (advection only).  !
      !==============================================!
      function advect_mass_fraction2(q,rho,limiter,zflux)

         use grid
         use vars, only:  w, rho_i
         use thuburn_mod
         use params

         implicit none

         ! In
         real, intent(in), dimension(:) :: q, rho
         logical, intent(in), optional :: limiter   ! Defaults to true

         ! Out
         real, dimension(nzm) :: advect_mass_fraction2
         real, dimension(nz), optional :: zflux

         real, dimension(1:nz) :: q_i
         integer :: k, kp
         real :: one_over_dz
         logical :: limiter2

         limiter2 = .true.
         if (present(limiter)) then
            limiter2 = limiter
         end if

         call thirdorder_adv(q,q_i)

         ! Apply a flux limiter
         if (limiter2) then
            call thuburn_vertical(q,q_i(2:nzm))
            !d_i = rho_i*q_i
         !else
         !   d_i = rho_i*q_i
         end if

         !if (present(zflux)) then
         !   zflux(1    ) = 0.
         !   zflux(2:nzm) = d_i(2:nzm)*w(2:nzm)
         !   zflux(nz ) = 0.
         !end if
         advect_mass_fraction2(:) = 0.

         ! Calculate the change in d
         advect_mass_fraction2(1) = -  0.5*w(  2     )*(q_i(2)-q_i(1)) * (1./adz(1)/dz)
         do k = 2,nzm-1
            kp = k + 1
            one_over_dz = (1./adz(k)/dz)
                  advect_mass_fraction2(k) = advect_mass_fraction2(k) + &
                    0.5 * (w(k )+w(kp))*(q_i(k ) - q_i(kp) ) * one_over_dz 
         end do
         advect_mass_fraction2(nzm) =  0.5*w( nzm)*(q_i(nzm)-q_i(nz)) * (1./adz(nzm)/dz)

      end function advect_mass_fraction2


      !===============================================!
      ! Subroutine thirdorder_adv                         !
      !                                               !
      ! Thirdorder obtains the one-dimensional upwind !
      ! third-order interpolation of the sss quantity !
      ! q onto the iss, sis, and ssi faces.  If the   !
      ! optional flag positive is set to true, then   !
      ! the face values are forced to be nonnegative. !
      !===============================================!
      ! difference to original thirdorder below is that 
      subroutine thirdorder_adv(q,q_i,positive)

         use vars, only : w
         use grid

         implicit none

         ! In
         real, intent(in), dimension(1:) :: q
         logical, intent(in), optional :: positive

         ! Out
         real, intent(out), dimension(1:nz) :: q_i   ! 

         integer :: k, ku, kc, kd
         real, parameter :: au = -1./6., ac = 5./6., ad = 1./3.
         logical :: positive2
         real :: temp

         if (present(positive)) then
            positive2 = positive
         else
            positive2 = .false.
         end if


         ! q_i
         ! this makes 1.order upstream in case of subsidence
         q_i(1) =  3./2.*q(1) - 1./2.*q(2)
         q_i(2) = 0.5*( q(1) + q(2) )
         do k = 3,nzm-1
              ku = k - int(0.5 + sign(1.5,w(k)))
              kc = k - int(0.5 + sign(0.5,w(k)))
              kd = k - int(0.5 - sign(0.5,w(k)))
              q_i(k) =        au*q(ku) + &
                              ac*q(kc) + &
                              ad*q(kd) 
         end do
         q_i(nzm) = 0.5*( q(nzm-1) + q(nzm) )
         ! this makes zero gradient 
         q_i(nz)  = q_i(nzm)

         ! Enforce nonnegativity
         if (positive2) then
            do k = 2,nzm
                temp = 0.5 + sign(0.5,q_i(k))
                q_i(k) = temp * q_i(k)
            end do
         end if

      end subroutine thirdorder_adv

      !===============================================!
      ! Subroutine thirdorder                         !
      !                                               !
      ! Thirdorder obtains the one-dimensional upwind !
      ! third-order interpolation of the sss quantity !
      ! q onto the iss, sis, and ssi faces.  If the   !
      ! optional flag positive is set to true, then   !
      ! the face values are forced to be nonnegative. !
      !===============================================!
      subroutine thirdorder(q,q_i,positive)

         use vars, only : w
         use grid

         implicit none

         ! In
         real, intent(in), dimension(1:) :: q
         logical, intent(in), optional :: positive

         ! Out
         real, intent(out), dimension(2:) :: q_i   ! 0:(nx+1),0:(ny+1),2:nzm

         integer :: k, ku, kc, kd
         real, parameter :: au = -1./6., ac = 5./6., ad = 1./3.
         logical :: positive2
         real :: temp

         if (present(positive)) then
            positive2 = positive
         else
            positive2 = .false.
         end if


         ! q_i
         q_i(2) = 0.5*( q(1) + q(2) )
         do k = 3,nzm-1
              ku = k - int(0.5 + sign(1.5,w(k)))
              kc = k - int(0.5 + sign(0.5,w(k)))
              kd = k - int(0.5 - sign(0.5,w(k)))
              q_i(k) =        au*q(ku) + &
                              ac*q(kc) + &
                              ad*q(kd) 
         end do
         q_i(nzm) = 0.5*( q(nzm-1) + q(nzm) )

         ! Enforce nonnegativity
         if (positive2) then
            do k = 2,nzm
                temp = 0.5 + sign(0.5,q_i(k))
                q_i(k) = temp * q_i(k)
            end do
         end if

      end subroutine thirdorder


      !=========================!
      ! Subroutine coriolis_add !
      !=========================!
      subroutine coriolis_add(SU, SV)

         use vars, only : u,v,ug,vg,rho
         use params
         use grid

         implicit none

         ! In/Out
         real, intent(inout), dimension(nzm) :: SU, SV   ! 1:nx,1:ny,1:nzm
      
         !                        |---- To average rho
         !                        |
         SU = SU + fcor * (v(1:nzm)-vg(1:nzm))
         SV = SV - fcor * (u(1:nzm)-ug(1:nzm))

      end subroutine coriolis_add

end module advect_mod

   
