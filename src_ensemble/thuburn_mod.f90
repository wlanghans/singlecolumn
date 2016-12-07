!==================================!
! Das Atmosphaerische Modell (DAM) !
! Copyright (c) David M. Romps     !
! http://www.romps.org             !
!==================================!

module thuburn_mod

   implicit none

   contains

      !================================================================!
      ! Subroutine thuburn_vertical                                    !
      !                                                                !
      ! Thuburn_vertical applies the multidimensional, conservative,   !
      ! shape-preserving flux limiter to the values of                 !
      ! q_w (a.k.a., q_ssi), where q_w is a *mass density*.            !
      !                                                                !
      ! @article{thuburn1996,                                          !
      !   author = {Thuburn, John},                                    !
      !   title = {{Multidimensional flux-limited advection schemes}}, !
      !   journal = {Journal of Computational Physics},                !
      !   year = {1996},                                               !
      !   volume = {123},                                              !
      !   pages = {74--83},                                            !
      !   number = {1},                                                !
      !   publisher = {Academic Press Professional,                    !
      !                 Inc. San Diego, CA, USA}                       !
      ! }                                                              !
      !                                                                !
      ! If the optional variables rho, rho_iss, rho_sis, and rho_ssi   !
      ! are provided as input, then this subroutine treats q as a mass !
      ! fraction and it bounds the mass fractions q_u, q_v, and q_w as !
      ! described in section 4 of thuburn1996.  Otherwise, q is        !
      ! treated as a density and the densities q_u, q_v, and q_w are   !
      ! bounded as described in section 6 of thuburn1996.              !
      !================================================================!
      subroutine thuburn_vertical(q,q_w)

         use grid
         use vars, only: w

         implicit none
         
         ! In
         real, intent(in), dimension(:) :: q   ! 1:nzm

         ! In/out
         real, intent(inout), dimension(2:) :: q_w   ! 2:nzm
      
         real, dimension(:), allocatable :: &
            q_w_in_min, q_w_in_max, &
            q_out_min, q_out_max, q_new_min, q_new_max, c_in_sum, c_out_sum
         real :: temp
         integer :: k, ku, kc, kd, km, kp, abflag
         !===================================================================!
         ! The last index of cz takes the value 1 (below) and 2 (above):     !
         !                                                                   !
         ! cz(:,1) are the Courant numbers on                                !
         !               the ssi faces 2:nzm below the sss cells 2: nzm      !
         ! cz(:,2) are the Courant numbers on                                !
         !               the ssi faces 2:nzm above the sss cells 1:(nzm-1)   !
         !===================================================================!
         real, dimension(:,:), allocatable :: cz   ! 2:nzm,1:2

         ! Check that array sizes are consistent
         if ( size(q_w) .ne. nzm-1 .or. &
              size(w  ) .ne. nz ) then
            print *,'Error in thuburn_vertical: dimensions of arrays are not correct'
            write(*,'(a,i3,a,i3,a,i3)') '   size(q_w) = ',size(q_w),' and nzm-1  = ',nzm-1
            write(*,'(a,i3,a,i3,a,i3)') '   size(w  ) = ',size(w  ),' and nz     = ',nz
            stop
         end if

         allocate( q_w_in_min(2:nzm), &
                   q_w_in_max(2:nzm), &
                   q_out_min (1:nzm), & 
                   q_out_max (1:nzm), & 
                   q_new_min (1:nzm), & 
                   q_new_max (1:nzm), & 
                   c_in_sum  (1:nzm), & 
                   c_out_sum (1:nzm) )
      
         allocate( cz(2:nzm,2) )


         ! Define the Courant numbers
         do k = 2,nzm
            cz(k,1) = abs(w(k)) * (dt * (1./adz(k)/dz) )   ! cz is below the sss cell
         end do
         do k = 2,nzm
            cz(k,2) = abs(w(k)) * (dt * (1./adz(k-1)/dz) )   ! cz is above the sss cell
         end do

   
         !================================!
         ! Define q_w_in_min(max) (32,33) !
         !================================!
   
   
         ! Define q_w_in_min and q_w_in_max
         do k = 2,nzm
            ku = k - int(0.5 + sign(0.5,w(k)))
            kd = k - int(0.5 - sign(0.5,w(k)))
            q_w_in_min(k) = min(q(kd),q(ku))
            q_w_in_max(k) = max(q(kd),q(ku))
         end do
   
   
         !===================!
         ! Bound q_w (34,35) !
         !===================!
   
   
         ! Bound q_w
         do k = 2,nzm
            q_w(k) = max(q_w_in_min(k),min(q_w_in_max(k),q_w(k)))
         end do
   

         !=======================!
         ! Define q_new_min (38) !
         !=======================!
   
   
         ! Set sss q_new_min(max) equal to the sss cell value
         q_new_min = q(1:nzm)
         q_new_max = q(1:nzm)

         ! Set sss q_new_min(max) to the min (max) of itself and inflowing q_w_in_min(max)
         do k = 2,nzm
            ! Set kd such that kd is the sss cell downwind of the k ssi face
            kd = k - int(0.5 - sign(0.5,w(k)))
            q_new_min(kd) = min(q_new_min(kd),q_w_in_min(k))
            q_new_max(kd) = max(q_new_max(kd),q_w_in_max(k))
         end do
   
   
         !==================================!
         ! Redefine q_w_in_min(max) (48,49) !
         !==================================!
               
   
         ! Define q_w_in_min(max) as the min (max) of q_up_min(max)_w and q_w
         do k = 2,nzm
            ku = k - int(0.5 + sign(0.5,w(k)))
            q_w_in_min(k) = min(q(ku),q_w(k))
            q_w_in_max(k) = max(q(ku),q_w(k))
         end do
   
   
         !===============================!
         ! Define c_in_sum and c_out_sum !
         !===============================!
   
   
         ! Add to c_in(out)_sum the Courant numbers from all faces
         c_in_sum(2:(nzm-1)) = cz(2:(nzm-1),1) + &  ! below
                               cz(3: nzm   ,2)      ! above
         c_in_sum(1        ) = cz(2        ,2)      ! above
         c_in_sum(   nzm   ) = cz(   nzm   ,1)      ! below
         c_out_sum = c_in_sum
   
   
         ! For c_in(out)_sum, add (subtract) each Courant number
         ! if the flux is in (out), otherwise subtract (add) it
         do k = 2,nzm-1
            kp = k + 1
            temp = sign(cz(k ,1),w(k )) - &   ! below
                   sign(cz(kp,2),w(kp))       ! above
            c_in_sum (k) = c_in_sum (k) + temp
            c_out_sum(k) = c_out_sum(k) - temp
         end do
         c_in_sum (1  ) = c_in_sum (1  ) - sign(cz(2  ,2),w(2  ))   ! above
         c_out_sum(1  ) = c_out_sum(1  ) + sign(cz(2  ,2),w(2  ))   ! above
         c_in_sum (nzm) = c_in_sum (nzm) + sign(cz(nzm,1),w(nzm))   ! below
         c_out_sum(nzm) = c_out_sum(nzm) - sign(cz(nzm,1),w(nzm))   ! below
   
         ! Mutiply c_in(out)_sum by 1/2
         c_in_sum  = 0.5*c_in_sum
         c_out_sum = 0.5*c_out_sum
   
   
         !===============================!
         ! Define q_out_min(max) (42,43) !
         !===============================!
   
   
         q_out_min = q(1:nzm) - q_new_max*(1. + c_in_sum - c_out_sum)
         q_out_max = q(1:nzm) - q_new_min*(1. + c_in_sum - c_out_sum)
   
         ! Add inflowing "cq" in the w direction
         do k = 2,nzm
            ! Set kd such that the sss cell kd is downwind of the ssi face k
            kd = k - int(0.5 - sign(0.5,w(k)))
            abflag = 1 + int(0.5 - sign(0.5,w(k)))   ! 1 for below, 2 for above
            q_out_min(kd) = q_out_min(kd) + cz(k,abflag) * q_w_in_max(k)
            q_out_max(kd) = q_out_max(kd) + cz(k,abflag) * q_w_in_min(k)
         end do

         where (c_out_sum .gt. 0.)
            q_out_min = q_out_min / c_out_sum
            q_out_max = q_out_max / c_out_sum
         end where
   
   
         !===================!
         ! Bound q_w (44,45) !
         !===================!
         
   
         ! Bound q_w
         do k = 2,nzm
            ! Set ku such that the sss cell ku is upwind of the ssi face k
            ku = k - int(0.5 + sign(0.5,w(k)))
            q_w(k) = min(q_out_max(ku),max(q_out_min(ku),q_w(k)))
         end do

      end subroutine thuburn_vertical

end module thuburn_mod
