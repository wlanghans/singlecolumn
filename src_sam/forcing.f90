
subroutine forcing()
 

use vars
use params
use advect_mod
use grid
implicit none

integer :: k

   if (docoriolis) call coriolis_add(tend_force_u,tend_force_v)

   if (doforcing) then

      tend_force_t      = cp * dthetadt
      tend_force_qt     = dqtdt
      tend_force_tke    = 0. 

   end if

if (dosequential) then

  t  = t   + dt*tend_force_t
  qt = qt  + dt*tend_force_qt
  u  = u   + dt*tend_force_u
  v  = v   + dt*tend_force_v
  tke= tke + dt*tend_force_tke

end if

end subroutine forcing

