
subroutine subsidence()
 

use vars
use params
use advect_mod
implicit none

integer :: k

   ! advect_rho needs to be called first since it sets rho
   ! on interfaces
   call advect_rho(rho)
   tend_sub_u   = advect_mass_fraction2(u,rho,limiter=.false.)
   tend_sub_v   = advect_mass_fraction2(v,rho,limiter=.false.)
   tend_sub_t   = advect_mass_fraction2(t,rho,limiter=dothuburn)
   tend_sub_qt  = advect_mass_fraction2(qt,rho,limiter=dothuburn)
   tend_sub_qp  = advect_mass_fraction2(qp,rho,limiter=dothuburn)
   tend_sub_tke = advect_mass_fraction2(tke,rho,limiter=dothuburn)

if (dosequential) then

  t=t + dt*tend_sub_t
  qt=qt + dt*tend_sub_qt
  qp=qp + dt*tend_sub_qp
  u=u + dt*tend_sub_u
  v=v + dt*tend_sub_v
  tke=tke + dt*tend_sub_tke

end if

end subroutine subsidence

