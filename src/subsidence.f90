
subroutine subsidence()
 

use vars
use params
use advect_mod
implicit none

real, dimension(nzm) :: rho_new
integer :: k

   ! advect_rho needs to be called first since it sets rho
   ! on interfaces
   tend_sub_rho     = advect_rho(rho)
   tend_sub_rho_u   = advect_mass_fraction2(u,rho,limiter=.false.)
   tend_sub_rho_v   = advect_mass_fraction2(v,rho,limiter=.false.)
   tend_sub_rho_thetav = advect_mass_fraction2(thetav,rho,limiter=dothuburn)
   tend_sub_rho_qv     = advect_mass_fraction2(qv,rho,limiter=dothuburn)
   tend_sub_rho_tke    = advect_mass_fraction2(tke,rho,limiter=dothuburn)
   ! now assume that density doesn't change to due mass continuity. 
   ! Only changes in water vapor mass then affect the density
   tend_sub_rho        = tend_sub_rho_qv
   

if (dosequential) then

  rho_new = rho + dt*tend_sub_rho_qv
  ! get u, v
  u=(rho*u + dt*tend_sub_rho_u)/rho_new
  v=(rho*v + dt*tend_sub_rho_v)/rho_new

! get thermo state
  thetav= (rho*thetav + dt * tend_sub_rho_thetav)/rho_new
  qv    = (rho*qv + dt * tend_sub_rho_qv)/rho_new
  qcl =0.0
  qci =0.0
  qpl =0.0
  qpi =0.0
  theta = thetav/(1.+epsv*qv)
  ! gas law
  pres  = p00 * (rgas * rho_new*thetav/(p00*100.))**(1./(1.-rgas/cp))
  tabs = theta * (pres/p00)**(rgas/cp)
  tke   = (rho*tke+dt*tend_sub_rho_tke)/rho_new

  ! get pressure on interfaces from pressure on mass levels and hydro balance
  presi(1) = pres(1)*(1. - ggr/cp/thetav(1)*(p00/pres(1))**(rgas/cp) * (zi(1)-z(1)) )**(cp/rgas)
  do k=2,nz
    presi(k) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k)-z(k)) )**(cp/rgas)
  end do

  rho=rho_new

end if

end subroutine subsidence

