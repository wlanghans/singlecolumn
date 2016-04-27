subroutine get_state()

use vars
use params
use grid

implicit none

integer :: k

! get u, v
u=rhou/rho
v=rhov/rho

! get thermo state
thetav= rhothv/rho
qv    = rhoqv/rho
qcl =0.0
qci =0.0
qpl =0.0
qpi =0.0
theta = thetav/(1.+epsv*qv)
! gas law
pres  = p00 * (rgas * rhothv/(p00*100.))**(1./(1.-rgas/cp))

tabs = theta * (pres/p00)**(rgas/cp) 

tke   = rhotke/rho

! get pressure on interfaces from pressure on mass levels and hydro balance
presi(1) = pres(1)*(1. - ggr/cp/thetav(1)*(p00/pres(1))**(rgas/cp) * (zi(1)-z(1)) )**(cp/rgas)
do k=2,nz
presi(k) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k)-z(k)) )**(cp/rgas)
end do

end subroutine get_state
