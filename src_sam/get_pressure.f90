subroutine get_pressure()
! computes pressure based on eos
! input: tabs, qv, thetav

use vars
use params
use grid

implicit none

integer :: k

! get pressure from gas law
pres = tabs * rho * rgas * (1. + qv * rv/rgas ) / 100.

! get pressure on interfaces from pressure on mass levels and hydro balance
presi(1) = pres(1)*(1. - ggr/cp/thetav(1)*(p00/pres(1))**(rgas/cp) * (zi(1)-z(1)) )**(cp/rgas)
do k=2,nz
presi(k) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k)-z(k)) )**(cp/rgas)
end do

end subroutine get_pressure
