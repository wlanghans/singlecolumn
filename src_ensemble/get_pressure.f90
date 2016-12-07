subroutine get_pressure()
! computes pressure based on eos
! input: tabs, qv, thetav

use vars
use params
use grid

implicit none

integer :: k


! get pressure on interfaces from pressure on mass levels and hydro balance
do k=1,nzm
 pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
 presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)
end do

end subroutine get_pressure
