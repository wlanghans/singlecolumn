
subroutine get_def2()

use vars, only : def2, u, v
use grid
implicit none
integer :: k

do k=1,nzm
  if (k.gt.1 .and. k.lt.nzm) then
    def2(k) =   0.25*( (u(k+1)-u(k))/(z(k+1) - z(k)) + (u(k)-u(k-1))/(z(k) - z(k-1)) )**2 + &
                0.25*( (v(k+1)-v(k))/(z(k+1) - z(k)) + (v(k)-v(k-1))/(z(k) - z(k-1)) )**2
  elseif (k.eq.1) then
    !def2(k) = 0.25*( (u(k+1)-u(k))/(z(k+1) - z(k)) )**2 + &
    !          0.25*( (v(k+1)-v(k))/(z(k+1) - z(k)) )**2 
    def2(k) = ( (u(k+1)-u(k))/(z(k+1) - z(k)) )**2 + &
              ( (v(k+1)-v(k))/(z(k+1) - z(k)) )**2 
  else
    !def2(k) = 0.25*( (u(k)-u(k-1))/(z(k) - z(k-1)) )**2 + &
    !          0.25*( (v(k)-v(k-1))/(z(k) - z(k-1)) )**2
   
    def2(k) = ( (u(k)-u(k-1))/(z(k) - z(k-1)) )**2 + &
              ( (v(k)-v(k-1))/(z(k) - z(k-1)) )**2
  end if
end do 

end subroutine get_def2


