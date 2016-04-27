subroutine get_dimlessdiff()

use vars
use params
use grid

implicit none

integer :: k

tkdimless=0.
do k=1,nzm

   if (tk(k).gt.0.0) then
     tkdimless(k)=0.5 * (adz(k)*dz)**2 / dt / (tk(k)+0.0001)
   else
     tkdimless(k)=-999.
   end if

end do

end subroutine get_dimlessdiff
