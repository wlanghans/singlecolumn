subroutine setgrid

use grid
use params

implicit none

integer k

if (doconstdz) then
  z(1) = 0.5*dz
  do k=2,nzm
     z(k)=z(k-1) + dz
  end do
else

! read in z 
write(*,*) 'doconstdz=.false. not implemented yet'
write(*,*) 'STOPPING'

stop

end if

if(.not.doconstdz) dz = 0.5*(z(1)+z(2))

do k=2,nzm
   adzw(k) = (z(k)-z(k-1))/dz
end do
adzw(1) = 1.
adzw(nz) = adzw(nzm)

adz(1) = 1.
do k=2,nzm-1
   adz(k) = 0.5*(z(k+1)-z(k-1))/dz
end do
adz(nzm) = adzw(nzm)
zi(1) = 0.
do k=2,nz
   zi(k) = zi(k-1) + adz(k-1)*dz
end do

end subroutine setgrid
