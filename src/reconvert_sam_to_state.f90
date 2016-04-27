subroutine reconvert_sam_to_state()

use vars
use params
use grid

implicit none

real :: oooprv(nzm)

! get mass fractions and air density
oooprv = 1./(1.+qv)

qv  =qv * oooprv
rho = rho  / (1.-qv)

qcl =qcl*oooprv
qci =qci*oooprv
qpl =qpl*oooprv
qpi =qpi*oooprv
tke =tke*oooprv

end subroutine reconvert_sam_to_state
