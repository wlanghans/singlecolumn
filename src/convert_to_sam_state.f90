subroutine convert_to_sam_state()

use vars
use params
use grid

implicit none

real :: oomqv(nzm)

! get mixing ratios and dry air density
oomqv = 1./(1.-qv)
rho = rho  * (1.-qv)

qv  =qv*oomqv
qcl =qcl*oomqv
qci =qci*oomqv
qpl =qpl*oomqv
qpi =qpi*oomqv
tke =tke *oomqv

t =  tabs & ! 
      + ggr/cp*z &
      - fac_cond * (qcl + qpl) &
      - fac_sub * (qci + qpi)



end subroutine convert_to_sam_state
