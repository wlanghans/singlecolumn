subroutine zero_stuff()

use vars
use params
use grid

implicit none

t0   = t
qt0  = qt
qp0  = qp
u0   = u
v0   = v
tke0 = tke

tend_sgs_qt      = 0.0
tend_sgs_qp      = 0.0
tend_sgs_t       = 0.0
tend_sgs_u       = 0.0
tend_sgs_v       = 0.0
tend_sgs_tke     = 0.0
tend_tke_buoy    = 0.0
tend_tke_shear   = 0.0 
tend_tke_diss    = 0.0

tend_force_qt    = 0.0
tend_force_t     = 0.0
tend_force_u     = 0.0
tend_force_v     = 0.0
tend_force_tke   = 0.0
tend_sub_qt      = 0.0
tend_sub_qp      = 0.0
tend_sub_t       = 0.0
tend_sub_u       = 0.0
tend_sub_v       = 0.0
tend_sub_tke     = 0.0
tend_rad_t       = 0.0


qt  = max(0.0,qt)
qp  = max(0.0,qp)
tke = max(0.0,tke)

tk          = 0.0
tkh         = 0.0
sumM        = 0.0
sumMt       = 0.0
sumMrt      = 0.0
sumMrp      = 0.0
sumMu       = 0.0
sumMv       = 0.0
sumMtke     = 0.0
tke_thvflx  = 0.0
radlwup     = 0.0
qcsgs_pdf   = 0.0
qcsgs_mf    = 0.0
cfrac_mf    = 0.0
cfrac_pdf   = 0.0
q1          = 0.0
qtgrad      = 0.0
thetaligrad = 0.0

Cm          = 0.0
Ctheta      = 0.0
Crv         = 0.0

end subroutine zero_stuff
