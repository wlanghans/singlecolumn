subroutine zero_stuff()

use vars
use params
use grid

implicit none

tend_sgs_rho         = 0.0
tend_sgs_rho_qv      = 0.0
tend_sgs_rho_thetav  = 0.0
tend_sgs_rho_u       = 0.0
tend_sgs_rho_v       = 0.0
tend_sgs_rho_tke     = 0.0
tend_rho_tke_buoy    = 0.0
tend_rho_tke_shear   = 0.0 
tend_rho_tke_diss    = 0.0

tend_force_rho        = 0.0
tend_force_rho_qv     = 0.0
tend_force_rho_thetav = 0.0
tend_force_rho_u      = 0.0
tend_force_rho_v      = 0.0
tend_force_rho_tke    = 0.0
tend_sub_rho          = 0.0
tend_sub_rho_qv       = 0.0
tend_sub_rho_thetav   = 0.0
tend_sub_rho_u        = 0.0
tend_sub_rho_v        = 0.0
tend_sub_rho_tke      = 0.0
tend_rad_rho_thetav = 0.0


rhoqv  = max(0.0,rhoqv)
rhotke = max(0.0,rhotke)

tk = 0.0
tkdimless = 0.0
tkh= 0.0
sumM = 0.0
sumMthetav=0.
sumMrv = 0.0
sumMqt = 0.0
sumMu  = 0.
sumMv  = 0.
sumMtke= 0.
tke_thvflx=0.
radlwup = 0.
qlsgs_ed = 0.
qlsgs_mf = 0.
thvflux_ed = 0.
thvflux_mf = 0.
q1 = 0.
qtgrad=0.
thetalgrad=0.

! set free slip and no transfer
Cm      = 0.0
Cthetav = 0.0
Ctheta  = 0.0
Crv     = 0.0

end subroutine zero_stuff
