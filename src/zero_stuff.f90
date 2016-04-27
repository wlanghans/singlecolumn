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

rhoqv  = max(0.0,rhoqv)
rhotke = max(0.0,rhotke)

tk = 0.0
tkdimless = 0.0
tkh= 0.0
sumM = 0.0
sumMthetav=0.
sumMrv = 0.0
sumMu  = 0.
sumMv  = 0.
sumMtke= 0.

! set free slip and no transfer
Cm      = 0.0
Cthetav = 0.0
Ctheta  = 0.0
Crv     = 0.0

sgs_thv_flux=0.0
sgs_sens_heat_flux=0.0
sgs_lat_heat_flux=0.0
taux=0.0
tauy=0.0
sgs_tke_flux=0.0

end subroutine zero_stuff
