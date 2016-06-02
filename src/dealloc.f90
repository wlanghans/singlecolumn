subroutine dealloc

use vars
use grid
use params

implicit none

deallocate(z,adz,zi,adzw)

deallocate(rhothv,rhoqv,rhou, rhov, rho, rhotke)
deallocate(u,v,t, theta, thetav, pres, tabs, &
                    qv, qcl, qpl, qci, qpi, tke, tkh, tk, def2, presi, & 
                    tkdimless, tabs0)

deallocate(tend_sgs_rho_qv, tend_sgs_rho, tend_sgs_rho_thetav, &
     tend_sgs_rho_u,tend_sgs_rho_v,tend_sgs_rho_tke, tend_sgs_rho_qt, &
     tend_rho_tke_buoy, tend_rho_tke_shear, tend_rho_tke_diss, & 
sumM, sumMrv, sumMqt, sumMu, sumMv, sumMtke, sumMthetav, sumMthetal, buoy_sgs, smix)

deallocate(sgs_sens_heat_flux, sgs_lat_heat_flux, sgs_thv_flux, taux, tauy, sgs_tke_flux)
deallocate(UPM,UPW,UPTHL,UPQT,UPQC,UPA,UPU,UPV,UPTHV,UPT,ENT,BUOY,DET,UPTHD)

end subroutine dealloc
