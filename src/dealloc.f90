subroutine dealloc

use vars
use grid
use params

implicit none

deallocate(z,adz,zi,adzw)

deallocate(rhothv,rhoqv,rhou, rhov, rho, rhotke)
deallocate(u,v,w,t, theta, thetav, pres, tabs, &
                    qv, qcl, qpl, qci, qpi, tke, tkh, tk, def2, presi, & 
                    tkdimless, tabs0, rho_i, ug, vg)

deallocate(tend_sgs_rho_qv, tend_sgs_rho, tend_sgs_rho_thetav, &
     tend_sgs_rho_u,tend_sgs_rho_v,tend_sgs_rho_tke, tend_sgs_rho_qt, &
     tend_rho_tke_buoy, tend_rho_tke_shear, tend_rho_tke_diss, & 
     sumM, sumMrv, sumMqt, sumMu, sumMv, sumMtke, sumMthetav, sumMthetal, buoy_sgs, smix,&
     tend_sub_rho, tend_sub_rho_qv,tend_rad_rho_thetav, &
     tend_sub_rho_thetav, tend_sub_rho_u, tend_sub_rho_v, tend_sub_rho_tke, &
     tend_force_rho, tend_force_rho_qv, tend_force_rho_u, tend_force_rho_v, tend_force_rho_tke, &
     tend_force_rho_thetav,dqvdt,dthetadt,radlwdn )

deallocate(sgs_sens_heat_flux, sgs_lat_heat_flux, sgs_thv_flux, taux, tauy, sgs_tke_flux)
deallocate(UPM,UPW,UPTHL,UPQT,UPQC,UPA,UPU,UPV,UPTHV,UPT,ENT,BUOY,DET,UPTHD)
deallocate(qlsgs_mf, qlsgs_ed, cfrac_mf, cfrac_ed, tke_thvflx)
deallocate(thvflux_ed,thvflux_mf,q1,qtgrad,thetalgrad)

end subroutine dealloc
