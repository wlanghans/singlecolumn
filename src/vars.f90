module vars

use grid
implicit none

! prognostic variables
real,allocatable,dimension(:), target ::  rhothv, rhoqv, rhou, rhov, rho, rhotke

! diagnostic variables
real,allocatable,dimension(:), target ::  u,v,t, theta, thetav, pres, tabs, tabs0, &
                    qv, qcl, qpl, qci, qpi, tke, tkh, tk, def2, tkdimless, rho_i, w, ug, vg, &
                    dqvdt, dthetadt
real,allocatable,dimension(:), target ::  presi

! tendencies
real, allocatable,dimension(:),target :: tend_sgs_rho_qv, tend_sgs_rho, tend_sgs_rho_thetav, &
     tend_sgs_rho_u,tend_sgs_rho_v,tend_sgs_rho_tke, tend_sgs_rho_qt, & 
     tend_rho_tke_buoy, tend_rho_tke_shear, tend_rho_tke_diss, tend_sub_rho, tend_sub_rho_qv, &
     tend_sub_rho_thetav, tend_sub_rho_u, tend_sub_rho_v, tend_sub_rho_tke, &
     tend_force_rho, tend_force_rho_qv, tend_force_rho_u, tend_force_rho_v, tend_force_rho_tke, &
     tend_force_rho_thetav, tend_rad_rho_thetav

! sgs variables
real vmag, ustar, r_s, thetav_s, theta_s
real, target :: pblh, wstar, Cthetav, Cm, Crv, Ctheta
real, allocatable,dimension(:), target :: sumM, sumMrv, sumMu, sumMv, sumMtke, sumMthetav, buoy_sgs, smix, sumMqt, sumMthetal
real,allocatable,dimension(:), target :: sgs_sens_heat_flux, sgs_lat_heat_flux, sgs_thv_flux, taux, tauy, sgs_tke_flux
real,allocatable,dimension(:), target :: qlsgs_mf, qlsgs_ed, cfrac_mf, cfrac_ed, tke_thvflx
real,allocatable,dimension(:), target :: thvflux_ed,thvflux_mf, q1, qtgrad, thetalgrad, radlwdn

real, parameter :: t00 = 300.   ! constant offset for sstxy 
real :: sst=0.                   ! perturbation from t00

! register functions:
real, external :: esatw,esati,dtesatw,dtesati,qsatw,qsati,dtqsatw,dtqsati
integer, external :: lenstr, bytes_in_rec

! plume properties
real, allocatable, dimension(:,:), target :: UPM,UPW,UPTHL,UPQT,UPQC,UPA,UPU,UPV,UPTHV,UPT,ENT,DET,BUOY,UPTHD



end module vars
