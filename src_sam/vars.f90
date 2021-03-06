module vars

use grid
implicit none

! prognostic variables
real,allocatable,dimension(:), target ::  t, qt, qp, u, v, tke, &
                                          t0, qt0, qp0, u0, v0, tke0

! diagnostic variables
real,allocatable,dimension(:), target ::  theta, thetav, thetar, thetali, pres, tabs, tabs0, &
                    qv, qcl, qpl, qci, qpi, qn, tkh, tk, def2, rho_i, w, ug, vg, &
                    dqtdt, dthetadt, rho, presi

! tendencies
real, allocatable,dimension(:),target :: tend_sgs_qt, tend_sgs_t, &
     tend_sgs_u,tend_sgs_v,tend_sgs_tke, tend_sgs_qp, & 
     tend_tke_buoy, tend_tke_shear, tend_tke_diss, tend_sub_qt, tend_sub_qp,&
     tend_sub_t, tend_sub_u, tend_sub_v, tend_sub_tke, &
     tend_force_qt, tend_force_u, tend_force_v, tend_force_tke, &
     tend_force_t, tend_rad_t

! sgs variables
real(8) vmag, r_s, t_s, theta_s, feddy
real, target :: pblh, wstar, Cm, Crv, Ctheta, lwp, tke_s, ustar
real, allocatable,dimension(:), target :: sumM, sumMu, sumMv, sumMtke, sumMt, sumMrt, sumMrp, buoy_sgs, smix, sumMthv,sumDEF2
real,allocatable,dimension(:), target :: sgs_t_flux, sgs_qt_flux, sgs_qp_flux, sgs_thv_flux, qt_flux_ed, qt_flux_mf, &
                                         t_flux_ed, t_flux_mf, taux, tauy, sgs_tke_flux, qp_flux_ed, qp_flux_mf
real,allocatable,dimension(:), target :: qcsgs_pdf, qcsgs_mf, tabs_mf, tke_mf, qisgs_pdf, qisgs_mf, qtsgs_mf, &
                                            cfrac_mf, frac_mf, cfrac_pdf, cfrac_tot, tke_thvflx, cthl, cqt, varwrt1
real,allocatable,dimension(:), target :: q1, sigmas, qtgrad, thetaligrad, radlwup, radlwdn, radswup, radswdn, radqrlw, radqrsw

real, parameter :: t00 = 300.   ! constant offset for sstxy 
real :: sst=0.                   ! perturbation from t00

! register functions:
real, external :: esatw,esati,dtesatw,dtesati,qsatw,qsati,dtqsatw,dtqsati
integer, external :: lenstr, bytes_in_rec

! plume properties
real, allocatable, dimension(:,:), target :: UPM,UPW,UPT,UPQT,UPQCL,UPQCI,UPA,UPU,UPV,UPTHV,UPTABS,ENT,DET,BUOY,UPTHD,UPCF



end module vars
