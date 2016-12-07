subroutine dealloc

use vars
use grid
use params

implicit none

deallocate(z,adz,zi,adzw)

deallocate(t,qt,qp,u,v,tke)
deallocate(t0,qt0,qp0,u0,v0,tke0)
deallocate(theta, thetav, thetar, thetali, pres, tabs, tabs0, &
                    qv, qcl, qpl, qci, qpi, qn, &
                    tkh, tk, def2, rho, &
                    rho_i, w, ug, vg,dqtdt, dthetadt, presi)

deallocate(tend_sgs_qt, tend_sgs_t,  &
     tend_sgs_u,tend_sgs_v,tend_sgs_tke, tend_sgs_qp, &
     tend_tke_buoy, tend_tke_shear, tend_tke_diss, &
     tend_sub_t, tend_sub_u, tend_sub_v, tend_sub_tke,tend_sub_qt,tend_sub_qp, &
     tend_force_qt,tend_force_u, tend_force_v, tend_force_tke, &
     tend_force_t,tend_rad_t)


deallocate(sumM, sumMu, sumMv, sumMrt,sumMrp,sumMt,sumMtke,sumMthv,sumDEF2, &
         buoy_sgs, smix)

deallocate(sgs_t_flux, sgs_qt_flux, sgs_thv_flux, qt_flux_mf, qt_flux_ed, t_flux_ed, t_flux_mf, taux, tauy, sgs_tke_flux,&
           sgs_qp_flux, qp_flux_ed, qp_flux_mf)

deallocate(qcsgs_mf, tabs_mf, tke_mf, qcsgs_pdf,qisgs_mf, qtsgs_mf, qisgs_pdf, cfrac_mf, frac_mf, cfrac_tot, &
     cfrac_pdf, tke_thvflx, cthl, cqt, varwrt1)
deallocate(q1,sigmas,thetaligrad,qtgrad,radlwup,radlwdn, radswup, radswdn, radqrlw, radqrsw)

deallocate(UPM,UPW,UPT,UPQT,UPQCL,UPQCI,&
         UPA,UPU,UPV,UPTHV,UPTABS,ENT,&
         BUOY,DET,UPTHD,UPCF)


end subroutine dealloc
