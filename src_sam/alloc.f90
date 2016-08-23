subroutine alloc

use vars
use grid
use params

implicit none

allocate(z(nzm),adz(nzm),zi(nz),adzw(nz))

allocate(t(nzm),qt(nzm),qp(nzm),u(nzm),v(nzm),tke(nzm))
allocate(t0(nzm),qt0(nzm),qp0(nzm),u0(nzm),v0(nzm),tke0(nzm))

allocate(theta(nzm), thetav(nzm), thetar(nzm), thetali(nzm), pres(nzm), tabs(nzm), tabs0(nzm), &
                    qv(nzm), qcl(nzm), qpl(nzm), qci(nzm), qpi(nzm), qn(nzm), &
                    tkh(nzm), tk(nzm), def2(nzm), rho(nzm),  &
                    rho_i(2:nzm), w(nz), ug(nzm), vg(nzm),dqtdt(nzm), dthetadt(nzm), presi(nz))

allocate(tend_sgs_qt(nzm), tend_sgs_t(nzm),  &
     tend_sgs_u(nzm),tend_sgs_v(nzm),tend_sgs_tke(nzm), tend_sgs_qp(nzm), &
     tend_tke_buoy(nzm), tend_tke_shear(nzm), tend_tke_diss(nzm), &
     tend_sub_t(nzm), tend_sub_u(nzm), tend_sub_v(nzm), tend_sub_tke(nzm), tend_sub_qt(nzm), &
     tend_force_qt(nzm),tend_force_u(nzm), tend_force_v(nzm), tend_force_tke(nzm),tend_sub_qp(nzm), &
     tend_force_t(nzm),tend_rad_t(nzm))


allocate(sumM(nz), sumMu(nz), sumMv(nz), sumMrt(nz), sumMrp(nz), sumMt(nz),sumMtke(nz), &
         buoy_sgs(nzm), smix(nzm))

allocate(sgs_t_flux(nz), sgs_qt_flux(nz), sgs_thv_flux(nz), qt_flux_ed(nz), qt_flux_mf(nz), &
        t_flux_ed(nz), t_flux_mf(nz), taux(nz), tauy(nz), sgs_tke_flux(nz), sgs_qp_flux(nz), &
        qp_flux_ed(nz), qp_flux_mf(nz))

allocate(qcsgs_mf(nz), tabs_mf(nz), qcsgs_pdf(nzm), qisgs_mf(nz), qtsgs_mf(nz), qisgs_pdf(nzm), &
                    cfrac_mf(nz), frac_mf(nz), cfrac_tot(nzm), cfrac_pdf(nzm), tke_thvflx(nzm),cthl(nzm), cqt(nzm),&
              varwrt1(nzm))
allocate(q1(nzm),sigmas(nzm),thetaligrad(nzm),qtgrad(nzm),radlwup(nz),radlwdn(nz), radswup(nz), radswdn(nz), radqrlw(nz)&
                 , radqrsw(nz))

allocate(UPM(nz,nup),UPW(nz,nup),UPT(nz,nup),UPQT(nz,nup),UPQCL(nz,nup),UPQCI(nz,nup),&
         UPA(nz,nup),UPU(nz,nup),UPV(nz,nup),UPTHV(nz,nup),UPTABS(nz,nup),ENT(nzm,nup),&
         DET(nzm,nup),UPTHD(nzm,nup),BUOY(nzm,nup),UPCF(nzm,nup))


end subroutine alloc
