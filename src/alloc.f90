subroutine alloc

use vars
use grid
use params

implicit none

allocate(z(nzm),adz(nzm),zi(nz),adzw(nz))

allocate(rhothv(nzm),rhoqv(nzm),rhou(nzm), rhov(nzm), rho(nzm), rhotke(nzm))
allocate(u(nzm),v(nzm),t(nzm), theta(nzm), thetav(nzm), pres(nzm), tabs(nzm), &
                    qv(nzm), qcl(nzm), qpl(nzm), qci(nzm), qpi(nzm), tke(nzm), &
                    tkh(nzm), tk(nzm), def2(nzm), presi(nz), tkdimless(nzm), tabs0(nzm), &
                    rho_i(2:nzm), w(nz), ug(nzm), vg(nzm))

allocate(tend_sgs_rho_qv(nzm), tend_sgs_rho(nzm), tend_sgs_rho_thetav(nzm), &
     tend_sgs_rho_u(nzm),tend_sgs_rho_v(nzm),tend_sgs_rho_tke(nzm), tend_sgs_rho_qt(nzm), &
     tend_rho_tke_buoy(nzm), tend_rho_tke_shear(nzm), tend_rho_tke_diss(nzm), & 
     sumM(nz), sumMrv(nz), sumMqt(nz),sumMthetal(nz),sumMu(nz), sumMv(nz), sumMtke(nz), &
     sumMthetav(nz), buoy_sgs(nzm), smix(nzm),                                          &
     tend_sub_rho(nzm), tend_sub_rho_qv(nzm), tend_rad_rho_thetav(nzm),&
     tend_sub_rho_thetav(nzm), tend_sub_rho_u(nzm), tend_sub_rho_v(nzm), tend_sub_rho_tke(nzm), & 
     tend_force_rho(nzm), tend_force_rho_qv(nzm), tend_force_rho_u(nzm), tend_force_rho_v(nzm), tend_force_rho_tke(nzm), &
     tend_force_rho_thetav(nzm), dqvdt(nzm), dthetadt(nzm))

allocate(sgs_sens_heat_flux(nzm), sgs_lat_heat_flux(nzm), sgs_thv_flux(nzm), taux(nzm), tauy(nzm), sgs_tke_flux(nzm), radlwup(nz))

allocate(UPM(nz,nup),UPW(nz,nup),UPTHL(nz,nup),UPQT(nz,nup),UPQC(nz,nup),UPA(nz,nup),UPU(nz,nup),UPV(nz,nup),&
          UPTHV(nz,nup),UPT(nz,nup))
allocate(qlsgs_mf(nz), qlsgs_ed(nzm), cfrac_mf(nz), cfrac_ed(nzm), tke_thvflx(nzm))
allocate(ENT(nzm,nup),BUOY(nzm,nup),DET(nzm,nup),UPTHD(nzm,nup))
allocate(thvflux_ed(nz),thvflux_mf(nz))
allocate(q1(nzm),thetalgrad(nzm),qtgrad(nzm))


end subroutine alloc
