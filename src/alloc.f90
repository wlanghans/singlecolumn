subroutine alloc

use vars
use grid
use params

implicit none

allocate(z(nzm),adz(nzm),zi(nz),adzw(nz))

allocate(rhothv(nzm),rhoqv(nzm),rhou(nzm), rhov(nzm), rho(nzm), rhotke(nzm))
allocate(u(nzm),v(nzm),t(nzm), theta(nzm), thetav(nzm), pres(nzm), tabs(nzm), &
                    qv(nzm), qcl(nzm), qpl(nzm), qci(nzm), qpi(nzm), tke(nzm), &
                    tkh(nzm), tk(nzm), def2(nzm), presi(nz), tkdimless(nzm), tabs0(nzm))

allocate(tend_sgs_rho_qv(nzm), tend_sgs_rho(nzm), tend_sgs_rho_thetav(nzm), &
     tend_sgs_rho_u(nzm),tend_sgs_rho_v(nzm),tend_sgs_rho_tke(nzm), tend_sgs_rho_qt(nzm), &
     tend_rho_tke_buoy(nzm), tend_rho_tke_shear(nzm), tend_rho_tke_diss(nzm), & 
sumM(nz), sumMrv(nz), sumMqt(nz),sumMthetal(nz),sumMu(nz), sumMv(nz), sumMtke(nz), sumMthetav(nz), buoy_sgs(nzm), smix(nzm))

allocate(sgs_sens_heat_flux(nzm), sgs_lat_heat_flux(nzm), sgs_thv_flux(nzm), taux(nzm), tauy(nzm), sgs_tke_flux(nzm))

allocate(UPM(nz,nup),UPW(nz,nup),UPTHL(nz,nup),UPQT(nz,nup),UPQC(nz,nup),UPA(nz,nup),UPU(nz,nup),UPV(nz,nup),&
          UPTHV(nz,nup),UPT(nz,nup))
allocate(ENT(nzm,nup),BUOY(nzm,nup),DET(nzm,nup),UPTHD(nzm,nup))


end subroutine alloc
