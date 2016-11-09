subroutine sgs_tendencies()

use vars
use params
use grid

implicit none

! local variables

real,dimension(nzm) :: a, b, c, d
integer :: k

!+++++++++++++++++++++++++++++++++++++++
! non-precip water mixing ratio qt 
!+++++++++++++++++++++++++++++++++++++++
if (doneuman.or.sfc_flx_fxd) then
  !Neuman
  call get_abcd(rho,qt,sumMrt,Crv,tkh,r_s,a,b,c,d,.true.,.true.,sgs_qt_flux (1))
else
  !Dirichlet
  call get_abcd(rho,qt,sumMrt,Crv,tkh,r_s,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d,nzm)

if (.not.(doneuman.or.sfc_flx_fxd)) then
  ! Diagnose implicit surface flux in kinematic units
  sgs_qt_flux (1) = - vmag * Crv * (betam*qt(1) + betap*d(1) - r_s )
else
  ! Nothing to do
  sgs_qt_flux (1) = sgs_qt_flux (1) 
end if
! Diagnose implicit atm. sgs fluxes in kinematic units
! sgs_qt_flux (2:nzm) = 0.5*(rho(1:nzm-1) + rho(2:nzm)) * (  &  
sgs_qt_flux (2:nzm) = 1.0 * (  &  
              (-1.) /adzw(2:nzm)/dz *         0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qt(2:nzm)-qt(1:nzm-1)) ) &
                           +(sumMrt(2:nzm) - (betap*d(2:nzm) + betam*qt(2:nzm)) * sumM(2:nzm) ))
qt_flux_ed(1)     = sgs_qt_flux (1)
qt_flux_ed(nz)    = 0.0
qt_flux_ed(2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *           &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qt(2:nzm)-qt(1:nzm-1)) )
qt_flux_mf(1)     = 0.
qt_flux_mf(nz)    = 0.
qt_flux_mf(2:nzm) = (sumMrt(2:nzm) - (betap*d(2:nzm) + betam*qt(2:nzm)) * sumM(2:nzm) )

! compute qt tendency
tend_sgs_qt = (d-qt)/dt    

! update qt
if (dosequential) qt = d 

!+++++++++++++++++++++++++++++++++++++++
! precip water mixing ratio qt 
!+++++++++++++++++++++++++++++++++++++++
if (doneuman.or.sfc_flx_fxd) then
  !Neuman
  call get_abcd(rho,qp,sumMrp,Crv,tkh,0.0,a,b,c,d,.false.,.true.,0.0)
else
  !Dirichlet
  call get_abcd(rho,qp,sumMrp,Crv,tkh,0.0,a,b,c,d,.false.,.false.)
end if
call tridiag(a,b,c,d,nzm)

sgs_qp_flux(1)  = 0.
sgs_qp_flux(nz) = 0.
sgs_qp_flux (2:nzm) =  &  
              (-1.) /adzw(2:nzm)/dz *         0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qp(2:nzm)-qp(1:nzm-1)) ) ! &
!                           +(sumMrp(2:nzm) - (betap*d(2:nzm) + betam*qp(2:nzm)) * sumM(2:nzm) )
qp_flux_ed(1)     = 0.
qp_flux_ed(nz)    = 0.0
qp_flux_ed(2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *           &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qp(2:nzm)-qp(1:nzm-1)) )
qp_flux_mf(1)     = 0.
qp_flux_mf(nz)    = 0.
qp_flux_mf(2:nzm) = 0. !(sumMrp(2:nzm) - (betap*d(2:nzm) + betam*qp(2:nzm)) * sumM(2:nzm) ) 

! compute qt tendency
tend_sgs_qp = (d-qp)/dt    

! update qt
if (dosequential) qp = d 


!+++++++++++++++++++++++++++++++++++++++
! hli 
!+++++++++++++++++++++++++++++++++++++++

if (dotlflux) then
  sgs_t_flux (1) = sgs_t_flux (1)/cp
  t = (t -ggr * z)/cp * (p00/pres)**(rgas/cp)
  t_s = t_s/cp * (p00/presi(1))**(rgas/cp)
end if



if (doneuman.or.sfc_flx_fxd) then 
  !Neuman
  call get_abcd(rho,t,sumMt,Ctheta,tkh,t_s,a,b,c,d,.true.,.true.,sgs_t_flux (1))
else
  !Dirichlet
  call get_abcd(rho,t,sumMt,Ctheta,tkh,t_s,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d,nzm)


if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface hli flux 
  if (.not.dotlflux) then
    sgs_t_flux (1)       = - vmag * Ctheta * (betam*t(1) + betap*d(1) - t_s )
  else
    sgs_t_flux (1)       = - cp * vmag * Ctheta * (betam*t(1) + betap*d(1) - t_s )
  end if
else
  ! Nothing to do but conversion 
  if (dotlflux) sgs_t_flux (1) = cp * sgs_t_flux (1)
end if

! Diagnose implicit surface buoyancy flux 
 sgs_thv_flux(1) = sgs_t_flux(1)/cp * (1.+epsv* qv(1)) + epsv * theta(1) * sgs_qt_flux(1)

if (dotlflux) then 
  t_s = t_s*cp*(presi(1)/p00)**(rgas/cp)
end if

! Diagnose atm. sgs fluxes 
sgs_t_flux (2:nzm)       = 1. * ( & 
            (-1.)/adzw(2:nzm)/dz * 0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(t(2:nzm)-t(1:nzm-1)) ) &
        + sumMt(2:nzm) - (betap*d(2:nzm) + betam*t(2:nzm)) * sumM(2:nzm))
if (dotlflux) sgs_t_flux (2:nzm) = cp * sgs_t_flux (2:nzm) 

t_flux_ed(1)     = sgs_t_flux (1)
t_flux_ed(nz)    = 0.0
t_flux_ed(2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *           &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(t(2:nzm)-t(1:nzm-1)) )
if (dotlflux) t_flux_ed(2:nzm) = cp * t_flux_ed(2:nzm) 
t_flux_mf(1)     = 0.
t_flux_mf(nz)    = 0.
t_flux_mf(2:nzm) = (sumMt(2:nzm) - (betap*d(2:nzm) + betam*t(2:nzm)) * sumM(2:nzm) )
if (dotlflux) t_flux_mf(2:nzm) = cp* t_flux_mf(2:nzm) 


if (dotlflux) then
  tend_sgs_t = cp * (d - t)/dt 
else
  tend_sgs_t = (d - t)/dt
end if

! get buoyancy flux from PDF scheme
do k=2,nzm
  sgs_thv_flux(k) = sgs_t_flux (k)
  if ((qp(k-1)+qp(k)).ne.0.0) sgs_thv_flux(k)=sgs_thv_flux(k) &
               + (lcond*(qpl(k-1)+qpl(k))+lsub*(qpi(k-1)+qpi(k)))/(qp(k-1)+qp(k)) * sgs_qp_flux(k) 
  sgs_thv_flux(k) = 0.5*(cthl(k-1)+cthl(k))   * (2.*p00/(pres(k-1) + pres(k)))**(rgas/cp) / cp &
                    * sgs_thv_flux(k) + 0.5*(cqt(k-1)+cqt(k)) * sgs_qt_flux(k)
end do



! update t
if (dosequential) then 
  if (.not.dotlflux) then 
     t = d
  else
     t = d*cp + ggr * z
  end if
end if

!+++++++++++++++++++++++++++++++++++++++
! u 
!+++++++++++++++++++++++++++++++++++++++
if (doneuman.or.sfc_flx_fxd) then 
  ! Neuman
  call get_abcd(rho,u,sumMu,Cm,tk,0.,a,b,c,d,.false.,.true.,taux(1))
else
  ! Dirichlet
  call get_abcd(rho,u,sumMu,Cm,tk,0.,a,b,c,d,.false.,.false.)
end if
call tridiag(a,b,c,d,nzm)

if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface momentum flux
  taux(1) = - Cm * vmag * (betap*d(1)+betam*u(1))
else
 ! Nothing to do
  taux(1) = taux(1) 
end if
!diagnose atm. fluxes, NOTE: no contribution from mass fluxes
taux(2:nzm) = -1./adzw(2:nzm)/dz * &
                            0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(u(2:nzm)-u(1:nzm-1)) )
!tendency
tend_sgs_u = (d - u)/dt

! update u
if (dosequential) u=d

!+++++++++++++++++++++++++++++++++++++++
! v (input has to be air density rho and u)
!+++++++++++++++++++++++++++++++++++++++
if (doneuman.or.sfc_flx_fxd) then 
  ! Neuman
  call get_abcd(rho,v,sumMv,Cm,tk,0.,a,b,c,d,.false.,.true.,tauy(1))
else
  ! Dirichlet
  call get_abcd(rho,v,sumMv,Cm,tk,0.,a,b,c,d,.false.,.false.)
end if
call tridiag(a,b,c,d,nzm)

if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface momentum flux
  tauy(1) = - Cm * vmag * (betap*d(1)+betam*v(1))
else
  tauy(1) = tauy(1) 
end if
!diagnose atm. fluxes, NOTE: no contribution from mass fluxes
tauy(2:nzm) = -1.0/adzw(2:nzm)/dz * &
                            0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(v(2:nzm)-v(1:nzm-1)) )
!tendency
tend_sgs_v = (d - v)/dt

! update v
if (dosequential) v=d

!+++++++++++++++++++++++++++++++++++++++
! tke transport
!+++++++++++++++++++++++++++++++++++++++
if (progtke) then
  tke_s = 3.75*ustar**2 + 0.2*wstar**2
  if (dotkedirichlet) then
    call get_abcd(rho,tke,sumMtke,2.* tk(1)/adz(1)/dz/vmag,tk,tke_s,a,b,c,d,.false.,.false.)
  else
    sgs_tke_flux(1) = 0.0   ! lower BC for tke transfer (questionable)
    call get_abcd(rho,tke,sumMtke,0.,tk,tke_s,a,b,c,d,.false.,.true.,sgs_tke_flux(1))
  end if
  call tridiag(a,b,c,d,nzm)

  if (dotkedirichlet) sgs_tke_flux(1) = - 2.*tk(1)/adz(1)/dz * (betap*d(1)+betam*tke(1))
  !diagnose atm. fluxes
  sgs_tke_flux(2:nzm) =    (-1.)/adzw(2:nzm)/dz * 0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(tke(2:nzm)-tke(1:nzm-1)) ) 

  ! tendency
  tend_sgs_tke = (d-tke)/dt

  ! update tke
  if (dosequential) tke = d
end if

end subroutine sgs_tendencies

