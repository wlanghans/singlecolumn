subroutine sgs_tendencies()

use vars
use params
use grid

implicit none

! local variables

real,dimension(nzm) :: a, b, c, d

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
call tridiag(a,b,c,d)

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

! compute qt tendency
tend_sgs_qt = (d-qt)/dt    

! update qt
if (dosequential) qt = d 

!+++++++++++++++++++++++++++++++++++++++
! precip water mixing ratio qt 
!+++++++++++++++++++++++++++++++++++++++
if (doneuman.or.sfc_flx_fxd) then
  !Neuman
  call get_abcd(rho,qp,sumMrp,Crv,tkh,0.0,a,b,c,d,.true.,.true.,0.0)
else
  !Dirichlet
  call get_abcd(rho,qp,sumMrp,Crv,tkh,0.0,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d)

! compute qt tendency
tend_sgs_qp = (d-qp)/dt    

! update qt
if (dosequential) qp = d 


!+++++++++++++++++++++++++++++++++++++++
! hli 
!+++++++++++++++++++++++++++++++++++++++
if (doneuman.or.sfc_flx_fxd) then 
  !Neuman
  call get_abcd(rho,t,sumMt,Ctheta,tkh,t_s,a,b,c,d,.true.,.true.,sgs_t_flux (1))
else
  !Dirichlet
  call get_abcd(rho,t,sumMt,Ctheta,tkh,t_s,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d)

if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface hli flux 
  sgs_t_flux (1)       = - vmag * Ctheta * (betam*t(1) + betap*d(1) - t_s )
! Diagnose implicit surface buoyancy flux 
  sgs_thv_flux(1) = sgs_t_flux(1)/cp * (1.+epsv* qv(1)) + epsv * theta(1) * sgs_qt_flux(1)
else
  ! Nothing to do
  sgs_t_flux (1)    =  sgs_t_flux (1)  
  sgs_thv_flux (1)  = sgs_thv_flux (1)  
end if
! Diagnose atm. sgs fluxes 
sgs_t_flux (2:nzm)       = 1. * ( & 
            (-1.)/adzw(2:nzm)/dz * 0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(t(2:nzm)-t(1:nzm-1)) ) &
        + sumMt(2:nzm) - (betap*d(2:nzm) + betam*t(2:nzm)) * sumM(2:nzm))

sgs_thv_flux(2:nzm) = 0.5*(tke_thvflx(1:nzm-1) + tke_thvflx(2:nzm))

tend_sgs_t = (d - t)/dt

! update t
if (dosequential) t = d

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
call tridiag(a,b,c,d)

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
call tridiag(a,b,c,d)

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
  if (dotkedirichlet) then
    call get_abcd(rho,tke,sumMtke,2.* tk(1)/adz(1)/dz/vmag,tk,0.,a,b,c,d,.false.,.false.)
  else
    sgs_tke_flux(1) = 0.0   ! lower BC for tke transfer (questionable)
    call get_abcd(rho,tke,sumMtke,0.,tk,0.,a,b,c,d,.false.,.true.,sgs_tke_flux(1))
  end if
  call tridiag(a,b,c,d)

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

