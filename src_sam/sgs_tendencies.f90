subroutine sgs_tendencies()

use vars
use params
use grid

implicit none

! local variables

real,dimension(nzm) :: a, b, c, d, tkhm, tkm, frac_mf_avg, qte, te, tp
integer :: k

!+++++++++++++++++++++++++++++++++++++++
! non-precip water mixing ratio qt 
!+++++++++++++++++++++++++++++++++++++++

! multiply Km by area fraction of environment
frac_mf_avg(1:nzm) = 0.5*(frac_mf(1:nzm)+frac_mf(2:nz))
if (doedmf.and..not.donoscale) then
  tkm(1:nzm)=(1.0-frac_mf_avg(1:nzm)) * tk(1:nzm) 
  tkhm(1:nzm)=(1.0-frac_mf_avg(1:nzm)) * tkh(1:nzm) 
else
  tkm=tk
  tkhm=tkh
end if

!if doenvedmf, then the EDMF equation has the environment values (not domain means)


if (doedmf.and.(.not.donoscale).and.doenvedmf) then
  do k=1,nzm
  if (frac_mf_avg(k).gt.0.0) then
    if (frac_mf(k+1).gt.0.0) then
      qte(k) = (qt(k) - frac_mf_avg(k)*0.5*(qtsgs_mf(k+1)+qtsgs_mf(k)))/(1.-frac_mf_avg(k))
    else
      qte(k) = (qt(k) - frac_mf_avg(k)*qtsgs_mf(k))/(1.-frac_mf_avg(k))
    end if
  else
    qte(k) = qt(k)
  end if
  end do
else
  qte = qt
end if

if (doneuman.or.sfc_flx_fxd) then
  !Neuman
  call get_abcd(rho,qte,sumMrt,Crv,tkhm,r_s,a,b,c,d,.true.,.true.,sgs_qt_flux (1))
else
  !Dirichlet
  call get_abcd(rho,qte,sumMrt,Crv,tkhm,r_s,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d,nzm)

if (.not.(doneuman.or.sfc_flx_fxd)) then
  ! Diagnose implicit surface flux in Wm2 for output
  sgs_qt_flux (1) = - (2.* (0.5*(rho(1)+rho(2) )) - 0.5*(rho(2)+rho(3))) * lcond * vmag * Crv * (betam*qte(1) + betap*d(1) - r_s )
else
  ! Nothing to do
  sgs_qt_flux (1) = (2.* (0.5*(rho(1)+rho(2) )) - 0.5*(rho(2)+rho(3))) * lcond * sgs_qt_flux (1) 
end if
  ! Diagnose implicit atm flux in Wm2 for output
sgs_qt_flux (2:nzm) =0.5*(rho(1:nzm-1)+rho(2:nzm)) * lcond  * (  &  
              (-1.) /adzw(2:nzm)/dz *         0.5*(tkhm(1:nzm-1) + tkhm(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qte(2:nzm)-qte(1:nzm-1)) ) &
                           +(sumMrt(2:nzm) - (betap*d(2:nzm) + betam*qte(2:nzm)) * sumM(2:nzm) ))
! and components in kinematic units
qt_flux_ed(1)     = sgs_qt_flux (1) / lcond/ (2.* (0.5*(rho(1)+rho(2) )) - 0.5*(rho(2)+rho(3)))
qt_flux_ed(nz)    = 0.0
qt_flux_ed(2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tkhm(1:nzm-1) + tkhm(2:nzm)) *           &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qte(2:nzm)-qte(1:nzm-1)) )
qt_flux_mf(1)     = 0.
qt_flux_mf(nz)    = 0.
qt_flux_mf(2:nzm) = (sumMrt(2:nzm) - (betap*d(2:nzm) + betam*qte(2:nzm)) * sumM(2:nzm) )

! compute qt tendency
tend_sgs_qt = (d-qte)/dt    

! update qt
if (dosequential) then 
  if (doedmf.and.(.not.donoscale).and.doenvedmf) then
    qt = qt + dt * tend_sgs_qt  ! tendency for environment and mean are the same
  else
    qt = d 
  end if
end if

!+++++++++++++++++++++++++++++++++++++++
! precip water mixing ratio qp 
!+++++++++++++++++++++++++++++++++++++++

if (doneuman.or.sfc_flx_fxd) then
  !Neuman
  call get_abcd(rho,qp,sumMrp,Crv,tkhm,0.0,a,b,c,d,.false.,.true.,0.0)
else
  !Dirichlet
  call get_abcd(rho,qp,sumMrp,Crv,tkhm,0.0,a,b,c,d,.false.,.false.)
end if
call tridiag(a,b,c,d,nzm)

  ! Diagnose implicit atm flux in Wm2 for output
sgs_qp_flux(1)  = 0.
sgs_qp_flux(nz) = 0.
sgs_qp_flux (2:nzm) =  0.5*(rho(1:nzm-1)+rho(2:nzm)) * lcond  * &  
              (-1.) /adzw(2:nzm)/dz *         0.5*(tkhm(1:nzm-1) + tkhm(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qp(2:nzm)-qp(1:nzm-1)) ) ! &
! and components in kinematic units
qp_flux_ed(1)     = 0.
qp_flux_ed(nz)    = 0.0
qp_flux_ed(2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tkhm(1:nzm-1) + tkhm(2:nzm)) *           &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qp(2:nzm)-qp(1:nzm-1)) )
qp_flux_mf(1)     = 0.
qp_flux_mf(nz)    = 0.
qp_flux_mf(2:nzm) = 0. !(sumMrp(2:nzm) - (betap*d(2:nzm) + betam*qp(2:nzm)) * sumM(2:nzm) ) 

! compute qt tendency
tend_sgs_qp = (d-qp)/dt    

! update qp
if (dosequential) qp = d 


!+++++++++++++++++++++++++++++++++++++++
! hli 
!+++++++++++++++++++++++++++++++++++++++

if (dotlflux) then
  sgs_t_flux (1) = sgs_t_flux (1)/cp
  t = (t -ggr * z)/cp * (p00/pres)**(rgas/cp)
  t_s = t_s/cp * (p00/presi(1))**(rgas/cp)

  if (doedmf.and.(.not.donoscale).and.doenvedmf) then

  do k=1,nzm
  if (frac_mf_avg(k).gt.0.0) then
  if (frac_mf(k+1).gt.0.0) then
    tp(k) = 0.5*(cp*tabs_mf(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k) +&
              cp*tabs_mf(k+1)  - lcond * qcsgs_mf(k+1) - lsub * qisgs_mf(k+1))/cp * (p00/pres(k))**(rgas/cp)
  else
    tp(k) = (cp*tabs_mf(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k))/cp * (p00/pres(k))**(rgas/cp)
  end if
  end if
  end do

  end if


else

  if (doedmf.and.(.not.donoscale).and.doenvedmf) then

  do k=1,nzm
  if (frac_mf_avg(k).gt.0.0) then
  if (frac_mf(k+1).gt.0.0) then
    tp(k) = 0.5*(cp*tabs_mf(k) + ggr*zi(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k) +&
              cp*tabs_mf(k+1) + ggr*zi(k+1) - lcond * qcsgs_mf(k+1) - lsub * qisgs_mf(k+1))
  else
    tp(k) = cp*tabs_mf(k) + ggr*z(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k)
  end if
  end if
  end do

  end if

end if


if (doedmf.and.(.not.donoscale).and.doenvedmf) then
do k=1,nzm
  if (frac_mf_avg(k).gt.0.0) then
    te(k) = (t(k)-frac_mf_avg(k)*tp(k))/(1.-frac_mf_avg(k))
  else
    te(k) = t(k)
  end if
end do
else
  te=t
end if

if (doneuman.or.sfc_flx_fxd) then 
  !Neuman
  call get_abcd(rho,te,sumMt,Ctheta,tkhm,t_s,a,b,c,d,.true.,.true.,sgs_t_flux (1))
else
  !Dirichlet
  call get_abcd(rho,te,sumMt,Ctheta,tkhm,t_s,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d,nzm)


if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface flux in W/m2 for output
  if (.not.dotlflux) then
    sgs_t_flux (1)       = - (2.* (0.5*(rho(1)+rho(2) )) - 0.5*(rho(2)+rho(3))) * &
                            vmag * Ctheta * (betam*te(1) + betap*d(1) - t_s )
  else
    sgs_t_flux (1)       = - (2.* (0.5*(rho(1)+rho(2) )) - 0.5*(rho(2)+rho(3))) * &
                            cp * vmag * Ctheta * (betam*te(1) + betap*d(1) - t_s )
  end if
else
  sgs_t_flux (1) = (2.* (0.5*(rho(1)+rho(2) )) - 0.5*(rho(2)+rho(3))) * sgs_t_flux (1)
  if (dotlflux) sgs_t_flux (1) = cp * sgs_t_flux (1)
end if

! Diagnose implicit surface buoyancy flux 
 sgs_thv_flux(1) = sgs_t_flux(1) * (1.+epsv* qv(1)) + &
 epsv * theta(1) * cp * sgs_qt_flux(1)/lcond

if (dotlflux) then 
  t_s = t_s*cp*(presi(1)/p00)**(rgas/cp)
end if

! Diagnose atm. sgs fluxes in Wm2 for output
sgs_t_flux (2:nzm)       = 0.5*(rho(1:nzm-1)+rho(2:nzm)) * ( & 
            (-1.)/adzw(2:nzm)/dz * 0.5*(tkhm(1:nzm-1) + tkhm(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(te(2:nzm)-te(1:nzm-1)) ) &
        + sumMt(2:nzm) - (betap*d(2:nzm) + betam*te(2:nzm)) * sumM(2:nzm))
if (dotlflux) sgs_t_flux (2:nzm) = cp * sgs_t_flux (2:nzm) 

! needs to be in kinematic units since used later
t_flux_ed(1)     = sgs_t_flux (1) / (2.* (0.5*(rho(1)+rho(2) )) - 0.5*(rho(2)+rho(3))) / cp
t_flux_ed(nz)    = 0.0
t_flux_ed(2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tkhm(1:nzm-1) + tkhm(2:nzm)) *           &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(te(2:nzm)-te(1:nzm-1)) ) / cp
if (dotlflux) t_flux_ed(2:nzm) = cp * t_flux_ed(2:nzm) 
t_flux_mf(1)     = 0.
t_flux_mf(nz)    = 0.
t_flux_mf(2:nzm) = (sumMt(2:nzm) - (betap*d(2:nzm) + betam*te(2:nzm)) * sumM(2:nzm) )/cp
if (dotlflux) t_flux_mf(2:nzm) = cp* t_flux_mf(2:nzm) 


if (dotlflux) then
  tend_sgs_t = cp * (d - te)/dt 
else
  tend_sgs_t = (d - te)/dt
end if

! get buoyancy flux from PDF scheme in Wm2 for output
do k=2,nzm
  sgs_thv_flux(k) = sgs_t_flux (k)
  if ((qp(k-1)+qp(k)).ne.0.0) sgs_thv_flux(k)=sgs_thv_flux(k) &
               + (lcond*(qpl(k-1)+qpl(k))+lsub*(qpi(k-1)+qpi(k)))/(qp(k-1)+qp(k)) * &
                      sgs_qp_flux(k)/lcond * (2.*p00/(pres(k-1) + pres(k)))**(rgas/cp)
  sgs_thv_flux(k) =  0.5*(cthl(k-1)+cthl(k)) * sgs_thv_flux(k) &
                         + cp * 0.5*(cqt(k-1)+cqt(k)) * sgs_qt_flux(k)/lcond
end do


! update t
if (dosequential) then 

  if (doedmf.and.(.not.donoscale).and.doenvedmf) then
       t = t + dt*tend_sgs_t ! environment and mean tendency are the same
  else
   if (.not.dotlflux) then 
       t = d
    else
       t = d*cp + ggr * z
    end if
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
  tke_s = (3.75*ustar**2 + 0.2*wstar**2)
  if (dotkedirichlet) then
    call get_abcd(rho,tke,sumMtke,2.* tk(1)/adz(1)/dz/vmag,tk,tke_s,a,b,c,d,.false.,.false.)
  else
    sgs_tke_flux(1) = 0.0   ! lower BC for tke transfer (questionable)
    call get_abcd(rho,tke,sumMtke,0.,tk,tke_s,a,b,c,d,.false.,.true.,sgs_tke_flux(1))
  end if
  call tridiag(a,b,c,d,nzm)

  if (dotkedirichlet) sgs_tke_flux(1) = - 2.*tk(1)/adz(1)/dz * (betap*d(1)+betam*tke(1)-tke_s)
  !diagnose atm. fluxes
  sgs_tke_flux(2:nzm) =    (-1.)/adzw(2:nzm)/dz * 0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(tke(2:nzm)-tke(1:nzm-1)) ) 

  ! tendency
  tend_sgs_tke = (d-tke)/dt

  ! update tke
  if (dosequential) tke = d
end if

end subroutine sgs_tendencies

