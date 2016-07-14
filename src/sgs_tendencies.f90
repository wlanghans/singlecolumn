subroutine sgs_tendencies()

use vars
use params
use grid

implicit none

! local variables

real,dimension(nzm) :: rhoold, a, b, c, d, qvspec

!qv holds mixing ratio
!rho holds rhod

! get actual density
rhoold = rho * (1.+qv)
! get specific humidity
qvspec  = qv/(1.+qv)

! water vapor mixing ratio rv (input has to be dry air density rho_d and vapor mixing ratio rv)
if (doneuman.or.sfc_flx_fxd) then
  !Neuman
  call get_abcd(rho,qv,sumMrv,Crv,tkh,r_s,a,b,c,d,.true.,.true.,sgs_lat_heat_flux (1))
else
  !Dirichlet
  call get_abcd(rho,qv,sumMrv,Crv,tkh,r_s,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d)

if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface flux in mass flux units
  !sgs_lat_heat_flux (1) = - vmag *rho(1)* Crv * (betam*qv(1) + betap*d(1) - r_s )
  sgs_lat_heat_flux (1) = - vmag * Crv * (betam*qv(1) + betap*d(1) - r_s )
else
  !sgs_lat_heat_flux (1) = sgs_lat_heat_flux (1) *rho(1) 
  sgs_lat_heat_flux (1) = sgs_lat_heat_flux (1) 
end if
! Diagnose implicit atm. sgs fluxes in mass flux units
!sgs_lat_heat_flux (2:nzm) = 0.5*(rho(1:nzm-1) + rho(2:nzm)) * (  &  
sgs_lat_heat_flux (2:nzm) = 1.0 * (  &  
              (-1.) /adzw(2:nzm)/dz *         0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qv(2:nzm)-qv(1:nzm-1)) ) &
                           +(sumMrv(2:nzm) - (betap*d(2:nzm) + betam*qv(2:nzm)) * sumM(2:nzm) ))

! compute rho and rhov tendency
tend_sgs_rho_qv = rho * (d-qv)/dt    ! =rho_d * (rvnew - rvold) /dt
tend_sgs_rho    = tend_sgs_rho_qv    ! and this gives the mass tendency

! update qv
! d is still mixing ratio
qv = d 

! WL ADD DIFFUSION OF ADDITIONAL WATER VARIABLE, to zero for now 
tend_sgs_rho_qt = 0.

! thetav (input has to be air density rho and thetav)
if (doneuman.or.sfc_flx_fxd) then 
  !Neuman
  call get_abcd(rhoold,thetav,sumMthetav,Cthetav,tkh,thetav_s,a,b,c,d,.true.,.true.,sgs_thv_flux (1))
else
  !Dirichlet
  call get_abcd(rhoold,thetav,sumMthetav,Cthetav,tkh,thetav_s,a,b,c,d,.true.,.false.)
end if
call tridiag(a,b,c,d)

if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface virt. pot. temp. flux 
  !sgs_thv_flux (1)       = - vmag * rhoold(1)*Cthetav * (betam*thetav(1) + betap*d(1) - thetav_s )
  sgs_thv_flux (1)       = - vmag * Cthetav * (betam*thetav(1) + betap*d(1) - thetav_s )
! Diagnose implicit surface sensible heat flux 
  !sgs_sens_heat_flux (1) =  rhoold(1)*( sgs_thv_flux (1)/rhoold(1) &
  !                  - epsv*thetav(1)/(1.+epsv*qvspec(1)) * sgs_lat_heat_flux (1)/rho(1) )/ (1.+epsv*qvspec(1))
  sgs_sens_heat_flux (1) =  ( sgs_thv_flux (1) &
                    - epsv*thetav(1)/(1.+epsv*qvspec(1)) * sgs_lat_heat_flux (1))/ (1.+epsv*qvspec(1))
else
  !sgs_thv_flux (1)       =  sgs_thv_flux (1)  * rhoold(1)
  !sgs_sens_heat_flux (1) = sgs_sens_heat_flux (1)  * rhoold(1)
  sgs_thv_flux (1)       =  sgs_thv_flux (1)  
  sgs_sens_heat_flux (1) = sgs_sens_heat_flux (1)  
end if
! Diagnose atm. sgs fluxes 
!sgs_thv_flux (2:nzm)       = 0.5*(rhoold(1:nzm-1) + rhoold(2:nzm)) * ( & 
sgs_thv_flux (2:nzm)       = 1. * ( & 
            (-1.)/adzw(2:nzm)/dz * 0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(thetav(2:nzm)-thetav(1:nzm-1)) ) &
        + sumMthetav(2:nzm) - (betap*d(2:nzm) + betam*thetav(2:nzm)) * sumM(2:nzm))
! WL don't even bother right now to compute it, set to zero
!sgs_sens_heat_flux (2:nzm) = 0.5*(rhoold(1:nzm-1)+rhoold(2:nzm) ) * ( 2. * sgs_thv_flux (2:nzm)/(rhoold(1:nzm-1)+rhoold(2:nzm)) &
sgs_sens_heat_flux (2:nzm) = 1. * ( sgs_thv_flux (2:nzm) &
                          - epsv*0.5*(thetav(1:nzm-1)+thetav(2:nzm))/(1.+epsv*0.5*(qvspec(1:nzm-1)+qvspec(2:nzm))) &
               * sgs_lat_heat_flux (2:nzm) )  / (1.+epsv * 0.5*(qvspec(1:nzm-1)+qvspec(2:nzm)))

thvflux_ed(1)     = sgs_thv_flux (1)
thvflux_ed(nz)    = 0.
thvflux_ed(2:nzm) = (-1.)/adzw(2:nzm)/dz * 0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(thetav(2:nzm)-thetav(1:nzm-1)) )
thvflux_mf(1)     = 0.
thvflux_mf(nz)    = 0.
thvflux_mf(2:nzm) = sumMthetav(2:nzm) - (betap*d(2:nzm) + betam*thetav(2:nzm)) * sumM(2:nzm)

tend_sgs_rho_thetav = ((rhoold + dt*tend_sgs_rho)*d - rhoold * thetav)/dt

! update thetav
thetav = d

! u (input has to be air density rho and u)
if (doneuman.or.sfc_flx_fxd) then 
  ! Neuman
  call get_abcd(rhoold,u,sumMu,Cm,tk,0.,a,b,c,d,.false.,.true.,taux(1))
else
  ! Dirichlet
  call get_abcd(rhoold,u,sumMu,Cm,tk,0.,a,b,c,d,.false.,.false.)
end if
call tridiag(a,b,c,d)

if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface momentum flux
!  taux(1) = - Cm * rhoold(1) * vmag * (betap*d(1)+betam*u(1))
  taux(1) = - Cm * vmag * (betap*d(1)+betam*u(1))
else
  !taux(1) = taux(1) * rhoold(1)
  taux(1) = taux(1) 
end if
!diagnose atm. fluxes, NOTE: no contribution from mass fluxes
!taux(2:nzm) = (-0.5)*(rhoold(1:nzm-1)+rhoold(2:nzm))/adzw(2:nzm)/dz * &
taux(2:nzm) = -1./adzw(2:nzm)/dz * &
                            0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(u(2:nzm)-u(1:nzm-1)) )
!tendency
tend_sgs_rho_u = ((rhoold+dt*tend_sgs_rho)*d - rhoold*u)/dt
! update u
u=d

! v (input has to be air density rho and u)
if (doneuman.or.sfc_flx_fxd) then 
  ! Neuman
  call get_abcd(rhoold,v,sumMv,Cm,tk,0.,a,b,c,d,.false.,.true.,tauy(1))
else
  ! Dirichlet
  call get_abcd(rhoold,v,sumMv,Cm,tk,0.,a,b,c,d,.false.,.false.)
end if
call tridiag(a,b,c,d)

if (.not.(doneuman.or.sfc_flx_fxd)) then
! Diagnose implicit surface momentum flux
  !tauy(1) = - Cm * rhoold(1) * vmag * (betap*d(1)+betam*v(1))
  tauy(1) = - Cm * vmag * (betap*d(1)+betam*v(1))
else
  !tauy(1) = tauy(1) * rhoold(1)
  tauy(1) = tauy(1) 
end if
!diagnose atm. fluxes, NOTE: no contribution from mass fluxes
!tauy(2:nzm) = (-0.5)*(rhoold(1:nzm-1)+rhoold(2:nzm))/adzw(2:nzm)/dz * &
tauy(2:nzm) = -1.0/adzw(2:nzm)/dz * &
                            0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(v(2:nzm)-v(1:nzm-1)) )
!tendency
tend_sgs_rho_v = ((rhoold+dt*tend_sgs_rho)*d - rhoold*v)/dt
! update v
v=d

if (progtke) then
  ! tke (input has to be dry air density rho_d and tke per mass of dry air) 
  ! use Neuman conditions right now, since Ctke unknown
  if (dotkedirichlet) then
    call get_abcd(rho,tke,sumMtke,tk(1)/adz(1)/dz/vmag,tk,0.,a,b,c,d,.false.,.false.)
  else
    sgs_tke_flux(1) = 0.0   ! lower BC for tke transfer (questionable)
    call get_abcd(rho,tke,sumMtke,0.,tk,0.,a,b,c,d,.false.,.true.,sgs_tke_flux(1))
  end if
  call tridiag(a,b,c,d)

  if (dotkedirichlet) sgs_tke_flux(1) = - tk(1)/adz(1)/dz * (betap*d(1)+betam*tke(1))
  !diagnose atm. fluxes
  !sgs_tke_flux(2:nzm) = 0.5*(rho(1:nzm-1)+rho(2:nzm))*(            &
  sgs_tke_flux(2:nzm) =    (-1.)/adzw(2:nzm)/dz * 0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(tke(2:nzm)-tke(1:nzm-1)) ) 

  ! tendency
  tend_sgs_rho_tke = rho * (d-tke)/dt
  ! update tke
  tke = d
end if

! update tabs
tabs=thetav/(1.+epsv*qv/(qv+1.))/(p00/pres)**(rgas/cp)

! update theta
theta = thetav/(1.+epsv*qv/(qv+1.))

end subroutine sgs_tendencies

