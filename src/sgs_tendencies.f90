subroutine sgs_tendencies()

use vars
use params
use grid

implicit none

! local variables

real,dimension(nzm) :: rhoold, a, b, c, d, qvspec

!qv holds mixing ratio
!rho holds rhod

rhoold = rho * (1.+qv)
qvspec  = qv/(1.+qv)

! water vapor mixing ratio rv (input has to be dry air density rho_d and vapor mixing ratio rv)
!if (.not.sfc_flx_fxd) then ! Dirichlet
!call get_abcd(rho,qv,sumMrv,Crv,tkh,r_s,a,b,c,d,sfc_flx_fxd)
!else ! Neuman
!sgs_lat_heat_flux (1) = fluxq0 ! explicit sfc flux
call get_abcd(rho,qv,sumMrv,Crv,tkh,r_s,a,b,c,d,.true.,sgs_lat_heat_flux (1))
!end if
call tridiag(a,b,c,d)
! Diagnose implicit surface flux in kinematic units
!sgs_lat_heat_flux (1) = - vmag *rho(1)* Crv * (betam*qv(1) + betap*d(1) - r_s )
!if (.not.sfc_flx_fxd) then
!  sgs_lat_heat_flux (1) = - vmag * Crv * (betam*qv(1) + betap*d(1) - r_s )
!end if
! Diagnose atm. sgs fluxes 
!sgs_lat_heat_flux (2:nzm) = -0.5*(rho(1:nzm-1) + rho(2:nzm)) /adzw(2:nzm)/dz * &
!                            0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
!                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qv(2:nzm)-qv(1:nzm-1)) )
sgs_lat_heat_flux (2:nzm) = - 1.0 /adzw(2:nzm)/dz * &
                            0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(qv(2:nzm)-qv(1:nzm-1)) ) &
                          + sumMrv(2:nzm) - 0.5 * (betap*(d(2:nzm) + d(1:nzm-1)) + betam*(qv(2:nzm)+qv(1:nzm-1))) * sumM(2:nzm)

! compute rho and rhov tendency
tend_sgs_rho_qv = rho * (d-qv)/dt    ! =rho_d * (rvnew - rvold) /dt
tend_sgs_rho    = tend_sgs_rho_qv

! update qv
! d is mixing ratio
qv = d 

! WL ADD DIFFUSION OF ADDITIONAL WATER VARIABLE, to zero for now 
tend_sgs_rho_qt = 0.

! thetav (input has to be air density rho and thetav)
!if (.not.sfc_flx_fxd) then ! Dirichlet
!call get_abcd(rhoold,thetav,sumMthetav,Cthetav,tkh,thetav_s,a,b,c,d,sfc_flx_fxd)
!else ! Neuman
!sgs_thv_flux (1)       = fluxt0*(1.+epsv * qvspec(1))+epsv*thetav(1)/(1.+epsv*qvspec(1)) * fluxq0 ! explicit surface flux
call get_abcd(rhoold,thetav,sumMthetav,Cthetav,tkh,thetav_s,a,b,c,d,.true.,sgs_thv_flux (1))
!end if
call tridiag(a,b,c,d)
tend_sgs_rho_thetav = ((rhoold + dt*tend_sgs_rho)*d - rhoold * thetav)/dt
! Diagnose implicit surface pot. temp. flux 
!sgs_thv_flux (1)       = - vmag * rhoold(1)*Cthetav * (betam*thetav(1) + betap*d(1) - thetav_s )
!if (.not.sfc_flx_fxd) then
!sgs_thv_flux (1)       = - vmag * Cthetav * (betam*thetav(1) + betap*d(1) - thetav_s )
!end if
! Diagnose implicit surface sensible heat flux 
!sgs_sens_heat_flux (1) =  rhoold(1)*( sgs_thv_flux (1)/rhoold(1) &
!                          - epsv*thetav(1)/(1.+epsv*qvspec(1)) * sgs_lat_heat_flux (1)/rho(1) ) &
!                                              / (1.+epsv * qvspec(1))
!sgs_sens_heat_flux (1) =  ( sgs_thv_flux (1) &
!                          - epsv*thetav(1)/(1.+epsv*qvspec(1)) * sgs_lat_heat_flux (1) ) &
!                                              / (1.+epsv * qvspec(1))
! Diagnose atm. sgs fluxes 
!sgs_thv_flux (2:nzm)       = -0.5*(rhoold(1:nzm-1) + rhoold(2:nzm)) /adzw(2:nzm)/dz * &
!                            0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
!                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(thetav(2:nzm)-thetav(1:nzm-1)) )
sgs_thv_flux (2:nzm)       = -1.0/adzw(2:nzm)/dz * &
                            0.5*(tkh(1:nzm-1) + tkh(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(thetav(2:nzm)-thetav(1:nzm-1)) ) &
                          + sumMthetav(2:nzm) - 0.5 * (betap*(d(2:nzm) + d(1:nzm-1)) + betam*(thetav(2:nzm)+thetav(1:nzm-1))) * sumM(2:nzm)
! WL don't even bother right now to compute it, set to zero
sgs_sens_heat_flux (2:nzm) = ( sgs_thv_flux (2:nzm) &
                          - epsv*0.5*(thetav(1:nzm-1)+thetav(2:nzm))/(1.+epsv*0.5*(qvspec(1:nzm-1)+qvspec(2:nzm))) &
                           * sgs_lat_heat_flux (2:nzm) )  / (1.+epsv * 0.5*(qvspec(1:nzm-1)+qvspec(2:nzm)))

! update thetav
thetav = d

! u (input has to be air density rho and u)
!if (.not.sfc_tau_fxd) then
!call get_abcd(rhoold,u,sumMu,Cm,tk,0.,a,b,c,d,sfc_tau_fxd)
!else
!taux(1)=-u(1)/vmag*tau0
call get_abcd(rhoold,u,sumMu,Cm,tk,0.,a,b,c,d,.true.,taux(1))
!end if
call tridiag(a,b,c,d)
tend_sgs_rho_u = ((rhoold+dt*tend_sgs_rho)*d - rhoold*u)/dt


!diagnose atm. fluxes
taux(2:nzm) = -1.0/adzw(2:nzm)/dz * &
                            0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(u(2:nzm)-u(1:nzm-1)) )
! update u
u=d

! v (input has to be air density rho and v)
!if (.not.sfc_tau_fxd) then
!call get_abcd(rhoold,v,sumMv,Cm,tk,0.,a,b,c,d,sfc_tau_fxd)
!else
!tauy(1)=-v(1)/vmag*tau0
call get_abcd(rhoold,v,sumMv,Cm,tk,0.,a,b,c,d,.true.,tauy(1))
!end if
call tridiag(a,b,c,d)
tend_sgs_rho_v = ((rhoold+dt*tend_sgs_rho)*d - rhoold*v)/dt


!diagnose atm. fluxes
tauy(2:nzm) = -1.0/adzw(2:nzm)/dz * &
                            0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(v(2:nzm)-v(1:nzm-1)) )
! update v
v=d

if (progtke) then
  ! tke (input has to be dry air density rho_d and tke per mass of dry air) 
  sgs_tke_flux(1) = 0.0
  call get_abcd(rho,tke,sumMtke,0.,tk,0.,a,b,c,d,.true.,sgs_tke_flux(1))
  call tridiag(a,b,c,d)
  tend_sgs_rho_tke = rho * (d-tke)/dt


  !diagnose atm. fluxes
  sgs_tke_flux(2:nzm) = -1.0/adzw(2:nzm)/dz * &
                            0.5*(tk(1:nzm-1) + tk(2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(tke(2:nzm)-tke(1:nzm-1)) ) &
                          + sumMtke(2:nzm) - 0.5 * (betap*(d(2:nzm) + d(1:nzm-1)) + betam*(tke(2:nzm)+tke(1:nzm-1))) * sumM(2:nzm)

  ! update tke
  tke = d
end if

! update tabs
tabs=thetav/(1.+epsv*qv/(qv+1.))/(p00/pres)**(rgas/cp)

! update theta
theta = thetav/(1.+epsv*qv/(qv+1.))

end subroutine sgs_tendencies

