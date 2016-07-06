subroutine surface()
	
use vars
use params
implicit none
	
real diag_ustar
integer i,j
real(8) buffer(2), buffer1(2)

real, parameter :: umin=1.0


 ! WL compute boundary values and wind speed which are needed independently of SFC_FLX_FXD
 vmag   = max(umin, sqrt(u(1) **2+v(1) **2))
 if(sst+t00.gt.271.) then
    r_s = salt_factor*qsatw(sst+t00,pres(1))
 else
    r_s = qsati(sst+t00,pres(1))
 end if
 thetav_s = (sst+t00)*(1.+epsv*r_s/(r_s+1.0)) * (p00/pres(1))**(rgas/cp)
 theta_s  = (sst+t00)* (p00/pres(1))**(rgas/cp)



if(.not.SFC_FLX_FXD) then

  if(OCEAN) then

         ! WL here requires mixing ratio at first level NOT specific humidity
         ! WL needs temperature, pot. temp, and virt. pot. temp at first model level
         ! NOTE: either a flux boundary condition or a Dirichlet can be used to solve for tendencies
         call oceflx( pres(1),u(1),v(1),tabs(1), &
                      qv(1)/(1.-qv(1)),thetav(1)/(1.+epsv* qv(1)),thetav(1),z(1),sst+t00, r_s, thetav_s)
         sgs_thv_flux(1) = sgs_sens_heat_flux(1) * (1.+epsv* qv(1)) + epsv * thetav(1)/(1.+epsv* qv(1)) * sgs_lat_heat_flux(1)
         Cthetav = -sgs_thv_flux(1) / (thetav(1) - thetav_s)/vmag

         if (SFC_CS_FXD) then
           ! reset transfer coefficients 
           Ctheta  = Cs_in
           Cthetav = Cs_in
           Crv     = Cs_in
           ! recompute fluxes in case explicit Neuman conditions are used
           sgs_sens_heat_flux(1)     = -Ctheta * vmag * (thetav(1)/(1.+epsv* qv(1)) - theta_s)
           sgs_lat_heat_flux(1) = -Crv  * vmag * (qv(1)/(1.-qv(1)) - r_s)
           sgs_thv_flux(1) = sgs_sens_heat_flux(1) * (1.+epsv* qv(1)) + epsv * thetav(1)/(1.+epsv* qv(1)) * sgs_lat_heat_flux(1)
         end if
                      
         ! WL re-compute Cm if tau is fixed
         if(SFC_TAU_FXD) then
         ! NOTE: if flux is prescribed, then we will use neuman conditions. Cm will actually not be used
            Cm = tau0 / vmag**2
            taux(1) = -tau0 * u(1)  / vmag 
            tauy(1) = -tau0 * v(1)  / vmag 
         elseif (SFC_CM_FXD) then
            Cm = Cm_in
            tau0 =  Cm * vmag**2
            taux(1) = -tau0 * u(1)  / vmag
            tauy(1) = -tau0 * v(1)  / vmag
         end if

  end if ! OCEAN

 ! WL land case not available yet
  if(LAND) then
     print *, 'LAND case is not available yet unless SFC_FLX_FXD = .True.'
     print *, 'STOPPING'
     stop
  end if ! LAND


else   ! IF SFC_FLX_FXD=True

  ! WL invert bulk transfer formulae to get drag coefficients C=-flux/dtheta/vmag, which are needed for implicit formulation
  ! WL coefficients might be negative just to ensure that fluxes equal their prescribed values; in this case the surface values are meaningless (that's ok)
  ! NOTE: if fluxes are prescribed, we will use neuman conditions. Cs will not be used
  sgs_sens_heat_flux(1)  = fluxt0 
  Ctheta = - fluxt0 /( thetav(1)/(1.+epsv* qv(1))  - theta_s) /vmag
  sgs_lat_heat_flux(1)  = fluxq0 
  Crv    = - fluxq0 / ( qv(1)/(1.-qv(1)) - r_s  ) / vmag 
  sgs_thv_flux(1) = sgs_sens_heat_flux(1) * (1.+epsv* qv(1)) + epsv * thetav(1)/(1.+epsv* qv(1)) * sgs_lat_heat_flux(1)
  Cthetav = -sgs_thv_flux(1) / (thetav(1) - thetav_s)/vmag


    if(.not.SFC_TAU_FXD) then
      if(OCEAN) z0 = 0.0001  ! for LAND z0 should be set in namelist (default z0=0.035)

      ustar = diag_ustar(z(1),  &
                ggr/tabs(1)* (fluxt0+epsv*tabs(1)*fluxq0),vmag,z0)
      tau0  = ustar**2
      Cm = tau0 / vmag**2
     
      if (SFC_CM_FXD) then
         Cm = Cm_in
         tau0 = Cm * vmag**2
      end if

    else
      Cm = tau0 / vmag**2
      ustar=sqrt(tau0)
    end if ! .not.SFC_TAU_FXD

    taux(1) = -tau0 * u(1)  / vmag 
    tauy(1) = -tau0 * v(1)  / vmag 

end if ! SFC_FLX_FXD


end




! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below 
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
! Code corrected 8th June 1999 (obukhov length was wrong way up,
! so now used as reciprocal of obukhov length)

      real function diag_ustar(z,bflx,wnd,z0)

      implicit none
      real, parameter      :: vonk =  0.4   ! von Karmans constant
      real, parameter      :: g    = 9.81   ! gravitational acceleration
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "
      real, parameter      :: eps  = 1.e-10 ! non-zero, small number

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: z0            ! momentum roughness height

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, rlmo, ustar

      lnz   = log(z/z0) 
      klnz  = vonk/lnz              
      c1    = 3.14159/2. - 3.*log(2.)

      ustar =  wnd*klnz
      if (bflx /= 0.0) then 
        do iterate=1,4
          rlmo   = -bflx * vonk/(ustar**3 + eps)   !reciprocal of
                                                   !obukhov length
          zeta  = z*rlmo
          if (zeta > 0.) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
          else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
            ustar = wnd*vonk/(lnz - psi1)
          end if
        end do
      end if

      diag_ustar = ustar

      return
      end function diag_ustar
! ----------------------------------------------------------------------

