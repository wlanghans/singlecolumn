
subroutine oceflx(pmidm1, ubot, vbot, tbot, rvbot, thbot, thvbot, zbot, ts, ssr, thvsfc)
                  

use params
use vars

implicit none	

!
! Compute ocean to atmosphere surface fluxes of sensible, latent heat
! and stress components:
!
! Assume:
!   1) Neutral 10m drag coeff: 
!         cdn = .0027/U10N + .000142 + .0000764 U10N
!   2) Neutral 10m stanton number: 
!         ctn = .0327 sqrt(cdn), unstable
!         ctn = .0180 sqrt(cdn), stable
!   3) Neutral 10m dalton number:  
!         cen = .0346 sqrt(cdn)
!
! Note:
!   1) here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   2) wind speeds should all be above a minimum speed (umin)
!
!--------------------------Code History---------------------------------
!
! Original:      Bill Large/M.Vertenstein, Sep. 1995 for CCM3.5
! Standardized:  L. Buja,     Feb 1996
! Reviewed:      B. Briegleb, March 1996
! Adopted for LES by Marat Khairoutdinov, July 1998
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
!                         
      real pmidm1  ! Bottom level pressure, mb
      real ubot    ! Bottom level u wind
      real vbot    ! Bottom level v wind
      real tbot    ! Bottom level temperature
      real rvbot   ! Bottom level mixing ratio
      real thvbot  ! Bottom level potential temperature 
      real thbot   ! Bottom level potential temperature 
      real zbot    ! Bottom level height above surface
      real ts      ! Surface temperature
      real thvsfc  ! WL Surface virtual potential temp (based on Dycore definition using qv (not rv))
      real ssr     ! WL Surface saturation mixing ratio

!                         
!---------------------------Local variables-----------------------------
!
      real tref    	  ! 2m reference temperature 
      real tstar          ! tstar
      real qstar          ! qstar, actually rstar
      real u10n           ! neutral 10 m wind speed over ocean
      real delt           ! potential T difference (K)
      real deltv          ! WL virtual potential T difference (K)
      real delrv          ! WL mixing ratio difference (kg/kg)
      real rdn            ! sqrt of neutral exchange coeff (momentum)
      real rhn            ! sqrt of neutral exchange coeff (heat)
      real ren            ! sqrt of neutral exchange coeff (tracers)          
      real rd             ! sqrt of exchange coeff (momentum)
      real rh             ! sqrt of exchange coeff (heat)
      real re             ! sqrt of exchange coeff (tracers)
      real hol            ! Ref hgt (10m) / monin-obukhov length
      real xsq            ! Temporary variable
      real xqq            ! Temporary variable
      real alz            ! ln(zbot/z10)
      real tau     ! Reference height stress
      real psimh          ! Stability funct at ref lev (momentum)
      real psixh          ! Stability funct at ref lev (heat & tracers) 
      real stable         ! Stability factor
      real rbot    	  ! Density at bottom model level
      real bn             ! exchange coef funct for interpolation
      real bh             ! exchange coef funct for interpolation
      real fac            ! interpolation factor
      real ln0            ! log factor for interpolation
      real ln3            ! log factor for interpolation
      
      real zref           ! 10m reference height
      real ztref          ! reference height for air temperature

      parameter  (zref =   10.0      , &
                  ztref=2.0 )



!--------------------------Statement functions--------------------------

      real psimhu         ! Unstable part of psimh
      real psixhu         ! Unstable part of psixh
      real cdn            ! Neutral drag coeff at bottom model level
      real xd             ! Dummy argument
      real Umps           ! Wind velocity (m/sec)

      cdn(Umps)  = 0.0027 / Umps + .000142 + .0000764 * Umps
      psimhu(xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) - 2.*atan(xd) + 1.571
      psixhu(xd) = 2. * log((1. + xd*xd)/2.)


!---------------------------------------------------------------
! Set up necessary variables
!---------------------------------------------------------------

         rbot = pmidm1*100.  / (rgas*tbot )
         delt  = thbot  - ts*(p00/pmidm1)**(rgas/cp)
         delrv = rvbot  - ssr  
         deltv = thvbot - (1.+epsv*ssr/(ssr+1.))*ts*(p00/pmidm1)**(rgas/cp)
         alz    = log(zbot /zref) 
!---------------------------------------------------------------
! First iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------

! Initial guess for roots of neutral exchange coefficients, 
! assume z/L=0. and u10n is approximated by vmag.
! Stable if (thbot > ts ).

         ! WL compute stability as theta difference between first level and surface
         stable = 0.5 + sign(0.5 , delt)
         rdn  = sqrt(cdn(vmag))
         rhn  = (1.-stable) * 0.0327 + stable * 0.018 
         ren  = 0.0346 

! Initial guess of ustar,tstar and qstar

         ustar = rdn*vmag
         tstar = rhn*delt
         qstar = ren*delrv

! Compute stability and evaluate all stability functions
! Stable if (thbot > ts or hol > 0 )

         hol = xkar*ggr*zbot*(tstar/thvbot+qstar/(1./0.606+rvbot))/ustar**2
         hol = sign( min(abs(hol),10.), hol )
         stable = 0.5 + sign(0.5 , hol)
         xsq   = max(sqrt(abs(1. - 16.*hol)) , 1.)
         xqq   = sqrt(xsq)
         psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
         psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)

! Shift 10m neutral wind speed using old rdn coefficient

         rd   = rdn / (1.+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn

! Update the neutral transfer coefficients at 10m and neutral stability

         rdn = sqrt(cdn(u10n))
         ren = 0.0346
         rhn = (1.-stable) * 0.0327 + stable * 0.018 

! Shift all coeffs to measurement height and stability

         rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
         re = ren / (1.+ren/xkar*(alz-psixh))

! Update ustar, tstar, qstar using updated, shifted coeffs 

         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delrv 

!---------------------------------------------------------------
! Second iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------

! Recompute stability & evaluate all stability functions  
! Stable if (thbot > ts or hol > 0 )
 
         hol = xkar*ggr*zbot*(tstar/thvbot+qstar/(1./0.606+rvbot))/ustar**2
         hol = sign( min(abs(hol),10.), hol )
         stable = 0.5 + sign(0.5 , hol)
         xsq   = max(sqrt(abs(1. - 16.*hol)) , 1.)
         xqq   = sqrt(xsq)
         psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
         psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)

! Shift 10m neutral wind speed using old rdn coefficient

         rd   = rdn / (1.+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn

! Update the neutral transfer coefficients at 10m and neutral stability

         rdn = sqrt(cdn(u10n))
         ren = 0.0346
         rhn = (1.-stable) * 0.0327 + stable * 0.018 

! Shift all coeffs to measurement height and stability

         rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
         re = ren / (1.+ren/xkar*(alz-psixh))

!---------------------------------------------------------------
! Compute the fluxes
!---------------------------------------------------------------

! Update ustar, tstar, qstar using updated, shifted coeffs 

         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delrv 

! Compute surface stress components and fluxes

         tau   =  rbot  * ustar * ustar 
         taux(1)  = -tau * ubot  / vmag 
         tauy(1)  = -tau * vbot  / vmag 

         sgs_sens_heat_flux(1)  = -tau * tstar / ustar / rbot
         sgs_lat_heat_flux(1)   = -tau * qstar / ustar / rbot


! WL get coefficients
         Ctheta  = rh * rd  
         Crv     = re * rd
         Cm      = ustar**2 / vmag**2

!---------------------------------------------------------------
! Following Geleyn(1988), interpolate ts to fixed height zref
!---------------------------------------------------------------

! Compute function of exchange coefficients. Assume that 
! cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
! 1/sqrt(cn )=1/rdn and sqrt(cm )/ch =1/rh 

         bn = xkar/rdn
         bh = xkar/rh

! Interpolation factor for stable and unstable cases

         ln0 = log(1.0 + (ztref/zbot )*(exp(bn) - 1.0))
         ln3 = log(1.0 + (ztref/zbot )*(exp(bn - bh) - 1.0))
         fac = (ln0 - ztref/zbot *(bn - bh))/bh * stable &
             + (ln0 - ln3)/bh * (1.-stable)
         fac = min(max(fac,0.),1.)

! Actual interpolation

         tref  = ts  + (tbot  - ts )*fac


end subroutine oceflx

 
