
subroutine sgscloudmf()

!input:
! t      ... liquid/frozen moist static energy
! presin ... pressure
! qt     ... total non-precipitating water
! qp     ... total precipitating water
! tke    ... turbulent kinetic energy

!output:
! cfrac_pdf ... cloud fraction from PDF scheme
! qn        ... total non-precip cloud condensate
! qcl       ... liquid cloud condensate
! qci       ... frozen cloud condensate
! tabs      ... absolute temperature
! theta     ... potential temperature
! thetav    ... virtual potential temperature
! thetar    ... density potential temperature


! following Bechtold 1992,1995 and Sommerica and Deardorff (1977)
! strategy:
! iterate over the following steps
! 1) get vertical gradients of hli and qt (extended here to include ice)
! 2) get variance of saturation deficit, sigmas, wherever tke.ne.0
! 3) iterate to get cloud condensate and tabs, etc. from PDF scheme

use vars
use params
use grid
use micro_params

implicit none

!local variables
integer k, kb, kc
real dtabs, an, bn, ap, bp, om, ag, omp, omn
real fac1,fac2
real fff,dfff,dqsat
real lstarn,dlstarn,lstarp,dlstarp, dtabsa, hlip, hlie, cfrac_mf2, qce, qie, qte, qtel, qne, tabse, qltot,qitot
integer niter
real :: lambdaf, alphaf, betaf, qsl, totheta, tl, sigmas
real,dimension (2) :: tabs1, tabs2

an = 1./(tbgmax-tbgmin) 
bn = tbgmin * an
ap = 1./(tprmax-tprmin) 
bp = tprmin * ap
fac1 = fac_cond+(1+bp)*fac_fus
fac2 = fac_fus*ap
ag = 1./(tgrmax-tgrmin) 


do k=1,nzm

totheta=(pres(k)/p00)**(rgas/cp)

if (k.eq.1.or.tke(k).eq.0.0) then

! no PDF scheme needed if no variability (or if in first layer)
! if condensate is present, then it comes from all-or-nothing adjustment
! and will not modified here
 tabs(k)=(t(k)-ggr*z(k))/cp + fac_cond*(qcl(k)+qpl(k)) + fac_fus*(qci(k)+qpi(k)) 
 thetali(k) = (tabs(k)- fac_cond*(qcl(k)) - fac_fus*(qci(k)) ) / totheta 
 thetav(k)  = tabs(k)/totheta * (1.+epsv*qv(k)) 
 ! density potential temperature will later be used in edmf to get buoyancy. 
 ! since qp is thought to be spread out evenly over grid-box it will not contribute to B
 thetar(k)  = thetav(k) - tabs(k)/totheta * qn(k)
 cfrac_pdf(k) = 0.
 if (qn(k).gt.0.0) then
   cfrac_tot(k) = 1.
 else 
   cfrac_tot(k) = 0.0
 end if

else

! get plume moist static energy
 hlip = 0.5*(cp*tabs_mf(k) + ggr*zi(k) - lcond * qcsgs_mf(k) - lsub * qisgs_mf(k) +&
            cp*tabs_mf(k+1) + ggr*zi(k+1) - lcond * qcsgs_mf(k+1) - lsub * qisgs_mf(k+1))
 cfrac_mf2 = 0.5*(cfrac_mf(k)+cfrac_mf(k+1))

! moist static energy in environment 
 hlie   = ((t(k)+lcond*qpl(k)+lsub*qpi(k)-cfrac_mf2*hlip))/(1.-cfrac_mf2)
! total water in environment
 qte = (qt(k) - cfrac_mf2*0.5*(qtsgs_mf(k)+qtsgs_mf(k+1)))/(1.-cfrac_mf2)
! first guess for tabs in environment outside of plumes: no condensate in environment
 tabse  = (hlie - ggr*z(k))/cp
! thetali in environment
 thetali(k) = tabs(k) / totheta 
 
 niter=0
 dtabs = 100.
 do while(abs(dtabs).gt.0.01.and.niter.lt.20)

 if (niter.gt.0) then 
   tabse= tabse + 0.5 * dtabs
   thetali(k) = tabs(k) / totheta - (lcond * qce - lsub*qie)/totheta
 end if

    ! environmental properties

    ! compute gradients of qt and thetali
    thetaligrad(k) = (thetali(k)-thetali(k-1))/(z(k)-z(k-1))
    qtgrad(k)      = (qte-qtel)/(z(k)-z(k-1))

    ! get liquid water temperature Tl
    tl     = totheta * thetali(k)
    ! get saturation mixing ratio at Tl
    qsl=qsatw(tl,pres(k))
    betaf = fac_cond*lcond/(tl**2*rv)
    lambdaf =  (1. + betaf*qsl )**(-1.)                                                       ! Bechtold Eq (14)
    alphaf = qsl * totheta * lcond/(tl**2*rv)
    ! use constant length scale of 250 m
    sigmas = 250.*(max(0.0,qtgrad(k)**2. + alphaf**2*thetaligrad(k)**2. - 2. * alphaf * thetaligrad(k)*qtgrad(k)))**(0.5)

    ! compute normalized saturation deficit
    q1(k)= min(10.,max(-10.,1000.*(qte-qsl)/(1000.*sigmas+1.e-7)))

    ! get cloud fraction and mean liquid water content
    cfrac_pdf(k) = min(1.0,max(0.,0.5 * (1. + erf(q1(k)/1.41))))                           ! Bechtold Eq (20)
    qne = max(0.0,sigmas * lambdaf  * (cfrac_pdf(k) * q1(k) + exp(-q1(k)**2/2.)/2.51 ))  ! Bechtold Eq (21)
 
    ! partition condensate into liquid and ice using temperature
    omn = max(0.,min(1.,(tabs(k)-tbgmin)*an))
    qce = qne*omn
    qie = qne*(1.-omn)

    ! get new estimates
    niter = niter+1
    dtabs = (hlie - ggr*z(k) + lcond * qce + lsub*qie )/cp - tabse

    end do ! end of iteration

    qltot = (1.-cfrac_mf2) * qce + 0.5 * cfrac_mf2 * (qcsgs_mf(k)+qcsgs_mf(k+1))
    qitot = (1.-cfrac_mf2) * qie + 0.5 * cfrac_mf2 * (qisgs_mf(k)+qisgs_mf(k+1))
    cfrac_tot(k) = cfrac_mf2 + (1.-cfrac_mf2) * cfrac_pdf(k)
    qn(k) = qltot + qitot
    qv(k) = qt(k) - qn(k)
    
    tabs(k) = (t(k)-ggr*z(k))/cp + fac_cond*(qpl(k)+qltot) + fac_sub*(qpi(k)+qitot) 
    thetav(k) = tabs(k)/totheta  * (1.+epsv*qv(k))
    thetar(k) = tabs(k)/totheta  * (1.+epsv*qv(k)-qn(k))

end if ! if tke is present

! store qte for computation of gradient at next level higher up
qtel=qte 

end do ! loop over k


!do k=1,nzm

! get buoyancy flux 
  ! get liquid water temperature Tl
!  tl     = totheta * thetali(k)
  ! get saturation mixing ratio at Tl
!  qsl=qsatw(tl,pres(k))
  ! unsaturated limit
!  bflx_cf = (1.+epsv*qt)*thetalflux + 0.61*thetali(k)*qtflux       ! Bechtold Eq (11) and Sommeria&Deardorff Eq(35)
  ! saturated limit
!  c0 = 1. + epsv * qsl - alphaf*lambdaf*thetal*(fac_cond/tl*(1.+0.61*qsl)-1.61)  ! Sommeria&Deardorff Eq(37)
!  bflx_c  = c0 * thetalflux + (c0*fac_cond/tl  -  1.) *thetal*qtflux       ! Sommeria&Deardorff Eq(40)

  ! linear combination
!  thetavflux = (1.-cfrac) * bflx_cf + cfrac * bflx_c

!end do


end subroutine sgscloudmf

