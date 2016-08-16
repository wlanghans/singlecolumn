
subroutine sgscloud()

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
real lstarn,dlstarn,lstarp,dlstarp
integer niter
real :: lambdaf, alphaf, betaf, qsl, totheta, tl
real :: tabs1, tabs2

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
! since grid-scale saturation will not occur, no condensate can be present
 qn(k)       = 0.
 qcl(k)      = 0.
 qci(k)      = 0.
 qv(k)       = qt(k)
 tabs1=(t(k)-ggr*z(k))/cp
! splitting qp into solid and liquid part
 tabs(k)=(tabs1+fac1*qp(k))/(1.+fac2*qp(k))
 thetali(k) = tabs2 / totheta 
 thetav(k)  = thetali(k) * (1.+epsv*qv(k)) 
 thetar(k)  = thetav(k)
 cfrac_pdf(k) = 0.

else

! first guess: no condensate
 tabs1=(t(k)-ggr*z(k))/cp
! second guess: splitting qp into solid and liquid part
 tabs2=(tabs1+fac1*qp(k))/(1.+fac2*qp(k))

! Warm cloud:
    if(tabs2.ge.tbgmax) then
      tabs(k)=tabs1+fac_cond*qp(k)
! Ice cloud:
    elseif(tabs2.le.tbgmin) then
      tabs(k)=tabs1+fac_sub*qp(k)
! Mixed-phase cloud:
    else
      om = an*tabs2-bn
      tabs(k)=tabs2
    endif

    ! first guess for
    ! frozen/liquid potential temperature thetali = theta*(1 - theta/t (Lv/cp*qcl + Ls/cp*qci )
    thetali(k) = tabs(k) / totheta 

    niter=0
    dtabs = 100.
    do while(abs(dtabs).gt.0.01.and.niter.lt.20)

    if (niter.gt.0) then 
      thetali(k) = thetali(k) + 0.5 * dtabs 
      tabs(k) = (t(k)-ggr*z(k))/cp + fac_cond*(qcl(k)+qpl(k)) + fac_sub*(qci(k)+qpi(k)) 
    end if

    ! compute gradients of qt and thetali
    thetaligrad(k) = (thetali(k)-thetali(k-1))/(z(k)-z(k-1))
    qtgrad(k)      = (qt(k)-qt(k-1))/(z(k)-z(k-1))

    ! get liquid water temperature Tl
    tl     = totheta * thetali(k)
    ! get saturation mixing ratio at Tl
    qsl=qsatw(tl,pres(k))
    betaf = fac_cond*lcond/(tl**2*rv)
    lambdaf =  (1. + betaf*qsl )**(-1.)                                                       ! Bechtold Eq (14)
    alphaf = qsl * totheta * lcond/(tl**2*rv)
    ! use constant length scale of 250 m
    sigmas(k) = 250.*(max(0.0,qtgrad(k)**2. + alphaf**2*thetaligrad(k)**2. - 2. * alphaf * thetaligrad(k)*qtgrad(k)))**(0.5)

    ! compute normalized saturation deficit
    q1(k)= min(10.,max(-10.,1000.*(qt(k)-qsl)/(1000.*sigmas(k)+1.e-7)))

    ! get cloud fraction and mean liquid water content
    cfrac_pdf(k) = min(1.0,max(0.,0.5 * (1. + erf(q1(k)/1.41))))                           ! Bechtold Eq (20)
    qn(k) = max(0.0,sigmas(k) * lambdaf  * (cfrac_pdf(k) * q1(k) + exp(-q1(k)**2/2.)/2.51 ))  ! Bechtold Eq (21)
    qv(k) = qt(k) - qn(k)
 
    ! partition condensate into liquid and ice using temperature
    omn = max(0.,min(1.,(tabs(k)-tbgmin)*an))
    qcl(k) = qn(k)*omn
    qci(k) = qn(k)*(1.-omn)

    ! partition precip into liquid and ice using temperature
    omp = max(0.,min(1.,(tabs(k)-tprmin)*ap))
    qpl(k) = qp(k)*omp
    qpi(k) = qp(k)*(1.-omp)
    
    ! get new estimates
    niter = niter+1
    dtabs      = tabs(k) / totheta  - (fac_cond*qcl(k) + fac_sub*qci(k))/totheta &
                 - thetali(k)

    end do ! end of iteration

    tabs(k) = (t(k)-ggr*z(k))/cp + fac_cond*(qcl(k)+qpl(k)) + fac_sub*(qci(k)+qpi(k)) 
    thetav(k) = tabs(k)/totheta  * (1.+epsv*qv(k))
    thetar(k) = tabs(k)/totheta  * (1.+epsv*qv(k)-qn(k)-qp(k))

end if ! if tke is present

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


end subroutine sgscloud

