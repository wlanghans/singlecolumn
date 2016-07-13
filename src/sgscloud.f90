
subroutine sgscloud(k,qtflux,thetalflux,varqt,varthetal,covarqtthetal,thetal,qt,presin, & ! input
                       cfrac,ql,thetavflux)                                           ! output
 

! following Bechtold 1992, JAS and Sommerica and Deardorff (1977)
! 1) get variance of saturation deficit, sigmas
! 2) get cloud fraction and mean liquid water content from PDF scheme
! 3) get buoyancy flux w'thetav' for partially cloudy grid box

use vars
use params
use grid
implicit none

integer, intent(in) :: k
real, intent(in) :: qtflux,thetalflux,varqt,varthetal,covarqtthetal,thetal,qt, presin
real, intent(out):: cfrac, ql, thetavflux


!local variables
real :: lambdaf, alphaf, betaf, qsl, totheta, tl, c0, bflx_c, bflx_cf, sigmas


!----------------------------------------------------
! preparations
!----------------------------------------------------
totheta=(presin/p00)**(rgas/cp)
tl     = totheta * thetal

! get saturation mixing ratio at tl
qsl=qsatw(tl,presin)
! convert to sat. spec. humidity
qsl=qsl/(qsl+1.)

!----------------------------------------------------
! get variance of saturation deficit
!----------------------------------------------------
betaf = fac_cond*lcond/(tl**2*rv)
lambdaf =  (1. + betaf*qsl )**(-1.)                                                       ! Bechtold Eq (14)
alphaf = qsl * totheta * lcond/(tl**2*rv)
sigmas = (max(0.0,varqt + alphaf**2*varthetal - 2. * alphaf * covarqtthetal))**(0.5)    ! Bechtold Eq (15)
q1(k)=min(10.,max(-10.,1000.*(qt-qsl)/(1000.*sigmas+1.e-9)))


!----------------------------------------------------
! get cloud fraction and mean liquid water content
!----------------------------------------------------
cfrac = 0.5 * (1. + erf(q1(k)/1.41))                            ! Bechtold Eq (20) and Sommeria 29b
ql    = sigmas * lambdaf  * (cfrac * q1(k) + exp(-(q1(k)**2)/2.)/2.51 )  ! Bechtold Eq (21)

!----------------------------------------------------
! get buoyancy flux
!----------------------------------------------------
! unsaturated limit
bflx_cf = (1.+0.61*qt)*thetalflux + 0.61*thetal*qtflux                    ! Bechtold Eq (11) and Sommeria&Deardorff Eq(35)
! saturated limit
c0 = 1. + 0.61 * qsl - alphaf*lambdaf*thetal*(fac_cond/tl*(1.+0.61*qsl)-1.61)  ! Sommeria&Deardorff Eq(37)
bflx_c  = c0 * thetalflux + (c0*fac_cond/tl  -  1.) *thetal*qtflux       ! Sommeria&Deardorff Eq(40)

thetavflux = (1.-cfrac) * bflx_cf + cfrac * bflx_c


end subroutine sgscloud

