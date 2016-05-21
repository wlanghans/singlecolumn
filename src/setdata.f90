subroutine setdata()

use vars
use params
use grid

implicit none

integer k
real zinvbottom,zinvtop,theta0,dthetadzinv,dthetadzfree,qv0,dqvdzinv,dthetadzpbl,dqvdzpbl,dqvdzfree
CHARACTER(LEN=100), PARAMETER :: FMT1 = "(A7,A9,A7,A9,A9,A9,A8)"
CHARACTER(LEN=100), PARAMETER :: FMT2 = "(f7.0,3x,f6.1,3x,f4.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f5.2)"

if (trim(case).eq.'MOENG07') then

zinvbottom  = 1000.
zinvtop     = 1150.
theta0      = 300.
qv0         = 0.e-03
dthetadzinv = 8./(zinvtop-zinvbottom)
dqvdzinv    = -6.e-03/400.
dthetadzfree= 3.e-03
dthetadzpbl = -.0e-03
presi(1)    = p00

do k = 1,nzm
   if (z(k).lt.zinvbottom) then
      theta(k) = theta0+z(k)*dthetadzpbl
      qv(k)    = qv0
   elseif (zinvbottom.le.z(k).and.zinvtop.gt.z(k)) then
      theta(k) = theta0+zinvbottom*dthetadzpbl + (z(k)-zinvbottom) * dthetadzinv 
         qv(k) = max(0.0,qv0 + (z(k)-zinvbottom)    * dqvdzinv)
   else
      theta(k) = theta0+zinvbottom*dthetadzpbl + dthetadzinv * (zinvtop-zinvbottom) + (z(k)-zinvtop)*dthetadzfree
         qv(k) = max(0.0,(qv0 + dqvdzinv * (zinvtop-zinvbottom))) * exp(-(z(k)-zinvtop)/3000.) 
   end if

   thetav(k)= theta(k)*(1.+epsv*qv(k))
   ! get pressure at cell center and at next face from hydrostatic eqn
   pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
   presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)
 
   ! get temp
   tabs(k) = thetav(k)/(1.+epsv*qv(k))*(pres(k)/p00)**(rgas/cp)
 
   ! get rho from gas law
   rho(k) = pres(k)/rgas/tabs(k)/(1.+epsv*qv(k))*100.
   
   u(k)     = 0.
   v(k)     = 0.
   tke(k)   = 0.
   qcl(k)   = 0.
   qpl(k)   = 0.
   qci(k)   = 0.
   qpi(k)   = 0.

   
end do 


elseif (trim(case).eq.'LCLMF') then

zinvbottom  = 520.
theta0      = 298.7
qv0         = 17.e-03
dthetadzfree= (302.4 - 298.7)/960.
dthetadzpbl = -.0e-03
dqvdzinv    = -6.e-03/960.
presi(1)    = p00

do k = 1,nzm
   if (z(k).lt.zinvbottom) then
      theta(k) = theta0
      qv(k)    = qv0
   else
      theta(k) = theta0+ (z(k)-zinvbottom) * dthetadzfree
      qv(k)    = qv0 + (z(k)-zinvbottom) * dqvdzinv
   end if

   thetav(k)= theta(k)*(1.+epsv*qv(k))
   ! get pressure at cell center and at next face from hydrostatic eqn
   pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
   presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)
 
   ! get temp
   tabs(k) = thetav(k)/(1.+epsv*qv(k))*(pres(k)/p00)**(rgas/cp)
 
   ! get rho from gas law
   rho(k) = pres(k)/rgas/tabs(k)/(1.+epsv*qv(k))*100.
   
   u(k)     = 0.
   v(k)     = 0.
   tke(k)   = 0.
   qcl(k)   = 0.
   qpl(k)   = 0.
   qci(k)   = 0.
   qpi(k)   = 0.

   
end do

elseif (trim(case).eq.'TEIXEIRA04'.or.trim(case).eq.'CHEINET03'.or.trim(case).eq.'TEST') then

zinvbottom  = 1400.
theta0      = 300.
qv0         = 0.e-03
dthetadzfree= 2.1e-03
dthetadzpbl = -.0e-03
presi(1)    = p00

do k = 1,nzm
   if (z(k).lt.zinvbottom) then
      theta(k) = theta0
   else
      theta(k) = theta0+ (z(k)-zinvbottom) * dthetadzfree
   end if
   qv(k)    = qv0

   thetav(k)= theta(k)*(1.+epsv*qv(k))
   ! get pressure at cell center and at next face from hydrostatic eqn
   pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
   presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)
 
   ! get temp
   tabs(k) = thetav(k)/(1.+epsv*qv(k))*(pres(k)/p00)**(rgas/cp)
 
   ! get rho from gas law
   rho(k) = pres(k)/rgas/tabs(k)/(1.+epsv*qv(k))*100.
   
   u(k)     = 0.
   v(k)     = 0.
   tke(k)   = 0.
   qcl(k)   = 0.
   qpl(k)   = 0.
   qci(k)   = 0.
   qpi(k)   = 0.

   
end do

elseif (trim(case).eq.'WITEK11') then

zinvbottom  = 1350.
theta0      = 300.
qv0         = 5.e-03
dthetadzfree= 2.e-03
dqvdzpbl    = - 3.7e-07
dqvdzfree   = - 9.4e-07
presi(1)    = p00

do k = 1,nzm
   if (z(k).lt.zinvbottom) then
      theta(k) = theta0
      qv(k) = qv0+ z(k) * dqvdzpbl
   else
      theta(k) = theta0+ (z(k)-zinvbottom) * dthetadzfree
      qv(k) = qv0+ zinvbottom * dqvdzpbl+ (z(k)-zinvbottom) * dqvdzfree
   end if

   thetav(k)= theta(k)*(1.+epsv*qv(k))
   ! get pressure at cell center and at next face from hydrostatic eqn
   pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
   presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)
 
   ! get temp
   tabs(k) = thetav(k)/(1.+epsv*qv(k))*(pres(k)/p00)**(rgas/cp)
 
   ! get rho from gas law
   rho(k) = pres(k)/rgas/tabs(k)/(1.+epsv*qv(k))*100.
   
   u(k)     = 0.
   v(k)     = 0.
   tke(k)   = 0.
   qcl(k)   = 0.
   qpl(k)   = 0.
   qci(k)   = 0.
   qpi(k)   = 0.

   
end do

elseif (trim(case).eq.'COS') then

zinvbottom  = 1500.
theta0      = 300.
qv0         = 0.e-03
dthetadzfree= 2.1e-03
dthetadzpbl = -.0e-03
presi(1)    = p00

do k = 1,nzm
   if (z(k).le.(zinvbottom+400.).and.z(k).ge.(zinvbottom-400.) ) then
      theta(k) = theta0+ max(0.,3.* cos( 0.5*pi * (z(k)-zinvbottom)/400. )  )
   else
      theta(k) = theta0
   end if
   qv(k)    = qv0

   thetav(k)= theta(k)*(1.+epsv*qv(k))
   ! get pressure at cell center and at next face from hydrostatic eqn
   pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
   presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)

   ! get temp
   tabs(k) = thetav(k)/(1.+epsv*qv(k))*(pres(k)/p00)**(rgas/cp)

   ! get rho from gas law
   rho(k) = pres(k)/rgas/tabs(k)/(1.+epsv*qv(k))*100.

   u(k)     = 0.
   v(k)     = 0.
   tke(k)   = 0.
   qcl(k)   = 0.
   qpl(k)   = 0.
   qci(k)   = 0.
   qpi(k)   = 0.


end do
 
else ! no such case name



write(*,*) 'CASE not implemented yet'
write(*,*) 'STOPPING'


end if

! prog variables
rho   = rho
rhothv=rho*thetav
rhoqv =rho*qv
rhotke=rho*tke
rhou  =rho*u
rhov  =rho*v
pblh = zi(2)
buoy_sgs=0.0
smix=0.

tabs0=tabs

!plume properties
 UPW=0.
 UPTHL=0.
 UPTHV=0.
 UPQT=0.
 UPA=0.
 UPU=0.
 UPV=0.
 UPQC=0.
 ENT=0.

 sgs_thv_flux=0.0
 sgs_sens_heat_flux=0.0
 sgs_lat_heat_flux=0.0
 taux=0.0
 tauy=0.0
 sgs_tke_flux=0.0
 tk=0.
 tkh=0.

 ! write initial condition
 open(unit=45,file='./'//trim(case)//'/'//trim(case)//'.init',&
           form='formatted', action='write')
 write(45,FMT1)'z','p','rho','T','th','thv','qv' 
 do k=1,nzm
   write (45,FMT2) z(k),pres(k),rho(k),tabs(k),theta(k),thetav(k),qv(k)*1000.
 end do
 close(45)


end subroutine setdata

