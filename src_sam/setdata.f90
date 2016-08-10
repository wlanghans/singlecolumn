subroutine setdata()

use vars
use params
use grid

implicit none

integer k
real zinvbottom,zinvtop,theta0,dthetadzinv,dthetadzfree,qv0,dqvdzinv,dthetadzpbl,dqvdzpbl,dqvdzfree
CHARACTER(LEN=100), PARAMETER :: FMT1 = "(A7,A9,A7,A9,A9,A9,A9,A8)"
CHARACTER(LEN=100), PARAMETER :: FMT2 = "(f7.0,3x,f6.1,3x,f4.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f5.2)"

!1) set profiles  of theta, qv, qcl, qci, qpl, qpi, u, v, ug, vg, tke, etc. for each case
!2) profiles of hli, qt, and qt are derived at the end (case independent)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++
! set profiles  of theta, qv, qcl, qci, qpl, qpi, u, v, ug, vg, tke, etc. for each case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

elseif (trim(case).eq.'DYCOMS') then


zinvbottom  = 840.
theta0      = 289.0
qv0         = 9.0e-03
presi(1)    = 1017.8 
fcor        = 2. * omega * sin(31.5 * pi/180.)

w(1)  = 0.
w(nz) = -3.75e-06 * zi(nz)
do k = 1,nzm
   if (z(k).le.zinvbottom) then
      theta(k) = theta0 
      qv(k) = qv0
   else
      theta(k) = 297.5 + (z(k)-zinvbottom)**(1./3.)
      qv(k) = 1.5e-03
   end if
   thetav(k)= theta(k)*(1.+epsv*qv(k))
   ! get pressure at cell center and at next face from hydrostatic eqn
   pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
   presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)
 
   ! get temp
   tabs(k) = thetav(k)/(1.+epsv*qv(k))*(pres(k)/p00)**(rgas/cp)
 
   ! get rho from gas law
   rho(k) = pres(k)/rgas/tabs(k)/(1.+epsv*qv(k))*100.
   
   w(k) = -3.75e-06 * zi(k)
   dthetadt(k) = 0.0
   dqtdt(k) = 0.0
   ug(k) = 7.
   vg(k) = -5.5  
   u(k)     = ug(k)
   v(k)     = vg(k)
   tke(k)   = 0.
   qcl(k)   = 0.
   qpl(k)   = 0.
   qci(k)   = 0.
   qpi(k)   = 0.

   
end do

elseif (trim(case).eq.'BOMEX') then


fcor        = 0.376e-04
zinvbottom  = 520.
theta0      = 298.7
qv0         = 17.0e-03
dthetadzpbl= (302.4-theta0)/(1480.-zinvbottom)
dthetadzinv= (308.2-302.4)/(2000.-1480.)
dthetadzfree=(311.85-308.2)/1000.
dqvdzpbl    = (16.3e-03-17.0e-03)/520.
dqvdzinv    = (10.7e-03-16.3e-03)/(1480.-520.)
dqvdzfree   = (4.2e-03-10.7e-03)/(2000.-1480.)
presi(1)    = 1015.

w(1)  = 0.
w(nz) = 0.
do k = 1,nzm
   if (z(k).lt.zinvbottom) then
      theta(k) = theta0 
      qv(k) = max(0.,qv0+ z(k) * dqvdzpbl)
   elseif (z(k).ge.zinvbottom.and.z(k).lt.1480. ) then
      theta(k) = theta0+ (z(k)-zinvbottom) * dthetadzpbl
      qv(k) = max(0.,qv0+ zinvbottom * dqvdzpbl &
                 + (z(k)-zinvbottom) * dqvdzinv)
   elseif (z(k).ge.1480..and.z(k).lt.2000. ) then
      theta(k) = theta0+ (1480.-zinvbottom)* dthetadzpbl+&
                 dthetadzinv*(z(k)-1480.)
      qv(k) = max(0.,qv0+ zinvbottom * dqvdzpbl &
                 + (1480.-zinvbottom) * dqvdzinv + &
                   (z(k)-1480.)*dqvdzfree )
   else
      theta(k) = theta0+ (1480.-zinvbottom)* dthetadzpbl+&
                 dthetadzinv*(2000.-1480.)              +&
                 dthetadzfree *(z(k)-2000.)
      qv(k) = max(0.,qv0+ zinvbottom * dqvdzpbl &
                 + (1480.-zinvbottom) * dqvdzinv + &
                   (2000.-1480.)*dqvdzfree       + &
                   (3.0e-3-4.2e-03)/1000.*(z(k)-2000.))
   end if
   thetav(k)= theta(k)*(1.+epsv*qv(k))
   ! get pressure at cell center and at next face from hydrostatic eqn
   pres(k) = presi(k)*(1. - ggr/cp/thetav(k)*(p00/presi(k))**(rgas/cp) * (z(k)-zi(k)) )**(cp/rgas)
   presi(k+1) = pres(k)*(1. - ggr/cp/thetav(k)*(p00/pres(k))**(rgas/cp) * (zi(k+1)-z(k)) )**(cp/rgas)
 
   ! get temp
   tabs(k) = thetav(k)/(1.+epsv*qv(k))*(pres(k)/p00)**(rgas/cp)
 
   ! get rho from gas law
   rho(k) = pres(k)/rgas/tabs(k)/(1.+epsv*qv(k))*100.
   
   if (zi(k).lt.1500.) then
     w(k) = -.65/100./1500. * zi(k)
   elseif (zi(k).ge.1500..and.zi(k).lt.2100.) then
     w(k) = -.65/100. + 0.65/100./(2100.-1500.)*(zi(k)-1500.)
   else
     w(k) = 0.
   end if
   if (z(k).lt.1500.) then
     dthetadt(k) = -2./(24.*3600.) 
   elseif (z(k).ge.1500..and.z(k).lt.3000.)  then
     dthetadt(k) = (-2. + (z(k)-1500.) * 2.0/1500.)/(24.*3600.)
   else
     dthetadt(k) = 0.0
   end if
   if (z(k).lt.300.) then
     dqtdt(k) = -1.2d-08
   elseif (z(k).ge.300. .and. z(k).lt.500.) then
     dqtdt(k) = -1.2d-08 + (z(k)-300.) * 1.2d-08/200.
   else
     dqtdt(k) = 0.0
   end if
   ug(k) = -10. + 1.8e-03*z(k)
   vg(k) = 0.   
   if (z(k).lt.700.) then
     u(k)     = -8.75
   elseif (z(k).ge.700..and.z(k).lt.3000.) then
     u(k)     = -8.75 + (z(k)-700.) * (8.75-4.61)/2300. 
   else
     u(k)     = 0.0
   end if
   v(k)     = 0.
   u(k)     = 0.
   tke(k)   = 0.
   qcl(k)   = 0.
   qpl(k)   = 0.
   qci(k)   = 0.
   qpi(k)   = 0.

   
end do

elseif (trim(case).eq.'LANGHANS') then

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
      qv(k) = max(0.,qv0+ z(k) * dqvdzpbl)
   else
      theta(k) = theta0+ (z(k)-zinvbottom) * dthetadzfree
      qv(k) = max(0.,qv0+ zinvbottom * dqvdzpbl+ (z(k)-zinvbottom) * dqvdzfree)
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
      qv(k) = max(0.,qv0+ z(k) * dqvdzpbl)
   else
      theta(k) = theta0+ (z(k)-zinvbottom) * dthetadzfree
      qv(k) = max(0.,qv0+ zinvbottom * dqvdzpbl+ (z(k)-zinvbottom) * dqvdzfree)
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

!++++++++++++++++++++++++++++++++++++++++++++++++++++++
! profiles of hli, qt, and qp are derived at the end
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
t=cp*tabs+ggr*z - fac_cond*cp * (qcl+qpl) - fac_sub*cp * (qci+qpi)
qn=qcl + qci
qt=qv + qn
qp=qpl+qpi
thetali = theta -(p00/pres)**(rgas/cp) * (fac_cond *qcl+fac_sub*qci) 
thetar=thetav - theta*qn

pblh = zi(2)
buoy_sgs=0.0
smix=0.

tabs0=tabs

!plume properties
UPM=0. 
UPW=0.
UPT=0.
UPQT=0.
UPQCL=0.
UPQCI=0.
UPA=0.
UPU=0.
UPV=0.
UPTHV=0.
UPTABS=0.
ENT=0.
DET=0.
BUOY=0.
UPTHD=0.
UPCF=0.

 sgs_t_flux=0.0
 sgs_qt_flux=0.0
 sgs_thv_flux=0.0
 taux=0.0
 tauy=0.0
 sgs_tke_flux=0.0
 tk=0.
 tkh=0.

 ! write initial condition
 open(unit=45,file='./'//trim(case)//'/'//trim(case)//'.init',&
           form='formatted', action='write')
 write(45,FMT1)'z','p','rho','T','th','thv','thli','qv' 
 do k=1,nzm
   write (45,FMT2) z(k),pres(k),rho(k),tabs(k),theta(k),thetav(k),thetali(k),qv(k)*1000.
 end do
 close(45)


end subroutine setdata

