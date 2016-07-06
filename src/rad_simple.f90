
subroutine rad_simple
 	
!	Simple Interactive Radiation
!  Coded by Ping Zhu for DycomsII LES intercomparison


use grid
use vars
use params
implicit none
	
real f0, xk, coef, coef1
integer i,k

!pzhu
real qq1,qq2,hfact,flux(nz),FTHRL(nzm),qt !bloss: flux(nz) rather than flux(nzm)
integer j,kk,itop,item3
!pzhu-end

!bloss
real deltaq(nzm), cpmassl(nzm), qzinf, qzeroz
real, parameter :: cp_spec = 1015. ! from DycomsII LES intercomparison specs


if(.not.dolongwave) return

coef=70.
coef1=22.
f0=3.75e-6
xk=85.

do k=1,nz
radlwdn(k) =0.
enddo

do k = 1,nzm
   !WL: remove density here since prog var is density times thetav
   cpmassl(k) = cp_spec*dz*adz(k) ! thermal mass of model layer.
end do

   ! search for inversion height (highest level where qt=8 g/kg) 
   !  and compute optical depth increments.
   itop = 1
   deltaq = 0.
   qzinf = 0. ! holds accumulated optical depth between z and domain top
   qzeroz = 0.! holds accumulated optical depth between surface and z
   do k = 1,nzm
      ! optical depth only includes that due to liquid water
      ! qlsgs_ed and qlsgs_mf need to be mixing ratios here too
      if(qcl(k)+qci(k)+qlsgs_ed(k)+qlsgs_mf(k).gt.0.) deltaq(k) = xk*rho(k)*(qcl(k)+qci(k)+qlsgs_ed(k)+qlsgs_mf(k))*adz(k)*dz

      ! inversion height is top of highest layer w/q>8g/kg.
      ! convert rt to qt
      qt=(qv(k)+qcl(k)+qci(k)+qlsgs_ed(k)+qlsgs_mf(k))/(1.+qv(k))
      if(qt.gt.0.008) itop=max(itop,k+1) ! note zi(k+1) is inversion hgt

      ! accumulate optical depth in qzinf
      qzinf = qzinf + deltaq(k) 
   end do
   ! note: qzinf now holds total optical depth of column (due to cloud).
   !       qzeroz is initialized to zero since first level is at surface.

   ! compute net upward longwave flux up to inversion height.
   do k = 1,itop
      flux(k) = coef*exp(-qzinf) + coef1*exp(-qzeroz)
      qzinf  = qzinf  - deltaq(k)
      qzeroz = qzeroz + deltaq(k)
   end do

   ! compute net upward longwave flux above inversion height.
   ! this includes correction for clearsky fluxes which balances
   ! the prescribed subsidence heating above the inversion.
   do k = itop+1,nzm
      flux(k) = coef*exp(-qzinf) + coef1*exp(-qzeroz) &
           !+ cp_spec*0.5*(rho(k)/(1.-qv(k))+rho(k-1)/(1.-qv(k-1)))*f0*(0.25*zi(k)+0.75*zi(itop)) &
           !                    *(zi(k)-zi(itop))**(1./3.)
           ! WL: changed to Stevens (2005), MWR, formulation
           + cp_spec*0.5*(rho(k)/(1.-qv(k))+rho(k-1)/(1.-qv(k-1)))*f0*(0.25*(zi(k)-zi(itop))**(4./3.) &
                               + zi(itop)*(zi(k)-zi(itop))**(1./3.))
      qzinf  = qzinf  - deltaq(k)
      qzeroz = qzeroz + deltaq(k)
   end do
   flux(nz) = coef*exp(-qzinf) + coef1*exp(-qzeroz) &
   !WL: linear interpolate density to model's top and adapt Stevens (2005), MWR, formulation
        + cp_spec*(rho(nzm)+ 0.5*adz(nzm)*dz* (rho(nzm)/(1.-qv(nzm))-rho(nzm-1)/(1.-qv(nzm-1)))/adzw(nzm)/dz)&
               *f0*(0.25*(zi(nz)-zi(itop))**(4./3.) + zi(itop)*(zi(nz)-zi(itop))**(1./3.))

   ! note that our flux differs from the specification in that it
   ! uses the local density (rather than the interface density) in
   ! computing the correction to the flux above the inversion. 
   ! Formulating things in this way ensures a rough balance between
   ! the subsidence heating and radiative cooling above the inversion.

   ! compute radiative heating as divergence of net upward lw flux.
   do k=1,nzm
      FTHRL(k)=-(flux(k+1)-flux(k))/cpmassl(k)
      radlwdn(k) = flux(k)
      tend_rad_rho_thetav(k) = (1.+epsv*qv(k)/(1.+qv(k))) * & 
      (p00/pres(k))**(rgas/cp)*FTHRL(k) 
   enddo
   

end 




