! WL this routine provides the fluxes sumM(nz), sumMrv(nz), sumMu(nz), sumMv(nz), sumMthetav(nz), sumMtke(nz)
! WL that result from the multiplume model provided by Kay Suselj and adapted by WL for the Dycore
subroutine edmf()

use vars
use grid
use params
implicit none

! ==============================================================================
! Multiplume stochastic EDMF model (some comments by Kay)
!  - does vertical  integration for the single horizontal grid-point (loop over horizontal points must be outside of this routine) 
!  - nup plumes initialized at the surface, integrated until the level where their vertical velocity becomes zero
!  - surface conditions - following Cheinet 2003
!  - plume integration - Suselj et al., 2013 and Suselj et al., 2014
!  - entrainment rate - stochastic process, following Poisson distribution
!
! output needed for tendencies (on faces)
! sumM = sum (a_i x w_i) [ms-1]
! sumMthetav = sum (a_i x w_i x thv_i) [Kms-1]
! sumMrv = sum (a_i x w_i x rv_i)  [ms-1]
! sumMu = sum (a_i x w_i x u_i)     = 0. [ms-1]^2
! sumMv = sum (a_i x w_i x v_i)     = 0. [ms-1]^2
! sumMtke = sum (a_i x w_i x tke_i) = 0. [ms-1]^3
!
! Input needed
! zw ... heights of the half levels
! p,u,v,thl,thv,qt (kts:kte) - profiles of grid-box mean variables 
! ust [ms-1],wthl [Kms-1],wqt [ms-1] - surface fluxes (ustar, sensible and latent heat)_
! WL with our implicit scheme wthv is unknown so far but, -- in line with the explicit computation of all mass fluxes) 
! WL -- let's just use the explicit fluxes for this purpose here
! pblh - boundary layer height
! ==============================================================================
! local variables:
 ! entrainment variables     
      REAl,DIMENSION(1:nzm,1:nup)    :: ENTf
      INTEGER,DIMENSION(1:nzm,1:nup) :: ENTi

      INTEGER :: K,I
      REAL :: wthv,wqt,wthl,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,sigmaTHV,zs, &
           pwmin,pwmax,wmin,wmax,wlv,wtv
      REAL :: B,B2,B3,QTn,THLn,THVn,QCn,Un,Vn,Wn2,EntEXP,EntW, hlp

! w parameters
      REAL,PARAMETER :: &
        &Wa=2./3., &
        &Wb=0.002,&
        &Wc=1.5

! entrainment parameters
      REAL,PARAMETER :: &
        & L0=100.,&
        & ENT0=0.1

!initialize plume properties
 UPW=0.
 UPTHL=0.
 UPT=0.
 UPTHV=0.
 UPQT=0.
 UPA=0.
 UPU=0.
 UPV=0.
 UPQC=0.
 ENT=0.

 
 ! surface fluxes
 ! sensible heat flux (K m s-1)
 wthl = sgs_sens_heat_flux(1)
 ! latent heat flux (m s-1)
 wqt  = sgs_lat_heat_flux(1) 
 ! virtual pot. temp flux
 wthv = sgs_thv_flux(1) 

 ! quit in case of a non-positive buoyancy flux
 if (wthv.le.0.0) return
 
 ! get entrainment rate
  do i=1,nup
    do k=1,nzm
      ENTf(k,i)=(zi(k+1)-zi(k))/L0
     enddo
   enddo

 ! get random Poisson number
   !call Poisson(1,nup,1,nzm,ENTf,ENTi)

 ! Calculate entrainment rate: Ent=Ent0/dz*P(dz/L0) for each layer             
    do i=1,nup
    do k=1,nzm
      ENT(k,i)= 1.e-03 ! real(ENTi(k,i))*Ent0/(zi(k+1)-zi(k))
    enddo
    enddo

! set initial conditions for updrafts
    zs=50.
    pwmin=0.5
    pwmax=3.

! see Lenschow et al. (1980), JAS
    wstar=max(0.,(ggr/thetav(1)*wthv*pblh)**(1./3.))
    qstar=wqt/wstar
    thstar=wthl/wstar 
    sigmaW=1.34*wstar*(zs/pblh)**(1./3.)*(1.-0.8*zs/pblh)
    sigmaQT=1.34*qstar*(zs/pblh)**(-1./3.)
    sigmaTH=1.34*thstar*(zs/pblh)**(-1./3.)
 
! now get sigmaTHV from linearization of pot. temp. and using a corr. coeff between qv and theta of 0.75 (see Sobjan 1991)
    sigmaTHV=sqrt(sigmaTH**2+epsv**2*(thetav(1)/(1.+epsv*qv(1)/(1.+qv(1))))**2*sigmaQT**2 &
        + 2.*epsv*(thetav(1)/(1.+epsv*qv(1)/(1.+qv(1))))*0.75*sigmaTH*sigmaQT)

    wmin=sigmaW*pwmin
    wmax=sigmaW*pwmax

    DO i=1,nup

       wlv=wmin+(wmax-wmin)/nup*(i-1)
       wtv=wmin+(wmax-wmin)/nup*i

       UPA(1,I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))
       !UPW(1,I)=0.5*(wlv+wtv)
       UPW(1,I)=  sigmaW/(UPA(1,I)*sqrt(2.*pi)) * (exp(-(wlv**2)/(2.*sigmaW**2)) - exp(-(wtv**2)/(2.*sigmaW**2))  ) 

       UPU(1,I)=u(1)
       UPV(1,I)=v(1)

       UPQC(1,I)=0.0
       ! specific humidity needed (will be convert back in the end)
       UPQT(1,I)=qv(1)/(1.+qv(1))+0.32*UPW(1,I)*sigmaQT/sigmaW
       ! according to cheinet the 0.58 is for thetav, hence thetav is initialized (instead of theta)
       UPTHV(1,I)=thetav(1)+0.58*UPW(1,I)*sigmaTHV/sigmaW
       UPTHL(1,I)=UPTHV(1,I)/(1.+epsv*UPQT(1,I))                      ! liquid water pot. temp = pot. temp since no condensate at sfc yet
       UPT(1,I) = UPTHL(1,I) * (pres(1)/p00)**(rgas/cp)
    ENDDO

    

        !    write(*,*) ' z=',zi(1),' w=',upw(1,1),' thv=',upthv(1,1),' thl=',upthl(1,1),' qt=',upqt(1,1)*1000.,' qc=',upqc(1,1) *1000.
        !    write(*,*) 'qv=',qv(1)/(1.+qv(1))*1000. 
! updraft integration     
    DO i=1,nup
       DO k=2,nzm

          EntExp=exp(-ENT(k-1,i)*(zi(k)-zi(k-1)))

          ! WL assume environment unsaturated: use qv instead of qtotal for now
          ! WL FIX later once the prog. water variables are defined in Dycore
          QTn=qv(k-1)/(1.+qv(k-1))*(1.-EntExp)+UPQT(k-1,i)*EntExp
          !QTn=qv(k)/(1.+qv(k))*(1.-EntExp)+UPQT(k-1,i)*EntExp
          ! WL assuming for now that the environment holds no liquid or solid water: theta_l=theta
          THLn=thetav(k-1)/(1.+epsv*qv(k-1)/(1.+qv(k-1)))& 
          !THLn=thetav(k)/(1.+epsv*qv(k)/(1.+qv(k)))& 
               *(1.-EntExp)+UPTHL(k-1,i)*EntExp
          !Un=u(k)*(1.-EntExp)+UPU(k-1,i)*EntExp
          !Vn=v(k)*(1.-EntExp)+UPV(k-1,i)*EntExp
          ! compute condensation
          call condensation_edmf(QTn,THLn,presi(k),THVn,QCn)
          !THVn = UPTHL(k-1,i)*(1.+epsv*UPQT(k-1,i))
          !QCn = 0.
          B=ggr*(0.5*(THVn+UPTHV(k-1,i))/thetav(k-1)-1.)
          !B=ggr*(0.5*(THVn+UPTHV(k-1,i))/thetav(k)-1.)
          EntW=exp(-2.*(Wb+Wc*ENT(k-1,i))*(zi(k)-zi(k-1)))
          Wn2=UPW(k-1,i)**2*EntW + 2.*Wa* B *(zi(k)-zi(k-1))   !+(1.-EntW)*Wa*B/(Wb+Wc*ENT(k-1,i))
          IF (Wn2 >0.) THEN
            UPW(k,i)=sqrt(Wn2)
            UPTHV(k,i)=THVn
            UPT  (k,i)=THVn/(1.+epsv*(QTn-QCn))*(presi(k)/p00)**(rgas/cp)
            UPTHL(k,i)=THLn
            UPQT(k,i)=QTn
            UPQC(k,i)=QCn
          !  UPU(k,i)=Un
          !  UPV(k,i)=Vn
            UPA(k,i)=UPA(k-1,i)
         ELSE
             EXIT
         END IF
       ENDDO
    ENDDO

! computing variables needed for tendency calculation


! mass flux is zero at surface and top
    sumM(1)       = 0.0
    sumMthetav(1) = 0.0
    sumMrv(1)     = 0.0
    sumMu(1)      = 0.0
    sumMv(1)      = 0.0
    sumMtke(1)    = 0.0
    sumM(nz)       = 0.0
    sumMthetav(nz) = 0.0
    sumMrv(nz)     = 0.0
    sumMu(nz)      = 0.0
    sumMv(nz)      = 0.0
    sumMtke(nz)    = 0.0

    DO k=2,nzm
      DO i=1,nup
        sumM(k)      =sumM(k)      +UPA(K,I)*UPW(K,I)
        ! WL and convert back to thetav
        sumMthetav(k)=sumMthetav(k)+UPA(K,i)*UPW(K,I)*                   &
        (UPTHL(K,I)+fac_cond*UPQC(k,i))*(1.+epsv*(UPQT(k,i)-UPQC(k,i)))
        ! WL and convert back to mixing ratio
        sumMrv(k)    =sumMrv(k)    +UPA(K,i)*UPW(K,I)*                   & 
        (UPQT(k,i)-UPQC(k,i))/(1.-(UPQT(k,i)-UPQC(k,i)))
        !sumMu(k)  =sumMu(k)+UPA(K,i)*UPW(K,I)*UPU(K,I)
        !sumMv(k)  =sumMv(k)+UPA(K,i)*UPW(K,I)*UPV(K,I)
      ENDDO
    ENDDO
 

end subroutine edmf

subroutine condensation_edmf(QT,THL,P,THV,QC)
!
! zero or one condensation for edmf: calculates THV and QC
!
use params
use vars, only : qsatw

implicit none
real,intent(in) :: QT,THL,P
real,intent(out):: THV,QC

integer :: niter,i
real :: diff,exn,t,qs,qcold

! number of iterations
niter=100
! minimum difference
diff=1.e-7

EXN=(P/p00)**(rgas/cp)
QC=0.

!print*, '+++++++++++++++++++++++++'
!print*, '+++++++++++++qc=0++++++++'
!print*, '+++++++++++++++++++++++++'
do i=1,NITER
T=EXN*(THL+lcond/cp*QC)
! WL get saturation mixing ratio
QS=qsatw(T,P)
! WL convert to sat. spec. humidity
QS=QS/(QS+1.)
QCOLD=QC
QC=max(0.5*QC+0.5*(QT-QS),0.)
!write(*,'(5A)') ' tabs','  ','   qs','  ','   qc'
!write(*,'(F5.1,2x,F5.2,2x,F5.2)'),T,qs*1000.,qc*1000.
if (abs(QC-QCOLD)<Diff) exit
!print*, '+++++++++++++++++++++++++'
enddo

T=exn*(THL+lcond/cp*QC)
QS=qsatw(T,P)
QS=QS/(QS+1.)
QC=max(QT-QS,0.)
THV=(THL+lcond/cp*QC)*(1.+epsv*(QT-QC));

end subroutine condensation_edmf

