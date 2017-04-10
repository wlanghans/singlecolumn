module params

implicit none

! physical constants
real, parameter :: cp = 1004.             ! Specific heat of air, J/kg/K
real, parameter :: ggr = 9.81             ! Gravity acceleration, m/s2
real, parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
real, parameter :: lfus = 0.3336e+06      ! Latent heat of fusion, J/kg
real, parameter :: lsub = 2.8440e+06      ! Latent heat of sublimation, J/kg
real, parameter :: rv = 461.              ! Gas constant for water vapor, J/kg/K
real, parameter :: rgas = 287.            ! Gas constant for dry air, J/kg/K
real, parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
real, parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
real, parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

real, parameter :: fac_cond = lcond/cp
real, parameter :: fac_fus = lfus/cp
real, parameter :: fac_sub = lsub/cp

real, parameter ::  pi = 3.141592653589793
real, parameter ::  epsv = rv/rgas-1.
real, parameter ::  p00 = 1.e3          ! in hPa
real, parameter ::  salt_factor = 1.
real, parameter ::  xkar = 0.4
real, parameter ::  omega = 7.29212e-5


! implicit weight for diffusion
real ::  betam = 0.0
real ::  betap 

! model parameters
real::   fluxt0 =0.  ! surface sensible flux, Km/s
real::   fluxq0 =0.  ! surface latent flux, m/s
real::   tau0   =0.  ! surface stress, m2/s2
real::   z0     =0.0035  ! roughness length
logical:: sfc_flx_fxd =.true. ! surface sensible flux is fixed
logical:: doenergyunit=.false.
logical:: sfc_tau_fxd =.true.! surface drag is fixed
logical :: sfc_cm_fxd=.false.
logical :: sfc_cs_fxd=.false.
real :: Cs_in = 1.e-2
real :: Cm_in = 1.e-3
logical :: doneuman = .false.
logical :: dotkedirichlet = .false.
logical:: ocean =.true.! 
logical:: land =.false.! 
real:: windshear = 0.

integer :: nup = 1

logical:: dosgs = .true.
logical:: dosgscloud = .true.
logical :: dovartrans=.false.
real :: lcld = 250.
real :: ctketau=0.5
logical:: doedmf = .false.
real :: Wc=0.5
real :: pwmin=1.2
logical:: donoplumesat=.false.
logical:: fixedeps=.true.
logical:: neggerseps=.false.
logical:: gregoryeps=.false.
real:: tauneggers=500.
logical:: randomeps=.false.
logical :: fixedfa=.true.
real   :: eps0=1.e-03
real   :: del0=1.e-03
logical:: dosmagor = .true.
logical:: dopblh = .false.
real :: fixedpblh = -100.
logical:: pblhfluxmin=.true.
logical:: pblhthgrad=.false.
logical:: doconsttk = .false.
logical:: dosurface = .true.
real   :: tkconst = 10.
logical :: doteixpbl=.false.
logical :: fixedtau=.true.
logical :: progtke=.true.
logical :: dosequential=.false.
logical :: dotkeles=.false.
logical :: dowitekpbl=.false.
logical :: dolanghanspbl=.false.
logical :: dosingleplume=.false.
real :: beta = 0.3
logical :: witekeps = .false.
logical :: dosubsidence = .false.
logical :: docoriolis = .false.
logical :: dothuburn = .true.
logical :: doforcing=.false.
real :: fcor = 0.376d-04
logical :: dolongwave=.false.
logical :: doshortwave=.false.
logical :: doradsimple=.true.
logical :: doradhomo=.false.
logical :: doseasons=.false.
logical :: doperpetual=.true.
logical :: dosolarconstant=.true.
logical :: notracegases=.false.
real    :: solar_constant = 413.98 ! solar constant (in W/m2)
real    :: zenith_angle = 50.5   ! zenith angle (in degrees)
real    :: nxco2 = 1         ! factor to modify co2 concentration
real    :: coszrs
real:: longitude0 = 0.    ! latitude of the domain's center 
real:: latitude0  = 0.    ! longitude of the domain's center 
logical :: dotlflux = .true.
logical :: dozerosigma=.false.
logical :: donoenvcloud=.false.




!output
logical:: snapshot_do = .true.
integer:: snapshot_start = 0
integer:: snapshot_period =1
integer:: snapshot_end = 1000
logical:: snapshot_as_double = .false.
character(1000) :: snapshot_fields = 'hli,u,v,w,th,thv,tabs,tke,p,rho,qt,qv,qn,qcl,qci,&
		qpl,qpi,qtflx,tflx,totbuoyflx,tkewthv,crv,ctheta,cm,q1,sigmas,cfrac_pdf,&
		cthl,cqt,varwrt1,tk,tkh,lmix,tend_mix_qt,tend_mix_t,tend_mix_tke,tend_buoy_tke,&
		tend_shear_tke,tend_diss_tke,tkeflx,pblh,B,upw,upthd,upqcl,upqci,upqt,upthv,ent,wstar,ustar,tend_rad_t,&
                radlwdn,radlwup,radqrlw,radswdn,radswup,radqrsw,thetaligrad,qtgrad,thl,lwp,a_mf,&
                tkemf,qtflx_ed,qtflx_mf,tflx_ed,tflx_mf,tke_s,cfrac_tot,cfrac_mf'


end module params
