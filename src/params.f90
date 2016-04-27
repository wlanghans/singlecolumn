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


! implicit weight for diffusion
real ::  betam = 0.0
real ::  betap 

! model parameters
real::   fluxt0 =0.  ! surface sensible flux, Km/s
real::   fluxq0 =0.  ! surface latent flux, m/s
real::   tau0   =0.  ! surface stress, m2/s2
real::   z0     =0.0035  ! roughness length
logical:: sfc_flx_fxd =.true. ! surface sensible flux is fixed
logical:: sfc_tau_fxd =.true.! surface drag is fixed
logical :: doneuman = .false.
logical:: ocean =.true.! 
logical:: land =.false.! 
real:: windshear = 0.

integer :: nup = 1

logical:: dosgs = .true.
logical:: doedmf = .false.
logical:: dosmagor = .true.
logical:: dopblh = .false.
logical:: pblhfluxmin=.true.
logical:: doconsttk = .false.
logical:: dosurface = .true.
logical:: dolargescale = .false.
real   :: tkconst = 10.
logical :: dolteix=.false.
logical :: fixedtau=.true.
logical :: progtke=.true.


!output
logical:: snapshot_do = .true.
integer:: snapshot_start = 0
integer:: snapshot_period =1
integer:: snapshot_end = 1000
logical:: snapshot_as_double = .false.
character(200) :: snapshot_fields = 'u,v,th,thv,tabs,tke,rho,qv,qcl,qci,vaporflx,heatflx,virtheatflx,cthetav,crv,ctheta,cm'


end module params
