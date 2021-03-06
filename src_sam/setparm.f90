subroutine setparm

use vars
use params
use grid

implicit none

integer :: ierr

NAMELIST /PARAMETERS/dz, dt, doconstdz, dosgs, dosmagor, doedmf, donoscale, doenvedmf,dosurface, &
                     fluxt0,fluxq0,tau0,nstop, &
                     betam, betap, z0, sfc_flx_fxd, sfc_tau_fxd, & 
                     ocean, land, nup, dopblh, windshear, &
                     snapshot_do, snapshot_start, snapshot_period, snapshot_end, & 
                     snapshot_as_double, snapshot_fields, doconsttk, tkconst, sst, &
                     doteixpbl,dowitekpbl,dolanghanspbl,pblhfluxmin,nzm, fixedtau, doneuman, fixedeps, eps0, Cs_in, Cm_in,&
                     sfc_cs_fxd, sfc_cm_fxd, neggerseps,gregoryeps, randomeps, del0, fixedfa, dosequential, &
                     dotkeles,dosingleplume,pblhthgrad,witekeps,dosgscloud, fcor, dosubsidence, docoriolis, dothuburn, doforcing, &
                     doshortwave, dolongwave, doradsimple, dotkedirichlet, donoplumesat, beta, tauneggers, fixedpblh, lcld, &
                     doenergyunit, dovartrans, dotlflux,pwmin,ctketau, neggersslooz, dozerosigma, Wc, donoenvcloud,alphaqt, alphathv

call getarg(1,path)
path =trim(adjustl(path))


open(55,file='./'//trim(path)//'/prm', status='old',form='formatted')
read (55,PARAMETERS,IOSTAT=ierr)
if (ierr.ne.0) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification in PARAMETERS namelist'
        stop
end if
close(55)

open(8,file='./'//trim(path)//'/CaseName',status='old',form='formatted')
read(8,'(a)') case
close (8)

write(*,*) 'Working on Case = ', trim(case) 

if (dosgs.and.count((/dotkeles,doteixpbl,dowitekpbl,dolanghanspbl,dosmagor,doconsttk/)).gt.1  ) then
  write(*,*) 'ERROR: only one eddy-diffusivity closure can be chosen'
  stop
end if
if (dosgs.and.doedmf.and.count((/witekeps,fixedeps,neggerseps,gregoryeps,randomeps/)).gt.1  ) then
  write(*,*) 'ERROR: only one entrainment option can be chosen.'
  stop
end if


if (dosgs.and.doedmf.and.randomeps ) then
  write(*,*) 'WARNING: random entrainment not implemented yet'
  write(*,*) 'Using Neggers entrainment eps~1/w instead'
  randomeps =.false.
  neggerseps=.true.
  gregoryeps=.false.
  fixedeps  = .false.
end if


if (doedmf.and.dosgs.and.fixedeps.and.eps0.lt.0.) then
   write(*,*) 'ERROR: eps0 has to be non-negative'
   stop
end if

if (doedmf.and.dosgs.and.del0.lt.0.) then
   write(*,*) 'ERROR: del0 has to be non-negative'
   stop
end if


!  Ensure parameters are set consistently
if ((betam+betap).ne.1.) betam=1.-betap
if (betam.gt.1.or.betam.lt.0.) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification of betap'
        stop
end if
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*) 'Solver for k>1: '
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*) 'Implicit weight is beta+=',betap
write(*,*) '  '
if (sfc_flx_fxd.and.sfc_cs_fxd) then
   write(*,*) 'Since sfc_flx is fixed, sfc_cs_fxd is set to .false.'
   write(*,*) '  '
   sfc_cs_fxd=.false.
end if
if (sfc_tau_fxd.and.sfc_cm_fxd) then
   write(*,*) 'Since tau is fixed, sfc_cm_fxd is set to .false.'
   write(*,*) '  '
   sfc_cm_fxd=.false.
end if
if (.not.dosurface) then
   write(*,*) 'dosurface = .false., thus zero all surface fluxes'
   write(*,*) 'fluxt0=0,fluxq0=0, and tau0=0'
   write(*,*) '  '
   fluxt0=0.0
   fluxq0=0.0
   tau0=0.0
   sfc_flx_fxd=.true.
   sfc_tau_fxd=.true.
end if
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*) 'Solver for k=1: '
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
if (sfc_flx_fxd.or.doneuman) then
   write(*,*) 'Neuman BC is used for sgs-tendencies of scalars'
   if (sfc_flx_fxd) then
     write(*,*) 'Explicit solve with fixed fluxes: ', fluxt0,' K m/s  and ', fluxq0, ' g/g m/s' 
   else
     write(*,*) 'Implicit weight is beta+=',betap
   end if
else
    write(*,*) 'Dirichlet BC is used for sgs-tendencies of scalars'
    write(*,*) 'Implicit weight is beta+=',betap
end if
if (sfc_tau_fxd.or.doneuman) then
   write(*,*) 'Neuman BC is used for sgs-tendencies of momentum'
   if (sfc_tau_fxd) then
     write(*,*) 'Explicit solve with fixed tau = ', tau0,' (m/s)**2 '
   else
     write(*,*) 'Implicit weight is beta+=',betap
   end if
else
    write(*,*) 'Dirichlet BC is used for sgs-tendencies of momentum'
    write(*,*) 'Implicit weight is beta+=',betap
end if

if (sfc_cs_fxd) then
   write(*,*) 'Cs is fixed at ', cs_in,', T_s =  ',sst+t00
end if
if (sfc_cm_fxd) then
   write(*,*) 'Cm is fixed at ', cm_in
end if

write(*,*) '  '
write(*,*) 'General other settings:  '
if (dosgs.and.doedmf) then 
   dopblh = .true.
end if
if (dopblh.and..not.dosgs) then
  write(*,*) 'Cant compute pblh without turbulence parameterization. Stopping'
  stop
end if
if (dosgs.and.(doteixpbl.or.dolanghanspbl).and..not.fixedtau) then
   dopblh=.true.
   !pblhfluxmin=.true.
   !pblhthgrad=.false.
   !write(*,*) 'Setting dopblh = .true. and pblhfluxmin = .true.'
end if
if (dosgs.and.dowitekpbl) then
  if (doedmf) then
    !dosingleplume=.true.
    fixedfa=.true.
  end if
  dopblh=.true.
  pblhthgrad=.true.
  pblhfluxmin=.false.
end if
if (doedmf.and.(.not.dosingleplume).and.(pwmin.lt.0..or.pwmin.gt.3.5)) then
  write(*,*) 'ERROR: pwmin has to be larger than 0 and smaller than 3.5'
  stop
end if
if (doconsttk) then
  if (dosgs) write(*,*) 'doconsttk=.true., a constant eddy-diffusivity K=',tkconst,' m2/s is used'
  progtke = .false.
elseif (dosmagor) then
  if (dosgs) write(*,*) 'Smagorinski closure is used to get K'
  progtke = .false.
elseif (dotkeles) then
  if (dosgs) write(*,*) 'TKE closure is used to get K'
  progtke = .true.
elseif (doteixpbl) then
  if (dosgs) write(*,*) 'Teixeira PBL TKE closure is used to get K'
  progtke = .true.
elseif (dowitekpbl) then
  if (dosgs) write(*,*) 'Witek PBL TKE closure is used to get K'
  progtke = .true.
elseif (dolanghanspbl) then
  if (dosgs) write(*,*) 'Langhans PBL TKE closure is used to get K '
  progtke = .true.
end if

if (land.eqv.ocean) then
  write(*,*) 'WARNING: locean == lland, thus setting ocean=true and land=false'
  ocean=.true.
  land =.false.
end if
if (snapshot_fields(1:1) .eq. '+') then
   snapshot_fields = 'hli,u,v,w,th,thv,tabs,tke,rho,p,qt,qv,qn,qcl,qci,&
                                qpl,qpi,qtflx,tflx,totbuoyflx,tkewthv,crv,ctheta,cm,q1,sigmas,cfrac_pdf,&
                                cthl,cqt,varwrt1,tk,tkh,lmix,tend_mix_qt,tend_mix_t,tend_mix_tke, &
                                tend_buoy_tke,tend_shear_tke,tend_diss_tke,pblh,B,upw,upthd,upqcl,upqci,upqt,upthv,ent,wstar,ustar,&
                                tend_rad_t,radlwdn,radlwup,radqrlw,radswdn,radswup,radqrsw,thetaligrad,&
                                qtgrad,thl,lwp,a_mf,tke_mf,qtflx_ed,qtflx_mf,tflx_ed,tflx_mf,tke_s,cfrac_tot,cfrac_mf,'&
                     //trim(snapshot_fields(2:len(snapshot_fields)))
end if

nz = nzm  + 1



! write namelist values out to file for documentation
      open(unit=55,file='./'//trim(path)//'/'//trim(case)//'.prm',&
            form='formatted', action='write')
      write (55,nml=PARAMETERS)
      write(55,*)
      close(55)


end subroutine setparm
