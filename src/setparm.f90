subroutine setparm

use vars
use params
use grid

implicit none

integer :: ierr

NAMELIST /PARAMETERS/dz, dt, doconstdz, dosgs, dosmagor, doedmf, dosurface, dolargescale, &
                     fluxt0,fluxq0,tau0,dt,nstop, dosmagor, doedmf, &
                     betam, betap, z0, sfc_flx_fxd, sfc_tau_fxd, & 
                     ocean, land, nup, dopblh, windshear, &
                     snapshot_do, snapshot_start, snapshot_period, snapshot_end, & 
                     snapshot_as_double, snapshot_fields, doconsttk, tkconst, sst, &
                     dolteix,pblhfluxmin,nzm, fixedtau

open(8,file='./CaseName',status='old',form='formatted')
read(8,'(a)') case
close (8)

open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')
read (55,PARAMETERS,IOSTAT=ierr)
if (ierr.ne.0) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification in PARAMETERS namelist'
        stop
end if
close(55)

path='./'//trim(case)//'/'

if (doconsttk) then
  write(*,*) 'doconsttk=.true., a constant eddy-diffusivity K=',tkconst,' m2/s is used'
end if


!  Ensure parameters are set consistently
if (dosgs.and.doedmf) then 
   dopblh = .true.
   pblhfluxmin = .true.
   write(*,*) 'Setting dopblh = .true. and pblhfluxmin = .true.'
end if
if (dopblh.and..not.dosgs) then
  write(*,*) 'Cant compute pblh without turbulence parameterization. Stopping'
  stop
end if
if (.not.dosurface) then
   write(*,*) 'dosurface = .false., thus zero all surface fluxes'
   write(*,*) 'fluxt0=0,fluxq0=0, and tau0=0'
   fluxt0=0.0
   fluxq0=0.0
   tau0=0.0
end if

if (dosmagor) then 
   if (dosgs) write(*,*) 'Smagorinski closure is used to get K'
   progtke = .false.
else
   if (dosgs) write(*,*) 'TKE closure is used to get K'
   progtke = .true.
end if

if (land) ocean=.false.
if ((betam+betap).ne.1.) betam=1.-betap
if (betam.gt.1.or.betam.lt.0.) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification of betap'
        stop
end if
if (snapshot_fields(1:1) .eq. '+') then
   snapshot_fields = 'u,v,th,thv,tabs,tke,rho,qv,qcl,qci,vaporflx,heatflx,virtheatflx,cthetav,crv,ctheta,cm,'&
                     //trim(snapshot_fields(2:len(snapshot_fields)))
end if

nz = nzm  + 1



! write namelist values out to file for documentation
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'.nml',&
            form='formatted', action='write')
      write (55,nml=PARAMETERS)
      write(55,*)
      close(55)


end subroutine setparm
