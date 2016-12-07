module snapshot

   use netcdf_mod, only: netcdf_data
   use params
   use grid
   
   implicit none

   type(netcdf_data) :: snap
   
   logical :: snapshot_dump
   
   public :: snapshot_init, snapshot_update, snapshot_close

   private :: snapshot_define
   
   contains   
   
      !============================!
      ! Subroutine snapshot_define !
      !============================!
      subroutine snapshot_define(ncsnap)
   
         use netcdf
         use netcdf_mod
         use params
         use vars
   
         implicit none
   
         ! In
         type(netcdf_data), target :: ncsnap
   
         integer :: zID, z2ID, npID, timeID
         integer :: num_fields, field_count, field
         integer :: k,i,j
         real, dimension(:), allocatable :: buffer
         character(50) :: field_name

         if (.not.snapshot_do) return
         
         ncsnap%title = 'Snapshot fields'
         ncsnap%file_path ='./'//trim(path)//'/'//trim(case)//'_snapshot.nc' 
   
         allocate( ncsnap%dimNames  (4), &
                   ncsnap%dimUnits  (4), &
                   ncsnap%dimLengths(4), &
                   ncsnap%dimIDs    (4), &
                   ncsnap%dimVarIDs (4), &
                   ncsnap%dimValues (3,max(nz,nup)) )
   
         ncsnap%dimNames = (/'z   ', 'z2  ', 'np  ', 'time' /)
         ncsnap%time = .true.   ! Time is one of the NetCDF file's dimensions
         ncsnap%zero_out = .false.   ! Do not set the data to zero after writing to the NetCDF file
         ncsnap%factor = 1.   ! Do not divide the data by any factor
                       !====================================================================!
                       ! The NetCDF routine assigns dimension ID's sequentially starting    !
         zID    = 1    ! with 1, so we know ahead of time that these will be the ID numbers !
         z2ID   = 2    !====================================================================!
         npID   = 3
         timeID = 4
         ncsnap%dimUnits = (/ 'm   ', 'm   ', '    ', 's   ' /)
         ncsnap%dimLengths = (/ nzm, nz, nup, nf90_unlimited /)
         ncsnap%dimValues(1,1:nzm)= (/(z(i) ,i=1,nzm)/)
         ncsnap%dimValues(2,1:nz) = (/(zi(i),i=1,nz)/)
         ncsnap%dimValues(3,1:nup) = (/(real(i),i=1,nup)/)
         ncsnap%timeCounter = 0
         if (snapshot_as_double) then
            ncsnap%as_double = .true.
         else
            ncsnap%as_double = .false.
         end if
         ncsnap%numWrites = 0
         
         field_count = 0
   
   
         !====================================================================================!
         ! Below are the variables to be included in the snapshot NetCDF file.                !
         ! Make sure that the number of fields written below matches num_fields.              !
         ! The dimensions argument is optional -- omit it if the field is zero dimensional.   !
         ! If the field is 1D, always associate or allocate it with dimensions (*:*,1:1,1:1). !
         ! If the field is 0D, always associate or allocate it with dimensions (1:1,1:1,1:1). !
         ! In the dimension ID list, the timeID must come last, if present, and the order     !
         !    of spatial dimensions must match that of the field.                             !
         !====================================================================================!
   

         num_fields = string_field_number(snapshot_fields)   ! This must match the number of fields below

         do field = 1,num_fields

            field_name = string_field(snapshot_fields,field)

            if (trim(field_name) .eq. 'u') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'u', &
                  ! Long name
                                           'Horizontal velocity along x-axis', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => u(1:nzm)
            end if

            if (trim(field_name) .eq. 'ug') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'ug', &
                  ! Long name
                                           'Geostrophic wind in x', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => ug(1:nzm)
            end if

            if (trim(field_name) .eq. 'vg') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'vg', &
                  ! Long name
                                           'Geostrophic wind in y', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => vg(1:nzm)
            end if
         
         
            if (trim(field_name) .eq. 'v') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'v', &
                  ! Long name
                                           'Horizontal velocity along y-axis', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => v(1:nzm)
            end if

            if (trim(field_name) .eq. 'tke') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tke', &
                  ! Long name
                                           'Turbulent kinetic energy; ED', &
                  ! Units
                                           'J/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tke(1:nzm)
            end if

            if (trim(field_name) .eq. 'tkemf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tkemf', &
                  ! Long name
                                           'Turbulent kinetic energy; MF', &
                  ! Units
                                           'J/kg', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tke_mf(1:nz)
            end if
         
            if (trim(field_name) .eq. 'thv') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'thv', &
                  ! Long name
                                           'Virtual Potential temperature', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => thetav(1:nzm)
            end if

            if (trim(field_name) .eq. 'thr') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'thr', &
                  ! Long name
                                           'Density Potential temperature', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => thetar(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'th') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'th', &
                  ! Long name
                                           'Potential temperature', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => theta(1:nzm)
            end if
            if (trim(field_name) .eq. 'thl') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'thl', &
                  ! Long name
                                           'Liquid water potential temperature', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => thetali(1:nzm)
            end if

            if (trim(field_name) .eq. 'rho') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'rho', &
                  ! Long name
                                           'Dry Air Density', &
                  ! Units
                                           'kg/m^3', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => rho(1:nzm)
            end if
         
         
            if (trim(field_name) .eq. 'tabs') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tabs', &
                  ! Long name
                                           'Temperature', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tabs(1:nzm)
            end if

            if (trim(field_name) .eq. 'hli') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'hli', &
                  ! Long name
                                           'Frozen/liquid moist static energy', &
                  ! Units
                                           'J/kg', &
                  ! Dimensions
                                           (/zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => t(1:nzm)
            end if

            if (trim(field_name) .eq. 'qcsgs_mf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qcsgs_mf', &
                  ! Long name
                                           'qc from plumes', &
                  ! Units
                                           'g/g', &
                  ! Dimensions
                                           (/z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qcsgs_mf(1:nz)
            end if

            if (trim(field_name) .eq. 'qcsgs_pdf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qcsgs_pdf', &
                  ! Long name
                                           'qc from PDF scheme', &
                  ! Units
                                           'g/g', &
                  ! Dimensions
                                           (/zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qcsgs_pdf(1:nzm)
            end if

            
         
            if (trim(field_name) .eq. 'tend_force_u') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_force_u', &
                  ! Long name
                                           'Coriolis forcing u', &
                  ! Units
                                           'kg/m3 m/s2', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_force_u(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'tend_force_v') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_force_v', &
                  ! Long name
                                           'Coriolis forcing v', &
                  ! Units
                                           'kg/m3 m/s2', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_force_v(1:nzm)
            end if

            if (trim(field_name) .eq. 'p') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'p', &
                  ! Long name
                                           'Pressure', &
                  ! Units
                                           'hPa', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => pres(1:nzm)
            end if
         
         
            if (trim(field_name) .eq. 'qv') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qv', &
                  ! Long name
                                           'Water vapor mixing ratio', &
                  ! Units
                                           'kg/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qv(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'qci') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qci', &
                  ! Long name
                                           'Non-precip cloud ice mixing ratio', &
                  ! Units
                                           'kg/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qci(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'qt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qt', &
                  ! Long name
                                           'Total water mixing ratio', &
                  ! Units
                                           'kg/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qt(1:nzm)
            end if

            if (trim(field_name) .eq. 'qn') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qn', &
                  ! Long name
                                           'Non-precip condensate mixing ratio', &
                  ! Units
                                           'kg/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qn(1:nzm)
            end if

            if (trim(field_name) .eq. 'qcl') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qcl', &
                  ! Long name
                                           'Non-precip cloud liquid mass faction', &
                  ! Units
                                           'kg/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qcl(1:nzm)
            end if

            if (trim(field_name) .eq. 'qpi') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qpi', &
                  ! Long name
                                           'Precip cloud ice mixing ratio', &
                  ! Units
                                           'kg/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qpi(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'qpl') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qpl', &
                  ! Long name
                                           'Precip cloud liquid mass faction', &
                  ! Units
                                           'kg/kg', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qpl(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'lmix') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'lmix', &
                  ! Long name
                                           'mixing length', &
                  ! Units
                                           'm', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => smix(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'presi') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'pi', &
                  ! Long name
                                           'Pressure on interfaces', &
                  ! Units
                                           'hPa', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => presi(1:nz)
            end if
         
            if (trim(field_name) .eq. 'tkh') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tkh', &
                  ! Long name
                                           'Eddy diffusivity for heat', &
                  ! Units
                                           'm2/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tkh(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'tk') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tk', &
                  ! Long name
                                           'Eddy diffusivity for momentum', &
                  ! Units
                                           'm2/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tk(1:nzm)
            end if
            
            if (trim(field_name) .eq. 'tend_rad_t') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_rad_t', &
                  ! Long name
                                           'hli tendency from radiation', &
                  ! Units
                                           'J/kg/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_rad_t(1:nzm)
            end if

            if (trim(field_name) .eq. 'tend_mix_tke') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_mix_tke', &
                  ! Long name
                                           'tke tendency from sgs mixing', &
                  ! Units
                                           'm2/s3', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_sgs_tke(1:nzm)
            end if

            if (trim(field_name) .eq. 'tend_shear_tke') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_shear_tke', &
                  ! Long name
                                           'tke tendency from shear production', &
                  ! Units
                                           'm2/s3', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_tke_shear(1:nzm)
            end if

            if (trim(field_name) .eq. 'tend_mix_qt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_mix_qt', &
                  ! Long name
                                           'Tot water mix ratio tend from mixing', &
                  ! Units
                                           '1/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_sgs_qt(1:nzm)
            end if

            if (trim(field_name) .eq. 'tend_buoy_tke') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_buoy_tke', &
                  ! Long name
                                           'tke tendency from buoy production', &
                  ! Units
                                           'm2/s3', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_tke_buoy(1:nzm)
            end if

            if (trim(field_name) .eq. 'tend_diss_tke') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_diss_tke', &
                  ! Long name
                                           'tke tendency from dissipation', &
                  ! Units
                                           'm2/s3', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_tke_diss(1:nzm)
            end if
         
         
         
            if (trim(field_name) .eq. 'tend_mix_t') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tend_mix_t', &
                  ! Long name
                                           'hli tendency from mixing', &
                  ! Units
                                           'J/kg/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tend_sgs_t(1:nzm)
            end if

            if (trim(field_name) .eq. 'n2') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'n2', &
                  ! Long name
                                           'Squared Brunt-Vaisala frequency', &
                  ! Units
                                           '1/s2', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => buoy_sgs(1:nzm)
            end if

            if (trim(field_name) .eq. 'dthetadt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'dtdt', &
                  ! Long name
                                           'large scale theta forcing', &
                  ! Units
                                           'K/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => dthetadt(1:nzm)
            end if
           
            if (trim(field_name) .eq. 'dqtdt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'dqtdt', &
                  ! Long name
                                           'large scale qt forcing', &
                  ! Units
                                           'g/g/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => dqtdt(1:nzm)
            end if
           
            if (trim(field_name) .eq. 'def2') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'def2', &
                  ! Long name
                                           'Squared deformation rate', &
                  ! Units
                                           '1/s2', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => def2(1:nzm)
            end if
           
            if (trim(field_name) .eq. 'sigmas') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'sigmas', &
                  ! Long name
                                           'std deviation of sat deficit', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sigmas(1:nzm)
            end if

            if (trim(field_name) .eq. 'varwrt1') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'varwrt1', &
                  ! Long name
                                           'additional variable for quick checks', &
                  ! Units
                                           'variable', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => varwrt1(1:nzm)
            end if

            if (trim(field_name) .eq. 'q1') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'q1', &
                  ! Long name
                                           'saturation deficit', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => q1(1:nzm)
            end if

            if (trim(field_name) .eq. 'cthl') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'cthl', &
                  ! Long name
                                           'buoyancy flux coeff. for wthl flux', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => cthl(1:nzm)
            end if

            if (trim(field_name) .eq. 'cqt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'cqt', &
                  ! Long name
                                           'buoyancy flux coeff. for wqt flux', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => cqt(1:nzm)
            end if

            if (trim(field_name) .eq. 'tkewthv') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tkewthv', &
                  ! Long name
                                           'kin. buoy flux used in TKE procution term', &
                  ! Units
                                           'K m/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tke_thvflx(1:nzm)
            end if
            if (trim(field_name) .eq. 'thetaligrad') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'thetaligrad', &
                  ! Long name
                                           'dthetali/dz', &
                  ! Units
                                           'K/m', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => thetaligrad(1:nzm)
            end if
            if (trim(field_name) .eq. 'qtgrad') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qtgrad', &
                  ! Long name
                                           'dqt/dz', &
                  ! Units
                                           'g/g/m', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qtgrad(1:nzm)
            end if


            if (trim(field_name) .eq. 'radqrsw') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'radqrsw', &
                  ! Long name
                                           'SW radiative heating rate', &
                  ! Units
                                           'K/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => radqrsw (1:nzm)
            end if
            if (trim(field_name) .eq. 'radswdn') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'radswdn', &
                  ! Long name
                                           'SW radiative flux down', &
                  ! Units
                                           'W/m2', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => radswdn (1:nz)
            end if
            if (trim(field_name) .eq. 'radswup') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'radswup', &
                  ! Long name
                                           'SW radiative flux up', &
                  ! Units
                                           'W/m2', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => radswup (1:nz)
            end if

            if (trim(field_name) .eq. 'radqrlw') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'radqrlw', &
                  ! Long name
                                           'LW radiative heating rate', &
                  ! Units
                                           'K/s', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => radqrlw (1:nzm)
            end if
            if (trim(field_name) .eq. 'radlwdn') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'radlwdn', &
                  ! Long name
                                           'LW radiative flux down', &
                  ! Units
                                           'W/m2', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => radlwdn (1:nz)
            end if
            if (trim(field_name) .eq. 'radlwup') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'radlwup', &
                  ! Long name
                                           'LW radiative flux up', &
                  ! Units
                                           'W/m2', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => radlwup (1:nz)
            end if
    
            if (trim(field_name) .eq. 'a_mf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'a_mf', &
                  ! Long name
                                           'Fractional area with plumes', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => frac_mf (1:nz)
            end if

            if (trim(field_name) .eq. 'cfrac_mf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'cfrac_mf', &
                  ! Long name
                                           'Cloud fraction from MF', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => cfrac_mf(1:nz)
            end if

            if (trim(field_name) .eq. 'cfrac_pdf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'cfrac_pdf', &
                  ! Long name
                                           'Cloud fraction in plume environment', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => cfrac_pdf (1:nzm)
            end if

            if (trim(field_name) .eq. 'cfrac_tot') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'cfrac_tot', &
                  ! Long name
                                           'Total cloud fraction (environ+plumes)', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ zID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => cfrac_tot(1:nzm)
            end if
         
            if (trim(field_name) .eq. 'qtflx') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qtflx', &
                  ! Long name
                                           'Total water mass flux; rho*L*w*qt', &
                  ! Units
                                           'W/m2', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sgs_qt_flux (1:nz)
            end if

            if (trim(field_name) .eq. 'tflx_ed') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tflx_ed', &
                  ! Long name
                                           'Thetal flux; ED', &
                  ! Units
                                           'K m/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => t_flux_ed (1:nz)
            end if

            if (trim(field_name) .eq. 'tflx_mf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tflx_mf', &
                  ! Long name
                                           'Thetal flux; MF', &
                  ! Units
                                           'K m/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => t_flux_mf (1:nz)
            end if

            if (trim(field_name) .eq. 'qtflx_ed') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qtflx_ed', &
                  ! Long name
                                           'Water mass flux; ED', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qt_flux_ed (1:nz)
            end if

            if (trim(field_name) .eq. 'qtflx_mf') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'qtflx_mf', &
                  ! Long name
                                           'Water mass flux; MF', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => qt_flux_mf (1:nz)
            end if
         
         
            if (trim(field_name) .eq. 'tflx') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tflx', &
                  ! Long name
                                           'Total thetal flux; rho*cp*w*thetal', &
                  ! Units
                                           'W/m2', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sgs_t_flux (1:nz)
            end if

            if (trim(field_name) .eq. 'taux') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'taux', &
                  ! Long name
                                           'x-momentumflux', &
                  ! Units
                                           'm/s m/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => taux (1:nz)
            end if

            if (trim(field_name) .eq. 'tauy') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tauy', &
                  ! Long name
                                           'y-momentumflux', &
                  ! Units
                                           'm/s m/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => tauy (1:nz)
            end if


            if (trim(field_name) .eq. 'tkeflx') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tkeflux', &
                  ! Long name
                                           'TKE flux', &
                  ! Units
                                           'm2/s2 m/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sgs_tke_flux (1:nz)
            end if

            if (trim(field_name) .eq. 'totbuoyflx') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'totbuoyflx', &
                  ! Long name
                                           'Total (ED+MF) buoyancy flux;rho*cp*w*thetav', &
                  ! Units
                                           'W/m2', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sgs_thv_flux (1:nz)
            end if
 
            if (trim(field_name) .eq. 'sumM') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'plume mass flux', &
                  ! Long name
                                           'sum of mass flux of all plumes', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sumM (1:nz)
            end if

            if (trim(field_name) .eq. 'sumMrt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'plume total moisture flux', &
                  ! Long name
                                           'sum of total moisture flux of all plumes', &
                  ! Units
                                           'g/g m/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sumMrt (1:nz)
            end if

            if (trim(field_name) .eq. 'sumMt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'plume hli flux', &
                  ! Long name
                                           'sum of hli fluxes of all plumes', &
                  ! Units
                                           'J/kg m/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => sumMt (1:nz)
            end if

            if (trim(field_name) .eq. 'upm') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upm', &
                  ! Long name
                                           'plume mass flux', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPM (1:nz,1:nup)
            end if

            if (trim(field_name) .eq. 'upw') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upw', &
                  ! Long name
                                           'plume velocity', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPW (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'ent') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'ent', &
                  ! Long name
                                           'plume frac. entrain. rate', &
                  ! Units
                                           '1/m', &
                  ! Dimensions
                                           (/ zID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => ENT (1:nzm,1:nup)
            end if
            if (trim(field_name) .eq. 'det') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'det', &
                  ! Long name
                                           'plume frac. detrain. rate', &
                  ! Units
                                           '1/m', &
                  ! Dimensions
                                           (/ zID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => DET (1:nzm,1:nup)
            end if

            if (trim(field_name) .eq. 'B') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'B', &
                  ! Long name
                                           'plume buoyancy', &
                  ! Units
                                           'm/s2', &
                  ! Dimensions
                                           (/ zID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => BUOY (1:nzm,1:nup)
            end if

            if (trim(field_name) .eq. 'upthd') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upthp', &
                  ! Long name
                                           'plume theta excess', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ zID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPTHD (1:nzm,1:nup)
            end if
 
            if (trim(field_name) .eq. 'w') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'w', &
                  ! Long name
                                           'subsidence velocity', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, timeID /) &
               )
               ncsnap%field(field_count)%data1 => w (1:nz)
            end if


            if (trim(field_name) .eq. 'uptabs') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upt', &
                  ! Long name
                                           'plume abs. temp.', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPTABS (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upthv') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upthv', &
                  ! Long name
                                           'plume virt. pot. temp.', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPTHV (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upa') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upa', &
                  ! Long name
                                           'plume area fraction', &
                  ! Units
                                           '', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPA (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upu') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upu', &
                  ! Long name
                                           'plume u velocity', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPU (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upv') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upv', &
                  ! Long name
                                           'plume v velocity', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPV (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upqcl') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upqcl', &
                  ! Long name
                                           'plume liquid water mass fraction', &
                  ! Units
                                           'g/g', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPQCL (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upqci') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upqci', &
                  ! Long name
                                           'plume frozen water mass fraction', &
                  ! Units
                                           'g/g', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPQCI (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upqt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upqt', &
                  ! Long name
                                           'plume total water mass fraction', &
                  ! Units
                                           'g/g', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPQT (1:nz,1:nup)
            end if
            if (trim(field_name) .eq. 'upt') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'upt', &
                  ! Long name
                                           'plume hli', &
                  ! Units
                                           'K', &
                  ! Dimensions
                                           (/ z2ID, npID, timeID /) &
               )
               ncsnap%field(field_count)%data2 => UPT (1:nz,1:nup)
            end if

            if (trim(field_name) .eq. 'ustar') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'ustar', &
                  ! Long name
                                           'characteristic u scale', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => ustar
            end if

            if (trim(field_name) .eq. 'wstar') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'wstar', &
                  ! Long name
                                           'characteristic w scale', &
                  ! Units
                                           'm/s', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => wstar
            end if

            if (trim(field_name) .eq. 'lwp') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'lwp', &
                  ! Long name
                                           'Liquid water path', &
                  ! Units
                                           'g/m2', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => lwp
            end if

            if (trim(field_name) .eq. 'tke_s') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'tke_s', &
                  ! Long name
                                           'TKE at surface', &
                  ! Units
                                           'm2/s2', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => tke_s
            end if

            if (trim(field_name) .eq. 'pblh') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'pblh', &
                  ! Long name
                                           'Boundary layer height', &
                  ! Units
                                           'm', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => pblh
            end if

            if (trim(field_name) .eq. 'ctheta') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'ctheta', &
                  ! Long name
                                           'Transfer coeff sensible heat flux', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => Ctheta
            end if


            if (trim(field_name) .eq. 'crv') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'crv', &
                  ! Long name
                                           'Transfer coeff latent heat flux', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => Crv
            end if

            if (trim(field_name) .eq. 'cm') then
               call fill_fields_snapshot( &
                  ! Short name
                                           'cm', &
                  ! Long name
                                           'Drag coefficient', &
                  ! Units
                                           '1', &
                  ! Dimensions
                                           (/ timeID /) &
               )
               ncsnap%field(field_count)%data0 => Cm
            end if


                   
         end do   ! Loop over field
   
         if (field_count .lt. num_fields) then
            print *,'Error: Number of fields in snapshot_define is ',field_count,'&
            , which is less than the expected ',num_fields,' fields'
            stop
         end if
      
         ! Set ncsnap%size_as_float
         call netcdf_set_size_as_float(ncsnap)

         write(*,'(a,i4)') ' Number of snapshot fields        = ',field_count
   
         contains
   
            !=================================!
            ! Subroutine fill_fields_snapshot !
            !=================================!
            subroutine fill_fields_snapshot(short_name, long_name, units, dims)
      
      
               implicit none
      
               ! In
               character(*) short_name
               character(*) long_name
               character(*) units
               integer, dimension(:), optional :: dims
      
               integer, dimension(:), allocatable :: locs   ! Holds the spatial ID's in dims
               integer :: i
      
               ! Update field counter
               field_count = field_count + 1
      
               ! If first call, allocate the appropiate number of fields
               if (field_count .eq. 1) then
                  allocate( ncsnap%field(num_fields) )
                  do i = 1,num_fields
                     nullify(ncsnap%field(i)%data2)
                     nullify(ncsnap%field(i)%data1)
                     nullify(ncsnap%field(i)%data0)
                  end do
               end if

               ! Check that there are not more than num_fields fields listed in snapshot_define
               if (field_count .gt. num_fields) then
                  print *,'Error: More than ',num_fields,' fields listed in snapshot_define'
                  stop
               end if
      
               ! Assign strings to the field
               ncsnap%field(field_count)%short_name = short_name
               ncsnap%field(field_count)%long_name = long_name
               ncsnap%field(field_count)%units = units
      
               if (present(dims)) then   ! Otherwise, leave dimIDs and location unallocated
      
                  allocate( ncsnap%field(field_count)%dimIDs(size(dims)) )
                  ncsnap%field(field_count)%dimIDs = dims
      
                  if (any(dims .eq. timeID)) then   ! If dims contains time, give location one less entry ...
      
                     ! Set the time flag to true
                     ncsnap%field(field_count)%time = .true.
      
                     ! Check that the timeID is the last in entry in dims
                     if (dims(size(dims)) .ne. timeID) then
                        print *,'Error: The timeID must come last in snapshot_define'
                        stop
                     end if
      
                     ! Set the size of locs, which holds spatial dimension IDs, to one less than dims
                     if (size(dims) .gt. 1) then
                        allocate(locs(size(dims)-1))
                        locs = dims(1:(size(dims)-1))
                     end if
      
                  else                              ! ... otherwise, give location the same size as dims
      
                     ! Set the time flag to false
                     ncsnap%field(field_count)%time = .false.
      
                     ! Set the size of locs to that of dims since dims only holds spatial dimension IDs
                     allocate(locs(size(dims)))
                     locs = dims
      
                  end if
      
               end if
      
      
               if (allocated(locs)) then   ! Assign spatial start values in location, otherwise leave location unallocated
      
                  allocate( ncsnap%field(field_count)%location(size(locs)) )
      
                  ! Assign start values according to dimension and subdomain
                  do i = 1,size(locs)
                     if (locs(i) .eq.  npID) then
                        ncsnap%field(field_count)%location(i) = 1
                     else if (locs(i) .eq.  zID) then
                        ncsnap%field(field_count)%location(i) = 1
                     else if (locs(i) .eq. z2ID) then
                        ncsnap%field(field_count)%location(i) = 1
                     else
                        write (*,*) 'Error: Invalid dimension ID in fill_fields_snapshot'
                        stop
                     end if
                  end do
      
               end if
      
               if (allocated(locs)) deallocate(locs)
      
            end subroutine fill_fields_snapshot

   
            !========================================================!
            ! Function string_field_number                           !
            !                                                        !
            ! This function returns the number of fields in string   !
            ! using sep (defaulted to ',') as the field separator.   !
            !========================================================!
            function string_field_number(string,sep)
         
               implicit none
         
               ! In
               character(*), intent(in) :: string
               character(*), intent(in), optional :: sep
         
               ! Out
               integer :: string_field_number
         
               character(100) :: sep2
               integer :: string_pos, string_field_pos, field_num
         
               sep2 = ','
               if (present(sep)) sep2 = sep
         
               string_field_number = 1
               do string_pos = 1,len(string)
                  if (string(string_pos:string_pos) .eq. sep2) then
                     string_field_number = string_field_number + 1
                  end if
               end do
         
            end function string_field_number


            !========================================================!
            ! Function string_field                                  !
            !                                                        !
            ! This function returns the field'th substring in string !
            ! using sep (defaulted to ',') as the field separator.   !
            !========================================================!
            function string_field(string,field,sep)
         
               implicit none
         
               ! In
               character(*), intent(in) :: string
               integer, intent(in) :: field
               character(*), intent(in), optional :: sep
         
               ! Out
               character(len(string)) :: string_field
         
               character(100) :: sep2
               integer :: string_pos, string_field_pos, field_num
         
               sep2 = ','
               if (present(sep)) sep2 = sep
         
               string_field = ''
         
               string_pos = 1
               string_field_pos = 1
               field_num = 1
               do while (string_pos .le. len(string) .and. field_num .le. field)
                  if (string(string_pos:string_pos) .eq. sep2) then
                     field_num = field_num + 1
                  else
                     if (field_num .eq. field) then
                        string_field(string_field_pos:string_field_pos) = string(string_pos:string_pos)
                        string_field_pos = string_field_pos + 1
                     end if
                  end if
                  string_pos = string_pos + 1
               end do
               
               if (string_field_pos .eq. 1) then
                   write(*,*) 'Error in string_field: There is no such field number'
                   stop
               end if
         
            end function string_field
   
      end subroutine snapshot_define

   
      !==========================!
      ! Subroutine snapshot_init !
      !==========================!
      subroutine snapshot_init()
      
         use grid
         use netcdf_mod
      
         implicit none

         if (snapshot_do) then
            call snapshot_define(snap)
            call netcdf_open(snap)
            if (snapshot_start .eq. 0) call netcdf_write_allproc(snap)
         end if 
      
      end subroutine snapshot_init    


      !============================!
      ! Subroutine snapshot_update !
      !============================!
      subroutine snapshot_update()

         use netcdf_mod, only: netcdf_write_allproc
         use vars
         use params
         use grid
      
         implicit none
         
         integer :: i, j, k

         snapshot_dump = .false.

         if ( snapshot_do               .and. &
              nstep .ge. snapshot_start .and. &
              nstep .le. snapshot_end   .and. &
              mod(nstep-snapshot_start,snapshot_period) .eq. 0 ) then
              call netcdf_write_allproc(snap)
              snapshot_dump = .true.
         end if

      end subroutine snapshot_update


      !===========================!
      ! Subroutine snapshot_close !
      !===========================!
      subroutine snapshot_close()

         use netcdf_mod, only: netcdf_close

         implicit none

         if (snapshot_do) call netcdf_close(snap)

      end subroutine snapshot_close

   
end module snapshot
