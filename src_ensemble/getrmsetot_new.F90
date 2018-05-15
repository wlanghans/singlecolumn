Program Parallel_Get_Rmse
  use mpi
  use netcdf
  USE IFPORT
  implicit none

  character(300) :: param, path_in,pathin,path_ncl
  character(60) :: runname
  character(200) :: namehlp,command, filein,fname_out

  integer :: my_id, Nproc, ierr, Nsetup, N_local_runs, N_max, ncount
  integer, allocatable, dimension(:) :: mpi_id, Nlocalvec,displace

  integer :: i,j,k,l,n,ca,m,o
  integer :: itmp(10), Ncombo
  integer :: run_start, run_end, status
  character(2) :: i_str,j_str,k_str,l_str,m_str,o_str
  real :: readparam(8),readparam2(8)
  integer, dimension(:), allocatable :: Recv_Request
  integer, dimension(:,:), allocatable :: mpi_status
  integer, dimension(:), allocatable :: mpi_status1
  integer ::  Send_Request


  integer,parameter :: Nparam=6
  !integer,parameter :: Ncases=1
  !character(10),dimension(Ncases) :: casename=(/"BOMEX"/)
  integer,parameter :: Ncases=2
  character(10),dimension(Ncases) :: casename=(/"BOMEX","CPBL2"/)
  integer,dimension(Nparam) :: Nrange
  real, allocatable,dimension(:,:,:,:,:,:) :: rmse
  real, allocatable,dimension(:,:) :: pvalue
  real, allocatable,dimension(:,:) :: rmse1Dloc,rmse1Dloc2
  real, allocatable,dimension(:,:,:) :: rmse1Dlocind
  real, allocatable,dimension(:,:,:) :: rmse1Dind
  real, allocatable,dimension(:,:) :: rmse1D,rmse1D2
  real, dimension(Ncases) :: rmseavg
  real, dimension(8,Ncases) :: rmseavgind
  real,dimension(8) :: refval,tmpvec

  integer :: output_ncid,dimids(Nparam), param1_varid,param2_varid,param3_varid,param4_varid,param5_varid,param6_varid
  character(50),dimension(Nparam):: dimnames
  integer, dimension(6) :: var_out_id



  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,ierr)

  if (Nparam.ne.6) then
     write(*,*) 'ERROR: Nparam has to be 4. Stopping.'
     call MPI_FINALIZE(ierr)
  end if

  allocate(mpi_id(Nproc-1))
  allocate(Nlocalvec(Nproc))
  allocate(Displace(Nproc))
  allocate(Recv_Request(Nproc-1))
  allocate(mpi_status(MPI_STATUS_SIZE,Nproc-1))
  allocate(mpi_status1(MPI_STATUS_SIZE))

  itmp(1) = 0
  do i = 0, Nproc-1
    if(i /= my_id) then
      itmp(1) = itmp(1) + 1
      mpi_id( itmp(1) ) = i
    endif
  enddo

  call getarg(1,param)
  path_in = trim(adjustl(param))

  Ncombo=1
  do i=1,Nparam

     write(i_str,'(i2)') i-1
     filein = trim(path_in)//'/param_'//trim(adjustl(i_str))
     open(77,file=trim(filein),status='old',form='formatted')
     k=0
     do while (.true.)
       if (k.eq.0) then 
          read (77,*,end=10) dimnames(i)
       else
          read (77,*,end=10) readparam(1)
       end if
       k=k+1
     end do
10   continue
     close(77)
     Nrange(i) = k-1
    
     Ncombo=Ncombo*Nrange(i)
  end do

  allocate(pvalue(Nparam,MAXVAL(Nrange)))

  do i=1,Nparam

     write(i_str,'(i2)') i-1
     filein = trim(path_in)//'/param_'//trim(adjustl(i_str))
     open(77,file=trim(filein),status='old',form='formatted')
     read (77,*) dimnames(i)
     k=1
     do while (k.le.Nrange(i))
       read (77,*) pvalue(i,k)
       k=k+1
     end do
     close(77)
  end do

 
  itmp(1) = Ncombo/Nproc

  Displace(1)=0
  do i=0,Nproc-1
  run_start = i*itmp(1)
  run_end   = run_start+itmp(1)-1
  if(i == Nproc-1) run_end = Ncombo-1
  Nlocalvec(i+1) = run_end - run_start+1
  if (i.ge.1) Displace(i+1) = Displace(i) + Nlocalvec(i) 
  end do


  run_start = my_id*itmp(1)
  run_end   = run_start+itmp(1)-1
  if(my_id == Nproc-1) run_end = Ncombo-1
  N_local_runs = run_end - run_start+1

  allocate(rmse1Dloc(Ncases+1,N_local_runs),rmse1Dloc2(Ncases+1,N_local_runs),rmse1Dlocind(8,Ncases,N_local_runs))
  rmse1Dloc=0.
  rmse1Dloc2=0.
  rmse1Dlocind=0.0

  ncount=1
  do n=run_start,run_end

    i = n/(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) 
    j = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))
    k = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(4)*Nrange(5)*Nrange(6))
    l = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(5)*Nrange(6))
    m = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6))-l*(Nrange(5)*Nrange(6)))/Nrange(6)
    o = n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6))-l*(Nrange(5)*Nrange(6))-m*Nrange(6)
    write(i_str,'(i2.2)') i
    write(j_str,'(i2.2)') j
    write(k_str,'(i2.2)') k
    write(l_str,'(i2.2)') l
    write(m_str,'(i2.2)') m
    write(o_str,'(i2.2)') o

    do ca=1,Ncases
      runname = trim(casename(ca))//'_'//trim(i_str)//'_'//trim(j_str)//'_'//trim(k_str)//'_'//trim(l_str)//'_'//trim(m_str)//'_'//trim(o_str)

      ! read rmsecase
      open(78,file=trim(path_in)//"/"//trim(runname)//"/"//'rmse_thvd_'//trim(casename(ca)),status='old',form='formatted')
      read (78,*) readparam(1:9)
      close(78)
      open(79,file=trim(path_in)//"/"//trim(runname)//"/"//'rmse_thvd_'//trim(casename(ca)),status='old',form='formatted')
      read (79,*) readparam2(1:9)
      close(79)
 
      tmpvec=readparam2(2:9)
      if (trim(casename(ca)).eq."BOMEX") then
        refval = (/300.,0.0,14.e-3,3.e-6,8.,150.,0.3,0.3/)
      else
        refval = (/300.,9.26e-6,5.e-3,0.0,20.,70.,0.3,0.3/)
      end if
      do m=1,8
        rmse1Dlocind(m,ca,ncount)=tmpvec(m) 
      end do

      rmse1Dloc(ca+1,ncount) = readparam(1)
      ! get rmse total of run 
      if (trim(casename(ca)).eq."BOMEX") then
        rmse1Dloc(1,ncount) = rmse1Dloc(1,ncount) + 0.5 * readparam(1)
      else
        rmse1Dloc(1,ncount) = rmse1Dloc(1,ncount) + 0.5 * readparam(1)
      end if

      rmse1Dloc2(ca+1,ncount) = readparam2(1)
      ! get rmse total of run 
      if (trim(casename(ca)).eq."BOMEX") then
        rmse1Dloc2(1,ncount) = rmse1Dloc2(1,ncount) + 0.5 * readparam2(1)
      else
        rmse1Dloc2(1,ncount) = rmse1Dloc2(1,ncount) + 0.5 * readparam2(1)
      end if
      
    
    end do

    ncount=ncount+1
  end do

  
  if (my_id.eq.0) then

  allocate(rmse1D(Ncases+1,Ncombo),rmse1D2(Ncases+1,Ncombo),rmse1Dind(8,Ncases,Ncombo))
  rmse1D=-999.
  rmse1D2=-999.
  rmse1Dind=-999.
 
  end if

  do i=1,Ncases+1
    call MPI_GATHERV(rmse1Dloc(i,:),N_local_runs,MPI_REAL8, &
                            rmse1D(i,:),Nlocalvec,Displace,MPI_REAL8, &
                            0,MPI_COMM_WORLD,ierr)
    call MPI_GATHERV(rmse1Dloc2(i,:),N_local_runs,MPI_REAL8, &
                            rmse1D2(i,:),Nlocalvec,Displace,MPI_REAL8, &
                            0,MPI_COMM_WORLD,ierr)
    if (i.le.Ncases) then
    do m=1,8
      call MPI_GATHERV(rmse1Dlocind(m,i,:),N_local_runs,MPI_REAL8, &
                            rmse1Dind(m,i,:),Nlocalvec,Displace,MPI_REAL8, &
                            0,MPI_COMM_WORLD,ierr)
    end do
    end if
  end do

  if (my_id.eq.0) then
 

  ! write netcdf file 
  fname_out =trim(path_in)//"/"//'RMSE.nc' 
  write(*,*) fname_out
  call check_nc( nf90_create(fname_out,nf90_clobber,output_ncid) )

  !  Define Dimensions
  call check_nc( nf90_def_dim(output_ncid,trim(dimnames(1)),Nrange(1),dimids(1)) )
  call check_nc( nf90_def_dim(output_ncid,trim(dimnames(2)),Nrange(2),dimids(2)) )
  call check_nc( nf90_def_dim(output_ncid,trim(dimnames(3)),Nrange(3),dimids(3)) )
  call check_nc( nf90_def_dim(output_ncid,trim(dimnames(4)),Nrange(4),dimids(4)) )
  call check_nc( nf90_def_dim(output_ncid,trim(dimnames(5)),Nrange(5),dimids(5)) )
  call check_nc( nf90_def_dim(output_ncid,trim(dimnames(6)),Nrange(6),dimids(6)) )

  !  Define Variables
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(1)),nf90_float,dimids(1),param1_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(2)),nf90_float,dimids(2),param2_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(3)),nf90_float,dimids(3),param3_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(4)),nf90_float,dimids(4),param4_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(5)),nf90_float,dimids(5),param5_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(6)),nf90_float,dimids(6),param6_varid) )


  !call check_nc( nf90_def_var( output_ncid,"RMSE_tot",nf90_float,dimids,var_out_id(1) ) )
  !call check_nc( nf90_def_var( output_ncid,"RMSE_"//trim(casename(1)),nf90_float,dimids,var_out_id(2) ) )
  !call check_nc( nf90_def_var( output_ncid,"RMSE_"//trim(casename(2)),nf90_float,dimids,var_out_id(3) ) )
  call check_nc( nf90_def_var( output_ncid,"RMSE_thvd_tot",nf90_float,dimids,var_out_id(4) ) )
  call check_nc( nf90_def_var( output_ncid,"RMSE_thvd_"//trim(casename(1)),nf90_float,dimids,var_out_id(5) ) )
  call check_nc( nf90_def_var( output_ncid,"RMSE_thvd_"//trim(casename(2)),nf90_float,dimids,var_out_id(6) ) )
  !call check_nc( nf90_def_var( output_ncid,"RMSE_"//trim(casename(3)),nf90_float,dimids,var_out_id(4) ) )

  call check_nc( nf90_enddef(output_ncid) )

  !  Put variables
  call check_nc( nf90_put_var(output_ncid, param1_varid, pvalue(1,1:Nrange(1)) ) )
  call check_nc( nf90_put_var(output_ncid, param2_varid, pvalue(2,1:Nrange(2)) ) )
  call check_nc( nf90_put_var(output_ncid, param3_varid, pvalue(3,1:Nrange(3)) ) )
  call check_nc( nf90_put_var(output_ncid, param4_varid, pvalue(4,1:Nrange(4)) ) )
  call check_nc( nf90_put_var(output_ncid, param5_varid, pvalue(5,1:Nrange(5)) ) )
  call check_nc( nf90_put_var(output_ncid, param6_varid, pvalue(6,1:Nrange(6)) ) )

  !fill parameter array
  allocate(rmse(Nrange(1),Nrange(2),Nrange(3),Nrange(4),Nrange(5),Nrange(6)))
  rmseavgind=0.0


  rmseavgind(1,1) = 0.298500898501127
  rmseavgind(2,1) = -999.
  rmseavgind(3,1) = 4.709051607813674e-4
  rmseavgind(4,1) = 2.291507888805414e-6
  rmseavgind(5,1) = 8.29152631298648
  rmseavgind(6,1) = 37.5370361230699
  rmseavgind(7,1) = 0.110751867145755
  rmseavgind(8,1) = 4.754294124718259e-4


  rmseavgind(1,2) = 0.142648591929376
  rmseavgind(2,2) = 9.579204282494366e-5
  rmseavgind(3,2) = 1.561694507888805e-4
  rmseavgind(4,2) = -999.
  rmseavgind(5,2) = 5.55608093415176
  rmseavgind(6,2) = 18.0438175345439
  rmseavgind(7,2) = 0.119297832160030
  rmseavgind(8,2) = 1.220929376408715e-4


  ! compute average error for each case and each parameter 
  write(*,*) 'Parameter-averaged rmse'
  do ca=1,Ncases
    do m=1,8
      !rmseavgind(m,ca) = sum(rmse1Dind(m,ca,:))/float(Ncombo) 
      !write(*,'(A,A,x,A,i1.1,x,f10.7)') 'case ',casename(ca),', parameter ', m, rmseavgind(m,ca)
      write(*,*) 'case ',casename(ca),', parameter ', m ,' ',rmseavgind(m,ca)
    end do
  end do

  ! compute error for each case by normalizing with the parameter averages
  rmse1D=0.0
  do ca=1,Ncases
    do m=1,6
      if ((trim(casename(ca)).eq."BOMEX".and.m.ne.2).or. &
      (trim(casename(ca)).ne."BOMEX".and.m.ne.4)) then
        write(*,*) trim(casename(ca)),' ',m
        rmse1D(1+ca,:) = rmse1D(1+ca,:) + 1./7. * rmse1Dind(m,ca,:) / rmseavgind(m,ca)
      end if
    end do
  end do

  ! compute mean case-specific error
  ! NOT NEEDED SINCE 1 PER DEFINITION
  !write(*,*) 'Case-averaged rmse (thvd NOT included)'
  !do ca=1,Ncases
  !  rmseavg(ca) = sum(rmse1D(1+ca,:))/float(Ncombo)
  !  write(*,*) 'case ',casename(ca),' ', rmseavg(ca)
  !end do
  ! compute total error from case-specific errors
  ! normalize each case by the mean case error (not needed since equal to 1)

  rmse1D(1,:) = 0.0
  do ca=1,Ncases
    rmse1D(1,:) = rmse1D(1,:) + 1./float(Ncases) * rmse1D(1+ca,:) 
  end do

  ! rearrange
  do ca=1,Ncases+1
  do n=0,Ncombo-1
    i = n/(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))
    j = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))
    k = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(4)*Nrange(5)*Nrange(6))
    l = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(5)*Nrange(6))
    m = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6))-l*(Nrange(5)*Nrange(6)))/Nrange(6)
    o = n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6))-l*(Nrange(5)*Nrange(6))-m*Nrange(6)    
    rmse(i+1,j+1,k+1,l+1,m+1,o+1) = rmse1D(ca,n+1)

    !if (ca.eq.1.and.pvalue(1,i+1).eq.1400. .and.pvalue(2,j+1).eq.0.8 .and.pvalue(3,k+1).eq.400. .and.pvalue(4,l+1).eq.0.4) then
    !  write(*,*) 'tot: ', rmse1D(1,n+1), ' BOMEX: ', rmse1D(2,n+1),' CPBL: ',rmse1D(3,n+1)
    !  do m=1,7
    !    write(*,*) m,' BOMEX: ', rmse1Dind(m,1,n+1),' CPBL: ',rmse1Dind(m,2,n+1) 
    !  end do
    !end if

  end do
 ! call check_nc( nf90_put_var(output_ncid, var_out_id(ca), rmse ) ) 
  end do
  deallocate(rmse1D)

  !! redo this for error including thvd
 
  ! compute error for each case by normalizing with the parameter averages
  !rmse1D2=0.0
  !do ca=1,Ncases
  !  do m=1,8
  !    if ((trim(casename(ca)).eq."BOMEX".and.m.ne.2).or. &
  !    (trim(casename(ca)).ne."BOMEX".and.m.ne.4)) then
  !      rmse1D2(1+ca,:) = rmse1D2(1+ca,:) + 1./7. * rmse1Dind(m,ca,:) / rmseavgind(m,ca)
  !    end if
  !  end do
!  end do

  ! compute mean case-specific error
  !write(*,*) 'Case-averaged rmse (thvd included)'
  !do ca=1,Ncases
  !  rmseavg(ca) = sum(rmse1D2(1+ca,:))/float(Ncombo)
  !  write(*,'(A,A,x,f10.7)') 'case ',casename(ca), rmseavg(ca)
  !end do
  ! compute total error from case-specific errors
  ! normalize each case by the mean case error (not needed since average is 1)
  !rmse1D2(1,:) = 0.0
  !do ca=1,Ncases
  !  rmse1D2(1,:) = rmse1D2(1,:) + 1./float(Ncases) * rmse1D2(1+ca,:) 
  !end do

  do ca=1,Ncases+1
  do n=0,Ncombo-1
    i = n/(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))
    j = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))
    k = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(4)*Nrange(5)*Nrange(6))
    l = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6)))/(Nrange(5)*Nrange(6))
    m = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6))-l*(Nrange(5)*Nrange(6)))/Nrange(6)
    o = n-i*(Nrange(2)*Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6))-j*(Nrange(3)*Nrange(4)*Nrange(5)*Nrange(6)) - k * (Nrange(4)*Nrange(5)*Nrange(6))-l*(Nrange(5)*Nrange(6))-m*Nrange(6)
    rmse(i+1,j+1,k+1,l+1,m+1,o+1) = rmse1D2(ca,n+1)
  end do
  call check_nc( nf90_put_var(output_ncid, var_out_id(3+ca), rmse ) ) 
  end do
  deallocate(rmse1D2)

  call check_nc( nf90_close(output_ncid) )


   
  end if
 
  call MPI_Finalize(ierr)


  if(my_id == 0) write(*,*) 'MPI finalized ...'

  
contains

  subroutine check_nc(status)
    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      call exit()
    endif
  end subroutine check_nc

  function char2int(c)
    implicit none

    character(*), intent(in) :: c

    integer :: char2int
    integer :: i
    character(100) :: c2

    c2 = trim(c)   ! Eliminate prefixed blanks
    char2int = 0
    do i = 1,len(trim(c2))
      char2int = 10*char2int + index('0123456789',c2(i:i)) - 1
    end do

  end function char2int


end Program Parallel_Get_Rmse

