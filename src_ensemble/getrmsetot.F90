Program Parallel_Get_Rmse
  use mpi
  use netcdf
  USE IFPORT
  implicit none

  character(300) :: param, path_in,pathin,path_ncl
  character(50) :: runname
  character(200) :: namehlp,command, filein,fname_out

  integer :: my_id, Nproc, ierr, Nsetup, N_local_runs, N_max, ncount
  integer, allocatable, dimension(:) :: mpi_id, Nlocalvec,displace

  integer :: i,j,k,l,n,ca
  integer :: itmp(10), Ncombo
  integer :: run_start, run_end, status
  character(2) :: i_str,j_str,k_str,l_str
  real :: readparam(6)
  integer, dimension(:), allocatable :: Recv_Request
  integer, dimension(:,:), allocatable :: mpi_status
  integer, dimension(:), allocatable :: mpi_status1
  integer ::  Send_Request


  integer,parameter :: Nparam=4
  integer,parameter :: Ncases=3
  character(10),dimension(Ncases) :: casename=(/"BOMEX","CPBL2","CPBL4"/)
  integer,dimension(Nparam) :: Nrange
  real, allocatable,dimension(:,:,:,:) :: rmse
  real, allocatable,dimension(:,:) :: pvalue
  real, allocatable,dimension(:,:) :: rmse1Dloc
  real, allocatable,dimension(:,:) :: rmse1D

  integer :: output_ncid,dimids(Nparam), param1_varid,param2_varid,param3_varid,param4_varid
  character(50),dimension(Nparam):: dimnames
  integer, dimension(4) :: var_out_id


  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,ierr)

  if (Ncases.ne.3.or.Nparam.ne.4) then
     write(*,*) 'ERROR: Ncases has to be 3 and Nparam has to be 4. Stopping.'
     call MPI_FINALIZE(ierr)
  end if

  allocate(mpi_id(Nproc-1))
  allocate(Nlocalvec(Nproc))
  allocate(Displace(Nproc))
  allocate(Recv_Request(Nproc-1))
  allocate( mpi_status(MPI_STATUS_SIZE,Nproc-1) )
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

  allocate(rmse1Dloc(Ncases+1,N_local_runs))
  rmse1Dloc=0.

  ncount=1
  do n=run_start,run_end

    i = n/(Nrange(2)*Nrange(3)*Nrange(4)) 
    j = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)))/(Nrange(3)*Nrange(4))
    k = (n-i*(Nrange(2)*Nrange(3)*Nrange(4))-j*(Nrange(3)*Nrange(4)))/(Nrange(4))
    l = n-i*(Nrange(2)*Nrange(3)*Nrange(4))-j*(Nrange(3)*Nrange(4)) - k * Nrange(4)
    write(i_str,'(i2.2)') i
    write(j_str,'(i2.2)') j
    write(k_str,'(i2.2)') k
    write(l_str,'(i2.2)') l

    do ca=1,Ncases
      runname = trim(casename(ca))//'_'//trim(i_str)//'_'//trim(j_str)//'_'//trim(k_str)//'_'//trim(l_str)

      ! read rmsecase
      open(78,file=trim(path_in)//"/"//trim(runname)//"/"//'rmse_'//trim(casename(ca)),status='old',form='formatted')
      read (78,*) readparam(1:6)
      close(78)

      rmse1Dloc(ca+1,ncount) = readparam(1)
      ! get rmse total of run 
      if (trim(casename(ca)).eq."BOMEX") then
        rmse1Dloc(1,ncount) = rmse1Dloc(1,ncount) + 0.5 * readparam(1)
      else
        rmse1Dloc(1,ncount) = rmse1Dloc(1,ncount) + 0.25 * readparam(1)
      end if
    
    end do

    ncount=ncount+1
  end do

  
  if (my_id.eq.0) then

  allocate(rmse1D(Ncases+1,Ncombo))
  rmse1D=-999.
 
  end if

  do i=1,Ncases+1
    call MPI_GATHERV(rmse1Dloc(i,:),N_local_runs,MPI_REAL8, &
                            rmse1D(i,:),Nlocalvec,Displace,MPI_REAL8, &
                            0,MPI_COMM_WORLD,ierr)
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

  !  Define Variables
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(1)),nf90_float,dimids(1),param1_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(2)),nf90_float,dimids(2),param2_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(3)),nf90_float,dimids(3),param3_varid) )
  call check_nc( nf90_def_var(output_ncid,trim(dimnames(4)),nf90_float,dimids(4),param4_varid) )


  call check_nc( nf90_def_var( output_ncid,"RMSE_tot",nf90_float,dimids,var_out_id(1) ) )
  call check_nc( nf90_def_var( output_ncid,"RMSE_"//trim(casename(1)),nf90_float,dimids,var_out_id(2) ) )
  call check_nc( nf90_def_var( output_ncid,"RMSE_"//trim(casename(2)),nf90_float,dimids,var_out_id(3) ) )
  call check_nc( nf90_def_var( output_ncid,"RMSE_"//trim(casename(3)),nf90_float,dimids,var_out_id(4) ) )

  call check_nc( nf90_enddef(output_ncid) )

  !  Put variables
  call check_nc( nf90_put_var(output_ncid, param1_varid, pvalue(1,1:Nrange(1)) ) )
  call check_nc( nf90_put_var(output_ncid, param2_varid, pvalue(2,1:Nrange(2)) ) )
  call check_nc( nf90_put_var(output_ncid, param3_varid, pvalue(3,1:Nrange(3)) ) )
  call check_nc( nf90_put_var(output_ncid, param4_varid, pvalue(4,1:Nrange(4)) ) )

  !fill parameter array
  allocate(rmse(Nrange(1),Nrange(2),Nrange(3),Nrange(4)))
  do ca=1,Ncases+1
  do n=0,Ncombo-1
    i = n/(Nrange(2)*Nrange(3)*Nrange(4))
    j = (n-i*(Nrange(2)*Nrange(3)*Nrange(4)))/(Nrange(3)*Nrange(4))
    k = (n-i*(Nrange(2)*Nrange(3)*Nrange(4))-j*(Nrange(3)*Nrange(4)))/(Nrange(4))
    l = n-i*(Nrange(2)*Nrange(3)*Nrange(4))-j*(Nrange(3)*Nrange(4)) - k *Nrange(4)
    rmse(i+1,j+1,k+1,l+1) = rmse1D(ca,n+1)
  end do
  call check_nc( nf90_put_var(output_ncid, var_out_id(ca), rmse ) ) 
  end do
  deallocate(rmse1D)

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

