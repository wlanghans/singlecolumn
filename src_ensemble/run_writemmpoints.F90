Program Parallel_Run_Ensemble
  use mpi
  use netcdf
  USE IFPORT
  implicit none

  character(300) :: param, path_in,pathin,path_ncl
  character(50),dimension(:),allocatable :: runname
  character(500) :: namehlp, casename,filein,fileinprev
  character(700) ::command, fileprev, filenow

  integer :: my_id, Nproc, ierr, Nruns, N_local_runs, N_max
  integer, allocatable, dimension(:) :: mpi_id

  integer :: i,j,k
  integer :: itmp(10)
  integer :: run_start, run_end, status

  logical :: fileexist, filenowexist
 
  real(8) :: t1,t2


  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,ierr)

  allocate(mpi_id(Nproc-1))

  itmp(1) = 0
  do i = 0, Nproc-1
    if(i /= my_id) then
      itmp(1) = itmp(1) + 1
      mpi_id( itmp(1) ) = i
    endif
  enddo

  call getarg(1,param)
  path_in = trim(adjustl(param))

  call getarg(2,param)
  Nruns = char2int( trim(adjustl(param)) )
  
  call getarg(3,param)
  path_ncl='"'//trim(adjustl(param))//'"'


  itmp(1) = Nruns/Nproc

  run_start = my_id*itmp(1)
  run_end   = run_start+itmp(1)-1
  if(my_id == Nproc-1) run_end = Nruns-1
  N_local_runs = run_end - run_start+1
  write(*,*)  'ID ',my_id, ' has ',N_local_runs,' runs'

  call MPI_Allreduce(N_local_runs,N_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

  allocate(runname(N_max))

  open(77,file=trim(path_in)//'/runlist',status='old',form='formatted')

  k=0
  do i=0,Nruns-1
      read(77,*) namehlp
      if (run_start.le.i.and.run_end.ge.i) then 
        k=k+1
        runname(k) = namehlp
      end if
  end do        

  close(77)

  if (k.ne.N_local_runs) then
      write(*,*) 'Something is wrong!',k,N_local_runs
      write(*,*) 'ID ',my_id, 'IDs:', run_start, ' IDe:',run_end 
      call MPI_FINALIZE(ierr)
  end if


  do i=1,N_local_runs
   namehlp = trim(runname(i))
   if (namehlp(1:5).eq."BOMEX"  ) then
      casename = '"BOMEX"'
      pathin='"'//trim(path_in)//"/"//trim(runname(i))//"/"//'"'
      filein='"BOMEX_snapshot.nc"'
   elseif (namehlp(1:5).eq."CPBL2"  ) then
      casename = '"CPBL2"'
      pathin='"'//trim(path_in)//"/"//trim(runname(i))//"/"//'"'
      filein='"WITEK11_snapshot.nc"'
   elseif (namehlp(1:5).eq."CPBL4"  ) then
      casename = '"CPBL4"'
      pathin='"'//trim(path_in)//"/"//trim(runname(i))//"/"//'"'
      filein='"WITEK11_snapshot.nc"'
   else
      write(*,*) 'ERROR: Case unknown. Stopping.'
      call MPI_FINALIZE(ierr)
   end if


   ! make sure previous run has completed yet
   if (i.gt.1) then
   INQUIRE(FILE=trim(fileprev),EXIST=fileexist)
   
   if (.not.fileexist) then

      t1  = MPI_WTIME() 
 
      do while (.not.fileexist)
         INQUIRE(FILE=trim(fileprev),EXIST=fileexist)
         t2 = MPI_WTIME() 
         if (t2-t1.gt.240.) then
            write(*,*) 'ERROR: Waited for 4 min but previous simulation still incomplete.Stopping.'
            call MPI_FINALIZE(ierr)
         end if 
      end do

   end if

   ! remove netcdf file from previous run
   write(*,*) 'Simulation ', trim(runname(i-1)), ' was successful'
   command=" rm -rf "//trim(path_in)//"/"//trim(runname(i-1))//"/"//trim(fileinprev)
   !status = system( trim(command))

   end if

   command="cd "//trim(path_in)//"/"//trim(runname(i))//"; ncl -Q "//trim(path_ncl)//"/writesaturationprofiles.ncl casename="//"'"//trim(casename)//"'"//" filein="//"'"//trim(filein)//"'"//" pathin="//"'"//trim(pathin)//"'"
   !command=trim(command)//"; rm -rf "//trim(path_in)//"/"//trim(runname(i))//"/"//trim(filein)

   filenow=trim(path_in)//"/"//trim(runname(i))//"/mminput.nc"
   !INQUIRE(FILE=trim(filenow),EXIST=filenowexist)
   !if (.not.filenowexist) then

   status = system( trim(command))
   !INQUIRE(FILE=trim(filenow),EXIST=filenowexist)
   !if (.not.filenowexist) then

   !   t1  = MPI_WTIME() 
 
   !   do while (.not.filenowexist)
   !      INQUIRE(FILE=trim(filenow),EXIST=filenowexist)
   !      t2 = MPI_WTIME() 
   !      if (t2-t1.gt.240.) then
   !         write(*,*) 'ERROR: Waited for 4 min to get postprocessing done but something is wrong. Stopping.'
   !         call MPI_FINALIZE(ierr)
   !      end if  
   !      command=" cd "//trim(path_in)//"/"//trim(runname(i))//"; ncl -Q "//trim(path_ncl)//"/getrmse.ncl benchdir="//"'"//trim(path_ncl)//"'"//" casename="//"'"//trim(casename)//"'"//" filein="//"'"//trim(filein)//"'"//" pathin="//"'"//trim(pathin)//"'"
   !      status = system( trim(command))
   !   end do

   !end if
 

   !else
 
   !  write(*,*) 'Simulation ', trim(runname(i)), ' was already completed earlier'

   !end if

!    if ( status .ne. 0 ) then 
!       write(*,*) 'ERROR: Simulation ', i,' /',N_local_runs,' ', trim(runname(i)), ' was unsuccessful'
!       call MPI_FINALIZE(ierr)
!    else
!       write(*,*) 'Simulation ', trim(runname(i)), ' was successful'
!    end if
!   command = "cd "//trim(path_in)//"/"//trim(runname(i))//"; ncl -Q "//trim(path_ncl)//"/getrmse.ncl casename="//"'"//trim(casename)//"'"//" filein="//"'"//trim(filein)//"'"//" pathin="//"'"//trim(pathin)//"'"
   !status = system(command)
   ! if ( status .ne. 0 ) then 
   !    write(*,*) 'ERROR: ncl postproc ', trim(runname(i)), ' was unsuccessful'
   !    write(*,*) command
   !    call MPI_FINALIZE(ierr)
   ! end if


   fileinprev=filein
   fileprev=filenow
  end do

  call MPI_Finalize(ierr)

  if(my_id == 0) write(*,*) 'MPI finalized ...'

  
contains

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


end Program Parallel_Run_Ensemble

