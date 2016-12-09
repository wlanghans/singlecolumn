Program Parallel_Run_Ensemble
  use mpi
  use netcdf
  USE IFPORT
  implicit none

  character(100) :: param, path_in 
  character(50),dimension(:),allocatable :: runname
  character(200) :: namehlp,command

  integer :: my_id, Nproc, ierr, Nruns, N_local_runs, N_max
   integer, allocatable, dimension(:) :: mpi_id

  integer :: i,k
  integer :: itmp(10)
  integer :: run_start, run_end, status


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
   command = "cd "//trim(path_in)//"; ./sc_sam_wl "//trim(runname(i))//" > "//trim(path_in)//"/"//trim(runname(i))//"/"//trim(runname(i))//".out"
   status = system( command)
    if ( status .ne. 0 ) then 
       write(*,*) 'ERROR: Simulation ', trim(runname(i)), ' was unsuccessful'
    end if
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
