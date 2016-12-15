Program Parallel_Get_Rmse
  use mpi
  use netcdf
  USE IFPORT
  implicit none

  character(300) :: param, path_in,pathin,path_ncl
  character(50),dimension(:),allocatable :: runname
  character(200) :: namehlp,command, casename,filein

  integer :: my_id, Nproc, ierr, Nsetup, N_local_runs, N_max
  integer, allocatable, dimension(:) :: mpi_id

  integer :: i,k
  integer :: itmp(10), Ncombo
  integer :: run_start, run_end, status
  character(2) :: i_str
  real :: readparam

  integer,parameter :: Nparam=4
  integer,dimension(Nparam) :: Nrange
  real, allocatable,dimension(:,:,:,:) :: pvalue


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

  Ncombo=1
  do i=1,Nparam

     write(i_str,'(i2)') i-1
     open(77,file=trim(path_in)//'/param_'//trim(adjustl(i_str)),status='old',form='formatted')
     k=0
     do while (.true.)
       read (77,*,end=10) readparam
       k=k+1
     end do
10   continue
     close(77)
     Nrange(i) = k
    
     Ncombo=Ncombo*Nrange(i)
  end do
 
  write(*,*) Ncombo
  write(*,*) Nrange(:)


  allocate(pvalue(Nrange(1),Nrange(2),Nrange(3),Nrange(4)))

 
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


end Program Parallel_Get_Rmse

