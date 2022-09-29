program grayscott

  implicit none

  integer :: n
  real(kind=8) :: tau, T, F, k
  real(kind=8), allocatable, dimension(:,:) :: Blu, Blv, Bru, Brv
  real(kind=8), allocatable, dimension(:) :: u,v,hu,hv

  real(kind=8) :: t_start, t_finish


  !read the input.dat file and import the correct variables
  open(unit=2,file="input.dat")

  read(2,*) n
  read(2,*) tau
  read(2,*) T
  read(2,*) F
  read(2,*) k

  close(2)

  write(*, '(A,X,I5)')   'n  = ', n
  write(*, '(A,X,F7.3)') 'tau= ', tau
  write(*, '(A,X,F8.1)') 'T  = ', T
  write(*, '(A,X,F9.5)') 'F  = ', F
  write(*, '(A,X,F9.5)') 'k  = ', k


  !allocate the appropriate matrices and populate them with the create_banded call
  allocate(Blu(2,n),Blv(2,n),Bru(2,n),Brv(2,n),hu(n),hv(n),u(n),v(n))

  call create_banded(n,tau,Blu,Blv,Bru,Brv,hu,hv)

  !call test_banded here to show the error between the bended_crank method and a random input vector
  call test_banded(n,tau,Blu,Blv,Bru,Brv,hu,hv)
  
  !call initial conditions for u and v (same as grayscott for non-banded case).
  call initial(n,u,v)

  call cpu_time(t_start)

  !call the new improved timestepping method and time the time it takes to run.
  call timestepping_improved(n,Blu,Blv,Bru,Brv,hu,hv,tau,T,F,k,u,v)

  call cpu_time(t_finish)

  !print out the middle value of u and the cpu time
  write(*, '(A,X,F15.12)') 'u_{n/2}  = ', u(n/2)
  write(*, '(A,X,F7.2)')  'cpu time = ', t_finish - t_start


  !save the solution to the file solution.dat and then deallocate the necessary matrices.
  call save_fields(n,u,v,'solution.dat')

  deallocate(Blu,Blv,Bru,Brv,hu,hv,u,v)


end program grayscott
