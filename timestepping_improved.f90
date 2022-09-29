subroutine timestepping_improved(n,Blu,Blv,Bru,Brv,hu,hv,tau,T,F,k,u,v)

  implicit none

  integer, intent(in) :: n
  real(kind=8), intent(in), dimension(2,n) :: Blu,Blv,Bru,Brv
  real(kind=8), intent(in), dimension(n) :: hu,hv
  real(kind=8), intent(in) :: tau, T, F, k
  real(kind=8), intent(inout), dimension(n) :: u,v

  !internal variables

  integer :: i

  !this is very similar to timestepping but with the new banded variables as the inputs. We repeat alg.1 with banded_crank here, with the new values of u and v overwriting the older ones. The heun call here is for the same function as before in timestepping.f90 (but this time without the annotations.
  do i = 0, (ceiling(T/tau)-1)
     !perform line 4 of alg.1 using the inproved and more efficient banded version of crank_nicholson
     call banded_crank(n,tau,Blu,Blv,Bru,Brv,hu,hv,u,v)
     !perform line 5 of alg.1 using the same heun as in timestepping.f90
     call heun(n,tau,F,k,u,v)
     !perform line 6 of alg.1 using the improved banded version of crank_nicholson
     call banded_crank(n,tau,Blu,Blv,Bru,Brv,hu,hv,u,v)
  end do

end subroutine timestepping_improved
