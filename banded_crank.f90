subroutine banded_crank(n,tau,Blu,Blv,Bru,Brv,hu,hv,u,v)

  implicit none

  integer, intent(in) :: n
  real(kind=8), intent(in) :: tau
  real(kind=8), intent(in), dimension(2,n) :: Blu, Blv, Bru, Brv
  real(kind=8), intent(in), dimension(n) :: hu, hv
  real(kind=8), intent(inout), dimension(n) :: u, v
  
  
  !internal variables
  
  integer :: info
  real(kind=8), external :: ddot
  real(kind=8), dimension(n) :: c, r, copy1
  real(kind=8), dimension(n-1) :: copy2
  real(kind=8) :: tauu, tauv, dot

  real(kind=8), parameter :: Du = 2d-5
  real(kind=8), parameter :: Dv = 1d-5

  !calculate tau_u and tau_v
  tauu = Du*tau*0.25_8*n**2
  tauv = Dv*tau*0.25_8*n**2

  !populate the vector c
  c = 0.0_8
  c(1) = 1.0_8
  c(n) = 1.0_8

  !calculate ru

  !calculate the matrix vector product B_ru*u and output the result in r
  call dsbmv('U',n,1,1.0_8,Bru,2,u,1,0.0_8,r,1)
  !calculate the scalar for the c vector in (14) and then add to the r
  dot = tauu*ddot(n,c,1,u,1)
  call daxpy(n,dot,c,1,r,1) 
  
  !copy the necessary parts of the matrix Blu to be able to solve the equation B_lu*uhat=r_u using r = r_u and outputing uhat as r
  call dcopy((n-1),Blu(1,2),2,copy2,1)
  call dcopy(n,Blu(2,1),2,copy1,1)
  call dptsv(n,1,copy1,copy2,r,n,info)

  !calculate u using the given equation u = uhat + h_u*(c^T*uhat), calculate c^T*uhat first then use daxpy to add the vectors together and then copy the result back to the output u
  dot = ddot(n,c,1,r,1)
  call daxpy(n,dot,hu,1,r,1)
  call dcopy(n,r,1,u,1)

  !calculate rv

  !calculate the matrix vector product B_rv*v and the output the result in r
  call dsbmv('U',n,1,1.0_8,Brv,2,v,1,0.0_8,r,1)
  !calculte the scalar for the c vector in (14) and then add to the r
  dot = tauv*ddot(n,c,1,v,1)
  call daxpy(n,dot,c,1,r,1)

  !copy the necessary parts of the matrix Blv to be able to solve the equation B_lv*vhat=r_v using r = r_u and outputting vhat as r
  call dcopy((n-1),Blv(1,2),2,copy2,1)
  call dcopy(n,Blv(2,1),2,copy1,1)
  call dptsv(n,1,copy1,copy2,r,n,info)

  !calculate v using the given equation v = vhat + h_v(c^T*vhat), calculate c^T*vhat first, then use daxpy to add the vectors together and then copy the result back to the output v
  dot = ddot(n,c,1,r,1)
  call daxpy(n,dot,hv,1,r,1)
  call dcopy(n,r,1,v,1)

end subroutine banded_crank
