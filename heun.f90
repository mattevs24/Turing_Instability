subroutine heun(n,tau,F,k,u,v)
  ! ====================================
  ! Implement your Heun method step here
  ! ====================================

  integer, intent(in) :: n
  real(kind=8), intent(in) :: tau, F, k
  real(kind=8), intent(inout), dimension(n) :: u,v

  !variables within the subroutine

  real(kind=8), dimension(n) :: Nu, Nv, ucopy, vcopy
  integer :: i

  !--------------------------------------------------------------------------------!

  !This is the same program as in partI and is commented properly in that directory!

  !--------------------------------------------------------------------------------!

  !calculate Nu and Nv of W0
  
  do i = 1,n
     Nu(i) = -u(i)*v(i)**2 + F - F*u(i)
     Nv(i) = u(i)*v(i)**2 - (F+k)*v(i)
  end do

  call dcopy(n,u,1,ucopy,1)
  call dcopy(n,v,1,vcopy,1)

  !(12.1, W1/2)
  call daxpy(n,tau,Nu,1,ucopy,1)
  call daxpy(n,tau,Nv,1,vcopy,1)

  !first part of the calc for W1
  call daxpy(n,tau/2.0_8,Nu,1,u,1)
  call daxpy(n,tau/2.0_8,Nv,1,v,1)

  !N(W(1/2))
  do i = 1,n
     Nu(i) = -ucopy(i)*vcopy(i)**2 + F - F*ucopy(i)
     Nv(i) = ucopy(i)*vcopy(i)**2 - (F+k)*vcopy(i)
  end do

  !(12.2, W1)
  call daxpy(n,tau/2.0_8,Nu,1,u,1)
  call daxpy(n,tau/2.0_8,Nv,1,v,1)

end subroutine heun

