subroutine initial(n,u,v)
  ! ====================================
  ! Gaussian initial conditions
  ! u = 1 - exp(-(x-1/2)^2/(2s^2)) / 2,
  ! v = exp(-(x-1/2)^2/(2s^2)) / 4,
  ! with s = 0.1
  ! ====================================

  implicit none
  ! Variables as per assignment
  integer, intent(in) :: n
  real(kind=8), dimension(n), intent(out) :: u,v
  
  integer :: i         ! loop index
  real(kind=8) :: h    ! spatial step size

  h = 1.0_8/n
  
  do i=1,n
     u(i) = 1.0_8 - 0.5_8*exp(-(i*h - 0.5_8)**2/(2.0_8*0.1_8**2))
     v(i) = 0.25_8*exp(-(i*h - 0.5_8)**2/(2.0_8*0.1_8**2))
  end do
  
end subroutine initial
