subroutine create_matrices(n,tau,Mlu,Mlv,Mru,Mrv)
  ! =============================================================
  ! A subroutine for assembling Crank-Nicolson matrices, Eq. (13)
  ! =============================================================

  implicit none
  ! Variables as per assignment
  integer, intent(in) :: n
  real(kind=8), intent(in) :: tau
  real(kind=8), dimension(n,n), intent(out) :: Mlu,Mlv,Mru,Mrv

  ! Diffusion coefficients
  real(kind=8), parameter :: Du = 2d-5
  real(kind=8), parameter :: Dv = 1d-5
  
  ! Loop index
  integer :: i
  ! Courant numbers
  real(kind=8) :: tau_u, tau_v

  
  ! n**2 = 1/h**2, where h is the spatial step size  
  tau_u = Du*tau*0.25_8*n**2
  tau_v = Dv*tau*0.25_8*n**2

  ! These matrices are mostly zero, so easy to enter zeros first
  Mlu = 0.0_8
  Mlv = 0.0_8
  Mru = 0.0_8
  Mrv = 0.0_8

  ! Enter the diagonal values
  do i=1,n
     Mlu(i,i) = 1.0_8 + tau_u*2
     Mlv(i,i) = 1.0_8 + tau_v*2
     Mru(i,i) = 1.0_8 - tau_u*2
     Mrv(i,i) = 1.0_8 - tau_v*2
  end do

  ! Sub- and super-diagonal values
  do i=1,n-1
     Mlu(i,i+1) = -tau_u
     Mlu(i+1,i) = Mlu(i,i+1)
     Mlv(i,i+1) = -tau_v
     Mlv(i+1,i) = Mlv(i,i+1)
     Mru(i,i+1) =  tau_u
     Mru(i+1,i) = Mru(i,i+1)
     Mrv(i,i+1) =  tau_v
     Mrv(i+1,i) = Mrv(i,i+1)
  end do

  ! Periodic boundary conditions
  Mlu(1,n) = -tau_u
  Mlu(n,1) = Mlu(1,n)
  Mlv(1,n) = -tau_v
  Mlv(n,1) = Mlv(1,n)
  Mru(1,n) =  tau_u
  Mru(n,1) = Mru(1,n)
  Mrv(1,n) =  tau_v
  Mrv(n,1) = Mrv(1,n)
  
end subroutine create_matrices
