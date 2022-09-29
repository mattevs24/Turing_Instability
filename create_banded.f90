subroutine create_banded(n,tau,Blu,Blv,Bru,Brv,hu,hv) 

  implicit none

  integer, intent(in) :: n
  real(kind=8), intent(in) :: tau
  real(kind=8), intent(out), dimension(2,n) :: Blu, Blv, Bru, Brv
  real(kind=8), intent(out), dimension(n) :: hu, hv

  !internal variables

  integer :: i, info
  real(kind=8) :: tauu, tauv, dot
  real(kind=8), dimension(n) :: c, copy1
  real(kind=8), dimension(n-1) :: copy2
  real(kind=8), external :: ddot
  
  !Diffusion coeffcients
  real(kind=8), parameter :: Du = 2d-5
  real(kind=8), parameter :: Dv = 1d-5


  !calculation of tau_u and tau_v
  tauu = Du*tau*0.25*n**2
  tauv = Dv*tau*0.25*n**2

  !calculating the matrices B** with an if statement to edit accordingly for the first and last entries
  do i = 1,n
     Blu(1,i) = -tauu
     Blv(1,i) = -tauv
     Bru(1,i) = tauu
     Brv(1,i) = tauv
     if ((i == 1).OR.(i == n)) then
        Blu(2,i) = 1 + 3*tauu
        Blv(2,i) = 1 + 3*tauv
        Bru(2,i) = 1 - 3*tauu
        Brv(2,i) = 1 - 3*tauv
        hu(i) = 1.0_8
     else
        Blu(2,i) = 1 + 2*tauu
        Blv(2,i) = 1 + 2*tauv
        Bru(2,i) = 1 - 2*tauu
        Brv(2,i) = 1 - 2*tauv
        hu(i) = 0.0_8
     end if
  end do

  !The start of the calculation of h_u and h_v

  !copies the vector hu (which is c) to hv and c itsself
  call dcopy(n,hu,1,hv,1)
  call dcopy(n,hu,1,c,1)

  !calculating the vectors h*
  !Bl* are spd by gergorins circle theorem, so used dptsv for this.

  !First copy each part of Blu to each copy vector to not overwrite Blu itself and then solve using dptsv
  call dcopy(n,Blu(2,1),2,copy1,1)
  call dcopy((n-1),Blu(1,2),2,copy2,1)
  call dptsv(n,1,copy1,copy2,hu,n,info)

  !first copy each part of Blv to each copy vector to not overwrite Blv itself and then solve using dptsv
  call dcopy(n,Blv(2,1),2,copy1,1)
  call dcopy((n-1),Blv(1,2),2,copy2,1)
  call dptsv(n,1,copy1,copy2,hv,n,info)

  !calculate the scalar for the vector hu and then scale hu by said product
  dot = tauu/(1-tauu*ddot(n,c,1,hu,1))
  call dscal(n,dot,hu,1)

  !calculate the scalar for the vector hv and then scale hv by said product
  dot = tauv/(1-tauv*ddot(n,c,1,hv,1))
  call dscal(n,dot,hv,1)

end subroutine create_banded
