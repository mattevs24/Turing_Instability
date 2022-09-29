subroutine test_banded(n,tau,Blu,Blv,Bru,Brv,hu,hv) 

  implicit none

  integer, intent(in) :: n
  real(kind=8), intent(in) :: tau
  real(kind=8), intent(in), dimension(n,n) :: Blu, Blv, Bru, Brv
  real(kind=8), intent(in), dimension(n) :: hu, hv

  !internal variables

  real(kind=8), dimension(n) :: utest, vtest, uin, vin, work
  integer, dimension(n) :: ipiv
  integer :: info
  real(kind=8), dimension(n,n) :: Mrcopy
  real(kind=8), dimension(n,n) :: Mlu,Mlv,Mru,Mrv

  !interface the function dnrm2
  interface
     function dnrm2(N,X,INCX) result(dnorm)
       integer, intent(in) :: N, INCX
       real(kind=8), dimension(N), intent(in) :: X
       real(kind=8) :: dnorm
     end function dnrm2
  end interface

  !create random vectors to test
  call random_number(utest)
  call random_number(vtest)


  !similarly to test_crank, input the matrices Mlu,Mlv,Mru and Mrv to be able to reverse the utest and vtest vectors to find the appropriate input vectors for banded_crank
  call create_matrices(n,tau,Mlu,Mlv,Mru,Mrv)

  call dsymv('U',n,1.0_8,Mlu,n,utest,1,0.0_8,uin,1)
  call dsymv('U',n,1.0_8,Mlv,n,vtest,1,0.0_8,vin,1)
  
  call dcopy(n**2,Mru,1,Mrcopy,1)
  call dsysv('U',n,1,Mrcopy,n,ipiv,uin,n,work,n,info)

  call dcopy(n**2,Mrv,1,Mrcopy,1)
  call dsysv('U',n,1,Mrcopy,n,ipiv,vin,n,work,n,info)

  !call banded crank and find the differnce in the methods

  call banded_crank(n,tau,Blu,Blv,Bru,Brv,hu,hv,uin,vin)

  !calculate the vectors differences to output them into uin and vin resp.
  call daxpy(n,-1.0_8,utest,1,uin,1)
  call daxpy(n,-1.0_8,vtest,1,vin,1)
  
  !print the norms of uin and vin
  print*, 'The u error is:', dnrm2(n,uin,1)
  print*, 'The v error is:', dnrm2(n,vin,1)


end subroutine test_banded
