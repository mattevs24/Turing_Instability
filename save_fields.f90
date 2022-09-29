subroutine save_fields(n,u,v,filename)
  ! ===========================================================
  ! Save solutions to disk.
  ! Writes a text file which contains the values of the vectors
  ! u and v at the nodal points.
  ! ===========================================================

  implicit none
  ! problem size
  integer, intent(in) :: n
  ! Velocity field $u,v$
  real(kind=8), dimension(n), intent(in) :: u,v
  ! Name of file to write to
  character(len=*), intent(in) :: filename
  ! Loop variable
  integer :: i
  ! File handle
  integer,parameter :: file_id = 99
                       ! avoid conflicts if file of parameters is not closed
  
  ! Write first component
  open(unit=file_id,file=trim(filename))
  do i=1, n
     write(file_id,'(E20.8e3," ")',advance='no') u(i)
  end do
  write(file_id,'("")') ! delimiter
  ! Write second component
  do i=1, n
     write(file_id,'(E20.8e3," ")',advance='no') v(i)
  end do
  write(file_id,'("")')
  close(file_id)

end subroutine save_fields
