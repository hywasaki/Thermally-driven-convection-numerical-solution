module tridiag_mod
  use iso_fortran_env, only: dp=>real64  
  implicit none
contains
!--------------------------------------------------------------------
! 1)  tri_dig   ---  tridiagonal solver
!--------------------------------------------------------------------
  subroutine tri_dig(a,b,c,r,x,n)
    use iso_fortran_env, only: dp=>real64  
    implicit none
    integer,           intent(in)  :: n
    real(dp),          intent(in)  :: a(n),b(n),c(n),r(n)
    real(dp),          intent(out) :: x(n)
    real(dp), allocatable :: beta(:), rho(:)
    integer :: j
    allocate(beta(n), rho(n))

    beta(1) = b(1);  if (abs(beta(1)) < 1.0e-12_dp) error stop "beta(1)=0"
    rho(1)  = r(1)

    do j = 2, n
       if (abs(beta(j-1)) < 1.0e-12_dp) error stop "beta(j-1)=0"
       rho(j)  = r(j) - a(j)*rho(j-1)/beta(j-1)
       beta(j) = b(j) - a(j)*c(j-1)/beta(j-1)
    end do
    if (abs(beta(n)) < 1.0e-16_dp) error stop "beta(n)=0"

    x(n) = rho(n)/beta(n)
    do j = n-1, 1, -1
       if (abs(beta(j)) < 1.0e-12_dp) error stop "beta(j)=0"
       x(j) = (rho(j) - c(j)*x(j+1))/beta(j)
    end do
    deallocate(beta, rho)
  end subroutine tri_dig
!--------------------------------------------------------------------
! 2)  check_nan   ---  detiene si encuentra NaN
!--------------------------------------------------------------------
  subroutine check_nan(arr, name, iter)
    real(dp), intent(in) :: arr(:)
    character(len=*), intent(in) :: name
    integer,       intent(in) :: iter
    if (any(.not.(arr == arr))) then
       write(*,*) "NaN detected in", trim(name), "at iteration", iter
       ! error stop
    end if
  end subroutine check_nan
!--------------------------------------------------------------------
! 3)  debug_coef  ---  imprime coeficientes seleccionados
!--------------------------------------------------------------------
  subroutine debug_coef(i, a, b, c, d, n)
    integer, intent(in) :: i, n
    real(dp), intent(in) :: a(n), b(n), c(n), d(n)
    integer :: j, mid
    print *, '=== Debug coef (iter=', i, ') ==='
    do j = 1, min(3,n)
       print '(I5,4F20.12)', j, a(j), b(j), c(j), d(j)
    end do
    if (n >= 7) then
       mid = n/2
       do j = mid-1, mid+1
          print '(I5,4F20.12)', j, a(j), b(j), c(j), d(j)
       end do
    end if
    do j = max(1,n-2), n
       print '(I5,4F20.12)', j, a(j), b(j), c(j), d(j)
    end do
    print *, '==============================='
  end subroutine debug_coef
!--------------------------------------------------------------------
! 4)  check_diag_dom  ---  comprueba dominancia diagonal
!--------------------------------------------------------------------
  subroutine check_diag_dom(tag,a,b,c,n,it)
    character(len=*),intent(in) :: tag
    integer, intent(in) :: n, it
    real(dp), intent(in) :: a(n), b(n), c(n)
    integer :: i, bad, shown
    bad = 0; shown = 0
    do i = 1, n
       if (abs(b(i)) < abs(a(i)) + abs(c(i))) then
          bad = bad + 1
          if (shown < 5) then
             write(*,'(A,": iter=",I0," fila=",I0)') trim(tag), it, i
             shown = shown + 1
          end if
       end if
    end do
    if (bad > 0) write(*,'("   total filas no dominantes: ",I0)') bad
  end subroutine check_diag_dom
!--------------------------------------------------------------------
   function diff2Cent(x) result (ddx)
            real(dp), intent(in):: x(:)
            real(dp) :: ddx(size(x)-2)
            integer :: im
            
            im = size(x)
            ddx(:) = x(3:im) - 2* x(2: im-1) + x(1:im-2)
   end function
end module tridiag_mod
