module precision
   implicit none
   integer, parameter :: dp = selected_real_kind(15,307)
end module precision


module tridiag_mod
   use precision
   implicit none
contains
   subroutine tri_dig(a,b,c,r,x,n)
      use precision
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in)  :: a(n),b(n),c(n),r(n)
      real(dp), intent(out) :: x(n)

      real(dp), allocatable :: beta(:), rho(:)
      integer :: j
      allocate(beta(n), rho(n))

      beta(1) = b(1)
      if (abs(beta(1)) < 1.0e-12_dp) then
         print *, "ERROR: beta(1) = 0 en tri_dig"
         stop
      end if

      rho(1) = r(1)

      do j = 2, n
         if (abs(beta(j-1)) < 1.0e-12_dp) then
            print *, "ERROR: beta(", j-1, ") = 0 en tri_dig"
            stop
         end if
         rho(j) = r(j) - a(j) / beta(j-1) * rho(j-1)
         beta(j) = b(j) - a(j) * c(j-1) / beta(j-1)
      end do

      if (abs(beta(n)) < 1.0e-16_dp) then
         print *, "ERROR: beta(n) = 0 en tri_dig"
         stop
      end if

      x(n) = rho(n) / beta(n)

      do j = n-1, 1, -1
         if (abs(beta(j)) < 1.0e-12_dp) then
            print *, "ERROR: beta(", j, ") = 0 en backward solve"
            stop
         end if
         x(j) = (rho(j) - c(j) * x(j+1)) / beta(j)
      end do

      deallocate(beta, rho)
   end subroutine tri_dig

   subroutine check_diag_dom(tag,a,b,c,n,it)
      character(len=*),intent(in) :: tag
      integer,intent(in) :: n,it
      real(dp),intent(in) :: a(n),b(n),c(n)
      integer :: i,bad,shown
      bad=0; shown=0
      do i=1,n
         if(abs(b(i))<abs(a(i))+abs(c(i)))then
            bad=bad+1
            if(shown<5)then
               write(*,'(A,": iter=",I0," fila=",I0)') trim(tag),it,i
               shown=shown+1
            end if
         end if
      end do
      if(bad>0) write(*,'("   total filas no dominantes: ",I0)') bad
   end subroutine check_diag_dom

   subroutine debug_coef(i, a, b, c, d, n)
      use precision
      implicit none
      integer, intent(in) :: i, n
      real(dp), intent(in) :: a(n), b(n), c(n), d(n)
      integer :: j, mid

      print *, '=== Debug coef (iteración =', i, ') ==='
      print *, '   j     |a(j)|                |b(j)|                |c(j)|                |d(j)|'

      ! First 3
      do j = 1, min(3,n)
         print '(I5, 4F20.12)', j, a(j), b(j), c(j), d(j)
      end do

      ! Middle 3
      if (n >= 7) then
         mid = n/2
         do j = mid-1, mid+1
            print '(I5, 4F20.12)', j, a(j), b(j), c(j), d(j)
         end do
      end if

      ! Last 3
      do j = max(1,n-2), n
         print '(I5, 4F20.12)', j, a(j), b(j), c(j), d(j)
      end do

      print *, '================================================================='
   end subroutine debug_coef
end module tridiag_mod

module nancheck_mod
  use precision
  implicit none
contains

  subroutine check_nan(arr, name, iter)
    real(dp), intent(in) :: arr(:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: iter

    if (any(.not.(arr == arr))) then
      write(*,*) "NaN detected in", trim(name), "at iteration", iter
      stop
    end if
  end subroutine check_nan

end module nancheck_mod
program main
  ! use iso_fortran_env, only: dp=>real64   
  use precision
  use tridiag_mod
  use nancheck_mod
  implicit none
  !---------------- parametros ------------------------------------------
  integer,parameter :: nx=2000, nt=2000000
  real(dp),parameter :: dx=0.5_dp, dt=4.0e-5_dp
  real(dp),parameter :: Le=20._dp, Sc=400._dp, Ra=2._dp, q=0.15_dp
  real(dp),parameter :: v0=1._dp/sqrt(2._dp)-0.0023_dp
  integer, parameter :: noise_ini = 1000, noise_end=1500

  !---- snapshots -------------------------------------------------------
  integer, parameter :: nsave   = 5              ! 'nsave' total saves
  integer, parameter :: stride = nt/(nsave-1)

  !---- Fields ----------------------------------------------------------
  !---- No dynamic arrays, hardcoding lenghts
  real(dp) :: C_old(nx),C_new(nx),T_old(nx),T_new(nx)
  real(dp) :: Cq_old(nx),Cq_new(nx),Tq_old(nx),Tq_new(nx)
  real(dp) :: Wq_old(nx),Wq_new(nx),Psiq_old(nx),Psiq_new(nx)

  real(dp) :: C_save(nsave,nx),T_save(nsave,nx),Cq_save(nsave,nx),Tq_save(nsave,nx)
  real(dp) :: Wq_save(nsave,nx),Psiq_save(nsave,nx)

  !----- coefficients ---------------------------------------
  real(dp):: a_con(nx),  b_con(nx),  c_con(nx),  d_con(nx)     ! Concentration_q
  real(dp):: a_temp(nx), b_temp(nx), c_temp(nx), d_temp(nx)    ! Temperature_q 
  real(dp):: a_vort(nx), b_vort(nx), c_vort(nx), d_vort(nx)    ! Vorticity
  real(dp):: a_str(nx),  b_str(nx),  c_str(nx),  d_str(nx)     ! Stream Function
  ! Arrays of size nx
  ! Arrays a and c are logically padded: a(1) = 0 and c(nx-2) = 0

  !!---------------- auxiliares ------------------------------------------
  integer :: i,save_cnt,prog_step
  real(dp) :: t0,tcur,pct
  integer  :: nseed,nxm,left,right
  integer, allocatable :: seed(:)
  real(dp) :: r


  !! Check what is declared as parameters
  write(*, *) "Reaction difussion" 
  write(*,*) "model thermally driven convection"

  write(*, *) "Mesh parameters"
  write(*,*) "Number of points"
  write(*, "(2(A8,I10, 3X))") "Space: ",nx,"Time: ", nt
  write(*,*) "Space between points"
  write(*, "(2(A8,F15.7, 3X))") "Space: ",dx,"Time: ", dt
  !! TODO: complete with other parameters (including snapshots).

  !======================================================================
  ! Initializing arrays
  !======================================================================
  C_old=0._dp 
  T_old=0._dp
  C_old(1:noise_ini)=1._dp
  T_old(1:noise_ini)=1._dp

   Cq_old = 0._dp
   Tq_old = 0._dp
   Wq_old = 0._dp
   Psiq_old = 0._dp


  !! This is only for seeding the random number generator
  call random_seed(size=nseed)
  allocate(seed(nseed)); 
  seed=12345; 
  call random_seed(put=seed); 
  deallocate(seed)

  !======================================================================
  ! Initial perturbation
  !======================================================================
  do i=noise_ini+1, noise_end
    call random_number(r); C_old(i) = 1.0e-3_dp * r
    call random_number(r); T_old(i) = 1.0e-3_dp * r
    call random_number(r); Cq_old(i) = 1.0e-3_dp * r
    call random_number(r); Tq_old(i) = 1.0e-3_dp * r
  enddo

  ! C0 = [1,1... random...0,0]
  ! T0 = [1,1... random...0,0]
  ! Cq = [0,0... random...0,0]
  ! Tq = [0,0... random...0,0]
  ! Wq = [0,0...0,0]
  ! Psiq = [0,0...,0,0]

  !======================================================================
  ! Diagonals (part 1)
  !======================================================================
  ! Initialazing arrays and adding padding for a and c  a(1) = 0 and c(nx-2) = 0
  a_con=0._dp ; b_con=1._dp;  c_con=0._dp;  d_con=0._dp
  a_temp=0._dp ;b_temp=1._dp; c_temp=0._dp; d_temp=0._dp
  a_vort=0._dp ;b_vort=1._dp; c_vort=0._dp; d_vort=0._dp
  a_str=0._dp ; b_str=1._dp;  c_str=0._dp ; d_str=0._dp

  ! Compute coefficients that do NOT  change in time
  ! Concentration (Con)
  a_con(2:nx-1)=-dt*(1/dx**2 - v0/(2*dx))
  c_con(2:nx-1)=-dt*(1/dx**2 + v0/(2*dx))

  ! Temperature (Temp)
   a_temp(2:nx-1)=-dt*(Le/dx**2 - v0/(2*dx))
   b_temp(2:nx-1)=1+dt*(Le*q**2+2*Le/dx**2)
   c_temp(2:nx-1)=-dt*(Le/dx**2 + v0/(2*dx))

  ! Vorticity (Vort)
   a_vort(2:nx-1)=-dt*(Sc/dx**2 - v0/(2*dx))
   b_vort(2:nx-1)=1+dt*(Sc*q**2+2*Sc/dx**2)
   c_vort(2:nx-1)=-dt*(Sc/dx**2 + v0/(2*dx))

  ! Stream function (Str)
   a_str(2:nx-1)= 1/dx**2
   b_str(2:nx-1)=-2/dx**2 - q**2
   c_str(2:nx-1)= 1/dx**2

  call cpu_time(t0)          ! Start timing the simulation
  prog_step = nt / 100       ! Print progress every 1% of total steps
  save_cnt  = 0              ! Counter for saving data 

  !======================================================================
  ! Time loop
  !======================================================================
  do i=0,nt
    !======================================================================
    ! Solving for flat front solutions C0 and T0 using Euler explicit methods
    ! C0
    C_new(2:nx-1) = C_old(2:nx-1)+ &
                  dt/dx**2*(C_old(3:nx)-2*C_old(2:nx-1)+C_old(1:nx-2))+ &
                  dt*(C_old(2:nx-1)**2)*(1-C_old(2:nx-1))+&
                  v0*dt/(2*dx)*(C_old(3:nx)-C_old(1:nx-2))
    ! T0
    T_new(2:nx-1)=T_old(2:nx-1)+&
                  Le*dt/dx**2*(T_old(3:nx)-2*T_old(2:nx-1)+T_old(1:nx-2))+ &
                  dt*(C_old(2:nx-1)**2)*(1-C_old(2:nx-1))+&
                  v0*dt/(2*dx)*(T_old(3:nx)-T_old(1:nx-2))
    ! Update
    C_old = C_new
    T_old = T_new  
    ! C0 and T0 Boundary conditions 
    ! C0(1)=1 = 1, C0(nx)=0 (DIRICHLET)
    C_old(1) = 1; C_old(nx) = 0

    ! dT/dx = 0 (Von Neumann) and T0(nx) =0 (DIRICHLET)
    T_old(1) = T_old(2); T_old(nx) = 0
    !======================================================================
    ! Diagonals (part 2)
    !======================================================================
    ! Compute coefficients that DO change in time
    ! Concentration (Con)
    b_con(2:nx-1)=1+dt*(q**2+2/dx**2 + C_old(2:nx-1)*(2-3*C_old(2:nx-1)))
    d_con(2:nx-1)=Cq_old(2:nx-1)+dt*q*Psiq_old(2:nx-1)/(2*dx)*(C_old(3:nx)-C_old(1:nx-2))

    ! Temperature (Temp)
    d_temp(2:nx-1)=Tq_old(2:nx-1)+dt*q*Psiq_old(2:nx-1)/(2*dx)*(T_old(3:nx)-T_old(1:nx-2))+ &
                dt*C_old(2:nx-1)*(2-3*C_old(2:nx-1))*Cq_old(2:nx-1)

    ! Vorticity (Vort)
    d_vort(2:nx-1)=Wq_old(2:nx-1)-dt*q*Ra*Sc*Tq_old(2:nx-1)

    ! Stream function (Str)
    d_str(2:nx-1)=Wq_old(2:nx-1)
    !======================================================================
    ! Boundary conditions
    ! Cq
    d_con(2) = d_con(2) - a_con(2); a_con(2) = 0._dp; d_con(1) = 1.0_dp! Cq(1) = 1
    c_con(nx-1) = 0.0_dp; d_con(nx)= 0.0_dp ! Cq(nx) = 0 

    ! Tq
    b_temp(2) = b_temp(2) + a_temp(2); a_temp(2) = 0.0_dp; c_temp(1) = -1.0_dp; d_temp(1) = 0.0_dp ! dT/dx (0) = 0
    c_temp(nx-1)=  0.0_dp; d_temp(nx)= 0.0_dp! Tq(nx) = 0 

    ! Voricity 
    a_vort(2)  = 0.0_dp; d_vort(1) = 0.0_dp ! Wq(0)  = 0 
    c_vort(nx-1) = 0.0_dp; d_vort(nx)= 0.0_dp ! Wq(nx) = 0 

    ! Stream function
    a_str(2)  = 0.0_dp; d_str(1) = 0.0_dp ! Psi(0) = 0 
    c_str(nx-1) = 0.0_dp; d_str(nx)= 0.0_dp ! Psi(nx) = 0 
   
   if (i == 1) then
      print *
      print *, "Coeficientes de concentración (Cq):"
      call debug_coef(i, a_con, b_con, c_con, d_con, nx)

      print *, "Coeficientes de temperatura (Tq):"
      call debug_coef(i, a_temp, b_temp, c_temp, d_temp, nx)

      print *, "Coeficientes de vorticidad (Wq):"
      call debug_coef(i, a_vort, b_vort, c_vort, d_vort, nx)

      print *, "Coeficientes de función de corriente (Psiq):"
      call debug_coef(i, a_str, b_str, c_str, d_str, nx)
      print *
   end if

    ! The arrays must be of the form (padding added for a and c):
    ! a_con=[0, 0, a2 ... a_n-1,0]  ; b_con=[1, b1... b_n-1, 1]        ; c_con=[0,c1...c_n-2, 0, 0]     ; d_con=[1, d1-a1, d2...d_n-1, 0]
    ! a_temp=[0, 0, a2 ... a_n-1,0] ; b_temp=[1, b1+a1, b2...b_n-1, 1] ; c_temp =[-1, c1,...c_n-2,0,0]  ; d_temp=[0,d1...d_n-1,0]
    ! a_vort=[0, 0, a2 ... a_n-1,0] ; b_vort=[1, b1... b_n-1, 1]       ; c_vort=[0,c1...c_n-2, 0, 0]    ; d_vort=[0, d1-a1, d2...d_n-1, 0]
    ! a_str=[0, 0, a2 ... a_n-1,0]  ; b_str=[1, b1... b_n-1, 1]        ; c_str=[0,c1...c_n-2, 0, 0]     ; d_str=[0, d1-a1, d2...d_n-1, 0]

    ! Check if matrix is diagonally dominant: |b(i)| ≥ |a(i)| + |c(i)|
    call check_diag_dom("Eq11 Cq",a_con,b_con,c_con,nx,i)
    call check_diag_dom("Eq12 Tq",a_temp,b_temp,c_temp,nx,i)
    call check_diag_dom("Eq13 wq",a_vort,b_vort,c_vort,nx,i)
    call check_diag_dom("Eq14 psi",a_str,b_str,c_str,nx,i)

    ! Solving the Tridiagonal sistem 
    call tri_dig(a_con(1:nx),  b_con(1:nx),  c_con(1:nx),  d_con(1:nx),  Cq_new(1:nx)  ,nx)
    call tri_dig(a_temp(1:nx), b_temp(1:nx), c_temp(1:nx), d_temp(1:nx), Tq_new(1:nx)  ,nx)
    call tri_dig(a_vort(1:nx), b_vort(1:nx), c_vort(1:nx), d_vort(1:nx), Wq_new(1:nx)  ,nx)
    call tri_dig(a_str(1:nx),  b_str(1:nx),  c_str(1:nx),  d_str(1:nx),  Psiq_new(1:nx),nx)
    
    ! Check for NaN values
    call check_nan(Cq_new,   "Cq_new", i)
    call check_nan(Tq_new,   "Tq_new", i)
    call check_nan(Wq_new,   "Wq_new", i)
    call check_nan(Psiq_new, "Psiq_new", i)
    
    ! Update
    Cq_old=Cq_new; Tq_old=Tq_new
    Wq_old=Wq_new; Psiq_old=Psiq_new
    !======================================================================
    !-- Print progress
    if(mod(i,prog_step)==0)then
      call cpu_time(tcur)
      pct=100._dp*real(i,dp)/real(nt,dp)
      write(*,'("Progreso: ",F5.1," %   t=",F8.1," s")') pct,tcur-t0
    end if

    !-- Save snapshots
    if (mod(i, stride) == 0 .and. save_cnt < nsave) then
      save_cnt = save_cnt + 1
      C_save(save_cnt,:)     = C_old
      T_save(save_cnt,:)     = T_old
      Cq_save(save_cnt,:)    = Cq_old
      Tq_save(save_cnt,:)    = Tq_old
      Wq_save(save_cnt,:)    = Wq_old
      Psiq_save(save_cnt,:)  = Psiq_old
      write(*,'("Snapshot t=",1PE12.4,"  iter=",I0)') dt*real(i,dp), i
    end if
  end do

  !======================================================================
  !  Saving files
  open(10,file='C_save.txt' ,status='replace')
  open(11,file='T_save.txt' ,status='replace')
  open(12,file='Cq_save.txt',status='replace')
  open(13,file='Tq_save.txt',status='replace')
  open(14,file='Wq_save.txt',status='replace')
  open(15,file='Psiq_save.txt',status='replace')

  do i = 1, save_cnt
    write(10,*) C_save (i,:)
    write(11,*) T_save (i,:)
    write(12,*) Cq_save(i,:)
    write(13,*) Tq_save(i,:)
    write(14,*) Wq_save(i,:)
    write(15,*) Psiq_save(i,:)
  end do
  close(10); close(11); close(12); close(13); close(14); close(15)
end program main
