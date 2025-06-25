program reaction_diffusion_convection
  use precision
  use tridiag_mod
  implicit none

  !-----------------------------------------------------------------
  ! Parameters
  integer, parameter :: nx     = 2000
  integer, parameter :: nxm2   = nx - 2
  
  integer, parameter :: nxm    = 1000
  integer, parameter :: itmax  = 2000000
  real(dp), parameter :: dx    = 0.5_dp
  real(dp), parameter :: dt    = 4.0e-4_dp
  real(dp), parameter :: sc    = 400.0_dp
  real(dp), parameter :: xle   = 20.0_dp
  real(dp), parameter :: ra    = 2.0_dp
  real(dp), parameter :: q     = 0.15_dp
  real(dp), parameter :: v0    = 1.0_dp / sqrt(2.0_dp)
  real(dp), parameter :: vadv  = v0 - 0.0023_dp

  !-----------------------------------------------------------------
  ! Fields
  real(dp) :: phi(nx),   phip(nx)
  real(dp) :: phi1(nx),  phi1p(nx)
  real(dp) :: th(nx),    thp(nx)
  real(dp) :: th1(nx),   th1p(nx)
  real(dp) :: ww(nx),    wwp(nx)
  real(dp) :: psi(nx)

  ! Arrays for tridiagonal solver
  real(dp) :: a(nxm2), b(nxm2), c(nxm2), rr(nxm2), u(nxm2)

  ! Local variables
  integer :: i, it, idx
  real(dp) :: dfpp, r
  character(len=32) :: fname
  !-----------
  real(dp) :: t0, tcur, pct
  integer  :: prog_step


  !-----------------------------------------------------------------
  ! Read initial φ and θ from external files
  open(unit=10, file='phi_init.txt', status='old', action='read')
  read(10, *) (phi(i), i = 1, nx)
  close(10)


  open(unit=11, file='th_init.txt', status='old', action='read')
  read(11, *) (th(i),  i = 1, nx)
  close(11)

  ! Initialize all other fields to zero
  phip  = 0.0_dp
  phi1  = 0.0_dp
  phi1p = 0.0_dp
  thp   = 0.0_dp
  th1   = 0.0_dp
  th1p  = 0.0_dp
  psi   = 0.0_dp
  ww    = 0.0_dp
  wwp   = 0.0_dp

  ! Seed random number generator
  call random_seed()

  ! Add noise around the front at i = nxm
  do i = nxm-40, nxm+40
    call random_number(r)
    phi1(i) = 1.0e-3_dp * r
    call random_number(r)
    th1(i)  = 1.0e-3_dp * r
  end do

  !-----------------------------------------------------------------
  ! Main time-stepping loop
  call cpu_time(t0)
  prog_step = itmax / 100

  do it = 1, itmax

    !---- explicit update for φ, φ'
    do i = 2, nx-1
      phip(i)=phi(i)+dt*(phi(i+1)+phi(i-1)-2.*phi(i))/dx**2+phi(i)**2*(1-phi(i))*dt
      phip(i)=phip(i)+dt*vadv*(phi(i+1)-phi(i-1))/(2.*dx)

      dfpp = phi(i)*(2.0_dp - 3.0_dp*phi(i))

      phi1p(i) = phi1(i) + dt*(phi1(i+1) + phi1(i-1) - 2.0_dp*phi1(i))/dx**2
      phi1p(i) = phi1p(i) + dt*vadv*(phi1(i+1) - phi1(i-1))/(2.0_dp*dx)
      phi1p(i) = phi1p(i) + dt*dfpp*phi1(i)
      phi1p(i) = phi1p(i) - dt*q**2*phi1(i)
      phi1p(i) = phi1p(i) - dt*q*psi(i)*(phi(i+1) - phi(i-1))/(2.0_dp*dx)

      thp(i)=th(i)+dt*xle*(th(i+1)+th(i-1)-2.*th(i))/dx**2+phi(i)**2*(1-phi(i))*dt
      thp(i)=thp(i)+dt*vadv*(th(i+1)-th(i-1))/(2.*dx)

      th1p(i) = th1(i) + dt*xle*(th1(i+1) + th1(i-1) - 2.0_dp*th1(i))/dx**2
      th1p(i) = th1p(i) + dt*vadv*(th1(i+1) - th1(i-1))/(2.0_dp*dx)
      th1p(i) = th1p(i) + dt*dfpp*phi1(i) - dt*xle*q**2*th1(i)
      th1p(i) = th1p(i) - dt*q*psi(i)*(th(i+1) - th(i-1))/(2.0_dp*dx)
    end do

    !---- enforce boundary conditions
    phi   = phip
    th    = thp
    phi1  = phi1p
    th1   = th1p

    phi(1)   = 1.0_dp;   phi(nx)   = 0.0_dp
    th(1)    = th(2);    th(nx)    = 0.0_dp
    phi1(1)  = 0.0_dp;   phi1(nx)  = 0.0_dp
    th1(1)   = th1(2);   th1(nx)   = 0.0_dp

    !-----------------------------------------------------------------
    ! Solve for vorticity ω with tridiagonal (implicit)
    do i = 2, nx-1
      rr(i-1) = -q*ra*th1(i) - ww(i)/(sc*dt)
      a(i-1)  = vadv/(sc*2.0_dp*dx) + 1.0_dp/dx**2
      b(i-1)  = -2.0_dp/dx**2 - q**2 - 1.0_dp/(sc*dt)
      c(i-1)  = 1.0_dp/dx**2 - vadv/(sc*2.0_dp*dx)
    end do
    a(1)      = 0.0_dp
    c(nxm2)   = 0.0_dp
    call tri_dig(a, b, c, rr, u, nxm2)
    do i = 2, nx-1
      wwp(i) = u(i-1)
    end do
    wwp(1) = 0.0_dp; wwp(nx) = 0.0_dp
    ww     = wwp

    !-----------------------------------------------------------------
    ! Solve for stream function ψ with tridiagonal
    do i = 2, nx-1
      rr(i-1) = ww(i)
      a(i-1)  = 1.0_dp/dx**2
      b(i-1)  = -2.0_dp/dx**2 - q**2
      c(i-1)  = 1.0_dp/dx**2
    end do
    a(1)      = 0.0_dp
    c(nxm2)   = 0.0_dp
    call tri_dig(a, b, c, rr, u, nxm2)
    do i = 2, nx-1
      psi(i) = u(i-1)
    end do
    psi(1)  = 0.0_dp; psi(nx) = 0.0_dp

    !-----------------------------------------------------------------
    ! Every 200 000 steps, write out φ, θ, ψ, φ₁, θ₁
    if (mod(it, 200000) == 0) then
      idx = it / 200000
      write(fname, '(A,I0,A)') 'output_', idx, '.dat'
      open(unit=20, file=trim(fname), status='replace', action='write')
      do i = 1, nx
        write(20,'(I6, 6(1X, ES12.5))') i, phi(i), th(i), psi(i), phi1(i), th1(i), ww(i)
      end do
      close(20)
    end if
    if (mod(it, prog_step) == 0) then
      call cpu_time(tcur)
      pct = 100.0_dp * real(it, dp) / real(itmax, dp)
      write(*,'("Progreso: ", F5.1, " %   t=", F8.1, " s")') pct, tcur - t0
    end if


  end do  ! end time loop

end program reaction_diffusion_convection
