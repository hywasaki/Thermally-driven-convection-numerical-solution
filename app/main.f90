program main
  use iso_fortran_env, only: dp=>real64
  implicit none
  
  !---------------- parametros ------------------------------------------
  integer,parameter :: nx=900, nt=9000000
  real(dp),parameter :: dx=0.5_dp, dt=4.0e-5_dp
  real(dp),parameter :: Le=20._dp, Sc=400._dp, Ra=2._dp, q=0.15_dp
  real(dp),parameter :: v0=1._dp/sqrt(2._dp)-0.0023_dp
  integer, parameter :: noise_ini = 400, noise_end=450

  !---- snapshots -------------------------------------------------------
  integer, parameter :: nsave   = 5              ! nº de archivos a guardar
  integer, parameter :: stride = nt/(nsave-1)

  !---- Fields ----------------------------------------------------------
  !---- No dynamic arrays, hardcoding lenghts
  real(dp) :: C_old(nx),C_new(nx),T_old(nx),T_new(nx)
  real(dp) :: Cq_old(nx),Cq_new(nx),Tq_old(nx),Tq_new(nx)
  real(dp) :: Wq_old(nx),Wq_new(nx),Psiq_old(nx),Psiq_new(nx)

  real(dp) :: C_save(nsave,nx),T_save(nsave,nx),Cq_save(nsave,nx),Tq_save(nsave,nx)
  real(dp) :: Wq_save(nsave,nx),Psiq_save(nsave,nx)

  !----- coefficients ---------------------------------------
  real(dp):: a_c(nx),b_c(nx),c_c(nx),d_c(nx)
  real(dp):: a_t(nx),b_t(nx),c_t(nx),d_t(nx)
  real(dp):: a_o(nx),b_o(nx),c_o(nx),d_o(nx)
  real(dp):: a_p(nx),b_p(nx),c_p(nx),d_p(nx)

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
  ! 1. Initializing arrays
  !======================================================================
  C_old=0._dp 
  T_old=0._dp
  C_old(1:noise_ini)=1._dp
  T_old(1:noise_ini)=1._dp

  !! This is only for seeding the random number generator
  call random_seed(size=nseed)
  allocate(seed(nseed)); 
  seed=12345; 
  call random_seed(put=seed); 
  deallocate(seed)

  !======================================================================
  ! 2. Initial perturbation
  !======================================================================
  do i=noise_ini+1, noise_end
    call random_number(r) 
    C_old(i) = 1.0e-3_dp * r
    call random_number(r) 
    T_old(i) = 1.0e-3_dp*r
  enddo

  !======================================================================
  ! 3. Diagonals
  !======================================================================


  
end program main
