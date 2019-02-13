module mod_ICDW_HMC
  implicit none
  integer,parameter :: DP=kind(1.0d0)
  real(DP),parameter :: PI=2.0_DP*acos(0.0_DP)
  integer,parameter :: NSITE=16
  integer :: mu(1:2,1:2),pn(-1:NSITE+2)

contains

subroutine set_unit_and_periodic (mu,pn)
  implicit none
  integer,intent(inout) :: mu(1:2,1:2),pn(-1:NSITE+2)
  integer :: i

  mu(1,1) = 1
  mu(1,2) = 0
  mu(2,1) = 0
  mu(2,2) = 1

  do i = 1, NSITE
    pn(i) = i
  end do
  pn(-1) = NSITE-1
  pn(0) = NSITE
  pn(NSITE+1) = 1
  pn(NSITE+2) = 2

  return
end subroutine

subroutine site_check (i,j)
  implicit none
  integer,intent(in) :: i,j

  if ((i<1 .or. NSITE<i) .or. (j<1 .or. NSITE<j)) &
    & stop 'site index is out of lattice (1 <= i,j <= NSITE).'

  return
end subroutine

function gauss_rand () result (rand_num)
  use mt19937
  implicit none
  real(DP) :: rr,tt,rand_num

  rr = grnd()
  tt = grnd()
  rr = 1.0_DP - rr
  tt = 1.0_DP - tt
  rand_num = sqrt(-2.0_DP*log(rr))*cos(2.0_DP*PI*tt)

  return
end function

subroutine generate_state_with_gaussian (s1,s0,var)
  implicit none
  real(DP),intent(inout) :: s1(:,:)
  real(DP),intent(in) :: s0(:,:),var
  integer :: i,j

  do j = 1, NSITE
  do i = 1, NSITE
    s1(i,j) = gauss_rand()*sqrt(var) + s0(i,j)
  end do
  end do

  return
end subroutine

subroutine generate_momentum_with_gaussian (p)
  implicit none
  real(DP),intent(inout) :: p(:,:)
  integer :: i,j

  do j = 1, NSITE
  do i = 1, NSITE
    p(i,j) = gauss_rand()
  end do
  end do

  return
end subroutine

function hamil (sx,sp,vx,vp,mchi,lam,nu,gg) result (h)
  implicit none
  real(DP),intent(in) :: sp(1:NSITE,1:NSITE),sx(1:NSITE,1:NSITE)
  real(DP),intent(in) :: vp(1:NSITE,1:NSITE,1:2),vx(1:NSITE,1:NSITE,1:2)
  real(DP),intent(in) :: mchi,gg,lam,nu
  real(DP) :: h,div_s,div_v,sum_v
  integer :: i,j

  h = 0.0_DP

  do j = 1, NSITE
  do i = 1, NSITE
    sum_v = sum_vec(vp,i,j)   ! calculation of vector momentum
    div_s = div_sca(sx,i,j)   ! calculation of scalar derivative

    h = h + 0.5_DP*sp(i,j)**2 - 0.5_DP*sx(i,j)*div_s + 0.5_DP*sum_v + lam*(sx(i,j)**2 - nu**2)**2

    sum_v = sum_vec(vx,i,j)
    div_v = div_vec(vx,i,j)

    h = h + 0.5_DP*sum_v + 0.5_DP*(div_v - gg*sx(i,j))**2/mchi**2

  end do
  end do

  return
end function

function div_sca (s,i,j) result (div_s)
  !
  ! take derivative for scalar field
  !
  implicit none
  real(DP),intent(in) :: s(:,:)
  integer,intent(in) :: i,j
  real(DP) :: div_s
  integer :: k

  call site_check (i,j)

  div_s = 0.0_DP
  do k = 1, 2
    div_s = div_s + s(pn(i+mu(k,1)),pn(j+mu(k,2))) - 2.0_DP*s(i,j) + s(pn(i-mu(k,1)),pn(j-mu(k,2)))
  end do

  return
end function

function div_vec (v,i,j) result (div_v)
  implicit none
  real(DP),intent(in) :: v(1:NSITE,1:NSITE,1:2)
  integer,intent(in) :: i,j
  real(DP) :: div_v
  integer :: k

  call site_check (i,j)

  div_v = 0.0_DP
  do k = 1, 2
    div_v = div_v + v(pn(i + mu(k,1)),pn(j + mu(k,2)),k) - v(pn(i - mu(k,1)),pn(j - mu(k,2)),k)
  end do

  div_v = 0.5_DP*div_v

  return
end function

function sum_vec (v,i,j) result (sum_v)
  implicit none
  real(DP),intent(in) :: v(1:NSITE,1:NSITE,1:2)
  integer,intent(in) :: i,j
  real(DP) :: sum_v
  integer :: k

  call site_check(i,j)

  sum_v = 0.0_DP
  do k = 1, 2
    sum_v = sum_v + v(i,j,k)**2
  end do

  return
end function

subroutine leapfrog_xpx_md (tau,nmd,sx,sp,vx,vp,mchi,lam,nu,gg)
  implicit none
  real(DP),intent(in) :: tau,mchi,lam,nu,gg
  integer,intent(in) :: nmd
  real(DP),intent(inout) :: vx(1:NSITE,1:NSITE,1:2),vp(1:NSITE,1:NSITE,1:2)
  real(DP),intent(inout) :: sx(1:NSITE,1:NSITE),sp(1:NSITE,1:NSITE)
  real(DP) :: dt,dt2
  integer :: i

  dt = tau/nmd
  dt2 = dt/2.0_DP

  call update_x(dt2,sx,sp,vx,vp,mchi,lam,nu,gg)
  do i = 1, nmd-1
    call update_p(dt,sx,sp,vx,vp,mchi,lam,nu,gg)
    call update_x(dt,sx,sp,vx,vp,mchi,lam,nu,gg)
  end do
  call update_p(dt,sx,sp,vx,vp,mchi,lam,nu,gg)
  call update_x(dt2,sx,sp,vx,vp,mchi,lam,nu,gg)

  return
end subroutine

subroutine update_x (dt,sx,sp,vx,vp,mchi,lam,nu,gg)
  implicit none
  real(DP),intent(inout) :: sx(1:NSITE,1:NSITE),vx(1:NSITE,1:NSITE,1:2)
  real(DP),intent(in) :: sp(1:NSITE,1:NSITE),vp(1:NSITE,1:NSITE,1:2)
  real(DP),intent(in) :: dt,mchi,lam,nu,gg

  sx(:,:) = sx(:,:) + sp(:,:)*dt
  vx(:,:,:) = vx(:,:,:) + vp(:,:,:)*dt

  return
end subroutine

subroutine update_p (dt,sx,sp,vx,vp,mchi,lam,nu,gg)
  implicit none
  real(DP),intent(inout) :: sp(1:NSITE,1:NSITE),vp(1:NSITE,1:NSITE,1:2)
  real(DP),intent(in) :: sx(1:NSITE,1:NSITE),vx(1:NSITE,1:NSITE,1:2)
  real(DP),intent(in) :: dt,mchi,lam,nu,gg
  real(DP) :: sf(1:NSITE,1:NSITE),vf(1:NSITE,1:NSITE,1:2)
  real(DP) :: div_s(1:NSITE,1:NSITE),div_v(1:NSITE,1:NSITE,1:2)
  integer :: i,j

  !
  ! update scalar field
  !
  sf(:,:) = get_scalar_force(sx,vx,mchi,lam,nu,gg)
  sp(:,:) = sp(:,:) + sf(:,:)*dt

  !
  ! update vector field
  !
  vf(:,:,:) = get_vector_force(sx,vx,mchi,lam,nu,gg)
  vp(:,:,:) = vp(:,:,:) + vf(:,:,:)*dt

  return
end subroutine

function get_scalar_force(sx,vx,mchi,lam,nu,gg) result(sf)
  implicit none
  real(DP),intent(in) :: sx(1:NSITE,1:NSITE),vx(1:NSITE,1:NSITE,1:2),mchi,lam,nu,gg
  real(DP) :: sf(1:NSITE,1:NSITE)
  real(DP) :: div_s,div_v
  integer :: i,j

  do i = 1, NSITE
  do j = 1, NSITE
    div_s = div_sca(sx,i,j)
    div_v = div_vec(vx,i,j)

    sf(i,j) = div_s - 4.0_DP*lam*sx(i,j)*(sx(i,j)**2 - nu**2) + gg*(div_v - gg*sx(i,j))/mchi**2

  end do
  end do

  return
end function

function get_vector_force(sx,vx,mchi,lam,nu,gg) result (vf)
  implicit none
  real(DP),intent(in) :: sx(1:NSITE,1:NSITE),vx(1:NSITE,1:NSITE,1:2),mchi,lam,nu,gg
  real(DP) :: div_v,vf(1:NSITE,1:NSITE,1:2)
  integer :: i,j,k,pni,pnj

  do k = 1, 2
  do j = 1, NSITE
  do i = 1, NSITE
    pni = pn(i + mu(k,1))
    pnj = pn(j + mu(k,2))
    div_v = div_vec(vx,pni,pnj)

    vf(i,j,k) = - vx(i,j,k) + 0.5_DP*(div_v - gg*sx(pni,pnj))/mchi**2

    pni = pn(i - mu(k,1))
    pnj = pn(j - mu(k,2))
    div_v = div_vec(vx,pni,pnj)

    vf(i,j,k) = + vf(i,j,k) - 0.5_DP*(div_v - gg*sx(pni,pnj))/mchi**2

  end do
  end do
  end do

  return
end function
end module

program ICDW_HMC
  use mod_ICDW_HMC
  use mt19937
  implicit none
  real(DP) :: sx0(1:NSITE,1:NSITE),sx1(1:NSITE,1:NSITE)
  real(DP) :: sp0(1:NSITE,1:NSITE),sp1(1:NSITE,1:NSITE)
  real(DP) :: vx0(1:NSITE,1:NSITE,1:2),vx1(1:NSITE,1:NSITE,1:2)
  real(DP) :: vp0(1:NSITE,1:NSITE,1:2),vp1(1:NSITE,1:NSITE,1:2)
  real(DP) :: mchi,lam,nu,gg

  real(DP) :: sx2(1:NSITE,1:NSITE),sp2(1:NSITE,1:NSITE)
  real(DP) :: vx2(1:NSITE,1:NSITE,1:2),vp2(1:NSITE,1:NSITE,1:2)

  integer :: NTHERM,NSKIP,NSAMPLE,seed

  real(DP) :: tau,pacc,rand_num
  integer :: nmd,iacc,itry

  integer :: istep,iout,i,j
  real(DP) :: rho,h0,h1,h2

  iout = 99
  open(iout,file="HMC_PARAM",status='old',form='formatted',action='read')
  read(iout,*) mchi,gg,lam,nu
  read(iout,*) NTHERM,NSKIP,NSAMPLE
  read(iout,*) nmd,seed
  close(iout)

  tau = 1.0_DP
  iacc = 0
  itry = 0

  call sgrnd(seed)
  call set_unit_and_periodic (mu,pn)

  sx0(:,:) = 0.0_DP
  vx0(:,:,:) = 0.0_DP
  call generate_state_with_gaussian(sx0,sx0,1.0_DP)
  call generate_state_with_gaussian(vx0(:,:,1),vx0(:,:,1),1.0_DP)
  call generate_state_with_gaussian(vx0(:,:,2),vx0(:,:,2),1.0_DP)

  do istep = 1, NTHERM + NSKIP*NSAMPLE
    !
    ! Generate initial momentum
    !
    call generate_momentum_with_gaussian(sp0)
    call generate_momentum_with_gaussian(vp0(:,:,1))
    call generate_momentum_with_gaussian(vp0(:,:,2))

    !
    ! Compute initial Hamiltonian
    !
    h0 = hamil(sx0,sp0,vx0,vp0,mchi,lam,nu,gg)
    write(*,'(ES24.15)') h0

    sx1(:,:) = sx0(:,:)
    sp1(:,:) = sp0(:,:)
    vx1(:,:,:) = vx0(:,:,:)
    vp1(:,:,:) = vp0(:,:,:)

    call leapfrog_xpx_md(tau,nmd,sx1,sp1,vx1,vp1,mchi,lam,nu,gg)

    sp1(:,:) = -sp1(:,:)
    vp1(:,:,:) = -vp1(:,:,:)

#ifdef _CHECK_REVERSE
    sx2(:,:) = sx1(:,:)
    sp2(:,:) = sp1(:,:)
    vx2(:,:,:) = vx1(:,:,:)
    vp2(:,:,:) = vp1(:,:,:)

    call leapfrog_xpx_md(tau,nmd,sx2,sp2,vx2,vp2,mchi,lam,nu,gg)
    sp2(:,:) = -sp2(:,:)
    vp2(:,:,:) = -vp2(:,:,:)
    h2 = hamil(sx2,sp2,vx2,vp2,mchi,lam,nu,gg)

    write(*,'("@",I10,6ES24.15)') &
      & istep,sp0(3,3)-sp2(3,3),sx0(16,3)-sx2(16,3), &
      & vp0(3,10,1)-vp2(3,10,1),vp0(1,9,2)-vp2(1,9,2), &
      & vx0(8,2,1)-vx2(8,2,1),vx0(12,7,2)-vx2(12,7,2)
#endif

    h1 = hamil(sx1,sp1,vx1,vp1,mchi,lam,nu,gg)
!    write(*,'(ES24.15)') h1

!    write(*,'(ES24.15)') h1-h0
!    rho = min(1.0_DP,exp(h0-h1))
!    rand_num = grnd()
!
!    itry = itry + 1
!    if (rand_num <= rho) then
!      iacc = iacc + 1
!      continue
!    else
!      sx1(:,:) = sx0(:,:)
!      vx1(:,:,:) = vx0(:,:,:)
!    end if
!
!    if (istep > NTHERM .and. mod(istep,NSKIP) == 0) then
!      write(*,'(I10,2ES24.15)') istep,sx1(8,8)
!    end if
!
!    sx0(:,:) = sx1(:,:)
!    vx0(:,:,:) = vx1(:,:,:)
!
  end do
!
!  pacc = real(iacc,kind=DP)/itry
!  write(*,'("# HMC Metropolis test statistics.")')
!  write(*,'("# dt= ",ES14.6," itry =",I10," iacc =",I10," Pacc =",F10.6)') &
!      & tau/nmd,itry,iacc,pacc*100
!  write(*,'("# mchi =",ES24.15," lam =",ES24.15)') mchi,lam
!  write(*,'("#  nu  =",ES24.15,"  g  =",ES24.15)') nu,gg
!
!  do i = 1, NSITE
!  do j = 1, NSITE
!    write(*,'(2I10,ES24.15)') i,j,sx1(i,j)
!  end do
!    write(*,*)
!  end do

  stop
end program
