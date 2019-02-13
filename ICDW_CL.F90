module mod_ICDW_2d_CL
  implicit none
  integer,parameter :: DP=kind(1.0d0)
  real(DP),parameter :: PI=2.0_DP*acos(0.0_DP)
  complex(DP),parameter :: zi=(0.0_DP,1.0_DP)
  integer,parameter :: NSITE=16
  integer :: mu(1:2,1:2),pn(0:NSITE+1)

contains

subroutine set_unit_and_periodic (mu,pn)
  implicit none
  integer,intent(inout) :: mu(1:2,1:2),pn(0:NSITE+1)
  integer :: i

  mu(1,1) = 1
  mu(1,2) = 0
  mu(2,1) = 0
  mu(2,2) = 1

  do i = 1, NSITE
    pn(i) = i
  end do
  pn(0) = NSITE
  pn(NSITE+1) = 1

  return
end subroutine

function gauss_rand (D,dt) result (rand_num)
  !
  ! Generate a Gaussiand random number via Box-Muller method
  ! mean = 0, variance = sqrt(2*D*dt)
  !
  use mt19937
  implicit none
  real(DP),intent(in) :: D,dt
  real(DP) :: rr,tt,rand_num

  rr = grnd()
  tt = grnd()
  rr = 1.0_DP - rr
  tt = 1.0_DP - tt
  rand_num = sqrt(-2.0_DP*log(rr))*cos(2.0_DP*PI*tt)
  rand_num = rand_num*sqrt(2.0_DP*D*dt)

  return
end function

function mu_sum_sca (s,i,j) result (sum_s)
  implicit none
  complex(DP),intent(in) :: s(1:NSITE,1:NSITE)
  integer,intent(in) :: i,j
  complex(DP) :: sum_s
  integer :: k

  sum_s = (0.0_DP,0.0_DP)
  do k = 1, 2
    sum_s = sum_s + s(pn(i + mu(k,1)),pn(j + mu(k,2))) + s(pn(i - mu(k,1)),pn(j - mu(k,2)))
  end do

  return
end function

subroutine update_config (phiz0,chiz0,phiz1,chiz1,mchi,lam,nu,gg,dt)
  implicit none
  complex(DP),intent(in) :: phiz0(1:NSITE,1:NSITE),chiz0(1:NSITE,1:NSITE)
  complex(DP),intent(inout) :: phiz1(1:NSITE,1:NSITE),chiz1(1:NSITE,1:NSITE)
  real(DP),intent(in) :: mchi,lam,nu,gg,dt
  complex(DP) :: dvph,dvch
  real(DP) :: wi
  integer :: ix,jy,ixp,jyp,ixm,jym,vec

  !
  ! Update configration of scalar field
  !
  do jy = 1, NSITE
  do ix = 1, NSITE
    dvph = (0.0_DP,0.0_DP)
    dvch = (0.0_DP,0.0_DP)
  do vec = 1, 2
    ixp = pn(ix + mu(vec,1))
    jyp = pn(jy + mu(vec,2))
    ixm = pn(ix - mu(vec,1))
    jym = pn(jy - mu(vec,2))

    dvph = dvph + phiz0(ixp,jyp) - 2.0_DP*phiz0(ix,jy) + phiz0(ixm,jym)
    dvch = dvch + chiz0(ixp,jyp) - 2.0_DP*chiz0(ix,jy) + chiz0(ixm,jym)
  end do
    wi = gauss_rand(1.0_DP,dt)
    phiz1(ix,jy) = phiz0(ix,jy) &
      & + (dvph - 4.0_DP*lam*phiz0(ix,jy)*(phiz0(ix,jy)**2-nu**2) + zi*gg*chiz0(ix,jy))*dt + wi

    wi = gauss_rand(1.0_DP,dt)
    chiz1(ix,jy) = chiz0(ix,jy) &
              & + (dvch - mchi**2*chiz0(ix,jy) + zi*gg*phiz0(ix,jy))*dt + wi
  end do
  end do

  return
end subroutine
end module

program ICDW_2d
  use mod_ICDW_2d_CL
  use mt19937
  implicit none
  complex(DP) :: phiz0(1:NSITE,1:NSITE),phiz1(1:NSITE,1:NSITE)
  complex(DP) :: chiz0(1:NSITE,1:NSITE),chiz1(1:NSITE,1:NSITE)

  real(DP) :: dt,mchi,lam,nu,gg

  integer :: NTHERM,NSKIP,NSAMPLE,seed
  integer :: tstep,ndo,iout,i,j,k

  !
  ! Read input parameter
  !
  iout = 99
  open(iout,file="CL_PARAM",status='old',form='formatted',action='read')
  read(iout,*) mchi,gg,lam,nu
  read(iout,*) NTHERM,NSKIP,NSAMPLE
  read(iout,*) dt,seed
  close(iout)

  !
  ! Set up of seed of pseudo-number generator
  !
  call sgrnd(seed)

  !
  ! Set up of unit two dimension vector "mu" and pn
  !
  call set_unit_and_periodic (mu,pn)

  !
  ! Set up of initial distribution
  !
  do j = 1, NSITE
  do i = 1, NSITE
    phiz0(i,j) = cmplx(grnd()-0.5_DP,0.0_DP)
    chiz0(i,j) = cmplx(grnd()-0.5_DP,0.0_DP)
  end do
  end do

  !
  ! complex Langevin dynamics part
  !
  do tstep = 1, NTHERM
    call update_config (phiz0,chiz0,phiz1,chiz1,mchi,lam,nu,gg,dt)
  end do

  do tstep = 1, NSAMPLE
  do ndo = 1, NSKIP
    call update_config (phiz0,chiz0,phiz1,chiz1,mchi,lam,nu,gg,dt)
  end do
    write(*,'(I10,6ES24.15)') tstep,phiz1(1,2),chiz1(3,6)
  end do

  write(*,'("# mchi =",ES24.15," lam =",ES24.15)') mchi,lam
  write(*,'("#  nu  =",ES24.15,"  g  =",ES24.15)') nu,gg

!  do j = 1, NSITE
!  do i = 1, NSITE
!    write(*,'(2I5,2ES24.15)') i,j,phiz0(i,j)
!  end do
!    write(*,*)
!  end do

  stop
end program
