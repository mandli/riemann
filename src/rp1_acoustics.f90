! =====================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Riemann solver for the acoustics equations in 1D.

! waves:     2
! equations: 2

! Conserved quantities:
!       1 pressure
!       2 velocity

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! On output, wave contains the waves,
!            s the speeds,
! 
!            amdq = A^- Delta q,
!            apdq = A^+ Delta q,
!                   the decomposition of the flux difference
!                       f(qr(i-1)) - f(ql(i))
!                   into leftgoing and rightgoing parts respectively.
! 

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    ! Input/Output
    integer, intent(in) :: maxm, meqn, mwaves, maux, mbc, mx
    real(kind=8), intent(inout) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: s(mwaves,1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: ql(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: qr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: auxl(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: auxr(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: apdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: amdq(meqn, 1-mbc:maxm+mbc)

    ! Local Storage
    integer :: i
    real(kind=8) :: delta(2), alpha(2)

    ! Common block Storage
    real(kind=8) :: rho, bulk, cc, zz

    ! density, bulk modulus, and sound speed, and impedence of medium:
    ! (should be set in setprob.f)
    common /cparam/ rho, bulk, cc, zz


    ! split the jump in q at each interface into waves
    ! find alpha(1) and alpha(2), the coefficients of the 2 eigenvectors:
    do i = 2-mbc, mx+mbc
        delta = ql(:, i) - qr(:, i - 1)

        alpha = [(-delta(1) + zz * delta(2)) / (2.0d0 * zz), &
                 ( delta(1) + zz * delta(2)) / (2.d0 * zz) ]
    
        ! Compute the waves - alpha * eigenvector
        wave(:, 1, i) = alpha(1) * [-zz, 1.d0]
        s(1, i) = -cc

        wave(:, 2, i) = alpha(2) * [zz, 1.d0]
        s(2, i) = cc
    end do

    ! Compute the left and right going fluctuations
    ! Note that in this case s(1) <= 0 and s(2) >= 0 always
    do i = 2 - mbc, mx + mbc
        amdq(:, i) = s(1, i) * wave(:, 1, i)
        apdq(:, i) = s(2, i) * wave(:, 2, i)
    end do

end subroutine rp1
