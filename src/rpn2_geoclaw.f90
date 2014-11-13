! ==============================================================================
!  Solve the normal Riemann problem for the 2D shallow water equations with
!   topography:
!       h_t + (hu)_x + (hv)_y = 0                     
!       (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x
!       (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y
!
! On input, ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!     or the y-direction if ixy=2.
!
!  Note that the i'th Riemann problem has left state qr(i-1,:)
!     and right state ql(i,:)
!  From the basic clawpack routines, this routine is called with
!     ql = qr
!
!  In order to support PyClaw use the values used from the modules are put into
!  the common block cparam defined as:
!
!    integer :: mcapa
!    real(kind=8) :: g, dry_tolerance, earth_radius, deg2rad
!    common /cparam/ mcapa, g, dry_tolerance, earth_radius, deg2rad
!
! ==============================================================================
subroutine rpn2(ixy, maxm, num_eqn, num_waves, num_aux, num_ghost, num_cells,  &
                ql, qr, auxl, auxr, fwave, s, amdq, apdq)

    implicit none

    ! Input
    integer, intent(in) :: ixy, maxm, num_eqn, num_waves, num_aux, num_ghost
    integer, intent(in) :: num_cells
    real(kind=8), intent(in out) :: ql(num_eqn, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in out) :: qr(num_eqn, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in out) :: auxl(num_aux, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in out) :: auxr(num_aux, 1 - num_ghost:maxm + num_ghost)

    ! Output
    real(kind=8), intent(in out) :: fwave(num_eqn, num_waves, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in out) :: s(num_waves, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in out) :: amdq(num_eqn, 1 - num_ghost:maxm + num_ghost)
    real(kind=8), intent(in out) :: apdq(num_eqn, 1 - num_ghost:maxm + num_ghost)

    ! Locals
    integer :: i, n, normal_index, transverse_index 
    real(kind=8) :: q_l(2), q_r(2), v_l, v_r, fwave_1d(num_eqn, num_waves)
    real(kind=8) :: s_1d(num_waves), amdq_1d(num_eqn), apdq_1d(num_eqn)

    ! Common block support - This is for compatibility with PyClaw
    integer :: mcapa
    real(kind=8) :: g, dry_tolerance, earth_radius, deg2rad
    common /cparam/ g, dry_tolerance, earth_radius, deg2rad, mcapa

    ! Initialize output
    fwave = 0.d0
    s = 0.d0
    amdq = 0.d0
    apdq = 0.d0

    ! Set direction indices
    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    end if

    ! Loop over grid cell edges
    do i = 2 - num_ghost, num_cells + num_ghost

        ! Skip problem if both states are dry
        if (qr(1, i - 1) <= dry_tolerance .and. ql(1, i) <= dry_tolerance) then
            cycle
        end if

        ! Extract transverse momentum state
        if (ql(1, i) > dry_tolerance) then
            v_r = ql(transverse_index, i) / ql(1, i)
        else
            v_r = 0.d0
        end if
        if (qr(1, i - 1) > dry_tolerance) then
            v_l = qr(transverse_index, i - 1) / qr(1, i - 1)
        else
            v_l = 0.d0
        end if

        ! Call 1D point-wise Riemann solver
        q_r = [ql(1, i), ql(normal_index, i)]
        q_l = [qr(1, i - 1), qr(normal_index, i - 1)]
        call rp1_ptwise(2, num_aux, 2, q_r, q_l, auxr(:, i - 1), auxl(:, i),   &
                        fwave_1d, s_1d, amdq_1d, apdq_1d)

        ! Unpack 1d solver results and handle transverse wave
        ! Note that fwave_1d(1, :) are the values of beta, the wave strength
        fwave(1, 1, i) = fwave_1d(1, 1)
        fwave(normal_index, 1, i) = fwave_1d(2, 1)
        fwave(transverse_index, 1, i) = s_1d(1)**2 * fwave_1d(1, 1)
        s(1, i) = s_1d(1)

        fwave(1, 3, i) = fwave_1d(1, 2)
        fwave(normal_index, 1, i) = fwave_1d(2, 2)
        fwave(transverse_index, 1, i) = s_1d(2)**2 * fwave_1d(1, 2)
        s(3, i) = s_1d(2)

        fwave(1, 2, i) = 0.d0
        fwave(normal_index, 2, i) = 0.d0
        fwave(transverse_index, 2  , i) = ql(normal_index, i) * v_r - qr(normal_index, i - 1) * v_l    &
                                            - fwave(transverse_index, 1, i)    &
                                            - fwave(transverse_index, 3, i)
        s(2, i) = 0.5d0 * (s(1, i) + s(3, i))

        ! Put fluctuations in the right spots
        amdq(1, i) = amdq_1d(1)
        apdq(1, i) = apdq_1d(1)
        amdq(normal_index, i) = amdq_1d(2)
        apdq(normal_index, i) = apdq_1d(2)

        ! Accumulate transverse values
        if (s(2, i) < 0.d0) then
            amdq(transverse_index, i) = amdq(transverse_index, i)              &
                                                 + fwave(transverse_index, 2, i)
        else if (s(2, i) > 0.d0) then
            apdq(transverse_index, i) = apdq(transverse_index, i)              &
                                                 + fwave(transverse_index, 2, i)
        else
            amdq(transverse_index, i) = amdq(transverse_index, i)              &
                                         + 0.5d0 * fwave(transverse_index, 2, i)
            apdq(transverse_index, i) = apdq(transverse_index, i)              &
                                         + 0.5d0 * fwave(transverse_index, 2, i)
        end if

        ! Handle lat-long coordinate transformation
        if (mcapa > 0) then
            if (ixy == 1) then
                s(:, i) = earth_radius * deg2rad * s(:, i)
                do n = 1, num_eqn
                    fwave(n, :, i) = earth_radius * deg2rad * fwave(n, :, i)
                end do
                amdq(:, i) = earth_radius * deg2rad * amdq(:, i)
                apdq(:, i) = earth_radius * deg2rad * apdq(:, i)
            else
                s(:, i) = earth_radius * deg2rad * auxl(3, i) * s(:, i)
                do n = 1, num_eqn
                    fwave(n, :, i) = earth_radius * deg2rad * auxl(3, i) * fwave(n, :, i)
                end do
                amdq(:, i) = earth_radius * deg2rad * auxl(3, i) * amdq(:, i)
                apdq(:, i) = earth_radius * deg2rad * auxl(3, i) * apdq(:, i)
            end if
        end if

    end do

end subroutine rpn2
