! ==============================================================================
!  Fwave based simple Riemann solver for 2D shallow water equations
!   - Handles bathymetry
! ==============================================================================
subroutine rpn2(ixy, maxm, meqn, mwaves, maux, mbc, mx, ql, qr, auxl, auxr,    &
                fwave, s, amdq, apdq)

    use geoclaw_module, only: g => grav, dry_tolerance
    ! use geoclaw_module, only: earth_radius, deg2rad
    ! use amr_module, only: mcapa

    implicit none

    ! Input/Output
    integer, intent(in) :: maxm, meqn, maux, mwaves, mbc, mx, ixy

    real(kind=8), intent(inout) :: fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: s(mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: ql(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: qr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: apdq(meqn,1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: amdq(meqn,1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: auxl(maux,1-mbc:maxm+mbc)
    real(kind=8), intent(inout) :: auxr(maux,1-mbc:maxm+mbc)

    ! Locals
    integer :: normal_index, transverse_index, i, mw
    logical :: dry_state_l, dry_state_r
    real(kind=8) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r, b_r, b_l
    real(kind=8) :: u_r, u_l, v_r, v_l, phi_r, phi_l, delta(3), beta(3)

    ! Initialize return values
    fwave = 0.d0
    s = 0.d0
    amdq = 0.d0
    apdq = 0.d0

    ! Set normal direction
    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    end if

    ! Primary grid-cell edge loop
    do i = 2 - mbc, mx + mbc

        ! Extract relevant states
        h_l = qr(1, i-1)
        h_r = ql(1, i)
        hu_l = qr(normal_index, i-1)
        hu_r = ql(normal_index, i)
        hv_l = qr(transverse_index, i-1)
        hv_r = ql(transverse_index, i)
        b_l = auxr(1, i-1)
        b_r = auxl(1, i)

        dry_state_r = .false.
        dry_state_l = .false.

        if (h_r > dry_tolerance) then
            u_r = hu_r / h_r
            v_r = hv_r / h_r
            phi_r = 0.5d0 * g * h_r**2 + hu_r**2 / h_r
        else
            dry_state_r = .true.
            h_r = 0.d0
            hu_r = 0.d0
            hv_r = 0.d0
            u_r = 0.d0
            v_r = 0.d0
            phi_r = 0.d0
        end if

        if (h_l > dry_tolerance) then
            u_l = hu_l / h_l
            v_l = hv_l / h_l
            phi_l = 0.5d0 * g * h_l**2 + hu_l**2 / h_l
        else
            dry_state_l = .true.
            h_l = 0.d0
            hu_l = 0.d0
            hv_l = 0.d0
            u_l = 0.d0
            v_l = 0.d0
            phi_l = 0.d0
        end if

        ! Skip this cell if both sides are dry
        if (dry_state_r .and. dry_state_l) then
            cycle
        end if

        ! Calculate speeds
        s(1, i) = u_l - sqrt(g * h_l)
        s(3, i) = u_r + sqrt(g * h_r)
        s(2, i) = 0.5d0 * (s(1, i) + s(3, i))

        ! Calculate flux differences, also take into account bathy jump
        delta(1) = hu_r - hu_l
        delta(2) = h_r * u_r * v_r - h_l * u_l * v_l
        delta(3) = phi_r - phi_l + g * 0.5d0 * (h_r + h_l) * (b_r - b_l)

        beta(1) = (s(3, i) * delta(1) - delta(3)) / (s(3, i) - s(1, i))
        beta(2) = (delta(3) - s(1, i) * delta(1)) / (s(3, i) - s(1, i))

        fwave(1, 1, i) = beta(1)
        fwave(normal_index, 1, i) = beta(1) * s(1, i)
        fwave(transverse_index, 1, i) = beta(1) * v_l

        fwave(1, 3, i) = beta(2)
        fwave(normal_index, 3, i) = beta(2) * s(3, i)
        fwave(transverse_index, 3, i) = beta(2) * v_r

        fwave(1, 2, i) = 0.d0
        fwave(normal_index, 2, i) = 0.d0
        fwave(transverse_index, 2, i) = delta(2) - fwave(transverse_index, 1, i) - fwave(transverse_index, 3, i)

    end do
    ! End of grid-cell edge loop

    ! Compute fluctuations
    do i = 2 - mbc, mx + mbc
        do mw = 1, mwaves
            if (s(mw, i) < 0.d0) then
                amdq(:, i) = amdq(:, i) + fwave(:, mw, i)
            else if (s(mw, i) > 0.d0) then
                apdq(:, i) = apdq(:, i) + fwave(:, mw, i)
            else
                amdq(:, i) = amdq(:, i) + 0.5d0 * fwave(:, mw, i)
                apdq(:, i) = apdq(:, i) + 0.5d0 * fwave(:, mw, i)
            end if
        end do
    end do

end subroutine rpn2
