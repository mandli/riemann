subroutine rpt2(ixy, imp, maxm, meqn, mwaves, maux, mbc, mx, ql, qr,           &
                aux1, aux2, aux3, asdq, bmasdq, bpasdq)

    use geoclaw_module, only: g => grav, dry_tolerance
!     use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

    implicit none

    ! Arguments
    integer, intent(in) :: ixy, maxm, meqn, maux, mwaves, mbc, mx, imp

    real(kind=8), intent(in) :: ql(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: qr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: asdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: aux1(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: aux2(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: aux3(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: bmasdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: bpasdq(meqn, 1-mbc:maxm+mbc)

    ! Locals
    integer :: i, normal_index, transverse_index, mw
    real(kind=8) :: h, hu, hv, u, v, b(3)
    real(kind=8) :: s(3), delta(3), beta(3), r(3,3)

    ! Initialize return values
    bmasdq = 0.d0
    bpasdq = 0.d0

    ! Set normal direction
    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    end if

    ! Primarly loop over cell edges
    do i = 2 - mbc, mx + mbc

        ! Skip this edge if both sides are dry
        if (qr(1, i - 1) < dry_tolerance .and. ql(1, i) < dry_tolerance) then
            cycle
        end if

        ! Extract states from correct side of edge
        if (imp == 1) then
            ! Left splitting
            h = qr(1, i - 1)
            hu = qr(normal_index, i - 1)
            hv = qr(transverse_index, i - 1)
            b(3) = aux3(1, i - 1)
            b(2) = aux2(1, i - 1)
            b(1) = aux1(1, i - 1)
        else
            ! Right splitting
            h = ql(1, i)
            hu = ql(normal_index, i)
            hv = ql(transverse_index, i)
            b(3) = aux3(1, i)
            b(2) = aux2(1, i)
            b(1) = aux1(1, i)
        end if

        if (h < dry_tolerance) then
            u = 0.d0
            v = 0.d0
        else
            u = hu / h
            v = hv / h
        end if

        ! If both transverse cells are higher than eta, skip this splitting
        if (h + b(2) < b(1) .and. h + b(2) < b(3)) then
            cycle
        end if

        ! Determine wave speeds
        s(1) = v - sqrt(g * h)
        s(2) = u
        s(3) = v + sqrt(g * h)

        ! Calculate asdq decomposition
        delta(1) = asdq(1, i)
        delta(2) = asdq(normal_index, i)
        delta(3) = asdq(transverse_index, i)
        beta(1) = s(3) * delta(1) / (s(3) - s(1)) - delta(3) / (s(3) - s(1))
        beta(2) = -s(2) * delta(1) + delta(2)
        beta(3) = delta(3) / (s(3) - s(1)) - s(1) * delta(1) / (s(3) - s(1))

        ! Eigenvectors
        r(1, :) = [1.d0, 0.d0, 1.d0]
        r(2, :) = [s(2), 1.d0, s(2)]
        r(3, :) = [s(1), 0.d0, s(3)]

        beta(3) = s(1) * delta(1) / (delta(3) - s(3) + s(1))
        beta(1) = (delta(3) - s(3) * beta(3)) / s(1)
        beta(2) = -(s(2) * beta(1) + s(2) * beta(3) - delta(2))

        ! Computate fluctuations
        do mw = 1, mwaves
            if (s(mw) < 0.d0 .and. h + b(2) >= b(1)) then
                bmasdq(1, i) = bmasdq(1, i) + s(mw) * beta(mw) * r(1, mw)
                bmasdq(normal_index, i) = bmasdq(normal_index, i)           &
                                        + s(mw) * beta(mw) * r(2, mw)
                bmasdq(transverse_index, i) = bmasdq(transverse_index, i)   &
                                    + s(mw) * beta(mw) * r(3, mw)
            else if (s(mw) > 0.d0 .and. h + b(2) >= b(3)) then
                bpasdq(1, i) = bpasdq(1, i) + s(mw) * beta(mw) * r(1, mw)
                bpasdq(normal_index, i) = bpasdq(normal_index, i)           &
                                        + s(mw) * beta(mw) * r(2, mw)
                bpasdq(transverse_index, i) = bpasdq(transverse_index, i)   &
                                    + s(mw) * beta(mw) * r(3, mw)
!             else
!                 bmasdq(1, i) = bmasdq(1, i) + 0.5d0 * s(mw) * beta(mw) * r(1, mw)
!                 bmasdq(normal_index, i) = bmasdq(normal_index, i)           &
!                                         + 0.5d0 * s(mw) * beta(mw) * r(normal_index, mw)
!                 bmasdq(transverse_index, i) = bmasdq(transverse_index, i)   &
!                                     + 0.5d0 * s(mw) * beta(mw) * r(transverse_index, mw)
!                 bpasdq(1, i) = bpasdq(1, i) + 0.5d0 * s(mw) * beta(mw) * r(1, mw)
!                 bpasdq(normal_index, i) = bpasdq(normal_index, i)           &
!                                         + 0.5d0 * s(mw) * beta(mw) * r(normal_index, mw)
!                 bpasdq(transverse_index, i) = bpasdq(transverse_index, i)   &
!                                     + 0.5d0 * s(mw) * beta(mw) * r(transverse_index, mw)
            end if
        end do

    end do
    ! End primary loop

end subroutine rpt2
