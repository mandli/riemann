subroutine rpn2(ixy, maxm, meqn, mwaves, maux, mbc, mx, ql, qr, auxl, auxr,    &
                fwave, s, amdq, apdq)
    !======================================================================
    !
    ! Solves normal Riemann problems for the 2D SHALLOW WATER equations
    !     with topography:
    !     #        h_t + (hu)_x + (hv)_y = 0                           #
    !     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
    !     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #
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
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                          !
    !      # This Riemann solver is for the shallow water equations.           !
    !                                                                          !
    !       It allows the user to easily select a Riemann solver in            !
    !       riemannsolvers_geo.f. this routine initializes all the variables   !
    !       for the shallow water equations, accounting for wet dry boundary   !
    !       dry cells, wave speeds etc.                                        !
    !                                                                          !
    !           David George, Vancouver WA, Feb. 2009                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use geoclaw_module, only: g => grav, drytol => dry_tolerance
    use geoclaw_module, only: earth_radius, deg2rad
    use amr_module, only: mcapa

    implicit none

    ! Input arguments
    integer, intent(in) :: ixy,maxm,meqn,mwaves,mbc,mx,maux
    real(kind=8), dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    real(kind=8), dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr

    ! Output arguments
    real(kind=8), dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: fwave
    real(kind=8), dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    real(kind=8), dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: apdq, amdq

    ! Counters
    integer :: m, i, mw

    ! Indexes
    integer :: normal_index, transverse_index

    ! Other
    real(kind=8) :: wall_mask(3)
    real(kind=8) :: fwave_local(3, 3), s_local(3)

    integer, parameter :: MAX_ITER

    ! State
    real(kind=8) :: h_R, h_L, hu_R, hu_L, u_R, u_L, hv_R, hv_L, v_R, v_L
    real(kind=8) :: phi_R, phi_L, b_R, b_L
    real(kind=8) :: h_star, h_star_test, h_star_HLL, s_L_test, s_R_test

    ! Riemann problem locals
    real(kind=8) :: s_L, s_R, s_Roe(2), s_einfeldt(2), u_hat, c_hat
    logical :: rare(2)

    real(kind=8) :: tw, dxdc
    double precision s1m,s2m

    ! Initialize output
    amdq = 0.d0
    apdq = 0.d0
    fwave = 0.d0
    s = 0.d0
    
    ! Set normal direction
    if (ixy == 1) then
        n_index = 2
        t_index = 3
    else
        n_index = 3
        t_index = 2
    endif

    ! Inform of bad Riemann problem
    ! ... need to do an any here...
    !-----------------------Initializing-----------------------------------
         ! !inform of a bad riemann problem from the start
         ! if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
         !    write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
         ! endif
         !          !zero (small) negative values if they exist
         ! if (qr(1,i-1).lt.0.d0) then
         !       qr(1,i-1)=0.d0
         !       qr(2,i-1)=0.d0
         !       qr(3,i-1)=0.d0
         ! endif

         ! if (ql(1,i).lt.0.d0) then
         !       ql(1,i)=0.d0
         !       ql(2,i)=0.d0
         !       ql(3,i)=0.d0
         ! endif

    ! ========================================================================
    ! Loop through Riemann problems
    ! ========================================================================
    do i=2-mbc,mx+mbc

        ! Skip problem if in a completely dry area
        if (qr(1, i - 1) <= dry_tolerance .and. ql(1, i) <= dry_tolerance) then
            cycle
        end if

        ! Extract states for this Reimann problem
        h_L = qr(1, i - 1)
        h_R = ql(1, i)
        hu_L = qr(n_index, i - 1)
        hu_R = ql(n_index, i)
        hv_L = qr(t_index, i - 1)
        hv_R = ql(t_index, i)

        b_L = auxr(1, i - 1)
        b_R = auxl(1, i)

        ! Check for wet/dry boundaries
        if (h_R > dry_tolerance) then
            u_R = hu_R / h_R
            v_R = hv_R / h_R
            phi_R = 0.5d0 * g * h_R**2 + hu_R**2 / h_R
        else
            hu_R = 0.d0
            hv_R = 0.d0
            u_R = 0.d0
            v_R = 0.d0
            phi_R = 0.d0
        end if

        if (h_L > dry_tolerance) then
            u_L = hu_L / h_L
            v_L = hv_L / h_L
            phi_L = 0.5d0 * g * h_L**2 + hu_L**2 / h_L
        else
            hu_L = 0.d0
            hv_L = 0.d0
            u_L = 0.d0
            v_L = 0.d0
            phi_L = 0.d0
        end if

        ! Initialize wall mask
        wall_mask = 1.d0

        ! If we are near a wet/dry boundary attempt to ascertain what kind of 
        ! problem we need to solve
        if (h_R <= dry_tolerance) then
            call riemanntype(h_L, h_L, u_L, -u_L, h_star, s_m(1), s_m(2),      &
                             rare(1), rare(2), 1, dry_tolerance, g)
            h_star_test = max(h_L, h_star)

            ! If true then right state should become ghost values so that we can
            ! solve a wall like boundary condition
            if (h_star_test + b_L < b_R) then
                b_R = h_star_test + b_L
                wall_mask(2:3) = 0.d0
                h_R = h_L
                hu_R = -hu_L
                b_R = b_L
                phi_R = phi_L
                u_R = -u_L
                v_R = v_L
            ! Otherwise we need to solve the full inundation problem
            else
                b_L = h_R + b_R
            end if

        else if (h_L <= dry_tolerance) then
            call riemanntype(h_R, h_R, -u_R, u_R, h_star, s_m(1), s_m(2),      &
                             rare(1), rare(2), 1, dry_tolerance, g)
            h_star_test = max(h_R, h_star)

            ! If true then left state should become ghost values so that we can
            ! solve a wall like boundary condition
            if (h_star_test + b_R < b_L) then
                b_L = h_star_test + b_R
                wall_mask(1:2) = 0.d0
                h_L = h_R
                hu_L = -hu_R
                b_L = b_R
                phi_L = phi_R
                u_L = -u_R
                v_L = v_R
            ! Otherwise we need to solve the full inundation problem
            else
                b_L = h_R + b_R
            end if
        end if

        ! ======================================================================
        !  Determine relevant speeds
        ! Left and right state speeds
        s_L = u_L - sqrt(g * h_L)
        s_R = u_R + sqrt(g * h_R)

        ! Roe-averages
        u_hat = (sqrt(g * h_L) * u_L + sqrt(g * h_R) * u_R)             &
                        / (sqrt(g * h_R) + sqrt(g * h_L))
        c_hat = sqrt(g * 0.5d0 * (h_R + h_L))
        s_Roe(1) = u_hat - c_hat
        s_Roe(2) = u_hat + c_hat

        ! Eindfeldt speeds
        s_E(1) = min(s_L, s_Roe(1))
        s_E(2) = max(s_R, s_Roe(2))

        ! Call specific Riemann solver
        call riemann_solver_func(MAX_ITER, 3, 3, h_L, h_R, hu_L, hu_R,        &
                                                 hv_L, hv_R, b_l, b_R,        &
                                                 u_L, u_R, v_L, v_R,          &
                                                 phi_L, phi_R,                &
                                                 s_E(1), s_E(2),              &
                                                 dry_tolerance, g,            &
                                                 s_local, fwave_local)

        ! Eliminate ghost fluxes for wall
        s_local = s_local * wall
        fwave_local(1, :) = fwave_local(1, :) * wall
        fwave_local(2, :) = fwave_local(2, :) * wall
        fwave_local(3, :) = fwave_local(3, :) * wall

        ! Store new waves and speeds in actuall output
        s(:, i) = s_local
        fwave(1, :, i) = fwave_local(1, :)
        fwave(normal_index, :, i) = fwave_local(normal_index, :)
        fwave(transverse_index, :, i) = fwave_local(transverse_index, :)

    end do
    ! == End of Riemann Solver Loop per grid cell ==============================

    ! ==========================================================================
    ! Handle capacity mapping from lat-long to physical space
    if (mcapa > 0) then
        if (ixy == 1) then
            dxdc = earth_radius * deg2rad
            do i = 2 - mbc, mx + mbc
                s(:, i) = dxdc * s(:, i)
                fwave(1, :, i) = dxdc * fwave(1, :, i)
                fwave(2, :, i) = dxdc * fwave(2, :, i)
                fwave(3, :, i) = dxdc * fwave(3, :, i)
            end do
        else
            do i = 2 - mbc, mx + mbc
                dxdc = earth_radius * cos(auxl(3, i)) * deg2rad
                s(:, i) = dxdc * s(:, i)
                fwave(1, :, i) = dxdc * fwave(1, :, i)
                fwave(2, :, i) = dxdc * fwave(2, :, i)
                fwave(3, :, i) = dxdc * fwave(3, :, i)
            end do
        end if
    end if

    ! ==========================================================================
    ! Compute fluctuations
    do i = 2 - mbc, mx + mbc
        do mw = 1, mwaves
            if (s(mw, i) < 0.d0) then
                amdq(1:3, i) = amdq(1:3, i) + fwave(1:3, mw, i)
            else if (s(mw, i) > 0.d0) then
                apdq(1:3, i) = apdq(1:3, i) + fwave(1:3, mw, i)
            else
                amdq(1:3, i) = amdq(1:3, i) + 0.5d0 * fwave(1:3, mw, i)
                apdq(1:3, i) = apdq(1:3, i) + 0.5d0 * fwave(1:3, mw, i)
            end if
        end do
    end do

end subroutine rpn2
