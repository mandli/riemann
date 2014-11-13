! ==============================================================================
!  Calculate the middle state of a ghost-cell wall problem, note that this is 
!   much more simple than the older riemanntype routine that used to be used.
! ==============================================================================
function calc_middle_state(h, u) result(h_middle)
    
    implicit none

    ! Input
    real(kind=8), intent(in) :: h, u

    ! Output
    real(kind=8) :: h_middle

    ! Locals
    integer :: n
    real(kind=8) :: h0, g_star, F0, dfdh, slope

    ! Constants
    integer, parameter :: MAX_ITERATIONS = 1

    ! Common block
    integer :: mcapa
    real(kind=8) :: g, dry_tolerance, earth_radius, deg2rad
    common /cparam/ g, dry_tolerance, earth_radius, deg2rad, mcapa

    h0 = h
    do n = 1, MAX_ITERATIONS
        g_star = sqrt(0.5d0 * g * (1.d0 / h0 + 1.d0 / h))
        F0 = 2.d0 * u + (h0 - h) * g_star + (h0 - h) * g_star
        dfdh = 2.d0 * g_star - g * (h0 - h) / (2.d0 * h0**2 * g_star)
        slope = 2.d0 * sqrt(h0) * dfdh
        h0 = (sqrt(h0) - F0 / slope)**2
    end do
    h_middle = h0

end function calc_middle_state

! ==============================================================================
!  Calculate the wave and speed structure using a simple fwave approach.
! ==============================================================================
subroutine riemann_fwave(h_l, h_r, hu_l, hu_r, phi_l, phi_r, b_l, b_r,         &
                         u_l, u_r, fwave, s)
    
    implicit none

    ! Input
    real(kind=8), intent(in) :: h_l, h_r, hu_l, hu_r, phi_l, phi_r
    real(kind=8), intent(in) :: b_l, b_r, u_l, u_r
    real(kind=8), intent(out) :: fwave(2,2), s(2)

    ! Local
    real(kind=8) :: s_l, s_r, u_hat, c_hat, delta(4), beta(2)
    real(kind=8) :: s_roe_l, s_roe_r, s_einfeldt_l, s_einfeldt_r

    real(kind=8) :: calc_middle_state

    ! Common block
    integer :: mcapa
    real(kind=8) :: g, dry_tolerance, earth_radius, deg2rad
    common /cparam/ g, dry_tolerance, earth_radius, deg2rad, mcapa

    ! Left/Right state speeds
    s_l = u_l - sqrt(g * h_l)
    s_r = u_r + sqrt(g * h_r)

    ! Roe speeds
    u_hat = (sqrt(g * h_l) * u_l + sqrt(g * h_r) * u_r)             &
                        / (sqrt(g * h_r) + sqrt(g * h_l))
    c_hat = sqrt(g * 0.5d0 * (h_r + h_l))
    s_roe_l = u_hat - c_hat
    s_roe_r = u_hat + c_hat

    ! Einfeldt speeds
    s_einfeldt_l = min(s_l, s_roe_l)
    s_einfeldt_r = max(s_r, s_roe_r)

    ! Calculate waves and speeds
    delta(1) = h_r - h_l
    delta(2) = hu_r - hu_l
    delta(3) = phi_r - phi_l
    delta(4) = -g * 0.5d0 * (h_r + h_l) * (b_r - b_l)

    beta(1) = (s_einfeldt_r * delta(2) - (delta(3) - delta(4))) / (s_einfeldt_r - s_einfeldt_l)
    beta(2) = ((delta(3) - delta(4)) - s_einfeldt_l * delta(2)) / (s_einfeldt_r - s_einfeldt_l)

    s(1) = s_einfeldt_l
    s(2) = s_einfeldt_r

    fwave(:, 1) = [1.d0, s(1)] * beta(1)
    fwave(:, 2) = [1.d0, s(2)] * beta(2)
    
end subroutine riemann_fwave

! ==============================================================================
!  
! ==============================================================================
subroutine rp1_ptwise(num_eqn, num_aux, num_waves, q_l, q_r, aux_l, aux_r,    &
                      fwave, s, amdq, apdq)

    implicit none

    ! Input Arguments
    integer, intent(in) :: num_eqn, num_aux, num_waves
    real(kind=8), intent(in out) :: q_l(num_eqn), q_r(num_eqn)
    real(kind=8), intent(in out) :: aux_l(num_aux), aux_r(num_aux)

    ! Output arguments
    real(kind=8), intent(in out) :: fwave(num_eqn, num_waves)
    real(kind=8), intent(in out) :: s(num_waves)
    real(kind=8), intent(in out) :: apdq(num_eqn), amdq(num_eqn)

    ! Locals
    integer :: mw
    real(kind=8) :: h_l, h_r, hu_l, hu_r, u_l, u_r, b_l, b_r, phi_l, phi_r
    real(kind=8) :: wall(2), fw(2, 2), sw(2), h_star

    real(kind=8) :: calc_middle_state

    ! Common block
    integer :: mcapa
    real(kind=8) :: g, dry_tolerance, earth_radius, deg2rad
    common /cparam/ g, dry_tolerance, earth_radius, deg2rad, mcapa

    fwave = 0.d0
    s = 0.d0
    amdq = 0.d0
    apdq = 0.d0

    ! Skip this Riemann problem if both states are dry
    if (q_l(1) < dry_tolerance .and. q_r(1) < dry_tolerance) then
        return
    end if

    ! Extract state
    h_l = q_l(1)
    h_r = q_r(1)
    hu_l = q_l(2)
    hu_r = q_r(2)
    b_l = aux_l(1)
    b_r = aux_r(1)

    ! Check for wet/dry boundary
    if (h_r > dry_tolerance) then
        u_r = hu_r / h_r
        phi_r = hu_r**2 / h_r**2 + 0.5d0 * g * h_r**2
    else
        h_r = 0.d0
        u_r = 0.d0
        phi_r = 0.d0
    end if
    if (h_l > dry_tolerance) then
        u_l = hu_l / h_l
        phi_l = hu_l**2 / h_l**2 + 0.5d0 * g * h_l**2
    else
        h_l = 0.d0
        u_l = 0.d0
        phi_l = 0.d0
    end if

    wall = 1.d0

    ! Check for type of Riemann problem, in particular whether to use a ghost
    ! cell wall approach if the wave does not crest the bathy jump
    if (h_l <= dry_tolerance) then
        stop
        if (u_r > 0.d0) then
            ! Need to check that middle state does not crest bathy jump
            h_star = calc_middle_state(h_r, u_r)
            if (h_star + h_l > b_r) then
                ! Inundation problem - Set right bathymetry to be equal to
                ! height of surface
                b_l = h_r + b_r
            else
                ! Set ghost cell wall problem
                h_l = h_r
                hu_l = -hu_r
                b_l = b_r
                u_l = -u_r
                phi_l = phi_r
                wall(2) = 0.d0
            end if
        else
            ! Set ghost cell wall problem
            h_l = h_r
            hu_l = -hu_r
            b_l = b_r
            u_l = -u_r
            phi_l = phi_r
                wall(2) = 0.d0
        end if

    else if (h_r <= dry_tolerance) then
        stop
        if (u_l > 0.d0) then
            ! Need to check that middle state does not crest bathy jump
            h_star = calc_middle_state(h_l, u_l)
            if (h_star + h_r > b_l) then
                ! Inundation problem - Set left bathymetry to be equal to
                ! height of surface
                b_r = h_l + b_l
            else
                ! Set ghost cell wall problem
                h_r = h_l
                hu_r = -hu_l
                b_r = b_l
                u_r = -u_l
                phi_r = phi_l
            end if
        else
            ! Set ghost cell wall problem
            h_r = h_l
            hu_r = -hu_l
            b_r = b_l
            u_r = -u_l
            phi_r = phi_l
        end if

    end if

    ! Calculate waves and speeds
    call riemann_fwave(h_l, h_r, hu_l, hu_r, phi_l, phi_r, b_l, b_r, u_l, u_r, &
                       fw, sw)


    ! Eliminate ghost fluxes for wall
    s = sw * wall
    do mw = 1, num_waves
        fwave(:, mw) = fw(:, mw) * wall
    end do

    do mw = 1, num_waves
        if (s(mw) < 0.d0) then
            amdq = amdq + fwave(:, mw)
        else if (s(mw) > 0.d0) then
            apdq = apdq + fwave(:, mw)
        else
            amdq = amdq + 0.5d0 * fwave(:, mw)
            apdq = apdq + 0.5d0 * fwave(:, mw)
        end if
    end do

end subroutine rp1_ptwise