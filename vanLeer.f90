subroutine vanleer


    ! This subroutine calculate the flluxes for van-leers fvs scheme. Note that
    ! this is not the MUSCL one my dear !

    use shared
    implicit none

    integer(kind=4) :: i, gs

    real(kind=8) :: rho, s1, s2

    gs = ng / 2.0d0


    do i = 1, total_mesh_points + ng


        ! Define some usefull properties.

        rho = q(i,1)


        ! All the matrices are multiplyed by the same coefficient, lets define.

        s1 = 1.0d0/4.0d0*rho*a_speed(i)*(1.0d0+mach(i))**2.0d0
        s2 = 1.0d0/4.0d0*rho*a_speed(i)*(1.0d0-mach(i))**2.0d0


        ! Lets calculate the positive fluxes

        f_pls(i,1) = s1 * 1.0d0
        f_pls(i,2) = s1 * ( 2.0d0*a_speed(i) ) / fgamma * ( (fgamma-1.0d0)/2.0d0*mach(i) + 1.0d0 )
        f_pls(i,3) = s1 * ( 2.0d0*a_speed(i)**2.0d0 ) / (fgamma**2.0d0-1.0d0) * ( (fgamma-1.0d0)/2.0d0*mach(i) + 1.0d0 )**2.0d0


        ! Now lets do the negative fluxes.

        f_min(i,1) = -s2 * 1.0d0
        f_min(i,2) = -s2 * ( 2.0d0*a_speed(i) ) / fgamma * ( (fgamma-1.0d0)/2.0d0*mach(i) - 1.0d0 )
        f_min(i,3) = -s2 * ( 2.0d0*a_speed(i)**2.0d0 ) / (fgamma**2.0d0-1.0d0) * ( (fgamma-1.0d0)/2.0d0*mach(i) - 1.0d0 )**2.0d0

    end do

end subroutine vanleer

subroutine vanleer_explicit_1stOrder

    use shared
    implicit none


    integer(kind=4) :: i
    integer(kind=4) :: gs = 0
    integer(kind=4) :: neq = 3

    real(kind=8), allocatable, dimension(:,:) :: q_aux

    gs = ng/2  !Ghosts in each side

    allocate(q_aux(total_mesh_points + ng,neq))

    q_aux(:,:) = q(:,:)


    do i = gs+1, total_mesh_points + gs

        q(i,1) = q_aux(i,1) - time_step*( (f_pls(i,1) - f_pls(i-1,1) + f_min(i+1,1) - f_min(i,1))/dx )
        q(i,2) = q_aux(i,2) - time_step*( (f_pls(i,2) - f_pls(i-1,2) + f_min(i+1,2) - f_min(i,2))/dx )
        q(i,3) = q_aux(i,3) - time_step*( (f_pls(i,3) - f_pls(i-1,3) + f_min(i+1,3) - f_min(i,3))/dx )

    end do


    deallocate(q_aux)

end subroutine vanleer_explicit_1stOrder

subroutine vanleer_explicit_2stOrder

    use shared
    implicit none


    integer(kind=4) :: i
    integer(kind=4) :: gs = 0
    integer(kind=4) :: neq = 3

    real(kind=8), allocatable, dimension(:,:) :: q_aux

    gs = ng/2  !Ghosts in each side

    allocate(q_aux(total_mesh_points + ng,neq))

    q_aux(:,:) = q(:,:)


    do i = gs+1, total_mesh_points + gs

        q(i,1) = q_aux(i,1) - time_step/(2.0d0*dx)*( (3.0d0*f_pls(i,1) - 4.0d0*f_pls(i-1,1) + f_pls(i-2,1) ) + (-3.0d0*f_min(i,1) + 4.0d0*f_min(i+1,1) - f_min(i+2,1)) )
        q(i,2) = q_aux(i,2) - time_step/(2.0d0*dx)*( (3.0d0*f_pls(i,2) - 4.0d0*f_pls(i-1,2) + f_pls(i-2,2) ) + (-3.0d0*f_min(i,2) + 4.0d0*f_min(i+1,2) - f_min(i+2,2)) )
        q(i,3) = q_aux(i,3) - time_step/(2.0d0*dx)*( (3.0d0*f_pls(i,3) - 4.0d0*f_pls(i-1,3) + f_pls(i-2,3) ) + (-3.0d0*f_min(i,3) + 4.0d0*f_min(i+1,3) - f_min(i+2,3)) )

    end do


    deallocate(q_aux)

end subroutine vanleer_explicit_2stOrder

