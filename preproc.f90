subroutine indat

    use shared
    implicit none


    ! Read namelist parameters declared in the module.

    namelist /PAR_physical/ p1, p4, fgamma, R_const, rho1, rho4, F_Cp, F_Cv
    namelist /PAR_geometry/ total_mesh_points, start_mesh_point, &
        final_mesh_point, print_step
    namelist /PAR_numeric/ time_step, iterations, scheme, dissip_omega
    namelist /PAR_dissip/ dissip_omega, dissp_scheme

    open(1,file='input.in')

    read(1,PAR_physical)
    read(1,PAR_geometry)
    read(1,PAR_numeric)
    read(1,PAR_dissip)


    close(1)

end subroutine indat

subroutine calcInitialProperties

    use shared
    implicit none


    ! Calculate the energies based on input variables.

    e1 = p1/(fgamma-1.0d0)
    e4 = p4/(fgamma-1.0d0)


    ! Lets see the results.

    write(*,*) ""
    write(*,*) "+--------------------------------------------------------------+"
    write(*,*) "|            +++ Initial properties are +++                    |"
    write(*,*) "+--------------------------------------------------------------+"

    write(*,'(A,F7.4)') "   p1   = ", p1 
    write(*,'(A,F7.4)') "   p4   = ", p4

    write(*,'(A,F7.4)') "   rho1 = ", rho1
    write(*,'(A,F7.4)') "   rho4 = ", rho4

    write(*,'(A,F7.4)') "   e1   = ", e1
    write(*,'(A,F7.4)') "   e4   = ", e4


end subroutine calcInitialProperties

subroutine createMeshPoints

    use shared
    implicit none

    integer(kind=4) :: i
    real(kind=8) :: length


    write(*,*) ""
    write(*,*) "+--------------------------------------------------------------+"
    write(*,*) "|                +++ Creating mesh points +++                  |"
    write(*,*) "+--------------------------------------------------------------+"


    length = (final_mesh_point - start_mesh_point)
    dx = length / real(total_mesh_points,8)


    
    mesh_points(1) = start_mesh_point


    ! Creating mesh points based on start point.

    do i = 2, total_mesh_points + ng
        mesh_points(i) = start_mesh_point + (real(i, 8)*dx)
    end do

    write(*,'(A,I7)')    "   Number of mesh points: ", total_mesh_points
    write(*,'(A,F10.4)') "   The mesh spacing dx = ", dx

end subroutine createMeshPoints

subroutine createInitialCondition

    use shared
    implicit none

    integer(kind=4) :: i


    write(*,*) ""
    write(*,*) "+--------------------------------------------------------------+"
    write(*,*) "|             +++ Creating initial condition +++               |"
    write(*,*) "+--------------------------------------------------------------+"


    ! Put all initial conditions in the properties vectors.
    
    do i = 1, total_mesh_points + ng
        
        if (mesh_points(i) <= 0.0d0) then

            q(i,1) = rho1
            q(i,2) = 0.0d0
            q(i,3) = e1

            f(i,1) = 0.0d0
            f(i,2) = 0.0d0
            f(i,3) = 0.0d0

        else if (mesh_points(i) > 0.0d0) then

            q(i,1) = rho4
            q(i,2) = 0.0d0
            q(i,3) = e4

            f(i,1) = 0.0d0
            f(i,2) = 0.0d0
            f(i,3) = 0.0d0

        end if

    end do

end subroutine createInitialCondition

subroutine sound_speed

    use shared
    implicit none

    integer(kind=4) :: i
    real(kind=8) :: p


    ! Simple as that my dear.

    do i = 1, total_mesh_points + ng

        p    = (fgamma-1.0d0)*( q(i,3) - 0.5d0 * (q(i,2)**2/q(i,1)) )
        a_speed(i) = sqrt(fgamma*p/q(i,1))

    end do

end subroutine sound_speed

subroutine machcalc

    use shared
    implicit none

    integer(kind=4) :: i
    real(kind=8) :: u
    
    do i = 1, total_mesh_points + ng

        u = q(i,2)/q(i,1)
        mach(i) = u/a_speed(i)

    end do

end subroutine machcalc

