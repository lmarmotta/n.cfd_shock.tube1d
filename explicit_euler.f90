subroutine explicitCentredSimple

    use shared
    implicit none

    integer(kind=4) :: i
    integer(kind=4) :: gs = 0
    integer(kind=4) :: neq = 3

    real(kind=8), allocatable, dimension(:,:) :: q_aux

    allocate(q_aux(total_mesh_points + ng,neq))

    q_aux(:,:) = q(:,:)


    ! Lets advance in time the solution.

    gs = ng/2  !Ghosts in each side

    do i = gs+1, total_mesh_points

        q(i,1) = q_aux(i,1) + time_step*((f(i+1,1)-f(i-1,1) + art_dissip(i,1)) / (2.0d0*dx) )
        q(i,2) = q_aux(i,2) + time_step*((f(i+1,2)-f(i-1,2) + art_dissip(i,2)) / (2.0d0*dx) )
        q(i,3) = q_aux(i,3) + time_step*((f(i+1,3)-f(i-1,3) + art_dissip(i,3)) / (2.0d0*dx) )

    end do


    deallocate(q_aux)

end subroutine explicitCentredSimple


