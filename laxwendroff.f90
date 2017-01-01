subroutine laxWendroff

    use shared
    implicit none

    integer(kind=4) :: i,ii,jj
    integer(kind=4) :: gs = 0

    real(kind=8) :: sigma1, sigma2

    real(kind=8), allocatable, dimension(:,:) :: q_aux
    real(kind=8), allocatable, dimension(:,:) :: a_plus_half
    real(kind=8), allocatable, dimension(:,:) :: a_minu_half
    real(kind=8), allocatable, dimension(:) :: mult_r1
    real(kind=8), allocatable, dimension(:) :: mult_r2
    real(kind=8), allocatable, dimension(:) :: difference1
    real(kind=8), allocatable, dimension(:) :: difference2


    allocate(q_aux(total_mesh_points + ng, 3))
    allocate(a_plus_half(3,3))
    allocate(a_minu_half(3,3))
    allocate(mult_r1(3))
    allocate(mult_r2(3))
    allocate(difference1(3))
    allocate(difference2(3))


    q_aux(:,:) = q(:,:)


    ! Lets prepare ou vectors for multiplication.

    gs = ng/2  ! Ghosts in each side.

    do i = gs+1, total_mesh_points-1


        ! Calculate the jacobians at middle points.

        do ii = 1,3
            do jj = 1,3
                a_plus_half(ii,jj) = ( a_j(ii,jj,i) + a_j(ii,jj,i+1) ) / 2.0d0
                a_minu_half(ii,jj) = ( a_j(ii,jj,i) + a_j(ii,jj,i-1) ) / 2.0d0
            end do
        end do


        ! Here Ill separate the matrix multiplication.

        difference1(1) = f(i+1,1) - f(i,1) !- art_dissip(i,1)
        difference1(2) = f(i+1,2) - f(i,2) !- art_dissip(i,2)
        difference1(3) = f(i+1,3) - f(i,3) !- art_dissip(i,3)

        difference2(1) = f(i,1) - f(i-1,1) !- art_dissip(i,1)
        difference2(2) = f(i,2) - f(i-1,2) !- art_dissip(i,2)
        difference2(3) = f(i,3) - f(i-1,3) !- art_dissip(i,3)


        ! In this vector we should have the result we need to proceed with
        ! time marching.

        mult_r1 = matmul(a_plus_half,difference1)
        mult_r2 = matmul(a_plus_half,difference2)


        ! Good variables to make things leaner.

        sigma1 = time_step/(2.0d0*dx) 
        sigma2 = ((time_step**2)/(2.0d0*dx**2) )


        ! Lets March the equations.

        q(i,1) = q_aux(i,1) - sigma1*(f(i+1,1) - f(i-1,1) - art_dissip(i,1)) + &
            sigma2*(mult_r1(1) - mult_r2(1))

        q(i,2) = q_aux(i,2) - sigma1*(f(i+1,2) - f(i-1,2) - art_dissip(i,2)) + &
            sigma2*(mult_r1(2) - mult_r2(2))

        q(i,3) = q_aux(i,3) - sigma1*(f(i+1,3) - f(i-1,3) - art_dissip(i,3)) + &
            sigma2*(mult_r1(3) - mult_r2(3))


    end do
    

    deallocate(a_plus_half)
    deallocate(a_minu_half)
    deallocate(mult_r1)
    deallocate(mult_r2)
    deallocate(difference1)
    deallocate(difference2)
    deallocate(q_aux)

end subroutine laxWendroff


