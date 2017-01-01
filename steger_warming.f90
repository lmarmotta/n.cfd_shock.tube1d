subroutine steger_warming

    use shared
    implicit none

    real(kind=8) :: m1, m2, m3, rho, H, u, p
    integer(kind=4) :: i, gs

    ! Lambda vectors.
    real(kind=8), allocatable, dimension(:,:) :: lmbd_pls
    real(kind=8), allocatable, dimension(:,:) :: lmbd_min

    gs = ng/2  !Ghosts in each side


    allocate(lmbd_pls(total_mesh_points + ng, 3))
    allocate(lmbd_min(total_mesh_points + ng, 3))

    
    lmbd_pls = 0.0d0
    lmbd_min = 0.0d0

    do i = 1, total_mesh_points + ng


        ! Bizurate some properties.

        u    = q(i,2) / q(i,1)


        ! Now lets proceed with our lambdas (positives).

        lmbd_pls(i,1) = ( (u - a_speed(i) ) + abs( u - a_speed(i) ) ) / 2.0d0
        lmbd_pls(i,2) =                 ( u + abs(u) )                / 2.0d0
        lmbd_pls(i,3) = ( (u + a_speed(i) ) + abs( u + a_speed(i) ) ) / 2.0d0


        ! Now lets proceed with our lambdas (negative).

        lmbd_min(i,1) = ( (u - a_speed(i) ) - abs( u - a_speed(i) ) ) / 2.0d0
        lmbd_min(i,2) =                 ( u - abs(u) )                / 2.0d0
        lmbd_min(i,3) = ( (u + a_speed(i) ) - abs( u + a_speed(i) ) ) / 2.0d0

    end do


    ! Now lets build the positive and negative fluxes.

    do i = 1, total_mesh_points + ng
        

        ! Bizurate some properties.

        rho  = q(i,1)
        u    = q(i,2) / q(i,1)
        p    = (fgamma-1.0d0)*( q(i,3) - 0.5d0 * rho*u*u )
        H    = (q(i,3) + p) / rho


        ! Vectors multiplication constants.

        m1 = lmbd_pls(i,1) * rho/(2.0d0*fgamma)
        m2 = lmbd_pls(i,2) * ( rho*(fgamma-1.0d0) )/ fgamma
        m3 = lmbd_pls(i,3) * rho/(2.0d0*fgamma)


        ! Build the positive fluxes.

        f_pls(i,1) =      m1 * 1.0d0          +      m2 * 1.0d0         +  m3 *     1.0d0
        f_pls(i,2) =  m1 * (u - a_speed(i))   +        m2 * u           +  m3 *  (u + a_speed(i))
        f_pls(i,3) = m1 * (H - u*a_speed(i))  +  m2 * 0.5d0*(u**2.0d0)  +  m3 * (H + u*a_speed(i))


        ! Vectors multiplication constants.

        m1 = lmbd_min(i,1) * rho/(2.0d0*fgamma)
        m2 = lmbd_min(i,2) * ( rho*(fgamma-1.0d0) )/ fgamma
        m3 = lmbd_min(i,3) * rho/(2.0d0*fgamma)


        ! Build the negative fluxes.

        f_min(i,1) =      m1 * 1.0d0          +      m2 * 1.0d0       +  m3 *     1.0d0
        f_min(i,2) =  m1 * (u - a_speed(i))   +        m2 * u         +  m3 *  (u + a_speed(i))
        f_min(i,3) = m1 * (H - u*a_speed(i))  +  m2 * 0.5d0*u**2.0d0  +  m3 * (H + u*a_speed(i))

    end do


    deallocate(lmbd_pls)
    deallocate(lmbd_min)

end subroutine steger_warming

subroutine stegerWarming_explicit_1stOrder

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

        q(i,1) = q_aux(i,1) - time_step*(f_pls(i,1) - f_pls(i-1,1) + f_min(i+1,1) - f_min(i,1))/dx 
        q(i,2) = q_aux(i,2) - time_step*(f_pls(i,2) - f_pls(i-1,2) + f_min(i+1,2) - f_min(i,2))/dx 
        q(i,3) = q_aux(i,3) - time_step*(f_pls(i,3) - f_pls(i-1,3) + f_min(i+1,3) - f_min(i,3))/dx 

    end do


    deallocate(q_aux)

end subroutine stegerWarming_explicit_1stOrder

subroutine stegerWarming_explicit_2stOrder

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

end subroutine stegerWarming_explicit_2stOrder

subroutine implicit_fvsfds_1stOrder
    
    use shared
    implicit none

    integer(kind=4) :: i,gs

    real(kind=8) :: rho
    real(kind=8) :: u  
    real(kind=8) :: p  
    real(kind=8) :: H  
    real(kind=8) :: h_dtdx

    ! Lambda vectors.
    real(kind=8), allocatable, dimension(:,:) :: lmbd_pls
    real(kind=8), allocatable, dimension(:,:) :: lmbd_min

    real(kind=8), dimension(3,3,total_mesh_points + ng) :: A_p
    real(kind=8), dimension(3,3,total_mesh_points + ng) :: A_m

    real(kind=8), dimension(3,3,total_mesh_points + ng) :: L_p
    real(kind=8), dimension(3,3,total_mesh_points + ng) :: L_m

    ! For the system.
    real(kind=8), dimension(3,3) :: identity

    real(kind=8), dimension(3,total_mesh_points)     :: RHS
    real(kind=8), dimension(3,total_mesh_points)     :: du
    real(kind=8), dimension(3,3,total_mesh_points)   :: maind
    real(kind=8), dimension(3,3,total_mesh_points)   :: lowerd
    real(kind=8), dimension(3,3,total_mesh_points)   :: upperd



    allocate(lmbd_pls(total_mesh_points + ng, 3))
    allocate(lmbd_min(total_mesh_points + ng, 3))


    h_dtdx = 0.5d0*(time_step/dx)


    gs = ng/2  !Ghosts in each side


    call identityCreate(3, identity)


    do i = 1, total_mesh_points + ng

        rho  = q(i,1)
        u    = q(i,2) / q(i,1)
        p    = (fgamma-1.0d0)*( q(i,3) - 0.5d0 * rho*u*u )
        H    = (q(i,3) + p) / rho

        call fvs_transformation_matrix(rho,a_speed(i),u,H,i)

    end do


    ! Calculate the labdas for our matrix.

    do i = 1, total_mesh_points + ng


        ! Bizurate some properties.

        u    = q(i,2) / q(i,1)


        ! Now lets proceed with our lambdas (positives).

        lmbd_pls(i,1) = ( (u - a_speed(i) ) + abs( u - a_speed(i) ) ) / 2.0d0
        lmbd_pls(i,2) =                 ( u + abs(u) )                / 2.0d0
        lmbd_pls(i,3) = ( (u + a_speed(i) ) + abs( u + a_speed(i) ) ) / 2.0d0


        ! Now lets proceed with our lambdas (negative).

        lmbd_min(i,1) = ( (u - a_speed(i) ) - abs( u - a_speed(i) ) ) / 2.0d0
        lmbd_min(i,2) =                 ( u - abs(u) )                / 2.0d0
        lmbd_min(i,3) = ( (u + a_speed(i) ) - abs( u + a_speed(i) ) ) / 2.0d0

    end do


    ! Calculate lambda matrix.

    L_p = 0.0d0
    L_m = 0.0d0

    do i = 1, total_mesh_points + ng

        
        L_p(1,1,i) = lmbd_pls(i,1)
        L_p(2,2,i) = lmbd_pls(i,2)
        L_p(3,3,i) = lmbd_pls(i,3)


        L_m(1,1,i) = lmbd_min(i,1)
        L_m(2,2,i) = lmbd_min(i,2)
        L_m(3,3,i) = lmbd_min(i,3)

    end do 


    ! Now lets play with our transformation matrix.

    do i = 1, total_mesh_points + ng

        A_p(:,:,i) = matmul(matmul(p_matrix(:,:,i),L_p(:,:,i)),p_matrix_inv(:,:,i))
        A_m(:,:,i) = matmul(matmul(p_matrix(:,:,i),L_m(:,:,i)),p_matrix_inv(:,:,i))

    end do


    ! Now we will start to build our system of equations. Here I use the
    ! complete operator that is not pre-factored. Ok... I was lazu enouth to
    ! pass my job to the computer !

    do i = 1, total_mesh_points

        RHS(1,i) = - time_step/dx * (f_pls(i+gs,1)-f_pls(i-1+gs,1) + f_min(i+1+gs,1)-f_min(i+gs,1))
        RHS(2,i) = - time_step/dx * (f_pls(i+gs,2)-f_pls(i-1+gs,2) + f_min(i+1+gs,2)-f_min(i+gs,2))
        RHS(3,i) = - time_step/dx * (f_pls(i+gs,3)-f_pls(i-1+gs,3) + f_min(i+1+gs,3)-f_min(i+gs,3))
        
    end do


    ! Main diagonal at first.

    do i = 1, total_mesh_points

        maind(:,:,i) = identity(:,:) + (time_step/dx)*(A_p(:,:,i+gs) - A_m(:,:,i+gs))

    end do


    ! Lower diagonal

    do i = 1, total_mesh_points

        lowerd(:,:,i) = -(time_step/dx)*A_p(:,:,i-1+gs)

    end do


    do i = 1, total_mesh_points

        upperd(:,:,i) = (time_step/dx)*A_m(:,:,i+1+gs)

    end do
    

    du = 0.0d0

    call cc299blktriad_opt(maind,lowerd,upperd,3,total_mesh_points,RHS,du)


    ! Its all done for now, lets now put our new time solution in place.

    do i = 1, total_mesh_points

        q(i+gs,1) = q(i+gs,1) + du(1,i)   
        q(i+gs,2) = q(i+gs,2) + du(2,i)   
        q(i+gs,3) = q(i+gs,3) + du(3,i)   

    end do

    deallocate(lmbd_pls)
    deallocate(lmbd_min)

end subroutine implicit_fvsfds_1stOrder

subroutine fvs_transformation_matrix(rho,a,u,h,i)

    ! Hirsch p.419.
    ! rho: Density.
    ! a  : Speed of sound.
    ! u  : u speed component.
    ! h  : entalpy

    use shared
    implicit none

    integer(kind=4) :: i

    real(kind=8) :: rho
    real(kind=8) :: a 
    real(kind=8) :: u
    real(kind=8) :: h  



    ! Calculating the normal matrix.

    p_matrix(1,1,i) =   1.0d0
    p_matrix(1,2,i) =   rho/2.0d0*a
    p_matrix(1,3,i) = - rho/2.0d0*a
     
    p_matrix(2,1,i) =    u
    p_matrix(2,2,i) =   (rho*(u+a))/(2.0d0*a)
    p_matrix(2,3,i) = - (rho*(u-a))/(2.0d0*a)
    
    p_matrix(3,1,i) =   (u**2.0d0)/2.0d0
    p_matrix(3,2,i) =   rho*(h+u*a)/(2.0d0*a)
    p_matrix(3,3,i) = - rho*(h-u*a)/(2.0d0*a)


    ! Inverting matrix to calc transformation

    call inv(p_matrix(:,:,i),p_matrix_inv(:,:,i),3)


end subroutine fvs_transformation_matrix
