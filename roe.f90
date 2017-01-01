subroutine roe_flux

    use shared
    implicit none

    integer(kind=4) :: i, gs
    real(kind=8) :: rhol, rhor, ul, ur, hl, hr, pl, pr

    ! Roe's averages vectors.
    real(kind=8) :: u_hat
    real(kind=8) :: h_hat
    real(kind=8) :: c_hat

    real(kind=8) :: rho_hat

    ! Eigenvalues.
    real(kind=8), dimension(3) :: lambda_hat

    ! Eigenvectors.
    ! r_hat(a,b,c):
    !     a = interfaces index.
    !     b = which one of the three eigenvectors.
    !     c = Element inside each eigenvector.
    real(kind=8), dimension(3,3) :: r_hat

    ! The intensities.
    real(kind=8), dimension(3) :: dv_hat

    ! Variation of properties.
    real(kind=8) :: dp
    real(kind=8) :: du
    real(kind=8) :: drho


    ! Ghosts in each side

    gs = ng/2  


    ! Calculate all the things for one face and step man, Fabio's tip !

    !do i = 1, total_mesh_points + (ng-1)
    do i = gs, total_mesh_points + gs


        ! Lets strt with the Roe averages man. Now I dont know if this will work si
        ! nce I am using averages points in ponints not in interfaces. If things go
        ! es wrong, this type of thing would be a good start during debugging. 
        ! Note: left  = property(i)
        !       right = property(i+1)


        ! Bizurate some properties.

        rhol = q(i,1)
        rhor = q(i+1,1)

        ul = q(i,2) / q(i,1)
        ur = q(i+1,2) / q(i+1,1)

        pl = press(i)
        pr = press(i+1)

        hl = (q(i,3) + pl) / rhol
        hr = (q(i+1,3) + pr) / rhor


        ! Build avgs..

        rho_hat = sqrt(rhol*rhor)

        u_hat = ( (sqrt(rhol)*ul) + (sqrt(rhor)*ur) ) / ( sqrt(rhol) + sqrt(rhor) )
        h_hat = ( (sqrt(rhol)*hl) + (sqrt(rhor)*hr) ) / ( sqrt(rhol) + sqrt(rhor) )
        c_hat = sqrt( (fgamma-1.0d0) * ( h_hat - 0.5d0*( u_hat**2.0d0) ) )


        ! We are almost done. Before we proceed with the intensity calculation,
        ! lets first build a very usefull dp vector.

        drho = rhor - rhol

        du = ur - ul
        dp = pr - pl


        ! Calculating the eigenvalues.

        lambda_hat(1) = u_hat - c_hat
        lambda_hat(2) = u_hat
        lambda_hat(3) = u_hat + c_hat


        ! Calculating the eigenvectors.

        r_hat(1,1) = 1.0d0
        r_hat(1,2) = u_hat - c_hat
        r_hat(1,3) = h_hat - ( u_hat*c_hat )

        r_hat(2,1) = 1.0d0
        r_hat(2,2) = u_hat
        r_hat(2,3) = 0.5d0*(u_hat**2.0d0)

        r_hat(3,1) = 1.0d0 
        r_hat(3,2) = u_hat + c_hat
        r_hat(3,3) = h_hat + ( u_hat*c_hat )


        ! Finally, lets calculate the intensities.

        dv_hat(1) = 0.5d0 * (dp - rho_hat*c_hat*du ) / c_hat**2.0d0
        dv_hat(2) = - (dp - c_hat**2.0d0*drho)/c_hat**2.0d0
        dv_hat(3) = 0.5d0 * (dp + rho_hat*c_hat*du ) / c_hat**2.0d0


        ! Now the fluxes...

        f(i,1) = 0.5d0 * ( f(i,1) + f(i+1,1) ) - 0.5d0 * (  &
            abs(lambda_hat(1)) * dv_hat(1) * r_hat(1,1) + &
            abs(lambda_hat(2)) * dv_hat(2) * r_hat(2,1) + &
            abs(lambda_hat(3)) * dv_hat(3) * r_hat(3,1)   )

        f(i,2) = 0.5d0 * ( f(i,2) + f(i+1,2) ) - 0.5d0 * (  &
            abs(lambda_hat(1)) * dv_hat(1) * r_hat(1,2) + &
            abs(lambda_hat(2)) * dv_hat(2) * r_hat(2,2) + &
            abs(lambda_hat(3)) * dv_hat(3) * r_hat(3,2)   ) 

        f(i,3) = 0.5d0 * ( f(i,3) + f(i+1,3) ) - 0.5d0 * (  &
            abs(lambda_hat(1)) * dv_hat(1) * r_hat(1,3) + &
            abs(lambda_hat(2)) * dv_hat(2) * r_hat(2,3) + &
            abs(lambda_hat(3)) * dv_hat(3) * r_hat(3,3)   )


    end do


    ! Lets separate the fluxes, this way we can use the same
    ! implicit time-stepping subroutine.

    if (scheme == 13) then

        do i = gs, total_mesh_points + gs

            f_pls(i,1) = f(i+1,1)
            f_pls(i,2) = f(i+1,2)
            f_pls(i,1) = f(i+1,3)

            f_min(i,1) = f(i,1)
            f_min(i,2) = f(i,2)
            f_min(i,1) = f(i,3)

        end do

    end if


end subroutine roe_flux

subroutine roe_explicit_1stOrder

    use shared
    implicit none

    integer(kind=4) :: i, gs
    integer(kind=4) :: neq = 3

    real(kind=8), allocatable, dimension(:,:) :: q_aux

    gs = ng/2  !Ghosts in each side

    allocate(q_aux(total_mesh_points + ng,neq))

    q_aux(:,:) = q(:,:)

    do i = gs+1, total_mesh_points + gs

        q(i,1) = q_aux(i,1) - (time_step/dx) * ( f(i,1) - f(i-1,1) ) 
        q(i,2) = q_aux(i,2) - (time_step/dx) * ( f(i,2) - f(i-1,2) ) 
        q(i,3) = q_aux(i,3) - (time_step/dx) * ( f(i,3) - f(i-1,3) ) 

    end do


    deallocate(q_aux)

end subroutine roe_explicit_1stOrder
