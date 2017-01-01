subroutine calcJacobian

    use shared
    implicit none

    integer(kind=4) :: k
    real(kind=8) :: u,rho,e


    ! Lets initialize the jacobian matrix.

    a_j(:,:,:) = 0.0d0

    do k = 1, total_mesh_points + ng

        u = q(k,2)/q(k,1)
        rho = q(k,1)
        e = q(k,3)

        a_j(1,1,k) = 0.0d0
        a_j(1,2,k) = 1.0d0
        a_j(1,3,k) = 0.0d0

        a_j(2,1,k) = 0.5d0*(fgamma-3.0d0) * u**2
        a_j(2,2,k) = (3.0d0-fgamma) * u
        a_j(2,3,k) = (fgamma-1.0d0)

        a_j(3,1,k) = - (fgamma*e*u)/rho + (fgamma-1.0d0)*u**3
        a_j(3,2,k) =   (fgamma*e)/rho - (3.0d0*(fgamma-1.0d0)*u**2)/2.0d0
        a_j(3,3,k) =    fgamma*u

    end do


end subroutine calcJacobian


