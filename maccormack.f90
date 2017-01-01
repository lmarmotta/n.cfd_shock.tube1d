subroutine maccormack_pred

    use shared
    implicit none

    integer(kind=4) :: i
    integer(kind=4) :: gs = 0
    real(kind=8) :: u, p, ei


    ! Lets now do the predictor step of the maccormak method.

    gs = ng/2

    do i = gs+1, total_mesh_points-1
        q_pred(i,1) = q(i,1) - (time_step / dx) * (f(i+1,1) - f(i,1) - art_dissip(i,1))
        q_pred(i,2) = q(i,2) - (time_step / dx) * (f(i+1,2) - f(i,2) - art_dissip(i,2))
        q_pred(i,3) = q(i,3) - (time_step / dx) * (f(i+1,3) - f(i,3) - art_dissip(i,3))
    end do


    ! We will also need the fluxes based on our predictor step variables.

    do i = gs+1, total_mesh_points-1
        

        ! Initial physical malabarisms.

        u    = q_pred(i,2) / q_pred(i,1)
        ei   = ( q_pred(i,3) - q_pred(i,1) * ( 0.5d0*( u**2 ) ) ) / q_pred(i,1)
        p    = (fgamma-1.0d0)*( q_pred(i,3) - 0.5d0*(q_pred(i,2)**2/q_pred(i,1)) )


        ! Lets do the fluxes.

        f_pred(i,1) = q_pred(i,1)*u
        f_pred(i,2) = ( q_pred(i,1)*(u**2) ) + p
        f_pred(i,3) = u * ( q_pred(i,3) + p )
     
    end do


    ! Now we need to extrapolate the fluxes for the ghosts.
    ! Initial ghosts gs = ng/2:

    do i = 1, gs
        f_pred(i,1) = f_pred(gs+1,1)
        f_pred(i,2) = f_pred(gs+1,2)
        f_pred(i,3) = f_pred(gs+1,3)
    end do


    ! for the final two ghosts.

    do i = total_mesh_points, total_mesh_points + ng
        f_pred(i,1) = f_pred(total_mesh_points,1)
        f_pred(i,2) = f_pred(total_mesh_points,2)
        f_pred(i,3) = f_pred(total_mesh_points,3)
    end do


end subroutine maccormack_pred

subroutine maccormack_corr

    use shared
    implicit none

    integer(kind=4) :: i
    integer(kind=4) :: gs = 0


    ! Lets now march all the equations using the values of the recently calcula
    ! ted predictor fluxes.

    gs = ng/2  ! Ghosts in each side

    do i = gs+1, total_mesh_points-1

        q(i,1) = 0.5d0 * ( q(i,1) + q_pred(i,1) - (time_step/dx) * &
            (f_pred(i,1) - f_pred(i-1,1) - art_dissip(i,1)) )
        q(i,2) = 0.5d0 * ( q(i,2) + q_pred(i,2) - (time_step/dx) * &
            (f_pred(i,2) - f_pred(i-1,2) - art_dissip(i,2)) )
        q(i,3) = 0.5d0 * ( q(i,3) + q_pred(i,3) - (time_step/dx) * &
            (f_pred(i,3) - f_pred(i-1,3) - art_dissip(i,3)) )
    end do


end subroutine maccormack_corr


