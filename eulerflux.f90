subroutine fluxCalculation

    use shared
    implicit none

    integer(kind=4) :: i, gs
    real(kind=8) :: u, p, ei


    gs = ng/2  !Ghosts in each side

    ! Calculate the fluxes for all the points.

    do i = 1, total_mesh_points + ng
        

        ! Initial physical malabarisms.

        u    = q(i,2) / q(i,1)
        ei   = ( q(i,3) - q(i,1) * ( 0.5d0*( u**2 ) ) ) / q(i,1)
        p    = (fgamma-1.0d0)*( q(i,3) - 0.5d0 * (q(i,2)**2/q(i,1)) )


        ! Lets do the fluxes.

        f(i,1) = q(i,1)*u
        f(i,2) = ( q(i,1)*(u**2) ) + p
        f(i,3) = u * ( q(i,3) + p )

        
        ! All we (and Jameson) need is pressure.

        press(i) = p

    end do


end subroutine fluxCalculation

subroutine artificialDissip

    use shared
    implicit none

    ! General loop variables.
    integer(kind=4) :: i, gs

    ! Jameson Artificial dissipation variables.
    real(kind=8), dimension(total_mesh_points + ng) :: eps2_ph, eps4_ph
    real(kind=8), dimension(total_mesh_points + ng) :: eta
    real(kind=8), dimension(total_mesh_points + ng,3) :: diff2
    real(kind=8), dimension(total_mesh_points + ng,3) :: diff4

    ! JST dissipation constant.
    real(kind=8) :: jst_k2
    real(kind=8) :: jst_k4

    jst_k2 = 1.0d0/4.0d0
    jst_k4 = 1.0d0/256.0d0


    gs = ng/2  !Ghosts in each side

    art_dissip(:,:) = 0.0d0
    art_dissip2d(:,:) = 0.0d0
    art_dissip4d(:,:) = 0.0d0


    ! Lets now add the artificial dissipation to the internal points. Since the
    ! artificial dissipation for the implicit scheme is added directly in its 
    ! own subroutine we should not run the artificial dissipation two times.

    ! Bare solution with no artificial dissipation.

    if (dissp_scheme == 0) then

        do i = gs+1, total_mesh_points
            q(i,1) = q(i,1) 
            q(i,2) = q(i,2)
            q(i,3) = q(i,3)
        end do


    ! Jameson's non-linear artificial dissipation.

    else if (dissp_scheme == 1) then


        ! Get dependencies for Jameson's artificial dissipation. As much ghosts 
        ! as possible... so, just one of then.

        do i = 1, total_mesh_points
            eta(i+gs) = dabs(press(i+1+gs) - 2.0d0*press(i+gs) + press(i-1+gs)) / &
                    (dabs(press(i+1+gs)) + 2.0d0*dabs(press(i+gs)) + dabs(press(i-1+gs)))
        end do


        do i = 1, total_mesh_points
            eps2_ph(i+gs) = jst_k2*max(eta(i+1+gs),eta(i+gs))
            eps4_ph(i+gs) = max(0.0d0,(jst_k4-eps2_ph(i-1+gs)))
        end do


        ! Build the artificial dissipation vector for all internal points.

        do i = 1, total_mesh_points


            ! This scratch of code deals with the full non-linaer jameson
            ! artificial dissipation.

            jst_d(i+gs,1) = (dx/time_step) * ( eps2_ph(i+1+gs)*( q(i+1+gs,1) - q(i+gs,1) ) - & 
                eps4_ph(i+1+gs) * ( q(i+2+gs,1) - 3.0d0*q(i+1+gs,1) + 3.0d0*q(i+gs,1) - q(i-1+gs,1) ) )

            jst_d(i+gs,2) = (dx/time_step) * ( eps2_ph(i+1+gs)*( q(i+1+gs,2) - q(i+gs,2) ) - & 
                eps4_ph(i+1+gs) * ( q(i+2+gs,2) - 3.0d0*q(i+1+gs,2) + 3.0d0*q(i+gs,2) - q(i-1+gs,2) ) )

            jst_d(i+gs,3) = (dx/time_step) * ( eps2_ph(i+1+gs)*( q(i+1+gs,3) - q(i+gs,3) ) - & 
                eps4_ph(i+1+gs) * ( q(i+2+gs,3) - 3.0d0*q(i+1+gs,3) + 3.0d0*q(i+gs,3) - q(i-1+gs,3) ) )


            ! Now lets take the second and fourth differences dissipation for
            ! our implicit scheme.

            diff2(i+gs,1) = dissip_omega*( q(i+1+gs,1) - q(i+gs,1) )
            diff2(i+gs,2) = dissip_omega*( q(i+1+gs,2) - q(i+gs,2) )
            diff2(i+gs,3) = dissip_omega*( q(i+1+gs,3) - q(i+gs,3) )

            diff4(i+gs,1) = dissip_omega*( q(i+2+gs,1) - 3.0d0*q(i+1+gs,1) + 3.0d0*q(i+gs,1) - q(i-1+gs,1) )
            diff4(i+gs,2) = dissip_omega*( q(i+2+gs,2) - 3.0d0*q(i+1+gs,2) + 3.0d0*q(i+gs,2) - q(i-1+gs,2) )
            diff4(i+gs,3) = dissip_omega*( q(i+2+gs,3) - 3.0d0*q(i+1+gs,3) + 3.0d0*q(i+gs,3) - q(i-1+gs,3) )

        end do


        ! Calculate the dissp for the non-linear mescled second and fourth
        ! derivatives.

        do i = 1, total_mesh_points
            art_dissip(i+gs,1) = jst_d(i+1+gs,1) - jst_d(i-1+gs,1)
            art_dissip(i+gs,2) = jst_d(i+1+gs,2) - jst_d(i-1+gs,2)
            art_dissip(i+gs,3) = jst_d(i+1+gs,3) - jst_d(i-1+gs,3)
        end do


        ! Building the second and fourth differences artificial dissipations.

        do i = 1, total_mesh_points

            art_dissip2d(i+gs,1) = diff2(i+1+gs,1) - diff2(i-1+gs,1)
            art_dissip2d(i+gs,2) = diff2(i+1+gs,2) - diff2(i-1+gs,2)
            art_dissip2d(i+gs,3) = diff2(i+1+gs,3) - diff2(i-1+gs,3)

            art_dissip4d(i+gs,1) = diff4(i+1+gs,1) - diff4(i-1+gs,1)
            art_dissip4d(i+gs,2) = diff4(i+1+gs,2) - diff4(i-1+gs,2)
            art_dissip4d(i+gs,3) = diff4(i+1+gs,3) - diff4(i-1+gs,3)

        end do


    ! Second difference linear artificial dissipation.

    else if (dissp_scheme == 2) then

        do i = 1, total_mesh_points

            jst_d(i+gs,1) = q(i+1+gs,1) - q(i+gs,1) 
            jst_d(i+gs,2) = q(i+1+gs,2) - q(i+gs,2) 
            jst_d(i+gs,3) = q(i+1+gs,3) - q(i+gs,3) 

        end do

        do i = 1, total_mesh_points

            art_dissip(i+gs,1) = dissip_omega*( jst_d(i+1+gs,1) - jst_d(i-1+gs,1) )
            art_dissip(i+gs,2) = dissip_omega*( jst_d(i+1+gs,2) - jst_d(i-1+gs,2) )
            art_dissip(i+gs,3) = dissip_omega*( jst_d(i+1+gs,3) - jst_d(i-1+gs,3) )
            
        end do


    ! Fourth difference linear artificial dissipation.

    else if (dissp_scheme == 3) then
        
        do i = 1, total_mesh_points
        
            jst_d(i+gs,1) = q(i+2+gs,1) - 3.0d0*q(i+1+gs,1) + 3.0d0*q(i+gs,1) - q(i-1+gs,1)
            jst_d(i+gs,2) = q(i+2+gs,2) - 3.0d0*q(i+1+gs,2) + 3.0d0*q(i+gs,2) - q(i-1+gs,2)
            jst_d(i+gs,3) = q(i+2+gs,3) - 3.0d0*q(i+1+gs,3) + 3.0d0*q(i+gs,3) - q(i-1+gs,3)


        end do

        do i = 1, total_mesh_points

            art_dissip(i+gs,1) = -dissip_omega*( jst_d(i+1+gs,1) - jst_d(i-1+gs,1) )
            art_dissip(i+gs,2) = -dissip_omega*( jst_d(i+1+gs,2) - jst_d(i-1+gs,2) )
            art_dissip(i+gs,3) = -dissip_omega*( jst_d(i+1+gs,3) - jst_d(i-1+gs,3) )

        end do


    end if

end subroutine artificialDissip

subroutine calc_pressure

    use shared
    implicit none

    integer(kind=4) :: i, gs
    real(kind=8) :: p


    ! Ghosts in each side

    gs = ng/2  

    ! Calculate the pressure.

    do i = 1, total_mesh_points + ng

        p = (fgamma-1.0d0)*( q(i,3) - 0.5d0 * (q(i,2)**2/q(i,1)) )

        press(i) = p

    end do

end subroutine calc_pressure
