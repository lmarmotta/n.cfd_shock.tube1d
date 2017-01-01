subroutine ausm_plus

    use shared
    implicit none


    integer(kind=4) :: i, gs
    real(kind=8) :: pl, pr, al, ar, a_half, beta, alfa, hl, hr, rhol, rhor
    real(kind=8) :: machl, machr, ul, ur, m_li, m_lip, m_lim

    real(kind=8) :: m_plus
    real(kind=8) :: p_plus

    real(kind=8) :: m_minus
    real(kind=8) :: p_minus

    real(kind=8), dimension(3) :: p
    real(kind=8), dimension(3) :: phil, phir


    gs = ng/2

    do i = gs, total_mesh_points + gs


        ! At first lets calculate the most basic properties.

        pl = (fgamma-1.0d0)*( q(i,3) - 0.5d0 * (q(i,2)**2/q(i,1)) )
        pr = (fgamma-1.0d0)*( q(i+1,3) - 0.5d0 * (q(i+1,2)**2/q(i+1,1)) )

        al = sqrt(fgamma*pl/q(i,1))
        ar = sqrt(fgamma*pr/q(i+1,1))

        ul = q(i,2)/q(i,1)
        ur = q(i+1,2)/q(i+1,1)

        rhol = q(i,1)
        rhor = q(i+1,1)

        hl = (q(i,3) + pl) / rhol
        hr = (q(i+1,3) + pr) / rhor


        ! Now lets calculate, by roe's avgs, the sound speed ate half point.

        a_half = sqrt(al*ar)


        ! Here there are the constants.

        alfa = 3.0d0/16.0d0
        beta = 1.0d0/8.0d0


        ! Lets start a bunch of if statments, so lets now define our mach number

        machl = ul/a_half
        machr = ur/a_half


        ! Let the conditional statments start !!  Get the lefts !

        if (abs(machl) >= 1.0d0) then

            m_plus = 0.5d0*(machl+abs(machl))
            p_plus = 0.5d0*(1.0d0+(machl/abs(machl)))

        else 

            m_plus = 0.25d0*(machl+1.0d0)**2.0d0 + beta*((machl**2.0d0)-1.0d0)**2.0d0
            p_plus = 0.25d0*(machl+1.0d0)**2.0d0*(2.0d0-machl) + alfa*machl*((machl**2.0d0)-1.0d0)**2.0d0


        end if


        ! Get the Rights !

        if (abs(machr) >= 1.0d0) then

            m_minus = 0.5d0*(machr-abs(machr))
            p_minus = 0.5d0*(1.0d0-(machr/abs(machr)))

        else 

            m_minus = - 0.25d0*(machr-1.0d0)**2.0d0 - beta*((machr**2.0d0)-1.0d0)**2.0d0
            p_minus =   0.25d0*(machr-1.0d0)**2.0d0*(2.0d0+machr) - alfa*machr*((machr**2.0d0)-1.0d0)**2.0d0

        end if


        ! Last steps of the method.

        m_li = m_plus + m_minus

        m_lip = 0.5d0*(m_li + abs(m_li))
        m_lim = 0.5d0*(m_li - abs(m_li))


        ! Define our p vector.

        p(1) = 0.0d0
        p(2) = p_plus*pl + p_minus*pr
        p(3) = 0.0d0


        ! Lets do the phi.

        phil(1) = q(i,1)
        phil(2) = q(i,2)
        phil(3) = q(i,1)*hl

        phir(1) = q(i+1,1)
        phir(2) = q(i+1,2)
        phir(3) = q(i+1,1)*hr


        ! Finally, the fluxes...

        f(i,1) = a_half*( m_lip*phil(1) + m_lim*phir(1) ) + p(1)
        f(i,2) = a_half*( m_lip*phil(2) + m_lim*phir(2) ) + p(2)
        f(i,3) = a_half*( m_lip*phil(3) + m_lim*phir(3) ) + p(3)

    end do

end subroutine ausm_plus

subroutine ausm_explicit_1stOrder

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

end subroutine  ausm_explicit_1stOrder
