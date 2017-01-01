subroutine harten_tvd_roeAvg(order)

    use shared
    implicit none


    integer(kind=4) :: gs, i, order

    real(kind=8) :: dtdx, signal

    real(kind=8) :: rhol,rhor,ul,ur,pl,pr,hl,hr,rho_hat,c1,c2

    real(kind=8), allocatable, dimension(:) :: u_hat
    real(kind=8), allocatable, dimension(:) :: h_hat
    real(kind=8), allocatable, dimension(:) :: c_hat

    real(kind=8), allocatable, dimension(:,:) :: v_k
    real(kind=8), allocatable, dimension(:,:,:) :: r
    real(kind=8), allocatable, dimension(:,:) :: alpha
    real(kind=8), allocatable, dimension(:,:) :: eps
    real(kind=8), allocatable, dimension(:,:) :: phi
    real(kind=8), allocatable, dimension(:,:) :: g_til
    real(kind=8), allocatable, dimension(:,:) :: s
    real(kind=8), allocatable, dimension(:,:) :: g
    real(kind=8), allocatable, dimension(:,:) :: gam
    real(kind=8), allocatable, dimension(:,:) :: hh

    ! Allocating Roe's avgs.
    allocate(u_hat(total_mesh_points + ng))
    allocate(h_hat(total_mesh_points + ng))
    allocate(c_hat(total_mesh_points + ng))

    ! Allocating auxiliary arrays for hartens scheme.
    allocate(v_k(3,total_mesh_points + ng))
    allocate(r(3,3,total_mesh_points + ng))
    allocate(alpha(3,total_mesh_points + ng))
    allocate(eps(3,total_mesh_points + ng))
    allocate(phi(3,total_mesh_points + ng))
    allocate(g_til(3,total_mesh_points + ng))
    allocate(s(3,total_mesh_points + ng))
    allocate(g(3,total_mesh_points + ng))
    allocate(gam(3,total_mesh_points + ng))
    allocate(hh(3,total_mesh_points + ng))


    gs = ng/2
    dtdx = time_step/dx


    ! Lets first define the Roe avgs...

    do i = 1, total_mesh_points + ng-1

       rhol = q(i,1)
       rhor = q(i+1,1)

       ul = q(i,2) / q(i,1)
       ur = q(i+1,2) / q(i+1,1)

       pl = press(i)
       pr = press(i+1)

       hl = (q(i,3) + pl) / rhol
       hr = (q(i+1,3) + pr) / rhor


       ! Build avgs..

       rho_hat = dsqrt(rhol*rhor)

       u_hat(i) = ( (dsqrt(rhol)*ul) + (dsqrt(rhor)*ur) ) / ( dsqrt(rhol) + dsqrt(rhor) )
       h_hat(i) = ( (dsqrt(rhol)*hl) + (dsqrt(rhor)*hr) ) / ( dsqrt(rhol) + dsqrt(rhor) )
       c_hat(i) = dsqrt( (fgamma-1.0d0) * ( h_hat(i) - 0.5d0*( u_hat(i)**2.0d0) ) )

    end do


    ! Now lets do the v_k.

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        v_k(1,i) = dtdx*(u_hat(i) - c_hat(i))
        v_k(2,i) = dtdx*(     u_hat(i)      )
        v_k(3,i) = dtdx*(u_hat(i) + c_hat(i))

    end do
    !$omp end do
    !$omp end parallel


    ! Now the R's.

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        r(1,1,i) = 1.0d0
        r(2,1,i) = u_hat(i) - c_hat(i)
        r(3,1,i) = h_hat(i) - u_hat(i)*c_hat(i)

        r(1,2,i) = 1.0d0
        r(2,2,i) = u_hat(i)
        r(3,2,i) = 0.5d0*u_hat(i)**2.0d0

        r(1,3,i) = 1.0d0
        r(2,3,i) = u_hat(i) + c_hat(i)
        r(3,3,i) = h_hat(i) + u_hat(i)*c_hat(i)

    end do
    !$omp end do
    !$omp end parallel


    ! At this part of the code I would calculate the eps but since I am lazy
    ! right now, I ll just make then fixed.

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        eps(1,i) = 0.005d0
        eps(2,i) = 0.0d0
        eps(3,i) = 0.0d0

    end do
    !$omp end do
    !$omp end parallel


    ! Now the alphas.

    do i = 1, total_mesh_points + ng-1

       c1 = ((fgamma-1.0d0)/(c_hat(i)**2.0d0)) * ( jump(q(i+1,3),q(i,3)) + 0.5d0 * u_hat(i)*jump(q(i+1,1),q(i,1)) - u_hat(i)*jump(q(i+1,2),q(i,2)) )
       c2 = (1.0d0/c_hat(i)) * ( jump(q(i+1,2),q(i,2)) - u_hat(i)*jump(q(i+1,1),q(i,1)) )

       alpha(1,i) = 0.5d0*( c1 - c2 )
       alpha(2,i) = jump(q(i+1,1),q(i,1)) - c1
       alpha(3,i) = 0.5d0*( c1 + c2 )

    end do


    ! For the g_til calculation we need our phis, take care with this vector
    ! because we will over write it later...

    phi = 0.0d0

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        phi(1,i) = psiHarten(v_k(1,i), eps(1,i))
        phi(2,i) = psiHarten(v_k(2,i), eps(2,i))
        phi(3,i) = psiHarten(v_k(3,i), eps(3,i))

    end do
    !$omp end do
    !$omp end parallel


    ! And finally lets calculate the g_tils

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

       g_til(1,i) = 0.5d0*(phi(1,i) - v_k(1,i)**2.0d0)*alpha(1,i) 
       g_til(2,i) = 0.5d0*(phi(2,i) - v_k(2,i)**2.0d0)*alpha(2,i) 
       g_til(3,i) = 0.5d0*(phi(3,i) - v_k(3,i)**2.0d0)*alpha(3,i) 

    end do
    !$omp end do
    !$omp end parallel


    ! Now the s's... 

    signal = 1.0d0

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        s(1,i) = sign(signal,g_til(1,i))
        s(2,i) = sign(signal,g_til(2,i))
        s(3,i) = sign(signal,g_til(3,i))

    end do
    !$omp end do
    !$omp end parallel


    ! Finally the g's !!

    do i = 2, total_mesh_points + ng-1

        if (order == 1) then

            g(1,i) = 0.0d0
            g(2,i) = 0.0d0
            g(3,i) = 0.0d0

        else if (order == 2) then

            g(1,i) = s(1,i) * max(0.0d0,min(abs(g_til(1,i)),g_til(1,i-1)*s(1,i)))
            g(2,i) = s(2,i) * max(0.0d0,min(abs(g_til(2,i)),g_til(2,i-1)*s(2,i)))
            g(3,i) = s(3,i) * max(0.0d0,min(abs(g_til(3,i)),g_til(3,i-1)*s(3,i)))

        end if

    end do


    ! Lets now calculate the gammas.

    gam = 0.0d0

    do i = 1, total_mesh_points + ng-1

            gam(1,i) = gamHarten(g(1,i), g(1,i+1), alpha(1,i))
            gam(2,i) = gamHarten(g(2,i), g(2,i+1), alpha(2,i))
            gam(3,i) = gamHarten(g(3,i), g(3,i+1), alpha(3,i))

    end do


    ! As I said, we will now re-build our phi vector to receive as arguments the
    ! v_k and gam summed... lets do it.

    phi = 0.0d0

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        phi(1,i) = psiHarten(v_k(1,i) + gam(1,i), eps(1,i))
        phi(2,i) = psiHarten(v_k(2,i) + gam(2,i), eps(2,i))
        phi(3,i) = psiHarten(v_k(3,i) + gam(3,i), eps(3,i))

    end do
    !$omp end do
    !$omp end parallel

    
    ! Lets now calculate our harten "dissipation". Note that in first order,
    ! the values of g will be zero.

    do i = 1, total_mesh_points + ng-1


        hh(1,i) = r(1,1,i)*( g(1,i) + g(1,i+1) - (phi(1,i)*alpha(1,i)) ) + &
                  r(1,2,i)*( g(2,i) + g(2,i+1) - (phi(2,i)*alpha(2,i)) ) + &
                  r(1,3,i)*( g(3,i) + g(3,i+1) - (phi(3,i)*alpha(3,i)) )

        hh(2,i) = r(2,1,i)*( g(1,i) + g(1,i+1) - (phi(1,i)*alpha(1,i)) ) + &
                  r(2,2,i)*( g(2,i) + g(2,i+1) - (phi(2,i)*alpha(2,i)) ) + &
                  r(2,3,i)*( g(3,i) + g(3,i+1) - (phi(3,i)*alpha(3,i)) )

        hh(3,i) = r(3,1,i)*( g(1,i) + g(1,i+1) - (phi(1,i)*alpha(1,i)) ) + &
                  r(3,2,i)*( g(2,i) + g(2,i+1) - (phi(2,i)*alpha(2,i)) ) + &
                  r(3,3,i)*( g(3,i) + g(3,i+1) - (phi(3,i)*alpha(3,i)) )

    end do


    ! We are almost done now.

    do i = 1, total_mesh_points + ng-1

        f_pls(i,1) = 0.5d0*( f(i,1) + f(i+1,1) + (1.0d0/dtdx)*hh(1,i) )
        f_pls(i,2) = 0.5d0*( f(i,2) + f(i+1,2) + (1.0d0/dtdx)*hh(2,i) )
        f_pls(i,3) = 0.5d0*( f(i,3) + f(i+1,3) + (1.0d0/dtdx)*hh(3,i) )

    end do

    deallocate(u_hat)
    deallocate(h_hat)
    deallocate(c_hat)
    deallocate(v_k)
    deallocate(r)
    deallocate(alpha)
    deallocate(eps)
    deallocate(phi)
    deallocate(g_til)
    deallocate(s)
    deallocate(g)
    deallocate(gam)
    deallocate(hh)

end subroutine harten_tvd_roeAvg

subroutine harten_tvd_avg(order)

    use shared
    implicit none


    integer(kind=4) :: gs, i, order

    real(kind=8) :: dtdx, signal

    real(kind=8) :: rhol,rhor,ul,ur,pl,pr,hl,hr,rho_hat,c1,c2

    real(kind=8), allocatable, dimension(:) :: u_hat
    real(kind=8), allocatable, dimension(:) :: h_hat
    real(kind=8), allocatable, dimension(:) :: c_hat

    real(kind=8), allocatable, dimension(:,:) :: v_k
    real(kind=8), allocatable, dimension(:,:,:) :: r
    real(kind=8), allocatable, dimension(:,:) :: alpha
    real(kind=8), allocatable, dimension(:,:) :: eps
    real(kind=8), allocatable, dimension(:,:) :: phi
    real(kind=8), allocatable, dimension(:,:) :: g_til
    real(kind=8), allocatable, dimension(:,:) :: s
    real(kind=8), allocatable, dimension(:,:) :: g
    real(kind=8), allocatable, dimension(:,:) :: gam
    real(kind=8), allocatable, dimension(:,:) :: hh

    ! Allocating Roe's avgs.
    allocate(u_hat(total_mesh_points + ng))
    allocate(h_hat(total_mesh_points + ng))
    allocate(c_hat(total_mesh_points + ng))

    ! Allocating auxiliary arrays for hartens scheme.
    allocate(v_k(3,total_mesh_points + ng))
    allocate(r(3,3,total_mesh_points + ng))
    allocate(alpha(3,total_mesh_points + ng))
    allocate(eps(3,total_mesh_points + ng))
    allocate(phi(3,total_mesh_points + ng))
    allocate(g_til(3,total_mesh_points + ng))
    allocate(s(3,total_mesh_points + ng))
    allocate(g(3,total_mesh_points + ng))
    allocate(gam(3,total_mesh_points + ng))
    allocate(hh(3,total_mesh_points + ng))


    gs = ng/2
    dtdx = time_step/dx


    ! Lets first define the Roe avgs...

    do i = 1, total_mesh_points + ng-1

       rhol = q(i,1)
       rhor = q(i+1,1)

       ul = q(i,2) / q(i,1)
       ur = q(i+1,2) / q(i+1,1)

       pl = press(i)
       pr = press(i+1)

       hl = (q(i,3) + pl) / rhol
       hr = (q(i+1,3) + pr) / rhor


       ! Build avgs..

       rho_hat = (rhol + rhor) / 2.0d0

       u_hat(i) = (ul + ur) / 2.0d0
       h_hat(i) = (hl + hr) / 2.0d0
       c_hat(i) = (a_speed(i) + a_speed(i+1))/2.0d0

    end do


    ! Now lets do the v_k.

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        v_k(1,i) = dtdx*(u_hat(i) - c_hat(i))
        v_k(2,i) = dtdx*(     u_hat(i)      )
        v_k(3,i) = dtdx*(u_hat(i) + c_hat(i))

    end do
    !$omp end do
    !$omp end parallel


    ! Now the R's.

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        r(1,1,i) = 1.0d0
        r(2,1,i) = u_hat(i) - c_hat(i)
        r(3,1,i) = h_hat(i) - u_hat(i)*c_hat(i)

        r(1,2,i) = 1.0d0
        r(2,2,i) = u_hat(i)
        r(3,2,i) = 0.5d0*u_hat(i)**2.0d0

        r(1,3,i) = 1.0d0
        r(2,3,i) = u_hat(i) + c_hat(i)
        r(3,3,i) = h_hat(i) + u_hat(i)*c_hat(i)

    end do
    !$omp end do
    !$omp end parallel


    ! At this part of the code I would calculate the eps but since I am lazy
    ! right now, I ll just make then fixed.

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        eps(1,i) = 0.005d0
        eps(2,i) = 0.0d0
        eps(3,i) = 0.0d0

    end do
    !$omp end do
    !$omp end parallel


    ! Now the alphas.

    do i = 1, total_mesh_points + ng-1

       c1 = ((fgamma-1.0d0)/(c_hat(i)**2.0d0)) * ( jump(q(i+1,3),q(i,3)) + 0.5d0 * u_hat(i)*jump(q(i+1,1),q(i,1)) - u_hat(i)*jump(q(i+1,2),q(i,2)) )
       c2 = (1.0d0/c_hat(i)) * ( jump(q(i+1,2),q(i,2)) - u_hat(i)*jump(q(i+1,1),q(i,1)) )

       alpha(1,i) = 0.5d0*( c1 - c2 )
       alpha(2,i) = jump(q(i+1,1),q(i,1)) - c1
       alpha(3,i) = 0.5d0*( c1 + c2 )

    end do


    ! For the g_til calculation we need our phis, take care with this vector
    ! because we will over write it later...

    phi = 0.0d0

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        phi(1,i) = psiHarten(v_k(1,i), eps(1,i))
        phi(2,i) = psiHarten(v_k(2,i), eps(2,i))
        phi(3,i) = psiHarten(v_k(3,i), eps(3,i))

    end do
    !$omp end do
    !$omp end parallel


    ! And finally lets calculate the g_tils

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

       g_til(1,i) = 0.5d0*(phi(1,i) - v_k(1,i)**2.0d0)*alpha(1,i) 
       g_til(2,i) = 0.5d0*(phi(2,i) - v_k(2,i)**2.0d0)*alpha(2,i) 
       g_til(3,i) = 0.5d0*(phi(3,i) - v_k(3,i)**2.0d0)*alpha(3,i) 

    end do
    !$omp end do
    !$omp end parallel


    ! Now the s's... 

    signal = 1.0d0

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        s(1,i) = sign(signal,g_til(1,i))
        s(2,i) = sign(signal,g_til(2,i))
        s(3,i) = sign(signal,g_til(3,i))

    end do
    !$omp end do
    !$omp end parallel


    ! Finally the g's !!

    do i = 2, total_mesh_points + ng-1

        if (order == 1) then

            g(1,i) = 0.0d0
            g(2,i) = 0.0d0
            g(3,i) = 0.0d0

        else if (order == 2) then

            g(1,i) = s(1,i) * max(0.0d0,min(abs(g_til(1,i)),g_til(1,i-1)*s(1,i)))
            g(2,i) = s(2,i) * max(0.0d0,min(abs(g_til(2,i)),g_til(2,i-1)*s(2,i)))
            g(3,i) = s(3,i) * max(0.0d0,min(abs(g_til(3,i)),g_til(3,i-1)*s(3,i)))

        end if

    end do


    ! Lets now calculate the gammas.

    gam = 0.0d0

    do i = 1, total_mesh_points + ng-1

            gam(1,i) = gamHarten(g(1,i), g(1,i+1), alpha(1,i))
            gam(2,i) = gamHarten(g(2,i), g(2,i+1), alpha(2,i))
            gam(3,i) = gamHarten(g(3,i), g(3,i+1), alpha(3,i))

    end do


    ! As I said, we will now re-build our phi vector to receive as arguments the
    ! v_k and gam summed... lets do it.

    phi = 0.0d0

    !$omp parallel
    !$omp do
    do i = 1, total_mesh_points + ng-1

        phi(1,i) = psiHarten(v_k(1,i) + gam(1,i), eps(1,i))
        phi(2,i) = psiHarten(v_k(2,i) + gam(2,i), eps(2,i))
        phi(3,i) = psiHarten(v_k(3,i) + gam(3,i), eps(3,i))

    end do
    !$omp end do
    !$omp end parallel

    
    ! Lets now calculate our harten "dissipation". Note that in first order,
    ! the values of g will be zero.

    do i = 1, total_mesh_points + ng-1


        hh(1,i) = r(1,1,i)*( g(1,i) + g(1,i+1) - (phi(1,i)*alpha(1,i)) ) + &
                  r(1,2,i)*( g(2,i) + g(2,i+1) - (phi(2,i)*alpha(2,i)) ) + &
                  r(1,3,i)*( g(3,i) + g(3,i+1) - (phi(3,i)*alpha(3,i)) )

        hh(2,i) = r(2,1,i)*( g(1,i) + g(1,i+1) - (phi(1,i)*alpha(1,i)) ) + &
                  r(2,2,i)*( g(2,i) + g(2,i+1) - (phi(2,i)*alpha(2,i)) ) + &
                  r(2,3,i)*( g(3,i) + g(3,i+1) - (phi(3,i)*alpha(3,i)) )

        hh(3,i) = r(3,1,i)*( g(1,i) + g(1,i+1) - (phi(1,i)*alpha(1,i)) ) + &
                  r(3,2,i)*( g(2,i) + g(2,i+1) - (phi(2,i)*alpha(2,i)) ) + &
                  r(3,3,i)*( g(3,i) + g(3,i+1) - (phi(3,i)*alpha(3,i)) )

    end do


    ! We are almost done now.

    do i = 1, total_mesh_points + ng-1

        f_pls(i,1) = 0.5d0*( f(i,1) + f(i+1,1) + (1.0d0/dtdx)*hh(1,i) )
        f_pls(i,2) = 0.5d0*( f(i,2) + f(i+1,2) + (1.0d0/dtdx)*hh(2,i) )
        f_pls(i,3) = 0.5d0*( f(i,3) + f(i+1,3) + (1.0d0/dtdx)*hh(3,i) )

    end do

    deallocate(u_hat)
    deallocate(h_hat)
    deallocate(c_hat)
    deallocate(v_k)
    deallocate(r)
    deallocate(alpha)
    deallocate(eps)
    deallocate(phi)
    deallocate(g_til)
    deallocate(s)
    deallocate(g)
    deallocate(gam)
    deallocate(hh)

end subroutine  harten_tvd_avg


subroutine harten_explicit_1stOrder

    use shared
    implicit none

    integer(kind=4) :: i, gs
    integer(kind=4) :: neq = 3

    real(kind=8), allocatable, dimension(:,:) :: q_aux

    gs = ng/2  !Ghosts in each side

    allocate(q_aux(total_mesh_points + ng,neq))

    q_aux(:,:) = q(:,:)

    do i = gs, total_mesh_points + gs

        q(i,1) = q_aux(i,1) - (time_step/dx) * ( f_pls(i,1) - f_pls(i-1,1) ) 
        q(i,2) = q_aux(i,2) - (time_step/dx) * ( f_pls(i,2) - f_pls(i-1,2) ) 
        q(i,3) = q_aux(i,3) - (time_step/dx) * ( f_pls(i,3) - f_pls(i-1,3) ) 

    end do


    deallocate(q_aux)

end subroutine  harten_explicit_1stOrder

