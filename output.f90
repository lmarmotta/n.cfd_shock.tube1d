subroutine output

    use shared 
    implicit none

    integer(kind=4) :: i,gs
    real(kind=8) :: p, u

    gs = ng/2  

    write(*,*) ""
    write(*,*) "+--------------------------------------------------------------+"
    write(*,*) "|               +++ Outputting solution +++                    |"
    write(*,*) "+--------------------------------------------------------------+"
    write(*,*) ""

    open(2,file='pressureOutput.out')
    open(3,file='densityOutput.out')
    open(4,file='velocityOutput.out')


    write(2,*) "mesh_points , rho , rho*u , energy" 


    do i = 1, total_mesh_points-1

        p    = (fgamma-1.0d0)*( q(i,3) - 0.5d0 * (q(i,2)**2/q(i,1)) )
        u = q(i,2)/q(i,1)

        write(2,'(5(f10.3))') mesh_points(i), p
        write(3,'(5(f10.3))') mesh_points(i), q(i,1)
        write(4,'(5(f10.3))') mesh_points(i), u

    end do


    close(2) 
    close(3) 
    close(4) 

end subroutine output

subroutine print_iter(actual)

    use shared
    implicit none

    integer(kind=4) :: actual
    real(kind=8) :: max_domain
    real(kind=8) :: now_time


    max_domain = maxval(q(:,1))

    now_time = actual*time_step


    if (actual == 1) then
        pux = print_step
    end if
    
    if (actual == 1) then

        write(*,'(A,I9,A,F10.4,A,F10.4,A,F10.4,A)') "| iter: ", actual, " | Max rho: ", max_domain, &
            " | Sol time: ", now_time, " [s]| Time per it: ", time_per_iter, " [s]"

    else if (actual == print_step) then

        write(*,'(A,I9,A,F10.4,A,F10.4,A,F10.4,A)') "| iter: ", actual, " | Max rho: ", max_domain, &
            " | Sol time: ", now_time, " [s]| Time per it: ", time_per_iter, " [s]"

        print_step = print_step + pux

    end if


end subroutine print_iter


subroutine stubex(p1,p4,r1,r4,u1,u4,tm,xl,xr,xd)

    implicit none

!     Computes the exact solution for the Sod's shock-tube problem.
!     See J.D. Anderson, Modern Compressible Flow (1984) for details.

!     ORIGIN: https://www.grc.nasa.gov/WWW/wind/valid/stube/stube.html
!     ADAPTED BY: Leonardo Motta Maia


      ! Input variables
      real(kind=8) :: gam,gm1,gp1,p4,p1,r4,r1,u1,tol,tm,xl,xr,xd

      ! Inside variables.
      integer(kind=4) :: nxp, n
      real(kind=8) :: rm,rx,px,ux,xx,xtail,xhead,xc,rmach4,rmach3,rmach2,rmach1
      real(kind=8) :: u4,u3,u2,up,a3,r3,p3,r2,a1,a4,p2p1,iterr,t2t1,r2r1,wsp
      real(kind=8) :: xs, p2



!c...Set constants.

      gam = 1.4
      gm1 = gam - 1.0
      gp1 = gam + 1.0

!c...Set initial states (non-dimensional).

      tol = 1.0E-05

!
!c...Compute acoustic velocities.

      a1 = sqrt( gam * p1 / r1 )
      a4 = sqrt( gam * p4 / r4 )

!c...Use a Newton-secant iteration to compute p2p1.

      call  sp2p1 ( gam, p1, a1, p4, a4, p2p1, iterr, tol )

!c...t2t1.

      t2t1 = p2p1 * ( gp1/gm1 + p2p1 ) / ( 1.0 + gp1 * p2p1 / gm1 )

!c...r2r1.

      r2r1 = ( 1.0 + gp1 * p2p1 / gm1 ) / ( gp1 / gm1 + p2p1 )

!c...W, shock-wave speed.

      wsp  = a1 * sqrt( gp1 * ( p2p1 - 1.0 ) / ( 2.0 * gam ) + 1.0 )

!c...Shock location.

      xs = xd + wsp * tm

!c...State 2.

      p2 = p2p1 * p1
      r2 = r2r1 * r1

!c...State 3.

      p3 = p2

!c...Isentropic between 3 and 4.

      r3 = r4 * ( p3 / p4 )**(1.0/gam)

      a3 = sqrt( gam * p3 / r3 )

!c...Speed of contact discontinuity.

      up = 2.0 * a4 * ( 1.0 - (p2/p4)**(0.5*gm1/gam) ) / gm1

      u2 = up
      u3 = up

!c...Mach numbers.

      rmach1 = u1 / sqrt( gam * p1 / r1 )
      rmach2 = u2 / sqrt( gam * p2 / r2 )
      rmach3 = u3 / sqrt( gam * p3 / r3 )
      rmach4 = u4 / sqrt( gam * p4 / r4 )

!c...Location of contact discontinuity.

      xc = xd + up * tm

!c...Location of expansion region.

      xhead = xd + ( u4 - a4 ) * tm
      xtail = xd + ( u3 - a3 ) * tm


      open ( unit= 7, file='apressure.out' )
      open ( unit= 8, file='adensity.out' )
      open ( unit= 9, file='aspeed.out' )
      open ( unit=10, file='amach.out' )

      write (  7, 100 ) xl, p4
      write (  8, 100 ) xl, r4
      write (  9, 100 ) xl, u4
      write ( 10, 100 ) xl, rmach4

      write (  7, 100 ) xhead, p4
      write (  8, 100 ) xhead, r4
      write (  9, 100 ) xhead, u4
      write ( 10, 100 ) xhead, rmach4

      nxp = 11

      do  n = 1, 11
        xx = xhead + ( xtail - xhead )  * n / ( nxp + 1.0 )
        ux = u4 + u3 * ( xx - xhead ) / ( xtail - xhead )
        px = p4 * ( 1.0 - 0.5 * gm1 * ( ux / a4 ) )**( 2.0 * gam / gm1 )
        rx = r4 * ( 1.0 - 0.5 * gm1 * ( ux / a4 ) )**( 2.0 / gm1 )
        rm = ux / sqrt( gam * px / rx )
        write ( 7, 100 ) xx, px
        write ( 8, 100 ) xx, rx
        write ( 9, 100 ) xx, ux
        write (10, 100 ) xx, rm
      enddo

      write (  7, 100 ) xtail, p3
      write (  8, 100 ) xtail, r3
      write (  9, 100 ) xtail, u3
      write ( 10, 100 ) xtail, rmach3

      write (  7, 100 ) xc, p3
      write (  8, 100 ) xc, r3
      write (  9, 100 ) xc, u3
      write ( 10, 100 ) xc, rmach3

      write (  7, 100 ) xc, p2
      write (  8, 100 ) xc, r2
      write (  9, 100 ) xc, u2
      write ( 10, 100 ) xc, rmach2

      write (  7, 100 ) xs, p2
      write (  8, 100 ) xs, r2
      write (  9, 100 ) xs, u2
      write ( 10, 100 ) xs, rmach2

      write (  7, 100 ) xs, p1
      write (  8, 100 ) xs, r1
      write (  9, 100 ) xs, u1
      write ( 10, 100 ) xs, rmach1

      write (  7, 100 ) xr, p1
      write (  8, 100 ) xr, r1
      write (  9, 100 ) xr, u1
      write ( 10, 100 ) xr, rmach1

!c...Format statements.

  100 format ( 3x, 2(2x,f12.6) )

!c-----------------------------------------------------------------------

      
end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sp2p1(gam,p1,a1,p4,a4,p2p1,iterr,tol) 

!c        Uses Newton-secant method to iterate on eqn 7.94 (Anderson, 
!c        1984) to fine p2p1 across moving shock wave.
!c
!c-----------------------------------------------------------------------
!

      implicit none

      ! Input variables.
      real(kind=8) :: gam, p1, a1, p4, a4, p2p1, iterr, tol 

      ! Insiders.
      real(kind=8) :: p2p1n,f,gm1,gp1,p2p1m,t1,t2,t3,fm
      integer(kind=4) :: iter,itmax

!c...Set some variables

      gm1 = gam - 1.0
      gp1 = gam + 1.0

!c...Initialize p2p1 for starting guess

      p2p1m = 0.25 * p4 / p1

      t1 = - 2.0 * gam / gm1

      t2 = gm1 * ( a1 / a4 ) * ( p2p1m - 1.0 )
      t3 = 2.0 * gam * ( 2.0 * gam + gp1 * ( p2p1m - 1.0 ) )
      fm = p4 / p1 - p2p1m * ( 1.0 - t2 / sqrt(t3) )**t1

!c...Perturb p2p1

      p2p1 = 0.95 * p2p1m

!c...Begin iteration

      iter  = 0
      itmax = 20

   10 continue

        iter = iter + 1

        t2 = gm1 * ( a1 / a4 ) * ( p2p1 - 1.0 )
        t3 = 2.0 * gam * ( 2.0 * gam + gp1 * ( p2p1 - 1.0 ) )

        f  = p4 / p1 - p2p1 * ( 1.0 - t2 / sqrt(t3) )**t1 

!        write (*,*) 'iter, p2p1, f: ', iter, p2p1, f

        if ( abs(f) .gt. tol  .and.  iter .lt. itmax ) then 
          p2p1n = p2p1 - f * ( p2p1 - p2p1m ) / ( f - fm )
          p2p1m = p2p1
          fm    = f
          p2p1  = p2p1n
          go to 10
        endif

!c...Check to see if maximum iterations reached

      iterr = 0
      if ( iter .lt. itmax )  iterr = 1

!c-----------------------------------------------------------------------

      return
end subroutine sp2p1 
