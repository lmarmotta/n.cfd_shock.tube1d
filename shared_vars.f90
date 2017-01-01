module shared        

    implicit none


!
! ** NAMELIST DECLARATIONS **
!
    ! Pressure in the start & end.
    real(kind=8) :: p1, p4             

    ! Gas constants.
    real(kind=8) :: fgamma, R_const    

    ! Densities.
    real(kind=8) :: rho1,rho4          

    ! More properties.
    real(kind=8) :: F_Cp, F_Cv         

    ! Dt
    real(kind=8) :: time_step, time_per_iter

    ! Number of time advances.
    integer(kind=4) :: iterations      

    ! Type of scheme.
    integer(kind=4) :: scheme          

    ! BWH artificial dissp. term.
    real(kind=8) :: dissip_omega       

    ! Artificial dissipation scheme.
    integer(kind=4) :: dissp_scheme    

    ! Print each x iterations.
    integer(kind=4) :: print_step,pux  

    ! Timer variable.
    real(kind=8) :: elapse_time


!
! ** MESH GEOMETRY INPUT **
!
    ! Number of ghosts.
    integer(kind=4) :: ng = 4

    ! Total internal points in the mesh.
    integer(kind=4) :: total_mesh_points

    ! Start point coordnate.
    real(kind=8) :: start_mesh_point

    ! Final point coordnate.
    real(kind=8) :: final_mesh_point


!
! ** GEOMETRICAL DEFINITIONS **
!
    ! Distance between two consecutive points (constant for this case)
    real(kind=8) :: dx
    
    ! THE mesh points.
    real(kind=8), allocatable, dimension(:) :: mesh_points


!
! ** VECTORS OF PROPERTIES **
!
    ! q(i,1) = rho
    ! q(i,2) = rho*u
    ! q(i,3) = energy
    real(kind=8), allocatable, dimension(:,:)   :: q

    ! f(i,1) = rho*u
    ! f(i,2) = rho*u^2 + p
    ! f(i,3) = rho*u(e+p/rho)
    real(kind=8), allocatable, dimension(:,:)   :: f

    ! Jacobian
    real(kind=8), allocatable, dimension(:,:,:) :: a_j

    ! MacCormack predition step.
    real(kind=8), allocatable, dimension(:,:)   :: q_pred
    real(kind=8), allocatable, dimension(:,:)   :: f_pred

    ! Pressure.
    real(kind=8), allocatable, dimension(:)     :: press

    ! Artificial dissipation vector.
    real(kind=8), allocatable, dimension(:,:)   :: art_dissip
    real(kind=8), allocatable, dimension(:,:)   :: art_dissip2d
    real(kind=8), allocatable, dimension(:,:)   :: art_dissip4d

    ! JSP dissipation.
    real(kind=8), allocatable, dimension(:,:)   :: jst_d

    ! Speed of sound stuff.
    real(kind=8), allocatable, dimension(:)     :: a_speed
    real(kind=8), allocatable, dimension(:)     :: mach

    ! Splited fluxes for prj.2
    real(kind=8), allocatable, dimension(:,:) :: f_pls 
    real(kind=8), allocatable, dimension(:,:) :: f_min

    ! Transformation matrix
    real(kind=8), allocatable, dimension(:,:,:) :: p_matrix
    real(kind=8), allocatable, dimension(:,:,:) :: p_matrix_inv



!
! ** PHYSICAL CALCULATED VARIABLES **
!
    ! Energy at the start of the tube.
    real(kind=8) :: e1 = 0.0d0

    ! Energy at the end of the tube.
    real(kind=8) :: e4 = 0.0d0


    contains

        !
        ! ** OPERATORS TO MAKE THINGS EASIER **
        !

        ! Forward 1st order differential operator.
        real(kind=8) function d_fwrd_1st(value_p0,value_p1,dx)

            implicit none

            ! Value at p(i+0)
            real(kind=8) :: value_p0

            ! Value at p(i+1)
            real(kind=8) :: value_p1

            ! dx (you know this one man !)
            real(kind=8) :: dx

            d_fwrd_1st = (value_p1 - value_p0) / dx

        end function  d_fwrd_1st


        ! Backward 1st order differential operator.
        real(kind=8) function d_bcwd_1st(value_p0,value_m1,dx)

            implicit none

            ! Value at p(i+0)
            real(kind=8) :: value_p0

            ! Value at p(i+1)
            real(kind=8) :: value_m1

            ! dx (you know this one man !)
            real(kind=8) :: dx

             d_bcwd_1st = (value_p0 - value_m1) / dx

        end function  d_bcwd_1st


        ! Difference operator proposed by S&W.
        ! Paper Steger-Warming 1981 Eq. 5.14

        real(kind=8) function diff_bcwd(value_p0,value_m1,value_m2,dx)

            implicit none

            ! Value at p(i+0)
            real(kind=8) :: value_p0

            ! Value at p(i-1)
            real(kind=8) :: value_m1

            ! Value at p(i-2)
            real(kind=8) :: value_m2

            ! dx (you know this one man !)
            real(kind=8) :: dx

            diff_bcwd = (3.0d0*value_p0 - 4.0d0*value_m1 + value_m2)/(2.0d0*dx)

        end function diff_bcwd


        ! Difference operator proposed by S&W.
        ! Paper Steger-Warming 1981 Eq. 5.14

        real(kind=8) function diff_frwd(value_p0,value_m1,value_m2,dx)

            implicit none

            ! Value at p(i+0)
            real(kind=8) :: value_p0

            ! Value at p(i-1)
            real(kind=8) :: value_m1

            ! Value at p(i-2)
            real(kind=8) :: value_m2

            ! dx (you know this one man !)
            real(kind=8) :: dx

            diff_frwd = (-3.0d0*value_p0 + 4.0d0*value_m1 - value_m2)/(2.0d0*dx)

        end function diff_frwd 


        ! Define the jump opertaor for harten scheme.

        real(kind=8) function jump(prop_p1, prop)

            implicit none

            ! Propertie in the advanced point.
            real(kind=8) :: prop_p1

            ! Propertie in the current point.
            real(kind=8) :: prop

            jump = prop_p1 - prop

        end function jump


        ! For Harten's TVD scheme.

        real(kind=8) function psiHarten(z, eps)

            implicit none

            real(kind=8) :: psi
            real(kind=8) :: z
            real(kind=8) :: eps

            psi = dabs(z)

            if (psi < eps) then 
                psi = (z**2.0d0 + eps**2.0d0) / (2.0d0 * eps)
            end if

            psiHarten = psi

        end function psiHarten


        ! This harten is very demanding man !

        real(kind=8) function gamHarten(g, gp1, alpha)

            implicit none

            real(kind=8) :: g, gp1, alpha

            if (alpha /= 0.0d0) then
                gamHarten = (gp1 - g) / alpha
            else
                gamHarten = 0.0d0
            end if

        end function gamHarten

end module shared


