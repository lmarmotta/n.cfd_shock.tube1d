! +---------------------------------------------------------------------------+! 
! |                                                                           |!
! |          SOD SHOCK TUBE - EULER EQUATIONS: CC299 PROJECT 01               |!
! |                                                                           |!
! +---------------------------------------------------------------------------+!
program euler

    use shared
    implicit none

    integer(kind=4) :: t = 0
    integer(kind=8) :: c1,c2,cr
    integer(kind=8) :: c11,c22,crr
    real(kind=8) :: rate, ratee


    call indat
    
    !
    ! ** ALLOCATE VARIABLES **
    ! 
    
    ! Points in the mesh in total (including ghosts).
    allocate(mesh_points(total_mesh_points + ng))

    ! Primitive vector of variables.
    allocate(q(total_mesh_points + ng,3))

    ! Vector of fluxes.
    allocate(f(total_mesh_points + ng,3))

    ! Jacobian (physyical one).
    allocate(a_j(3,3,total_mesh_points + ng))

    ! For MacCormak method, the cloned predictive step vectors.
    allocate(q_pred(total_mesh_points + ng,3))
    allocate(f_pred(total_mesh_points + ng,3))

    ! Pressure vector (mainly used for JST dissipation).
    allocate(press(total_mesh_points + ng))

    ! JST artificial dissipation vector.
    allocate(jst_d(total_mesh_points + ng,3))

    ! General artificial dissipation vector (all dissip schemes).
    allocate(art_dissip(total_mesh_points + ng,3))

    ! Fourth and second difference artificial dissipation (implicit schemes).
    allocate(art_dissip2d(total_mesh_points + ng,3))
    allocate(art_dissip4d(total_mesh_points + ng,3))

    ! Speed of sound for our flux based schemes.
    allocate(a_speed(total_mesh_points + ng))
    allocate(mach(total_mesh_points + ng))

    ! Splited fluxes for the FVS methods.
    allocate(f_pls(total_mesh_points + ng, 3))
    allocate(f_min(total_mesh_points + ng, 3))

    ! Transformation matrices
    allocate(p_matrix(3,3,total_mesh_points + ng))
    allocate(p_matrix_inv(3,3,total_mesh_points + ng))

    !
    ! ** NOW THE PARTY BEGINS **
    !


    ! Creating initial conditions for the case.

    call calcInitialProperties
    call createMeshPoints
    call createInitialCondition


    ! Print out some useful information about the configuration.

    write(*,*) "" 
    if (scheme == 1) then
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Explicit Euler (t) centered (space)."
    else if (scheme == 2) then
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Lax-Wendroff. "
    else if (scheme == 3) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with MacCormak."
    else if (scheme == 4) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Implicit Beam-Warming."
    else if (scheme == 5) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with 1st order Steger-Warming scheme."
    else if (scheme == 6) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with 2st order Steger-Warming scheme."
    else if (scheme == 7) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with 1st order Van-Leer scheme."
    else if (scheme == 8) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with 2nd order Van-Leer scheme."
    else if (scheme == 9) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with 1st order Roe scheme."
    else if (scheme == 10) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with AUSM scheme."
    else if (scheme == 11) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Implicit Steger-Warming scheme."
    else if (scheme == 12) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Implicit VanLeer scheme."
    else if (scheme == 13) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Harten TVD 1st order scheme."
    else if (scheme == 14) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Harten TVD 2nd order scheme."
    else if (scheme == 15) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Harten TVD 1st order scheme."
    else if (scheme == 16) then 
        write(*,'(A,I7,A)') " Running: ", iterations, " iterations with Harten TVD 2nd order scheme."
    end if
    write(*,*) "" 
    

    ! Calling analytical solutin

    call stubex(p4,p1,rho4,rho1,0.0d0,0.0d0, &
        iterations*time_step,start_mesh_point,final_mesh_point,0.0d0)


    !  Start main time step loop.

    elapse_time = 0.0d0
    time_per_iter = 0.0d0

    ! Time measurements stuff.

    call system_clock(count_rate=cr)
    call system_clock(count_rate=crr)

    rate = REAL(cr)
    ratee = REAL(crr)

    call system_clock(c1)


    !
    !    +++ MAIN ITERATION LOOP +++
    !

    do t = 1, iterations


        if (scheme == 1) then

            call system_clock(c11)

            call fluxCalculation
            call artificialDissip
            call explicitCentredSimple

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 2) then

            call system_clock(c11)

            call fluxCalculation
            call calcJacobian
            call artificialDissip
            call laxWendroff

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 3) then

            call system_clock(c11)

            call fluxCalculation
            call artificialDissip
            call maccormack_pred
            call artificialDissip
            call maccormack_corr

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 4) then

            call system_clock(c11)

            call fluxCalculation
            call calcJacobian
            call artificialDissip
            call implicitBeamWarming

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 5) then

            call system_clock(c11)

            call sound_speed
            call artificialDissip
            call steger_warming
            call stegerWarming_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 6) then

            call system_clock(c11)

            call sound_speed
            call artificialDissip
            call steger_warming
            call stegerWarming_explicit_2stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 7) then

            call system_clock(c11)

            call sound_speed
            call machcalc
            call vanleer
            call vanleer_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 8) then

            call system_clock(c11)

            call sound_speed
            call machcalc
            call vanleer
            call vanleer_explicit_2stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 9) then

            call system_clock(c11)

            call calc_pressure
            call fluxCalculation
            call roe_flux
            call roe_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 10) then

            call system_clock(c11)

            call ausm_plus
            call ausm_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 11) then

            call system_clock(c11)

            call sound_speed
            call steger_warming
            call implicit_fvsfds_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 12) then

            call system_clock(c11)

            call sound_speed
            call machcalc
            call vanleer
            call implicit_fvsfds_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 13) then

            call system_clock(c11)

            call calc_pressure
            call fluxCalculation
            call harten_tvd_roeAvg(1)
            call harten_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 14) then

            call system_clock(c11)

            call calc_pressure
            call fluxCalculation
            call harten_tvd_roeAvg(2)
            call harten_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 15) then

            call system_clock(c11)

            call calc_pressure
            call fluxCalculation
            call sound_speed
            call harten_tvd_avg(1)
            call harten_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        else if (scheme == 16) then

            call system_clock(c11)

            call calc_pressure
            call fluxCalculation
            call sound_speed
            call harten_tvd_avg(2)
            call harten_explicit_1stOrder

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/ratee

        end if 

        call print_iter(t)

    end do

    call system_clock(c2)

    elapse_time = real(c2 - c1,kind=8)/rate


    write(*,*) ""
    write(*,'(A,F10.6,A)') "| Total elapse time: ",elapse_time," [s]| "


    call output


    ! Deallocate variables.

    deallocate(mesh_points)
    deallocate(     q     )
    deallocate(     f     )
    deallocate(    a_j    )
    deallocate(  q_pred   )
    deallocate(  f_pred   )
    deallocate(   press   )
    deallocate(   jst_d   )
    deallocate(art_dissip )

    deallocate(art_dissip2d)
    deallocate(art_dissip4d)

    deallocate(  a_speed  )
    deallocate(   mach    )

    deallocate(p_matrix)
    deallocate(p_matrix_inv)

end program euler


