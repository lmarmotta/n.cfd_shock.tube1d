subroutine implicitBeamWarming

    ! This subroutine prepares advance the primitive variables vectors q(i,eq)
    ! in time using implicit Beam-Warming scheme. The cc299tridiag solver is c
    ! alled every time step to solve the system and bring back the n+1 solutio
    ! n. If a sigular matrix is created, the solver will stop and youll go back 
    ! to gdb.

    use shared
    implicit none


    ! Loop variables.
    integer(kind=4) :: i,ii,jj,gs

    ! Our scalar dx/dt
    real(kind=8) :: dtdx

    ! Product between primitives and jacobian.
    real(kind=8), allocatable, dimension(:,:) :: a_jxq

    ! Identity matrix.
    real(kind=8), dimension(3,3) :: identity

    ! RHS side of the system (Ax=B | B = RHS)
    real(kind=8), allocatable, dimension(:,:) :: RHS

    ! Main diagonal of our system. 
    real(kind=8), allocatable, dimension(:,:,:) :: maind

    ! Lower diagonal of our system. 
    real(kind=8), allocatable, dimension(:,:,:) :: lowerd

    ! Upper diagonal of our system. 
    real(kind=8), allocatable, dimension(:,:,:) :: upperd

    ! Advanced in time solution.
    real(kind=8), allocatable, dimension(:,:) :: unp1


    ! Ghosts in each side.

    gs = ng/2  


    ! The dt/dx constant (making things cleaner.

    dtdx = time_step/dx


    ! We will need an identity later in this episode.

    call identityCreate(3, identity)


    ! Allocating the product between the jacobian (3x3) and the primitive varia
    ! bles q(i,eq) (3x1) resulting in a_jxq (3x1).

    allocate(a_jxq(3, total_mesh_points + ng))


    ! Lets now allocate our RHS following the calling index ordering, in this 
    ! case we will need just internal points so.

    allocate(RHS(3, total_mesh_points))


    ! Lets now allocate our main matrix, this matrix has the same size of the nu
    ! mber of internal points.

    allocate(maind(3, 3, total_mesh_points))


    ! Allocate lower diagonal.

    allocate(lowerd(3, 3, total_mesh_points))


    ! Allocate the upper diagonal.

    allocate(upperd(3, 3, total_mesh_points))


    ! Allocate the n+1 vector.

    allocate(unp1(3,total_mesh_points))


    ! We will use a lot the product between the jacobian matrix a_j and the pr
    ! imitive vectors, so now we will create a vector a_jxq(ii,i) for every
    ! point in the mesh so things will ge easier. Lets loop through all.

    a_jxq = 0.0d0

    do i = 1, total_mesh_points 

        a_jxq(:,i+gs) = matmul(a_j(:,:,i+gs),q(i+gs,:))

    end do


    ! Lets now create our RHS(B) for the internal points. Note that although the
    ! RHS is receiving the correct index for the formation of the matrix, our
    ! properties variables are not !  So we need to perform some basic index ar
    ! itmetics on primitives and flux vectors seeking correct index allignement.

    RHS = 0.0d0

    do i = 1, total_mesh_points

        RHS(1,i) = q(i+gs,1) - 0.5d0*dtdx*( f(i+1+gs,1) - f(i-1+gs,1) ) + 0.25d0*dtdx*(a_jxq(1,i+1+gs) - a_jxq(1,i-1+gs)) + art_dissip(i+gs,1)/8.0d0
        RHS(2,i) = q(i+gs,2) - 0.5d0*dtdx*( f(i+1+gs,2) - f(i-1+gs,2) ) + 0.25d0*dtdx*(a_jxq(2,i+1+gs) - a_jxq(2,i-1+gs)) + art_dissip(i+gs,2)/8.0d0
        RHS(3,i) = q(i+gs,3) - 0.5d0*dtdx*( f(i+1+gs,3) - f(i-1+gs,3) ) + 0.25d0*dtdx*(a_jxq(3,i+1+gs) - a_jxq(3,i-1+gs)) + art_dissip(i+gs,3)/8.0d0

    end do


    ! With our RHS created, we can now proceed to create our lower main and
    ! upper diagonals. Lets start with the main one.

    do i = 1, total_mesh_points

        do ii = 1, 3
        do jj = 1, 3
        
            maind(ii,jj,i) = identity(ii,jj) 

        end do
        end do
        
    end do


    ! Now lets create the lower diagonal of the system, always paying atention 
    ! to our indexes alligment. Note that the lowerd vector is allocated in the
    ! same size as the maind matrix. This is required by the tridiagonal solver.

    do i = 1, total_mesh_points

        do ii = 1, 3
        do jj = 1, 3

            lowerd(ii,jj,i) = - ( 0.25d0*dtdx*a_j(ii,jj,i-1+gs) )

        end do
        end do
            
    end do


    ! Following the style of the lower diagonal lets create our upper diagonal.

    do i = 1, total_mesh_points

        do ii = 1, 3
        do jj = 1, 3

            upperd(ii,jj,i) = ( 0.25d0*dtdx*a_j(ii,jj,i+1+gs) ) 

        end do
        end do
            
    end do


    ! Ok folks... now we have all the vectors needed by our system solver, so,
    ! lets call it !. Note that there are two different solvers in this code.
    ! One of the (...)_opt uses matmul for the matrix multiplication, the othe
    ! er one makes extensive use of lapacks subroutines. The results show that
    ! lapack is slower so... lets use the simpler matmul.

    unp1 = 0.0d0   ! Cleanning the answer vector.

    !call cc299blktriad(maind,lowerd,upperd,3,total_mesh_points,RHS,unp1)
    call cc299blktriad_opt(maind,lowerd,upperd,3,total_mesh_points,RHS,unp1)


    ! Now we can give update our solution primitive vector with the brand new
    ! solution.

    do i = 1, total_mesh_points
        q(i+gs,1) = unp1(1,i)
        q(i+gs,2) = unp1(2,i)
        q(i+gs,3) = unp1(3,i)
    end do


    deallocate(  a_jxq  )
    deallocate(   RHS   )
    deallocate(  maind  )
    deallocate( lowerd  )
    deallocate( upperd  )
    deallocate(  unp1   )

end subroutine implicitBeamWarming

subroutine identityCreate(sizeOfMatrix, array)

      implicit none

      integer(kind=4) :: i, j, sizeOfMatrix
      real(kind=8), dimension(sizeOfMatrix,sizeOfMatrix) :: array

      forall(i = 1:sizeOfMatrix, j = 1:sizeOfMatrix) array(i,j) = (i/j)*(j/i) 

end subroutine identityCreate

subroutine cc299blktriad(maind,lower,upper,id,md,xb,x)

     !|     B(1)    C(1)             | | xb(1)  |
     !| A(2)  B(2)    C(2)           | |        |
     !|   A(3)  B(3)                 | |        |
     !|     .     .                  |*|        | = B[1:mb*3]
     !|       .     .                | |        |
     !|               .       C(mb-1)| |        |
     !|         A(mb)  B(mb)         | |xb(n*id)|

    ! id = inner matrices dimension.
    ! md = number matrices.
    ! maind = main diagonal of matrices  format: maind(id,id,md)
    ! lower = lower diagonal of matrices format: lower(id,id,2:md)
    ! upper = upper diagonal of matrices format: maind(id,id,md-1)
    ! xb    = B vector in Ax=B           format: xb(md*id)
    ! x     = x vector in Ax=B           format: x(md*id)
    !

    implicit none

    ! +++ Inputs +++
    !
    ! Scalar input variables.
    integer(kind=4) :: id,md

    ! Main diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: maind

    ! Lower diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: lower

    ! Upper diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: upper

    ! Vector of equalties B in Ax=B.
    real(kind=8), dimension(id,md)    :: xb

    ! Vector of answers x in Ax=B.
    real(kind=8), dimension(id,md)    :: x

    ! ++ Inside variables ++
    !
    ! Scalar variables
    integer(kind=4) :: i,ii,jj

    ! Array of gamma coefficients.
    real(kind=8), dimension(id,id,md) :: gamm

    ! Array of beta coefficients.
    real(kind=8), dimension(id,md) :: beta

    ! Auxiliary arrays.
    real(kind=8), allocatable, dimension(:,:) :: aux_copy
    real(kind=8), allocatable, dimension(:,:) :: aux_mult
    real(kind=8), allocatable, dimension(:,:) :: aux_summ
    real(kind=8), allocatable, dimension(:)   :: aux_dumm

    allocate(aux_mult(id,id))
    allocate(aux_summ(id,id))
    allocate(aux_copy(id,id))
    allocate(aux_dumm(id))


    !--------------------------------------------------------------------------!
    !                  Step 1: BLOCK TRIANGULARIZATION                         !
    !--------------------------------------------------------------------------!


    ! Zero out the auxiliar vectors.

    beta = 0.0d0
    gamm = 0.0d0


    ! Lets first get our first gamma.

    call dlacpy('A',id,id,maind(:,:,1),id,aux_copy,id)

    call inv(aux_copy,aux_copy,id)

    call dgemm('N','N',id,id,id,1.0d0,aux_copy,id, & 
        upper(:,:,1),id,1.0d0,gamm(:,:,1),id)


    ! Now that we have our first gamma, lets get the rest of then.

    do i = 2, md-1

        aux_mult = 0.0d0
        aux_summ = 0.0d0

        call dgemm('N','N',id,id,id,1.0d0,lower(:,:,i),id, &
            gamm(:,:,i-1),id,1.0d0,aux_mult,id)

        do jj = 1, id
            do ii = 1, id
                aux_summ(ii,jj) = maind(ii,jj,i) - aux_mult(ii,jj)
            end do
        end do

        call inv(aux_summ,aux_summ,id)

        call dgemm('N','N',id,id,id,1.0d0,aux_summ,id, &
            upper(:,:,i),id,1.0d0,gamm(:,:,i),id)
            
    end do

    
    ! Now that we have our gammas, lets get the betas, starting from the first
    ! ones. Note that now the calls done by the Lapack library will get a bit
    ! more complicated so lets use matmul...

    aux_copy = 0.0d0

    call dlacpy('A',id,id,maind(:,:,1),id,aux_copy,id)

    call inv(aux_copy,aux_copy,id)

    beta(:,1) = matmul(aux_copy,xb(:,1))


    ! We now have our first beta, lets get the rest.

    do i = 2, md

        aux_mult = 0.0d0
        aux_summ = 0.0d0

        call dgemm('N','N',id,id,id,1.0d0,lower(:,:,i),id, &
            gamm(:,:,i-1),id,1.0d0,aux_mult(:,:),id)

        !aux_mult(:,:) = matmul(lower(:,:,i),gamm(:,:,i-1))

        do jj = 1, id
            do ii = 1, id
                aux_summ(ii,jj) = maind(ii,jj,i) - aux_mult(ii,jj)
            end do
        end do

        call inv(aux_summ,aux_summ,id)

        aux_dumm(:) = xb(:,i) - matmul(lower(:,:,i),beta(:,i-1))

        beta(:,i) = matmul(aux_summ(:,:),aux_dumm(:))

    end do


    !--------------------------------------------------------------------------!
    !                  Step 2: BACKWARD SWEEP                                  !
    !--------------------------------------------------------------------------!


    ! How cool is that, lets start build our solution vector... iupiiii!

    x = 0.0d0

    x(:,md) = beta(:,md)

    do i = md-1,1,-1

        aux_dumm(:) = matmul(gamm(:,:,i),x(:,i+1))

        do ii = 1, id
            x(ii,i) = beta(ii,i) - aux_dumm(ii)
        end do

    end do


    deallocate(aux_mult)
    deallocate(aux_summ)
    deallocate(aux_copy)
    deallocate(aux_dumm)

end subroutine cc299blktriad

subroutine inv(A,A_inv,m)

  Implicit none
  integer :: m
  real(kind=8), dimension(m,m)::A, A_inv
  real(kind=8),dimension(m)::WORK
  integer,dimension(m)::IPIV
  integer info

  A_inv = A

  call DGETRF(M,M,A_inv,M,IPIV,info)

  if (info /=  0) then
    write(*,*)"DGETRF: Failed during matrix factorization"
    stop
  end if

  call DGETRI(M,A_inv,M,IPIV,WORK,M,info)

  if (info /=  0) then
   write(*,*)"DGETRI: Failed during matrix inversion."
   stop
  end if

end subroutine inv

subroutine cc299blktriad_opt(maind,lower,upper,id,md,xb,x)

     !|     B(1)    C(1)             | | xb(1)  |
     !| A(2)  B(2)    C(2)           | |        |
     !|   A(3)  B(3)                 | |        |
     !|     .     .                  |*|        | = B[1:mb*3]
     !|       .     .                | |        |
     !|               .       C(mb-1)| |        |
     !|         A(mb)  B(mb)         | |xb(n*id)|

    ! id = inner matrices dimension.
    ! md = number matrices.
    ! maind = main diagonal of matrices  format: maind(id,id,md)
    ! lower = lower diagonal of matrices format: lower(id,id,2:md)
    ! upper = upper diagonal of matrices format: maind(id,id,md-1)
    ! xb    = B vector in Ax=B           format: xb(md*id)
    ! x     = x vector in Ax=B           format: x(md*id)
    !

    implicit none

    ! +++ Inputs +++
    !
    ! Scalar input variables.
    integer(kind=4) :: id,md

    ! Main diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: maind

    ! Lower diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: lower

    ! Upper diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: upper

    ! Vector of equalties B in Ax=B.
    real(kind=8), dimension(id,md)    :: xb

    ! Vector of answers x in Ax=B.
    real(kind=8), dimension(id,md)    :: x

    ! ++ Inside variables ++
    !
    ! Scalar variables
    integer(kind=4) :: i,ii,jj

    ! Array of gamma coefficients.
    real(kind=8), dimension(id,id,md) :: gamm

    ! Array of beta coefficients.
    real(kind=8), dimension(id,md) :: beta

    ! Auxiliary arrays.
    real(kind=8), allocatable, dimension(:,:) :: aux_copy
    real(kind=8), allocatable, dimension(:,:) :: aux_mult
    real(kind=8), allocatable, dimension(:,:) :: aux_summ
    real(kind=8), allocatable, dimension(:)   :: aux_dumm

    allocate(aux_mult(id,id))
    allocate(aux_summ(id,id))
    allocate(aux_copy(id,id))
    allocate(aux_dumm(id))


    !--------------------------------------------------------------------------!
    !                  Step 1: BLOCK TRIANGULARIZATION                         !
    !--------------------------------------------------------------------------!


    ! Zero out the auxiliar vectors.

    beta = 0.0d0
    gamm = 0.0d0


    ! Lets first get our first gamma.

    aux_copy = maind(:,:,1)


    call inv(aux_copy,aux_copy,id)


    gamm(:,:,1) = matmul(aux_copy,upper(:,:,1))


    ! Now that we have our first gamma, lets get the rest of then.

    do i = 2, md-1

        aux_mult = 0.0d0
        aux_summ = 0.0d0

        aux_mult = matmul(lower(:,:,i),gamm(:,:,i-1))

        do jj = 1, id
            do ii = 1, id
                aux_summ(ii,jj) = maind(ii,jj,i) - aux_mult(ii,jj)
            end do
        end do

        call inv(aux_summ,aux_summ,id)

        gamm(:,:,i) = matmul(aux_summ,upper(:,:,i))
            
    end do

    
    ! Now that we have our gammas, lets get the betas, starting from the first
    ! ones. Note that now the calls done by the Lapack library will get a bit
    ! more complicated so lets use matmul...

    aux_copy = 0.0d0

    aux_copy = maind(:,:,1)

    call inv(aux_copy,aux_copy,id)

    beta(:,1) = matmul(aux_copy,xb(:,1))


    ! We now have our first beta, lets get the rest.

    do i = 2, md

        aux_mult = 0.0d0
        aux_summ = 0.0d0

        aux_mult(:,:) = matmul(lower(:,:,i),gamm(:,:,i-1))

        do jj = 1, id
            do ii = 1, id
                aux_summ(ii,jj) = maind(ii,jj,i) - aux_mult(ii,jj)
            end do
        end do

        call inv(aux_summ,aux_summ,id)

        aux_dumm(:) = xb(:,i) - matmul(lower(:,:,i),beta(:,i-1))

        beta(:,i) = matmul(aux_summ(:,:),aux_dumm(:))

    end do


    !--------------------------------------------------------------------------!
    !                  Step 2: BACKWARD SWEEP                                  !
    !--------------------------------------------------------------------------!


    ! How cool is that, lets start build our solution vector... iupiiii!

    x = 0.0d0

    x(:,md) = beta(:,md)

    do i = md-1,1,-1

        aux_dumm(:) = matmul(gamm(:,:,i),x(:,i+1))

        do ii = 1, id
            x(ii,i) = beta(ii,i) - aux_dumm(ii)
        end do

    end do


    deallocate(aux_mult)
    deallocate(aux_summ)
    deallocate(aux_copy)
    deallocate(aux_dumm)

end subroutine cc299blktriad_opt


