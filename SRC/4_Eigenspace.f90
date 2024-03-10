    MODULE Eigenspace
    !**********************************************************************************************
    !******This module is used to assembly the stability matrix and calculate the eigenvalue*******
    !******Only has one subroutine: Calculate_Eigen*******
    !**********************************************************************************************

    IMPLICIT NONE

    CONTAINS

    !**********************************************************************************************
    !*****************************************Calculate_Eigen***************************************
    !**********************************************************************************************
    SUBROUTINE Calculate_Eigen()

    USE Common_Data
    USE Calculate_Gradients

    IMPLICIT NONE

    REAL(p2) :: len1,len2,len3,len4,s                           !The length of four faces and the area of the cell
    REAL(p2) :: nx1,ny1,nx2,ny2,nx3,ny3,nx4,ny4                 !The unit normal vector of for faces
    REAL(p2) :: rho_,u_,v_,q_                                   !The local state

    INTEGER  :: col,row                                         !The number of row and column
    INTEGER  :: n,i,j,ii,jj
    INTEGER  :: N_Matrix,INFO,LDA,LDVL,LDVR,LWORK               !Parameters needed by DGEEV

    REAL(p2),DIMENSION(4,4) :: alphaL1_1,alphaL1_2,alphaL1_3,alphaL1_4,alphaL1_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: alphaL2_1,alphaL2_2,alphaL2_3,alphaL2_4,alphaL2_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: alphaL3_1,alphaL3_2,alphaL3_3,alphaL3_4,alphaL3_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: alphaL4_1,alphaL4_2,alphaL4_3,alphaL4_4,alphaL4_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: alphaR1_1,alphaR1_2,alphaR1_3,alphaR1_4,alphaR1_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: alphaR2_1,alphaR2_2,alphaR2_3,alphaR2_4,alphaR2_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: alphaR3_1,alphaR3_2,alphaR3_3,alphaR3_4,alphaR3_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: alphaR4_1,alphaR4_2,alphaR4_3,alphaR4_4,alphaR4_5!Parameters matrix of the stability matrix
    REAL(p2),DIMENSION(4,4) :: A1,A2,A3,A4                      !Gradient matrix
    REAL(p2),DIMENSION(4,4) :: B1,B2,B3,B4                      !Gradient matrix
    REAL(p2),DIMENSION(4,4) :: C1,C2,C3,C4                      !Gradient matrix
    REAL(p2),DIMENSION(4,4) :: D1,D2,D3,D4                      !Gradient matrix
    REAL(p2),DIMENSION(4,4) :: L1,L2,L3,L4                      !Gradient matrix
    REAL(p2),DIMENSION(4,4) :: R1,R2,R3,R4                      !Gradient matrix
    REAL(p2),DIMENSION(4,4) :: dW_dU                            !Transformation matrix between conservative primitive variables
    REAL(p2),DIMENSION(4,4) :: E                                !unit matrix
    REAL(p2),DIMENSION(4,4) :: Ln1,Ln2,Ln3,Ln4
    REAL(p2),DIMENSION(4,4) :: Ln1_L,Ln2_L,Ln3_L,Ln4_L
    REAL(p2),DIMENSION(4,4) :: Ln1_R,Ln2_R,Ln3_R,Ln4_R

    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: Matrix               !Stability matrix
    REAL(p2),DIMENSION(:)  ,ALLOCATABLE :: WORK                 !Parameters needed by DGEEV

    row = jd-1
    col = id-1

    ALLOCATE(Matrix(4*row*col,4*col*row))

    Matrix = zero
    dW_dU  = zero

    A1 = zero; A2 = zero; A3 = zero; A4 = zero
    B1 = zero; B2 = zero; B3 = zero; B4 = zero
    C1 = zero; C2 = zero; C3 = zero; C4 = zero
    D1 = zero; D2 = zero; D3 = zero; D4 = zero
    L1 = zero; L2 = zero; L3 = zero; L4 = zero
    R1 = zero; R2 = zero; R3 = zero; R4 = zero

    E  = zero
    Ln1 = zero; Ln2 = zero; Ln3 = zero; Ln4 = zero
    Ln1_L = zero; Ln2_L = zero; Ln3_L = zero; Ln4_L = zero
    Ln1_R = zero; Ln2_R = zero; Ln3_R = zero; Ln4_R = zero

    DO i = 1,4
        E(i,i) = one
    END DO

    DO j = 1,row
        DO i = 1,col
            n = i+(j-1)*col

            !********************Calculate the information of the current cell*********************
            nx1  = nxx(i+1,j)
            ny1  = nxy(i+1,j)
            nx2  = nyx(i,j+1)
            ny2  = nyy(i,j+1)
            nx3  = nxx(i,j)
            ny3  = nxy(i,j)
            nx4  = nyx(i,j)
            ny4  = nyy(i,j)
            len1 = SQRT(sxx(i+1,j+0)**2+sxy(i+1,j+0)**2)
            len2 = SQRT(syx(i+0,j+1)**2+syy(i+0,j+1)**2)
            len3 = SQRT(sxx(i+0,j+0)**2+sxy(i+0,j+0)**2)
            len4 = SQRT(syx(i+0,j+0)**2+syy(i+0,j+0)**2)
            s    = vol(i,j)

            !************************Calculate the gradients of four faces*************************
            !The gradient of the variables on the left side
            CALL Gradientofflux_L(i,j,1,nx1,ny1,L1,alphaL1_1,alphaL1_2,alphaL1_3,alphaL1_4,alphaL1_5,Ln1_L) !Right
            CALL Gradientofflux_L(i,j,2,nx2,ny2,L2,alphaL2_1,alphaL2_2,alphaL2_3,alphaL2_4,alphaL2_5,Ln2_L) !Upper
            CALL Gradientofflux_L(i,j,3,nx3,ny3,L3,alphaL3_1,alphaL3_2,alphaL3_3,alphaL3_4,alphaL3_5,Ln3_L) !Left
            CALL Gradientofflux_L(i,j,4,nx4,ny4,L4,alphaL4_1,alphaL4_2,alphaL4_3,alphaL4_4,alphaL4_5,Ln4_L) !Lower

            !The gradient of the variables on the right side
            CALL Gradientofflux_R(i,j,1,nx1,ny1,R1,alphaR1_1,alphaR1_2,alphaR1_3,alphaR1_4,alphaR1_5,Ln1_R) !Right
            CALL Gradientofflux_R(i,j,2,nx2,ny2,R2,alphaR2_1,alphaR2_2,alphaR2_3,alphaR2_4,alphaR2_5,Ln2_R) !Upper
            CALL Gradientofflux_R(i,j,3,nx3,ny3,R3,alphaR3_1,alphaR3_2,alphaR3_3,alphaR3_4,alphaR3_5,Ln3_R) !Left
            CALL Gradientofflux_R(i,j,4,nx4,ny4,R4,alphaR4_1,alphaR4_2,alphaR4_3,alphaR4_4,alphaR4_5,Ln4_R) !Lower

            SELECT CASE(Reconstruction_Variables)
            CASE(1)           !Reconstruct the primitive variable

                rho_ = rho(i,j)
                u_   = u(i,j)
                v_   = v(i,j)
                q_   = u_**2+v_**2

                dW_dU = RESHAPE([one  , -u_/rho_ , -v_/rho_ , (gamma-one)*q_/two,&
                    zero , one/rho_ , zero     , -(gamma-one)*u_,&
                    zero , zero     , one/rho_ , -(gamma-one)*v_,&
                    zero , zero     , zero     , gamma-one] ,[4,4])

                A1 = ( MATMUL(L1,alphaL1_3) + MATMUL(R1,alphaR1_4) )*len1/s
                A2 = ( MATMUL(L2,alphaL2_3) + MATMUL(R2,alphaR2_4) )*len2/s
                A3 = ( MATMUL(R3,alphaR3_3) + MATMUL(L3,alphaL3_4) )*len3/s
                A4 = ( MATMUL(R4,alphaR4_3) + MATMUL(L4,alphaL4_4) )*len4/s

                B1 = ( len1*( MATMUL(L1,alphaL1_4)+MATMUL(R1,alphaR1_3) ) &
                    +  len3*( MATMUL(L3,alphaL3_5)+MATMUL(R3,alphaR3_2) ) )/s
                B2 = ( len2*( MATMUL(L2,alphaL2_4)+MATMUL(R2,alphaR2_3) ) &
                    +  len4*( MATMUL(L4,alphaL4_5)+MATMUL(R4,alphaR4_2) ) )/s
                B3 = ( len3*( MATMUL(R3,alphaR3_4)+MATMUL(L3,alphaL3_3) ) &
                    +  len1*( MATMUL(R1,alphaR1_5)+MATMUL(L1,alphaL1_2) ) )/s
                B4 = ( len4*( MATMUL(R4,alphaR4_4)+MATMUL(L4,alphaL4_3) ) &
                    +  len2*( MATMUL(R2,alphaR2_5)+MATMUL(L2,alphaL2_2) ) )/s

                C1 = ( len1*( MATMUL(L1,alphaL1_5)+MATMUL(R1,alphaR1_2) ) &
                    +  len3*MATMUL(R3,alphaR3_1) )/s
                C2 = ( len2*( MATMUL(L2,alphaL2_5)+MATMUL(R2,alphaR2_2) ) &
                    +  len4*MATMUL(R4,alphaR4_1) )/s
                C3 = ( len3*( MATMUL(R3,alphaR3_5)+MATMUL(L3,alphaL3_2) ) &
                    +  len1*MATMUL(L1,alphaL1_1) )/s
                C4 = ( len4*( MATMUL(R4,alphaR4_5)+MATMUL(L4,alphaL4_2) ) &
                    +  len2*MATMUL(L2,alphaL2_1) )/s

                D1 = MATMUL(R1,alphaR1_1)*len1/s
                D2 = MATMUL(R2,alphaR2_1)*len2/s
                D3 = MATMUL(L3,alphaL3_1)*len3/s
                D4 = MATMUL(L4,alphaL4_1)*len4/s

                A1 = MATMUL(dW_DU,A1); A2 = MATMUL(dW_DU,A2); A3 = MATMUL(dW_DU,A3); A4 = MATMUL(dW_DU,A4)
                B1 = MATMUL(dW_DU,B1); B2 = MATMUL(dW_DU,B2); B3 = MATMUL(dW_DU,B3); B4 = MATMUL(dW_DU,B4)
                C1 = MATMUL(dW_DU,C1); C2 = MATMUL(dW_DU,C2); C3 = MATMUL(dW_DU,C3); C4 = MATMUL(dW_DU,C4)
                D1 = MATMUL(dW_DU,D1); D2 = MATMUL(dW_DU,D2); D3 = MATMUL(dW_DU,D3); D4 = MATMUL(dW_DU,D4)
            CASE(2)           !Reconstruct the characteristic variable
                DO ii = 1,4
                    DO jj = 1,4
                        IF ((Ln1_L(ii,jj)/=Ln1_R(ii,jj)) .OR. (Ln2_L(ii,jj)/=Ln2_R(ii,jj)) .OR.&
                            (Ln3_L(ii,jj)/=Ln3_R(ii,jj)) .OR. (Ln4_L(ii,jj)/=Ln4_R(ii,jj))) THEN
                            WRITE(*,*)" "
                            WRITE(*,*)"  The characteristic matrices of the left side and right side are unequal"
                            PAUSE
                            STOP
                        ELSE
                            Ln1(ii,jj) = Ln1_L(ii,jj); Ln2(ii,jj) = Ln2_L(ii,jj)
                            Ln3(ii,jj) = Ln3_L(ii,jj); Ln4(ii,jj) = Ln4_L(ii,jj)
                        END IF
                    END DO
                END DO

                A1 = MATMUL( ( MATMUL(L1,alphaL1_3) + MATMUL(R1,alphaR1_4) ),Ln1 )*len1/s
                A2 = MATMUL( ( MATMUL(L2,alphaL2_3) + MATMUL(R2,alphaR2_4) ),Ln2 )*len2/s
                A3 = MATMUL( ( MATMUL(R3,alphaR3_3) + MATMUL(L3,alphaL3_4) ),Ln3 )*len3/s
                A4 = MATMUL( ( MATMUL(R4,alphaR4_3) + MATMUL(L4,alphaL4_4) ),Ln4 )*len4/s

                B1 =  MATMUL( ( MATMUL(L1,alphaL1_4)+MATMUL(R1,alphaR1_3) ),Ln1 )*len1/s &
                    + MATMUL( ( MATMUL(L3,alphaL3_5)+MATMUL(R3,alphaR3_2) ),Ln3 )*len3/s
                B2 =  MATMUL( ( MATMUL(L2,alphaL2_4)+MATMUL(R2,alphaR2_3) ),Ln2 )*len2/s &
                    + MATMUL( ( MATMUL(L4,alphaL4_5)+MATMUL(R4,alphaR4_2) ),Ln4 )*len4/s
                B3 =  MATMUL( ( MATMUL(R3,alphaR3_4)+MATMUL(L3,alphaL3_3) ),Ln3 )*len3/s &
                    + MATMUL( ( MATMUL(R1,alphaR1_5)+MATMUL(L1,alphaL1_2) ),Ln1 )*len1/s
                B4 =  MATMUL( ( MATMUL(R4,alphaR4_4)+MATMUL(L4,alphaL4_3) ),Ln4 )*len4/s &
                    + MATMUL( ( MATMUL(R2,alphaR2_5)+MATMUL(L2,alphaL2_2) ),Ln2 )*len2/s

                C1 = MATMUL( ( MATMUL(L1,alphaL1_5)+MATMUL(R1,alphaR1_2) ),Ln1 )*len1/s &
                    +MATMUL( MATMUL(R3,alphaR3_1),Ln3 )*len3/s
                C2 = MATMUL( ( MATMUL(L2,alphaL2_5)+MATMUL(R2,alphaR2_2) ),Ln2 )*len2/s &
                    +MATMUL( MATMUL(R4,alphaR4_1),Ln4 )*len4/s
                C3 = MATMUL( ( MATMUL(R3,alphaR3_5)+MATMUL(L3,alphaL3_2) ),Ln3 )*len3/s &
                    +MATMUL( MATMUL(L1,alphaL1_1),Ln1 )*len1/s
                C4 = MATMUL( ( MATMUL(R4,alphaR4_5)+MATMUL(L4,alphaL4_2) ),Ln4 )*len4/s &
                    +MATMUL( MATMUL(L2,alphaL2_1),Ln2 )*len2/s

                D1 = MATMUL( MATMUL(R1,alphaR1_1),Ln1 )*len1/s
                D2 = MATMUL( MATMUL(R2,alphaR2_1),Ln2 )*len2/s
                D3 = MATMUL( MATMUL(L3,alphaL3_1),Ln3 )*len3/s
                D4 = MATMUL( MATMUL(L4,alphaL4_1),Ln4 )*len4/s
            END SELECT

            !****************************Assembly the stability matrix*****************************
            Matrix(4*n-3:4*n,4*n-3:4*n)= -(A1+A2+A3+A4)

            IF (i >= 2) THEN
                Matrix(4*n-3:4*n,4*n-7:4*n-4) = -B3
            END IF

            IF (i <= col-1) THEN
                Matrix(4*n-3:4*n,4*n+1:4*n+4) = -B1
            END IF

            IF (j >= 2) THEN
                Matrix(4*n-3:4*n,4*n-4*col-3:4*n-4*col) = -B4
            END IF

            IF (j <= row-1) THEN
                Matrix(4*n-3:4*n,4*n+4*col-3:4*n+4*col) = -B2
            END IF

            IF (i >= 3) THEN
                Matrix(4*n-3:4*n,4*n-11:4*n-8) = -C3
            END IF

            IF (i <= col-2) THEN
                Matrix(4*n-3:4*n,4*n+5:4*n+8) = -C1
            END IF

            IF (j >= 3) THEN
                Matrix(4*n-3:4*n,4*n-4*2*col-3:4*n-4*2*col) = -C4
            END IF

            IF (j <= row-2) THEN
                Matrix(4*n-3:4*n,4*n+4*2*col-3:4*n+4*2*col) = -C2
            END IF

            IF (i >= 4) THEN
                Matrix(4*n-3:4*n,4*n-15:4*n-12) = -D3
            END IF

            IF (i <= col-3) THEN
                Matrix(4*n-3:4*n,4*n+9:4*n+12) = -D1
            END IF

            IF (j >= 4) THEN
                Matrix(4*n-3:4*n,4*n-4*3*col-3:4*n-4*3*col) = -D4
            END IF

            IF (j <= row-3) THEN
                Matrix(4*n-3:4*n,4*n+4*3*col-3:4*n+4*3*col) = -D2
            END IF

        END DO
    END DO

    !**********************Calculate the eigenvalue and eigenvector by DGEEV***********************
    N_Matrix = 4*col*row

    LDA  = N_Matrix
    LDVL = N_Matrix
    LDVR = N_Matrix

    ALLOCATE(WR(N_Matrix),WI(N_Matrix),WORK(4*N_Matrix),VectorL(LDVL,N_Matrix),VectorR(LDVR,N_Matrix))

    LWORK = 4*N_Matrix
    CALL DGEEV("Vectors","Vectors",N_Matrix,Matrix,LDA,WR,WI,VectorL,LDVL,VectorR,LDVR,WORK,LWORK,INFO)

    DEALLOCATE(WORK,Matrix)

    END SUBROUTINE Calculate_Eigen

    END MODULE Eigenspace