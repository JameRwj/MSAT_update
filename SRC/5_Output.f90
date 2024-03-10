    MODULE Output
    !**********************************************************************************************
    !******************This module is used to output the result of the computation*****************
    !******************Contains Output_Matrix, Output_1D, Output_Initialization********************
    !**********************************************************************************************

    IMPLICIT NONE
    
    CONTAINS
    
    !**********************************************************************************************
    !****************************************Output_Matrix*****************************************
    !**********************************************************************************************
    SUBROUTINE Output_Matrix()
    !********************This subroutine outputs the result of Matrix analysis*********************

    USE Common_Data,ONLY:p2,WR,WI,id,jd,VectorR,VectorL,zero

    IMPLICIT NONE

    INTEGER :: N_Matrix,i,j,k,n,m,ii,jj
    REAL(p2):: max_eigvalue

    CHARACTER*50 :: fileName
    CHARACTER*50 :: number

    REAL(p2),DIMENSION(:),ALLOCATABLE :: eigVector_rho,eigVector_u,eigVector_v,eigVector_p
    INTEGER ,DIMENSION(:),ALLOCATABLE :: position

    ALLOCATE(eigVector_rho((id-1)*(jd-1)),eigVector_u((id-1)*(jd-1)),eigVector_v((id-1)*(jd-1)),&
             eigVector_p((id-1)*(jd-1)),position(4*(id-1)*(jd-1)))

    N_Matrix = SIZE(WR,1)

    !****************************Output the scatter of all eigenvalues*****************************
    OPEN(115,FILE="Scatter.plt",FORM="FORMATTED")
    DO i = 1,N_Matrix
        WRITE(115,*)WR(i),WI(i)
    END DO
    CLOSE(115)
    
    !********************Determine the maximum eigenvalues and their locations*********************
    Max_eigvalue = MAXVAL(WR)
    k = 0
    position = zero
    DO i = 1,N_Matrix
        IF (WR(i) == Max_eigvalue) THEN
            k = k + 1
            position(k) = i
        END IF
    END DO
    
    !**************************Write the maximum eigenvalues on the screen*************************
    WRITE(*,"(A80)")"=====================================Result====================================="
    Write(*,"(A26)")"The maximal eigenvalue is:"
    DO i = 1,k
        WRITE(*,"(A26,F16.8,A4,F16.8,A2)")" ",WR(position(i)),"+",WI(position(i)),"i"
    END DO

    !**********************Out put the eigenvalues of the maximum eigenvalues**********************
    IF (Max_eigvalue > zero) THEN
        DO i = 1,k
            DO j = 1,N_Matrix,4
                n = 1 + FLOOR(j/4.0)
                eigVector_rho(n) = VectorR(j  ,position(i))
                eigVector_u(n)   = VectorR(j+1,position(i))
                eigVector_v(n)   = VectorR(j+2,position(i))
                eigVector_p(n)   = VectorR(j+3,position(i))
            END DO

            WRITE(number,*)i

            fileName = "Eigenvector_rho_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(115,FILE=TRIM(ADJUSTL(fileName)),FORM="FORMATTED")
            WRITE(115,*)' variables="x","y","<greek>r</greek>" '
            WRITE(115,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            fileName = "Eigenvector_u_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(116,FILE=TRIM(ADJUSTL(fileName)),FORM="FORMATTED")
            WRITE(116,*)' variables="x","y","u" '
            WRITE(116,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            fileName = "Eigenvector_v_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(117,FILE=TRIM(ADJUSTL(fileName)),FORM="FORMATTED")
            WRITE(117,*)' variables="x","y","v" '
            WRITE(117,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            fileName = "Eigenvector_p_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(118,FILE=fileName,FORM="FORMATTED")
            WRITE(118,*)' variables="x","y","p" '
            WRITE(118,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            m = 1
            DO ii = 1,jd-1
                DO jj = 1,id-1
                    WRITE(115,*)jj,ii,eigVector_rho(m)
                    WRITE(116,*)jj,ii,eigVector_u(m)
                    WRITE(117,*)jj,ii,eigVector_v(m)
                    WRITE(118,*)jj,ii,eigVector_p(m)
                    m = m+1
                END DO
            END DO

            CLOSE(115)
            CLOSE(116)
            CLOSE(117)
            CLOSE(118)
        END DO
    END IF

    DEALLOCATE(eigVector_rho,eigVector_u,eigVector_v,eigVector_p,position)
    END SUBROUTINE Output_Matrix
    
    
    
    !**********************************************************************************************
    !******************************************output_1D*******************************************
    !**********************************************************************************************
    SUBROUTINE output_1D()
    !**************************This subroutine outputs the 1D computation**************************
    USE Common_Data,ONLY:half,x_1D,id,rho_1D,u_1D,p_1D,gamma,one

    IMPLICIT NONE

    INTEGER :: i

    OPEN(101,FILE="Result_1D.plt")
    WRITE(101,"(1X,A)") "VARIABLES= ""X"", ""Density"", ""Velocity"", ""Pressure"",""Total Energy"",""Mach Number"",""Mass flux"""
    WRITE(101,"(1X,A)") "ZONE DATAPACKING=POINT"
    DO i=1,id-1
        WRITE(101,"(7E20.12)") half*(x_1D(i)+x_1D(i+1)),&
            rho_1D(i),&
            u_1D(i),  &
            p_1D(i),  &
            p_1D(i)/(gamma-one)+half*rho_1D(i)*u_1D(i)*u_1D(i),&
            u_1D(i)/SQRT(gamma*p_1D(i)/rho_1D(i)),&
            rho_1D(i)*u_1D(i)
    END DO
    CLOSE(101)

    END SUBROUTINE output_1D
    
    
    
    !**********************************************************************************************
    !************************************Output_Initialization*************************************
    !**********************************************************************************************
    SUBROUTINE Output_Initialization()
    !*********************This subroutine outputs the result of initialization*********************
    USE Common_Data,ONLY : x,y,rho,u,v,p,id,jd
    
    IMPLICIT NONE
    
    INTEGER :: i,j
    
    OPEN(113,FILE="Result_FlowField.plt",FORM="FORMATTED")
    WRITE(113,*) "TITLE             = result"
    WRITE(113,*) "VARIABLES=   ", "x   ","y   ","rho   ","u   ","v   ","p   "
    WRITE(113,*) "ZONE T=","ZONE1"
    WRITE(113,*) "STRANDID=0, SOLUTIONTIME=0"
    WRITE(113,*) "I= ",id," J= ",jd," K= ",1
    WRITE(113,*) "DATAPACKING=BLOCK"
    WRITE(113,*) "VARLOCATION=([3-6]=CELLCENTERED)"
    
    DO j=1,jd
        DO i = 1,id
            WRITE(113,*)x(i,j)
        END DO
    END DO
    DO j=1,jd
        DO i = 1,id
            WRITE(113,*)y(i,j)
        END DO
    END DO
    DO j=1,jd-1
        DO i = 1,id-1
            WRITE(113,*)rho(i,j)
        END DO
    END DO
    DO j=1,jd-1
        DO i = 1,id-1
            WRITE(113,*)u(i,j)
        END DO
    END DO
    DO j=1,jd-1
        DO i = 1,id-1
            WRITE(113,*)v(i,j)
        END DO
    END DO
    DO j=1,jd-1
        DO i = 1,id-1
            WRITE(113,*)p(i,j)
        END DO
    END DO

    CLOSE(113)
    
    END SUBROUTINE Output_Initialization

    END MODULE Output