    MODULE Calculate_Gradients
    !**********************************************************************************************
    !***************This module is used to calculate the gradient matrix, including****************
    !***************Gradientofflux_L and Gradientofflux_R******************************************
    !**********************************************************************************************
    USE Common_Data,ONLY:p2,rho,u,v,p,zero,two,Reconstruction_Method
    USE Reconstruction
    USE Calculate_Flux

    IMPLICIT NONE

    CONTAINS

    !**********************************************************************************************
    !***************************************Gradientofflux_L***************************************
    !**********************************************************************************************
    SUBROUTINE Gradientofflux_L(i,j,Face,nx,ny,L,alphaL_1,alphaL_2,alphaL_3,alphaL_4,alphaL_5,Ln)
    !*****This subroutine calculate the gradient of the variables on the left side of the face*****
    IMPLICIT NONE

    INTEGER i,j,Face,m,n
    REAL(p2),DIMENSION(4,4) :: L,alphaL_1,alphaL_2,alphaL_3,alphaL_4,alphaL_5  !Gradient matrix and coefficient matrix
    REAL(p2),DIMENSION(4)   :: flux1_plus ,flux2_plus ,flux3_plus ,flux4_plus  !Numerical flux
    REAL(p2),DIMENSION(4)   :: flux1_minus,flux2_minus,flux3_minus,flux4_minus !Numerical flux
    REAL(p2),DIMENSION(4)   :: UL,UR,UL0,UR0                !Variables on both sides of the current face
    REAL(p2),DIMENSION(4,4) :: Ln,Rn

    REAL(p2),DIMENSION(4,5) :: alphaL,alphaR

    REAL(p2) :: delta,deltaL,deltaR
    REAL(p2) :: nx,ny

    alphaL_1 = zero
    alphaL_2 = zero
    alphaL_3 = zero
    alphaL_4 = zero
    alphaL_5 = zero
    Ln = zero
    Rn = zero
    
    delta = 1.0E-7_p2

    !**********************************Reconstruct the variables***********************************    
    CALL Reconstruct(i,j,Face,nx,ny,alphaL,UL,alphaR,UR,Ln,Rn)
    UL0 = UL
    UR0 = UR

    deltaL = delta
    deltaR = zero

    !************************************Calculate the gradient***********************************
    !The first item plus delta
    UL = UL0
    UR = UR0

    UL(1) = UL(1)+deltaL
    UR(1) = UR(1)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux1_plus)

    !The first item subtracts delta
    UL = UL0
    UR = UR0

    UL(1) = UL(1)-deltaL
    UR(1) = UR(1)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux1_minus)

    !The second item plus delta
    UL = UL0
    UR = UR0

    UL(2) = UL(2)+deltaL
    UR(2) = UR(2)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux2_plus)

    !The second item subtracts delta
    UL = UL0
    UR = UR0

    UL(2) = UL(2)-deltaL
    UR(2) = UR(2)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux2_minus)

    !The third item plus delta
    UL = UL0
    UR = UR0

    UL(3) = UL(3)+deltaL
    UR(3) = UR(3)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux3_plus)

    !The third item subtracts delta
    UL = UL0
    UR = UR0

    UL(3) = UL(3)-deltaL
    UR(3) = UR(3)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux3_minus)

    !The forth item plus delta
    UL = UL0
    UR = UR0

    UL(4) = UL(4)+deltaL
    UR(4) = UR(4)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux4_plus)

    !The forth item subtracts delta
    UL = UL0
    UR = UR0

    UL(4) = UL(4)-deltaL
    UR(4) = UR(4)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux4_minus)

    !*********************************Calculate the gradient matrix********************************
    L(:,1) = (flux1_plus-flux1_minus)/two/delta
    L(:,2) = (flux2_plus-flux2_minus)/two/delta
    L(:,3) = (flux3_plus-flux3_minus)/two/delta
    L(:,4) = (flux4_plus-flux4_minus)/two/delta

    IF ((Face == 3) .OR. (Face == 4)) THEN
        DO m = 1,4
            DO n = 1,4
                L(m,n) = -L(m,n)!This is because the direction of the unit normal vector of the left and lower faces
            END DO
        END DO
    END IF

    !*******************************Calculate the coefficient matrix*******************************
    DO n = 1,4
        alphaL_1(n,n) = alphaL(n,1)
        alphaL_2(n,n) = alphaL(n,2)
        alphaL_3(n,n) = alphaL(n,3)
        alphaL_4(n,n) = alphaL(n,4)
        alphaL_5(n,n) = alphaL(n,5)
    END DO

    END SUBROUTINE Gradientofflux_L



    !**********************************************************************************************
    !***************************************Gradientofflux_R***************************************
    !**********************************************************************************************
    SUBROUTINE Gradientofflux_R(i,j,Face,nx,ny,R,alphaR_1,alphaR_2,alphaR_3,alphaR_4,alphaR_5,Ln)
    !*****This subroutine calculate the gradient of the variables on the left right of the face****

    IMPLICIT NONE

    INTEGER i,j,Face,m,n
    REAL(p2),DIMENSION(4,4) :: R,alphaR_1,alphaR_2,alphaR_3,alphaR_4,alphaR_5  !Gradient matrix and coefficient matrix
    REAL(p2),DIMENSION(4)   :: flux1_plus ,flux2_plus ,flux3_plus ,flux4_plus  !Numerical flux
    REAL(p2),DIMENSION(4)   :: flux1_minus,flux2_minus,flux3_minus,flux4_minus !Numerical flux
    REAL(p2),DIMENSION(4)   :: UL,UR,UL0,UR0             !Variables on both sides of the current face
    REAL(p2),DIMENSION(4,5) :: alphaL,alphaR
    REAL(p2),DIMENSION(4,4) :: Ln,Rn
    
    REAL(p2) :: delta,deltaL,deltaR
    REAL(p2) :: nx,ny

    alphaR_1 = zero
    alphaR_2 = zero
    alphaR_3 = zero
    alphaR_4 = zero
    alphaR_5 = zero
    Ln = zero
    Rn = zero
    
    delta = 1.0E-7_p2

    !**********************************Reconstruct the variables***********************************
    
    CALL Reconstruct(i,j,Face,nx,ny,alphaL,UL,alphaR,UR,Ln,Rn)
    UL0 = UL
    UR0 = UR
    
    deltaL = zero
    deltaR = delta

    !************************************Calculate the gradient***********************************
    !The first item plus delta
    UL = UL0
    UR = UR0

    UL(1) = UL(1)+deltaL
    UR(1) = UR(1)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux1_plus)

    !The first item subtracts delta
    UL = UL0
    UR = UR0

    UL(1) = UL(1)-deltaL
    UR(1) = UR(1)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux1_minus)

    !The second item plus delta
    UL = UL0
    UR = UR0

    UL(2) = UL(2)+deltaL
    UR(2) = UR(2)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux2_plus)

    !The second item subtracts delta
    UL = UL0
    UR = UR0

    UL(2) = UL(2)-deltaL
    UR(2) = UR(2)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux2_minus)

    !The third item plus delta
    UL = UL0
    UR = UR0

    UL(3) = UL(3)+deltaL
    UR(3) = UR(3)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux3_plus)

    !The third item subtracts delta
    UL = UL0
    UR = UR0

    UL(3) = UL(3)-deltaL
    UR(3) = UR(3)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux3_minus)

    !The forth item plus delta
    UL = UL0
    UR = UR0

    UL(4) = UL(4)+deltaL
    UR(4) = UR(4)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux4_plus)

    !The forth item subtracts delta
    UL = UL0
    UR = UR0

    UL(4) = UL(4)-deltaL
    UR(4) = UR(4)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,Rn,flux4_minus)

    !*********************************Calculate the gradient matrix********************************
    R(:,1) = (flux1_plus-flux1_minus)/two/delta
    R(:,2) = (flux2_plus-flux2_minus)/two/delta
    R(:,3) = (flux3_plus-flux3_minus)/two/delta
    R(:,4) = (flux4_plus-flux4_minus)/two/delta

    IF ((Face == 3) .OR. (Face == 4)) THEN
        DO m = 1,4
            DO n = 1,4
                R(m,n) = -R(m,n)!This is because the direction of the unit normal vector of the left and lower faces
            END DO
        END DO
    END IF

    !*******************************Calculate the coefficient matrix*******************************
    DO n = 1,4
        alphaR_1(n,n) = alphaR(n,1)
        alphaR_2(n,n) = alphaR(n,2)
        alphaR_3(n,n) = alphaR(n,3)
        alphaR_4(n,n) = alphaR(n,4)
        alphaR_5(n,n) = alphaR(n,5)
    END DO

    END SUBROUTINE Gradientofflux_R
    END MODULE Calculate_Gradients