    MODULE Reconstruction
    !**********************************************************************************************
    !********This module reconstruct the variables and calculate the coefficient, including********
    !********Reconstruct_L ,Reconstruct_R, MUSCL_L, MUSCL_R, ROUND_L, ROUND_R**********************
    !**********************************************************************************************

    USE Common_Data,ONLY: Reconstruction_Method,Reconstruction_Variables,&
        p2,id,jd,ghostLayers,zero,one,two,half,rho,u,v,p,gamma
    USE Reconstruction_1D,ONLY: Limiter_function

    CONTAINS

    !**********************************************************************************************
    !******************************************Reconstruct*****************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruct(i,j,Face,nx,ny,alphaL,UUL,alphaR,UUR,Ln,Rn)

    IMPLICIT NONE

    INTEGER ::i,j,Face,n,m
    REAL(p2)::nx,ny

    REAL(p2) :: rhoL,uL,vL,pL
    REAL(p2) :: rhoR,uR,vR,pR
    REAL(p2) :: rho_L,u_L,v_L,p_L
    REAL(p2) :: rho_R,u_R,v_R,p_R
    REAL(p2) :: C1L,C2L,C3L,C4L
    REAL(p2) :: C1R,C2R,C3R,C4R

    REAL(p2),DIMENSION(4)   :: UUL,UUR
    REAL(p2),DIMENSION(5)   :: alphaL1,alphaL2,alphaL3,alphaL4
    REAL(p2),DIMENSION(5)   :: alphaR1,alphaR2,alphaR3,alphaR4
    REAL(p2),DIMENSION(4,5) :: alphaL,alphaR
    REAL(p2),DIMENSION(4,4) :: Ln,Rn
    REAL(p2),DIMENSION(4,1) :: UU,CL,CR,UU_L,UU_R,C

    REAL(p2),DIMENSION(:,:)  ,ALLOCATABLE :: U1,U2,U3,U4

    Ln = zero
    Rn = zero
    
    SELECT CASE(Reconstruction_Variables)
    CASE(1)
        CALL Reconstruct_L(i,j,rho,Face,alphaL1,rhoL)
        CALL Reconstruct_L(i,j,u,  Face,alphaL2,uL)
        CALL Reconstruct_L(i,j,v,  Face,alphaL3,vL)
        CALL Reconstruct_L(i,j,p,  Face,alphaL4,pL)

        CALL Reconstruct_R(i,j,rho,Face,alphaR1,rhoR)
        CALL Reconstruct_R(i,j,u,  Face,alphaR2,uR)
        CALL Reconstruct_R(i,j,v,  Face,alphaR3,vR)
        CALL Reconstruct_R(i,j,p,  Face,alphaR4,pR)
        
        UUL(1) = rhoL
        UUL(2) = uL
        UUL(3) = vL
        UUL(4) = pL

        UUR(1) = rhoR
        UUR(2) = uR
        UUR(3) = vR
        UUR(4) = pR
    CASE(2)
        
        ALLOCATE(U1(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
        ALLOCATE(U2(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
        ALLOCATE(U3(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
        ALLOCATE(U4(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))

        SELECT CASE(Face)
        CASE(1)
            rho_L = rho(i,j)
            u_L   = u  (i,j)
            v_L   = v  (i,j)
            p_L   = p  (i,j)

            rho_R = rho(i+1,j)
            u_R   = u  (i+1,j)
            v_R   = v  (i+1,j)
            p_R   = p  (i+1,j)
        CASE(2)
            rho_L = rho(i,j)
            u_L   = u  (i,j)
            v_L   = v  (i,j)
            p_L   = p  (i,j)

            rho_R = rho(i,j+1)
            u_R   = u  (i,j+1)
            v_R   = v  (i,j+1)
            p_R   = p  (i,j+1)
        CASE(3)
            rho_L = rho(i-1,j)
            u_L   = u  (i-1,j)
            v_L   = v  (i-1,j)
            p_L   = p  (i-1,j)

            rho_R = rho(i,j)
            u_R   = u  (i,j)
            v_R   = v  (i,j)
            p_R   = p  (i,j)
        CASE(4)
            rho_L = rho(i,j-1)
            u_L   = u  (i,j-1)
            v_L   = v  (i,j-1)
            p_L   = p  (i,j-1)

            rho_R = rho(i,j)
            u_R   = u  (i,j)
            v_R   = v  (i,j)
            p_R   = p  (i,j)
        END SELECT

        CALL CalculateLw(rho_L,u_L,v_L,p_L,rho_R,u_R,v_R,p_R,Ln,Rn,nx,ny)

        DO n = 1-ghostLayers,jd-1+ghostLayers
            DO m = 1-ghostLayers,id-1+ghostLayers
                UU(1,1) = rho(m,n)
                UU(2,1) = rho(m,n)*u(m,n)
                UU(3,1) = rho(m,n)*v(m,n)
                UU(4,1) = p(m,n)/(gamma-one) + half*rho(m,n)*( u(m,n)*u(m,n)+v(m,n)*v(m,n) )

                C = MATMUL(Ln,UU)
                
                U1(m,n) = C(1,1)
                U2(m,n) = C(2,1)
                U3(m,n) = C(3,1)
                U4(m,n) = C(4,1)
            END DO
        END DO

        CALL Reconstruct_L(i,j,U1,Face,alphaL1,C1L)
        CALL Reconstruct_L(i,j,U2,Face,alphaL2,C2L)
        CALL Reconstruct_L(i,j,U3,Face,alphaL3,C3L)
        CALL Reconstruct_L(i,j,U4,Face,alphaL4,C4L)

        CALL Reconstruct_R(i,j,U1,Face,alphaR1,C1R)
        CALL Reconstruct_R(i,j,U2,Face,alphaR2,C2R)
        CALL Reconstruct_R(i,j,U3,Face,alphaR3,C3R)
        CALL Reconstruct_R(i,j,U4,Face,alphaR4,C4R)
        
        UUL(1) = C1L
        UUL(2) = C2L
        UUL(3) = C3L
        UUL(4) = C4L
        
        UUR(1) = C1R
        UUR(2) = C2R
        UUR(3) = C3R
        UUR(4) = C4R
        
        DEALLOCATE(U1,U2,U3,U4)
        
    END SELECT

    DO n = 1,5
        alphaL(1,n) = alphaL1(n)
        alphaL(2,n) = alphaL2(n)
        alphaL(3,n) = alphaL3(n)
        alphaL(4,n) = alphaL4(n)

        alphaR(1,n) = alphaR1(n)
        alphaR(2,n) = alphaR2(n)
        alphaR(3,n) = alphaR3(n)
        alphaR(4,n) = alphaR4(n)
    END DO

    END SUBROUTINE Reconstruct
    
    
    !**********************************************************************************************
    !******************************************CalculateLw*****************************************
    !**********************************************************************************************
    SUBROUTINE CalculateLw(rhoL_,uL_,vL_,pL_,rhoR_,uR_,vR_,pR_,Ln,Rn,nx,ny)
    IMPLICIT NONE
    
    REAL(p2)::rhoL_,uL_,vL_,pL_,rhoR_,uR_,vR_,pR_,nx,ny
    REAL(p2)::qn,ql,lx,ly,l,q
    REAL(p2)::Ln(4,4),Rn(4,4)
    REAL(p2)::rhoFace,uFace,vFace,hFace,cFace,EL,ER,hL,hR
    
    EL = pL_/(gamma-one)+half*rhoL_*(uL_*uL_+vL_*vL_)
    ER = pR_/(gamma-one)+half*rhoR_*(uR_*uR_+vR_*vR_)

    hL = (gamma/(gamma-one))*pL_/rhoL_+(uL_*uL_+vL_*vL_)/two
    hR = (gamma/(gamma-one))*pR_/rhoR_+(uR_*uR_+vR_*vR_)/two
    
    rhoFace = SQRT(rhoL_*rhoR_)
    hFace   = (SQRT(rhoL_)*hL+SQRT(rhoR_)*hR) / (SQRT(rhoL_)+SQRT(rhoR_))
    uFace   = (SQRT(rhoL_)*uL_+SQRT(rhoR_)*uR_) / (SQRT(rhoL_)+SQRT(rhoR_))
    vFace   = (SQRT(rhoL_)*vL_+SQRT(rhoR_)*vR_) / (SQRT(rhoL_)+SQRT(rhoR_))
    cFace   = SQRT( (gamma-one) *( hFace-half*(uFace*uFace+vFace*vFace) ) )
    
    lx = -ny
    ly = nx
    qn = uFace*nx+vFace*ny
    ql = uFace*lx+vFace*ly
    
    q = sqrt(uFace*uFace+vFace*vFace)
    
    Ln(1,1) = half*( (gamma-one)/two/cFace/cFace*q*q+qn/cFace )
    Ln(1,2) = -half*( (gamma-one)/cFace/cFace*uFace+nx/cFace )
    Ln(1,3) = -half*( (gamma-one)/cFace/cFace*vFace+ny/cFace )
    Ln(1,4) = (gamma-one)/two/cFace/cFace
    Ln(2,1) = one-(gamma-one)/two/cFace/cFace*q*q
    Ln(2,2) = (gamma-one)/cFace/cFace*uFace
    Ln(2,3) = (gamma-one)/cFace/cFace*vFace
    Ln(2,4) = -(gamma-one)/cFace/cFace
    Ln(3,1) = half*((gamma-one)/two/cFace/cFace*q*q-qn/cFace)
    Ln(3,2) = -half*( (gamma-one)/cFace/cFace*uFace-nx/cFace )
    Ln(3,3) = -half*( (gamma-one)/cFace/cFace*vFace-ny/cFace )
    Ln(3,4) = (gamma-one)/two/cFace/cFace
    Ln(4,1) = -ql
    Ln(4,2) = lx
    Ln(4,3) = ly
    Ln(4,4) = zero
    
    Rn(1,1) = one
    Rn(1,2) = one
    Rn(1,3) = one
    Rn(1,4) = zero
    Rn(2,1) = uFace-cFace*nx
    Rn(2,2) = uFace
    Rn(2,3) = uFace+cFace*nx
    Rn(2,4) = lx
    Rn(3,1) = vFace-cFace*ny
    Rn(3,2) = vFace
    Rn(3,3) = vFace+cFace*ny
    Rn(3,4) = ly
    Rn(4,1) = hFace-qn*cFace
    Rn(4,2) = q*q/two
    Rn(4,3) = hFace+qn*cFace
    Rn(4,4) = ql
    
    END SUBROUTINE CalculateLw

    
    
    !**********************************************************************************************
    !*****************************************Reconstruct_L****************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruct_L(i,j,U,Face,alphaL,UL)
    !********************Reconstruct the variables on the left side of the face********************

    IMPLICIT NONE

    INTEGER i,j,Face
    REAL(p2) :: UL
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)
    REAL(p2),DIMENSION(5) :: alphaL

    alphaL = zero

    !*********Reconstruct the variables on the left side of the face by different methods**********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_L(i,j,U,Face,alphaL,UL)
    CASE(2)
        CALL ROUND_L(i,j,U,Face,alphaL,UL)
    CASE(3)
        CALL WENO5_L(i,j,U,Face,alphaL,UL)
    END SELECT

    END SUBROUTINE Reconstruct_L



    !**********************************************************************************************
    !*****************************************Reconstruct_R****************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruct_R(i,j,U,Face,alphaR,UR)
    !********************Reconstruct the variables on the left side of the face********************

    IMPLICIT NONE

    INTEGER i,j,Face
    REAL(p2) :: UR
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)
    REAL(p2),DIMENSION(5) :: alphaR

    alphaR = zero

    !*********Reconstruct the variables on the left side of the face by different methods**********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_R(i,j,U,Face,alphaR,UR)
    CASE(2)
        CALL ROUND_R(i,j,U,Face,alphaR,UR)
    CASE(3)
        CALL WENO5_R(i,j,U,Face,alphaR,UR)
    END SELECT

    END SUBROUTINE Reconstruct_R



    !**********************************************************************************************
    !*******************************************MUSCL_L********************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_L(i,j,U,Face,alphaL,UL)
    !***************Reconstruct the variables on the left side of the face by MUSCL****************

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: UL
    REAL(p2) :: delta_plus,delta_minus,delta,fai
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)
    REAL(p2),DIMENSION(5) :: alphaL

    SELECT CASE(Face)
    CASE(1)
        delta_plus = U(i+1,j)-U(i,j)
        delta_minus= U(i,j)  -U(i-1,j)
    CASE(2)
        delta_plus = U(i,j+1)-U(i,j)
        delta_minus= U(i,j)  -U(i,j-1)
    CASE(3)
        delta_plus = U(i,j)  -U(i-1,j)
        delta_minus= U(i-1,j)-U(i-2,j)
    CASE(4)
        delta_plus = U(i,j)  -U(i,j-1)
        delta_minus= U(i,j-1)-U(i,j-2)
    END SELECT

    CALL Limiter_function(delta_plus,delta_minus,delta,fai)

    alphaL = zero

    SELECT CASE(Face)
    CASE(1)
        UL = U(i,j)+delta
    CASE(2)
        UL = U(i,j)+delta
    CASE(3)
        UL = U(i-1,j)+delta
    CASE(4)
        UL = U(i,j-1)+delta
    END SELECT

    alphaL(3) = one+fai
    alphaL(2) = -fai

    IF (delta_minus==zero)THEN
        alphaL(3) = one
        alphaL(2) = zero     !If faiL = 0, the reconstruction becomes first-order
    END IF

    END SUBROUTINE MUSCL_L



    !**********************************************************************************************
    !*******************************************MUSCL_R********************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_R(i,j,U,Face,alphaR,UR)
    !***************Reconstruct the variables on the left side of the face by MUSCL****************

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: UR
    REAL(p2) :: delta_plus,delta_minus,delta,fai
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)
    REAL(p2),DIMENSION(5) :: alphaR

    SELECT CASE(Face)
    CASE(1)
        delta_plus = U(i+2,j)-U(i+1,j)
        delta_minus= U(i+1,j)-U(i,j)
    CASE(2)
        delta_plus = U(i,j+2)-U(i,j+1)
        delta_minus= U(i,j+1)-U(i,j)
    CASE(3)
        delta_plus = U(i+1,j)-U(i,j)
        delta_minus= U(i,j)  -U(i-1,j)
    CASE(4)
        delta_plus = U(i,j+1)-U(i,j)
        delta_minus= U(i,j)  -U(i,j-1)
    END SELECT

    CALL Limiter_function(delta_minus,delta_plus,delta,fai)

    alphaR = zero

    SELECT CASE(Face)
    CASE(1)
        UR = U(i+1,j)-delta
    CASE(2)
        UR = U(i,j+1)-delta
    CASE(3)
        UR = U(i,j)-delta
    CASE(4)
        UR = U(i,j)-delta
    END SELECT

    alphaR(3) = one+fai
    alphaR(2) = -fai

    IF (delta_plus==zero)THEN
        alphaR(3) = one
        alphaR(2) = zero     !If faiL = 0, the reconstruction becomes first-order
    END IF

    END SUBROUTINE MUSCL_R



    !**********************************************************************************************
    !*******************************************ROUND_L********************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_L(i,j,U,Face,alphaL,UL)
    !***************Reconstruct the variables on the left side of the face by ROUND****************

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: u1,u2,u3
    REAL(p2) :: faiL,UL
    REAL(p2) :: faiL_bar,UL_bar
    REAL(p2) :: omega_0,omega_1
    REAL(p2) :: gamma_0,gamma_1,lambda_1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: delta

    REAL(p2),DIMENSION(5) :: alphaL

    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    delta = 1.0E-16

    SELECT CASE(Face)
    CASE(1)
        u1 = U(i+1,j)
        u2 = U(i-0,j)
        u3 = U(i-1,j)
    CASE(2)
        u1 = U(i,j+1)
        u2 = U(i,j-0)
        u3 = U(i,j-1)
    CASE(3)
        u1 = U(i-0,j)
        u2 = U(i-1,j)
        u3 = U(i-2,j)
    CASE(4)
        u1 = U(i,j-0)
        u2 = U(i,j-1)
        u3 = U(i,j-2)
    END SELECT

    alphaL = zero

    faiL_bar = (u2-u3)/(u1-u3)

    gamma_0  = 1100.0_p2
    gamma_1  = 800.0_p2
    lambda_1 = 0.15_p2

    omega_0 = one/(one+gamma_0*(faiL_bar-one)**4)**2
    omega_1 = one/(one+gamma_1*(faiL_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega_0 + two*faiL_bar*(one-omega_0)
    Temp2 = two*faiL_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega_1 + (lambda_1*faiL_bar-lambda_1+one)*(one-omega_1)
    Temp4 = lambda_1*faiL_bar-lambda_1+one

    IF ((faiL_bar > zero) .AND. (faiL_bar <= half)) THEN
        IF (Temp1 <= Temp2) THEN
            alphaL(4) = one/3.0_p2*omega_0
            alphaL(3) = two-7.0_p2/6.0_p2*omega_0
            alphaL(2) = 5.0_p2/6.0_p2*omega_0-one
        ELSE
            alphaL(4) = zero
            alphaL(3) = two
            alphaL(2) = -one
        END IF
    ELSEIF ((faiL_bar > half) .AND. (faiL_bar <= one)) THEN
        IF (Temp3 <= Temp4) THEN
            alphaL(4) = one/3.0_p2*omega_1 + (one-omega_1)*(one-lambda_1)
            alphaL(3) = 5.0_p2/6.0_p2*omega_1 + lambda_1*(one-omega_1)
            alphaL(2) = -one/6.0_p2*omega_1
        ELSE
            alphaL(4) = one-lambda_1
            alphaL(3) = lambda_1
            alphaL(2) = zero
        END IF
    ELSE
        alphaL(4) = zero
        alphaL(3) = one
        alphaL(2) = zero
    END IF

    IF (u1-u3==0.0) THEN
        alphaL(4) = zero
        alphaL(3) = one
        alphaL(2) = zero
    END IF

    UL = alphaL(4)*u1+alphaL(3)*u2+alphaL(2)*u3

    END SUBROUTINE ROUND_L



    !**********************************************************************************************
    !*******************************************ROUND_R********************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_R(i,j,U,Face,alphaR,UR)
    !***************Reconstruct the variables on the right side of the face by ROUND***************

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: u1,u2,u3
    REAL(p2) :: faiR,UR
    REAL(p2) :: faiR_bar,UR_bar
    REAL(p2) :: omega_0,omega_1
    REAL(p2) :: gamma_0,gamma_1,lambda_1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: delta

    REAL(p2),DIMENSION(5) :: alphaR

    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    delta = 1.0E-16

    SELECT CASE(Face)
    CASE(1)
        u1 = U(i+0,j)
        u2 = U(i+1,j)
        u3 = U(i+2,j)
    CASE(2)
        u1 = U(i,j+0)
        u2 = U(i,j+1)
        u3 = U(i,j+2)
    CASE(3)
        u1 = U(i-1,j)
        u2 = U(i+0,j)
        u3 = U(i+1,j)
    CASE(4)
        u1 = U(i,j-1)
        u2 = U(i,j+0)
        u3 = U(i,j+1)
    END SELECT

    alphaR = zero

    faiR_bar = (u2-u3)/(u1-u3)

    gamma_0  = 1100.0_p2
    gamma_1  = 800.0_p2
    lambda_1 = 0.15_p2

    omega_0 = one/(one+gamma_0*(faiR_bar-one)**4)**2
    omega_1 = one/(one+gamma_1*(faiR_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega_0 + two*faiR_bar*(one-omega_0)
    Temp2 = two*faiR_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega_1 + (lambda_1*faiR_bar-lambda_1+one)*(one-omega_1)
    Temp4 = lambda_1*faiR_bar-lambda_1+one

    IF ((faiR_bar > zero) .AND. (faiR_bar <= half)) THEN
        IF (Temp1 <= Temp2) THEN
            alphaR(4) = one/3.0_p2*omega_0
            alphaR(3) = two-7.0_p2/6.0_p2*omega_0
            alphaR(2) = 5.0_p2/6.0_p2*omega_0-one
        ELSE
            alphaR(4) = zero
            alphaR(3) = two
            alphaR(2) = -one
        END IF
    ELSEIF ((faiR_bar > half) .AND. (faiR_bar <= one)) THEN
        IF (Temp3 <= Temp4) THEN
            alphaR(4) = one/3.0_p2*omega_1 + (one-omega_1)*(one-lambda_1)
            alphaR(3) = 5.0_p2/6.0_p2*omega_1 + lambda_1*(one-omega_1)
            alphaR(2) = -one/6.0_p2*omega_1
        ELSE
            alphaR(4) = one-lambda_1
            alphaR(3) = lambda_1
            alphaR(2) = zero
        END IF
    ELSE
        alphaR(4) = zero
        alphaR(3) = one
        alphaR(2) = zero
    END IF

    IF (u1-u3==0.0_p2) THEN
        alphaR(4) = zero
        alphaR(3) = one
        alphaR(2) = zero
    END IF

    UR = alphaR(4)*u1+alphaR(3)*u2+alphaR(2)*u3

    END SUBROUTINE ROUND_R


    !**********************************************************************************************
    !*******************************************WENO5_L********************************************
    !**********************************************************************************************
    SUBROUTINE WENO5_L(i,j,U,Face,alphaL,UL)
    !***************Reconstruct the variables on the right side of the face by WENO5L**************
    INTEGER :: i,j,Face
    REAL(p2) :: UL
    REAL(p2),DIMENSION(5) :: alphaL

    REAL(p2) :: C30,C31,C32
    REAL(p2) :: v0,v1,v2
    REAL(p2) :: IS0,IS1,IS2
    REAL(p2) :: tao5,alpha0,alpha1,alpha2
    REAL(p2) :: weight0,weight1,weight2

    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    C30 = 0.1_p2
    C31 = 0.6_p2
    C32 = 0.3_p2

    SELECT CASE(Face)
    CASE(1)
        v0 = one/3.0_p2*U(i-2,j) - 7.0_p2/6.0_p2*U(i-1,j) + 11.0_p2/6.0_p2*U(i,j)
        v1 =-one/6.0_p2*U(i-1,j) + 5.0_p2/6.0_p2*U(i,j)   + one/3.0_p2    *U(i+1,j)
        v2 = one/3.0_p2*U(i,j)   + 5.0_p2/6.0_p2*U(i+1,j) - one/6.0_p2    *U(i+2,j)

        IS0 = 13.0_p2/12.0_p2*( U(i-2,j)-two*U(i-1,j)+       U(i,j) )**2 &
            + one/4.0_p2  *( U(i-2,j)-4.0_p2*U(i-1,j)+3.0_p2*U(i,j) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i-1,j)-two*U(i,j)  +U(i+1,j) )**2 &
            + one/4.0_p2  *( U(i-1,j)                -U(i+1,j) )**2
        IS2 = 13.0_p2/12.0_p2*( U(i,j)  -two*U(i+1,j)+U(i+2,j) )**2 &
            + one/4.0_p2  *( 3.0_p2*U(i,j)-4.0_p2*U(i+1,j)+U(i+2,j) )**2
    CASE(2)
        v0 = one/3.0_p2*U(i,j-2) - 7.0_p2/6.0_p2*U(i,j-1) + 11.0_p2/6.0_p2*U(i,j)
        v1 =-one/6.0_p2*U(i,j-1) + 5.0_p2/6.0_p2*U(i,j)   + one/3.0_p2 *U(i,j+1)
        v2 = one/3.0_p2*U(i,j)   + 5.0_p2/6.0_p2*U(i,j+1) - one/6.0_p2 *U(i,j+2)

        IS0 = 13.0_p2/12.0_p2*( U(i,j-2)-two*U(i,j-1)+       U(i,j) )**2 &
            + one/4.0_p2  *( U(i,j-2)-4.0_p2*U(i,j-1)+3.0_p2*U(i,j) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i,j-1)-two*U(i,j)  +U(i,j+1) )**2 &
            + one/4.0_p2  *( U(i,j-1)                -U(i,j+1))**2
        IS2 = 13.0_p2/12.0_p2*( U(i,j)  -two*U(i,j+1)+U(i,j+2) )**2 &
            + one/4.0_p2  *( 3.0_p2*U(i,j)-4.0_p2*U(i,j+1)+U(i,j+2) )**2
    CASE(3)
        v0 = one/3.0_p2*U(i-3,j) - 7.0_p2/6.0_p2*U(i-2,j) + 11.0_p2/6.0_p2*U(i-1,j)
        v1 =-one/6.0_p2*U(i-2,j) + 5.0_p2/6.0_p2*U(i-1,j) + one/3.0_p2 *U(i,j)
        v2 = one/3.0_p2*U(i-1,j) + 5.0_p2/6.0_p2*U(i,j)   - one/6.0_p2 *U(i+1,j)

        IS0 = 13.0_p2/12.0_p2*( U(i-3,j)-   two*U(i-2,j)+       U(i-1,j) )**2 &
            + one/4.0_p2     *( U(i-3,j)-4.0_p2*U(i-2,j)+3.0_p2*U(i-1,j) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i-2,j)-two*U(i-1,j) +U(i,j) )**2 &
            + one/4.0_p2  *( U(i-2,j)                 -U(i,j) )**2
        IS2 = 13.0_p2/12.0_p2*( U(i-1,j)-two*U(i,j)+U(i+1,j) )**2 &
            + one/4.0_p2  *( 3.0_p2*U(i-1,j)-4.0_p2*U(i,j)+U(i+1,j) )**2
    CASE(4)
        v0 = one/3.0_p2*U(i,j-3) - 7.0_p2/6.0_p2*U(i,j-2) + 11.0_p2/6.0_p2*U(i,j-1)
        v1 =-one/6.0_p2*U(i,j-2) + 5.0_p2/6.0_p2*U(i,j-1) + one/3.0_p2 *U(i,j)
        v2 = one/3.0_p2*U(i,j-1) + 5.0_p2/6.0_p2*U(i,j)   - one/6.0_p2 *U(i,j+1)

        IS0 = 13.0_p2/12.0_p2*( U(i,j-3)-two*U(i,j-2)+       U(i,j-1) )**2 &
            + one/4.0_p2  *( U(i,j-3)-4.0_p2*U(i,j-2)+3.0_p2*U(i,j-1) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i,j-2)-two*U(i,j-1) +U(i,j) )**2 &
            + one/4.0_p2  *( U(i,j-2)                 -U(i,j) )**2
        IS2 = 13.0_p2/12.0_p2*( U(i,j-1)-two*U(i,j)+U(i,j+1) )**2 &
            + one/4.0_p2  *( 3.0_p2*U(i,j-1)-4.0_p2*U(i,j)+U(i,j+1) )**2
    END SELECT

    tao5 = ABS(IS0-IS2)

    alpha0 = C30*(1+tao5/(is0+1.0E-15_p2))
    alpha1 = C31*(1+tao5/(is1+1.0E-15_p2))
    alpha2 = C32*(1+tao5/(is2+1.0E-15_p2))

    weight0 = alpha0/(alpha0+alpha1+alpha2)
    weight1 = alpha1/(alpha0+alpha1+alpha2)
    weight2 = alpha2/(alpha0+alpha1+alpha2)

    alphaL(1) = one/3.0_p2*weight0
    alphaL(2) =-one/6.0_p2*(7.0_p2*weight0+weight1)
    alphaL(3) = one/6.0_p2*(11.0_p2*weight0+5.0_p2*weight1+two*weight2)
    alphaL(4) = one/6.0_p2*(two*weight1+5.0_p2*weight2)
    alphaL(5) =-one/6.0_p2*weight2

    UL = weight0*v0+weight1*v1+weight2*v2

    END SUBROUTINE WENO5_L


    !**********************************************************************************************
    !*******************************************WENO5_R********************************************
    !**********************************************************************************************
    SUBROUTINE WENO5_R(i,j,U,Face,alphaR,UR)
    !***************Reconstruct the variables on the right side of the face by WENO5R**************
    INTEGER :: i,j,Face
    REAL(p2) :: UR
    REAL(p2),DIMENSION(5) :: alphaR

    REAL(p2) :: C30,C31,C32
    REAL(p2) :: v0,v1,v2
    REAL(p2) :: IS0,IS1,IS2
    REAL(p2) :: tao5,alpha0,alpha1,alpha2
    REAL(p2) :: weight0,weight1,weight2

    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    c30 = 0.1_p2
    c31 = 0.6_p2
    c32 = 0.3_p2

    SELECT CASE(Face)
    CASE(1)
        v0 = one/3.0_p2*U(i+3,j) - 7.0_p2/6.0_p2*U(i+2,j) + 11.0_p2/6.0_p2*U(i+1,j)
        v1 =-one/6.0_p2*U(i+2,j) + 5.0_p2/6.0_p2*U(i+1,j) + one/3.0_p2 *U(i,j)
        v2 = one/3.0_p2*U(i+1,j) + 5.0_p2/6.0_p2*U(i,j)   - one/6.0_p2 *U(i-1,j)

        IS0 = 13.0_p2/12.0_p2*( U(i+3,j)-two*U(i+2,j)+       U(i+1,j) )**2 &
            + one/4.0_p2  *( U(i+3,j)-4.0_p2*U(i+2,j)+3.0_p2*U(i+1,j) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i+2,j)-two*U(i+1,j)+U(i,j) )**2 &
            + one/4.0_p2  *( U(i+2,j)                -U(i,j) )**2
        IS2 = 13.0_p2/12.0_p2*( U(i+1,j)       -two*U(i,j)+U(i-1,j) )**2 &
            + one/4.0_p2  *( 3.0_p2*U(i+1,j)-4.0_p2*U(i,j)+U(i-1,j) )**2
    CASE(2)
        v0 = one/3.0_p2*U(i,j+3) - 7.0_p2/6.0_p2*U(i,j+2) + 11.0_p2/6.0_p2*U(i,j+1)
        v1 =-one/6.0_p2*U(i,j+2) + 5.0_p2/6.0_p2*U(i,j+1) + one/3.0_p2 *U(i,j)
        v2 = one/3.0_p2*U(i,j+1) + 5.0_p2/6.0_p2*U(i,j)   - one/6.0_p2 *U(i,j-1)

        IS0 = 13.0_p2/12.0_p2*( U(i,j+3)-two*U(i,j+2)+       U(i,j+1) )**2 &
            + one/4.0_p2  *( U(i,j+3)-4.0_p2*U(i,j+2)+3.0_p2*U(i,j+1) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i,j+2)-two*U(i,j+1)+U(i,j) )**2 &
            + one/4.0_p2  *( U(i,j+2)                -U(i,j) )**2
        IS2 = 13.0_p2/12.0_p2*( U(i,j+1)     -two*U(i,j)+U(i,j-1) )**2 &
            + one/4.0_p2*( 3.0_p2*U(i,j+1)-4.0_p2*U(i,j)+U(i,j-1) )**2
    CASE(3)
        v0 = one/3.0_p2*U(i+2,j) - 7.0_p2/6.0_p2*U(i+1,j) + 11.0_p2/6.0_p2*U(i,j)
        v1 =-one/6.0_p2*U(i+1,j) + 5.0_p2/6.0_p2*U(i,j)   + one/3.0_p2    *U(i-1,j)
        v2 = one/3.0_p2*U(i,j)   + 5.0_p2/6.0_p2*U(i-1,j) - one/6.0_p2    *U(i-2,j)

        IS0 = 13.0_p2/12.0_p2*( U(i+2,j)-two*U(i+1,j)+       U(i,j) )**2 &
            + one/4.0_p2  *( U(i+2,j)-4.0_p2*U(i+1,j)+3.0_p2*U(i,j) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i+1,j)-two*U(i,j) +U(i-1,j) )**2 &
            + one/4.0_p2     *( U(i+1,j)            -U(i-1,j) )**2
        IS2 = 13.0_p2/12.0_p2*( U(i,j)     -two*U(i-1,j)+U(i-2,j) )**2 &
            + one/4.0_p2*( 3.0_p2*U(i,j)-4.0_p2*U(i-1,j)+U(i-2,j) )**2
    CASE(4)
        v0 = one/3.0_p2*U(i,j+2) - 7.0_p2/6.0_p2*U(i,j+1) + 11.0_p2/6.0_p2*U(i,j)
        v1 =-one/6.0_p2*U(i,j+1) + 5.0_p2/6.0_p2*U(i,j)   + one/3.0_p2 *U(i,j-1)
        v2 = one/3.0_p2*U(i,j)   + 5.0_p2/6.0_p2*U(i,j-1) - one/6.0_p2 *U(i,j-2)

        IS0 = 13.0_p2/12.0_p2*( U(i,j+2)-two*U(i,j+1)+       U(i,j) )**2 &
            + one/4.0_p2  *( U(i,j+2)-4.0_p2*U(i,j+1)+3.0_p2*U(i,j) )**2
        IS1 = 13.0_p2/12.0_p2*( U(i,j+1)-two*U(i,j) +U(i,j-1) )**2 &
            + one/4.0_p2     *( U(i,j+1)            -U(i,j-1) )**2
        IS2 = 13.0_p2/12.0_p2*( U(i,j)          -two*U(i,j-1)+U(i,j-2) )**2 &
            + one/4.0_p2*     ( 3.0_p2*U(i,j)-4.0_p2*U(i,j-1)+U(i,j-2) )**2
    END SELECT

    tao5 = ABS(IS0-IS2)

    alpha0 = C30*(1+tao5/(is0+1.0E-15_p2))
    alpha1 = C31*(1+tao5/(is1+1.0E-15_p2))
    alpha2 = C32*(1+tao5/(is2+1.0E-15_p2))

    weight0 = alpha0/(alpha0+alpha1+alpha2)
    weight1 = alpha1/(alpha0+alpha1+alpha2)
    weight2 = alpha2/(alpha0+alpha1+alpha2)

    alphaR(1) = one/3.0_p2*weight0
    alphaR(2) =-one/6.0_p2*(7.0_p2*weight0+weight1)
    alphaR(3) = one/6.0_p2*(11.0_p2*weight0+5.0_p2*weight1+two*weight2)
    alphaR(4) = one/6.0_p2*(two*weight1+5.0_p2*weight2)
    alphaR(5) =-one/6.0_p2*weight2

    UR = weight0*v0+weight1*v1+weight2*v2

    END SUBROUTINE WENO5_R

    END MODULE Reconstruction