    MODULE Reconstruction_1D
    !**********************************************************************************************
    !***This Module is used to reconstructed the variables on both sides of the face, including:***
    !***Reconstruction_L, Reconstruction_R, MUSCL_L, MUSCL_R***************************************
    !**********************************************************************************************
    CONTAINS

    SUBROUTINE Reconstruct_1D(i,variableL,variableR)
    USE Common_Data ,ONLY:p2,ghostLayers,id,Reconstruction_Variables,rho_1D,u_1D,p_1D
    USE Variable_Conversion_1D
    IMPLICIT NONE

    INTEGER :: i,ii
    REAL(p2):: rhoL,uL,pL
    REAL(p2):: rhoR,uR,pR
    REAL(p2):: C1L,C2L,C3L,C1R,C2R,C3R
    REAL(p2):: rhoL_,uL_,pL_
    REAL(p2):: rhoR_,uR_,pR_
    
    REAL(p2),DIMENSION(3)   :: variableL,variableR
    REAL(p2),DIMENSION(3,3) :: Lw,Rw
    REAL(p2),DIMENSION(3,1) :: W,C,CL,CR,UUL,UUR

    REAL(p2),DIMENSION(:),ALLOCATABLE :: rho_u,rho_e,C1,C2,C3

    ALLOCATE(rho_u(1-ghostLayers:id+ghostLayers),rho_e(1-ghostLayers:id+ghostLayers))
    ALLOCATE(C1(1-ghostLayers:id+ghostLayers),C2(1-ghostLayers:id+ghostLayers),C3(1-ghostLayers:id+ghostLayers))

    SELECT CASE(Reconstruction_Variables)
    CASE(1)
        CALL Reconstruction_L(i,rhoL,rho_1D)
        CALL Reconstruction_L(i,uL,u_1D)
        CALL Reconstruction_L(i,pL,p_1D)
        CALL Reconstruction_R(i,rhoR,rho_1D)
        CALL Reconstruction_R(i,uR,u_1D)
        CALL Reconstruction_R(i,pR,p_1D)
    CASE(2)
        rhoL_ = rho_1D(i-1)
        uL_   = u_1D  (i-1)
        pL_   = p_1D  (i-1)

        rhoR_ = rho_1D(i)
        uR_   = u_1D  (i)
        pR_   = p_1D  (i)
        
        CALL CalculateLw_1D(rhoL_,uL_,pL_,rhoR_,uR_,pR_,Lw,Rw)
        
        DO ii = 1-ghostLayers, id+ghostLayers-1
            W(1,1)=rho_1D(ii)
            W(2,1)=u_1D(ii)
            W(3,1)=p_1D(ii)
            
            C = MATMUL(Lw,W)
            
            C1(ii) = C(1,1)
            C2(ii) = C(2,1)
            C3(ii) = C(3,1)
        END DO
        
        CALL Reconstruction_L(i,C1L,C1)
        CALL Reconstruction_L(i,C2L,C2)
        CALL Reconstruction_L(i,C3L,C3)
        CALL Reconstruction_R(i,C1R,C1)
        CALL Reconstruction_R(i,C2R,C2)
        CALL Reconstruction_R(i,C3R,C3)
        
        CL(1,1)=C1L
        CL(2,1)=C2L
        CL(3,1)=C3L
        
        CR(1,1)=C1R
        CR(2,1)=C2R
        CR(3,1)=C3R
        
        UUL = MATMUL(Rw,CL)
        UUR = MATMUL(Rw,CR)
        
        rhoL = UUL(1,1)
        uL   = UUL(2,1)
        pL   = UUL(3,1)
        
        rhoR = UUR(1,1)
        uR   = UUR(2,1)
        pR   = UUR(3,1)        
    END SELECT

    variableL(1) = rhoL
    variableL(2) = uL
    variableL(3) = pL

    variableR(1) = rhoR
    variableR(2) = uR
    variableR(3) = pR
    
    DEALLOCATE(C1,C2,C3)
    
    END SUBROUTINE Reconstruct_1D
    
    
    
    !**********************************************************************************************
    !***************************************Reconstruction_L***************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruction_L(i,variableL,variable)
    !********************Reconstruct the variables on the left side of the face********************
    USE Common_Data ,ONLY: Reconstruction_Method,p2,ghostLayers,id,Riemann_solver

    IMPLICIT NONE

    INTEGER i
    REAL(p2) :: variableL,variable(1-ghostLayers:id+ghostLayers)

    !*********Reconstruct the variables on the left side of the face by different methods**********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_L_1D(i,variableL,variable)
    CASE(2)
        CALL ROUND_L_1D(i,variableL,variable)
    CASE(3)
        IF (Riemann_Solver == 1 .OR. Riemann_Solver == 7 .OR. Riemann_Solver == 8) THEN
            CALL WENO5_L_1D(i,variableL,variable)
        ELSE
            CALL MUSCL_L_1D(i,variableL,variable)
        END IF
    END SELECT

    END SUBROUTINE Reconstruction_L



    !**********************************************************************************************
    !***************************************Reconstruction_R***************************************
    !**********************************************************************************************
    SUBROUTINE Reconstruction_R(i,variableR,variable)
    !*******************Reconstruct the variables on the right side of the face********************
    USE Common_Data ,ONLY: Reconstruction_Method,p2,ghostLayers,id,Riemann_solver

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableR,variable(1-ghostLayers:id+ghostLayers)

    !*********Reconstruct the variables on the right side of the face by different methods*********
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_R_1D(i,variableR,variable)
    CASE(2)
        CALL ROUND_R_1D(i,variableR,variable)
    CASE(3)
        IF (Riemann_Solver == 1 .OR. Riemann_Solver == 7 .OR. Riemann_Solver == 8) THEN
            CALL WENO5_R_1D(i,variableR,variable)
        ELSE
            CALL MUSCL_R_1D(i,variableR,variable)
        END IF
    END SELECT

    END SUBROUTINE Reconstruction_R



    !**********************************************************************************************
    !******************************************MUSCL_L_1D******************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_L_1D(i,variableL,variable)
    !***************Reconstruct the variables on the left side of the face by MUSCL****************

    USE Common_Data ,ONLY: p2,id,ghostLayers

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2) :: variableL,fai
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta_plus=variable(i)-variable(i-1)
    delta_minus=variable(i-1)-variable(i-2)
    CALL Limiter_function(delta_plus,delta_minus,delta,fai)
    variableL=variable(i-1)+delta

    END SUBROUTINE MUSCL_L_1D



    !**********************************************************************************************
    !******************************************MUSCL_R_1D******************************************
    !**********************************************************************************************
    SUBROUTINE MUSCL_R_1D(i,variableR,variable)
    !**************Reconstruct the variables on the right side of the face by MUSCL****************

    USE Common_Data ,ONLY: p2,id,ghostLayers

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2) :: variableR,fai
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta_plus=variable(i+1)-variable(i)
    delta_minus=variable(i)-variable(i-1)
    CALL Limiter_function(delta_minus,delta_plus,delta,fai)
    variableR=variable(i)-delta

    END SUBROUTINE MUSCL_R_1D



    !**********************************************************************************************
    !******************************************ROUND_L_1D******************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_L_1D(i,variableL,variable)
    !***************Reconstruct the variables on the left side of the face by ROUND****************

    USE Common_Data ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableL,variableL1
    REAL(p2) :: delta
    REAL(p2) :: faiL_bar,variableL_bar
    REAL(p2) :: omega0,omega1
    REAL(p2) :: gamma0,gamma1,lambda1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: alpha1,alpha2,alpha3

    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta = 1.0E-16

    gamma0  = 1100.0_p2
    gamma1  = 800.0_p2
    lambda1 = 0.15_p2

    faiL_bar = (variable(i-1)-variable(i-2))/(variable(i)-variable(i-2))

    omega0 = one/(one+gamma0*(faiL_bar-one)**4)**2
    omega1 = one/(one+gamma1*(faiL_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega0 + two*faiL_bar*(one-omega0)
    Temp2 = two*faiL_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega1 + (lambda1*faiL_bar-lambda1+one)*(one-omega1)
    Temp4 = lambda1*faiL_bar-lambda1+one

    IF ((faiL_bar > zero) .AND. (faiL_bar <= half)) THEN
        variableL_bar = MIN(Temp1,Temp2)
    ELSEIF ((faiL_bar > half) .AND. (faiL_bar <= one)) THEN
        variableL_bar = MIN(Temp3,Temp4)
    ELSE
        variableL_bar = faiL_bar
    END IF

    variableL = variableL_bar*(variable(i)-variable(i-2))+variable(i-2)

    IF (variable(i)-variable(i-2) == 0) THEN
        variableL = variable(i-1)
    END IF

    END SUBROUTINE ROUND_L_1D



    !**********************************************************************************************
    !******************************************ROUND_R_1D******************************************
    !**********************************************************************************************
    SUBROUTINE ROUND_R_1D(i,variableR,variable)
    !***************Reconstruct the variables on the right side of the face by ROUND***************

    USE Common_Data ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableR
    REAL(p2) :: delta
    REAL(p2) :: faiR_bar,variableR_bar
    REAL(p2) :: omega0,omega1
    REAL(p2) :: gamma0,gamma1,lambda1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)
    REAL(p2) :: alpha1,alpha2,alpha3

    delta = 1.0E-16

    gamma0  = 1100.0_p2
    gamma1  = 800.0_p2
    lambda1 = 0.15_p2

    faiR_bar = (variable(i)-variable(i+1))/(variable(i-1)-variable(i+1))

    omega0 = one/(one+gamma0*(faiR_bar-one)**4)**2
    omega1 = one/(one+gamma1*(faiR_bar-one)**4)**2

    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega0 + two*faiR_bar*(one-omega0)
    Temp2 = two*faiR_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega1 + (lambda1*faiR_bar-lambda1+one)*(one-omega1)
    Temp4 = lambda1*faiR_bar-lambda1+one

    IF ((faiR_bar > 0.0) .AND. (faiR_bar <= half)) THEN
        variableR_bar = MIN(Temp1,Temp2)
    ELSEIF ((faiR_bar > half) .AND. (faiR_bar <= one)) THEN
        variableR_bar = MIN(Temp3,Temp4)
    ELSE
        variableR_bar = faiR_bar
    END IF

    variableR = variableR_bar*(variable(i-1)-variable(i+1))+variable(i+1)

    IF (variable(i-1)-variable(i+1)==0)THEN
        variableR = variable(i)
    END IF

    END SUBROUTINE ROUND_R_1D


    !**********************************************************************************************
    !******************************************WENO5_L_1D******************************************
    !**********************************************************************************************
    SUBROUTINE WENO5_L_1D(i,variableL,variable)
    USE Common_Data ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableL0,variableL1,variableL2
    REAL(p2) :: IS0,IS1,IS2
    REAL(p2) :: C30,C31,C32
    REAL(p2) :: alpha0,alpha1,alpha2
    REAL(p2) :: weight0,weight1,weight2
    REAL(p2) :: tao5
    REAL(p2) :: variableL
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    C30=0.1_p2
    C31=0.6_p2
    C32=0.3_p2

    variableL0 =  one/3.0_p2*variable(i-3) -7.0_p2/6.0_p2*variable(i-2) + 11.0_p2/6.0_p2*variable(i-1)
    variableL1 = -one/6.0_p2*variable(i-2) +5.0_p2/6.0_p2*variable(i-1) + one/3.0_p2*variable(i-0)
    variableL2 =  one/3.0_p2*variable(i-1) +5.0_p2/6.0_p2*variable(i-0) - one/6.0_p2*variable(i+1)

    IS0 = 13.0_p2/12.0_p2*(variable(i-3)-two   *variable(i-2)+       variable(i-1))**2 &
        + one/4.0_p2*     (variable(i-3)-4.0_p2*variable(i-2)+3.0_p2*variable(i-1))**2
    IS1 = 13.0_p2/12.0_p2*(variable(i-2)-two*variable(i-1)+variable(i-0))**2 &
        + one/4.0_p2*     (variable(i-2)                  -variable(i-0))**2
    IS2 = 13.0_p2/12.0_p2*(variable(i-1)-two*variable(i-0)+variable(i+1))**2 &
        +one/4.0_p2*      (3.0_p2*variable(i-1)-4.0_p2*variable(i-0)+variable(i+1))**2

    tao5 = ABS(IS0-IS2)

    alpha0=C30*(1+tao5/(IS0+1.0E-15))
    alpha1=C31*(1+tao5/(IS1+1.0E-15))
    alpha2=C32*(1+tao5/(IS2+1.0E-15))

    weight0=alpha0/(alpha0+alpha1+alpha2)
    weight1=alpha1/(alpha0+alpha1+alpha2)
    weight2=alpha2/(alpha0+alpha1+alpha2)

    variableL=weight0*variableL0+weight1*variableL1+weight2*variableL2
    END SUBROUTINE WENO5_L_1D


    !**********************************************************************************************
    !******************************************WENO5_R_1D******************************************
    !**********************************************************************************************
    SUBROUTINE WENO5_R_1D(i,variableR,variable)
    USE Common_Data ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableR0,variableR1,variableR2
    REAL(p2) :: IS0,IS1,IS2
    REAL(p2) :: C30,C31,C32
    REAL(p2) :: alpha0,alpha1,alpha2
    REAL(p2) :: weight0,weight1,weight2
    REAL(p2) :: tao5
    REAL(p2) :: variableR
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    C30=0.1_p2
    C31=0.6_p2
    C32=0.3_p2


    variableR0 =  one/3.0_p2*variable(i+2) -7.0_p2/6.0_p2*variable(i+1) +11.0_p2/6.0_p2*variable(i+0)
    variableR1 = -one/6.0_p2*variable(i+1) +5.0_p2/6.0_p2*variable(i+0) + one/3.0_p2*variable(i-1)
    variableR2 =  one/3.0_p2*variable(i+0) +5.0_p2/6.0_p2*variable(i-1) - one/6.0_p2*variable(i-2)

    IS0 = 13.0_p2/12.0_p2*(variable(i+2) -two  *variable(i+1)+       variable(i+0))**2 &
        + one/4.0_p2     *(variable(i+2)-4.0_p2*variable(i+1)+3.0_p2*variable(i+0))**2
    IS1 = 13.0_p2/12.0_p2*(variable(i+1)-two*variable(i+0)+variable(i-1))**2 &
        + one/4.0_p2     *(variable(i+1)                  -variable(i-1))**2
    IS2 = 13.0_p2/12.0_p2*(variable(i+0)-two*variable(i-1)+variable(i-2))**2 &
        + one/4.0_p2     *(3.0_p2*variable(i+0)-4.0_p2*variable(i-1)+variable(i-2))**2

    tao5 = ABS(IS0-IS2)

    alpha0=C30*(1+tao5/(IS0+1.0E-15))
    alpha1=C31*(1+tao5/(IS1+1.0E-15))
    alpha2=C32*(1+tao5/(IS2+1.0E-15))


    weight0=alpha0/(alpha0+alpha1+alpha2)
    weight1=alpha1/(alpha0+alpha1+alpha2)
    weight2=alpha2/(alpha0+alpha1+alpha2)

    variableR=weight0*variableR0+weight1*variableR1+weight2*variableR2

    END SUBROUTINE WENO5_R_1D

    !**********************************************************************************************
    !***************************************Limiter_function***************************************
    !**********************************************************************************************
    SUBROUTINE Limiter_function(delta_plus,delta_minus,delta,fai)
    !**This subroutine includes the limiter function used in MUSCl, including Superbee, van Leer,**
    !**van Albada, minmod, and the limiter function proposed by Xi Deng****************************

    USE Common_Data,ONLY:p2,half,zero,one,Limiter,two

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2):: r,psi,fai

    r=delta_plus/delta_minus

    IF(ABS(r) .GT. 1.0E+16)THEN
        r=SIGN(one,r)*1.0E+16
    END IF
    IF(delta_minus .EQ. zero)THEN
        r=SIGN(one,delta_plus)*1.0E+16
    END IF
    IF((delta_plus .EQ. zero))THEN
        r=zero
    END IF

    SELECT CASE(Limiter)
    CASE(0)		!no limiter
        delta = zero
        psi   = zero
    CASE(1)		!Superbee
        psi   = MAX(MIN(two*r,one),MIN(r,two))
        psi   = MAX(zero,psi)
        delta = half*psi*delta_minus
    CASE(2)		!van Leer
        psi   = (r+ABS(r))/(one+ABS(r))
        delta = half*psi*delta_minus
    CASE(3)		!van Albada
        psi   = (r*r+r)/(one+r*r)
        delta = half*psi*delta_minus
    CASE(4)		!Minmod
        psi   = MAX(MIN(r,one),zero)
        delta = half*psi*delta_minus
    CASE(5)     !From the paper of Deng Xi in JCP
        IF (r >= zero) THEN
            psi   = (two*r+two*two*r**2)/(one+two*r+3.0_p2*r**2)
        ELSE
            psi = zero
        END IF
        delta = half*psi*delta_minus
    END SELECT
    fai = half*psi

    END SUBROUTINE Limiter_function

    END MODULE Reconstruction_1D