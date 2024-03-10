    MODULE Calculate_Flux
    !**********************************************************************************************
    !****This Module calculate the numerical flux, including subroutines of Calculateflux, Roe,**** 
    !*****HLLC, HLL, van Leer, AUSM+, SLAU, HLLE, and HLLEM****************************************
    !**********************************************************************************************
    USE Common_Data,ONLY:p2,Riemann_solver,zero,half,one,two,gamma,Reconstruction_Variables

    IMPLICIT NONE

    CONTAINS

    
    
    !**********************************************************************************************
    !****************************************Calculateflux*****************************************
    !**********************************************************************************************
    SUBROUTINE Calculateflux(UL,UR,nx,ny,Rn,flux)
    !**********************This subroutine is used to select Riemann solvers***********************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: UL,UR,flux
    REAL(p2) :: nx,ny
    REAL(p2) :: rho_L,u_L,v_L,p_L
    REAL(p2) :: rho_R,u_R,v_R,p_R
    REAL(p2),DIMENSION(4,4) :: Rn
    REAL(p2),DIMENSION(4,1) :: CL,CR,UU_L,UU_R

    IF (Reconstruction_Variables == 2) THEN
        CL(1,1) = UL(1)
        CL(2,1) = UL(2)
        CL(3,1) = UL(3)
        CL(4,1) = UL(4)
        
        CR(1,1) = UR(1)
        CR(2,1) = UR(2)
        CR(3,1) = UR(3)
        CR(4,1) = UR(4)
        
        UU_L = MATMUL(Rn,CL)
        UU_R = MATMUL(Rn,CR)
        
        rho_L = UU_L(1,1)
        u_L   = UU_L(2,1)/UU_L(1,1)
        v_L   = UU_L(3,1)/UU_L(1,1)
        p_L   = (gamma-one) * ( UU_L(4,1)-half*rho_L*( u_L*u_L+v_L*v_L ) )
        
        rho_R = UU_R(1,1)
        u_R   = UU_R(2,1)/UU_R(1,1)
        v_R   = UU_R(3,1)/UU_R(1,1)   
        p_R   = (gamma-one) * (UU_R(4,1)-half*rho_R*( u_R*u_R+v_R*v_R ) )
        
        UL(1) = rho_L
        UL(2) = u_L
        UL(3) = v_L
        UL(4) = p_L

        UR(1) = rho_R
        UR(2) = u_R
        UR(3) = v_R
        UR(4) = p_R
    END IF
    
    SELECT CASE(Riemann_solver)
    CASE(1)
        CALL Roe(UL,UR,nx,ny,flux)
    CASE(2)
        CALL HLLC(UL,UR,nx,ny,flux)
    CASE(3)
        CALL HLL(UL,UR,nx,ny,flux)
    CASE(4)
        CALL van_Leer(UL,UR,nx,ny,flux)
    CASE(5)
        CALL AUSMplus(UL,UR,nx,ny,flux)
    CASE(6)
        CALL SLAU(UL,UR,nx,ny,flux)
    CASE(7)
        CALL HLLE(UL,UR,nx,ny,flux)
    CASE(8)
        CALL HLLEM(UL,UR,nx,ny,flux)

    END SELECT

    END SUBROUTINE Calculateflux

    
    
    !**********************************************************************************************
    !**********************************************Roe*********************************************
    !**********************************************************************************************
    SUBROUTINE Roe(variableL,variableR,nx,ny,flux)
    !***************************This subroutine contains the Roe solver****************************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: R1,R2,R3,R4
    REAL(p2),DIMENSION(4) :: qL,qR,FluxL,FluxR,diss
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL,tnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR,tnR
    REAL(p2) :: rhoFace,cFace,uFace,vFace,hFace,qnFace,tnFace
    REAL(p2) :: lamda1,lamda2,lamda3,lamda4,ee,ee2
    REAL(p2) :: alfa1,alfa2,alfa3,alfa4

    REAL(p2):: nx,ny,tx,ty

    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)

    eL = pL/(gamma-one)+half*rhoL*(uL**2+vL**2)
    eR = pR/(gamma-one)+half*rhoR*(uR**2+vR**2)

    hL = (EL+pL)/rhoL
    hR = (ER+pR)/rhoR

    rhoFace = SQRT(rhoL*rhoR)
    uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
    vFace = (SQRT(rhoL)*vL+SQRT(rhoR)*vR)/(SQRT(rhoL)+SQRT(rhoR))
    hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
    cFace = SQRT( (gamma-one)*(hFace-half*(uFace**2+vFace**2)) )

    tx = -ny
    ty = nx

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny
    qnFace = (SQRT(rhoL)*qnL+SQRT(rhoR)*qnR)/(SQRT(rhoL)+SQRT(rhoR))

    tnL = uL*tx+vL*ty
    tnR = uR*tx+vR*ty
    tnFace = (SQRT(rhoL)*tnL+SQRT(rhoR)*tnR)/(SQRT(rhoL)+SQRT(rhoR))

    lamda1 = abs(qnFace-cFace)
    lamda2 = abs(qnFace)
    lamda3 = abs(qnFace)
    lamda4 = abs(qnFace+cFace)

    !ee = 0.2_p2
    !if( lamda1<ee )  lamda1 = half*(lamda1*lamda1/ee+ee)
    !if( lamda4<ee )  lamda4 = half*(lamda4*lamda4/ee+ee)

    !ee2 = 0.9_p2*cFace
    !if( lamda2<ee2 )  lamda2 = half*(lamda2*lamda2/ee2+ee2)

    !ee2 = 2.0_p2*cFace
    !if( lamda3<ee2 )  lamda3 = half*(lamda3*lamda3/ee2+ee2)

    alfa1 = (pR-pL-rhoFace*cFace*(qnR-qnL))/(2.0_p2*cFace*cFace)
    alfa2 = rhoR-rhoL-(pR-pL)/cFace/cFace
    alfa3 = rhoFace*(tnR-tnL)/cFace
    alfa4 = (pR-pL+rhoFace*cFace*(qnR-qnL))/(2.0_p2*cFace*cFace)

    R1(1) = one
    R1(2) = uFace-cFace*nx
    R1(3) = vFace-cFace*ny
    R1(4) = hFace-qnFace*cFace

    R2(1) = one
    R2(2) = uFace
    R2(3) = vFace
    R2(4) = half*(uFace**2+vFace**2)

    R3(1) = zero
    R3(2) = cFace*tx
    R3(3) = cFace*ty
    R3(4) = cFace*tnFace

    R4(1) = one
    R4(2) = uFace+cFace*nx
    R4(3) = vFace+cFace*ny
    R4(4) = hFace+qnFace*cFace

    diss(:) = lamda1*alfa1*R1(:)+lamda2*alfa2*R2(:)+lamda3*alfa3*R3(:)+lamda4*alfa4*R4(:)

    FluxL(1) = rhoL*qnL
    FluxL(2) = rhoL*qnL*uL+nx*pL
    FluxL(3) = rhoL*qnL*vL+ny*pL
    FluxL(4) = rhoL*qnL*hL

    FluxR(1) = rhoR*qnR
    FluxR(2) = rhoR*qnR*uR+nx*pR
    FluxR(3) = rhoR*qnR*vR+ny*pR
    FluxR(4) = rhoR*qnR*hR

    flux(:) = half*(FluxL(:)+FluxR(:)-diss(:))
    END SUBROUTINE Roe

    
    
    !**********************************************************************************************
    !*********************************************HLLC*********************************************
    !**********************************************************************************************
    SUBROUTINE HLLC(variableL,variableR,nx,ny,flux)
    !***************************This subroutine contains the HLLC solver***************************
    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,FluxStarL,FluxStarR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    REAL(p2) :: cFace,uFace,vFace,hFace,qnFace,rhoFace
    REAL(p2) :: rhoStarL,rhoStarR,eStarL,eStarR
    REAL(p2) :: sL,sR,sStar,pStar
    REAL(p2) :: nx,ny
    REAL(p2) :: alfaL,alfaR


    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)

    eL = pL/(gamma-one)/rhoL+half*(uL**2+vL**2)
    eR = pR/(gamma-one)/rhoR+half*(uR**2+vR**2)

    HL = gamma/(gamma-one)*pL/rhoL+half*(uL**2+vL**2)
    HR = gamma/(gamma-one)*pR/rhoR+half*(uR**2+vR**2)

    rhoFace = SQRT(rhoL*rhoR)
    uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
    vFace = (SQRT(rhoL)*vL+SQRT(rhoR)*vR)/(SQRT(rhoL)+SQRT(rhoR))
    hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
    cFace = SQRT( (gamma-one)*(hFace-half*(uFace**2+vFace**2)) )

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny
    qnFace = (SQRT(rhoL)*qnL+SQRT(rhoR)*qnR)/(SQRT(rhoL)+SQRT(rhoR))

    FluxL(1) = rhoL*qnL
    FluxL(2) = rhoL*qnL*uL+nx*pL
    FluxL(3) = rhoL*qnL*vL+ny*pL
    FluxL(4) = rhoL*qnL*hL

    FluxR(1) = rhoR*qnR
    FluxR(2) = rhoR*qnR*uR+nx*pR
    FluxR(3) = rhoR*qnR*vR+ny*pR
    FluxR(4) = rhoR*qnR*hR

    sL = min(qnR-cR,qnL-cL)
    sR = max(qnL+cL,qnR+cR)

    alfaL = rhoL*(sL-qnL)
    alfaR = rhoR*(sR-qnR)

    sStar = (alfaR*qnR-alfaL*qnL+pL-pR)/(alfaR-alfaL)

    rhoStarL = alfaL/(sL-sStar)
    rhoStarR = alfaR/(sR-sStar)

    eStarL = eL+(sStar-qnL)*(sStar+pL/alfaL)
    eStarR = eR+(sStar-qnR)*(sStar+pR/alfaR)

    pStar = (alfaR*pL-alfaL*pR-alfaL*alfaR*(qnL-qnR))/(alfaR-alfaL)

    FluxStarL(1) = rhoStarL*sStar
    FluxStarL(2) = rhoStarL*sStar*(uL+nx*(sStar-qnL))+nx*pStar
    FluxStarL(3) = rhoStarL*sStar*(vL+ny*(sStar-qnL))+ny*pStar
    FluxStarL(4) = sStar*(rhoStarL*eStarL+pStar)

    FluxStarR(1) = rhoStarR*sStar
    FluxStarR(2) = rhoStarR*sStar*(uR+nx*(sStar-qnR))+nx*pStar
    FluxStarR(3) = rhoStarR*sStar*(vR+ny*(sStar-qnR))+ny*pStar
    FluxStarR(4) = sStar*(rhoStarR*eStarR+pStar)

    IF(sL .GE. zero)THEN
        Flux(:) = FluxL
    ENDIF

    IF( (sL .LE. zero) .and. (sStar .GE. zero))THEN
        Flux(:) = FluxStarL
    ENDIF

    IF( (sStar .LE. zero) .and. (sR .GE. zero))THEN
        Flux(:) = FluxStarR
    ENDIF

    IF( sR .LE. zero)THEN
        Flux(:) = FluxR
    END IF

    END SUBROUTINE HLLC

    
    
    !**********************************************************************************************
    !**********************************************HLL*********************************************
    !**********************************************************************************************
    SUBROUTINE HLL(variableL,variableR,nx,ny,flux)
    !***************************This subroutine contains the HLL solver****************************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,qL,qR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    REAL(p2) :: sL,sR
    REAL(p2) :: nx,ny


    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)

    eL = pL/(gamma-one)/rhoL+half*(uL**2+vL**2)
    eR = pR/(gamma-one)/rhoR+half*(uR**2+vR**2)

    HL = gamma/(gamma-one)*pL/rhoL+half*(uL**2+vL**2)
    HR = gamma/(gamma-one)*pR/rhoR+half*(uR**2+vR**2)

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny

    FluxL(1) = rhoL*qnL
    FluxL(2) = rhoL*qnL*uL+nx*pL
    FluxL(3) = rhoL*qnL*vL+ny*pL
    FluxL(4) = rhoL*qnL*hL

    FluxR(1) = rhoR*qnR
    FluxR(2) = rhoR*qnR*uR+nx*pR
    FluxR(3) = rhoR*qnR*vR+ny*pR
    FluxR(4) = rhoR*qnR*hR

    sL = MIN(qnR-cR,qnL-cL)
    sR = MAX(qnL+cL,qnR+cR)

    qL(1) = rhoL
    qL(2) = rhoL*uL
    qL(3) = rhoL*vL
    qL(4) = rhoL*eL

    qR(1) = rhoR
    qR(2) = rhoR*uR
    qR(3) = rhoR*vR
    qR(4) = rhoR*eR

    IF(sL .GE. zero)THEN
        Flux(:) = FluxL
    ENDIF

    IF( (sL .LE. zero) .and. (sR .GE. zero))THEN
        Flux(:) = sR/(sR-sL)*FluxL(:)-sL/(sR-sL)*FluxR(:)+sR*sL/(sR-sL)*(qR(:)-qL(:))
    ENDIF

    IF( sR .LE. zero)THEN
        Flux(:) = FluxR
    END IF

    END SUBROUTINE HLL



    !**********************************************************************************************
    !********************************************van_Leer******************************************
    !**********************************************************************************************
    SUBROUTINE van_Leer(variableL,variableR,nx,ny,flux)
    !*************************This subroutine contains the van Leer solver*************************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,qL,qR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL,maL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR,maR
    REAL(p2) :: sL,sR
    REAL(p2) :: nx,ny


    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)

    eL = pL/(gamma-one)/rhoL+half*(uL**2+vL**2)
    eR = pR/(gamma-one)/rhoR+half*(uR**2+vR**2)

    HL = gamma/(gamma-one)*pL/rhoL+half*(uL**2+vL**2)
    HR = gamma/(gamma-one)*pR/rhoR+half*(uR**2+vR**2)

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny

    maL = qnL/cL
    maR = qnR/cR

    IF(maL .GE. one)THEN
        fluxL(1) = rhoL*qnL
        fluxL(2) = rhoL*qnL*uL+pL*nx
        fluxL(3) = rhoL*qnL*vL+pL*ny
        fluxL(4) = (pL*gamma/(gamma-one)+rhoL*(uL*uL+vL*vL)*half)*qnL
    ELSE IF(maL .LE. -one)THEN
        fluxL(1) = zero
        fluxL(2) = zero
        fluxL(3) = zero
        fluxL(4) = zero
    ELSE
        fluxL(1) = 0.25_p2*rhoL*cL*(MaL+one)*(MaL+one)
        fluxL(2) = fluxL(1)*(nx*(-qnL+two*cL)/gamma+uL)
        fluxL(3) = fluxL(1)*(ny*(-qnL+two*cL)/gamma+vL)
        fluxL(4) = fluxL(1)*half*(((gamma-one)*qnL+two*cL)**2/(gamma*gamma-one)&
            +(uL*uL+vL*vL-qnL*qnL))
    END IF

    IF(maR .GE. one)THEN
        fluxR(1) = zero
        fluxR(2) = zero
        fluxR(3) = zero
        fluxR(4) = zero
    ELSE IF(maR .LE. -one)THEN
        fluxR(1) = rhoR*qnR
        fluxR(2) = rhoR*qnR*uR+pR*nx
        fluxR(3) = rhoR*qnR*vR+pR*ny
        fluxR(4) = (pR*gamma/(gamma-one)+rhoR*(uR*uR+vR*vR)*half)*qnR
    ELSE
        fluxR(1) = -0.25_p2*rhoR*cR*(MaR-one)*(MaR-one)
        fluxR(2) = fluxR(1)*(nx*(-qnR-two*cR)/gamma+uR)
        fluxR(3) = fluxR(1)*(ny*(-qnR-two*cR)/gamma+vR)
        fluxR(4) = fluxR(1)*half*(((gamma-one)*qnR-two*cR)**2/(gamma*gamma-one)+(uR*uR+vR*vR-qnR*qnR))
    END IF

    Flux(:) = fluxL(:)+fluxR(:)

    END SUBROUTINE van_Leer

    
    
    !**********************************************************************************************
    !********************************************AUSMplus******************************************
    !**********************************************************************************************
    SUBROUTINE AUSMplus(variableL,variableR,nx,ny,flux)
    !**************************This subroutine contains the AUSM+ solver***************************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2) :: rhoL,uL,vL,pL,hL,cL,maL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,hR,cR,maR,qnR
    REAL(p2) :: cCriticalL,cCriticalR,cFace,maFace
    REAL(p2) :: maPlus,maMinus
    REAL(p2) :: pPlus,pMinus
    REAL(p2) ::sL,sR
    REAL(p2) ::nx,ny


    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)

    hL = (gamma/(gamma-one))*pL/rhoL+uL*uL/two+vL*vL/two
    hR = (gamma/(gamma-one))*pR/rhoR+uR*uR/two+vR*vR/two

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny

    cCriticalL = SQRT(two*(gamma-one)/(gamma+one)*hL)
    cCriticalR = SQRT(two*(gamma-one)/(gamma+one)*hR)

    cL = (cCriticalL*cCriticalL)/MAX(cCriticalL,ABS(qnL))
    cR = (cCriticalR*cCriticalR)/MAX(cCriticalR,ABS(qnR))

    cFace = half*(cL+cR)

    maL = qnL/cFace
    maR = qnR/cFace

    IF(ABS(maL) .LE. one)THEN
        maPlus = 0.25_p2*(maL+one)**2+one/8.0_p2*(maL**2-one)**2
    ELSE
        maPlus = half*(maL+ABS(maL))
    END IF

    IF(ABS(maR) .LE. one)THEN
        maMinus = -0.25_p2*(maR-one)**2-one/8.0_p2*(maR**2-one)**2
    ELSE
        maMinus = half*(maR-ABS(maR))
    END IF

    IF(ABS(maL) .LE. one)THEN
        pPlus = 0.25_p2*(maL+one)**2*(two-maL)+3.0_p2/16.0_p2*maL*(maL**2-one)**2
    ELSE
        pPlus = half*(one+SIGN(one,maL))
    END IF

    IF(ABS(maR) .LE. one)THEN
        pMinus = 0.25_p2*(maR-one)**2*(two+maR)-3.0_p2/16.0_p2*maR*(maR**2-one)**2
    ELSE
        pMinus = half*(one-SIGN(one,maR))
    END IF

    maFace = maPlus+maMinus

    IF(maFace .GT. zero)THEN
        Flux(1) = maFace*rhoL*cFace
        Flux(2) = maFace*rhoL*uL*cFace+(pPlus*pL+pMinus*pR)*nx
        Flux(3) = maFace*rhoL*vL*cFace+(pPlus*pL+pMinus*pR)*ny
        Flux(4) = maFace*rhoL*hL*cFace
    ELSE
        Flux(1) = maFace*rhoR*cFace
        Flux(2) = maFace*rhoR*uR*cFace+(pPlus*pL+pMinus*pR)*nx
        Flux(3) = maFace*rhoR*vR*cFace+(pPlus*pL+pMinus*pR)*ny
        Flux(4) = maFace*rhoR*hR*cFace
    END IF

    END SUBROUTINE AUSMplus

    
    
    !**********************************************************************************************
    !**********************************************SLAU********************************************
    !**********************************************************************************************
    SUBROUTINE SLAU(variableL,variableR,nx,ny,flux)
    !**************************This subroutine contains the SLAU solver****************************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,flux_plus,flux_minus
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    REAL(p2) :: cFace,delP,delRho,maL,maR,alphaL,alphaR,betaL,betaR,vtFace,Mcap,Xi,Vnabs,VnabsL,VnabsR,fnG,pBar,mass
    REAL(p2) :: nx,ny

    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)
    cFace = half*(cL+cR)
    hL = (gamma/(gamma-one))*pL/rhoL+(uL*uL+vL*vL)/two
    hR = (gamma/(gamma-one))*pR/rhoR+(uR*uR+vR*vR)/two
    delP = pR-pL
    delRho = rhoR-rhoL

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny

    maL = qnL/cFace
    maR = qnR/cFace

    alphaL = MAX(zero, one-FLOOR(ABS(maL)))
    alphaR = MAX(zero, one-FLOOR(ABS(maR)))

    betaL = (one-alphaL)*half*(one+SIGN(one,maL)) + (alphaL)*0.25_p2*(two-maL)*((maL+one)**2)
    betaR = (one-alphaR)*half*(one-SIGN(one,maR)) + (alphaR)*0.25_p2*(two+maR)*((maR-one)**2)

    vtface = SQRT(half*((uL*uL) + (vL*vL) + (uR*uR) + (vR*vR)))
    Mcap   = MIN(one, vtface/cFace)
    Xi     = (one - Mcap)**2

    Vnabs = (rhoL *ABS(qnL) + rhoR*ABS(qnR))/(rhoL + rhoR)

    fnG = -one*MAX(MIN(maL,zero),-one)*MIN(MAX(maR,zero),one)

    pbar = half*((pL+pR) + (betaL-betaR)*(pL-pR) + (one-xi)*(betaL+betaR-one)*(pL+pR))

    VnabsL = (one - fnG)*Vnabs + fnG*ABS(qnL)
    VnabsR = (one - fnG)*Vnabs + fnG*ABS(qnR)

    mass = half*((rhoL*(qnL+VnabsL) + rhoR*(qnR-VnabsR)) - (Xi*delp/cFace))

    flux_plus(1) = half*(mass + ABS(mass))
    flux_plus(2) = flux_plus(1)*uL
    flux_plus(3) = flux_plus(1)*vL
    flux_plus(4) = flux_plus(1)*hL

    flux_minus(1) = half*(mass - ABS(mass))
    flux_minus(2) = flux_minus(1)*uR
    flux_minus(3) = flux_minus(1)*vR
    flux_minus(4) = flux_minus(1)*hR


    Flux(:) = flux_plus+flux_minus
    Flux(2) = Flux(2)+pbar*nx
    Flux(3) = Flux(3)+pbar*ny

    END SUBROUTINE SLAU

    
    
    !**********************************************************************************************
    !*********************************************HLLE*********************************************
    !**********************************************************************************************
    SUBROUTINE HLLE(variableL,variableR,nx,ny,flux)
    !**************************This subroutine contains the HLLE solver****************************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,qL,qR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    REAL(p2) :: rhoFace,cFace,uFace,vFace,hFace,qnFace
    REAL(p2) :: sL,sR,sStar
    REAL(p2) :: nx,ny


    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)

    EL = pL/(gamma-one)+half*rhoL*(uL*uL+vL*vL)
    ER = pR/(gamma-one)+half*rhoR*(uR*uR+vR*vR)

    hL = (EL+pL)/rhoL
    hR = (ER+pR)/rhoR

    qL(1) = rhoL
    qL(2) = rhoL*uL
    qL(3) = rhoL*vL
    qL(4) = EL

    qR(1) = rhoR
    qR(2) = rhoR*uR
    qR(3) = rhoR*vR
    qR(4) = ER

    rhoFace = SQRT(rhoL*rhoR)
    uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
    vFace = (SQRT(rhoL)*vL+SQRT(rhoR)*vR)/(SQRT(rhoL)+SQRT(rhoR))
    hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
    cFace = SQRT( (gamma-one)*(hFace-half*(uFace*uFace+vFace*vFace)) )

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny
    qnFace = (SQRT(rhoL)*qnL+SQRT(rhoR)*qnR)/(SQRT(rhoL)+SQRT(rhoR))

    sL = MIN(zero,qnL-cL,qnFace-cFace)
    sR = MAX(zero,qnR+cR,qnFace+cFace)

    FluxL(1) = rhoL*qnL
    FluxL(2) = rhoL*qnL*uL+nx*pL
    FluxL(3) = rhoL*qnL*vL+ny*pL
    FluxL(4) = rhoL*qnL*hL

    FluxR(1) = rhoR*qnR
    FluxR(2) = rhoR*qnR*uR+nx*pR
    FluxR(3) = rhoR*qnR*vR+ny*pR
    FluxR(4) = rhoR*qnR*hR

    Flux(:) = sR*FluxL/(sR-sL)-sL*FluxR/(sR-sL)+sL*sR*(qR-qL)/(sR-sL)

    END SUBROUTINE HLLE

    
    
    !**********************************************************************************************
    !*********************************************HLLEM********************************************
    !**********************************************************************************************
    SUBROUTINE HLLEM(variableL,variableR,nx,ny,flux)
    !**************************This subroutine contains the HLLEM solver***************************

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,qL,qR,R1,R2
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    REAL(p2) :: rhoFace,cFace,uFace,vFace,hFace,qnFace
    REAL(p2) :: sL,sR
    REAL(p2) :: alfa1,alfa2
    REAL(p2) :: Del
    REAL(p2) :: nx,ny

    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL = SQRT(gamma*pL/rhoL)
    cR = SQRT(gamma*pR/rhoR)

    EL = pL/(gamma-one)+half*rhoL*(uL*uL+vL*vL)
    ER = pR/(gamma-one)+half*rhoR*(uR*uR+vR*vR)

    hL = (EL+pL)/rhoL
    hR = (ER+pR)/rhoR

    qL(1) = rhoL
    qL(2) = rhoL*uL
    qL(3) = rhoL*vL
    qL(4) = EL

    qR(1) = rhoR
    qR(2) = rhoR*uR
    qR(3) = rhoR*vR
    qR(4) = ER

    rhoFace = SQRT(rhoL*rhoR)
    uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
    vFace = (SQRT(rhoL)*vL+SQRT(rhoR)*vR)/(SQRT(rhoL)+SQRT(rhoR))
    hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
    cFace = SQRT( (gamma-one)*(hFace-half*(uFace*uFace+vFace*vFace)) )

    qnL = uL*nx+vL*ny
    qnR = uR*nx+vR*ny
    qnFace = (SQRT(rhoL)*qnL+SQRT(rhoR)*qnR)/(SQRT(rhoL)+SQRT(rhoR))

    sL = MIN(zero,qnL-cL,qnFace-cFace)
    sR = MAX(zero,qnR+cR,qnFace+cFace)

    Del = cFace/(ABS(qnFace)+cFace)

    alfa1 = rhoR-rhoL-(pR-pL)/cFace/cFace
    alfa2 = rhoFace

    R1(1) = one
    R1(2) = uFace
    R1(3) = vFace
    R1(4) = half*(uFace*uFace+vFace*vFace)

    R2(1) = zero
    R2(2) = uR-uL-(qnR-qnL)*nx
    R2(3) = vR-vL-(qnR-qnL)*ny
    R2(4) = uFace*(uR-uL)+vFace*(vR-vL)-qnFace*(qnR-qnL)

    FluxL(1) = rhoL*qnL
    FluxL(2) = rhoL*qnL*uL+nx*pL
    FluxL(3) = rhoL*qnL*vL+ny*pL
    FluxL(4) = rhoL*qnL*hL

    FluxR(1) = rhoR*qnR
    FluxR(2) = rhoR*qnR*uR+nx*pR
    FluxR(3) = rhoR*qnR*vR+ny*pR
    FluxR(4) = rhoR*qnR*hR

    Flux(:) = sR*FluxL/(sR-sL)-sL*FluxR/(sR-sL)+sL*sR*(qR-qL-Del*alfa1*R1-Del*alfa2*R2)/(sR-sL)

    END SUBROUTINE HLLEM

    END MODULE