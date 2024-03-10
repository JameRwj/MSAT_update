    MODULE Riemann_Solver_1D
    !**********************************************************************************************
    !This Module contains 8 Riemann solvers: Roe, HLLC, HLL, van Leer, AUSM+, SLAU, HLLE, and HLLEM
    !**********************************************************************************************

    CONTAINS

    !**********************************************************************************************
    !********************************************Roe_1D********************************************
    !**********************************************************************************************
    SUBROUTINE Roe_1D()
    !***************************This subroutine contains the Roe solver****************************
    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,EL               !Variables on the left side of the face
    REAL(p2):: rhoR,uR,pR,hR,cR,ER               !Variables on the right side of the face
    REAL(p2):: rhoFace,cFace,uFace,hFace         !Roe average
    REAL(p2):: sL,sR                             !Wave speed
    REAL(p2):: alfa1,alfa2,alfa3                 !Wave strength
    REAL(p2):: lamda1,lamda2,lamda3              !Averaged eigenvalues

    REAL(p2):: xFlux(3,id)                       !Numerical flux
    REAL(p2):: R1(3),R2(3),R3(3),qL(3),qR(3),FL(3),FR(3)
    REAL(p2):: variableL(3),variableR(3)
    
    DO i=1,id

        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        EL = pL/(gamma-one)+half*rhoL*uL*uL
        ER = pR/(gamma-one)+half*rhoR*uR*uR

        hL = (gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR = (gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        cL = SQRT(gamma*pL/rhoL)
        cR = SQRT(gamma*pR/rhoR)

        rhoFace = sqrt(rhoL*rhoR)
        hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
        uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
        cFace = SQRT((gamma-one)*(hFace-half*uFace**2))

        alfa1 = (pR-pL-rhoFace*cFace*(uR-uL))/two/cFace/cFace
        alfa2 = rhoR-rhoL-(pR-pL)/cFace/cFace
        alfa3 = (pR-pL+rhoFace*cFace*(uR-uL))/two/cFace/cFace

        qL(1) = rhoL
        qL(2) = rhoL*uL
        qL(3) = EL

        qR(1) = rhoR
        qR(2) = rhoR*uR
        qR(3) = ER

        R1(1) = one
        R1(2) = uFace - cFace
        R1(3) = hFace - uFace*cFace

        R2(1)=  one
        R2(2)=  uFace
        R2(3)=  half*uFace*uFace

        R3(1) = one
        R3(2) = uFace + cFace
        R3(3) = hFace + uFace*cFace

        FL(1) = rhoL*uL
        FL(2) = rhoL*uL*uL+pL
        FL(3) = uL*(EL+pL)

        FR(1) = rhoR*uR
        FR(2) = rhoR*uR*uR+pR
        FR(3) = uR*(ER+pR)

        lamda1 = ABS(uFace-cFace)
        lamda2 = ABS(uFace)
        lamda3 = ABS(uFace+cFace)

        xFlux(:,i) = half*(FL(:)+FR(:)-lamda1*alfa1*R1(:)-lamda2*alfa2*R2(:)-lamda3*alfa3*R3(:))

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i=1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    END SUBROUTINE Roe_1D



    !**********************************************************************************************
    !********************************************HLLC_1D*******************************************
    !**********************************************************************************************
    SUBROUTINE HLLC_1D()
    !***************************This subroutine contains the HLLC solver***************************
    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: rhoL,uL,pL,cL,eL,hL                         !Variables on the left side of the face
    REAL(p2):: rhoR,uR,pR,cR,eR,hR                         !Variables on the right side of the face
    REAL(p2):: rhoFace,uFace,cFace,hFace                   !Roe average
    REAL(p2):: sL,sR,sStar                                 !Wave speed
    REAL(p2):: qL(3),qR(3),qStarL(3),qStarR(3)
    REAL(p2):: fluxL(3),fluxR(3),fluxStarL(3),fluxStarR(3)
    REAL(p2):: xFlux(3,id)
    REAL(p2):: variableL(3),variableR(3)

    DO i=1,id
        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        cL = SQRT(gamma*pL/rhoL)
        cR = SQRT(gamma*pR/rhoR)

        eL = pL/(gamma-one)+half*rhoL*uL*uL
        eR = pR/(gamma-one)+half*rhoR*uR*uR

        hL = (EL+pL)/rhoL
        hR = (ER+pR)/rhoR

        hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
        uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
        cFace = SQRT((gamma-one)*(hFace-half*uFace*uFace))

        qL(1) = rhoL
        qL(2) = rhoL*uL
        qL(3) = eL

        qR(1) = rhoR
        qR(2) = rhoR*uR
        qR(3) = eR

        !sL = min(uL-cL,uFace-cFace)
        !sR = max(uR+cR,uFace+cFace)

        sL = min(uL-cL,uR-cR)
        sR = max(uL+cL,uR+cR)

        sStar = ( pR-pL+rhoL*uL*(sL-uL)-rhoR*uR*(sR-uR) )/( rhoL*(sL-uL)-rhoR*(sR-uR) )

        fluxL(1) = rhoL*uL
        fluxL(2) = rhoL*uL*uL+pL
        fluxL(3) = rhoL*uL*hL

        fluxR(1) = rhoR*uR
        fluxR(2) = rhoR*uR*uR+pR
        fluxR(3) = rhoR*uR*hR

        qStarL(1) = rhoL*(sL-uL)/(sL-sStar)
        qStarL(2) = rhoL*(sL-uL)/(sL-sStar)*sStar
        qStarL(3) = rhoL*(sL-uL)/(sL-sStar)*(eL/rhoL+(sStar-uL)*(sStar+pL/rhoL/(sL-uL)))

        qStarR(1) = rhoR*(sR-uR)/(sR-sStar)
        qStarR(2) = rhoR*(sR-uR)/(sR-sStar)*sStar
        qStarR(3) = rhoR*(sR-uR)/(sR-sStar)*(eR/rhoR+(sStar-uR)*(sStar+pR/rhoR/(sR-uR)))

        fluxStarL(:) = fluxL(:)+sL*(qStarL(:)-qL(:))
        fluxStarR(:) = fluxR(:)+sR*(qStarR(:)-qR(:))

        IF(sL >= 0)THEN
            xFlux(:,i) = fluxL(:)
        ELSEIF(sL <= zero .AND. sStar >= zero)THEN
            xFlux(:,i)=fluxStarL(:)
        ELSEIF(sStar <= zero .AND. sR >= zero)THEN
            xFlux(:,i) = fluxStarR(:)
        ELSEIF(sR <= zero)THEN
            xFlux(:,i) = fluxR(:)
        ENDIF

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i=1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    end subroutine HLLC_1D



    !**********************************************************************************************
    !********************************************HLL_1D********************************************
    !**********************************************************************************************
    SUBROUTINE HLL_1D()
    !***************************This subroutine contains the HLL solver****************************
    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,EL         !Variables on the left side of the face
    REAL(p2):: rhoR,uR,pR,hR,cR,ER         !Variables on the right side of the face
    REAL(p2):: sL,sR                       !Wave speed
    REAL(p2):: xFlux(3,id)                 !Numerical flux
    REAL(p2):: qL(3),qR(3)                 !Conservative variables on both sides of the face
    REAL(p2):: FL(3),FR(3)                 !Numerical flux on both sides of the face
    REAL(p2):: variableL(3),variableR(3)

    DO i=1,id
        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        EL = pL/(gamma-one)+half*rhoL*uL*uL
        ER = pR/(gamma-one)+half*rhoR*uR*uR

        hL = (gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR = (gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        cL = SQRT(gamma*pL/rhoL)
        cR = SQRT(gamma*pR/rhoR)

        qL(1) = rhoL
        qL(2) = rhoL*uL
        qL(3) = EL

        qR(1) = rhoR
        qR(2) = rhoR*uR
        qR(3) = ER

        FL(1) = rhoL*uL
        FL(2) = rhoL*uL*uL+pL
        FL(3) = uL*(EL+pL)

        FR(1) = rhoR*uR
        FR(2) = rhoR*uR*uR+pR
        FR(3) = uR*(ER+pR)

        sR = MAX(zero,uL+cL,uR+cR)
        sL = MIN(zero,uR-cR,uL-cL)

        xFlux(:,i) = sR/(sR-sL)*FL-sL/(sR-sL)*FR+sR*sL/(sR-sL)*(qR-qL)

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i=1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    END SUBROUTINE HLL_1D



    !**********************************************************************************************
    !******************************************van_Leer_1D*****************************************
    !**********************************************************************************************
    SUBROUTINE van_Leer_1D()
    !*************************This subroutine contains the van Leer solver*************************
    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,maL,maPlus
    REAL(p2):: rhoR,uR,pR,hR,cR,maR,maMinus
    REAL(p2):: maN
    REAL(p2):: fluxL(3),fluxR(3)
    REAL(p2):: xFlux(3,id)
    REAL(p2):: variableL(3),variableR(3)

    DO i=1,id
        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        cL  = SQRT(gamma*pL/rhoL)
        maL = uL/cL
        IF(maL .GE. one)THEN
            maPlus = maL
        ELSE IF(maL .LE. -one)THEN
            maPlus = zero
        ELSE
            maPlus = 0.25_p2*(maL+one)**2
        END IF

        cR  = SQRT(gamma*pR/rhoR)
        maR = uR/cR
        IF(maR .GE. one)THEN
            maMinus = zero
        ELSE IF(maR .LE. -one)THEN
            maMinus = maR
        ELSE
            maMinus = -0.25_p2*(maR-one)**2
        END IF

        maN = maPlus+maMinus

        IF(maL .GE. one)THEN
            fluxL(1) = rhoL*uL
            fluxL(2) = rhoL*uL*uL+pL
            fluxL(3) = (pL*gamma/(gamma-one)+half*rhoL*(uL*uL))*uL
        ELSE IF(maL .LE. -one)THEN
            fluxL(1) = zero
            fluxL(2) = zero
            fluxL(3) = zero
        ELSE
            fluxL(1) = 0.25_p2*rhoL*cL*(MaL+one)*(MaL+one)
            fluxL(2) = fluxL(1)*((-uL+two*cL)/gamma+uL)
            fluxL(3) = fluxL(1)*half*(((gamma-one)*uL+two*cL)**2/(gamma*gamma-one))
        END IF

        IF(maR .GE. one)THEN
            fluxR(1) = zero
            fluxR(2) = zero
            fluxR(3) = zero
        ELSE IF(maR .LE. -one)THEN
            fluxR(1) = rhoR*uR
            fluxR(2) = rhoR*uR*uR+pR
            fluxR(3) = (pR*gamma/(gamma-one)+half*rhoR*(uR*uR))*uR
        ELSE
            fluxR(1) = -0.25_p2*rhoR*cR*(MaR-one)*(MaR-one)
            fluxR(2) = fluxR(1)*((-uR-two*cR)/gamma+uR)
            fluxR(3) = fluxR(1)*half*(((gamma-one)*uR-two*cR)**2/(gamma*gamma-one))
        END IF

        xFlux(:,i) = fluxL(:)+fluxR(:)

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i=1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    END SUBROUTINE van_Leer_1D



    !**********************************************************************************************
    !******************************************AUSMplus_1D*****************************************
    !**********************************************************************************************
    SUBROUTINE AUSMplus_1D()
    !**************************This subroutine contains the AUSM+ solver***************************

    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,maL,maR
    REAL(p2):: rhoR,uR,pR,hR,cR
    REAL(p2):: cCriticalL,cCriticalR,cFace,maFace
    REAL(p2):: maPlus,maMinus
    REAL(p2):: pPlus,pMinus
    REAL(p2):: sL,sR,del,fp
    REAL(p2):: xFlux(3,id)
    REAL(p2):: variableL(3),variableR(3)

    DO i=1,id

        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        hL = (gamma/(gamma-one))*pL/rhoL+uL*uL/two
        hR = (gamma/(gamma-one))*pR/rhoR+uR*uR/two

        cCriticalL = SQRT(two*(gamma-one)/(gamma+one)*hL)
        cCriticalR = SQRT(two*(gamma-one)/(gamma+one)*hR)

        cL = (cCriticalL*cCriticalL)/MAX(cCriticalL,ABS(uL))
        cR = (cCriticalR*cCriticalR)/MAX(cCriticalR,ABS(uR))

        cFace = MIN(cL,cR)

        maL = uL/cFace
        maR = uR/cFace

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
            xFlux(1,i) = maFace*rhoL*cFace
            xFlux(2,i) = maFace*rhoL*uL*cFace+(pPlus*pL+pMinus*pR)
            xFlux(3,i) = maFace*rhoL*hL*cFace
        ELSE
            xFlux(1,i) = maFace*rhoR*cFace
            xFlux(2,i) = maFace*rhoR*uR*cFace+(pPlus*pL+pMinus*pR)
            xFlux(3,i) = maFace*rhoR*hR*cFace
        END IF

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i = 1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    END SUBROUTINE AUSMplus_1D



    !**********************************************************************************************
    !********************************************SLAU_1D*******************************************
    !**********************************************************************************************
    SUBROUTINE SLAU_1D()
    !**************************This subroutine contains the SLAU solver****************************

    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,EL
    REAL(p2):: rhoR,uR,pR,hR,cR,ER
    REAL(p2):: cFace,delP,delRho,maL,maR,alphaL,alphaR
    REAL(p2):: betaL,betaR,vtFace,Mcap,Xi,Vnabs,VnabsL
    REAL(p2):: VnabsR,fnG,pbar,mass
    REAL(p2):: xFlux(3,id)
    REAL(p2):: flux_plus(3),flux_minus(3)
    REAL(p2):: variableL(3),variableR(3)

    DO i=1,id
        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        EL = pL/(gamma-one)+half*rhoL*uL*uL
        ER = pR/(gamma-one)+half*rhoR*uR*uR

        cL = SQRT(gamma*pL/rhoL)
        cR = SQRT(gamma*pR/rhoR)
        cFace = half*(cL+cR)

        hL = (gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR = (gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        delP = pR-pL
        delRho = rhoR-rhoL

        maL = uL/cFace
        maR = uR/cFace

        alphaL = MAX(zero, one-FLOOR(ABS(maL)))
        alphaR = MAX(zero, one-FLOOR(ABS(maR)))

        betaL = (one-alphaL)*half*(one+SIGN(one,maL)) + (alphaL)*0.25*(two-maL)*((maL+one)**2)
        betaR = (one-alphaR)*half*(one-SIGN(one,maR)) + (alphaR)*0.25*(two+maR)*((maR-one)**2)

        vtface = SQRT(half*((uL*uL)+ (uR*uR)))
        Mcap   = MIN(one, vtface/cFace)
        Xi     = (one - Mcap)**2

        Vnabs = (rhoL *ABS(uL) + rhoR*ABS(uR))/(rhoL + rhoR)

        fnG = -one*MAX(MIN(maL,zero),-one)*MIN(MAX(maR,zero),one)

        pbar = half*((pL+pR) + (betaL-betaR)*(pL-pR) + (one-xi)*(betaL+betaR-one)*(pL+pR))

        VnabsL = (one - fnG)*Vnabs + fnG*ABS(uL)
        VnabsR = (one - fnG)*Vnabs + fnG*ABS(uR)

        mass = half*((rhoL*(uL+VnabsL) + rhoR*(uR-VnabsR)) - (Xi*delp/cFace))

        flux_plus(1) = half*(mass + ABS(mass))
        flux_plus(2) = flux_plus(1)*uL
        flux_plus(3) = flux_plus(1)*hL

        flux_minus(1) = half*(mass - ABS(mass))
        flux_minus(2) = flux_minus(1)*uR
        flux_minus(3) = flux_minus(1)*hR

        xFlux(:,i) = flux_plus+flux_minus
        xFlux(2,i) = xFlux(2,i)+pbar

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i = 1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    END SUBROUTINE SLAU_1D



    !**********************************************************************************************
    !********************************************HLLE_1D*******************************************
    !**********************************************************************************************
    SUBROUTINE HLLE_1D()
    !**************************This subroutine contains the HLLE solver****************************

    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,EL
    REAL(p2):: rhoR,uR,pR,hR,cR,ER
    REAL(p2):: cFace,uFace,hFace
    REAL(p2):: sL,sR
    REAL(p2):: xFlux(3,id)
    REAL(p2):: qL(3),qR(3),FL(3),FR(3)
    REAL(p2):: variableL(3),variableR(3)

    DO i=1,id
        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        EL = pL/(gamma-one)+half*rhoL*uL*uL
        ER = pR/(gamma-one)+half*rhoR*uR*uR

        hL = (gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR = (gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        cL = SQRT(gamma*pL/rhoL)
        cR = SQRT(gamma*pR/rhoR)

        EL = pL/(gamma-one)+half*rhoL*uL*uL
        ER = pR/(gamma-one)+half*rhoR*uR*uR

        hL = (gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR = (gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        cL = SQRT(gamma*pL/rhoL)
        cR = SQRT(gamma*pR/rhoR)

        hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
        uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
        cFace = SQRT((gamma-one)*(hFace-half*uFace*uFace))

        qL(1) = rhoL
        qL(2) = rhoL*uL
        qL(3) = EL

        qR(1) = rhoR
        qR(2) = rhoR*uR
        qR(3) = ER

        FL(1) = rhoL*uL
        FL(2) = rhoL*uL*uL+pL
        FL(3) = uL*(EL+pL)

        FR(1) = rhoR*uR
        FR(2) = rhoR*uR*uR+pR
        FR(3) = uR*(ER+pR)

        sR = MAX(zero,uFace+cFace,uR+cR)
        sL = MIN(zero,uFace-cFace,uL-cL)

        xFlux(:,i) = sR/(sR-sL)*FL-sL/(sR-sL)*FR+sR*sL/(sR-sL)*(qR-qL)

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i = 1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    END SUBROUTINE HLLE_1D



    !**********************************************************************************************
    !*******************************************HLLEM_1D*******************************************
    !**********************************************************************************************
    SUBROUTINE HLLEM_1D()
    !**************************This subroutine contains the HLLEM solver***************************

    USE Common_Data,ONLY:gamma,id,one,two,half,p2,zero,rho_1D,u_1D,p_1D,rhs
    USE Reconstruction_1D

    IMPLICIT NONE

    INTEGER i
    REAL(p2)::rhoL,uL,pL,hL,cL,EL
    REAL(p2)::rhoR,uR,pR,hR,cR,ER
    REAL(p2)::cFace,uFace,hFace
    REAL(p2)::sL,sR
    REAL(p2)::del,alfa2
    REAL(p2):: xFlux(3,id)
    REAL(p2):: R2(3),qL(3),qR(3),FL(3),FR(3)
    REAL(p2):: variableL(3),variableR(3)

    DO i=1,id
        !***************Reconstruct the variables on the both side of every interface**************
        CALL Reconstruct_1D(i,variableL,variableR)

        rhoL = variableL(1)
        uL   = variableL(2)
        pL   = variableL(3)

        rhoR = variableR(1)
        uR   = variableR(2)
        pR   = variableR(3)

        !*******************************Calculate the numerical flux*******************************
        EL=pL/(gamma-one)+half*rhoL*uL*uL
        ER=pR/(gamma-one)+half*rhoR*uR*uR

        hL=(gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR=(gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        cL = SQRT(gamma*pL/rhoL)
        cR = SQRT(gamma*pR/rhoR)

        hFace = (SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
        uFace = (SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
        cFace = SQRT((gamma-one)*(hFace-half*uFace*uFace))

        qL(1) = rhoL
        qL(2) = rhoL*uL
        qL(3) = EL

        qR(1) = rhoR
        qR(2) = rhoR*uR
        qR(3) = ER

        R2(1) = one
        R2(2) = uFace
        R2(3) = half*uFace*uFace

        FL(1) = rhoL*uL
        FL(2) = rhoL*uL*uL+pL
        FL(3) = uL*(EL+pL)

        FR(1) = rhoR*uR
        FR(2) = rhoR*uR*uR+pR
        FR(3) = uR*(ER+pR)

        sR = uFace+cFace
        sL = uFace-cFace

        del = cFace/(ABS(uFace)+cFace)

        alfa2 = rhoR-rhoL-(pR-pL)/cFace/cFace

        IF(sL >= zero)THEN
            xFlux(:,i) = FL
        ELSEIF(sL<zero .AND. sR>zero)THEN
            xFlux(:,i) = sR/(sR-sL)*FL-sL/(sR-sL)*FR+sR*sL/(sR-sL)*(qR-qL-del*alfa2*R2)
        ELSEIF(sR <= zero)THEN
            xFlux(:,i) = FR
        ENDIF

    END DO

    !*********************************Calculate the right-end term*********************************
    DO i = 1,id-1
        rhs(:,i) = rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    END SUBROUTINE HLLEM_1D

    END MODULE Riemann_Solver_1D