    MODULE Shock_1D
    !**********************************************************************************************
    !***This Module performs the 1D computation to obtain the stable 1D flow field, including two**
    !***subroutine: OneDShock and Initialization_1D************************************************
    !**********************************************************************************************
    USE Common_Data,ONLY:p2,Limiter,Riemann_Solver,Initialization_Method,id,x,ghostLayers,&
                         zero,epsilon,x_1D,u_1D,rho_1D,p_1D,deltaT,qq,dqq,rhs,qq_0,residual,&
                         resid,stageCoefficient,DeltaX,one,Mainf,gamma,two,factor,ncount,&
                         MaxCount,dt,nst,half
    IMPLICIT NONE
    
    CONTAINS 
    
    !**********************************************************************************************
    !*******************************************OneDShock******************************************
    !**********************************************************************************************
    SUBROUTINE OneDShock()
    !******************This subroutine is the mian program of the 1D computation*******************
    
    USE Temporal_Discretization_1D
    USE Output
    
    IMPLICIT NONE

    INTEGER :: i,n

    WRITE(*,"(A80)")"========================Performing the 1D computation...========================"
    
    !******************************Allocate memory for 1D variables********************************
    ALLOCATE(x_1D(id),deltaT(id-1),DeltaX(id-1))
    ALLOCATE(rho_1D(1-ghostLayers:id-1+ghostLayers))
    ALLOCATE(u_1D  (1-ghostLayers:id-1+ghostLayers))
    ALLOCATE(p_1D  (1-ghostLayers:id-1+ghostLayers))
    ALLOCATE(qq(3,1:id-1),dqq(3,1:id-1),rhs(3,1:id-1),qq_0(3,1:id-1))
    ALLOCATE(residual(1,3),resid(3))
    
    !**********************************Coefficient of Runge-Kutta**********************************
    stageCoefficient(1)=0.1481_p2
    stageCoefficient(2)=0.4000_p2
    stageCoefficient(3)=1.0000_p2
    
    DO i = 1,id
        x_1D(i) = x(i,1)
    END DO

    DO i = 1,id-1
        DeltaX(i) = x_1D(i+1)-x_1D(i)
    END DO
    
    CALL Initialization_1D()
    
    !*Parameter of controling the modification of Wenjia Xie. If factor = 1, the fix is acting.****
    !*This modification is from "On numerical instabilities of Godunov-type schemes for strong*****
    !*shocks". If there is only one point within the numerical shock structure and epsilon < 0.4,**
    !*the modification is used, since the shock instability will also occour in the 1D computation*
    !*in such conditions.**************************************************************************
    IF (Riemann_Solver == 1 .OR. Riemann_Solver == 7 .OR. Riemann_Solver == 8) THEN
        factor = 1
    ELSE
        factor = 0
    END IF
    
    CALL Boundary()                     !Set boundaries
    
    ncount=0
    DO WHILE (ncount <= maxCount)
        ncount = ncount+1

        CALL CalTimeStep()              !Compute the time step
        CALL Multistage_Runge_Kutta()   !Use multistage Runge Kutta

        !********************Compute the residual of the computation and output********************
        resid=zero
        DO i=1,id-1
            DO n=1,3
                resid(n)=resid(n)+(dqq(n,i)/dt)**2
            END DO
        END DO

        IF (MOD(ncount-1,nst)==0)THEN
            WRITE(*,202) ncount-1
202         FORMAT(10("-"), " Step ", I7, 1X, " >>", 1X,10("-"))

            WRITE(*,203) LOG10(SQRT(SUM(resid(:))))
203         FORMAT(4X,"Residual=",F10.4)
        END IF

        IF(SUM(resid) /= SUM(resid))THEN
            PAUSE
            STOP
        END IF

        IF (MOD(ncount-1,nst)==0)THEN
            IF(ncount == 1)then
                OPEN(112,FILE="Residual_1D.plt",FORM="FORMATTED")
            ELSE
                OPEN(112,FILE="Residual_1D.plt",FORM="FORMATTED",POSITION="APPEND")
            ENDIF

            WRITE(112,"(1X,I8,4(E15.6))") ncount,log10(SQRT(resid(1)*resid(1)))
        END IF
        CLOSE(112)
    END DO
    
    CALL Output_1D()                    !Output the 1D result

    DEALLOCATE(x_1D,deltaT,DeltaX)
    DEALLOCATE(qq,dqq,rhs,qq_0)
    DEALLOCATE(residual,resid)

    WRITE(*,"(A80)")"=========================The 1D computation is finished========================="
    WRITE(*,*)
    
    END SUBROUTINE OneDShock

    
    
    !**********************************************************************************************
    !***************************************Initialization_1D**************************************
    !**********************************************************************************************
    SUBROUTINE Initialization_1D()
    !*******************This subroutine is used to initialize the 1D flow field********************
    IMPLICIT NONE

    INTEGER i,ii

    REAL(p2)::rhoL,rhoR,rhoM,uL,uR,uM,pL,pR,pM
    REAL(p2)::deltaU,deltaP

    !*******************************Define the state of the upstream*******************************
    rhoL = one
    uL   = one
    pL   = one/gamma/mainf/mainf

    !*********************Calculate the state of the downstream by shock theory********************
    rhoR = one/( two/( (gamma+one)*mainf*mainf )+(gamma-one)/(gamma+one) )
    uR   = two/( (gamma+one)*mainf*mainf )+(gamma-one)/(gamma+one)
    pR   = ( two*gamma*mainf*mainf/(gamma+one)-(gamma-one)/(gamma+one) )/gamma/mainf/mainf

    !*****Calculate the state of the numerical shock structure by Rankine Hugoniot conditions******
    deltaU = one-(one-epsilon)/SQRT( one+epsilon*(mainf*mainf-1.0)/(one+(gamma-one)*mainf*mainf/two) )&
                              /SQRT( one+epsilon*(mainf*mainf-one)/( one-two*gamma*mainf*mainf/(gamma-one) ) )
    deltaP = epsilon/SQRT( one + (one-epsilon)*(gamma+one) / (gamma-one)*(mainf*mainf-one)/mainf/mainf )

    rhoM = (one-epsilon)*rhoL+epsilon*rhoR
    uM   = (one-deltaU)*uL+deltaU*uR
    pM   = (one-deltaP)*pL+deltaP*pR
    
    !**********************************Initialize the initial flow*********************************
    ii = FLOOR(0.5*id)
    DO i = 1,ii-1
        rho_1D(i) = rhoL
        u_1D(i)   = uL
        p_1D(i)   = pL
    END DO

    DO i = ii+1,id-1
        rho_1D(i) = rhoR
        u_1D(i)   = uR
        p_1D(i)   = pR
    END DO

    rho_1D(ii) = rhoM
    u_1D(ii)   = uM  
    p_1D(ii)   = pM  
    END SUBROUTINE Initialization_1D
    
    END MODULE Shock_1D