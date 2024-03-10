    MODULE Initialization
    !**********************************************************************************************
    !****This Module is used to initialized the initial flow field, containg three subroutine: ****
    !****Rankine_Hugoniot，Computation_1D, and Read_File*******************************************
    !**********************************************************************************************
    USE Common_Data,ONLY:rho,u,v,p,id,jd,epsilon,Mainf,zero,one,two,id,jd,ghostLayers,p2,gamma

    IMPLICIT NONE

    CONTAINS

    !**********************************************************************************************
    !**************************************Rankine_Hugoniot****************************************
    !**********************************************************************************************
    SUBROUTINE Rankine_Hugoniot()
    !******This subroutine uses the Rankine Hugoniot condition to initialize the initial flow******
    !************************Only used in the 2D steady normal shock problem***********************

    IMPLICIT NONE

    INTEGER :: i,j,ii

    REAL(p2) :: rhoL,uL,vL,pL          !State of upstream
    REAL(p2) :: rhoR,uR,vR,pR          !State of downstream
    REAL(p2) :: rhoM,uM,vM,pM          !State of the numerical structure
    REAL(p2) :: epsilonU,epsilonP

    !*******************************Define the state of the upstream*******************************
    rhoL = one
    uL   = one
    vL   = zero
    pL   = one/gamma/mainf/mainf
    
    !*********************Calculate the state of the downstream by shock theory********************
    rhoR = one/(two/( (gamma+one)*mainf*mainf )+(gamma-one)/(gamma+one))
    uR   = two/((gamma+one)*mainf*mainf)+(gamma-one)/(gamma+one)
    vR   = zero
    pR   = (two*gamma*mainf*mainf/(gamma+one)-(gamma-one)/(gamma+one))/gamma/mainf/mainf

    !*****Calculate the state of the numerical shock structure by Rankine Hugoniot conditions******
    epsilonU = one-(one-epsilon)/SQRT( one+epsilon*(mainf*mainf-one) / (one+(gamma-one)*mainf*mainf/two) )&
                                /SQRT( one+epsilon*(mainf*mainf-one) / ( one-two*gamma*mainf*mainf/(gamma-one) ) )
    epsilonP = epsilon/SQRT(one+(one-epsilon)*(gamma+one)/(gamma-one)*(mainf*mainf-one)/mainf/mainf)

    rhoM = (one-epsilon)*rhoL+epsilon*rhoR
    uM   = (one-epsilonU)*uL+epsilonU*uR
    vM   = zero
    pM   = (one-epsilonP)*pL+epsilonP*pR

    !**********************************Initialize the initial flow*********************************
    ii = FLOOR(0.5*id)               !Location of the numerical shock structure
    
    DO i = 1-ghostLayers,ii-1
        DO j = 1-ghostLayers,jd-1+ghostLayers
            rho(i,j) = rhoL
            u(i,j)   = uL
            v(i,j)   = vL
            p(i,j)   = pL
        END DO
    END DO

    DO j = 1-ghostLayers,jd-1+ghostLayers
        rho(ii,j) = rhoM
        u(ii,j)   = uM
        v(ii,j)   = vM
        p(ii,j)   = pM
    END DO

    DO i = ii+1,id-1+ghostLayers
        DO j = 1-ghostLayers,jd-1+ghostLayers
            rho(i,j) = rhoR
            u(i,j)   = uR
            v(i,j)   = vR
            p(i,j)   = pR
        END DO
    END DO

    END SUBROUTINE Rankine_Hugoniot

    
    
    !**********************************************************************************************
    !****************************************Computation_1D****************************************
    !**********************************************************************************************
    SUBROUTINE Computation_1D()
    !***This subroutine initializes the initial flow by projecting the 1D result to the 2D plane***
    !************************Only used in the 2D steady normal shock problem***********************
    USE Common_Data,ONLY:rho_1D,u_1D,p_1D
    USE Shock_1D,ONLY:OneDShock

    IMPLICIT NONE

    INTEGER :: i,j,ii
    REAL(p2) :: rhoL,uL,vL,pL
    REAL(p2) :: rhoR,uR,vR,pR

    CALL OneDShock()

    !*******************************Define the state of the upstream*******************************
    rhoL = one
    uL   = one
    vL   = zero
    pL   = one/gamma/mainf/mainf
    
    !*********************Calculate the state of the downstream by shock theory********************
    rhoR = one/(two/( (gamma+one)*mainf*mainf )+(gamma-one)/(gamma+one))
    uR   = two/((gamma+one)*mainf*mainf)+(gamma-one)/(gamma+one)
    vR   = zero
    pR   = (two*gamma*mainf*mainf/(gamma+one)-(gamma-one)/(gamma+one))/gamma/mainf/mainf

    DO i=1,id-1
        rho(i,1)= rho_1D(i)
        u(i,1)  = u_1D(i)
        p(i,1)  = p_1D(i)
    END DO

    !*****************************Reduce the effects of perturbations******************************
    ii=FLOOR(0.5*id)
    DO i=1,ii-1
        IF (ABS(rho(i,1)-rhoL)<1.0E-7) THEN
            rho(i,1)=rhoL
        END IF
        IF (ABS(u(i,1)-uL)<1.0E-7) THEN
            u(i,1)=uL
        END IF
        IF (ABS(p(i,1)-pL)<1.0E-7) THEN
            p(i,1)=pL
        END IF
    END DO
    DO i=ii+1,id-1
        IF (ABS(rho(i,1)-rhoR)<1.0E-7) THEN
            rho(i,1)=rhoR
        END IF
        IF (ABS(u(i,1)-uR)<1.0E-7) THEN
            u(i,1)=uR
        END IF
        IF (ABS(p(i,1)-pR)<1.0E-7) THEN
            p(i,1)=pR
        END IF
    END DO
    
    !**********************************Project onto the 2D plane***********************************
    DO i=1-ghostLayers,0
        rho(i,1)=rho(1,1)
        u(i,1)=u(1,1)
        p(i,1)=p(1,1)
    END DO
    DO i=id,id-1+ghostLayers
        rho(i,1)=rho(id-1,1)
        u(i,1)=u(id-1,1)
        p(i,1)=p(id-1,1)
    END DO
    DO j=1-ghostLayers,jd-1+ghostLayers
        DO i=1-ghostLayers,id-1+ghostLayers
            rho(i,j)=rho(i,1)
            u(i,j)=  u(i,1)
            p(i,j)=  p(i,1)
            v(i,j)=  zero
        END DO
    END DO

    DEALLOCATE(rho_1D,u_1D,p_1D)

    END SUBROUTINE Computation_1D

    
    
    !**********************************************************************************************
    !*******************************************Read_File******************************************
    !**********************************************************************************************
    SUBROUTINE Read_File()
    !*********This subroutine initialize the initial flow by reading the initial flow file*********
    USE Common_Data,ONLY:rho,u,v,p

    IMPLICIT NONE

    INTEGER :: i,j

    !*************************Open the initial flow file and read the data*************************
    OPEN(113,FILE="InitialFlow_rho.dat")
    OPEN(114,FILE="InitialFlow_u.dat")
    OPEN(115,FILE="InitialFlow_v.dat")
    OPEN(116,FILE="InitialFlow_p.dat")
    DO i = 1,id-1
        DO j = 1,jd-1
            READ(113,*)rho(i,j)
            READ(114,*)u(i,j)
            READ(115,*)v(i,j)
            READ(116,*)p(i,j)
        END DO
    END DO

    !************************************Assign the ghost cells************************************
    DO i = 1,id-1
        DO j = 1-ghostLayers,0
            rho(i,j) = rho(i,1)
            u(i,j)   = u(i,1)
            v(i,j)   = v(i,1)
            p(i,j)   = p(i,1)
        END DO
    END DO
    DO i = 1,id-1
        DO j = jd,jd+ghostLayers-1
            rho(i,j) = rho(i,jd-1)
            u(i,j)   = u(i,jd-1)
            v(i,j)   = v(i,jd-1)
            p(i,j)   = p(i,jd-1)
        END DO
    END DO
    DO i = 1-ghostLayers,0
        DO j = 1,jd-1
            rho(i,j) = rho(1,j)
            u(i,j)   = u(1,j)
            v(i,j)   = v(1,j)
            p(i,j)   = p(1,j)
        END DO
    END DO
    DO i = id,id+ghostLayers-1
        DO j = 1,jd-1
            rho(i,j) = rho(id-1,j)
            u(i,j)   = u(id-1,j)
            v(i,j)   = v(id-1,j)
            p(i,j)   = p(id-1,j)
        END DO
    END DO

    END SUBROUTINE Read_File

    END MODULE Initialization