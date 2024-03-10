    MODULE Variable_Conversion_1D
    !**********************************************************************************************
    !*This Module performs the conversion between conservative and primitive variables, including**
    !*Primitive_to_Conservative and Conservative_to_Primitive two subroutine***********************
    !**********************************************************************************************

    IMPLICIT NONE
    
    CONTAINS

    !**********************************************************************************************
    !***********************************Primitive_to_Conservative**********************************
    !**********************************************************************************************
    SUBROUTINE Primitive_to_Conservative()
    !*********************This subroutine calculate the conservative variables*********************

    USE Common_Data,ONLY:id,p2,gamma,one,half,rho_1D,u_1D,p_1D,qq,factor

    IMPLICIT NONE

    INTEGER i,ii

    DO i=1,id-1
        qq(1,i)=rho_1D(i)
        qq(2,i)=rho_1D(i)*u_1D(i)
        qq(3,i)=p_1D(i)/(gamma-one)+half*rho_1D(i)*u_1D(i)*u_1D(i)
    END DO

    !*******************************The modification from Wenjia Xie*******************************
    !******If there is only one point within the numerical shock structure and epsilon < 0.4,******
    !*the modification is used, since the shock instability will also occour in the 1D computation*
    !*in such conditions.**************************************************************************
    !******Fix the mass flux of the cells behind the shock equal to that in front of the shock*****
    ii=FLOOR((id)*half)
    IF (factor==1) THEN
        qq(2,ii+1)=qq(2,ii-1)          
    END IF

    END SUBROUTINE Primitive_to_Conservative

    
    
    !**********************************************************************************************
    !***********************************Conservative_to_Primitive**********************************
    !**********************************************************************************************
    SUBROUTINE Conservative_to_Primitive()
    !***********************This subroutine calculate the primitive variables**********************

    USE Common_Data,ONLY:p2,id,gamma,one,half,rho_1D,u_1D,p_1D,qq

    IMPLICIT NONE

    INTEGER i

    DO i=1,id-1
        rho_1D(i) = qq(1,i)
        u_1D(i)   = qq(2,i)/rho_1D(i)
        p_1D(i)   = (gamma-one)*(qq(3,i)-half*rho_1D(i)*u_1D(i)*u_1D(i))
    END DO

    !*****************************Determine if the flow variable is NAN****************************
    DO i=1,id-1
        IF(rho_1D(i) /= rho_1D(i))THEN
            PAUSE
        END IF
        IF(u_1D(i) /= u_1D(i))THEN
            PAUSE
        END IF
        IF(p_1D(i) /= p_1D(i))THEN
            PAUSE
        END IF
    END DO

    END SUBROUTINE Conservative_to_Primitive
    
    !**********************************************************************************************
    !*****************************************CalculateLw_1D***************************************
    !**********************************************************************************************
    SUBROUTINE CalculateLw_1D(rhoL_,uL_,pL_,rhoR_,uR_,pR_,Lw,Rw)
    USE Common_Data,ONLY:p2,gamma,one,half,two,zero
    IMPLICIT NONE
    
    REAL(p2)::rhoL_,uL_,pL_,rhoR_,uR_,pR_
    REAL(p2)::Lw(3,3),Rw(3,3)
    REAL(p2)::rhoFace,uFace,hFace,cFace,EL,ER,hL,hR
    
    EL = pL_/(gamma-one)+half*rhoL_*uL_*uL_
    ER = pR_/(gamma-one)+half*rhoR_*uR_*uR_

    hL = (gamma/(gamma-one))*pL_/rhoL_+(uL_*uL_)/two
    hR = (gamma/(gamma-one))*pR_/rhoR_+(uR_*uR_)/two
    
    rhoFace = SQRT(rhoL_*rhoR_)
    hFace   = (SQRT(rhoL_)*hL+SQRT(rhoR_)*hR) / (SQRT(rhoL_)+SQRT(rhoR_))
    uFace   = (SQRT(rhoL_)*uL_+SQRT(rhoR_)*uR_) / (SQRT(rhoL_)+SQRT(rhoR_))
    cFace   = SQRT((gamma-one)*(hFace-half*uFace**2))
    
    Lw(1,1) = zero
    Lw(1,2) = one
    Lw(1,3) = -one/rhoFace/cFace
    Lw(2,1) = one
    Lw(2,2) = zero
    Lw(2,3) = -one/cFace/cFace
    Lw(3,1) = zero
    Lw(3,2) = one
    Lw(3,3) = one/rhoFace/cFace
    
    Rw(1,1) = -rhoFace/two/cFace
    Rw(1,2) = one
    Rw(1,3) = rhoFace/two/cFace
    Rw(2,1) = half
    Rw(2,2) = zero
    Rw(2,3) = half
    Rw(3,1) = -rhoFace*cFace/two
    Rw(3,2) = zero
    Rw(3,3) = rhoFace*cFace/two
    
    END SUBROUTINE CalculateLw_1D
    
    END MODULE Variable_Conversion_1D