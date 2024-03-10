    MODULE Temporal_Discretization_1D
    !**********************************************************************************************
    !*****This Module performs the temporal discretization of the 1D computation, including two****
    !*****subroutine: CalTimeStep and Multistage_Runge_Kutta***************************************
    !**********************************************************************************************
    IMPLICIT NONE

    CONTAINS

    !**********************************************************************************************
    !******************************************CalTimeStep*****************************************
    !**********************************************************************************************
    SUBROUTINE CalTimeStep()
    !**********************This subroutine is used to calculate the time step**********************

    USE Common_Data,ONLY:p2,id,gamma,p_1D,rho_1D,u_1D,x_1D,deltaT,dt,cfl_1D,deltaX

    IMPLICIT NONE

    INTEGER ::i
    REAL(p2)::maxspeed

    DO i=1,id-1
        maxspeed = ABS(u_1D(i))+SQRT(gamma*p_1D(i)/rho_1D(i))!The maximum speed in every cell (u+c)
        deltaT(i)= cfl_1D*deltaX(i)/maxspeed                 !Calculate the time step of every cell
    END DO

    !********The time step of the computation is the minimum of the time steps in all cells********
    dt = deltaT(1)
    DO i = 1,id-1
        IF (deltaT(i) <= dt) THEN
            dt=deltaT(i)
        END IF
    END DO

    END SUBROUTINE CalTimeStep



    !**********************************************************************************************
    !************************************Multistage_Runge_Kutta************************************
    !**********************************************************************************************
    SUBROUTINE Multistage_Runge_Kutta()
    !***********This subroutine uses the multistage Runge Kutta method to solve the PDEs***********

    USE Common_Data,ONLY:id,zero,stageCoefficient,dqq,qq,dt,deltaX,qq_0,rhs,Riemann_solver
    USE Variable_Conversion_1D
    USE Riemann_Solver_1D

    IMPLICIT NONE

    INTEGER i,m

    CALL Primitive_to_Conservative()                          !Calculate the conservative veriables
                                               
                                               
    qq_0=qq                                    !Conservative variables of last step
    !                                          
    DO m = 1,3                                 
        CALL Primitive_to_Conservative()                      !Calculate the conservative veriables
                                               
        rhs=zero                               !Initialize the right-end term of PDEs

        !***********************************Choose Riemann solver**********************************
        SELECT CASE(Riemann_solver)
        CASE(1)
            CALL Roe_1D()                      !Roe
        CASE(2)                                
            CALL HLLC_1D()                     !HLLC
        CASE(3)                                
            CALL HLL_1D()                      !HLL
        CASE(4)                                
            CALL van_Leer_1D()                 !van Leer
        CASE(5)                                
            CALL AUSMplus_1D()                 !AUSM+
        CASE(6)                                
            CALL SLAU_1D()                     !SLAU
        CASE(7)                                
            CALL HLLE_1D()                     !HLLE
        CASE(8)                                
            CALL HLLEM_1D()                    !HLLEM
        END SELECT

        dqq=zero

        !********************Solve the PDEs with multistage Runge Kutta method********************
        DO i=1,id-1
            dqq(:,i)=dt/deltaX(i)*stageCoefficient(m)*rhs(:,i)
        END DO

        qq=qq_0+dqq

        CALL Conservative_to_Primitive()                         !Calculate the primitive variables
                                               
        CALL Boundary()                        !Update the boundary
    END DO

    END SUBROUTINE Multistage_Runge_Kutta

    
    
    !**********************************************************************************************
    !*******************************************Boundary*******************************************
    !**********************************************************************************************
    SUBROUTINE Boundary()
    !************************This subroutine is used to update the boundary************************
    !************************The boundary conditions are given by ghost cells**********************
    USE Common_Data,ONLY: id,gamma,mainf,p2,ghostlayers,one,rho_1D,u_1D,p_1D

    IMPLICIT NONE

    INTEGER :: i

    !*******************************Left boundary: the upstream state******************************
    DO i=0,-1,1-ghostLayers
        rho_1D(i) = one
        u_1D(i)   = one
        p_1D(i)   = one/gamma/mainf/mainf
    END DO

    !*****************************Right boundary: the downstream state*****************************
    DO i=id,id+ghostlayers-1
        rho_1D(i) = rho_1D(id-1)
        u_1D(i)   = one/rho_1D(id)
        p_1D(i)   = p_1D(id-1)
    END DO

    END SUBROUTINE Boundary

    END MODULE Temporal_Discretization_1D