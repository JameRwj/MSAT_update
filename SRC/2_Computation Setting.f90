    MODULE Computation_Setting

    !**********************************************************************************************
    !****This Module contains three subroutines: Read_Grid, Read_Settings, and Information_Write****
    !**********************************************************************************************

    IMPLICIT NONE

    CONTAINS

    !**********************************************************************************************
    !******************************************Read_Grid*******************************************
    !**********************************************************************************************
    SUBROUTINE Read_Grid()
    !***********This subroutine is used to read the grid file and compute its information**********
    USE Common_Data,ONLY:id,jd,x,y,sxx,sxy,sxMod,nxx,nxy,syx,syy,syMod,nyx,nyy,vol,half,GridName

    IMPLICIT NONE

    INTEGER :: ni,nj,i,j

    !*******************Open the grid file and read the coordinate of the node*********************
    OPEN(104,FILE=TRIM(ADJUSTL(GridName)))
    READ(104,*) ni,nj
    DO i = 1,id
        DO j = 1,jd
            READ(104,*)x(i,j),y(i,j)
        END DO
    END DO
    CLOSE(104)

    !******************************Calculate the unit normal vector********************************
    DO i = 1,id
        DO j = 1,jd-1
            sxx(i,j) = y(i,j+1)-y(i,j)
            sxy(i,j) = x(i,j)-x(i,j+1)
            sxMod(i,j)=SQRT(sxx(i,j)**2+sxy(i,j)**2)
            nxx(i,j)=sxx(i,j)/sxMod(i,j)
            nxy(i,j)=sxy(i,j)/sxMod(i,j)
        END DO
    END DO

    DO j = 1,jd
        DO i = 1,id-1
            syx(i,j) = y(i,j)-y(i+1,j)
            syy(i,j) = x(i+1,j)-x(i,j)
            syMod(i,j)=SQRT(syx(i,j)**2+syy(i,j)**2)
            nyx(i,j)=syx(i,j)/syMod(i,j)
            nyy(i,j)=syy(i,j)/syMod(i,j)
        END DO
    END DO

    !*************************Calculate the volume of the control volume***************************
    DO j=1,jd-1
        DO i=1,id-1
            vol(i,j)=half*( (x(i+1,j+1)-x(i,j)) * (y(i,j+1)-y(i+1,j))&
                           -(x(i,j+1)-x(i+1,j)) * (y(i+1,j+1)-y(i,j)) )
        END DO
    END DO

    END SUBROUTINE Read_Grid



    !**********************************************************************************************
    !****************************************Read_Settings******************************************
    !**********************************************************************************************
    SUBROUTINE Read_Settings()
    !*************************This subroutine is used to read setting file*************************
    USE Common_Data,ONLY:Limiter,Riemann_Solver,epsilon,Mainf,Initialization_Method,&
                         TestCase,Reconstruction_Method,maxCount,Reconstruction_Variables

    IMPLICIT NONE

    OPEN(103,FILE="Settings.dat")
    READ(103,*) Reconstruction_Method
    !1- MUSCL; 2-ROUND; 3-WENO5
    READ(103,*) Reconstruction_Variables
    !1- Primitive Variables; 2-Characteristic Variables
    READ(103,*) Limiter
    !Works only when the MUSCL approach is used
    !0-No limiter(First order); 1-Superbee limiter; 2-van Leer limiter;
    !3-van Albada limiter;      4-minmod limiter;   5-limiter of Deng Xi
    READ(103,*) Riemann_Solver
    !1-Roe; 2-HLLC; 3-HLL; 4-van Leer; 5-AUSMplus; 6-SLAU; 7-HLLE; 8-HLLEM
    READ(103,*) TestCase
    !1-two dimensional steady normal shock; 2-other test case
    READ(103,*) epsilon                    !Numerical shock structure
    READ(103,*) Mainf                      !Mach number
    READ(103,*) Initialization_Method
    !Initialized by Rankine-Hugoniot conditions; 2-Initialized by 1D computation
    READ(103,*) maxCount                   !Maximum iteration steps of the 1D computation
    !Only works when the 2D normal shock problem is used and the initial flow field is
    !initialized by the 1D computation
    CLOSE(103)

    END SUBROUTINE Read_Settings



    !**********************************************************************************************
    !**************************************Information_Write***************************************
    !**********************************************************************************************
    SUBROUTINE Information_Write()
    !****************This subroutine is used to write the information on the screen****************
    USE Common_Data,ONLY : Limiter,Riemann_Solver,Mainf,epsilon,id,jd,Initialization_Method,&
        TestCase,Reconstruction_Method,Reconstruction_Variables

    IMPLICIT NONE

    CHARACTER*50 :: CharTemp1,CharTemp2
    CHARACTER*100:: CharGrid,CharLimiter,CharSolver,CharInitialization,CharTestCase,CharReconstruction,CharVariable
    INTEGER :: key,error

    WRITE(CharTemp1,*)id-1
    WRITE(CharTemp2,*)jd-1
    CharGrid = TRIM(ADJUSTL(TRIM(ADJUSTL(CharTemp1))//"*"//TRIM(ADJUSTL(CharTemp2))))

    error = 0

    SELECT CASE (Reconstruction_Method)
    CASE(1)
        CharReconstruction = "MUSCL"
    CASE(2)
        CharReconstruction = "ROUND"
	CASE(3)
        CharReconstruction = "WENO5"
	CASE DEFAULT
        CharReconstruction = "Warning! The reconstruction method is wrong"
        error = 1
    END SELECT

    SELECT CASE (Reconstruction_Variables)
    CASE(1)
        CharVariable = "Primitive Variable"
    CASE(2)
        CharVariable = "Characteristic Variable"
	CASE DEFAULT
        CharReconstruction = "Warning! The reconstruction variable is wrong"
        error = 1
    END SELECT
    
    SELECT CASE (Limiter)
    CASE (0)
        CharLimiter = "No Limiter (First Order)"
    CASE (1)
        CharLimiter = "Superbee Limiter"
    CASE (2)
        CharLimiter = "van Leer Limiter"
    CASE (3)
        CharLimiter = "van Albada Limiter"
    CASE (4)
        CharLimiter = "minmod Limiter"
    CASE (5)
        CharLimiter = "Limiter of Deng Xi"
        CASE DEFAULT
        CharLimiter = "Warning! The limiter is wrong"
        error = 1
    END SELECT

    SELECT CASE (Riemann_Solver)
    CASE (1)
        CharSolver = "Roe Solver"
    CASE (2)
        CharSolver = "HLLC Solver"
    CASE (3)
        CharSolver = "HLL Solver"
    CASE (4)
        CharSolver = "van Leer Solver"
    CASE (5)
        CharSolver = "AUSMplus Solver"
    CASE (6)
        CharSolver = "SLAU Solver"
    CASE (7)
        CharSolver = "HLLE Solver"
    CASE (8)
        CharSolver = "HLLEM Solver"
        CASE DEFAULT
        CharSolver = "Warning! The solver is wrong"
        error = 1
    END SELECT

    SELECT CASE (TestCase)
    CASE(1)
        CharTestCase = "2D steady normal shock"
        SELECT CASE (Initialization_Method)
        CASE(1)
            CharInitialization = "the Rankine-Hugoniot conditions"
        CASE(2)
            CharInitialization = "the 1D computation"
            CASE DEFAULT
            CharInitialization = "Warning! The initialization method is wrong"
            error = 1
        END SELECT
    CASE(2)
        CharTestCase = "From the file in the InitialFlow floder"
        CASE DEFAULT
        CharTestCase = "Warning! The test case is wrong"
        error = 1
    END SELECT

    WRITE(*,*)
    WRITE(*,"(A80)")"===============================Computation Settings==============================="
    WRITE(*,"(A35,A60)")"The computational grid is:         " , CharGrid
    WRITE(*,*)
    WRITE(*,"(A35,A60)")"The reconstruction method is:      " , CharReconstruction
    WRITE(*,*)
    IF (Reconstruction_Method == 1) THEN
        WRITE(*,"(A35,A60)")"The used limiter is:               " , CharLimiter
        WRITE(*,*)
    END IF
    WRITE(*,"(A35,A60)")"The reconstruction variable is:    " , CharVariable
    WRITE(*,*)
    WRITE(*,"(A35,A60)")"The used Riemann solver is:        " , CharSolver
    WRITE(*,*)
    WRITE(*,"(A35,A60)")"The test case is:                  " , CharTestCase
    IF (TestCase == 1) THEN
        WRITE(*,*)
        WRITE(*,"(A35,F5.2)")"The Mach number is:                " , Mainf
        WRITE(*,*)
        WRITE(*,"(A35,F4.2)")"The numerical shock position is:   " , epsilon
        Write(*,*)
        WRITE(*,"(A35,A60)")"The flow field is initialized by:   " , CharInitialization
        WRITE(*,*)
    END IF

    IF (error==1) THEN
        WRITE(*,*)"There is at least one error in Setting.dat. Please check."
        PAUSE
        STOP
    END IF

    WRITE(*,"(A80)")"=====================Start the calculation?(Yes-1,No-else)======================"
    READ(*,*)key

    IF (key == 1) THEN
        WRITE(*,*)"The computation is ongoing, please wait......"
    ELSE
        STOP
    END IF

    END SUBROUTINE Information_Write

    END MODULE Computation_Setting