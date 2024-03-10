    MODULE Common_Data
    !**********************************************************************************************
    !********************This module contains the common data of the computation*******************
    !**********************************************************************************************

    IMPLICIT NONE

    INTEGER,PARAMETER :: p2 = SELECTED_REAL_KIND(P=15) !P=15 Double precision; p=6 is single precision

    !*****************************Information of the computational grid****************************
    INTEGER,PARAMETER :: ghostLayers = 3               !Number of ghost cells. At least 2 for MUSCL and ROUND
    CHARACTER*50 :: GridName                           !Name of the grid file
    INTEGER :: id,jd                                   !Number of nodes in two direction
    REAL(p2),ALLOCATABLE :: x(:,:),y(:,:)              !Coordinates of the node
    REAL(p2),ALLOCATABLE :: sxx(:,:),sxy(:,:)          !Face vertors
    REAL(p2),ALLOCATABLE :: syx(:,:),syy(:,:)          !Face vertors
    REAL(p2),ALLOCATABLE :: vol(:,:)                   !Volume of the control volume
    REAL(p2),ALLOCATABLE :: sxMod(:,:),syMod(:,:)      !Area of the cell face
    REAL(p2),ALLOCATABLE :: nxx(:,:),nxy(:,:)          !unit normal vector
    REAL(p2),ALLOCATABLE :: nyx(:,:),nyy(:,:)          !unit normal vector

    !**********************************Definition of the analysis**********************************
    INTEGER :: Reconstruction_Method                   !Number of the way to reconstruct the variables
    INTEGER :: Reconstruction_Variables                !The reconstructed variables
    INTEGER :: Limiter                                 !Number of the limiter function used in MUSCL
    INTEGER :: Riemann_Solver                          !Number of the Riemann solver
    INTEGER :: TestCase                                !The test case analyzed in this software
    INTEGER :: Initialization_Method                   !The way to initialize the initial flow field
    !                                                   Only work when the 2D steady normal shock problem is used
    REAL(p2):: epsilon,Mainf                           !Numerical shock structure and Mach number

    !*********************************Definition of some constants*********************************
    REAL(p2),PARAMETER ::  zero = 0.0_p2
    REAL(p2),PARAMETER ::   one = 1.0_p2
    REAL(p2),PARAMETER ::   two = 2.0_p2
    REAL(p2),PARAMETER ::  half = 0.5_p2
    REAL(p2),PARAMETER :: gamma = 1.4_p2

    !*****************************Primitive variables of the flow field****************************
    REAL(p2),ALLOCATABLE :: rho(:,:),u(:,:),v(:,:),p(:,:)

    !******************************Eigenspace of the stability matrix******************************
    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: VectorR,VectorL
    REAL(p2),DIMENSION(:)  ,ALLOCATABLE :: WR,WI

    !*******************************Parameter of the 1D computation********************************
    REAL(p2),PARAMETER:: cfl_1D=0.2_p2                !CFL
    INTEGER,PARAMETER :: nst=1000                     !Steps to output the result

    INTEGER :: factor                                 !Parameter of controling the modification of Wenjia Xie 
    !The modification is from "On numerical instabilities of Godunov-type schemes for strong shocks"
    INTEGER :: maxCount                               !Maximum steps of 1D computation
    INTEGER :: ncount                                 !Current step
    
    REAL(p2) :: stageCoefficient(3)                   !Coefficient of Runge-Kutta
    REAL(p2) :: dt                                    !Time step

    REAL(p2),ALLOCATABLE :: x_1D(:)                   !Coordinates of the node
    REAL(p2),ALLOCATABLE :: rho_1D(:),u_1D(:),p_1D(:) !Variables of the flow field
    REAL(p2),ALLOCATABLE :: DeltaX(:),deltaT(:)       !Arrays to store space and time steps
    REAL(p2),ALLOCATABLE :: qq(:,:),dqq(:,:),rhs(:,:),qq_0(:,:)
    REAL(p2),ALLOCATABLE :: residual(:,:),resid(:)

    END MODULE Common_Data