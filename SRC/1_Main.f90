    PROGRAM Main
    
    !**********************************************************************************************
    !***************************This is the main program of this software**************************
    !**********************************************************************************************
    
    USE Common_Data
    USE Computation_Setting
    USE Initialization
    USE Eigenspace
    USE Output
    
    IMPLICIT NONE
    
    !************************Read the name of the grid file from the screen************************
    WRITE(*,*)
    WRITE(*,"(A60)")"Please enter the name of the grid file (without suffix):          "
    READ(*,*)GridName
    GridName = TRIM(ADJUSTL(GridName))//".dat"
    OPEN(104,FILE=TRIM(ADJUSTL(GridName)))
    READ(104,*) id,jd               !Read the number of the grid node
    CLOSE(104)
    
    !********************************Allocate memory for variables*********************************
    ALLOCATE(u  (1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
    ALLOCATE(v  (1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
    ALLOCATE(rho(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
    ALLOCATE(p  (1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1))
    ALLOCATE(x(1:id,1:jd),y(1:id,1:jd))  
    ALLOCATE(sxx(1:id,1:jd-1),sxy(1:id,1:jd-1),syx(1:id-1,1:jd),syy(1:id-1,1:jd),vol(1:id-1,1:jd-1))
    ALLOCATE(nxx(1:id,1:jd-1),nxy(1:id,1:jd-1),nyx(1:id-1,1:jd),nyy(1:id-1,1:jd))
	ALLOCATE(sxMod(1:id,1:jd-1),syMod(1:id-1,1:jd))
    
    CALL Read_Grid()                     !Read and compute the information of the computational grid           
    
    CALL Read_Settings()                 !Read the setting of the computation
    
    CALL Information_write()             !Write the information of the computation on the screen
    
    !*******************************Choose the initialization method*******************************
    !If the test case is the 2D steady normal shock, there are two ways to initialized the initial 
    !flow field: Rankine Hugoniot condition and 1D computation. If the computation uses other test
    !cases, the initialization can only use the initial flow file
    SELECT CASE (TestCase)
    CASE(1)
        SELECT CASE (Initialization_Method)
        CASE(1)
            CALL Rankine_Hugoniot()
        CASE(2)
            CALL Computation_1D()
        END SELECT
    CASE(2)
        CALL Read_File()
    END SELECT
    
    CALL Output_Initialization()
    
    CALL Calculate_Eigen()                !Assembly the stability matrix and calculate the eigenvalue
    
    CALL Output_Matrix()                 !Output the result
    
    DEALLOCATE(rho,u,v,p)
    DEALLOCATE(x,y,sxx,sxy,syx,syy,vol,nxx,nxy,nyx,nyy,sxMod,syMod)
    DEALLOCATE(WR,WI,VectorL,VectorR)
    
    WRITE(*,*)
    WRITE(*,"(A80)")"The computation is finished, corresponding files are outputted in Results folder."
    WRITE(*,*)
    PAUSE
    
    END PROGRAM Main