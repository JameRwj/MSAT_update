# MSAT_update

This version of MSAT provides an extension of [MSAT ](https://www.sciencedirect.com/science/article/pii/S2352711023002625)by extending the reconstruction stencil from three points to five points. In this way, the updated MSAT can quantitatively analyze shock instability problems for higher-order schemes while still maintaining the functionality of the initial version of the MSAT.

## Authors

- Weijie Ren <rwj@nudt.edu.cn>

- Wenjia Xie <xiewenjia@nudt.edu.cn>

- Ye Zhang <zhangye17a@nudt.edu.cn>

- Hang Yu <yuhang18@gfkd.edu.cn>

- Zhengyu Tian <tianzhengyu_kd@163.com>

## Usage

Following the next step to use MSAT:

1.  Prepare the Fortran complier and LAPACK on the computer.

2.  Clone the directory with <https://github.com/JameRwj/MSAT_SoftwareX.git>.

3.  Generate the computational grid and export to _Grid.dat._

4.  Determine the computation setting in _Settings.dat._

5.  Compile and run the MSAT.

6.  Enter the name of grid file (without extension).

7.  View the calculation results and analyze the processing data.

## Instructions for the compilation

MSAT is written in **Fortran** and is tested by three compilers:[ Intel Visual Studio](https://visualstudio.microsoft.com/zh-hans/), [Ifortran (Intel® Fortran Compiler)](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/fortran-compiler.html), and [gfortran (GNU Fortran)](https://gcc.gnu.org/fortran/). Here are some instructions for the compilation with different compilers.

1.  MSAT is tested by Visual Studio on the Windows system. The Visual Studio (Visual Studio 2022) must be equipped with [Intel® oneAPI Toolkits](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) (2023.2) (**Intel oneAPI Base Toolkit** for Math Kernel Library and **Intel oneAPI HPC Toolkit** for Intel Fortran Compiler). Before performing the compilation of MSAT, users can create a new project and add ".f90" type files in the source code to this project. Then, some necessary settings must be made:

    - The platform is set to **x64**

    - Open the Intel Math Kernel Library (Set **Parallel** in _Project properties/Fortran/Libraries_)

    Finally, MSAT can be compiled with Visual Studio.

2.  MSAT is also tested by Ifortran and gfortran compiler on the Linux system. The compilation process is consistent regardless of which compiler is used. **Ifortran**/**gfortran** and **LAPACK** must be installed before compiling. **Cmake** and **make** tools make the compilation process simple. To use these two tools, users need to prepare \*CMakeLists.txt \*just as follows:

    ```markup
    cmake_minimum_required(VERSION 3.16)

    set(MAKE_Fortran_COMPILER   /usr/bin/gfortran)

    enable_language(Fortran)

    Project(MSAT)

    file(GLOB_RECURSE SRC_FILES ./*.f90)

    add_executable(MSAT ${SRC_FILES})
    find_package(LAPACK)

    target_link_libraries(MSAT LAPACK::LAPACK)
    ```

    The first line declares the minimum version of cmake and the second line declares the path of the compiler. If there are no special needs, the rest can be left unchanged. Then, MSAT can be compiled by the commands "cmake.", "make", "./MSAT" in turn.

## Settings.dat

The main settings of the MSAT are stored in _Settings.dat_, including:

- The reconstruction method used in MAST: 1. MUSCL; 2. ROUND; 3. fifth-order WENO.

- The reconstruction variables: 1. primitive variables; 2. characteristic variables.

- The limiter function used in the MUSCL approach. The correspondence between the number and the limiter function is: 1. Superbee limiter; 2. van Leer limiter; 3. van Albada limiter; 4. minmod limiter; 5. the limiter proposed by [Xi Deng](https://doi.org/10.1016/j.jcp.2023.112052) ; And it will be first-order accurate if choose 0.

- The Riemann solver used in MSAT. The correspondence between the number and the Riemann solver is: 1. Roe solver; 2. HLLC solver; 3. HLL solver; 4. van Leer solver; 5. AUSM+ solver; 6. SLAU solver; 7. HLLE solver; 8. HLLEM solver.

- Test case: 1. 2D normal shock 2. other test cases

- Numerical shock structure $\varepsilon$ of the 2D normal shock, which is between 0 and 1.

- Mach number of the 2D normal shock.

- The method of Initializing the 2D normal shock. If it is 1, the 2D flow field is initialized by the Rankine-Hugoniot conditions. And the 2D flow field will be initialized by projecting the steady flow field from 1D computation onto the 2D domain.

- The iteration steps of the 1D computation.

## Grid.dat

The coordinates of the grid nodes are stored in the grid file, which is in the ".dat" type. Note that the origin of the coordinates is the lower-left corner of the grid. There are two parts in the grid file. The first part specifies the number of grid nodes in the _x_ and _y_ dimensions. The second part is divided into three columns: the first column contains the _x_ coordinates of each grid point, the second column contains the _y_ coordinates, and the third column displays the _z_ values, which should be 0 due to the two-dimensional nature of the computation.

## InitialFlow

These files contains the initial flow field, including _InitialFlow_rho.dat_, _InitialFlow_u.dat_, _InitialFlow_v.dat_, and _InitialFlow_p.dat_. Those files will be used to initialize the flow field if other test cases are used. Note that other programs should obtain the four files.

## Results

The results of MSAT include the flow field, result and residual of the 1D computation (if have), scatters of all eigenvalues, and the eigenvectors of the most unstable eigenvalue if it exceeds 0. All the results are ".plt" type and can be opened by Tecplot.
