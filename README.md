# Quasi-two-dimensional (Q2D) packed-bed reactor model
A Quasi-two-dimensional code for modeling a catalytic packed-bed membrane reactor 

## Installation
Q2D model can be run using the pre-existing executables or by compiling the code on the user’s local machine from the source. All operating parameters needed to run the model are specified in the `input.dat' file.

## Pre-existing compiled versions
1. Download the OS-specific compiled versions of the Q2D model from the `run' folder.
2. Copy the `input.dat' and chemistry mechanism files in the same folder.
3. To run the model, doubleclick or type `./packbedQ2D' in the terminal.

## Compiling from the source code

#### Installation requirements
1. Compiled Cantera code (Works with Cantera 3.1)
2. Python, Boost, and Sundials - The Q2D PackBed code does not use Python or Sundials explicitly, but Cantera needs this software.
3. Microsoft Visual Studio (Required on Windows)

#### Compilation on Windows
1. Open the file 'SConstruct' from the `Compile' folder.
2. Change the directory paths to the local directory paths where Cantera and Boost suites are installed.
3. Open the MSVC command prompt in the Packbed folder and run 'scons'. The code should compile without any error. 
4. Run the command 'packbedQ2D.exe' to run the program.
5. The code can also be run by double-clicking on packbedQ2D.exe.

#### Compilation on Linux
1. Open the file 'Makefile' from the `compile' folder.
2. Change the directory paths to the local directory paths where Cantera and Boost suites are installed.
3. Open the terminal in the same folder and run 'make all'. The code should compile without any error.
4. Run the command './output' to run the program.

## Input and output files
1. The Q2D model reads input data from the file 'input.dat'. 
2. The solution is saved in output files 'molefractions_T_P_V.csv' and 'solution_T_P_V.csv'. 
3. The file 'molefractions_T_P_V.csv' stores mole fractions of gas-phase species, density, pressure, temperature, and surface coverages.
   T, P, and V correspond to the inlet temperature, pressure, and velocity(or flow rate).
4. The file 'solution_T_P_V.csv' stores the variables used in the solution vector 
   (i.e. mass fractions of gas-phase species, density, pressure, temperature, and surface coverages). 
5. The CSV file starting with "Err_" is saved in case of code failure. It saves the last successful steady-state solution.
6. The input file (input.dat) contains a parameter called 'RESTART' to specify if the user wants to use the previous solution as the initial guess. 
7. If RESTART is set to 1, then the model reads the CSV file named "savedSolution.csv" as the initial guess. The file "savedSolution.csv" gets saved every time the model is run successfully.
8. The files 'steadysolve_n_*.csv' store the solution before `n' th refinement. 

## Citation
Please cite this package along with the research articles, when used in a scholarly work.
1. `Computationally efficient and robust models of non-ideal thermodynamics, gas-phase kinetics and heterogeneous catalysis in chemical reactors', G. Kogekar, Doctoral dissertation, 2021
2. `Quasi-two-dimensional model of catalytic packed-bed membrane reactors', G. Kogekar, H. Zhu, F. Goldsmith, R.J. Kee, In preparation, 2025

## License
Q2D packed-bed model is released under the MIT license; see LICENSE for details.
