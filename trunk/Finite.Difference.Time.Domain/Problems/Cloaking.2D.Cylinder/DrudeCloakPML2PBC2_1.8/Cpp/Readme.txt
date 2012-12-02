1. Compilation
   a. Windows
   b. Linux
2. Simulation
3. Post-Processing

1. Compilation
a. Windows.
To compile on windows, open the Drude_2D_DNG_VS2010.sln file in Visual Studio 2010 (VC++ 2010 Express is recommended). Build and run the solution.
b. Linux.
Makefile with rules is provided. To compile, type 'make'. To run, simply type 'make run'. Default configuration is release and default platform is x86.
To compile/run in debug configuration add dbg=1 to compile or run command. For example, 'make dbg=1' or 'make run dbg=1'.
To compile/run a 64 bit binary, add m64bit=1 to compile or run command. For example, 'make m64bit=1' or 'make run 64bit=1'.
Configuration and platform options can be combined, e.g., 'make run dbg=1 m64bit=1'

2. Simulation
Once simulation is started, fields will be stored in 'FieldData' folder. This folder with be deleted if already present before simulation. A new folder will be created for a new simulation.
Memory and hard-disk space required will be displayed every time a simulation runs.

3. Post-Processing
A Matlab script PostProcess.m is provided that will read in simulation parameters stored in the FieldData folder, do all the post-processing and display results.
