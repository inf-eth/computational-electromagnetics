2D FDTD TMz simulation using C++ in a lossy medium.
- Source files are located in 'src' directory.
- On windows you can compile using either VS6 or VS2010. Go into the appropriate directory (FDTD2DTMzVS 6 or 10). Executibles are placed in respective project directories under 'debug' or 'release' folders depending on configuration.
- On linux just run compile. Output executible is placed in 'bin/debug'.
- Field data is saved in the 'FieldData' folder.
- Either run matlab script named 'FieldPlotBin.m' or 'FieldPlotBin.exe' if matlab isn't available to view the simulation.

NOTE: 'sprintf_s()' may not work in VC6, use 'sprintf()' instead.