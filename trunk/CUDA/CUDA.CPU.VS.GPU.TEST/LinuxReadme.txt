=============== Compilation instructions: ====================
1. Modify Makefile for dependencies. These can be additional
.cpp, .cu or any other code files.
2. ROOTDIR should point to $(NVIDIA_SDK_INSTALL_ROOT)/C/src
3. Specify libraries for linking. Check Makefiles for examples
in SDK for more info.
4. Run 'make'
==============================================================

# Example Makefile for additional src and include directories.
# Taken from FDTD3D example.
=============================================================
# Add source files here
EXECUTABLE	:= FDTD3d
# CUDA source files
CUFILES     := src/FDTD3dGPU.cu
# CUDA dependencies
CU_DEPS     := inc/FDTD3dGPUKernel.cuh
# C/C++ source files (compiled with gcc / c++)
SRCDIR		:= src/
CCFILES		:= FDTD3d.cpp FDTD3dReference.cpp
INCLUDES    += -Iinc/
=============================================================

=======================================================================================================
Manual compilation is no longer needed since Makefile is now included. Following is for reference only.
=======================================================================================================
- Run './compilenvcc' to invoke nvidia C compiler (nvcc) and generate object file(s)
- Run './link' to generate binary file called 'template'.
- Run './template' to execute the binary.

NOTE:
- Nvidia compiler is assumed to be located at /usr/local/cuda/bin/nvcc. If your installation directory
is different, you need to change the 'compilenvcc' to reflect this.
- Nvidia GPU Computing SDK is assumed to be located at '/opt/NVIDIA_GPU_Compute_SDK_4.1'. If this is 
different then you need to change the paths in both 'compilenvcc' and 'link' files.
- If you get a 'cannot find -lcutil' error, then go to 'NVIDIA_GPU_Computing_SDK/CUDALibraries/common/lib'
folder and rename libcutil_i386.a or libcutil_i386_x64 to libcutil.a. 'cp libcutil_i386 libcutil.a' works well.
=======================================================================================================

