- Run './compilenvcc' to invoke nvidia C compiler (nvcc) and generate object file(s)
- Run './link' to generate binary file called 'template'.
- Run './template' to execute the binary.

NOTE:
- Nvidia compiler is assumed to be located at /usr/local/cuda/bin/nvcc. If your installation directory is different, you need to change the 'compilenvcc' to reflect this.
- Nvidia GPU Computing SDK is assumed to be located at '/opt/NVIDIA_GPU_Compute_SDK_4.1'. If this is different then you need to change the paths in both 'compilenvcc' and 'link' files.
- If you get a 'cannot find -lcutil' error, then go to 'NVIDIA_GPU_Computing_SDK/CUDALibraries/common/lib' folder and rename libcutil_i386.a or libcutil_i386_x64 to libcutil.a. 'cp libcutil_i386 libcutil.a' works well.


