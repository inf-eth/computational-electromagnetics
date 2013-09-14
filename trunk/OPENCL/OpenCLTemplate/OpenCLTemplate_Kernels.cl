#define PRECISION double

#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

// Kernel function.
__kernel void OpenCLTemplateKernel(__global PRECISION *input, __global PRECISION *output, const PRECISION Multiplier)
{
	unsigned int i = get_global_id(0);

	output[i] = Multiplier * input[i];
}
