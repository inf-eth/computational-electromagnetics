#include <CUDATemplate.hpp>

// CUDATemplate Kernel.
template <unsigned int BlockX, unsigned int BlockY> __global__ void CUDATemplateKernel(const PRECISION* input, PRECISION* output, PRECISION Multiplier)
{
	const unsigned int i = BlockX*blockIdx.x+threadIdx.x;

	output[i] = Multiplier * input[i];
}
