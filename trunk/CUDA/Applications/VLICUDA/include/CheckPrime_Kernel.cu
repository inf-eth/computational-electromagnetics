#ifndef __CHECK_PRIME_KERNEL
#define __CHECK_PRIME_KERNEL

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void CheckPrime_Kernel(int A, int B)
{
	// Because of the simplicity of this tutorial, we are going to assume that
	// every block has 256 threads. Each thread simply multiplies two numbers,
	// and then stores the result.

	// The grid of blocks is 128 blocks long.

	int tid = (blockIdx.y * 128 * 256) + blockIdx.x * 256 + threadIdx.x;	// This gives every thread a unique ID.
	// By no coincidence, we'll be using this thread ID to determine which data elements to multiply.

	//pResult[tid] = pDataA[tid] * pDataB[tid];		// Each thread only multiplies one data element.
	//pResult[tid] = pDataA[tid] * pDataB[tid] / 12.34567;
	//pResult[tid] = sqrt(pDataA[tid] * pDataB[tid] / 12.34567);
	//pResult[tid] = sqrt(pDataA[tid] * pDataB[tid] / 12.34567) * sin(pDataA[tid]);

}
#endif
