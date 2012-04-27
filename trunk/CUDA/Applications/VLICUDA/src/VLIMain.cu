#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <VLI.hpp>
#include <Timer.h>
#include <iostream>

using std::cout;
using std::endl;

__global__ void CheckPrime_Kernel(int A, int B, CVLI C, bool Check)
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

int main(int argc, char * argv[])
{
	unsigned int hTimer;

	CUT_DEVICE_INIT(argc, argv);
	CUT_SAFE_CALL(cutCreateTimer(&hTimer));

	//dim3 blockGridRows(blockGridWidth, blockGridHeight);
	//dim3 threadBlockRows(256, 1);

	dim3 blockGridRows(4096, 4096, 1);
	dim3 threadBlockRows(16, 16, 1);

	CUT_SAFE_CALL(cutResetTimer(hTimer));
	CUT_SAFE_CALL(cutStartTimer(hTimer));

	int A=0, B=0;
	bool Check = false;
	CVLI C;
	__int64 tStart, tEnd;

	tStart = GetTimeus64();
	for (unsigned int i=0; i<2; i++)
	{
		CheckPrime_Kernel<<<blockGridRows, threadBlockRows>>>(A, B, C, Check);
		CUT_CHECK_ERROR("multiplyNumbersGPU() execution failed\n");
		CUDA_SAFE_CALL(cudaThreadSynchronize());
	}
	tEnd = GetTimeus64();
	cout << "Time taken for kernel execution: " << ((double)(tEnd-tStart))/(1000000.) << " seconds." << endl;

	CUT_SAFE_CALL(cutStopTimer(hTimer));
	double gpuTime = cutGetTimerValue(hTimer);
	cout << "GPU time: " << gpuTime*0.001 << " seconds." << endl;

	CUT_SAFE_CALL(cutDeleteTimer(hTimer));
	CUT_EXIT(argc, argv);
}
