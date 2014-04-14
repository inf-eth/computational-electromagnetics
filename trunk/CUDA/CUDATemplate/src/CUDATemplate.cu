#include <CUDATemplate.hpp>
#include <iostream>
#include <CUDATemplate_Kernels.cu>
#include <helper_cuda.h>
#include <helper_functions.h> // helper functions for SDK examples
using namespace std;

CUDATemplate::CUDATemplate(unsigned int pWidth, PRECISION pMultiplier): Width(pWidth), Multiplier(pMultiplier), input(NULL), output(NULL), tStart(0LL), tEnd(0LL), tDelta(0LL), tPaused(true)
{
	// Printing simulation size.
	cout << "Width      = " << Width << endl;
	cout << "Multiplier = " << Multiplier << endl;
}

// Allocate memory for data arrays.
int CUDATemplate::AllocateMemoryCPU()
{
	input = new PRECISION[Width];
	output = new PRECISION[Width];

	return 0;
}

// Initialise CPU data.
int CUDATemplate::InitialiseCPU()
{
	for (unsigned int i=0; i<Width; i++)
		input[i] = (PRECISION)(rand()%100);

	return 0;
}

int CUDATemplate::AllocateMemoryGPU()
{
	// Device memory allocation
	checkCudaErrors(cudaMalloc((void **)&d_input, sizeof(PRECISION)*Width));
	checkCudaErrors(cudaMalloc((void **)&d_output, sizeof(PRECISION)*Width));

	return 0;
}

int CUDATemplate::CopyDataCPUtoGPU()
{
	checkCudaErrors(cudaMemcpy(d_input, input, sizeof(PRECISION)*Width, cudaMemcpyHostToDevice));

	return 0;
}

int CUDATemplate::RunSimulationCPU()
{
	for (unsigned int i=0; i<Width; i++)
		output[i] = Multiplier * input[i];

	// Display results.
	cout << "Multiplier is " << Multiplier << endl;
	cout << "Input array is: " << endl;
	SafeCall(DisplayArray(Width, input), "Error: Displaying input array.");
	cout << "Output array is: " << endl;
	SafeCall(DisplayArray(Width, output), "Error: Displaying output array.");

	return 0;
}

int CUDATemplate::RunSimulationGPU()
{
	// Total local threads in a block. Can be thought of as Block dimensions.
	const unsigned int ThreadsX = 256;
	const unsigned int ThreadsY = 1;
	
	// Total blocks in simulation grid. Can be thought of as no. of blocks in grid.
	// Size should be divisible by 256.
	unsigned int BlocksX = Width/ThreadsX;
	unsigned int BlocksY = 1;

	// Kernel parameters.
	dim3 Blocks(BlocksX, BlocksY);
	dim3 Threads(ThreadsX, ThreadsY);

	cout << "Simulation (GPU) started..." << endl;

	cout << "Block dimensions: " << ThreadsX << "x" << ThreadsY << endl;
	cout << "Grid dimensions: " << BlocksX << "x" << BlocksY << endl;

	StopWatchInterface *Timer = 0;
	sdkCreateTimer(&Timer);
	sdkResetTimer(&Timer);

	sdkStartTimer(&Timer);
	// Kernel call.
	CUDATemplateKernel <ThreadsX, ThreadsY> <<<Blocks, Threads>>>(d_input, d_output, Multiplier);
	// Error checking.
	getLastCudaError("Kernel execution failed");
	sdkStopTimer(&Timer);

	checkCudaErrors(cudaMemcpy(output, d_output, sizeof(PRECISION)*Width, cudaMemcpyDeviceToHost));
	// Display results.
	cout << "Multiplier is " << Multiplier << endl;
	cout << "Input array is: " << endl;
	SafeCall(DisplayArray(Width, input), "Error: Displaying input array.");
	cout << "Output array is: " << endl;
	SafeCall(DisplayArray(Width, output), "Error: Displaying output array.");

	cout << "\r" << "kernel execution time = " << sdkGetTimerValue(&Timer) << " ms." << endl;
	sdkDeleteTimer(&Timer);

	return 0;
}

int CUDATemplate::CompleteRunCPU()
{
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");
	SafeCall(RunSimulationCPU(), "Error: Running on CPU.");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");

	return 0;
}

int CUDATemplate::CompleteRunGPU()
{
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");
	SafeCall(AllocateMemoryGPU(), "Error: Allocating memory on GPU.");
	SafeCall(CopyDataCPUtoGPU(), "Error: Copying data from CPU to GPU.");
	SafeCall(RunSimulationGPU(), "Error: Running on GPU.");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");
	SafeCall(CleanupGPU(), "Error: Cleaning up CPU.");

	return 0;
}

// Display array.
int CUDATemplate::DisplayArray(const unsigned int Size, PRECISION* Array)
{
	for (unsigned int i=0; i<Size; i++)
		cout << Array[i] << " ";
	cout << endl;

	return 0;
}

// Timing.
void CUDATemplate::StartTimer()
{
	if (tPaused == true)
	{
		tStart = GetTimeus64();
		tPaused = false;
	}
}
void CUDATemplate::StopTimer()
{
	if (tPaused == false)
	{
		tEnd = GetTimeus64();
		tDelta += tEnd - tStart;
		tStart = tEnd;
		tPaused = true;
	}
}
void CUDATemplate::ResetTimer()
{
	if (tPaused == true)
		tStart = tEnd;
	else
		tStart = GetTimeus64();

	tDelta = 0UL;
}
double CUDATemplate::GetElapsedTime()
{
	if (tPaused == false)
		tEnd = GetTimeus64();

	return ((double)(tEnd-tStart+tDelta))/(1000000.);
}

int CUDATemplate::SafeCall(int Status, const char *Error)
{
	if (Status != 0)
	{
		if (Error!=NULL) cout << Error << endl;
		exit(Status);
	}
	return Status;
}

template<typename T> void DeleteArray(T *&ptr)
{
	if (ptr != NULL)
	{
		delete[] ptr;
		ptr = NULL;
	}
}

int CUDATemplate::CleanupCPU()
{
	// Host cleanup.
	DeleteArray(input);
	DeleteArray(output);

	return 0;
}

int CUDATemplate::CleanupGPU()
{
	// Device cleanup.
	checkCudaErrors(cudaFree(d_input));
	checkCudaErrors(cudaFree(d_output));
	cudaDeviceReset();

	return 0;
}

CUDATemplate::~CUDATemplate ()
{
	// Cleanup.
	DeleteArray(input);
	DeleteArray(output);
}
