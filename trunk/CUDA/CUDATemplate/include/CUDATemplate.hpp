#ifndef CUDATemplate_H_
#define CUDATemplate_H_

#define PRECISION double

#include <Timer.h>

class CUDATemplate
{
private:
	// Size of problem domain
	const unsigned int Width;

	// Host data arrays.
	PRECISION* input;
	PRECISION* output;

	// Scalar.
	const PRECISION Multiplier;

	// Device data arrays.
	PRECISION* d_input;
	PRECISION* d_output;

	// Timer variables.
	__int64 tStart;
	__int64 tEnd;
	__int64 tDelta;
	bool tPaused;

public:
	CUDATemplate(unsigned int=256U, const PRECISION=2.0);	// Default width of problem is 256 and multiplier is 2.

	// Memory allocation and initialisation.
	int AllocateMemoryCPU();
	int InitialiseCPU();
	int AllocateMemoryGPU();
	int CopyDataCPUtoGPU();

	// Simulations.
	int RunSimulationGPU(); // RunCUDAKernels
	int RunSimulationCPU();

	// Complete run encapsulating all the sub-functions.
	int CompleteRunCPU();
	int CompleteRunGPU();

	// Timing.
	void StartTimer();
	void StopTimer();
	void ResetTimer();
	PRECISION GetElapsedTime();

	int DisplayArray(const unsigned int, PRECISION*);

	int SafeCall(int, const char []=NULL);

	int CleanupCPU();
	int CleanupGPU();
	~CUDATemplate();
};
#endif // #ifndef CUDATemplate_H_
