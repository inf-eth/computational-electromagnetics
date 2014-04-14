#include <CUDATemplate.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	CUDATemplate CUDATemplateSim(/*Size=*/256U, /*Multiplier=*/2.0);

	// ================== GPU Simulation ================
	CUDATemplateSim.StartTimer();
	CUDATemplateSim.CompleteRunGPU(); // Complete GPU run.
	CUDATemplateSim.StopTimer();
	cout << "Total time taken = " << CUDATemplateSim.GetElapsedTime() << " seconds." << endl;

	return 0;
}
