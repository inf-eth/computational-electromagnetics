#include <CUDATemplate.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	CFDTD1DDNG FDTD1DDNGSim(/*Size=*/4U*1024U, /*SourceLocation=*/10U, /*SnapshotInterval=*/16U, /*SourceChoice=*/1U);

	// ================== GPU Simulation ================
	FDTD1DDNGSim.StartTimer();
	FDTD1DDNGSim.CompleteRunGPU(true);
	FDTD1DDNGSim.StopTimer();
	cout << "Time taken = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;

	// ================== CPU Simulation ================
	/*FDTD1DDNGSim.StartTimer();
	FDTD1DDNGSim.CompleteRunCPU(true);
	FDTD1DDNGSim.StopTimer();
	cout << "Time taken = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;*/

	return 0;
}
