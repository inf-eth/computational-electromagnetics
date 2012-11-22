#include <FDTD2DDNG.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	CFDTD2DDNG FDTD2DDNGSim(/*Size=*/4U*1024U, /*SourceLocation=*/10U, /*SnapshotInterval=*/16U, /*SourceChoice=*/1U);

	// ================== GPU Simulation ================
	FDTD2DDNGSim.StartTimer();
	FDTD2DDNGSim.CompleteRunGPU(true);
	FDTD2DDNGSim.StopTimer();
	cout << "Time taken = " << FDTD2DDNGSim.GetElapsedTime() << " seconds." << endl;
	/*
	// ================== CPU Simulation ================
	FDTD1DDNGSim.StartTimer();
	FDTD1DDNGSim.CompleteRunCPU(true);
	FDTD1DDNGSim.StopTimer();
	cout << "Time taken = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;
	*/
	return 0;
}
