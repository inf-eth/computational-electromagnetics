#include <FDTD2DDNG.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	CFDTD2DDNG FDTD2DDNGSim(
							512,	// I
							512,	// J
							64,		// PMLw
							4*512,	// MaxTime
							1,		// Snapshot resolution
							10,		// Snapshot interval
							1,		// Source choice
							1,		// Source is plane wave?
							256,	// Source location X
							64+5);	// Source location Y

	// ================== GPU Simulation ================
	FDTD2DDNGSim.StartTimer();
	FDTD2DDNGSim.CompleteRunGPU(true); // Save field snapshots to hard disk?
	FDTD2DDNGSim.StopTimer();
	cout << "Total time taken for GPU run = " << FDTD2DDNGSim.GetElapsedTime() << " seconds." << endl;

	// ================== CPU Simulation ================
	/*FDTD1DDNGSim.StartTimer();
	FDTD1DDNGSim.CompleteRunCPU(true);
	FDTD1DDNGSim.StopTimer();
	cout << "Total time taken for CPU run = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;*/

	return 0;
}
