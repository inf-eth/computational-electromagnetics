#include <FDTD2DDNG.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	CFDTD2DDNG FDTD2DDNGSim(
							512,	// I
							512,	// J
							64,		// PMLw
							4*256,	// MaxTime
							1,		// Snapshot resolution
							10,		// Snapshot interval
							1,		// Source choice
							1,		// Source is plane wave?
							256,	// Source location X
							64+5);	// Source location Y

	// ================== GPU Simulation ================
	FDTD2DDNGSim.StartTimer();
	FDTD2DDNGSim.CompleteRunGPU(true);	// Save field snapshot data to hard disk?
	FDTD2DDNGSim.StopTimer();
	cout << "Total time taken for GPU run = " << FDTD2DDNGSim.GetElapsedTime() << " seconds." << endl;

	// ================== CPU Simulation ================
	/*FDTD2DDNGSim.StartTimer();
	FDTD2DDNGSim.CompleteRunCPU(true);
	FDTD2DDNGSim.StopTimer();
	cout << "Total time taken for CPU run = " << FDTD2DDNGSim.GetElapsedTime() << " seconds." << endl;*/

	return 0;
}
