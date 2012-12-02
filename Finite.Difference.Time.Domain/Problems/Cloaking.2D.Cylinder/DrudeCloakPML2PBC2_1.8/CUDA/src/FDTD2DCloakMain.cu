#include <FDTD2DCloak.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	CFDTD2DCloak FDTD2DCloakSim(
							512,	// I
							512,	// J
							64,		// PMLw
							8*512,	// MaxTime
							1,		// Snapshot resolution
							32,		// Snapshot interval
							2,		// Source choice
							1,		// Source is plane wave?
							256,	// Source location X
							64+5);	// Source location Y

	// ================== GPU Simulation ================
	FDTD2DCloakSim.StartTimer();
	FDTD2DCloakSim.CompleteRunGPU(true); // Save field snapshots to hard disk?
	FDTD2DCloakSim.StopTimer();
	cout << "Total time taken for GPU run = " << FDTD2DCloakSim.GetElapsedTime() << " seconds." << endl;

	// ================== CPU Simulation ================
	/*FDTD1DDNGSim.StartTimer();
	FDTD1DDNGSim.CompleteRunCPU(true);
	FDTD1DDNGSim.StopTimer();
	cout << "Total time taken for CPU run = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;*/

	return 0;
}
