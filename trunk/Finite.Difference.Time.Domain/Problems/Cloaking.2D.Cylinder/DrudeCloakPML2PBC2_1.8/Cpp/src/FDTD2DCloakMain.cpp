#include <FDTD2DCloak.h>
#include <iostream>
using namespace std;

int main()
{
	CFDTD2DCloak FDTD2DCloakSim(
							256,	// I
							256,	// J
							64,		// PMLw
							8*256,	// MaxTime
							1,		// Snapshot resolution
							12,		// Snapshot interval
							2,		// Source choice
							1,		// Source is plane wave?
							128,	// Source location X
							64+5);	// Source location Y

	// ================== CPU Simulation ================
	FDTD2DCloakSim.StartTimer();
	FDTD2DCloakSim.CompleteRunCPU(true);
	FDTD2DCloakSim.StopTimer();
	cout << "Time taken = " << FDTD2DCloakSim.GetElapsedTime() << " seconds." << endl;

	return 0;
}
