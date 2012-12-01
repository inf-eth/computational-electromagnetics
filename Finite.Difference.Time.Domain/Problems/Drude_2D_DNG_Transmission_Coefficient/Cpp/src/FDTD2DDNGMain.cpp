#include <FDTD2DDNG.h>
#include <iostream>
using namespace std;

int main()
{
	CFDTD2DDNG FDTD2DDNGSim(
							256,	// I
							256,	// J
							64,		// PMLw
							4*256,	// MaxTime
							1,		// Snapshot resolution
							1,		// Snapshot interval
							1,		// Source choice
							1,		// Source is plane wave?
							128,	// Source location X
							64+5);	// Source location Y

	// ================== CPU Simulation ================
	FDTD2DDNGSim.StartTimer();
	FDTD2DDNGSim.CompleteRunCPU(true);
	FDTD2DDNGSim.StopTimer();
	cout << "Time taken = " << FDTD2DDNGSim.GetElapsedTime() << " seconds." << endl;

	return 0;
}
