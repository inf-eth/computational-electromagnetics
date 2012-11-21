#include <FDTD2DDNG.h>
#include <iostream>
using namespace std;

int main()
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
							50+5);	// Source location Y

	// ================== CPU Simulation ================
	FDTD2DDNGSim.StartTimer();
	FDTD2DDNGSim.CompleteRunCPU(true);
	FDTD2DDNGSim.StopTimer();
	cout << "Time taken = " << FDTD2DDNGSim.GetElapsedTime() << " seconds." << endl;

	return 0;
}
