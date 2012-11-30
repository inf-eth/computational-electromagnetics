#include <FDTD2DDNG.hpp>
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char * argv[])
{
	fstream Temporal2D;
	Temporal2D.open("Temporal2D.txt", ios::out);
	unsigned int NN[] = {25, 50, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
	for (int i=0; i<17; i++)
	{
		CFDTD2DDNG FDTD2DDNGSim(
								512,	// I
								512,	// J
								0,		// PMLw
								NN[i],	// MaxTime
								1,		// Snapshot resolution
								128,		// Snapshot interval
								1,		// Source choice
								1,		// Source is plane wave?
								256,	// Source location X
								64+5);	// Source location Y

		// ================== GPU Simulation ================
		FDTD2DDNGSim.StartTimer();
		FDTD2DDNGSim.CompleteRunGPU(false); // Save field snapshots to hard disk?
		FDTD2DDNGSim.StopTimer();
		cout << "Total time taken for GPU run = " << FDTD2DDNGSim.GetElapsedTime() << " seconds." << endl;
		Temporal2D << FDTD2DDNGSim.GetElapsedTime() << " ";
	}
	// ================== CPU Simulation ================
	/*FDTD1DDNGSim.StartTimer();
	FDTD1DDNGSim.CompleteRunCPU(true);
	FDTD1DDNGSim.StopTimer();
	cout << "Total time taken for CPU run = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;*/
	Temporal2D.close();
	return 0;
}
