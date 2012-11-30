#include <FDTD1DDNG.hpp>
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char * argv[])
{
	fstream Temporal1D;
	Temporal1D.open("Temporal1D.txt", ios::out);
	unsigned int NN[] = {100, 250, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000};
	for (int i=0; i<17; i++)
	{
		CFDTD1DDNG FDTD1DDNGSim(/*Size=*/65536U, NN[i], /*SourceLocation=*/10U, /*SnapshotInterval=*/128U, /*SourceChoice=*/1U);

		// ================== GPU Simulation ================
		FDTD1DDNGSim.ResetTimer();
		FDTD1DDNGSim.StartTimer();
		FDTD1DDNGSim.CompleteRunGPU(false);
		FDTD1DDNGSim.StopTimer();
		cout << "Time taken = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;

		// ================== CPU Simulation ================
		/*FDTD1DDNGSim.StartTimer();
		FDTD1DDNGSim.CompleteRunCPU(true);
		FDTD1DDNGSim.StopTimer();
		cout << "Time taken = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;*/
		Temporal1D << FDTD1DDNGSim.GetElapsedTime() << " ";
	}
	Temporal1D.close();
	return 0;
}
