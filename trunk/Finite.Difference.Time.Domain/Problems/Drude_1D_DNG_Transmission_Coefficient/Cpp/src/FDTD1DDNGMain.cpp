#include <FDTD1DDNG.h>
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
	fstream Temporal1D;
	Temporal1D.open("Temporal1D.txt", ios::out);
	unsigned int NN[] = {100, 250, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000};
	for (int i=0; i<17; i++)
	{
		CFDTD1DDNG TestSim(/*Size=*/65536U, NN[i], /*SourceLocation=*/10U, /*SnapshotInterval=*/128U, /*SourceChoice=*/1U);
		TestSim.ResetTimer();
		TestSim.StartTimer();
		cout << "Memory required for simulation = " << TestSim.SimSize() << " bytes (" << (double)TestSim.SimSize()/1024UL << "kB/" << (double)TestSim.SimSize()/1024UL/1024UL << "MB)." << endl;
		cout << "HDD space required for data storage = " << TestSim.HDDSpace() << " bytes (" << (double)TestSim.HDDSpace()/1024UL << "kB/" << (double)TestSim.HDDSpace()/1024UL/1024UL << "MB)." << endl;
		TestSim.AllocateMemoryCPU();
		TestSim.InitialiseCPU();
		TestSim.DryRunCPU();
		TestSim.InitialiseExHyCPU();
		TestSim.RunSimulationCPU(false);
		TestSim.StopTimer();
		cout << "Time taken = " << TestSim.GetElapsedTime() << " seconds." << endl;
		Temporal1D << TestSim.GetElapsedTime() << " ";
	}
	cout << "Exiting..." << endl;
	Temporal1D.close();
	return 0;
}
