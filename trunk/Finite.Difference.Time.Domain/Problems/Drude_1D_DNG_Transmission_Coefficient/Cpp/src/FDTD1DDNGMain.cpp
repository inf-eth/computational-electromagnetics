#include <FDTD1DDNG.h>
#include <iostream>
using namespace std;

int main()
{
	CFDTD1DDNG TestSim(/*Size=*/4U*1024U, /*SourceLocation=*/10U, /*SnapshotInterval=*/16U, /*SourceChoice=*/1U);
	TestSim.StartTimer();
	cout << "Memory required for simulation = " << TestSim.SimSize() << " bytes (" << (double)TestSim.SimSize()/1024UL << "kB/" << (double)TestSim.SimSize()/1024UL/1024UL << "MB)." << endl;
	cout << "HDD space required for data storage = " << TestSim.HDDSpace() << " bytes (" << (double)TestSim.HDDSpace()/1024UL << "kB/" << (double)TestSim.HDDSpace()/1024UL/1024UL << "MB)." << endl;
	TestSim.AllocateMemoryCPU();
	TestSim.InitialiseCPU();
	TestSim.DryRunCPU();
	TestSim.InitialiseExHyCPU();
	TestSim.RunSimulationCPU(true);
	TestSim.StopTimer();
	cout << "Time taken = " << TestSim.GetElapsedTime() << " seconds." << endl;
	cout << "Exiting..." << endl;

	return 0;
}
