#include "FDTD1DDNG.h"
#include <iostream>
using namespace std;

int main()
{
	cout << "Testing simulation..." << endl;
	CFDTD1DDNG TestSim;
	TestSim.StartTimer();
	cout << "Memory required for simulation = " << TestSim.SimSize() << " bytes (" << (double)TestSim.SimSize()/1024 << "kB/" << (double)TestSim.SimSize()/1024/1024 << "MB)." << endl;
	TestSim.AllocateMemoryCPU();
	TestSim.InitialiseCPU();
	TestSim.DryRunCPU();
	TestSim.StopTimer();
	cout << "Exiting..." << endl;

	return 0;
}
