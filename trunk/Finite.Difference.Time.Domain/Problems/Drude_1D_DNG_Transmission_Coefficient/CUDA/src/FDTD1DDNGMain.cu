#include <FDTD1DDNG.hpp>
#include <cutil.h>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	CFDTD1DDNG FDTD1DDNGSim(/*Size=*/4U*1024U, /*SourceLocation=*/10U, /*SnapshotInterval=*/16U, /*SourceChoice=*/1U);
	/*
	// =================== Initialization ===================
	FDTD1DDNGSim.StartTimer();
	if (FDTD1DDNGSim.Initialize() == 1)
		return 1;

	CUT_DEVICE_INIT(argc, argv);
	
	if (FDTD1DDNGSim.initializeFDTD1DDNGKernel() == 1)
		return 1;

	FDTD1DDNGSim.DisplaySimulationParameters();

	FDTD1DDNGSim.StopTimer();
	std::cout << "Initialization elapsed time (sec): " << FDTD1DDNGSim.GetElapsedTime() << std::endl;

	// ================== Simulation ========================
	FDTD1DDNGSim.StartTimer();
	if (FDTD1DDNGSim.runFDTD1DDNGKernels(true) == 1)
		return 1;
	FDTD1DDNGSim.StopTimer();
	std::cout << "Simulation total elapsed time (sec): " << FDTD1DDNGSim.GetElapsedTime() << std::endl;

	// ================== Clean up ======================
	FDTD1DDNGSim.StartTimer();
	if (FDTD1DDNGSim.CleanupGPU() == 1)
		return 1;
	
	FDTD1DDNGSim.StopTimer();
	std::cout << "Cleanup total elapsed time (sec): " << FDTD1DDNGSim.GetElapsedTime() << std::endl;
	*/
	/*// ================== CPU Simulation ================
	FDTD1DDNGSim.StartClock ();
	if (FDTD1DDNGSim.RunSimulationCPU (false) == 1)
		return 1;
	FDTD1DDNGSim.StopClock ();
	std::cout << "CPU Simulation time (sec): " << FDTD1DDNGSim.GetElapsedTime () << std::endl;
	*/
	FDTD1DDNGSim.StartTimer();
	cout << "Memory required for simulation = " << FDTD1DDNGSim.SimSize() << " bytes (" << (double)FDTD1DDNGSim.SimSize()/1024UL << "kB/" << (double)FDTD1DDNGSim.SimSize()/1024UL/1024UL << "MB)." << endl;
	cout << "HDD space required for data storage = " << FDTD1DDNGSim.HDDSpace() << " bytes (" << (double)FDTD1DDNGSim.HDDSpace()/1024UL << "kB/" << (double)FDTD1DDNGSim.HDDSpace()/1024UL/1024UL << "MB)." << endl;
	FDTD1DDNGSim.AllocateMemoryCPU();
	FDTD1DDNGSim.InitialiseCPU();

	// ==== GPU Simulation ====
	FDTD1DDNGSim.AllocateMemoryGPU();
	FDTD1DDNGSim.CopyDataCPUtoGPU();
	FDTD1DDNGSim.DryRunGPU();
	FDTD1DDNGSim.InitialiseExHyCPU();
	FDTD1DDNGSim.CopyExHyCPUtoGPU();
	FDTD1DDNGSim.RunSimulationGPU();

	// ==== CPU Simulation ====
	/*FDTD1DDNGSim.DryRunCPU();
	FDTD1DDNGSim.InitialiseExHyCPU();
	FDTD1DDNGSim.RunSimulationCPU();*/

	FDTD1DDNGSim.StopTimer();
	cout << "Time taken = " << FDTD1DDNGSim.GetElapsedTime() << " seconds." << endl;

	CUT_EXIT(argc, argv);
	return 0;
}
