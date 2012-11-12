#include <FDTD1DDNG.hpp>

int main(int argc, char * argv[])
{
	CFDTD1DDNG FDTD1DDNGSim;
	
	// =================== Initialization ===================
	FDTD1DDNGSim.StartClock ();	
	if (FDTD1DDNGSim.Initialize () == 1)
		return 1;
		
	CUT_DEVICE_INIT(argc, argv);
	
	if (FDTD1DDNGSim.initializeFDTD1DDNGKernel () == 1)
		return 1;
	
	FDTD1DDNGSim.DisplaySimulationParameters ();

	FDTD1DDNGSim.StopClock ();
	std::cout << "Initialization elapsed time (sec): " << FDTD1DDNGSim.GetElapsedTime () << std::endl;

	// ================== Simulation ========================
	FDTD1DDNGSim.StartClock ();
	if (FDTD1DDNGSim.runFDTD1DDNGKernels (true) == 1)
		return 1;
	FDTD1DDNGSim.StopClock ();
	std::cout << "Simulation total elapsed time (sec): " << FDTD1DDNGSim.GetElapsedTime () << std::endl;

	// ================== Clean up ======================
	FDTD1DDNGSim.StartClock ();
	if (FDTD1DDNGSim.Cleanup () == 1)
		return 1;
	
	FDTD1DDNGSim.StopClock ();
	std::cout << "Cleanup total elapsed time (sec): " << FDTD1DDNGSim.GetElapsedTime () << std::endl;
		
	/*// ================== CPU Simulation ================
	FDTD1DDNGSim.StartClock ();
	if (FDTD1DDNGSim.RunSimulationCPU (false) == 1)
		return 1;
	FDTD1DDNGSim.StopClock ();
	std::cout << "CPU Simulation time (sec): " << FDTD1DDNGSim.GetElapsedTime () << std::endl;
	*/
	CUT_EXIT(argc, argv);
	return 0;
}
