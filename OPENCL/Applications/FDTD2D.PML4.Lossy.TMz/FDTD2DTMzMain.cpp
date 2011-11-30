#include "FDTD2D.hpp"
int main(int argc, char * argv[])
{
	CFDTD2D FDTD2DSim;
	
	// =================== Initialization ===================
	FDTD2DSim.StartClock ();	
	if (FDTD2DSim.Initialize () == 1)
		return 1;
	
	if (FDTD2DSim.initializeCL () == 1)
		return 1;
	
	if (FDTD2DSim.initializeFDTD2DKernel () == 1)
		return 1;
	
	FDTD2DSim.DisplaySimulationParameters ();

	FDTD2DSim.StopClock ();
	std::cout << "Initialization elapsed time (sec): " << FDTD2DSim.GetElapsedTime () << std::endl;

	// ================== Simulation ========================
	FDTD2DSim.StartClock ();
	if (FDTD2DSim.runCLFDTD2DKernels (true) == 1)
		return 1;
	FDTD2DSim.StopClock ();
	std::cout << "Simulation total elapsed time (sec): " << FDTD2DSim.GetElapsedTime () << std::endl;

	// ================== Clean up ======================
	if (FDTD2DSim.CleanupCL () == 1)
		return 1;

	return 0;
}
