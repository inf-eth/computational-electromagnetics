#include <FDTD2D.hpp>

int main(int argc, char * argv[])
{
	CFDTD2D FDTD2DSim;
	
	// =================== Initialization ===================
	FDTD2DSim.StartClock ();	
	if (FDTD2DSim.Initialize () == 1)
		return 1;
		
	CUT_DEVICE_INIT(argc, argv);
	
	if (FDTD2DSim.initializeFDTD2DKernel () == 1)
		return 1;
	
	FDTD2DSim.DisplaySimulationParameters ();

	FDTD2DSim.StopClock ();
	std::cout << "Initialization elapsed time (sec): " << FDTD2DSim.GetElapsedTime () << std::endl;

	// ================== Simulation ========================
	FDTD2DSim.StartClock ();
	if (FDTD2DSim.runFDTD2DKernels (true) == 1)
		return 1;
	FDTD2DSim.StopClock ();
	std::cout << "Simulation total elapsed time (sec): " << FDTD2DSim.GetElapsedTime () << std::endl;

	// ================== Clean up ======================
	FDTD2DSim.StartClock ();
	if (FDTD2DSim.Cleanup () == 1)
		return 1;
	
	FDTD2DSim.StopClock ();
	std::cout << "Cleanup total elapsed time (sec): " << FDTD2DSim.GetElapsedTime () << std::endl;
	CUT_EXIT(argc, argv);
	
	return 0;
}
