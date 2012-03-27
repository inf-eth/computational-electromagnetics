#include "FDTD2D.h"
int main (int argc, char **argv)
{
	CFDTD2D FDTDSim;
	FDTDSim.StartClock ();
	FDTDSim.Initialize ();
	FDTDSim.StopClock ();
	std::cout << "Initialization Elapsed Time (sec): " << FDTDSim.GetElapsedTime () << std::endl;

	FDTDSim.StartClock ();
	FDTDSim.RunSimulation (false);
	FDTDSim.StopClock ();
	std::cout << "Simulation Elapsed Time (sec): " << FDTDSim.GetElapsedTime () << std::endl;

	return 0;
}
