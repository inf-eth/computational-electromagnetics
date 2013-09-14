#include <OpenCLTemplate.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
	COpenCLTemplate OpenCLTemplateSim(/*Size=*/256U, /*Multiplier=*/2.0);

	// ================== Simulation ================
	OpenCLTemplateSim.StartTimer();
	OpenCLTemplateSim.CompleteRun(); // Complete GPU run.
	OpenCLTemplateSim.StopTimer();
	cout << "Total time taken = " << OpenCLTemplateSim.GetElapsedTime() << " seconds." << endl;

	return 0;
}
